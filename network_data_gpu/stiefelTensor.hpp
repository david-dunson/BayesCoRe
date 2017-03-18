#include "tensorOps.hpp"

class stiefelTensor {
 public:
  afarray A;
  int r, n, p;
  afarray posSet, zeroSet;
  afarray U;
  afarray D;
  afarray Z;
  afarray tau;

  bool fastMode;

  std::vector<float> trace_D;

  stiefelTensor(afarray _A, int _r, int _fastMode = 0) {
    A = _A.as(u8);
    n = A.dims(0);
    p = A.dims(2);
    r = _r;

    fastMode = _fastMode == 1;

    initialize();
  }

  void initialize() {
    // initialize stiefel factor
    {
      afarray avgA = (sum(A, 2)).as(f32) / (float)p;
      afarray Di, Vt;

      avgA(avgA > (1 - 1E-5)) = (1 - 1E-5);
      avgA(avgA < 1E-5) = 1E-5;

      avgA = qnorm(avgA);
      avgA.eval();

      if (anyTrue<bool>(isNaN(avgA))) {
        af_print(avgA(isNaN(avgA)));
        exit(0);
      }

      // cout << (min(min(avgA))).scalar<float>() << endl;
      // cout << (max(max(avgA))).scalar<float>() << endl;

      svd(U, Di, Vt, avgA);

      U = U.cols(0, r - 1);
    }

    // latent probit variale

    afarray m = A.as(f32) - 0.5;
    afarray zeros = constant(0, A.dims());
    Z = A.as(f32) - 0.5;

    // make the diagonal 2 and upper matrix 3, so that they don't get updated
    {
      afarray upper3 = upper(constant(3, dim4(n, n)));
      for (int i = 0; i < p; i++) {
        // make upper matrix 3
        A.slice(i) += upper3;
        // make diagonal 2
        // afarray diag_adjust = -diag(A.slice(i)) + 2;
        // A.slice(i) += diag(diag_adjust, 0, false);
        A.slice(i) = setDiag(A.slice(i), constant(2, n));
      }
      af::sync();
      A(A >= 3) = 3;

      // location index set
      posSet = A == 1;
      zeroSet = A == 0;
      // diagonalSet = A == 2;
      // upperSet = A == 3;
    }
    // variance matrix
    tau = constant(1000.0, dim4(r, r));
    // tau_copy = tile(tau, dim4(1, 1, p));

    Z = gpu_rtruncnorm_zero(m, posSet, zeroSet);

    for (int i = 0; i < p; i++) {
      Z.slice(i) = setDiag(Z.slice(i), randn(n) * sqrt(2));
    }
    af::sync();
    Z = symmetrize(Z);

    updateD();
    updateZ();

    // prepare trace vector
    trace_D.resize(0);
  }

  // symmetrize matrix in a tensor
  afarray symmetrize(const afarray X) {
    afarray Y = X;
    // for (int i = 0; i < X.dims(2); i++) {
    gfor(seq i, X.dims(2)) {
      afarray L = lower(X(span, span, i));
      afarray U = upper(X(span, span, i));
      Y(span, span, i) = X(span, span, i) - U + transpose(L);
    }
    af::sync();

    return Y;
  }

  // set diagonal value to matrix (matrix only)
  afarray setDiag(const afarray X, const afarray diagvec) {
    afarray Y = X;
    afarray diag_adjust = -diag(X) + diagvec;
    Y = X + diag(diag_adjust, 0, false);
    return Y;
  }

  void updateD() {
    afarray UZU = tensorMulAtBA(U, Z);

    // for (int i = 0; i < p; i++)
    // {
    //   UZU.slice(i) = matmul(transpose(U), Z.slice(i), U);
    // }

    // off diagonal part
    afarray v = 1.0 / (0.5 + 1.0 / tau);
    afarray m = tile(v, dim4(1, 1, p)) * (UZU / 2.0);

    afarray v_d = 1.0 / (1.0 + 1.0 / diag(tau));
    v = setDiag(v, v_d);
    // diagonal part
    for (int i = 0; i < p; i++) {
      afarray uzu_d = diag(UZU.slice(i));
      afarray m_d = v_d * uzu_d;
      m.slice(i) = setDiag(m.slice(i), m_d);
    }

    af::sync();

    D = randn(dim4(r, r, p)) * tile(sqrt(v), dim4(1, 1, p)) + m;
    D = symmetrize(D);
    af::sync();

    if (anyTrue<bool>(isNaN(D))) {
      cout << "D has NaN" << endl;
      exit(0);
    }
  }

  void updateZ() {
    afarray zeros = constant(0, dim4(n, n));

    // ver 1: large memory
    /*
    afarray UDU = tensorMulAtBA(transpose(U), D);

    Z = gpu_rtruncnorm_zero(UDU, posSet, zeroSet);

    for (int i = 0; i < p; i++) {
      Z.slice(i) = setDiag(Z.slice(i), randn(n) * sqrt(2));
    }
    Z = symmetrize(Z);
    */

    // ver2 memory efficient
    for (int i = 0; i < p; i++) {
      afarray Z_local = Z.slice(i);
      afarray UDU = matmul(U, D.slice(i), transpose(U));

      afarray posSetLocal = posSet.slice(i);
      afarray zeroSetLocal = zeroSet.slice(i);

      Z_local = gpu_rtruncnorm_zero(UDU, posSetLocal, zeroSetLocal);
      // Z_local(posSetLocal) =
      //     gpu_rtruncnorm(UDU(posSetLocal), 1.0, zeros(posSetLocal), true);
      // Z_local(zeroSetLocal) =
      //     gpu_rtruncnorm(UDU(zeroSetLocal), 1.0, zeros(zeroSetLocal),
      // false);
      // Z_local.eval();

      Z_local = setDiag(Z_local, randn(n) * sqrt(2) + diag(UDU));
      Z.slice(i) = Z_local;
    }
    af::sync();

    Z = symmetrize(Z);
    af::sync();

    if (anyTrue<bool>(isNaN(Z))) {
      cout << "Z has NaN" << endl;

      exit(0);
    }
  }

  void updateZfast() {
    afarray zeros = constant(0, dim4(n, n));

    // ver 1: large memory, run fast

    afarray UDU = tensorMulAtBA(transpose(U), D);
    af::sync();

    int batchSize = 15;
    int batch = p / batchSize;
    int remain = p % batchSize;
    for (int b = 0; b < (batch + (remain > 0)); b++) {
      int startIdx = b * batchSize;
      int endIdx = b * batchSize + batchSize - 1;
      if (endIdx > (p - 1)) endIdx = p - 1;
      Z.slices(startIdx, endIdx) = gpu_rtruncnorm_zero(
          UDU.slices(startIdx, endIdx), posSet.slices(startIdx, endIdx),
          zeroSet.slices(startIdx, endIdx));
    }

    while (anyTrue<bool>(isNaN(Z))) {
      afarray nanSet = isNaN(Z);

      auto nanCountDev = sum(sum(sum(nanSet.as(f32))));

      int nanCount = (int)nanCountDev.scalar<float>();
      cout << "nan count: " << nanCount << endl;
      if (nanCount < 10) {
        Z(nanSet) = 0;
      } else {
        // need fix
        cout << "h" << endl;
        Z(nanSet) =
            gpu_rtruncnorm_zero(UDU(nanSet), A(nanSet) == 1, A(nanSet) == 0);
      }
    }

    for (int i = 0; i < p; i++) {
      // Z.slices(i) =
      // gpu_rtruncnorm_zero(UDU.slice(i), posSet.slice(i), zeroSet.slice(i));

      Z.slice(i) = setDiag(Z.slice(i), randn(n) * sqrt(2));
    }

    af::sync();

    Z = symmetrize(Z);
    af::sync();

    if (anyTrue<bool>(isNaN(Z))) {
      cout << "Z has NaN" << endl;

      exit(0);
    }
  }

  void Run(int steps) {
    af_device_gc();
    for (int j = 0; j < steps; j++) {
      updateD();

      if (fastMode) {
        updateZfast();
      } else {
        updateZ();
      }
      std::vector<float> curD = conv_to_std_vector<float>(lower(D));
      trace_D.insert(trace_D.end(), curD.begin(), curD.end());

      cout << j << endl;
    }
  }
};
