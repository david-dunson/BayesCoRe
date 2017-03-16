class stiefelTensor
{
public:
  afarray A;
  int r, n, p;
  afarray posSet, zeroSet, upperSet, diagonalSet;
  afarray U;
  afarray D;
  afarray Z;
  afarray tau, tau_copy;

  stiefelTensor(afarray _A, int _r)
  {
    A = _A;
    n = A.dims(0);
    p = A.dims(2);
    r = _r;
    initialize();
  }

  void initialize()
  {
    // make the diagonal 2 and upper matrix 3, so that they don't get updated
    afarray upper3 = upper(constant(3, dim4(n, n)));
    for (int i = 0; i < p; i++)
    {
      // make upper matrix 3
      A.slice(i) += upper3;
      // make diagonal 2
      // afarray diag_adjust = -diag(A.slice(i)) + 2;
      // A.slice(i) += diag(diag_adjust, 0, false);
      A.slice(i) = setDiag(A.slice(i), constant(2, n));
    }
    A(A >= 3) = 3;

    // location index set
    posSet = A == 1;
    zeroSet = A == 0;
    diagonalSet = A == 2;
    upperSet = A == 3;

    // variance matrix
    tau = constant(1000.0, dim4(r, r));
    tau_copy = tile(tau, dim4(1, 1, p));

    // initialize stiefel factor
    afarray avgA = sum(A, 2) / (float)p;
    afarray Di, Vt;
    svd(U, Di, Vt, avgA);
    U = U.cols(0, r - 1);

    // latent probit variale
    afarray m = A - 0.5;
    afarray zeros = constant(0, A.dims());

    Z = A - 0.5;

    Z(posSet) = gpu_rtruncnorm(m(posSet), 1.0, zeros(posSet), true);
    Z(zeroSet) = gpu_rtruncnorm(m(zeroSet), 1.0, zeros(zeroSet), false);
    for (int i = 0; i < p; i++)
    {
      Z.slice(i) = setDiag(Z.slice(i), randn(n) * sqrt(2));
    }
    Z = symmetrize(Z);

    updateD();
    cout << "here" << endl;
    updateZ();
  }

  // symmetrize matrix in a tensor
  afarray symmetrize(const afarray X)
  {
    afarray Y = X;
    // for (int i = 0; i < X.dims(2); i++) {
    gfor(seq i, X.dims(2))
    {
      afarray L = lower(X(span, span, i));
      afarray U = upper(X(span, span, i));
      Y(span, span, i) = X(span, span, i) - U + transpose(L);
    }
    return Y;
  }

  // set diagonal value to matrix (matrix only)
  afarray
  setDiag(const afarray X, const afarray diagvec)
  {
    afarray Y = X;
    afarray diag_adjust = -diag(X) + diagvec;
    Y = X + diag(diag_adjust, 0, false);
    return Y;
  }

  void updateD()
  {

    afarray UZU = moddims(Z, dim4(n, n * p));
    UZU = matmul(transpose(U), UZU);
    UZU = moddims(UZU, dim4(r, n, p));

    {
      afarray temp = constant(0, dim4(n, r, p));
      gfor(seq i, p)
      {
        temp(span, span, i) = transpose(UZU(span, span, i));
      }
      UZU = temp;
    }

    UZU = moddims(UZU, dim4(n, r * p));
    UZU = matmul(transpose(U), UZU);
    UZU = moddims(UZU, dim4(r, r, p));

    // afarray UZU = constant(0, dim4(r, r, p));

    // for (int i = 0; i < p; i++)
    // {
    //   UZU.slice(i) = matmul(transpose(U), Z.slice(i), U);
    // }

    // off diagonal part
    afarray v = 1.0 / (0.5 + 1.0 / tau_copy);
    afarray m = v * (UZU / 2.0);

    // diagonal part
    for (int i = 0; i < p; i++)
    {

      afarray uzu_d = diag(UZU.slice(i));
      afarray v_d = 1.0 / (1.0 + 1.0 / diag(tau));
      afarray m_d = v_d * uzu_d;

      v.slice(i) = setDiag(v.slice(i), v_d);
      m.slice(i) = setDiag(m.slice(i), m_d);
    }

    D = randn(dim4(r, r, p)) * sqrt(v) + m;
    D = symmetrize(D);
  }

  void updateZ()
  {
    afarray zeros = constant(0, dim4(n, n));

    cout << "wha" << endl;

    for (int i = 0; i < p; i++)
    {

      afarray Z_local = Z.slice(i);
      afarray UDU = matmul(U, D.slice(i), transpose(U));

      afarray posSetLocal = posSet.slice(i);
      afarray zeroSetLocal = zeroSet.slice(i);

      cout << "wha1" << endl;

      Z_local(posSetLocal) = gpu_rtruncnorm(UDU(posSetLocal), 1.0, zeros(posSetLocal), true);
      Z_local(zeroSetLocal) = gpu_rtruncnorm(UDU(zeroSetLocal), 1.0, zeros(zeroSetLocal), false);

      Z_local.eval();
      cout << "wha2" << endl;

      Z_local = setDiag(Z_local, randn(n) * sqrt(2) + diag(UDU));
      Z.slice(i) = Z_local;
    }
    cout << "wha" << endl;
    Z = symmetrize(Z);
  }

  void Run(int steps)
  {
    for (int j = 0; j < steps; j++)
    {
      updateD();
      updateZ();

      cout << j << endl;
    }
  }
};
