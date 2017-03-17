afarray tensorMulAtBA(afarray A, afarray B) {

  int n = B.dims(0);
  int p = B.dims(2);
  int r = A.dims(1);

  afarray AtBA = moddims(B, dim4(n, n * p));
  AtBA = matmul(transpose(A), AtBA);

  AtBA = moddims(AtBA, dim4(r, n, p));

  {
    afarray temp = constant(0, dim4(n, r, p));
    gfor(seq i, p) { temp(span, span, i) = transpose(AtBA(span, span, i)); }
    AtBA = temp;
  }

  AtBA = moddims(AtBA, dim4(n, r * p));
  AtBA = matmul(transpose(A), AtBA);
  AtBA = moddims(AtBA, dim4(r, r, p));

  return AtBA;
}