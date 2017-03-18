afarray qnorm(afarray p, afarray mu, afarray sigma) {
  bool lower_tail = true;
  bool log_p = false;

  afarray r, val;

  val = p;

  // R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
  val(p == 0) = -af::Inf;
  val(p == 1) = af::Inf;
  val(sigma < 0) = af::NaN;
  val(sigma == 0) = mu(sigma == 0);

  afarray q0 = p - 0.5;

  afarray fabsq425 = abs(q0) <= .425;

  if (anyTrue<bool>(fabsq425)) { /* 0.075 <= p <= 0.925 */
    // cout << "t" << endl;

    afarray q = q0(fabsq425);
    r = .180625 - q * q;
    val(fabsq425) =
        q *
        (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r +
              67265.770927008700853) *
                 r +
             45921.953931549871457) *
                r +
            13731.693765509461125) *
               r +
           1971.5909503065514427) *
              r +
          133.14166789178437745) *
             r +
         3.387132872796366608) /
        (((((((r * 5226.495278852854561 + 28729.085735721942674) * r +
              39307.89580009271061) *
                 r +
             21213.794301586595867) *
                r +
            5394.1960214247511077) *
               r +
           687.1870074920579083) *
              r +
          42.313330701600911252) *
             r +
         1.);
    // cout << "t" << endl;
  }

  if (anyTrue<bool>(!fabsq425)) { /* closer than 0.075 from {0,1} boundary */
    // cout << "f" << endl;

    afarray q = q0(!fabsq425);
    afarray pL = p(!fabsq425);
    /* r = min(p, 1-p) < 0.075 */
    r = q;
    if (anyTrue<bool>(q > 0)) r(q > 0) = 1 - pL(q > 0); /* 1-p */
    if (anyTrue<bool>(q < 0)) r(q < 0) = pL(q < 0);     /* = R_DT_Iv(p) ^=  p */

    r = sqrt(-log(r));

    afarray valNOTfabsq425 = val(!fabsq425);

    /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

    afarray rL5bool = r <= 5.;
    if (anyTrue<bool>(
            rL5bool)) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */

      afarray rL5 = r(rL5bool);
      rL5 += -1.6;
      valNOTfabsq425(rL5bool) =
          (((((((rL5 * 7.7454501427834140764e-4 + .0227238449892691845833) *
                    rL5 +
                .24178072517745061177) *
                   rL5 +
               1.27045825245236838258) *
                  rL5 +
              3.64784832476320460504) *
                 rL5 +
             5.7694972214606914055) *
                rL5 +
            4.6303378461565452959) *
               rL5 +
           1.42343711074968357734) /
          (((((((rL5 * 1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    rL5 +
                .0151986665636164571966) *
                   rL5 +
               .14810397642748007459) *
                  rL5 +
              .68976733498510000455) *
                 rL5 +
             1.6763848301838038494) *
                rL5 +
            2.05319162663775882187) *
               rL5 +
           1.);
    }
    if (anyTrue<bool>(!rL5bool)) { /* very close to  0 or 1 */
      afarray rL5 = r(!rL5bool);
      rL5 += -5.;
      valNOTfabsq425(!rL5bool) =
          (((((((rL5 * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) *
                    rL5 +
                .0012426609473880784386) *
                   rL5 +
               .026532189526576123093) *
                  rL5 +
              .29656057182850489123) *
                 rL5 +
             1.7848265399172913358) *
                rL5 +
            5.4637849111641143699) *
               rL5 +
           6.6579046435011037772) /
          (((((((rL5 * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    rL5 +
                1.8463183175100546818e-5) *
                   rL5 +
               7.868691311456132591e-4) *
                  rL5 +
              .0148753612908506148525) *
                 rL5 +
             .13692988092273580531) *
                rL5 +
            .59983220655588793769) *
               rL5 +
           1.);
    }

    if (anyTrue<bool>(q < 0)) valNOTfabsq425(q < 0) *= -1;

    val(!fabsq425) = valNOTfabsq425;

    /* return (q >= 0.)? r : -r ;*/
    // cout << "f" << endl;
  }

  return mu + sigma * val;
}

afarray qnorm(afarray p, afarray mu) {
  return qnorm(p, mu, constant(1, p.dims()));
}

afarray qnorm(afarray p) {
  return qnorm(p, constant(0, p.dims()), constant(1, p.dims()));
}