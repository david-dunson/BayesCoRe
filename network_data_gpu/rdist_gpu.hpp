template<class T>
afarray gpu_rexp(dim4 dim, T alpha) {
	return -log(randu(dim)) / alpha;
}

afarray normalize_prob2(afarray a) {
	int n = a.dims(2);
	return a / tile(sum(a, 2), dim4(1, 1, n));
}

afarray normalize_prob3(afarray a) {
	int n = a.dims(3);
	return a / tile(sum(a, 3), dim4(1, 1, 1, n));
}

afarray softmax2(afarray a) {
	int n = a.dims(2);
	afarray b = a - tile(max(a, 2), dim4(1, 1, n));
	afarray exp_b = exp(b);
	return exp_b;
}

afarray softmax3(afarray a) {
	int n = a.dims(3);
	afarray b = a - tile(max(a, 3), dim4(1, 1, 1, n));
	afarray exp_b = exp(b);
	return exp_b;
}

afarray gpu_multinomial2(afarray prob, int n1, int n2) {
	afarray accu_prob = accum(prob, 2);
	int n = prob.dims(2);

	accu_prob = tile(randu(n1, n2), dim4(1, 1, n)) < accu_prob;
	// 0 1 2 2  4    <- 1.5 => 3 , which would be wrong
	afarray C = constant(NaN, n1, n2);
	for (int kappa = 0; kappa < n; ++kappa) {
		afarray select = (accu_prob(span, span, kappa) && isNaN(C));
		if (anyTrue<bool>(select))
			C(select) = kappa;
	}

	return C.as(u32);
}

afarray gpu_multinomial3(afarray prob, int n1, int n2 ,int p) {
	afarray accu_prob = accum(prob, 3);
	int n = prob.dims(3);

	accu_prob = tile(randu(n1, n2, p), dim4(1, 1, 1, n)) < accu_prob;
	afarray C = constant(NaN, n1, n2, p);
	for (int kappa = 0; kappa < n; ++kappa) {
		afarray select = (accu_prob(span, span, span, kappa) && isNaN(C));
		if (anyTrue<bool>(select))
			C(select) = kappa;
	}

	return C.as(u32);
}

afarray gpu_rand_int(dim4 dim, int k) {
	return floor(randu(dim) * k);
}

afarray gpu_rbern(afarray prob) {
	return (randu(prob.dims()) < prob).as(u32);
}

afarray gpu_pnorm(afarray q) {
	return 0.5f + erf(q / sqrt(2.0f)) / 2;
}

afarray gpu_rtruncated_std_lb_aux(dim4 n, afarray alpha, afarray lb) {
	afarray z = gpu_rexp(n, alpha) + lb;
	afarray condition = (lb < alpha);
	afarray delta = condition.as(f32) * exp(-(alpha - z) * (alpha - z) / 2.0)
			+ (!condition).as(f32)
					* exp(
							-(alpha - z) * (alpha - z) / 2.0
									+ (lb - alpha) * (lb - alpha) / 2.0);
//	delta.eval();
	afarray u = randu(n);
//	u.eval();
	afarray not_satisfied = u > delta;
	int n1 = sum<int>(not_satisfied);
	if (n1 > 0) {
		z(not_satisfied) = gpu_rtruncated_std_lb_aux(n1, alpha(not_satisfied),
				lb(not_satisfied));
	}
	return z;
}

afarray gpu_rtruncnorm_std_lb(dim4 dim, afarray lb) {
	afarray alpha = (lb + sqrt(lb * lb + 4)) / 2.0;
	afarray delta = constant(1, dim);
	afarray u = constant(1.5, dim);
	afarray z = constant(0, dim);
	afarray not_satisfied = u > delta;
//	int n = sum<int>(not_satisfied);
	z = gpu_rtruncated_std_lb_aux(dim, alpha, lb);
	return z;
}

afarray gpu_rtruncnorm_std_ub(dim4 dim, afarray ub) {
	return -gpu_rtruncnorm_std_lb(dim, -ub);
}

afarray gpu_rtruncnorm(afarray mean, float sigma, afarray bound, bool lower_bound) {
	afarray bound1 = (bound - mean) / sigma;
	afarray choice1 = bound1 > 6;
	if (anyTrue<bool>(choice1))
		bound1(choice1) = 6;
	afarray choice2 = bound1 < -6;
	if (anyTrue<bool>(choice2))
		bound1(choice2) = -6;
	dim4 n = mean.dims();
	afarray r(n);
	if (lower_bound) {
		r = gpu_rtruncnorm_std_lb(n, bound1);
	} else {
		r = gpu_rtruncnorm_std_ub(n, bound1);
	}
	assert(anyTrue<bool>(!isNaN(r)));
	return r * sigma + mean;

}
