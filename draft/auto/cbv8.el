(TeX-add-style-hook
 "cbv8"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt" "fleqn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("natbib" "round") ("ulem" "normalem")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "color"
    "natbib"
    "rotating"
    "graphicx"
    "subcaption"
    "float"
    "bbm"
    "amsthm"
    "amsmath"
    "amssymb"
    "mathrsfs"
    "ulem"
    "tikz"
    "algorithm"
    "algpseudocode"
    "array")
   (TeX-add-symbols
    '("KL" 2)
    '("mc" 1)
    '("bb" 1)
    '("alex" 1)
    '("leo" 1)
    "xbeta"
    "xtheta"
    "sgamma"
    "core"
    "be"
    "ee"
    "Binom"
    "No"
    "PG"
    "IG"
    "Ga"
    "Bern"
    "U"
    "Poi"
    "NB"
    "cov"
    "var"
    "diag"
    "Diag"
    "1"
    "bigO"
    "dt"
    "mass"
    "hess")
   (LaTeX-add-labels
    "eq:sharp"
    "EQ:Rel_Dens_Motivation"
    "eq:DD"
    "eq:ID"
    "fig:two_distances"
    "EQ:relaxedDensityPosMeasure"
    "fig:gaussian_inequality"
    "SEC:Zero_Measure_Methods"
    "EQ:relaxedDensityZeroMeasure"
    "fig:torus"
    "SEC:Positive_measure_theory"
    "EQ:Expectation_Positive_Measure_Constraint"
    "EQ:Expectation_Positive_Measure_Relaxed"
    "THM:positive_measure_approximation_error"
    "THM:Positive_measure_convergence_rate"
    "SEC:Zero_measure_theory"
    "TABLE:Equality_constraints_examples"
    "THM:RCP_construction"
    "EQ:Constrained_rcp"
    "THM:Relaxed_Expectation_Convergence_Measure_Zero"
    "hamiltonian"
    "leap-frog"
    "eq:hessian_extrinsic"
    "sphere_examples"
    "ordered_dp_prior"
    "dirichlet"
    "network_model_basis"
    "network_model"
    "APP:Positive_Measure_Convergence_Proofs"
    "EQ:Expectation_Identity_Positive_Measure"
    "table_circle")
   (LaTeX-add-bibliographies
    "reference")
   (LaTeX-add-amsthm-newtheorems
    "theorem"
    "lemma"
    "corollary"
    "remark"
    "example"))
 :latex)

