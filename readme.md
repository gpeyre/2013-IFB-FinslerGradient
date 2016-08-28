This Matlab package contains the source code to reproduce the figure of the articles:

G. Nardi, G. Peyré, F-X. Vialard. [Geodesics on Shape Spaces with Bounded Variation and Sobolev Metrics](https://hal.archives-ouvertes.fr/hal-00952672/). SIAM Journal on Imaging Sciences, 9(1), pp. 238–274, 2016.

G. Charpiat, G. Nardi, G. Peyré, F-X. Vialard. [Finsler Steepest Descent with Applications to Piecewise-regular Curve Evolution](https://hal.archives-ouvertes.fr/hal-00849885/). Interfaces and Free Boundaries, 18(1), pp. 1–44, 2016.

![Example of piecewise-rigid flow](pr-flow.png)

* load_curve.m loads the different curves used for the tests.

* load_rigidification_operator.m defines all operators associated with a given curve in order to calculate the Finsler gradient. This function defines the tangent and normal vectors to the curve and all penalties introduced in the paper.

* compute_rkhs_energy.m takes as input a given target curve  and an initial curve  and  calculates the energy  intrduced in the paper. The two parameters needed for such a function are the standard deviation of the Gaussian kernels,  σ and δ.

* perform_rigidification_linprog.m takes as input the current curve and the gradient of the energy at such a curve and calculates the Finsler gradient by solving a linear programming with Mosek.  Two models are implemented :
    -	Finsler gradient associated with the penalty for piecewise rigid motions which is influenced by the parameter ρ;
    -	Finsler gradient associated with the penalty for piecewise similarity motions which is influenced by the parameters ρ and λ .

* test_man.m, test_horse.m, and test_harm.m reproduce the Finsler evolutions of the paper (figure 2 and 5).
In each of these function we can choose options.gradtype as ‘sobolev’ or ‘finsler’ in order to perform a Sobolev or Finsler evolution. In the case of the Finsler evolution the parameter option.constraints allows one to choose between the piecewise rigid Finsler gradient (options.constraints=1) and the piecewise similarity Finsler gradient (options.conatrints=2).

* compare_parameter_rho.m, compare_parameter_lambda.m, and  compare_gradients.m reproduce the evolutions in figure 1, 3, and 4, respectively.

Copyright (c) 2014 G. Nardi, G. Peyré
