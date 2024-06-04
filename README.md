# AbunNB

The R package **AbunNB** implements the penalized empirical likelihood (EL) inference for abundance via EM algorithms. Under the zero-truncated negative binomial regression model, this package provides three main functions to make statistical inferences.

+ abun_nb(): used to calculate the maximum (penalized) EL estimators of the abundance, regression coefficients, and the dispersion parameter.

- abun_nb_ci(): used to calculate the (penalized) EL ratio confidence interval of abundance.

+ summary(): used to summarize the maximum (penalized) EL estimation results.


# Usage

In the **R** software, the following codes are used to install the package:

install.packages("devtools")

library(devtools)

install_github("ecnuliuyang/AbunNB")



As an example, the following codes are provided to show the usage of the main functions in this packag: 

library(AbunNB)

example(abun_nb_ci)


# Reference
Ji, Y., and Liu, Y. (2024). A penalized empirical likelihood approach for estimating population sizes under the zero-truncated negative binomial regression model. *Biometrical Journal*, submitted.


#

For questions, comments or remarks about the code please contact Y. Liu at this email address <liuyangecnu@163.com>.
