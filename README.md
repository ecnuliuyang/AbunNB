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

example(abun_nb)



# Case study 

To reproduce analyses in the case study in Ji and Liu (2024), we present the corresponding code as follows. 

+ Code to summarize the frequencie of the number of captures and sex information (Table 3 in Ji and Liu (2024))
  
library(AbunNB)
dat <- data.frame(y = rowSums(blackbear[,-9]), sex = blackbear[,9])
tab3 <- aggregate(dat$y, by = list(dat$sex), table)$x
names(tab3) <- c("Male", "Female")
tab3


- Code to plot the subfigure A of Figure 1 in Ji and Liu (2024)
  
fy <- table(dat$y)
y <- as.numeric(names(fy))
r <- y[-1] * fy[-1]/fy[-length(fy)]

par(mar = c(4.8,5.2,2.4,2.4))
plot(y[-length(y)], r, type = "p", xlab = "y",
     ylab = expression((y+1)*f[y+1]/f[y]), cex.lab = 1.2, yaxt = "n",
     main = "A: Ratio plot")
fit <- lm(r ~ y[-length(y)])
abline(fit, lty = 2, lwd = 1.5)
axis(2, seq(2,14,2))


+ Code to calculate the maximum (penalized) empirical likelihood estimates and standard errors of abundance under the negative binomial regression model
mele <- abun_nb(numCap=dat$y, x=dat$sex)
mpele <- abun_nb(numCap=dat$y, x=dat$sex, method = "PEL")
round(mele@N)
round(abun_nb_se(mele)$se_N,1)
round(mpele@N)
round(abun_nb_se(mpele)$se_N,1)


- Code to calculate the Chao lower bound
  
round(mpele@Nchao)


+ Code to calculate the ratio regression estimate with weights (3.2) in Rocchetti et al. (2011)
  
fit <- lm(log(r)~y[-length(y)], weights = 1/(1/fy[-length(fy)] + 1/fy[-1]))
f0_ratioReg <- fy[1]*exp(-fit$coefficients[1])
N_ratioReg <- length(dat$y) + f0_ratioReg
round(N_ratioReg)
### Standard error; see (2.3) in Rocchetti et al. (2011; AoAS)
varN <- length(dat$y)*f0_ratioReg/N_ratioReg + exp(-2*fit$coefficients[1]) *
  fy[1] * (fy[1]*summary(fit)$coefficients[1,2]^2+1)
round(sqrt(varN),1)



- Code to calculate the maximum penalized empirical likelihood estimate and standard error of the dispersion parameter
round(mpele@k,1)
round(abun_nb_se(mpele)$se_k,1)


+ Code to calculate the (penalized) empirical likelihood ratio confidence intervals of abundance under the negative binomial regression model
round(abun_nb_ci(mele))
round(abun_nb_ci(mpele))


- Code to plot the subfigure B of Figure 1 in Ji and Liu (2024)
Ns <- seq(50, 550, length = 200)
elr <- NULL
pelr <- NULL
for (i in 1:200) {
  elr <- c(elr, -2*( abun_nb(dat$y, dat$sex, N0 = Ns[i])@loglikelihood -
                     mele@loglikelihood))
  pelr <- c(pelr,- 2*( abun_nb(dat$y, dat$sex, N0 = Ns[i],
                               method = "PEL")@loglikelihood -
                         mpele@loglikelihood))
}

par(mar = c(4.8,5.2,2.4,2.4))
plot(Ns, elr, xlab = "N", ylab = "Log-likelihood ratio function",
     type="l", cex.lab = 1.2, lwd = 1.5, main = "B: Log-EL ratio curve")
lines(Ns,pelr, col=2, lty=2, lwd = 1.5)
legend("bottomright", c("With penalty", "Without penalty"), col = c(2, 1),
       lty = c(2,1), cex = 1, lwd = 1.5)




# Reference
Ji, Y., and Liu, Y. (2024). A penalized empirical likelihood approach for estimating population sizes under the negative binomial regression model. *Journal of Computational and Applied Mathematics*, submitted.

Rocchetti, I., Bunge, J., and B$\ddot o$hning, D. (2011). Populationsize estimation based upon ratios of recapture probabilities. *Annals of Applied Statistics*, **5**, 1512–1533.

#

For questions, comments or remarks about the code please contact Y. Liu at this email address <liuyangecnu@163.com>.
