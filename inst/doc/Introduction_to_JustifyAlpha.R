## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#if(!require(JustifyAlpha)){devtools::install_github("Lakens/JustifyAlpha")}
library(JustifyAlpha)
library(pwr)
library(ggplot2)


## ----minimize, fig.width = 6--------------------------------------------------
res1 <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power",
              error = "minimize",
              costT1T2 = 1,
              priorH1H0 = 1,
              verbose = FALSE, 
              printplot = TRUE)

res1$alpha
res1$beta
res1$errorrate

## ----echo=T, results='hide', fig.show='hide'----------------------------------
# Note that printing plots is suppressed with rmarkdown here and it is simply used to generate the data.
resplot1 <- optimal_alpha(power_function = "pwr::pwr.t.test(d = 1, n = 3, sig.level = x, type = 'two.sample', alternative = 'two.sided')$power", error = "minimize", costT1T2 = 1, printplot = TRUE)
resplot2 <- optimal_alpha(power_function = "pwr::pwr.t.test(d = 1, n = 10, sig.level = x, type = 'two.sample', alternative = 'two.sided')$power", error = "minimize", costT1T2 = 1, printplot = TRUE)
resplot3 <- optimal_alpha(power_function = "pwr::pwr.t.test(d = 1, n = 30, sig.level = x, type = 'two.sample', alternative = 'two.sided')$power", error = "minimize", costT1T2 = 1, printplot = TRUE)

## ----mudgeplot, fig.width = 6-------------------------------------------------
plot_data <- rbind(resplot1$plot_data, resplot2$plot_data, resplot3$plot_data)
plot_data$n <- as.factor(rep(c(3, 10, 30), each = 9999))

w_c_alpha_plot <- ggplot(data=plot_data, aes(x=alpha_list, y=w_c_list)) +
  geom_line(size = 1.3, aes(linetype = n)) +
  geom_point(aes(x = resplot1$alpha, y = (1 * resplot1$alpha + 1 * (resplot1$beta)) / (1 + 1)), color="red", size = 3) +
  geom_point(aes(x = resplot2$alpha, y = (1 * resplot2$alpha + 1 * (resplot2$beta)) / (1 + 1)), color="red", size = 3) +
  geom_point(aes(x = resplot3$alpha, y = (1 * resplot3$alpha + 1 * (resplot3$beta)) / (1 + 1)), color="red", size = 3) +
  theme_minimal(base_size = 16) +
  scale_x_continuous("alpha", seq(0,1,0.1)) +
  scale_y_continuous("weighted combined error rate", seq(0,1,0.1), limits = c(0,1))
w_c_alpha_plot

## ----balance, fig.width = 6---------------------------------------------------
res2 <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power",
              error = "balance",
              costT1T2 = 1,
              priorH1H0 = 1,
              verbose = FALSE, 
              printplot = TRUE)

res2$alpha
res2$beta
res2$errorrate

## ----cost minimize, fig.width = 6---------------------------------------------
res3 <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power",
              error = "minimize",
              costT1T2 = 4,
              priorH1H0 = 1,
              verbose = FALSE, 
              printplot = TRUE)

res3$alpha
res3$beta
res3$errorrate

## ----cost balance, fig.width = 6----------------------------------------------
res4 <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power",
              error = "balance",
              costT1T2 = 4,
              priorH1H0 = 1,
              verbose = FALSE, 
              printplot = TRUE)

res4$alpha
res4$beta
res4$errorrate

## ----prior, fig.width = 6-----------------------------------------------------
res5 <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power",
              error = "minimize",
              costT1T2 = 1,
              priorH1H0 = 4,
              verbose = FALSE, 
              printplot = TRUE)

res5$alpha
res5$beta
res5$errorrate

## ----optimal_sample-----------------------------------------------------------
res6 <- optimal_sample(power_function = "pwr.t.test(d=0.5, n = sample_n, sig.level = x, type='two.sample', alternative='two.sided')$power", 
                       error = "minimize", 
                       errorgoal = 0.05, 
                       costT1T2 = 1,
                       priorH1H0 = 1)

res6

## ----eval = FALSE-------------------------------------------------------------
# res <- optimal_alpha(power_function = "TOSTER::powerTOSTtwo(alpha=x, N=200, low_eqbound_d=-0.4, high_eqbound_d=0.4) ",
#               error = "minimize",
#               costT1T2 = 1,
#               priorH1H0 = 1)
# 
# res$alpha
# res$beta
# res$errorrate
# res$plot

## ----eval = FALSE-------------------------------------------------------------
# res <- optimal_alpha(power_function = "Superpower::ANOVA_exact( (Superpower::ANOVA_design(design = '2b', n = 64, mu = c(0, 0.5), sd = 1, plot = FALSE)), alpha_level = x, verbose = FALSE)$main_results$power/100",
#               error = "minimize",
#               costT1T2 = 1,
#               priorH1H0 = 1)
# 
# res$alpha
# res$beta
# res$errorrate

## ----p-plot, fig.width=7, fig.height=5, echo=FALSE, message=FALSE, warning=FALSE, results='hide', fig.cap="*P*-value distributions for a two-sided independent *t*-test with N = 150 and d = 0.5 (black curve) or d = 0 (horizontal dashed line) which illustrates how *p*-values just below 0.05 can be more likely when there is no effect than when there is an effect."----
N <- 150
d <- 0.5
p <- 0.05
p_upper <- 0.05 + 0.00000001
p_lower <- 0.00000001
ymax <- 10 # Maximum value y-scale (only for p-curve)

# Calculations
se <- sqrt(2 / N) # standard error
ncp <- (d * sqrt(N / 2)) # Calculate non-centrality parameter d

# p-value function
pdf2_t <- function(p) 0.5 * dt(qt(p / 2, 2 * N - 2, 0), 2 * N - 2, ncp) / dt(qt(p / 2, 2 * N - 2, 0), 2 * N - 2, 0) + dt(qt(1 - p / 2, 2 * N - 2, 0), 2 * N - 2, ncp) / dt(qt(1 - p / 2, 2 * N - 2, 0), 2 * N - 2, 0)
plot(-10,
  xlab = "P-value", ylab = "", axes = FALSE,
  main = "P-value distribution for d = 0.5 and N = 150", xlim = c(0, .1), ylim = c(0, ymax), cex.lab = 1.5, cex.main = 1.5, cex.sub = 1.5
)
axis(side = 1, at = seq(0, 1, 0.01), labels = seq(0, 1, 0.01), cex.axis = 1.5)
# cord.x <- c(p_lower,seq(p_lower,p_upper,0.001),p_upper)
# cord.y <- c(0,pdf2_t(seq(p_lower, p_upper, 0.001)),0)
# polygon(cord.x,cord.y,col=rgb(0.5, 0.5, 0.5,0.5))
curve(pdf2_t, 0, .1, n = 1000, col = "black", lwd = 3, add = TRUE)
ncp <- (0 * sqrt(N / 2)) # Calculate non-centrality parameter d
curve(pdf2_t, 0, 1, n = 1000, col = "black", lty = 2, lwd = 3, add = TRUE)
abline(v = 0.05, col = "black", lty = 3, lwd = 3)


## ----p-bf, fig.width=7, fig.height=7------------------------------------------

n1 <- 100
n2 <- 100

loops <- seq(from = 0, to = 3, by = 0.001)
p <- numeric(length(loops))
bf <- numeric(length(loops))
d <- numeric(length(loops))
tval <- numeric(length(loops))
i <- 0

for(t in loops){
  i <- i+1
  bf[i] <- exp(BayesFactor::ttest.tstat(t, n1, n2)$bf)
  p[i] <- 2*pt(t, ((n1 + n2) - 2), lower=FALSE)
  tval[i] <- t
  d[i] <- t * sqrt((1/n1)+(1/n2))
}

plot(p, bf, type="l", lty=1, lwd=2, log = "y")
abline(v = seq(0,1,0.1), h = c(0, 1/10, 1/3, 1, 3, 10), col = "gray", lty = 1)


## -----------------------------------------------------------------------------
p[which.min(abs(bf-1))]

## -----------------------------------------------------------------------------
res8 <- ttestEvidence("lindley", 100, 100)[[1]]
res8

