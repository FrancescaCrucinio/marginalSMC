library(dlm)
library(ggplot2)
library(resample)
library(ggpubr)
set.seed(1234)

# data
Time.step <- 10
sigmaX <- 1
sigmaY <- 1
x <- rep(0, times = Time.step)
y <- rep(0, times = Time.step)
x[1] <- rnorm(n = 1, 0, sigmaX)
y[1] <- rnorm(n = 1, x[1], sigmaY)
for(i in 2:Time.step){
  x[i] <- rnorm(n = 1, x[i-1], sigmaX)
  y[i] <- rnorm(n = 1, x[i], sigmaY)
}
# Kalman filter
mod <- dlm(FF = 1, V = sigmaY, GG = 1, W = sigmaX, m0 = 0, C0 = sigmaX)
res <- dlmFilter(y, mod, debug = FALSE, simplify = FALSE)

true_mean <- res$m[Time.step+1]
true_variance <- c(res$U.C[[Time.step+1]]^2 * res$D.C[Time.step+1,]^2)
ub <- true_mean + sqrt(true_variance)
lb <- true_mean - sqrt(true_variance)
true_probability <- pnorm(ub, mean = true_mean, sd = sqrt(true_variance)) -
  pnorm(lb, mean = true_mean, sd = sqrt(true_variance))


Nparticles <- c(100, 500, 1000, 5000, 10000)
Nrep <- 100

mean_fa_apf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
mean_apf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
mean_mapf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)

indicator_fa_apf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
indicator_apf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
indicator_mapf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)

for (n in 1:length(Nparticles)) {
  print(n)
  for (i in 1:Nrep) {
    randomness <- rnorm(Nparticles[n])
    uniforms_resampling <- sort(runif(Nparticles[n]))
    # FA-APF
    x_past <- rep(0, times = Nparticles[n])
    opt_var <- (1/sigmaX + 1/sigmaY)^(-1)
    opt_mean <- (y[1]/sigmaY)*opt_var*rep(1, times = Nparticles[n])
    x_fa_apf <- opt_mean + sqrt(opt_var)*randomness
    W <- W_apf(x_fa_apf, y[1], sigmaY, x_past, sigmaX, opt_mean, opt_var)
    twisting <- dnorm(x_fa_apf, mean = y[2], sd = sqrt(sigmaX + sigmaY))
    ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
    x_fa_apf <- x_fa_apf[ancestors]
    # APF
    x_past <- rep(0, times = Nparticles[n])
    proposal_var <- 2*opt_var
    proposal_mean <- (y[1]/sigmaY)*opt_var*rep(1, times = Nparticles[n])
    x_apf <- proposal_mean + sqrt(proposal_var)*randomness
    W <- W_apf(x_apf, y[1], sigmaY, x_past, sigmaX, proposal_mean, proposal_var)
    twisting <- optimal_auxiliary_function(x_apf, y[2], proposal_mean, sigmaX, sigmaY, proposal_var)
    ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
    x_apf <- x_apf[ancestors]
    # MAPF
    x_past <- rep(0, times = Nparticles[n])
    proposal_var <- 2*opt_var
    proposal_mean <- (y[1]/sigmaY)*opt_var*rep(1, times = Nparticles[n])
    x_mapf <- proposal_mean + sqrt(proposal_var)*randomness
    W <- W_apf(x_mapf, y[1], sigmaY, x_past, sigmaX, proposal_mean, proposal_var)
    twisting <- optimal_auxiliary_function(x_mapf, y[2], proposal_mean, sigmaX, sigmaY, proposal_var)
    W_twisted <- W*twisting/sum(W*twisting)
    ancestors <- multinomial_resample(W_twisted, Nparticles[n], uniforms_resampling)
    x_mapf_resampled <- x_mapf[ancestors]
    for (t in 2:(Time.step-1)) {
      randomness <- rnorm(Nparticles[n])
      uniforms_resampling <- sort(runif(Nparticles[n]))
      # FA-APF
      x_past <- x_fa_apf
      opt_mean <- (x_fa_apf/sigmaX + y[t]/sigmaY)*opt_var
      x_fa_apf <- opt_mean + sqrt(opt_var)*randomness
      W <- W_apf(x_fa_apf, y[t], sigmaY, x_past, sigmaX, opt_mean, opt_var)
      twisting <- dnorm(x_fa_apf, mean = y[t+1], sd = sqrt(sigmaX + sigmaY))/
        dnorm(x_past, mean = y[t], sd = sqrt(sigmaX + sigmaY))
      ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
      x_fa_apf <- x_fa_apf[ancestors]
      # APF
      x_past <- x_apf
      proposal_mean <- (x_apf/sigmaX + y[t]/sigmaY)*proposal_var
      x_apf <- proposal_mean + sqrt(proposal_var)*randomness
      proposal_mean_future <- (x_apf/sigmaX + y[t+1]/sigmaY)*proposal_var
      W <- W_apf(x_apf, y[t], sigmaY, x_past, sigmaX, proposal_mean, proposal_var)
      twisting <- optimal_auxiliary_function(x_apf, y[t+1], proposal_mean_future, sigmaX, sigmaY, proposal_var)/
        optimal_auxiliary_function(x_past, y[t], proposal_mean, sigmaX, sigmaY, proposal_var)
      ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
      x_apf <- x_apf[ancestors]
      # MAPF
      x_past <- x_mapf
      proposal_mean <- (x_mapf_resampled/sigmaX + y[t]/sigmaY)*opt_var
      x_mapf <- proposal_mean + sqrt(proposal_var)*randomness
      proposal_mean_future <- (x_mapf/sigmaX + y[t+1]/sigmaY)*opt_var
      W_twisted <- W_mapf(x_mapf, x_past, y[t], y[t+1], proposal_mean, proposal_mean_future, sigmaX, sigmaY, proposal_var, W_twisted)
      ancestors <- multinomial_resample(W_twisted, Nparticles[n], uniforms_resampling)
      x_mapf_resampled <- x_mapf[ancestors]
    }
    # last time step
    # FA-APF
    randomness <- rnorm(Nparticles[n])
    x_past <- x_fa_apf
    opt_mean <- (x_fa_apf/sigmaX + y[Time.step]/sigmaY)*opt_var
    x_fa_apf <- opt_mean + sqrt(opt_var)*randomness
    W <- W_apf(x_fa_apf, y[Time.step], sigmaY, x_past, sigmaX, opt_mean, opt_var)
    twisting <- 1/dnorm(x_past, mean = y[Time.step], sd = sqrt(sigmaX + sigmaY))
    ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
    x_fa_apf <- x_fa_apf[ancestors]
    # APF
    x_past <- x_apf
    proposal_mean <- (x_apf/sigmaX + y[Time.step]/sigmaY)*proposal_var
    x_apf <- proposal_mean + sqrt(proposal_var)*randomness
    W <- W_apf(x_past, y[Time.step], sigmaY, x_past, sigmaX, proposal_mean, proposal_var)
    twisting <- 1/optimal_auxiliary_function(x_past, y[Time.step], proposal_mean, sigmaX, sigmaY, proposal_var)
    ancestors <- multinomial_resample(W*twisting/sum(W*twisting), Nparticles[n], uniforms_resampling)
    x_apf <- x_apf[ancestors]
    # MAPF
    x_past <- x_mapf
    proposal_mean <- (x_mapf_resampled/sigmaX + y[Time.step]/sigmaY)*opt_var
    x_mapf <- proposal_mean + sqrt(proposal_var)*randomness
    W_twisted <- W_mapf_last(x_mapf, x_past, y[Time.step], proposal_mean, sigmaX, sigmaY, proposal_var, W_twisted)
    ancestors <- multinomial_resample(W_twisted, Nparticles[n], uniforms_resampling)
    x_mapf_resampled <- x_mapf[ancestors]

    mean_fa_apf[i, n] <- mean(x_fa_apf)
    mean_apf[i, n] <- mean(x_apf)
    mean_mapf[i, n] <- mean(x_mapf_resampled)

    indicator_fa_apf[i, n] <- sum((x_fa_apf <= ub) & (x_fa_apf >= lb))/Nparticles[n]
    indicator_apf[i, n] <- sum((x_apf <= ub) & (x_apf >= lb))/Nparticles[n]
    indicator_mapf[i, n] <- sum((x_mapf_resampled <= ub) & (x_mapf_resampled >= lb))/Nparticles[n]
  }
}

cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
# means
df_mean <- data.frame("var" = c(colVars(mean_fa_apf), colVars(mean_apf), colVars(mean_mapf)))
df_mean$algo <- rep(c("FA-APF", "APF", "MAPF"), each = length(Nparticles))
df_mean$algo <- factor(df_mean$algo, levels = c("FA-APF", "MAPF", "APF"))
df_mean$N <- rep(Nparticles, times = 3)

ggplot(data = df_mean, aes(x = N, y = var, group = algo, fill = algo, colour = algo)) +
  geom_line(aes(linetype = algo), size = 2) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none")
# ggsave("lgssm_auxiliary_mean.pdf", width = 12, height = 8, dpi = 300)

# indicator
df_indicator <- data.frame("var" = c(colVars(indicator_fa_apf), colVars(indicator_apf), colVars(indicator_mapf)))
df_indicator$algo <- rep(c("FA-APF", "APF", "MAPF"), each = length(Nparticles))
df_indicator$algo <- factor(df_mean$algo, levels = c("FA-APF", "MAPF", "APF"))
df_indicator$N <- rep(Nparticles, times = 3)

ggplot(data = df_indicator, aes(x = N, y = var, group = algo, fill = algo, colour = algo)) +
  geom_line(aes(linetype = algo), size = 2) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none")
# ggsave("lgssm_auxiliary_indicator.pdf", width = 12, height = 8, dpi = 300)


p <- ggplot(data = df_indicator, aes(x = N, y = var, group = algo, fill = algo, colour = algo)) +
  geom_line(aes(linetype = algo), size = 2) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=10),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="bottom",
        legend.key.size = unit(2, 'cm'))
my_legend <- get_legend(p)
legend_p <- as_ggplot(my_legend)
# ggsave("lgssm_auxiliary_legend.pdf", width = 10, height = 0.8, dpi = 300)
