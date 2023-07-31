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

mean_bpf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
mean_pf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
mean_mpf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)

indicator_bpf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
indicator_pf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)
indicator_mpf <- matrix(0, ncol = length(Nparticles), nrow = Nrep)

for (n in 1:length(Nparticles)) {
  print(n)
  for (i in 1:Nrep) {
    randomness <- rnorm(Nparticles[n])
    uniforms_resampling <- sort(runif(Nparticles[n]))
    # BPF
    x_bpf <- sqrt(sigmaX)*randomness
    lW <- lW_bpf(x_bpf, y[1], sigmaY)
    ancestors <- multinomial_resample(weights(lW), Nparticles[n], uniforms_resampling)
    x_bpf <- x_bpf[ancestors]
    # PF
    x_past <- rep(0, times = Nparticles[n])
    opt_var <- (1/sigmaX + 1/sigmaY)^(-1)
    opt_mean <- (y[1]/sigmaY)*opt_var
    x_pf <- opt_mean + sqrt(opt_var)*randomness
    lW <- lW_pf(x_pf, y[1], sigmaY, x_past, opt_mean, opt_var)
    ancestors <- multinomial_resample(weights(lW), Nparticles[n], uniforms_resampling)
    x_pf <- x_pf[ancestors]
    # MPF
    x_past <- rep(0, times = Nparticles[n])
    opt_var <- (1/sigmaX + 1/sigmaY)^(-1)
    opt_mean <- (y[1]/sigmaY)*opt_var
    x_mpf <- opt_mean + sqrt(opt_var)*randomness
    lW <- lW_pf(x_mpf, y[1], sigmaY, x_past, opt_mean, opt_var)
    W_mpf <- weights(lW)
    ancestors <- multinomial_resample(W_mpf, Nparticles[n], uniforms_resampling)
    x_mpf_resampled <- x_mpf[ancestors]
    for (t in 2:Time.step) {
      randomness <- rnorm(Nparticles[n])
      uniforms_resampling <- sort(runif(Nparticles[n]))
      # BPF
      x_bpf <- x_bpf + sqrt(sigmaX)*randomness
      lW <- lW_bpf(x_bpf, y[t], sigmaY)
      ancestors <- multinomial_resample(weights(lW), Nparticles[n], uniforms_resampling)
      x_bpf <- x_bpf[ancestors]
      # PF
      x_past <- x_pf
      opt_var <- (1/sigmaX + 1/sigmaY)^(-1)
      opt_mean <- (x_pf/sigmaX + y[t]/sigmaY)*opt_var
      x_pf <- opt_mean + sqrt(opt_var)*randomness
      lW <- lW_pf(x_pf, y[t], sigmaY, x_past, opt_mean, opt_var)
      ancestors <- multinomial_resample(weights(lW), Nparticles[n], uniforms_resampling)
      x_pf <- x_pf[ancestors]
      # MPF
      x_past <- x_mpf
      opt_var <- (1/sigmaX + 1/sigmaY)^(-1)
      opt_mean <- (x_mpf_resampled/sigmaX + y[t]/sigmaY)*opt_var
      x_mpf <- opt_mean + sqrt(opt_var)*randomness
      lW <- lW_mpf(x_mpf, y[t], sigmaY, x_past, opt_mean, opt_var, W_mpf)
      W_mpf <- weights(lW)
      ancestors <- multinomial_resample(W_mpf, Nparticles[n], uniforms_resampling)
      x_mpf_resampled <- x_mpf[ancestors]
    }
    mean_bpf[i, n] <- mean(x_bpf)
    mean_pf[i, n] <- mean(x_pf)
    mean_mpf[i, n] <- mean(x_mpf_resampled)

    indicator_bpf[i, n] <- sum((x_bpf <= ub) & (x_bpf >= lb))/Nparticles[n]
    indicator_pf[i, n] <- sum((x_pf <= ub) & (x_pf >= lb))/Nparticles[n]
    indicator_mpf[i, n] <- sum((x_mpf_resampled <= ub) & (x_mpf_resampled >= lb))/Nparticles[n]
  }
}

cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
# means
df_mean <- data.frame("var" = c(colVars(mean_bpf), colVars(mean_pf), colVars(mean_mpf)))
df_mean$algo <- rep(c("BPF", "PF", "MPF"), each = length(Nparticles))
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
# ggsave("lgssm_simple_mean.pdf", width = 12, height = 8, dpi = 300)

# indicator
df_indicator <- data.frame("var" = c(colVars(indicator_bpf), colVars(indicator_pf), colVars(indicator_mpf)))
df_indicator$algo <- rep(c("BPF", "PF", "MPF"), each = length(Nparticles))
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
# ggsave("lgssm_simple_indicator.pdf", width = 12, height = 8, dpi = 300)


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
# ggsave("lgssm_simple_legend.pdf", width = 10, height = 0.8, dpi = 300)

