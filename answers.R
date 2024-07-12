# 1 #
library(ggplot2)
theme_set(theme_classic())
subcjt <- subset(data, Continent %in% continents)
labeled <- subset(subcjt, Country %in% countries)
ggplot(subcjt, aes(x = GDP, y = HCI)) +
  geom_point(aes(color = Continent), size = 2, alpha = 0.8) +
  geom_label(data = labeled, aes(label = Country), size = 2.5, alpha = 0) +
  scale_x_log10() +
  labs(title = "Human Capital Index versus GDP per capita",
       x = "GDP per capita adjusted (international dollars, IMF, 2023)",
       y = "Human Capital Index (World Bank, 2020)")

# 2 #
library(ggplot2)
library(dplyr)
library(readr)
data <- read_csv("/data/master.csv")
filtered_data <- data %>%
 filter(year == 2002, age == '55-74 years')
ggplot(filtered_data, aes(x = sex, y = `suicides/100k pop`, fill = sex)) +
 geom_boxplot() +
 labs(title = "Suicides per 100k population in 2002 for age group 55-74 years",
 x = "Sex",
 y = "Suicides per 100k population") +
 theme_minimal()

# 3 #
library(ggplot2)
theme_set(theme_light())
data <- subset(electricity, COUNTRY %in% countries & YEAR >= 2015 & PRODUCT == "Renewables")
data$DATE <- as.Date(paste(1, data$MONTH, data$YEAR), format = "%d %m %Y") 
data$share <- as.numeric(data$share) * 100

ggplot(data) +
  geom_line(aes(DATE, share, color = COUNTRY), linewidth = 1) +
  ylim(c(0, 100)) +
  labs(title="Monthly share of renewable sources in electricity production",
       y = "RENEWABLE SOURCES (% of total production)")

# 4 #
set.seed(seed)
probabilities <- signals / sum(signals)
for (realization in 1:num_realizations) {
  emitted_signals <- sample(signals, num_circuits, replace = TRUE, prob = probabilities)
  if (2 %in% emitted_signals && !(1 %in% emitted_signals)) {
    warning_and_not_off_count <- warning_and_not_off_count + 1
  }
  if (!(1 %in% emitted_signals)) {
    not_off_count <- not_off_count + 1
  }
}
proportion <- warning_and_not_off_count / not_off_count

# 5 #
set.seed(seed)
prop <- vector(length = r)
for(j in 1:r) {
  prop[j] <- 0
  for(i in 1:m) {
    z <- rnorm(n + 1)
    t <- z[1] / sqrt(sum(z[-1]^2) / n)
    prop[j] <- prop[j] + (t <= value)
  }
}
prop <- prop / m
solution <- abs(pt(value, n) - mean(prop)) * 100

# 6 #
set.seed(seed)
samples <- replicate(num_samples, sum(rexp(n, rate = lambda)))
proportion_greater_90 <- mean(samples > 90)
exact_prob <- 1 - pgamma(90, shape = shape, scale = scale)
difference <- abs(proportion_greater_90 - exact_prob) * 100
rounded_difference <- round(difference, 4)
rounded_difference

# 7 #
library(stats4)
logLik <- function(theta) {
  n <- length(observations)
  logL <- n * log(theta) + n * theta * log(a) - (theta + 1) * sum(log(observations))
  return(-logL)
}
mle_est <- mle(logLik, start = list(theta = 3.6))
theta_hat <- coef(mle_est)
quantile_75 <- a * (1 / (1 - p))^(1 / theta_hat)
quantile_75_true <- a * (1 / (1 - p))^(1 / theta_true)
absolute_deviation <- abs(quantile_75 - quantile_75_true)
rounded_absolute_deviation <- round(absolute_deviation, 4)
rounded_absolute_deviation

# 8 #
set.seed(seed)
data <- sample(sample_data, n)
a <- qchisq((1 - gamma) / 2, n - 1)
b <- qchisq((1 + gamma) / 2, n - 1)
ci_width <- (n - 1) * var(data) * (1 / a - 1 / b)
library(pracma)
constraints <- function(x) {
  F1 <- pchisq(x[2], n - 1) - pchisq(x[1], n - 1) - gamma
  F2 <- dchisq(x[2], n + 3) - dchisq(x[1], n + 3)
  return(c(F1, F2))
}
root <- fsolve(constraints, c(a, b))$x
min_ci_width <- (n - 1) * var(data) * (1 / root[1] - 1 / root[2])
solution <- ci_width - min_ci_width

# 9 #
set.seed(seed)
for (i in 1:m) {
  sample_H0 <- rpois(n, lambda0)
  sample_H1 <- rpois(n, lambda1)
  mean_H0 <- mean(sample_H0)
  mean_H1 <- mean(sample_H1)
  if (mean_H0 > k) {
    type_I_error <- type_I_error + 1
  }
  if (mean_H1 <= k) {
    type_II_error <- type_II_error + 1
  }
}
type_I_error_freq <- type_I_error / m
type_II_error_freq <- type_II_error / m
error_ratio <- type_II_error_freq / type_I_error_freq

# 10 #
F_X <- function(x) {
  if (x < a) {
    return(0)
  } else if (x < (a + b) / 2) {
    return(2 * ((x - a) / (b - a))^2)
  } else if (x < b) {
    return(1 - 2 * ((b - x) / (b - a))^2)
  } else {
    return(1)
  }
}
breaks <- seq(a, b, length.out = k + 1)
observed_factor <- cut(observed, breaks = breaks, include.lowest = TRUE, right = FALSE)
observed_freq <- table(observed_factor)
expected_freq <- numeric(k)
for (i in 1:k) {
  expected_freq[i] <- n * (F_X(breaks[i + 1]) - F_X(breaks[i]))
}
chisq_test <- chisq.test(observed_freq, p = expected_freq / sum(expected_freq))
