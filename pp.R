
###
library(ggplot2)
library(dplyr)
library(readr)
data <- read_csv("/data/Paises_PIB_ICH.csv")
# Filter data for countries in Europe and Americas
filtered_data <- data %>%
 filter(Continent %in% c("Europe", "Americas"))
ggplot(filtered_data, aes(x = GDP, y = HCI, color = Continent)) +
 geom_point() +
 scale_x_log10() +
 geom_text(data = filtered_data %>% filter(Country %in% c("Lithuania",
"Iceland", "United States", "Saint Lucia")),
 aes(label = Country), vjust = -1, hjust = 1) +
 labs(title = "Human Capital Index vs GDP per capita (2023)",
 x = "GDP per capita (log scale)",
 y = "Human Capital Index") +
 theme_minimal()

###
library(ggplot2)
library(dplyr)
library(readr)
data <- read_csv("/data/master.csv")
# Filter data for the year 2002 and age group 55-74 years
filtered_data <- data %>%
 filter(year == 2002, age == '55-74 years')
ggplot(filtered_data, aes(x = sex, y = `suicides/100k pop`, fill = sex)) +
 geom_boxplot() +
 labs(title = "Suicides per 100k population in 2002 for age group 55-74
years",
 x = "Sex",
 y = "Suicides per 100k population") +
 theme_minimal()

###
set.seed(seed)

ligado <- 0
aviso <- 0

probs <- 1:sinais * 2 / (sinais * (sinais + 1))

for (i in 1:simulacoes) {
  replica <- sample(sinais, size = circuitos, replace = TRUE, prob = probs)
  minimo <- min(replica)
  if (minimo > 1) {
    ligado <- ligado + 1
    if (minimo == 2) aviso <- aviso + 1
  }
}
solution <- aviso / ligado

###
set.seed(seed)

falsa_aceitacao <- 0
falsa_rejeicao <- 0

for(i in 1:m) {
  falsa_rejeicao <- falsa_rejeicao + (mean(rpois(n, lambda0)) > k)
  falsa_aceitacao <- falsa_aceitacao + (mean(rpois(n, lambda1)) <= k)
}

solution <- falsa_aceitacao / falsa_rejeicao

###
delta <- (b - a) / k
lim.classes <- a + delta * 0:k
classes <- cut(amostra, breaks = lim.classes)
freq.obs <- table(classes)

Fd <- function(x) {
  res <- ifelse(x < (a + b) / 2,
                2 * ((x -a) / (b - a))^2,
                1 - 2 * ((b - x) / (b - a))^2)
  return(res)
}
prob.H0 <- diff(Fd(lim.classes))

solution <- chisq.test(freq.obs, p = prob.H0)$p.value

#####################
#####################
subcjt <- subset(dados, Continent %in% continentes)
etiquetados <- subset(subcjt, Country %in% paises)

ggplot(subcjt, aes(x = GDP, y = HCI)) +
  geom_point(aes(color = Continent), size = 2, alpha = 0.8) +
  geom_label(data = etiquetados, aes(label = Country), size = 2.5, alpha = 0) +
  scale_x_log10() +
  labs(title = "Human Capital Index versus GDP per capita",
       x = "GDP per capita adjusted (international dollars, IMF, 2023)",
       y = "Human Capital Index (World Bank, 2020)")

library(ggplot2)
theme_set(theme_minimal())


#####################
subcjt <- subset(dados, year == ano & age == idade)

ggplot(subcjt) +
  geom_boxplot(aes(sex, `suicides/100k pop`)) + 
  labs(title = paste("Rate of suicides for age group", idade, "in", ano))

library(ggplot2)
theme_set(theme_light())


#####################
dados <- subset(eletro, COUNTRY %in% paises & YEAR >= 2015 & PRODUCT == "Renewables")
dados$DATE <- as.Date(paste(1, dados$MONTH, dados$YEAR), format = "%d %m %Y") 
dados$share <- as.numeric(dados$share) * 100

ggplot(dados) +
  geom_line(aes(DATE, share, color = COUNTRY), linewidth = 1) +
  ylim(c(0, 100)) +
  labs(title="Monthly share of renewable sources in electricity production",
       y = "RENEWABLE SOURCES (% of total production)")


#####################
set.seed(seed)

ligado <- 0
aviso <- 0

probs <- 1:sinais * 2 / (sinais * (sinais + 1))

for (i in 1:simulacoes) {
  replica <- sample(sinais, size = circuitos, replace = TRUE, prob = probs)
  minimo <- min(replica)
  if (minimo > 1) {
    ligado <- ligado + 1
    if (minimo == 2) aviso <- aviso + 1
  }
}
solution <- aviso / ligado


#####################
set.seed(seed)
prop <- vector(length = r)

for(j in 1:r) {
  prop[j] <- 0
  for(i in 1:m) {
    z <- rnorm(n + 1)
    t <- z[1] / sqrt(sum(z[-1]^2) / n)
    prop[j] <- prop[j] + (t <= valor)
  }
}
prop <- prop / m

solution <- abs(pt(valor, n) - mean(prop)) * 100


#####################
set.seed(seed)
y <- vector(length = 1000)

for (i in 1:1000) {
  y[i] <- sum(rexp(n, rate = 1 / a))
}
sim <- sum(y > valor) / 1000

ext <- 1 - pgamma(valor, n, rate = 1 / a)

solution <- abs(sim - ext) * 100


#####################
library(stats4)

soma_logdata <- sum(log(data))
n <- length(data)

fun <- function(theta) {
  res <- -n * log(theta) + (theta + 1) * soma_logdata - n * theta * log(a)
  return(res)
}

result <- mle(minuslog = fun, start = list(theta = valor_inicial))
est_theta <- as.vector(result@coef)

quantil <- function(p, theta) {
  res <- a * (1 - p)^(-1 / theta)
}

estimativa <- quantil(p, est_theta)
valor_exato <- quantil(p, theta)

solution <- abs(estimativa - valor_exato)

#####################
set.seed(seed)
dados <- sample(amostra, n)

a <- qchisq((1 - gama) / 2, n - 1)
b <- qchisq((1 + gama) / 2, n - 1)
amp_ic <- (n - 1) * var(dados) * (1 / a - 1 / b)

library(pracma)

restricoes <- function(x) {
  F1 <- pchisq(x[2], n - 1) - pchisq(x[1], n - 1) - gama
  F2 <- dchisq(x[2], n + 3) - dchisq(x[1], n + 3)
  return(c(F1, F2))
}

raiz <- fsolve(restricoes, c(a, b))$x

amp_ic_min <- (n - 1) * var(dados) * (1 / raiz[1] - 1 / raiz[2])

solution <- amp_ic - amp_ic_min

#####################
set.seed(seed)

falsa_aceitacao <- 0
falsa_rejeicao <- 0

for(i in 1:m) {
  falsa_rejeicao <- falsa_rejeicao + (mean(rpois(n, lambda0)) > k)
  falsa_aceitacao <- falsa_aceitacao + (mean(rpois(n, lambda1)) <= k)
}

solution <- falsa_aceitacao / falsa_rejeicao

#####################
delta <- (b - a) / k
lim.classes <- a + delta * 0:k
classes <- cut(amostra, breaks = lim.classes)
freq.obs <- table(classes)

Fd <- function(x) {
  res <- ifelse(x < (a + b) / 2,
                2 * ((x -a) / (b - a))^2,
                1 - 2 * ((b - x) / (b - a))^2)
  return(res)
}
prob.H0 <- diff(Fd(lim.classes))

solution <- chisq.test(freq.obs, p = prob.H0)$p.value
