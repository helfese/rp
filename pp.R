
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
