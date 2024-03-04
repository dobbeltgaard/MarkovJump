#################################
### Defect evolution forecast ###
#################################
library(Matrix)
library(ggplot2)
library(tidyr)

forecast.state.probs <- function(states, data, class, track, position, pars, exo.cols, span) {
  span <- 1:(span * 365) # time span
  m <- length(states) # number of states
  dp <- data[data$Track == track, ] # define working data
  
  idx.From.To <- dp$From <= position & dp$To >= position # position within defined track ranges
  idx.To <- which.min(abs(dp$To - position)) # closest point in track ranges
  if (sum(idx.From.To) > 0) { # check if position is within defined track ranges
    idx <- which(idx.From.To)[1] # find defect position
  } else {
    idx <- idx.To # find defect position
  }
  # print(dp[idx, ])
  exo <- exp(sum(pars[m:length(pars)] * dp[idx, exo.cols])) # compute exogenous contribution
  lambda <- exp(pars[1:(m - 1)]) * exo # transform base transition rates to natural domain
  
  A <- construct.A(m, lambda, A_param = "A3") # transition rate matrix
  
  Pt <- states == class # compute init state probs
  Ptu <- matrix(NaN, ncol = m, nrow = length(span)) # initialize state distribution
  for (i in 1:length(span)) {
    Ptu[i, ] <- as.matrix(Pt %*% Matrix::expm(A * span[i] / 365))
  }
  
  
  ##### MAKE GGPLOT #####
  dpp <- as.data.frame(Ptu)
  colnames(dpp) <- states
  dpp$time <- span
  text.size <- 12
  
  dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
  
  p <- ggplot(data = dpp_long, aes(x = time, y = Probability, color = Column)) +
    geom_line(size = 0.75) +
    theme(
      text = element_text(size = text.size, family = "serif"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_line(color = "lightgray"),
      legend.key = element_blank()
    ) +
    annotate("text",
             x = max(dpp_long$time), y = 1, size = text.size / 2, label = str_to_title(track),
             family = "serif", vjust = 1, hjust = 1
    ) +
    annotate("text",
             x = max(dpp_long$time), y = 0.925, size = text.size / 2.5, label = paste0("Class ", class, " at ", position, "km"),
             family = "serif", vjust = 1, hjust = 1
    ) +
    xlab("Time [days]") +
    ylab("Probability") +
    labs(color = "Classes")
  
  return(list(Ptu, p))
}


forecast.state.probs2 <- function(m, states, z, class, pars, span, position, track) {
  span <- 1:(span * 365) # time span
  m <- length(states) # number of states
  
  exo <- exp(sum(pars[m:length(pars)] * z)) # compute exogenous contribution
  lambda <- exp(pars[1:(m - 1)]) * exo # transform base transition rates to natural domain
  
  A <- construct.A(m, lambda, A_param = "A3") # transition rate matrix
  
  Pt <- states == class # compute init state probs
  Ptu <- matrix(NaN, ncol = m, nrow = length(span)) # initialize state distribution
  for (i in 1:length(span)) {
    Ptu[i, ] <- as.matrix(Pt %*% Matrix::expm(A * span[i] / 365))
  }
  
  
  ##### MAKE GGPLOT #####
  dpp <- as.data.frame(Ptu)
  colnames(dpp) <- states
  dpp$time <- span
  text.size <- 16
  
  dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
  
  p <- ggplot(data = dpp_long, aes(x = time, y = Probability, color = Column)) +
    geom_line(size = 0.75) +
    theme(
      text = element_text(size = text.size, family = "serif"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_line(color = "lightgray"),
      legend.position = c(0.6, 0.8),
      legend.direction = "horizontal", # Set legend direction to horizontal
      legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
      legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
    ) +
    # scale_y_continuous(breaks = 1:5) +
    annotate("text",
             x = max(dpp_long$time), y = 1, size = text.size / 2.5, label = str_to_title(get_full_name(track)),
             family = "serif", vjust = 1, hjust = 1
    ) +
    annotate("text",
             x = max(dpp_long$time), y = 0.925, size = text.size / 3, label = paste0(class, " propagation at ", position, "km"),
             family = "serif", vjust = 1, hjust = 1
    ) +
    xlab("Time [days]") +
    ylab("Probability") +
    labs(color = "Classes")
  
  return(list(Ptu, p))
}