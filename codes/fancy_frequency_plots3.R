fancy_freq_plots3 <- function(index, FREQ, my_col = "gold") {
  #this function generates frequency plots of X+Y experiments
  #INPUTS: index = index of profile; FREQ = vector of frequencies of each outcome;
  #my_col = color of frequency bars

  #normalize outcome frequencies
  tot.number = sum(FREQ)
  FREQ = FREQ / sum(FREQ)
  entropy = shannon.entropy(FREQ)
  
  #the lines below define a graphical representation of the 17 cases
  
  D = 1.5
  M = 5
  d = 4.4
  
  cases = rbind(
    c(M - d, M - d, M - d, 0),
    c(M, M - D, M - D, 0),
    c(M, M + D, M + D, 0),
    c(M, M + D, M, 0),
    c(M, M, M + D, 0),
    c(M, M - D, M, 0),
    c(M, M, M - D, 0),
    c(M, M + D, M + 2 * D, 0),
    c(M, M + 2 * D, M + D, 0),
    c(M, M - D, M - 2 * D, 0),
    c(M, M - 2 * D, M - D, 0),
    c(M, M - D, M + D, 0),
    c(M, M + D, M - D, 0),
    c(M, M + 2 * D, M - D, 0),
    c(M, M - D, M + 2 * D, 0),
    c(M, M - 2 * D, M + D, 0),
    c(M, M + D, M - 2 * D, 0)
  )
  
  
  #P = cases[index, ]
  Mu = cases[index, ]
  
  additive = Mu[1] + Mu[2] - Mu[1] + Mu[3] - Mu[1]
  
  UPB = max(c(Mu, additive)) + 0.5
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  barplot(
    Mu,
    adj = 0,
    names.arg = c("0", "X", "Y", "X+Y"),
    main = paste(
      "case=",
      as.character(index),
      "\nn tot=",
      as.character(tot.number),
      "\nentropy=",
      as.character(round(entropy, 2))
    ),
    cex.names = 1.5,
    ylim = c(0, UPB),
    yaxt = 'n',
    ann = FALSE,
    col = "gray40",
    font = 2
  )
  
  segments(0, additive, 3.8, additive, col = "red", lwd = 3.0)
  
  axis(
    3,
    xaxp = c(3.8, 4.8, 2),
    labels = c(0, 0.5, 1),
    at = seq(3.8, 4.8, by = 0.5),
    cex.axis =
      0.8,
    font = 2
  )
  
  w = .15
  
  
  if (index == 1)
    w = 0.05
  
  BREAKS = sort(unique(c(Mu, additive, UPB)))
  BREAKSl = BREAKS[-c(1, length(BREAKS))] - w
  BREAKSu = BREAKS[-c(1, length(BREAKS))] + w
  BREAKS = sort(c(BREAKS[1], BREAKSl, BREAKSu, BREAKS[length(BREAKS)]))
  
  for (h in (2:(length(BREAKS)))) {
    col = "gray90"
    if (h %% 2 != 0)
      col = "gray80"
    rect(
      3.8,
      BREAKS[h - 1],
      4.8,
      BREAKS[h],
      lwd = 1,
      border = "white",
      col = col
    )
  }
  
  #lightgreen
  #gold
  
  for (h in (2:(length(BREAKS)))) {
    rect(3.8,
         BREAKS[h - 1],
         3.802 + FREQ[h - 1],
         BREAKS[h],
         col = my_col,
         border = "gray40")
  }
  
}


shannon.entropy <- function(p)
{
  p.norm <- p[p > 0]
  - sum(log2(p.norm) * p.norm) / log2(length(p))
}
