fancy.freq.plots3 <- function(index, FREQ, my_col = "gold") {
  #this function generates frequency plots of X+Y experiments
  #first input is the index; second parameters is a vector of frequencies;
  #third parameter is a number of genes
  #fourth parameter is an entropy
  
  tot.number = sum(FREQ)
  FREQ = FREQ / sum(FREQ)
  entropy = shannon.entropy(FREQ)
  
  D = 1.5
  M = 5
  d = 4.4
  EX = rbind(
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
  
  
  
  P = EX[index, ]
  
  ADD = P[1] + P[2] - P[1] + P[3] - P[1]
  
  UPB = max(c(P, ADD)) + 0.5
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  barplot(
    P,
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
  
  segments(0, ADD, 3.8, ADD, col = "red", lwd = 3.0)
  
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
  
  BREAKS = sort(unique(c(P, ADD, UPB)))
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
  #  if (sum(is.na(p)) == 0)
  #    return(NA)
  p.norm <- p[p > 0]
  - sum(log2(p.norm) * p.norm) / log2(length(p))
}
