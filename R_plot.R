args <- commandArgs(trailingOnly = T)

pwd <- toString(args[4])
mask <- toString(args[5])
name <- toString(args[6])

"%+%" <- function(...){
  paste0(...)
}
name_1 <- pwd %+% "/LD_" %+% name %+% ".txt"
LD_txt <- read.csv(name_1, sep = '\t', 
                 stringsAsFactors = F)
LD_txt <- LD_txt[order(LD_txt$Len, decreasing = F),] 
x <- LD_txt$Len
y <- LD_txt$LD

wind <- as.integer(round(length(x)/10, 0))
if (wind == 0) {wind <- 2}

name_plot <- pwd %+% "/LD_plot_" %+% name %+% "_" %+% mask %+% "_" %+% wind  %+% ".png"

png(file= name_plot, width=1500, height=900, res=120)
plot(x, y, col=grey(.7),  main = "LD plot of " %+% name %+% ". " %+% "Mask = " %+% mask %+% ". Window = " %+% wind, xlab ="Distance", ylab = "LD", ylim=c(0, max(y)))
grid()
f <- rep(1/wind, wind)
y_lag <- filter(y, f, sides=1)
lines(x, y_lag, col="red")
dev.off()
