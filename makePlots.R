## Make plots for the presentation since ggplot2 is not installed in ShareLatex
library("ggplot2")
load("plots.Rdata")

pdf("eda1.pdf")
eda1
dev.off()

pdf("eda2.pdf")
eda2
dev.off()

pdf("p1.pdf")
p1
dev.off()

pdf("p2.pdf")
p2
dev.off()

pdf("p3.pdf")
p3
dev.off()

pdf("p4.pdf")
p4
dev.off()

pdf("p5.pdf")
p5
dev.off()

pdf("p6.pdf")
p6
dev.off()

