x = read.table("./plot/barplot")
title = read.table("./plot/title")

pdf("./plot/sen.pdf")
q=c(3,6);
r=c(2,5)

maxy=max(x[,q]) + 3

barx = barplot(t(x[,q]), beside=TRUE, col=c(2,3), names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Sensitivity (%)", main = title[1,1])
text(x = barx, y = t(x[,q]), label = t(x[,r]), pos = 3, cex = 1.0, col = c(2,3))
axis(2, seq(0,maxy,by=2))

legend(5.0, maxy, c("Scallop", "StringTie"), fill=c(2,3));

dev.off()

pdf("./plot/pre.pdf")
q=c(4,7);
r=c(2,5)

maxy=max(x[,q]) + 10

barx = barplot(t(x[,q]), beside=TRUE, col=c(2,3), names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Precision (%)", main = title[1,1])
#text(x = barx, y = t(x[,q]), label = t(x[,r]), pos = 3, cex = 1.0, col = c(2,3))
axis(2, seq(0,maxy,by=10))

legend(5.0, maxy, c("Scallop", "StringTie"), fill=c(2,3));

dev.off()
