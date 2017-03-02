x = read.table("./plot/barplot")
title = read.table("./plot/title")

pdf("./plot/sen.pdf",width=8,height=7)
q = c(3,6,9,12,15,18,21,24)
cc = c(2,2,2,3,3,3,4,4)
#r=c(2,5,8)

maxy=max(x[,q]) + 2

#barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Sensitivity (%)", main = title[1,1])
barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Sensitivity (%)", grid = TRUE)
axis(2, seq(0,maxy,by=2))

legend(5.0, maxy, c("Scallop(TopHat, STAR, HISAT)", "StringTie(TopHat, STAR, HISAT)", "TransComb(TopHat, STAR)"), fill=c(2,3,4));
dev.off()

pdf("./plot/pre.pdf",width=8,height=7)
q=c(4,7,10,13,16,19,22,25)

maxy=max(x[,q]) + 10

#barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Precision (%)", main = title[1,1])
barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Precision (%)", grid = TRUE)
axis(2, seq(0,maxy,by=10))

legend(5.0, maxy, c("Scallop(TopHat, STAR, HISAT)", "StringTie(TopHat, STAR, HISAT)", "TransComb(TopHat, STAR)"), fill=c(2,3,4));

dev.off()
