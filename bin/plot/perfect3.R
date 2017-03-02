x = read.table("./plot/barplot")
#title = read.table("./plot/title")

pdf("./plot/perfect.pdf",width=8,height=7)
q = c(4,7,10)
cc = c(2,3,4)
r=c(2,5,8)

maxy=max(x[,q]) + 10

#barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Sensitivity (%)", main = title[1,1])
barx = barplot(t(x[,q]), beside=TRUE, col=cc, names.arg=t(x[,1]), ylim=c(0,maxy), ylab = "Precision (%)", grid = TRUE)
text(x = barx, y = t(x[,q]), label = t(x[,r]), pos = 3, cex = 1.0, col = c(2,3,4))
axis(2, seq(0,maxy,by=10))

legend(3.0, maxy, c("Scallop", "Scallop + StringTie", "Scallop + StringTie + TransComb"), fill=c(2,3,4));
dev.off()
