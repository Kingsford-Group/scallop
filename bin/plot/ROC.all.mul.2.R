pdf("plot/ROC.pdf");

title = read.table("plot/title");

data1 = read.table("plot/roc1");
data2 = read.table("plot/roc2");
data4 = read.table("plot/roc4");
data5 = read.table("plot/roc5");

dataA = read.table("plot/rocA");
dataB = read.table("plot/rocB");
dataD = read.table("plot/rocD");
dataE = read.table("plot/rocE");

xmax = max(max(data1[,16]), max(data2[,16]), max(data4[,16]), max(data5[,16]));
ymax = max(max(data1[,13]), max(data2[,13]), max(data4[,13]), max(data5[,13]));
xmin = min(min(data1[,16]), min(data2[,16]), min(data4[,16]), min(data5[,16]));
ymin = min(min(data1[,13]), min(data2[,13]), min(data4[,13]), min(data5[,13]));

plot(-1, -1, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Precision (%)", ylab = "Sensitivity (%)", main = title[1, 1]);

grid();

lines(data1[,16], data1[,13], col = 2, lwd = 2, lty = 1);
lines(data2[,16], data2[,13], col = 3, lwd = 2, lty = 1);
lines(data4[,16], data4[,13], col = 2, lwd = 2, lty = 2);
lines(data5[,16], data5[,13], col = 3, lwd = 2, lty = 2);

points(dataA[,16], dataA[,13], col = 2, lwd = 2, pch = 19);
points(dataB[,16], dataB[,13], col = 3, lwd = 2, pch = 19);
points(dataD[,16], dataD[,13], col = 2, lwd = 2, pch = 1);
points(dataE[,16], dataE[,13], col = 3, lwd = 2, pch = 1);

legend(xmin + (xmax - xmin) * 0.7, ymax, c("All Transcripts", "Multiple Exon"), lwd = c(2, 2), pch = c(19, 1), lty = c(1, 2));
legend(xmin, ymin + (ymax - ymin) * 0.2, c("Scallop", "StringTie"), lwd = c(2, 2), col = c(2,3), lty = c(1, 1));


dev.off();
