pdf("plot/ROC.pdf", width = 8, height = 6.5);
title = read.table("plot/title");

data1 = read.table("plot/roc1");
data2 = read.table("plot/roc2");
data3 = read.table("plot/roc3");
data4 = read.table("plot/roc4");
data5 = read.table("plot/roc5");
data6 = read.table("plot/roc6");
data7 = read.table("plot/roc7");
data8 = read.table("plot/roc8");


dataA = read.table("plot/rocA");
dataB = read.table("plot/rocB");
dataC = read.table("plot/rocC");
dataD = read.table("plot/rocD");
dataE = read.table("plot/rocE");
dataF = read.table("plot/rocF");
dataG = read.table("plot/rocG");
dataH = read.table("plot/rocH");


xmax = max(max(data1[,16]), max(data2[,16]), max(data3[,16]), max(data4[,16]), max(data5[,16]), max(data6[,16]), max(data7[,16]), max(data8[,16]) );
ymax = max(max(data1[,13]), max(data2[,13]), max(data3[,13]), max(data4[,13]), max(data5[,13]), max(data6[,13]), max(data7[,13]), max(data8[,13]) );
xmin = min(min(data1[,16]), min(data2[,16]), min(data3[,16]), min(data4[,16]), min(data5[,16]), min(data6[,16]), min(data7[,16]), min(data8[,16]) );
ymin = min(min(data1[,13]), min(data2[,13]), min(data3[,13]), min(data4[,13]), min(data5[,13]), min(data6[,13]), min(data7[,13]), min(data8[,13]) );

#plot(-1, -1, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Precision (%)", ylab = "Sensitivity (%)", main = title[1, 1]);
plot(-1, -1, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Precision (%)", ylab = "Sensitivity (%)");

grid();

lines(data1[,16], data1[,13], lwd = 2, col = 2, lty = 3);
lines(data2[,16], data2[,13], lwd = 2, col = 3, lty = 3);
lines(data3[,16], data3[,13], lwd = 2, col = 4, lty = 3);
lines(data4[,16], data4[,13], lwd = 2, col = 2, lty = 5);
lines(data5[,16], data5[,13], lwd = 2, col = 3, lty = 5);
lines(data6[,16], data6[,13], lwd = 2, col = 4, lty = 5);
lines(data7[,16], data7[,13], lwd = 2, col = 2, lty = 1);
lines(data8[,16], data8[,13], lwd = 2, col = 3, lty = 1);


points(dataA[,16], dataA[,13], lwd = 2, col = 2, pch = 2);
points(dataB[,16], dataB[,13], lwd = 2, col = 3, pch = 2);
points(dataC[,16], dataC[,13], lwd = 2, col = 4, pch = 2);
points(dataD[,16], dataD[,13], lwd = 2, col = 2, pch = 1);
points(dataE[,16], dataE[,13], lwd = 2, col = 3, pch = 1);
points(dataF[,16], dataF[,13], lwd = 2, col = 4, pch = 1);
points(dataG[,16], dataG[,13], lwd = 2, col = 2, pch = 19);
points(dataH[,16], dataH[,13], lwd = 2, col = 3, pch = 19);

legend(xmin + (xmax - xmin) * 0.7, ymax, c("HISAT", "STAR", "TopHat"), lwd = c(2, 2, 2), pch = c(19, 1, 2), lty = c(1, 5, 3));
legend(xmin, ymin + (ymax - ymin) * 0.2, c("Scallop", "StringTie", "TransComb"), lwd = c(2, 2, 2), col = c(2,3,4), lty = c(1, 1, 1));

dev.off();
