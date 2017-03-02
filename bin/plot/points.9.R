x = read.table("./plot/barplot")
title = read.table("./plot/title")

minx = min(x[,c(4,7,10,13,16,19,22,25)]);
maxx = max(x[,c(4,7,10,13,16,19,22,25)]);
miny = min(x[,c(2,5,8,11,14,17,20,23)]);
maxy = max(x[,c(2,5,8,11,14,17,20,23)]);

pdf("fig.pdf");
for (k in seq(1, length(x[,1])))
{
	xx = x[k,];
	tt = xx[1, 1];
	yy = xx[1, seq(2, length(xx[1,]))];
	mm = t(matrix(yy, nrow = 3));
	plot(-1, -1, xlim = c(minx, maxx), ylim = c(miny, maxy), main = tt, xlab = "Precision", ylab = "Number of Correct Transcripts");
	points(mm[c(1,2,3), c(3,1)], col = 2, pch = c(15,16,17));
	points(mm[c(4,5,6), c(3,1)], col = 3, pch = c(15,16,17));
	points(mm[c(7,8), c(3,1)], col = 4, pch = c(15,16));
}
dev.off();
