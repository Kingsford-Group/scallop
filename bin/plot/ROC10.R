draw10roc = function(file)
{
	x = read.table(file);

	m = length(x[,1]);
	n = length(x[1,]);
	cols = c(2,2,2,3,3,3,4,4);
	pchs = c(15,16,17,15,16,17,15,16);
	p = seq(1, n / 2) * 2 - 0;
	q = seq(1, n / 2) * 2 - 1;

	xmin = min(x[,p]);
	xmax = max(x[,p]);
	ymin = min(x[,q]);
	ymax = max(x[,q]);
	#xmin = 0;
	#xmax = 100;
	#ymin = 0;
	#ymax = 12;

	plot(-1, -1, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Precision", ylab = "Sensitivity", main = file);
	#plot(-1, -1, xlim = c(xmin, xmax), ylim = c(ymin, ymax));

	for (k in seq(1, m))
	{
		a = x[k, p];
		b = x[k, q];
		points(t(a), t(b), col = cols[k], pch = pchs[k]);
		lines(t(a), t(b), col = cols[k]);
	}

	legend(xmax - 15, ymax, c("Scallop", "StringTie", "TransComb"), col = c(2,3,4), lty = c(1,1,1));
	legend(xmax - 15, ymax - 1.5, c("TopHat", "STAR", "HISAT"), pch = c(15,16,17));
}

pdf("plot/ROC10.pdf");
#par(mar=c(5,5,2,2) + 0.1);
#par(mfrow=c(2,5)) 

draw10roc("plot/GSM981244");
draw10roc("plot/GSM981256");
draw10roc("plot/GSM984609");
draw10roc("plot/SRR307903");
draw10roc("plot/SRR307911");
draw10roc("plot/SRR315323");
draw10roc("plot/SRR315334");
draw10roc("plot/SRR387661");
draw10roc("plot/SRR534307");
draw10roc("plot/SRR545723");

dev.off();
