auc_vector = function(v)
{
	m = length(v) / 2;
	x1 = 0;
	y1 = v[1];
	auc = 0;
	for(k in 1:m)
	{
		x2 = v[k * 2 - 0];
		y2 = v[k * 2 - 1];
		auc = auc + (y1 + y2) * 0.5 * (x2 - x1);
		x1 = x2;
		y1 = y2;
	}
	return(auc);
}

auc_file = function(file)
{
	data = as.matrix(read.table(file));
	m = length(data[,1]);

	auc = c();
	for(k in 1:m)
	{
		auc[k] = auc_vector(data[k,]);
	}
	return(auc);
}

m = 8;
n = 10;
auc = matrix(nrow = m, ncol = n);
auc[,1] = auc_file("plot/GSM981244");
auc[,2] = auc_file("plot/GSM981256");
auc[,3] = auc_file("plot/GSM984609");
auc[,4] = auc_file("plot/SRR307903");
auc[,5] = auc_file("plot/SRR307911");
auc[,6] = auc_file("plot/SRR315323");
auc[,7] = auc_file("plot/SRR315334");
auc[,8] = auc_file("plot/SRR387661");
auc[,9] = auc_file("plot/SRR534307");
auc[,10] = auc_file("plot/SRR545723");

pdf("plot/AUC10.pdf");

cols = c(2,2,2,3,3,3,4,4);
pchs = c(15,16,17,15,16,17,15,16);
p = seq(1, n / 2) * 2 - 0;
q = seq(1, n / 2) * 2 - 1;

xmin = min(x[,p]);
xmax = max(x[,p]);
ymin = min(x[,q]);
ymax = max(x[,q]);

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

dev.off();
