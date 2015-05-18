#Please note that data folder should contain only data files
#Datasets should be in tab separated matrix format, each row containing expression of that and each column indicating that sample 
#There should be the only two unique column names in a dataset (each referring to one condition e.g. cancer, healthy) and these should be the same across all datasets 
#The first column should contain gene names
#Enter the two locations below:

data.loc = "your data folder path along with the data folder name"
res.loc = "your result folder path along with the result folder name"

#You do not need to enter anything after this line of code, run the rest of the code as it is!

extr.nth.elm=function(x, n)
return(unlist(x)[n])

data.files = dir(data.loc)

genes = NULL
conditions = NULL
for(k in 1:length(data.files))
{
	data = read.table(paste(data.loc, data.files[k], sep="/"), sep="\t", header=T, row.names = 1)
	data = as.matrix(data)
	genes = unique(c(genes, rownames(data)))
	conditions  = c(conditions, unique(colnames(data)))
}

conditions = sapply(strsplit(conditions, "[.]"), extr.nth.elm, 1)
conditions = unique(conditions)
if(length(conditions)!=2)
stop("one or more dataset(s) contain do not contain exactly two conditions")
genes = sort(unique(genes))

m2 = matrix(rep(genes, length(genes)), length(genes), length(genes))
m1 = t(m2)
ltm = lower.tri(m1)
n1 = m1[ltm]
n2 = m2[ltm]
pair.names = paste(n1, n2, sep="_")
rm(m1, m2, n1, n2)

write.table(genes, paste(res.loc, "gene.list.txt", sep="/"), row.names=F, sep="\t")
write.table(pair.names, paste(res.loc, "gene.pairs.list.txt", sep="/"), row.names=F, sep="\t")

cnt = matrix(rep(0, length(genes)*length(genes)), length(genes), length(genes) )
for(k in 1:length(data.files))
{
	data = read.table(paste(data.loc, data.files[k], sep="/"), sep="\t", header=T)
	data=as.matrix(data)
	data.genes = data[,1]
	data = data[,-1]
	hdrs = sapply(strsplit(colnames(data), "[.]"), extr.nth.elm, 1)
	colnames(data) = hdrs

	new.data = matrix(NA, nrow = length(genes), ncol = length(data[1,]))
	colnames(new.data) = colnames(data)
	index = match(genes, data.genes)
	new.data = data[index, ]
	rownames(new.data) = genes

	con1.data = new.data[, which(colnames(new.data) == conditions[1])]
	con2.data = new.data[, which(colnames(new.data) == conditions[2])]

	con1.data = t(apply(con1.data, 2, as.numeric))
	con2.data = t(apply(con2.data, 2, as.numeric))

	con1.mat = rcorr(con1.data)
	con2.mat = rcorr(con2.data)

	rm(data, new.data, con1.data, con2.data)

	con1.cor = con1.mat[[1]]
	con1.n = con1.mat[[2]]
	rm(con1.mat)
	
	ind.na.con1 = which(con1.n<4)

	con2.cor = con2.mat[[1]]
	con2.n = con2.mat[[2]]
	rm(con2.mat)

	ind.na.con2 = which(con2.n<4)
	
	con1.z = atanh(con1.cor)
	con2.z = atanh(con2.cor)
	diff.z = con1.cor - con2.cor
	var.z = (1/(con1.n - 3)) + (1/(con2.n - 3))
	cnt = cnt + 1
	
	ind.na = unique(c(which(is.na(diff.z)), ind.na.con1, ind.na.con2))
	if(length(ind.na) > 1){
	diff.z[ind.na] = 0
	var.z[ind.na] = 0
	cnt[ind.na] = cnt[ind.na] - 1 }
	
	if(k == 1){
	s.diff.z = diff.z
	s.var.z = var.z}
	
	rm(diff.z, var.z)
	
	if(k>1)
	s.diff.z = s.diff.z + diff.z
	s.var.z = s.var.z + var.z
	
}
q.diff.z = s.diff.z / sqrt(s.var.z)
rm(s.diff.z, s.var.z)
p.diff.z = 2*pnorm(abs(q.diff.z), lower.tail = F)

colnames(q.diff.z) = genes
rownames(q.diff.z) = genes
colnames(p.diff.z) = genes
rownames(p.diff.z) = genes
colnames(cnt) = genes
rownames(cnt) = genes

write.table(q.diff.z, paste(res.loc, "q_scores.txt", sep = "/"), sep = "\t")
write.table(p.diff.z, paste(res.loc, "p_values.txt", sep = "/"), sep = "\t")
write.table(cnt, paste(res.loc, "counts.txt", sep = "/"), sep = "\t")
