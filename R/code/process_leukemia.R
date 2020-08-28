#加载文件
#file_path = "/Users/liguangping/Desktop/R/data/GSE9476_RAW";
#file_path = "/Users/brown/Desktop/R/data/GSE9476_RAW";
file_path = "/Users/brown/paper/R/data/GSE9476_RAW";
out_path = "/Users/brown/paper/R/output_files/"
setwd(file_path);
library(simpleaffy);
raw_data = read.affy(covdesc = 'covdesc.txt', path = '.')

#质量控制
Data.qc <- qc(raw_data);
plot(Data.qc)

my_filter = Data.qc@scale.factors <= 3
#t = Data.qc@average.background <= 50
#my_filter = my_filter & t
Ra = ratios(Data.qc)
t = Ra[,'gapdh3/gapdh5']<=1.25
my_filter = my_filter & t
t = Ra[,'actin3/actin5']<=3
my_filter = my_filter & t
raw_data = raw_data[,my_filter]
Data.qc <- qc(raw_data);
plot(Data.qc)

library(affyPLM);
Pset <- fitPLM(raw_data);
RLE(Pset);
NUSE(Pset);
Pset_stat = NUSE(Pset,type="stats")
my_filter = Pset_stat[1,] < 1.05
raw_data = raw_data[,my_filter]

Pset <- fitPLM(raw_data);
NUSE(Pset);

library(affy);
Data.deg <- AffyRNAdeg(raw_data);
plotAffyRNAdeg(Data.deg)

#归一化并分组
eset = mas5(raw_data)
ill = get.array.subset(eset, 'state', 'leukemia')
normal = get.array.subset(eset, 'state', 'normal')

bone_marrow_ill = get.array.subset(ill, 'type', 'bone_marrow')
bone_marrow_normal = get.array.subset(normal, 'type', 'bone_marrow')

peripheral_blood_ill = get.array.subset(ill, 'type', 'peripheral_blood')
peripheral_blood_normal = get.array.subset(normal, 'type', 'peripheral_blood')

#保存

##保存bone_marrow_ill
ncol(bone_marrow_ill) 
file_name = paste(out_path, 'leukemia_bone_marrow_ill_7mas5.csv',sep = "")
mat = exprs(bone_marrow_ill)
write.csv(mat, file_name)

##保存bone_marrow_normal
ncol(bone_marrow_normal) 
file_name = paste(out_path, 'leukemia_bone_marrow_normal_10mas5.csv',sep = "")
mat = exprs(bone_marrow_normal)
write.csv(mat, file_name)
