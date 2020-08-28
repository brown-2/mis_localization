#加载文件
#file_path = "/Users/brown/paper/R/data/GSE121248_RAW(hepatitis cavinma)/normal"
#file_path = "/Users/brown/paper/R/data/GSE121248_RAW(hepatitis cavinma)/ill"
file_path = "/Users/brown/paper/R/data/GSE121248_RAW(all)"
out_path = "/Users/brown/paper/R/output_files/"
setwd(file_path);
library(simpleaffy);
raw_data = read.affy(covdesc = 'covdesc.txt', path = '.')
#raw_data = ReadAffy()
#质量控制
Data.qc <- qc(raw_data);
plot(Data.qc)

my_filter = Data.qc@scale.factors <= 3
t = Data.qc@average.background <= 50
my_filter = my_filter & t
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
ill = get.array.subset(eset, 'state', 'Tumor')
normal = get.array.subset(eset, 'state', 'Normal')

#保存

##保存ill
ncol(ill)
file_name = paste(out_path, 'hepatitis_ill_21mas5.csv',sep = "")
mat = exprs(ill)
write.csv(mat, file_name)

##保存normal
ncol(normal) 
file_name = paste(out_path, 'hepatitis_normal_6mas5.csv',sep = "")
mat = exprs(normal)
write.csv(mat, file_name)
