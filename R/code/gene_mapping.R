file_path = "/Users/brown/paper/R/data/GSE27567_RAW";
setwd(file_path)
raw_data = ReadAffy();

library(affy);
eset = rma(raw_data)
eset_mas5 = mas5(raw_data)
mat = exprs(eset_mas5)
write.csv(mat, 'eset_matrix_23mas5.csv')

library(annotate)
affy_IDs = featureNames(eset_mas5)
an = annotation(eset_mas5)
affydb = annPkgName(an, type = "db")
library(affydb, character.only = TRUE)

library("hgu133plus2.db")
#EntrezID = select(hgu133plus2.db, keys = affy_IDs, columns = 'ENTREZID')
UniprotID = select(hgu133plus2.db, keys = affy_IDs, columns = 'UNIPROT')
#write.csv(EntrezID, 'EntrezID.csv')
write.csv(UniprotID, 'probe2UniprotAC_file.csv')#目标文件
