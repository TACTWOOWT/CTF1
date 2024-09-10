
library(limma)
library(estimate)
library(reshape2)
library(ggpubr)


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=data[rowMeans(data)>0,]
data=avereps(data)
#删除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]


out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#运行estimate包
filterCommonGenes(input.f="uniq.symbol.txt",                   
                                output.f="commonGenes.gct",                   
                                id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct",              
                      output.ds="estimateScore.gct")


sameSample=intersect(row.names(exp), row.names(score))
exp=exp[sameSample,"Type",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(score, exp)
rt$Type=factor(rt$Type, levels=c("Low", "High"))
#将合并后的数据转换为ggplot2的输入文件
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "scoreType", "Score")