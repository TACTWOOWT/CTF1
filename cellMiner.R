#清空当前工作环境
rm(list = ls())

library(readxl)
rt1<- read_excel(path = "DTP_NCI60_ZSCORE.xlsx", skip = 7)
colnames(rt1)<- rt1[1,]
rt1<-rt1[-1,-c(67,68)]

#筛选药物标准
table(rt1$`FDA status`)
rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
rt1 <- rt1[,-c(1, 3:6)]
write.table(rt1, file = "drug.txt",sep = "\t",row.names = F,quote = F)

#基因表达数据的准备
rt2<- read_excel(path = "RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(rt2)<-rt2[1,]
rt2<-rt2[-1,-c(2:6)]
write.table(rt2,file = "geneExp.txt",sep = "\t",row.names = F,quote = F)

#药物敏感性分析
rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")

library(impute)
library(limma)

#读取药物输入文件
rt <- read.table("drug.txt",sep="\t",header=T,check.names=F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
drug <- rt[,2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)

mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug)

#读取表达输入文件
exp <- read.table("geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]

#提取特定基因表达
gene <- read.table("gene.txt",sep="\t",header=F,check.names=F)
genelist <- as.vector(gene[,1])
genelist


genelist <- gsub(" ","",genelist)
genelist <- intersect(genelist,row.names(exp))
exp <- exp[genelist,]

#药物敏感性计算
outTab <- data.frame()

for(Gene in row.names(exp)){
    x <- as.numeric(exp[Gene,])
    #对药物循环
    for(Drug in row.names(drug)){
        y <- as.numeric(drug[Drug,])
        corT <- cor.test(x,y,method="pearson")
        cor <- corT$estimate
        pvalue <- corT$p.value
        if(pvalue < 0.01){
            outVector <- cbind(Gene,Drug,cor,pvalue)
            outTab <- rbind(outTab,outVector)
        }
    }
}


outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab, file="drugCor.txt", sep="\t", row.names=F, quote=F)


library(ggplot2)
library(ggpubr)

#散点图
plotList_1 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
    corPlotNum=nrow(outTab)
}

for(i in 1:corPlotNum){
    Gene <- outTab[i,1]
    Drug <- outTab[i,2]
    x <- as.numeric(exp[Gene,])
    y <- as.numeric(drug[Drug,])
    cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
    pvalue=0
    if(as.numeric(outTab[i,4])<0.001){
        pvalue="p<0.001"
    }else{
        pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
    }
    df1 <- as.data.frame(cbind(x,y))
    p1=ggplot(data = df1, aes(x = x, y = y))+
        geom_point(size=1)+
        stat_smooth(method="lm",se=FALSE, formula=y~x)+
        labs(x="Expression",y="IC50",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
        theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
        theme_bw()
    plotList_1[[i]]=p1
}

#箱线图
plotList_2 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
    corPlotNum=nrow(outTab)
}

for(i in 1:corPlotNum){
    Gene <- outTab[i,1]
    Drug <- outTab[i,2]
    x <- as.numeric(exp[Gene,])
    y <- as.numeric(drug[Drug,])
    df1 <- as.data.frame(cbind(x,y))
    colnames(df1)[2] <- "IC50"
    df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
    compaired <- list(c("low", "high"))
    p1 <- ggboxplot(df1,
                    x = "group", y = "IC50",
                    fill = "group", palette = c("#00AFBB", "#E7B800"),
                    add = "jitter", size = 0.5,
                    xlab = paste0("The_expression_of_", Gene),
                    ylab = paste0("IC50_of_", Drug)) +
        stat_compare_means(comparisons = compaired,
                           method = "wilcox.test",   #设置统计方法
                           symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                            symbols = c("***", "**", "*", "ns")))
    plotList_2[[i]]=p1
}


nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)
ggarrange(plotlist=plotList_2,nrow=nrow,ncol=ncol)