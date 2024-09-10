load('TCGA_XXX.RData')
tcga.data=log(dt+1)
load('TCGA_XXX.RData')
gtex.data=log(dt+1)
load('Clinical.RData')
gene.type=read.csv('GeneTag.txt',stringsAsFactors = F,row.names = 1,sep = '\t',check.names = F)
gene.type1=gene.type[which(gene.type$TYPE=='protein_coding'),]
nor.data=tcga.data[,grep('-11',colnames(tcga.data))]
comsamples=intersect(c(paste0(rownames(clinical),'-01'),paste0(rownames(clinical),'-03'),paste0(rownames(clinical),'-06')),colnames(tcga.data))
tum.data=tcga.data[,comsamples]
clin.cut=clinical[substr(comsamples,1,12),]
comgenes=intersect(intersect(rownames(tcga.data),rownames(gtex.data)),gene.type1$SYMBOL)

ttt1= data_tum1 %>%
  rownames_to_column('Samples') %>%
  pivot_longer(cols = 2:(ncol(data_tum1)+1),names_to='Celltype',values_to='Values')
datas1=data.frame(Type=ttt1$Celltype,Values=ttt1$Values,Group=rep('Tumor',nrow(ttt1)))

ttt0= GTEx.data %>%
  rownames_to_column('Samples') %>%
  pivot_longer(cols = 2:(ncol(GTEx.data)+1),names_to='Celltype',values_to='Values')
datas.GETx=data.frame(Type=ttt0$Celltype,Values=ttt0$Values,Group=rep('Normal',nrow(ttt0)))

ttt0= TCGA_nor %>%
  rownames_to_column('Samples') %>%
  pivot_longer(cols = 2:(ncol(TCGA_nor)+1),names_to='Celltype',values_to='Values')
datas.nor=data.frame(Type=ttt0$Celltype,Values=ttt0$Values,Group=rep('Normal',nrow(ttt0)))

datas.final=as.data.frame(rbind(datas1,datas.nor,datas.GETx))

pvalues <- sapply(unique(datas.final$Type), function(x) {
  res <- wilcox.test(Values ~ Group, data = datas.final[which(datas.final$Type==x),])
  res$p.value
})

pv <- data.frame(Type = unique(datas.final$Type), pvalue = pvalues)
pv$sigcode <- cut(as.numeric(pv$pvalue), c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '', ''))

p11=ggboxplot(datas.final1,x='Type',y='Values',color='Group',palette = "nejm",size = 0.5,shape=16,
              add = "jitter",add.params=list(size=0.1),bxp.errorbar =T,width=0.55,outlier.shape=NA,ggtheme=theme_bw())
p11=p11+geom_text(aes(Type, y=max(datas.final1$Values)*1.1,label=pv1$sigcode),data=pv1, inherit.aes=F,size=3) + xlab(NULL)+ylab("Genes")

