load('TCGA_GBM.RData')
tcga.data=log(dt+1)
load('Clinical.RData')
comsamples=intersect(c(paste0(rownames(clinical),'-01'),paste0(rownames(clinical),'-03'),paste0(rownames(clinical),'-06')),colnames(tcga.data))
gs=c('XX','XX')
data=tcga.data[gs,comsamples]
library(ggcorrplot)
cor1 = cor(data)
aaa1=pheatmap(cor1,show_colnames = F)
inds11=aaa1$tree_row$order
corr_table=cor(dt[,inds11], use = "p")
p1=ggcorrplot(corr_table, method = 'circle',type = "lower", show.diag = T,lab=T ,hc.order = T,colors = c(mypal[2],'white',mypal[1]))

