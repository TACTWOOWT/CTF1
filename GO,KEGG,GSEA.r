
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GOplot)


读入数据rt=read.table("data.txt",sep="\t",check.names=F)
转换为ENTREZID
entrezIDs = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
使用entrezIDs gene<- entrezIDs$ENTREZID



##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #
##可视化
#条形图
pdf(file="GO-barplot.pdf",width = 10,height = 15)
barplot(go, drop = TRUE, showCategory =10,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
##气泡图
pdf(file="GO-bubble.pdf",width = 10,height = 15)
dotplot(go,showCategory = 10,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         

##可视化
##条形图
pdf(file="KEGG-barplot.pdf",width = 10,height = 13)
barplot(kk, drop = TRUE, showCategory = 15,label_format=100)
dev.off()
##气泡图
pdf(file="KEGG-bubble.pdf",width = 10,height = 13)
dotplot(kk, showCategory = 15,label_format=100)
dev.off()


##### gsea富集 ####
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,                  
                                             organism     = organism, #人hsa 鼠mmu                   
                                             pvalueCutoff = 0.25)  #实际为padj阈值可调整 
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez,                              
                                                 OrgDb=OrgDb,                             
                                                 keyType='ENTREZID')#转化id               

GO_kk_entrez <- gseGO(geneList     = geneList,               
                                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"               
                                      OrgDb        = OrgDb,#人类org.Hs.eg.db 鼠org.Mm.eg.db               
                                      keyType      = "ENTREZID",              
                                      pvalueCutoff = 0.25)   #实际为padj阈值可调整

GO_kk <- DOSE::setReadable(GO_kk_entrez,                            
                                              OrgDb=OrgDb,                           
                                              keyType='ENTREZID')#转化id 
save(KEGG_kk_entrez, GO_kk_entrez, file = "GSEA_result.RData")

##选取富集结果
kk_gse <- KEGG_kk
kk_gse_entrez <- KEGG_kk_entrez

###条件筛选 
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]

#### 经典的GSEA图 
up_gsea$Description
i=2
gseap1 <- gseaplot2(kk_gse,                    
                                 up_gsea$ID[i],#富集的ID编号                    
                                 title = up_gsea$Description[i],#标题                   
                                 color = "red", #GSEA线条颜色                    
                                 base_size = 20,#基础字体大小                   
                                 rel_heights = c(1.5, 0.5, 1),#副图的相对高度                    
                                 subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图                    
                                 ES_geom = "line", #enrichment score用线还是用点"dot"                    
                                 pvalue_table = T) #显示pvalue等信息
ggsave(gseap1, filename = 'GSEA_up_1.pdf', width =10, height =8)





