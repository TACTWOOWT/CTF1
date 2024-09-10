genes=c('CD274','CTLA4','HAVCR2','LAG3','PDCD1','PDCD1LG2','TIGIT','SIGLEC15')
cor.result=rbind()
for (i in 1:nrow(can)) {
  if (can[i,]!='Nor') {
    if (length(grep('GTE',can[i,]))!=1) {
      data=getdatas_more(can[i,],genes)
      dat.gene=getdatas_cli(can[i,],clinical,clinType,dnaName)
      data=as.data.frame(data)
      dat.gene=as.data.frame(dat.gene)
      clinical1=clinical[unique(substr(colnames(data),1,12)),]
      comid=get_cli_data(as.character(clinType),data,can[i,],clinical1)
      data.tum=t(data[,comid])
      dat.gene=dat.gene[,comid]
      dat.gene=as.data.frame(dat.gene)
      #print(dim(dat.gene))
      data.tum=as.data.frame(data.tum)
      ######
      gs=as.numeric(dat.gene)
      cl <- makeCluster(THREAD_NUM)
      # save(gs,file = paste0(BASE_IMG_DIR,'fun_Pancancer_analysis/',fileNamedata1,'/gs.RData'))
      # cor11 = tryCatch({
      #   cor11=get_cor_pva(data.tum = data.tum,gs = gs,cl = cl,path = paste0(BASE_IMG_DIR,'fun_Pancancer_analysis/',fileNamedata1,'/gs.RData'))
      # }, error = function(e) {
      #   print(e)
      # }, finally = {
      #   tryCatch(
      #     stopCluster(cl))
      # }
      # )
      cor11=t(apply(data.tum,2,function(x){
        dd  <- cor.test(as.numeric(x),gs,method="spearman")
        res <- cbind(cor=round(dd$estimate,3),p.value=dd$p.value)
      }))
      cor.result=rbind(cor.result,cbind(Type=rep(can[i,],nrow(cor11)),Genes=rownames(cor11),Cor=cor11[,1],pvalue=cor11[,2]))
    }
  }
}
cor.result=as.data.frame(cor.result)
cor.result$Cor=as.numeric(cor.result$Cor)
cor.result$pvalue=as.numeric(cor.result$pvalue)
cor.result$pvalue[which(cor.result$pvalue==0)]=2.2e-116

cor.result$pstar <- ifelse(cor.result$pvalue < 0.05,ifelse(cor.result$pvalue < 0.01,"**","*"),"")
p1=ggplot(cor.result, aes(Genes,Type)) + 
  geom_tile(aes(fill = Cor),size=1)+
  scale_fill_gradient2(low = mypal[1],mid = "white",high = mypal[2])+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1,size = 10),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))

