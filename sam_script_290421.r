test=read.csv("../OneDrive - Imperial College London/yeast_10x/WILD_2021/counts_norm_by_strain_t20_cont_filt.csv",stringsAsFactors = F)
row.names(test)=test[,1]
test=test[,-1]
#
test1=test[which(rowSums(abs(log2(test) > 2)) != 0),]
#
test2=test1[-grep("NCRNA",row.names(test1)),]
#
a=pheatmap(log2(test2),cluster_cols = F,breaks=seq(-7.5,7.5,1),col=colorRampPalette(c("dark green","grey","dark red"))(15),clustering_method = "ward.D2",cutree_rows = 10,show_rownames = F)
b=cutree(a$tree_row,k=10)
BIG[which(BIG[,1] %in% names(b)[which(b==5)]),c("common.BIG","Annot")]
boxplot(log2(test2),outline = F,las=2,cex.names=0.1)
#
#
#
#hit=BIG[grep("^rpl",BIG[,"common.BIG"]),1]
hit=row.names(GOtable)[which(GOfinal[,99]==1)]
hit=row.names(GOtable)[which(GOfinal[,100]==1)]
hit=names(b)[which(b==10)]
plot_gene(gene=hit,dat=test,yli=c(-5,5))
#




plot_gene=function(gene=NULL,dat,yli=c(-5,15))
{
if(length(gene) > 1)
{
 dat=colMeans(dat[which(row.names(dat) %in% gene),])
 dat=rbind(dat,dat)
 gene=1
}
dat=log2(dat)
par(mfrow=c(1,2),mar=c(4,2,4,1))
wild=c(758,759,762,838)
my.col=c("dark red","dark green","dark blue","orange")
plot(c(0,20,60),dat[gene,c(grep(paste("WT","_20_cont",sep=''),colnames(dat)),grep("WT.*L",colnames(dat)))],type="b",ylim=yli,ylab="fold change",xlab="Time",main="Low")
abline(h=dat[gene,grep(paste("WT","_60_cont",sep=''),colnames(dat))],lty=2,col="dark red")
for(i in 1:4)
{
  points(c(0,20,60),dat[gene,c(grep(paste(wild[i],"_20_cont",sep=''),colnames(dat)),grep(paste(wild[i],".*L",sep=''),colnames(dat)))],type="b",col=my.col[i],pch=20)
 abline(h=dat[gene,grep(paste(wild[i],"_60_cont",sep=''),colnames(dat))],lty=2,col=my.col[i])
}
abline(h=0)
legend(x="topleft",legend = gene,bty="n")
abline(h=dat[gene,grep(paste("WT","_60_cont",sep=''),colnames(dat))],lty=2)
plot(c(0,20,60),dat[gene,c(grep(paste("WT","_20_cont",sep=''),colnames(dat)),grep("WT.*M",colnames(dat)))],type="b",ylim=yli,ylab="fold change",xlab="Time",main="Medium")
for(i in 1:4)
{
 points(c(0,20,60),dat[gene,c(grep(paste(wild[i],"_20_cont",sep=''),colnames(dat)),grep(paste(wild[i],".*M",sep=''),colnames(dat)))],type="b",col=my.col[i],pch=20)
 abline(h=dat[gene,grep(paste(wild[i],"_60_cont",sep=''),colnames(dat))],lty=2,col=my.col[i])
}
abline(h=0)
legend(x="topright",legend = wild,text.col = my.col,bty="n")
#print(BIG[which(BIG[,1]==gene),c("common.BIG","Annot")])
}




