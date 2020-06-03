##将学出的embendding初始化：
breast_emb=read.table("data/breast_out.txt",header=FALSE)
load("data/breast_gene.Rdata")
t1=breast_emb[order(breast_emb[,1]),]
rownames(t1)=t1[,1]
t1=t1[,-1]
breast_emb=t1
rm(t1)
breast_mean=apply(breast_emb,2,mean)
breast_var=apply(breast_emb,2,var)
breast_zscore <- matrix(data=0,nr=length(breast_gene),nc=128)
for(i in 1:128)
{for(j in 1:length(breast_gene))
{z=(breast_emb[j,i]-breast_mean[i])/breast_var[i]
breast_zscore[j,i]=z
}
}



breast_list <- list()
for(i in 1:128)
{gene <- breast_gene
score <- breast_zscore[,i]
groupi <- data.frame(gene,score) 
breast_list[[i]] <- groupi
rm(groupi)
}

#将gene赋给breast_list[[i]]行名
for(i in 1:128)
{
  #Groups_list[[i]]=Groups_list[[i]][order(Groups_list[[i]]$gene),]
  rownames(breast_list[[i]])=breast_gene
}


#rm(t_value,breast_list1),调参确定最佳阈值
t_value <- rep(0 , times =128, length.out =128, each = 1)
#lung-500,prostate-1,breast-5
for(i in 1:128)
{t_value[i]=(breast_var[i]*breast_var[i])/(1+breast_var[i]*breast_var[i])*5*breast_var[i]}

breast_list1=list()
for(i in 1:128)
{breast_list1[[i]]=breast_list[[i]][which(breast_list[[i]][,2]>t_value[i]),]}