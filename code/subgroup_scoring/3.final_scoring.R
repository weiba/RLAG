#rm(list=ls())
load('data/final_result_breast.Rdata')
load('data/gene_sample_number.Rdata')
load('data/breast_sample_genelist.Rdata')
final_result_breast=as.data.frame(final_result_breast)
final_result_breast[,2]=as.numeric(as.character(final_result_breast[,2]))
x <- seq(0, 1, by = 0.1)


finalscores=function()
{
  n=1
  new_breast_final_score=list()
  
  rownames(final_result_breast)=final_result_breast[,1]
  for(i in 1:length(breast_sample_genelist))
  {
    temp1=final_result_breast[breast_sample_genelist[[i]],]
    score=c()
    for(k in 1:length(breast_sample_genelist[[i]]))
    {score[k]=alpha*temp1[k,2]+(1-alpha)*log10((gene_sample_number[[i]][k]+1)/length(breast_sample_genelist))}
    temp1[,2]=score
    new_breast_final_score[[n]]=temp1
    n=n+1
   
  }
  
  

  load("data/breast_gene.Rdata")
  gene=breast_gene
  final_breast_score=c()
  for(i in 1:length(new_breast_final_score))
  {
    final_breast_score=rbind(final_breast_score,new_breast_final_score[[i]])
  } 
  
  
  total_gene=c()
  for(i in 1:length(gene)) 
  {
    gene_score=max(as.numeric(final_breast_score[which(as.character(final_breast_score[,1])==gene[i]),2]))
    total_gene=rbind(total_gene,cbind(gene[i],gene_score))      
    
  }
  total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
  breast_final_result=total_gene1
  
}



for(j in 1:length(x))
{ 
  alpha=x[j]
  breast_final_result=finalscores()
  t1=paste0("data/breast_",x[j],"_results.csv")
  write.csv(breast_final_result,file = t1,row.names = FALSE)
  
}

