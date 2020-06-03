load('data/Breast.Rdata')#Prostate/Lung/Breast
load("data/breastw_dis.Rdata")
load("data/breast_gene.Rdata")
load('data/breast_list1.Rdata')
colnames(patMutMatrix)=toupper(colnames(patMutMatrix))

main_function=function(breastw_dis)
{
  influ=breastw_dis
  tt=normalize_grph_freq(influ,patMutMatrix)
  ent=compute_entrophy(tt$total_graph,tt$frequency)
  final_result_breast=combine_sort_new(ent$subgroup,influ)
  final_result_breast
}

normalize_grph_freq=function(total_graph,patMatrix)
{ 
  gen_fr=colSums(patMatrix)
  gene_frequence=gen_fr/nrow(patMatrix)
  gene_frequence=cbind(names(gene_frequence),gene_frequence)
  list(total_graph=total_graph,frequency=gene_frequence)
}

compute_entrophy=function(total_graph,frequency)
{
  group=list()
  for(i in 1:length(breast_list1))
  {
    gene_num <- row.names(breast_list1[[i]]) 
    group_i <- rep(i , times =nrow(breast_list1[[i]]), length.out =nrow(breast_list1[[i]]), each = 1)
    group[[i]] <- data.frame(Gene=gene_num,Groups=group_i)
    group[[i]] <- as.matrix(group[[i]])
  }
  H_G=c()
  for(j in 1:length(group))###calculate the entrophy for each sub-group
  { 
    
    sub_graph=total_graph[group[[j]][,1],group[[j]][,1]]
    i_entropy=c()
    
    for(k in 1:nrow(sub_graph))
    {
      n=length(which(sub_graph[k,]!=0))
      #wij=sub_graph[k,]
      wij=sub_graph[k,which(sub_graph[k,]!=0)]
      yj=frequency[colnames(sub_graph)[which(sub_graph[k,]!=0)],]
      if(is.null(dim(yj))==TRUE)
      {Wi=sum(wij*log(as.numeric(frequency[row.names(sub_graph)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_graph)[k],2]))}
      else
      {Wi=sum(wij*log(as.numeric(frequency[row.names(sub_graph)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_graph)[k],2]))}
      
     
      Pi=n/(length(which(sub_graph!=0))/2)
      Ei=Wi*Pi*log2(1/Pi)
      if(n==0||Pi==0)
      {Ei=0}
      i_entropy=rbind(i_entropy,entropy=Ei)
    } 
    group[[j]]=cbind(group[[j]],entrophy=i_entropy)
    H_k=sum(as.numeric(group[[j]][,3]))
    H_G=rbind(H_G,cbind(as.character(group[[j]][1,2]),H_k))
  }
  list(subgroup=group,H_G=H_G)
}

combine_sort_new=function(finalresult,ECC_matrix)
{
  finalresult1=c()
  for(i in 1:length(finalresult))
  {
    finalresult1=rbind(finalresult1,finalresult[[i]])
  }
  
  
  ######below is weight of the compartment edges
  comp=c(1:length(finalresult))
  comp=as.character(comp)
  com=c()
  
  for(i in 1:length(comp))
  {
    gen_nam=finalresult1[which(comp[i]==finalresult1[,2]),1]
    com_graph=ECC_matrix[gen_nam,gen_nam]
    com_value=length(which(com_graph!=0))/length(which(ECC_matrix!=0))
    com=rbind(com,cbind(as.character(comp[i]),com_value))
  }
  #########
  for(j in 1:nrow(finalresult1))
  {
    finalresult1[j,3]=as.numeric(com[which(finalresult1[j,2]==com[,1]),2])*as.numeric(finalresult1[j,3])
  }
  gene=unique(finalresult1[,1])
  total_gene=c()
  for(i in 1:length(gene)) 
  {
    G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))
    total_gene=rbind(total_gene,cbind(gene[i],G_v))      
    ##128组中gene最大的评分（gene+最大评分）
  }
  total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]     ##gene最大评分降序排列，第一列为gene序号，第二列为最大评分
  return(total_gene1)
}