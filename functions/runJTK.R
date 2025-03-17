
runJTK=function(output.list, group.by, ...){
  if(!group.by %in% colnames(output.list[[1]])){
    stop('group.by must be a column name of output.list[[1]]!')
  }
  output=matrix(nrow=length(output.list), ncol=nrow(output.list[[1]]))
  rownames(output)=names(output.list)
  colnames(output)=apply(output.list[[1]][,-ncol(output.list[[1]])],1, function(x){paste(x, collapse = '_')})
  for(i in 1:nrow(output)){
    output[i,]=output.list[[i]][,ncol(output.list[[i]])]
  }
  output.meta=output.list[[1]][,-ncol(output.list[[1]])]
  
  JTK=vector(mode = "list", length=length(unique(output.meta[,group.by])))
  names(JTK)=unique(output.meta[,group.by])
  
  for(i in 1:length(JTK)){
    group.i=names(JTK)[i]
    output.i=output[,output.meta[,group.by]==group.i]
    output.meta.i=output.meta[output.meta[,group.by]==group.i,]
    output.i=data.frame(gene=names(output.list),output.i)
    
    write.csv(output.i, file='temp.csv',row.names = FALSE)
    JTK.i=meta2d(infile = 'temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = output.meta.i$ZT, cycMethod = c("JTK"),
                 maxper = 24, minper = 8, combinePvalue = 'bonferroni')$meta
    file.remove('temp.csv')
    JTK[[i]]=cbind(group=group.i, JTK.i)
  }
  return(JTK)
}



runJTK.FUN=function(output.list, group.by, cell.meta, FUN, ...){
  if(!group.by %in% colnames(output.list[[1]])){
    stop('group.by must be a column name of output.list[[1]]!')
  }
  output=matrix(nrow=length(output.list), ncol=nrow(output.list[[1]]))
  rownames(output)=names(output.list)
  colnames(output)=apply(output.list[[1]][,cell.meta],1, function(x){paste(x, collapse = '_')})
  for(i in 1:nrow(output)){
    output[i,]=output.list[[i]][,FUN]
  }
  output.meta=output.list[[1]][,cell.meta]
  
  JTK=vector(mode = "list", length=length(unique(output.meta[,group.by])))
  names(JTK)=unique(output.meta[,group.by])
  
  for(i in 1:length(JTK)){
    group.i=names(JTK)[i]
    output.i=output[,output.meta[,group.by]==group.i]
    if(FUN=='cv' | FUN=='nonzero.mean' | FUN =='burst.size'){
      output.i[is.nan(output.i)]=0
    } else if(FUN=='gini'){
      output.i[is.nan(output.i)]=1
    }
    output.meta.i=output.meta[output.meta[,group.by]==group.i,]
    output.i=data.frame(gene=names(output.list),output.i)
    
    write.csv(output.i, file='temp.csv',row.names = FALSE)
    JTK.i=meta2d(infile = 'temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = output.meta.i$ZT, cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni')$meta
    file.remove('temp.csv')
    JTK[[i]]=cbind(group=group.i, JTK.i)
  }
  return(JTK)
}