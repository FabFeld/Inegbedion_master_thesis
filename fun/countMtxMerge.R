countMtxMerge = function(cmtx.ls){
  mergeMtx = function(df1,df2){merge(df1,df2, by=0, all = T)}
  mtx = Reduce(mergeMtx, cmtx.ls)
  row.names(mtx) <- mtx$Row.names
  return(mtx[,-1])
}