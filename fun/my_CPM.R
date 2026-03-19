# This custom function is used, because cpm() is from edgeR package, which does not to be installed this way


## CPM transformation function 
my_CPM = function(mtx) {
  cpm_mtx = apply(mtx, 2, function(x) (x/sum(x)) * 1e6)
  return(as.data.frame(cpm_mtx))
}
