importCountMtx = function(file){
  if(any(grepl(".txt",file))==T) { # import txt
    tmp1 = file[grepl(".txt",file)]
    count.matrix = read.table(file = tmp1, row.names = 1, header = T)
    if(ncol(count.matrix) <= 1) {
      count.matrix = read.delim2(file = tmp1, row.names = 1, header = T)
      if(ncol(count.matrix) <= 1) {
        count.matrix = read.delim(file = tmp1, row.names = 1, header = T)
      }
    }
  } else { # import csv
    tmp1 = file[grepl(".csv",file)]
    count.matrix = read.csv2(tmp1, header = T, row.names = 1, fill = F)
    if(ncol(count.matrix) <= 1) {
      count.matrix = read.csv(tmp1, header = T, row.names = 1, fill = F)
    }
  }
  return(count.matrix)
}