for(year in 1990:2016){
  destfile<-paste0('us.congress.',year,'.votes.txt')
  if(!file.exists(destfile)){
    next
  }else{
    test1<-read.table(destfile)
    outfile<-paste0('us.congress.',year,'.bills.txt')
    write.table(t(test1),file=outfile,col.names = F,row.names = F)
  }
  
}