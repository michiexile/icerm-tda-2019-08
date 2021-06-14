#setwd('D:\\ATMCS2020\\voting_data\\GCC_politics_voting_19902016')
#year<-1990
#lambda<-0
#pdf('plot_bill.pdf',height=4,width=8)

pdf('plot_bill_sq.pdf',height=5,width=5)
for(year in c(1990,1998,2006)){
  #2006 L2 better
  #1990 L1 betterset
  if(year == 1993){next}
  if(year == 1995){next}
  if(year == 2003){next}
  if(year == 2005){next}
  if(year == 2007){next}
  if(year == 2009){next}
  if(year == 2010){next}
  cat('\n year=',year)
  for(lambda in c(0,0.5,1)){
  cat('\n lambda(alpha)=',lambda)   
  #Read party data
  party<-read.csv(paste0('D:\\ATMCS2020\\politics_data\\us.congress.',year,'.keys.csv'),header = F)
  party<-as.factor(party[,3])
  
  for(i in 1:20){
    destfile<-paste0('us.congress.',year,'.bills_CircularCoordinates_',lambda,'_',i,'.txt')
    if(!file.exists(destfile)){
      break
    }
    N_cocycle<-i
  }
  
  df<-as.data.frame(matrix(NA,nrow = 0,ncol=7))
  df_plot <- df
  
  #Calculate %Democrat
  destbill<-paste0('us.congress.',year,'.bills.txt')
  tmpbill<-read.table(destbill)
  if(year==1998){
    rmv<-which(apply(tmpbill, 2, var)==0)
    tmpbill[1,rmv]<-1
  }
  pcabill <- prcomp(tmpbill, center = TRUE, scale. = TRUE, retx = T)
  
  perDemocrat<-c()
  perRepublican<-c()
  
  for(k in 1:NROW(tmpbill)){
    n_Democrat<-sum( as.numeric(party=='D' & tmpbill[k,]==1) )
    n_Republican<-sum( as.numeric(party=='R' & tmpbill[k,]==1) )
    perDemocrat<-c(perDemocrat,n_Democrat/sum(party=='D'))
    perRepublican<-c(perRepublican,n_Republican/sum(party=='R'))
  }

  
  for(j in 0:N_cocycle){
    destfile<-paste0('us.congress.',year,'.bills_CircularCoordinates_',lambda,'_',j,'.txt')
    tmp<-read.table(destfile)
    N_total<-NROW(tmp)
    
    #destbill<-paste0('us.congress.',year,'.bills.txt')
    #tmpbill<-read.table(destbill)
    
    
    #plot(pcabill$x[,1:2])
    df<-rbind(df,cbind(1:N_total,tmp,j,pcabill$x[,1],pcabill$x[,2],perDemocrat,perRepublican))
    
  }
  colnames(df)<-c('bill','GCC','cocycle','PCA1','PCA2','Democrat','Republican')
  df$cocycle<-as.factor(df$cocycle)
  
  write.table(df,paste0('us.congress.',year,'.bills.plotting'))
  library(ggplot2)
  sharedTitle<-paste0('Year ',year,' GCC (mod 1) with penalty=',1-lambda,'*L^1+',lambda,'*L^2')
  
  gp1<-ggplot()+ggtitle(sharedTitle)+
    geom_point(data=df,aes(x=bill,y=GCC,colour=cocycle),size=.5)+
    #geom_col(data=df,aes(x=voter,y=GCC,fill=cocycle),position='stack')+
    #geom_point(data=as.data.frame(matrix(NA,nrow=N_total)),aes(x=1:N_total,y=rep(0,N_total),colour=party),size=3)+
    theme_bw()+xlim(0,range(df$bill)[2])+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
    theme(plot.title = element_text(hjust = 0.5))
  #+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))
  
  df1<-cbind(1:N_total,0)
  for(k in unique(df$cocycle)){
    df1[,2]<-df1[,2]+df[df$cocycle==k,]$GCC
  }
  df1<-as.data.frame(df1)
  colnames(df1)<-c('voter','GCC')
  #gp2<-ggplot()+#geom_line(data=df1,aes(x=voter,y=GCC))+
  #  geom_point(data=as.data.frame(matrix(NA,nrow=N_total)),aes(x=1:N_total,y=df1$GCC%%1,colour=df1%cocycle),size=.3)+
  #  theme_bw()+xlim(0,450)+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
  #  theme(plot.title = element_text(hjust = 0.5))+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))
  
  gp3<-ggplot()+ggtitle(sharedTitle)+
    geom_density(aes(x=df1$GCC%%1,y=..scaled..),size=1)+geom_density(data=df,aes(x=GCC,y=..scaled..,group=cocycle,color=cocycle),adjust=2)+
    theme_bw()+xlim(0,1)+ylim(-0.05,1.05)+xlab('GCC')+ylab('density')+
    theme(plot.title = element_text(hjust = 0.5))
  
  hsv = colorspace::HSV(seq(0, 360, length=13)[-13], 1, 1)
  
  gp4<-ggplot()+
    geom_point(data=df,aes(x=PCA1,y=PCA2,colour=Democrat),size=2)+ggtitle(sharedTitle)+
    #geom_col(data=df,aes(x=voter,y=GCC,fill=cocycle),position='stack')+
    #geom_point(data=as.data.frame(matrix(NA,nrow=N_total)),aes(x=1:N_total,y=rep(0,N_total),colour=party),size=3)+
    theme_bw()+#+xlim(0,range(df$bill)[2])+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
    theme(plot.title = element_text(hjust = 0.5))+scale_colour_gradient(low='red',high='blue',limits=c(0,1),"D%")
  #+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))

  gp5<-ggplot()+#geom_line(data=df1,aes(x=voter,y=GCC))+
    geom_point(data=df,aes(x=PCA1,y=PCA2,colour=GCC),size=2)+ggtitle(sharedTitle)+
    theme_bw()+#xlim(0,450)+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
    theme(plot.title = element_text(hjust = 0.5))+scale_colour_gradientn(colours=rainbow(100),limits=c(0,1))#+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))
  
  gp6<-ggplot()+#geom_line(data=df1,aes(x=voter,y=GCC))+
    geom_point(data=df,aes(x=Democrat,y=Republican, colour=GCC),size=2)+ggtitle(sharedTitle)+
    theme_bw()+xlim(-0.05,1.05)+ylim(-0.05,1.05)+xlab('% of all Democrats favored')+ylab('% of all Republicans favored')+
    theme(plot.title = element_text(hjust = 0.5))+scale_colour_gradientn(colours=rainbow(100),limits=c(0,1))#+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))
  
    
  #gridExtra::grid.arrange(gp1,gp2,gp3,ncol=1,top=paste0('Year ',year,' GCC (mod 1) with penalty=',1-lambda,'*L^1+',lambda,'*L^2'))
  print(gp1)
  #print(gp2)
  print(gp3)
  print(gp4)
  print(gp5)
  print(gp6)
  
  
  }
}
dev.off()