setwd('/media/hrluo/ALL/icerm-tda-2019-08/2021-06/SECTION4_2')
#year<-1990
#lambda<-0
library(clusterCrit)
pdf('plot_voting_Separate.pdf',height=4,width=8)
for(year in c(1990,1998,2006)){
  #2006 L2 better
  #1990 L1 better
  if(!(year %in% c(1990,1998,2006)) ){next}

  cat('\n year=',year)
  for(lambda in c(0,0.5,1)){
  cat('\n lambda(alpha)=',lambda)   
  #Read party data
  party<-read.csv(paste0('us.congress.',year,'.keys.csv'),header = F)
  party<-as.factor(party[,3])
  
  N_cocycle<-0
  for(i in 0:20){
    destfile<-paste0('us.congress.',year,'.votes_CircularCoordinates_',lambda,'_',i,'.txt')
    if(!file.exists(destfile)){
      break
    }
    N_cocycle<-i
  }
  
  df<-as.data.frame(matrix(NA,nrow = 0,ncol=2))
  for(j in 0:N_cocycle){
  #Pick only the longest 1-cocycle with persistence greater than 1 manually
  #if(year==1990){year_leading_vec=c(8,1,10)}
  #if(year==1998){year_leading_vec=c(5,2,4)}
  #if(year==2006){year_leading_vec=c(2,4,3)}
  #for(j in year_leading_vec){
    destfile<-paste0('us.congress.',year,'.votes_CircularCoordinates_',lambda,'_',j,'.txt')
    tmp<-read.table(destfile)
    N_total<-NROW(tmp)
    df<-rbind(df,cbind(1:N_total,tmp,j,party))
  }
  colnames(df)<-c('voter','GCC','cocycle','party')
  df$cocycle<-as.factor(df$cocycle)
  #df_GCC_mod1<-cbind( abs(df$GCC),abs(df$GCC+1) )
  #df$GCC<-apply(df_GCC_mod1, 1, FUN=min)
  df$GCC<-df$GCC%%1
  df$GCC_x<-cos(df$GCC*2*pi)
  df$GCC_y<-sin(df$GCC*2*pi)
  
  library(ggplot2)
  #intIdx <- clusterCrit::intCriteria(cbind(df$GCC_x,df$GCC_y),as.integer(df$party),"all")
  #extIdx <- clusterCrit::extCriteria()
  
  #sharedTitle<-paste0('Year ',year,' GCC (mod 1) with penalty=',1-lambda,'*L^1+',lambda,'*L^2','\n DBI=',round(intIdx$davies_bouldin,3),' CHI=',round(intIdx$calinski_harabasz,3),' SIL=',round(intIdx$silhouette,3))
  sharedTitle<-'Single'
  gp1<-ggplot()+geom_point(data=df,aes(x=voter,y=GCC,shape=cocycle,colour=party, fill=party),size=3,pch=21)+
    #geom_col(data=df,aes(x=voter,y=GCC,fill=cocycle),position='stack')+
    #geom_point(data=as.data.frame(matrix(NA,nrow=N_total)),aes(x=1:N_total,y=rep(0,N_total),colour=party, fill=party),size=3)+
    theme_bw()+xlim(0,450)+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
    theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))+scale_fill_manual(values = c("D" = "blue", "I" = "green","R"="red"))+ggtitle(sharedTitle)
  
  df1<-cbind(1:N_total,0)
  for(k in unique(df$cocycle)){
    df1[,2]<-df1[,2]+df[df$cocycle==k,]$GCC
  }
  df1<-as.data.frame(df1)
  colnames(df1)<-c('voter','GCC')
  
  df1$GCC<-df1$GCC%%1
  df1$GCC_x<-cos(df1$GCC*2*pi)
  df1$GCC_y<-sin(df1$GCC*2*pi)
  df1$party<-party
  library(ggplot2)
  gc()
  intIdx1 <- clusterCrit::intCriteria(cbind(df1$GCC_x,df1$GCC_y),as.integer(df1$party),"all")
  #extIdx <- clusterCrit::extCriteria()
  sharedTitle<-paste0('Year ',year,' GCC (mod 1) with penalty=',1-lambda,'*L^1+',lambda,'*L^2','\n DBI=',round(intIdx1$davies_bouldin,3),' CHI=',round(intIdx1$calinski_harabasz,3),' TAU=',round(intIdx1$tau,3))
  
  
  gp2<-ggplot()+#geom_line(data=df1,aes(x=voter,y=GCC))+
    geom_point(data=as.data.frame(matrix(NA,nrow=N_total)),aes(x=1:N_total,y=df1$GCC%%1,colour=party, fill=party),size=3,pch=21)+
    theme_bw()+xlim(0,450)+ylim(-0.05,1.05)+xlab('voter')+ylab('GCC')+
    theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))+scale_colour_manual(values = c("D" = "blue", "I" = "green","R"="red"))+scale_fill_manual(values = c("D" = "blue", "I" = "green","R"="red"))+ggtitle(sharedTitle)
  
  gp3<-ggplot()+geom_density(aes(x=df1$GCC%%1,y=..scaled..),size=1)+geom_density(data=df,aes(x=GCC,y=..scaled..,group=cocycle,color=cocycle),adjust=2)+
    theme_bw()+xlim(0,1)+ylim(-0.05,1.05)+xlab('GCC')+ylab('density')+
    theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))+ggtitle(sharedTitle)
    
  #gridExtra::grid.arrange(gp1,gp2,gp3,ncol=1,top=paste0('Year ',year,' GCC (mod 1) with penalty=',1-lambda,'*L^1+',lambda,'*L^2'))
  print(gp1)
  print(gp2)
  print(gp3)
  
  }
}
dev.off()