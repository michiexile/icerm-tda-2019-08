library(RcppCNPy)
library(ggplot2)
pdf('collection_5_plot.pdf',height=2,width=6)
for (alpha in c(0,0.5,1)){
        fnam<-paste0("collection_5_tshort_CircularCoordinates_",alpha,"_0.txt")
        fmat <- read.table(fnam)
        df   <- as.data.frame(cbind(1:NROW(fmat),fmat))
        colnames(df) <- c('X','Y')
        gf   <- ggplot(df)+geom_point(data=df,aes(x=X,y=Y),fill='gray',size=3,pch=21)+
                ggtitle(paste0('Circular coordinates with ',(1-alpha),'*L^1+',(alpha),'*L^2'))+
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
                xlab('Index')+ylab('Circular coordinates')+theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))
        print(gf)
}
dev.off()

pdf('collection_4_plot.pdf',height=2,width=6)
for (alpha in c(0,0.5,1)){
        fnam<-paste0("collection_4_tshort_CircularCoordinates_",alpha,"_0.txt")
        fmat <- read.table(fnam)
        df   <- as.data.frame(cbind(1:NROW(fmat),fmat))
        colnames(df) <- c('X','Y')
        gf   <- ggplot(df)+geom_point(data=df,aes(x=X,y=Y),fill='gray',size=3,pch=21)+
                ggtitle(paste0('Circular coordinates with ',(1-alpha),'*L^1+',(alpha),'*L^2'))+
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
                xlab('Index')+ylab('Circular coordinates')+theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))
        print(gf)
}
dev.off()

pdf('collection_3_plot.pdf',height=2,width=6)
for (alpha in c(0,0.5,1)){
        fnam<-paste0("collection_3_tshort_CircularCoordinates_",alpha,"_0.txt")
        fmat <- read.table(fnam)
        df   <- as.data.frame(cbind(1:NROW(fmat),fmat))
        colnames(df) <- c('X','Y')
        gf   <- ggplot(df)+geom_point(data=df,aes(x=X,y=Y),fill='gray',size=3,pch=21)+
                ggtitle(paste0('Circular coordinates with ',(1-alpha),'*L^1+',(alpha),'*L^2'))+
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
                xlab('Index')+ylab('Circular coordinates')+theme(axis.title=element_text(size=16),plot.title = element_text(hjust = 0.5, size=18, face="bold"))
        print(gf)
}
dev.off()

