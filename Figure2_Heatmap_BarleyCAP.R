rm(list=ls())
#setwd("~/Desktop")  #Set WD
heatmap<-read.csv("Heatmap_R_input.csv", header=T)

library(RColorBrewer)
library(gplots)
##create a red color palette
ylred <- brewer.pal(9, "Reds")

hits_matrix<-data.matrix(heatmap)
hits=hits_matrix[,-1]

rownames(hits) <- heatmap$QTL

#For assigning colors to the CHR legend
#rlab= c(rep("white",2), rep("white",8), rep("white",6), rep("white",1), rep("white",25))

#Defining vector for adding spaces/lines in heatmap (see function "rowsep")
row<-c(0,2,5,7,8,12)
#col <-c(1,2,3,4,5,6,7,8,9,10,11)
#Creating vectors for legends
#names<-c("1H","2H","3H","4H","6H","UNK") #legend names
#col.legend<-c("blue","dodgerblue","blue","dodgerblue","blue","dodgerblue") #legend #colors 

pdf("Heatmap-QTLs.pdf", width=8, height=7) #Size of the graphics output

#set margins of plots
par(mar=c(2,4,4,4)) # bottom, left, top, right
#Creating a heatmap
heatmap.2(hits,
           # dendrogram control
           Rowv = NULL,
           Colv=TRUE,
           dendrogram = "none",
           scale ="none",
          #image plot
          revC = identical("Colv", "Rowv"),
           na.rm=TRUE,
           col=ylred,
          # block sepration
          rowsep=row,
          #colsep = col,
          sepcolor="black",
          sepwidth=c(0.0001,0.0001),
          trace="none",
          # Row/Column Labeling
          margins = c(8,14),
          #RowSideColors=rlab,
          labRow = heatmap$QTL,
          labCol = NULL,
          srtRow = 0,
          srtCol = 60,
          adjRow = c(0,NA),
          offsetRow = 0.00,
          offsetCol = 0.75,
          # color key + density info 
          key = TRUE,
           keysize = 1,  
           key.title = "Color key",
            key.xlab = "-log(P) value", 
            key.ylab = "Frequency",
           #plot labels,
           xlab = "Association mapping panels",
           ylab = "QTLs")
dev.off()



