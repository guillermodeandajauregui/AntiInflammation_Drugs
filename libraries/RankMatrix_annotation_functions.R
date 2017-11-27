#######################################
#Functions for annotating RankedMatrix
#######################################

#######################################
#Libraries
#######################################
library("data.table")
library("annotate")
library("hgu133a.db")
#######################################
#ReadMatrix 
#Reads a Ranked Matrix
#Rows are probesets
#Cols are samples
#######################################
ReadMatrix<- function(FilePath, DF=TRUE){
  
  RankMatrix = fread(FilePath, data.table = FALSE)
  rownames(RankMatrix)<-RankMatrix[,1]
  colnames(RankMatrix)<-RankMatrix[1,]
  RankMatrix<-RankMatrix[-1,-1]
  
  if(DF == TRUE){
    return(as.data.frame(RankMatrix))
  }
  else{
    return(as.matrix(RankMatrix))
  }
}

#######################################
#AnnotateAggregator
#Aggregates a Ranked Matrix 
#With Affy probes as rownames
#into gene symbols
#Based on a chip database
#Either Max or Median aggregated
#######################################

AnnotateAggregator <- function(RankMatrix, Chip="hgu133a.db", Method="Median"){
  x = as.data.frame(RankMatrix) #make sure RankMatrix behaves like DataFrame
  genesymbols<-data.frame() #will contain genesymbols
  genesymbols<-(getSYMBOL(as.character(rownames(x)), 
                          Chip)
                ) #takes symbols from Affymetrix database; if no gene symbol in annotation, returns NA
  sym<-as.data.frame(genesymbols)
  x$GS<-genesymbols #will add a column with gene symbols ; if no gene symbol in annotation, returns NA
  x$GS<-ifelse(is.na(x$GS), 
               as.vector(rownames(x)), 
               as.vector(x$GS)
               ) #replace NA with AffyID
  if(Method=="Median"){
    x.med <- aggregate(. ~ GS, data = x, median) #Agregar valores de Genesymbol por mediana#
    rownames(x.med)<-x.med$GS
    y.med<-x.med[,-1]
    return(y.med)
  }
  else{
    if(Method=="Max"){
      x.max <- aggregate(. ~ GS, data = x, max) #Agregar valores de Genesymbol por mediana#
      rownames(x.max)<-x.max$GS
      y.max<-x.max[,-1]
      return(y.max)
    }
  }
}

#######################################
#CustomAggregator
#Aggregates a Ranked Matrix 
#With Affy probes as rownames
#into gene symbols
#Based on a Custom Annotation
#Either Max or Median aggregated
#######################################

CustomAggregator <- function(RankMatrix, CustomAnnotation, Method="Median"){
  x = as.data.frame(RankMatrix) #make sure RankMatrix behaves like DataFrame
  genesymbols<-data.frame()
  #CustomAnnotation must be a data frame with column "AffyID" and column "Symbol"
  genesymbols<-cbind(CustomAnnotation$AffyID[match(x = rownames(x), 
                                                         table = Merged_Annotation_Data$AffyID)],
                     Merged_Annotation_Data$Symbol[match(x = rownames(x), 
                                                         table = Merged_Annotation_Data$AffyID)]
                    )
  sym<-as.data.frame(genesymbols, 
                     stringsAsFactors = FALSE)
  rownames(sym)<-sym[,1]
  sym<-sym[,-1, drop = FALSE]
  rownames(x)<-rownames(sym)
  x$GS<-sym$V2
  x$GS<-ifelse(is.na(x$GS), 
               as.vector(rownames(x)), 
               as.vector(x$GS)
  ) #replace NA with AffyID
  if(Method=="Median"){
    x.med <- aggregate(. ~ GS, data = x, median) #Agregar valores de Genesymbol por mediana#
    rownames(x.med)<-x.med$GS
    y.med<-x.med[,-1]
    return(y.med)
  }
  else{
    if(Method=="Max"){
      x.max <- aggregate(. ~ GS, data = x, max) #Agregar valores de Genesymbol por mediana#
      rownames(x.max)<-x.max$GS
      y.max<-x.max[,-1]
      return(y.max)
    }
  }
}

#######################################
#ReRank
#Takes an Aggregated Rank matrix 
# ReRanks the new aggregated values
#Wrapper for apply Rank over columns
#######################################
ReRank <-function(RankMatrix, ties.method="first"){
  y= apply(X = RankMatrix, 
        MARGIN = 2, 
        FUN = rank, 
        ties.method=ties.method)
  return(y)
}
