qpcr <- function(dat, numsamp=12, numgene=1, rep=2, col=3, genenames=NA, treatment=NA, savegraphs=F, drawheatmap=T){
  getData <- function(){
    #numsamp=number of samples per repetition, numgene=number of genes assessed, rep=the repetitions of each gene assessed, col=the collumn of the original data table containing the cf value
    c=1
    g=1
    t=1
    s=0
    skip=0
    if(numsamp>12){
      skip=12-(numsamp%%12)
    }
    else if(numsamp<12){
      skip=12-numsamp
    }
    out = matrix(NA, numsamp, numgene*rep)
    for(i in dat$Well){
      if(s>0){
        s=s-1;
      }
      else{
        out[c,g] = dat[t,col]
        c=c+1
        if(c>=numsamp+1){
          c=1
          g=g+1
          s=skip
        }
      }
      t=t+1
    }
    out
  }
  meancf <- function(dat){
    out=matrix(NA, nrow(dat),1)
    c=1
    while(c<=nrow(dat)){
      out[c,1]=mean(dat[c,1:ncol(dat)], na.rm=TRUE)
      c=c+1
    }
    out
  }
  processdata<-function(data){
    c=0
    genes=matrix(NA,nrow(data),numgene)
    while(c<numgene){
      #str(extract[,(c*rep+1):(c*rep+rep)])
      genes[,c+1]<-meancf(data[,(c*rep+1):(c*rep+rep)])
      c=c+1
    }
    genes
  }
  qPCRHeatmap <- function(dat){
    require(lattice)
    dataMatrix = matrix(dat$Cq, 12, 8)
    mat2 = dataMatrix[,ncol(dataMatrix):1]
    graph = levelplot(mat2, col.regions=heat.colors(96, 1), main="Plate CF Values", scales=list(draw=F), xlab="", ylab="")
    if(savegraphs){
      png("img/heatmap.png")
      print(graph)
      dev.off()
    }
    else{
      graph
    }
  }
  if(drawheatmap){
    qPCRHeatmap(dat)
  }
  extract = getData()
  out = processdata(extract)
  if(!is.na(genenames)){
    colnames(out)<-genenames
  }
  out=as.data.frame(out)
  if(!is.na(treatment)){
    if(length(treatment)!=numsamp){
      message("Could not apply treatment names to table, because number of treatments does not equal the number of samples")
    }
    else{
      out=cbind(treatment, out)
    }
  }
  c=1
  dir.create("img", showWarnings = FALSE)
  for(i in genenames){
    v=out[[c+1]]
    if(savegraphs){
      png(paste("img/",genenames[c],".png",sep=""))
      boxplot(v~out$treat)
      title(i)
      c=c+1
      dev.off()
    }
    else{
      boxplot(v~out$treat)
      title(i)
      c=c+1
    }
  }
  write.csv(out, "qpcrCfMeans.csv")
  out
}
extract = qpcr(dat, 14,2,2,3,c("GAPDH", "UCP-1"),savegraphs=T,treatment=c("WT","WT","WT","WT","IL-4","IL-4","IL-4","Thio","Thio","Thio","IL-4+Thio","IL-4+Thio","IL-4+Thio","IL-4+Thio"))