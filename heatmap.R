qPCRHeatmap <- function(dat){
  require(lattice)
  dataMatrix = matrix(dat$Cq, 12, 8)
  mat2 = dataMatrix[,ncol(dataMatrix):1]
  png("heatmap.png")
  print(levelplot(mat2, col.regions=heat.colors(96, 1), main="Plate CF Values", scales=list(draw=F), xlab="", ylab=""))
  dev.off()
}
qPCRHeatmap(dat)