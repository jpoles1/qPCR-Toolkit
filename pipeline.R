sel = dat[dat$Fluor=="FAM",]
#genes=c("GAPDH","YM-1","ARG-1","Fizz-1")
#genes=c("UCP-1","PD-L2","Gata-6","MR1")
names=c("WT","WT","WT","WT","IL4-C","IL4-C","IRF-4 KO","IRF-4 KO","IRF-4 KO")
extract2 = qpcr(sel, 9,4,2,3,genes,savegraphs=F,treatment=names)
extract = cbind(extract1, extract2[,2:5])
#pdf("data.pdf", width=18,height=10)
#datplot=qplot(treatment, value,data=pdat,geom=c("boxplot","jitter"), facets=.~variable, fill=treatment)
#print(datplot)
#dev.off()c
pdat = cbind(extract$treatment, as.data.frame(2^(-1*abs(extract[,3:9]-extract[,2]))))
names(pdat)[1] = "Treatment"
pmelt = melt(pdat)
qplot(Treatment, value,data=pmelt,geom=c("point"), facets=.~variable, fill=Treatment)
w=pdat[pdat$Treatment=="WT",]
wtavg=colMeans(w[2:8], na.rm=T)
minwt = pdat[5:9,]
minwt = minwt/wtavg
plotdat=cbind(minwt[,1], paste(names[5:9], 1:5),log10(minwt[,2:8]/wtavg))
names(plotdat)[1] = "Treatment"
names(plotdat)[2] = "TreatID"
plotdatmelt= melt(plotdat)
ggplot(data=plotdatmelt, aes(x=TreatID, y=value, fill=Treatment))+geom_bar(stat="identity")+facet_grid(.~variable)+scale_x_discrete(breaks=NULL, name="")+scale_y_continuous(name="Log Fold Change")
