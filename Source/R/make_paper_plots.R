######################################################################
## script for making the plots from reactive vac and hotspots paper ##
## (has some extras)                                                ##
######################################################################
source("Source/R/plots_from_python_functions.R")

#####################################
## load all the runs in global env ##
#####################################
load.all.current.runs()

##########################################
## first all things with 2-patch models ##
##########################################

## Making the 2-pop subplots for the paper figure
pdf("Plots/subplot_uncon_nolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()
pdf("Plots/subplot_mild_nolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()
pdf("Plots/subplot_high_nolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

## now we will make the full plots for larger megaplot
pdf("Plots/pctbetter_alpha0p2.pdf",width=11,height=6)
make.pct.better.plot(high$fs,high$timings.save,
                     high$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,
                     plot.epi.curves = F,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p01.pdf",width=11,height=6)
make.pct.better.plot(mild$fs,mild$timings.save,
                     mild$null.run,pct.vac.plot.lim=1,
                     plot.vac.thresh = F,plot.epi.curves = F,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p0.pdf",width=11,height=6)
make.pct.better.plot(fs.con = uncon$fs,
                     timings.save = uncon$timings.save,
                     uncon$null.run,pct.vac.plot.lim=1,
                     plot.vac.thresh = FALSE,
                     plot.epi.curves = FALSE,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p2_curves.pdf",width=11,height=6)
make.pct.better.plot(high$fs,timings.save = high$timings.save,
                     high$null.run,pct.vac.plot.lim=1,
                     plot.vac.thresh = TRUE,
                     plot.epi.curves = TRUE,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p01_curves.pdf",width=11,height=6)
make.pct.better.plot(mild$fs,mild$timings.save,
                     mild$null.run,pct.vac.plot.lim=1,
                     plot.vac.thresh = TRUE,
                     plot.epi.curves = TRUE,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p0_curves.pdf",width=11,height=6)
make.pct.better.plot(uncon$fs,uncon$timings.save,uncon$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = TRUE,
                     plot.epi.curves = TRUE,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()


###########################################################################
## Now lets compare only targeting the hotspot to targeting the coldspot ##
###########################################################################
pdf("Plots/pctbetter_alpha0p0_rest_targeted.pdf",width=11,height=6)
make.pct.better.plot.restricted(uncon$fs,uncon$timings.save,uncon$null.run,
                                pvacs=seq(0,1,length=41),
                                ptiles=seq(0,1,length=41),
                                Rs=expand.grid(seq(1.5,2.5,length=5),
                                    seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()


pdf("Plots/pctbetter_alpha0p01_rest_targeted.pdf",width=11,height=6)
make.pct.better.plot.restricted(mild$fs,mild$timings.save,mild$null.run,
                                pvacs=seq(0,1,length=41),
                                ptiles=seq(0,1,length=41),
                                Rs=expand.grid(seq(1.5,2.5,length=5),
                                    seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p2_rest_targeted.pdf",width=11,height=6)
make.pct.better.plot.restricted(high$fs,high$timings.save,high$null.run,
                                pvacs=seq(0,1,length=41),
                                ptiles=seq(0,1,length=41),
                                Rs=expand.grid(seq(1.5,2.5,length=5),
                                    seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()


#####################################################
## Epidemic Plots with decision points highlighted ##
#####################################################

run.num <- 12
fs.con <- mild$fs
timings <- mild$timings.save[run.num,]
cis <- apply(mild$null.run[[run.num]][,9:10],2,diff)
cis <- cis[1:max(timings),]
stat.mat <- get.stat.mats(fs.con=mild$fs,
                          num.vac.strat = 3,
                          zlims=c(0,1),
                          run.num=run.num)


#quartz("",height=4,width=5)
pdf("Plots/realtimecurves_mild12.pdf",height=4,width=5)
par(mfrow=c(2,1),mar=c(0.1,4,0.1,1),oma=c(3,1,3,0),
    mgp=c(3,.5,0),tck=-.05,mgp=c(1.5,.5,0))
plot(cis[,1],type="l",axes=FALSE,col="black",ylab="Incident Cases \n (hotspot)")
abline(v=seq(0,300,by=20),col=AddAlpha("grey",.1))
axis(2,cex.axis=.85)
axis(3,labels=seq(0,1,length=41)[c(1,3,11,21,31,39,41)],at=round(timings)[c(1,3,11,21,31,39,41)],cex.axis=.6,las=2)
#axis(3,labels=seq(0,1,length=41),at=round(timings),cex.axis=.6,las=2)

## get the borders
red.border <- timings[range(which(!is.na(stat.mat$stat1[,run.num])))]
blue.border <- timings[range(which(!is.na(stat.mat$stat2[,run.num])))]
green.border <- timings[range(which(!is.na(stat.mat$stat3[,run.num])))]
red.green.midpoint <- (red.border[2]+green.border[1])/2
green.blue.midpoint <- (green.border[2] + blue.border[1])/2
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,green.blue.midpoint,green.blue.midpoint,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
polygon(x=c(green.blue.midpoint,1000,1000,
            green.blue.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Blues")[3],.2),
        border=FALSE)
text(par("usr")[2]*.8,par("usr")[4]*.9,"Hotspot")
plot(cis[,2],type="l",axes=FALSE,col="black",ylab="Incident Cases \n (non-hotspot)")
axis(2,cex.axis=.85)
abline(v=seq(0,300,by=20),col=AddAlpha("grey",.1))
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,green.blue.midpoint,green.blue.midpoint,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
polygon(x=c(green.blue.midpoint,1000,1000,
            green.blue.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Blues")[3],.2),
        border=FALSE)
text(par("usr")[2]*.8,par("usr")[4]*.9,"Non-hotspot")
axis(1,outer=F)
mtext("Epidemic Day",side=1,outer=T,line=1.5,at=.55)
mtext("Percent of Uncontrolled Global Epidemic Elapsed",side=3,outer=T,line=1.5,at=.55)
dev.off()

## for unconnected

run.num <- 12
fs.con <- uncon$fs
timings <- uncon$timings.save[run.num,]
cis <- apply(uncon$null.run[[run.num]][,9:10],2,diff)
cis <- cis[1:max(timings),]
stat.mat <- get.stat.mats(fs.con=uncon$fs,
                          num.vac.strat = 3,
                          zlims=c(0,1),
                          run.num=run.num)


#quartz("",height=4,width=5)
pdf("Plots/realtimecurves_uncon12.pdf",height=4,width=5)
par(mfrow=c(2,1),mar=c(0.1,4,0.1,1),oma=c(3,1,3,0),
    mgp=c(3,.5,0),tck=-.05,mgp=c(1.5,.5,0))
plot(cis[,1],type="l",axes =FALSE,col="black",ylab="Incident Cases \n (hotspot)")
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
axis(2,cex.axis=.85)
axis(3,labels=seq(0,1,length=41)[c(1,3,21,31,39,41)],
     at=round(timings)[c(1,3,21,31,39,41)],cex.axis=.6,las=2)
#axis(3,labels=seq(0,1,length=41),at=round(timings),cex.axis=.6,las=2)

## get the borders
red.border <- c(0,.01) # so that it shows up
blue.border <- timings[range(which(!is.na(stat.mat$stat2[,run.num])))]
green.border <- timings[range(which(!is.na(stat.mat$stat3[,run.num])))]
red.green.midpoint <- green.border[1]/2
blue.green.midpoint <- (green.border[2] + blue.border[1])/2
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,blue.green.midpoint,blue.green.midpoint,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
polygon(x=c(blue.green.midpoint,1000,1000,
            blue.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Blues")[3],.2),
        border=FALSE)
text(par("usr")[2]*.8,par("usr")[4]*.9,"Hotspot")
plot(cis[,2],type="l",axes =FALSE,col="black",ylab="Incident Cases \n (non-hotspot)")
axis(2,cex.axis=0.85)
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,blue.green.midpoint,blue.green.midpoint,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
polygon(x=c(blue.green.midpoint,1000,1000,
            blue.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Blues")[3],.2),
        border=FALSE)

text(par("usr")[2]*.8,par("usr")[4]*.9,"Non-hotspot")
axis(1,outer=F,cex.axis=.85)
mtext("Epidemic Days",side=1,outer=T,line=1.5,at=.55)
mtext("Percent of Uncontrolled Global Epidemic Elapsed",side=3,outer=T,line=1.5,at=.55)
dev.off()

## for highly connected
run.num <- 12
fs.con <- high$fs
timings <- high$timings.save[run.num,]
cis <- apply(high$null.run[[run.num]][,9:10],2,diff)
cis <- cis[1:max(timings),]
stat.mat <- get.stat.mats(fs.con=high$fs,
                          num.vac.strat = 3,
                          zlims=c(0,1),
                          run.num=run.num)


quartz("",height=4,width=5)
pdf("Plots/realtimecurves_high12.pdf",height=4,width=5)
par(mfrow=c(2,1),mar=c(0.1,4,0.1,1),oma=c(3,1,3,0),
    mgp=c(3,.5,0),tck=-.05,mgp=c(1.5,.5,0))
plot(cis[,1],type="l",axes=FALSE,col="black",ylab="Incident Cases \n (hotspot)")
axis(2,cex.axis=.85)
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
axis(3,labels=seq(0,1,length=41)[c(1,3,11,21,31,39,41)],at=round(timings)[c(1,3,11,21,31,39,41)],cex.axis=.6,las=2)
#axis(3,labels=seq(0,1,length=41),at=round(timings),cex.axis=.6,las=2)

## get the borders
red.border <- timings[range(which(!is.na(stat.mat$stat1[,run.num])))]
#blue.border <- timings[range(which(!is.na(stat.mat$stat2[,run.num])))]
green.border <- timings[range(which(!is.na(stat.mat$stat3[,run.num])))]

## interpolating between red and grreen so th transition is smooth
red.green.midpoint <- (red.border[2] + green.border[1])/2
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,1000,1000,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
text(par("usr")[2]*.8,par("usr")[4]*.9,"Hotspot")
plot(cis[,2],type="l",axes=FALSE,col="black",ylab="Incident Cases \n (non-hotspot)")
axis(2,cex.axis=.85)
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,1000,1000,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)

text(par("usr")[2]*.8,par("usr")[4]*.9,"Non-hotspot")
axis(1,outer=F)
mtext("Epidemic Days",side=1,outer=T,line=1.5,at=.55)
mtext("Percent of Uncontrolled Global Epidemic Elapsed",
      side=3,outer=T,line=1.5,at=.55)
dev.off()


## for 5 patch system
## for highly connected
run.num <- 12
fs.con.HCCCC<- high.HCCCC$fs
timings.HCCCC <- high.HCCCC$timings.save[run.num,]
cis.HCCCC <- apply(high.HCCCC$null.run[[run.num]][,21:25],2,diff)
cis.HCCCC <- cis.HCCCC[1:max(timings.HCCCC),]
stat.mat.HCCCC <- get.stat.mats(fs.con=high.HCCCC$fs,
                               num.vac.strat = 3,
                               zlims=c(0,1),
                               run.num=run.num)


run.num <- 12
fs.con <- high$fs
timings <- high$timings.save[run.num,]
cis <- apply(high$null.run[[run.num]][,9:10],2,diff)
cis <- cis[1:max(timings),]
stat.mat <- get.stat.mats(fs.con=high$fs,
                          num.vac.strat = 3,
                          zlims=c(0,1),
                          run.num=run.num)



quartz("",height=4,width=5)
pdf("Plots/realtimecurves_high12_5pop.pdf",height=4,width=5)
par(mfrow=c(2,1),mar=c(0.1,4,0.1,1),oma=c(3,1,3,0),
    mgp=c(3,.5,0),tck=-.05,mgp=c(1.5,.5,0))
plot(cis[,1],type="l",axes=FALSE,col=AddAlpha("grey",.5),ylab="Incident Cases \n (hotspot)")
lines(cis.HCCCC[,1],col="black")

axis(2,cex.axis=.85)
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
axis(3,labels=seq(0,1,length=41)[c(1,3,11,21,31,39,41)],at=round(timings.HCCCC)[c(1,3,11,21,31,39,41)],cex.axis=.6,las=2)
axis(3,labels=seq(0,1,length=41)[c(1,3,11,21,31,39,41)],
     at=round(timings)[c(1,3,11,21,31,39,41)],
     cex.axis=.6,las=2,col=AddAlpha("grey",.5))

#axis(3,labels=seq(0,1,length=41),at=round(timings.HCCCC),cex.axis=.6,las=2)

## get the borders
red.border <- timings.HCCCC[range(which(!is.na(stat.mat.HCCCC$stat1[,run.num])))]
green.border <- timings.HCCCC[range(which(!is.na(stat.mat.HCCCC$stat3[,run.num])))]

## interpolating between red and grreen so th transition is smooth
red.green.midpoint <- (red.border[2] + green.border[1])/2

## get it for the two patch system
red.border.2p <- timings[range(which(!is.na(stat.mat$stat1[,run.num])))]
green.border.2p <- timings[range(which(!is.na(stat.mat$stat3[,run.num])))]

## interpolating between red and grreen so th transition is smooth
red.green.midpoint.2p <- (red.border.2p[2] + green.border.2p[1])/2

abline(v=red.green.midpoint,col="black",lty=2,lwd=.5)
abline(v=red.green.midpoint.2p,col=AddAlpha("grey",.5),lty=2,lwd=.5)

polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,1000,1000,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)
text(par("usr")[2]*.8,par("usr")[4]*.9,"Hotspot")
plot(cis[,2],type="l",axes=FALSE,col=AddAlpha("grey",.5),ylab="Incident Cases \n (non-hotspot)")
lines(cis.HCCCC[,2],col="black")

axis(2,cex.axis=.85)
abline(v=seq(0,600,by=20),col=AddAlpha("grey",.1))
polygon(x=c(red.border[1],red.green.midpoint,red.green.midpoint,red.border[1]),
        y=c(0,0,1e5,1e5),col=AddAlpha(brewer.pal(3,"Reds")[3],.2),
        border=FALSE)
polygon(x=c(red.green.midpoint,1000,1000,
            red.green.midpoint),y=c(0,0,1e5,1e5),
        col=AddAlpha(brewer.pal(3,"Greens")[3],.2),border=FALSE)

text(par("usr")[2]*.8,par("usr")[4]*.9,"Non-hotspot")
axis(1,outer=F)

abline(v=red.green.midpoint,col="black",lty=2,lwd=.5)
abline(v=red.green.midpoint.2p,col=AddAlpha("grey",.5),lty=2,lwd=.5)

mtext("Epidemic Days",side=1,outer=T,line=1.5,at=.55)
mtext("Percent of Uncontrolled Global Epidemic Elapsed",
      side=3,outer=T,line=1.5,at=.55)
dev.off()


#########################################
## Asymetric connectivity explorations ##
## 2-pops                              ##
#########################################

## load asym runs
## first for Rhot = 1.5, Rnot=0.94
asym.Rh1p5.Rn0p94.a0p01 <- load.numpy.asym.plots(timestamp=138843811778,
                                                 alpha.string = "0.01",
                                                 pvacs=seq(0,1,length=41),
                                                  ptiles=seq(0,1,length=41))

asym.Rh1p5.Rn0p94.a0p2 <- load.numpy.asym.plots(timestamp=138844427882,
                                                alpha.string = "0.2",
                                                pvacs=seq(0,1,length=41),
                                                  ptiles=seq(0,1,length=41))


## then for Rhot = 2.5, Rnot=1.125
asym.Rh2p5.Rn1p125.a0p2 <- load.numpy.asym.plots(timestamp=138842619746,
                                                 alpha.string = "0.2",
                                                 pvacs=seq(0,1,length=41),
                                                 ptiles=seq(0,1,length=41))

asym.Rh2p5.Rn1p125.a0p01 <- load.numpy.asym.plots(timestamp=138842152458,
                                                  alpha.string = "0.01",
                                                  pvacs=seq(0,1,length=41),
                                                  ptiles=seq(0,1,length=41))

pdf("Plots/asymplot_Rh1p5.Rn0p94_a0p2.pdf")
make.threshold.contours(asym.Rh1p5.Rn0p94.a0p2)
dev.off()

pdf("Plots/asymplot_Rh1p5.Rn0p94_a0p01.pdf")
make.threshold.contours(asym.Rh1p5.Rn0p94.a0p01)
dev.off()

pdf("Plots/asymplot_Rh2p5_Rn1p125_a0p01.pdf")
make.threshold.contours(asym.Rh2p5.Rn1p125.a0p01)
dev.off()

pdf("Plots/asymplot_Rh2p5_Rn1p125_a0p2.pdf")
make.threshold.contours(asym.Rh2p5.Rn1p125.a0p2)
dev.off()

#############################
## Now 3 patch systems fun ##
#############################

## make subplots for paper
pdf("Plots/subplot_uncon_HCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon.HCC,labs=FALSE,
                      num.vac.strat = 4,
                      ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41),
                      )
dev.off()
pdf("Plots/subplot_mild_HCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild.HCC,
                      labs=FALSE,
                      num.vac.strat = 4,
                      ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))

dev.off()
pdf("Plots/subplot_high_HCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high.HCC,labs=FALSE,
                      ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41),
                      num.vac.strat = 4)
dev.off()

## ## make more plots for larger megaplot
## pdf("Plots/pctbetter_alpha0p2_HCH.pdf",width=11,height=6)
## make.pct.better.plot(high.HCH$fs,high.HCH$timings.save,
##                      high.HCH$null.run,pct.vac.plot.lim=1,
##                      plot.vac.thresh = F,plot.epi.curves = F,
##                      spotconfig="HCH",
##                      pvacs=seq(0,1,length=21),
##                      ptiles=seq(0,1,length=21),
##                      num.vac.strat = 4,
##                      plot.legend = FALSE)
## dev.off()
## pdf("Plots/pctbetter_alpha0p01_HCH.pdf",width=11,height=6)
## make.pct.better.plot(mild.HCH$fs,mild.HCH$timings.save,
##                      mild.HCH$null.run,pct.vac.plot.lim=1,
##                      plot.vac.thresh = F,plot.epi.curves = F,
##                      spotconfig="HCH",
##                      pvacs=seq(0,1,length=21),
##                      ptiles=seq(0,1,length=21))
## dev.off()
## pdf("Plots/pctbetter_alpha0p0_HCH.pdf",width=11,height=6)
## make.pct.better.plot(uncon.HCH$fs,uncon.HCH$timings.save,
##                      uncon.HCH$null.run,pct.vac.plot.lim=1,
##                      plot.vac.thresh = F,plot.epi.curves = F,
##                      spotconfig="HCH",
##                      pvacs=seq(0,1,length=21),
##                      ptiles=seq(0,1,length=21))
## dev.off()

pdf("Plots/pctbetter_alpha0p2_HCC.pdf",width=11,height=6)
make.pct.better.plot(high.HCC$fs,high.HCC$timings.save,
                     high.HCC$null.run,pct.vac.plot.lim=1,
                     plot.vac.thresh = F,plot.epi.curves = F,
                     spotconfig="HCC",
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     num.vac.strat=4,
                     plot.legend = FALSE)
dev.off()
pdf("Plots/pctbetter_alpha0p01_HCC.pdf",width=11,height=6)
make.pct.better.plot(mild.HCC$fs,mild.HCC$timings.save,
                     mild.HCC$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,plot.epi.curves = F,
                     spotconfig="HCC",
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     num.vac.strat=4,
                     plot.legend = FALSE)
dev.off()

pdf("Plots/pctbetter_alpha0p0_HCC.pdf",width=11,height=6)
make.pct.better.plot(uncon.HCC$fs,uncon.HCC$timings.save,
                     uncon.HCC$null.run,
                     pct.vac.plot.lim=1,plot.vac.thresh = F,
                     plot.epi.curves = F,
                     spotconfig="HCC",
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     num.vac.strat=4,
                     plot.legend = FALSE)
dev.off()


max.pct.dif <- numeric(5)
for (i in 21:25){
    stat.mats <- get.stat.mats(fs.con=high.HCC$fs,
                               num.vac.strat = 4,
                               zlims=c(0,1),
                               run.num=22)
    max.pct.dif[i] <- sum(stat.mats$stat1,na.rm=T)
}

##########################
## 4 patch simulations  ##
##########################

## make subplots for paper
pdf("Plots/subplot_uncon_HCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon.HCCC,labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()
pdf("Plots/subplot_mild_HCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild.HCCC,labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()
pdf("Plots/subplot_high_HCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high.HCCC,labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()


max.pct.dif <- numeric(5)
for (i in 21:25){
    stat.mats <- get.stat.mats(fs.con=high.HCCC$fs,
                               num.vac.strat = 4,
                               zlims=c(0,1),
                               run.num=22)
    max.pct.dif[i] <- sum(stat.mats$stat1,na.rm=T)
}

##############################
## 5 population simulations ##
##############################

## make subplots for paper
pdf("Plots/subplot_uncon_HCCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon.HCCCC,labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()
pdf("Plots/subplot_mild_HCCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild.HCCCC,
                      labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()
pdf("Plots/subplot_high_HCCCCnolab.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high.HCCCC,labs=FALSE,
                      pvacs=seq(0,1,length=41),num.vac.strat = 4)
dev.off()

max.pct.dif <- numeric(5)
for (i in 21:25){
    stat.mats <- get.stat.mats(fs.con=high.HCCCC$fs,
                               num.vac.strat = 4,
                               zlims=c(0,1),
                               run.num=22)
    max.pct.dif[i] <- sum(stat.mats$stat1,na.rm=T)
}



##################################
## Exploring 2 pop optimization ##
##################################

opt.ids <- c(138842068531,138842076887,138842078015)

## matrices represent fraction of vaccine optimally allocated in HS
opt.uncon <- load.numpy.opt(timestamp=opt.ids[1],alpha.string="0.0")
opt.mild <- load.numpy.opt(timestamp=opt.ids[2],alpha.string="0.01")
opt.high <- load.numpy.opt(timestamp=opt.ids[3],alpha.string="0.2")

pdf("Plots/opt_2pop_uncon.pdf",width=11,heigh=6)
plot.2pop.opt(opt.uncon)
dev.off()
pdf("Plots/opt_2pop_mild.pdf",width=11,heigh=6)
plot.2pop.opt(opt.mild)
dev.off()
pdf("Plots/opt_2pop_high.pdf",width=11,heigh=6)
plot.2pop.opt(opt.high)
dev.off()


####################################################################
## Create subplot version to match some of the other main figures ##
####################################################################

pdf("Plots/opt_2pop_uncon_sub.pdf",height=4,width=4.5)
plot.2pop.opt.subplots(opt.uncon)
dev.off()
pdf("Plots/opt_2pop_mild_sub.pdf",height=4,width=4.5)
plot.2pop.opt.subplots(opt.mild)
dev.off()
pdf("Plots/opt_2pop_high_sub.pdf",width=4.5,heigh=4)
plot.2pop.opt.subplots(opt.high)
dev.off()



###########################################
## Plots of proportion of cases averted  ##
###########################################
#HC.ids <- c(138657171722,138659127137,138661343024)
## load 2 pops data
uncon <- load.numpy.megaplot(timestamp=HC.ids[1],
                             alpha.string = "0.0",
                             ptiles=seq(0,1,length=41),
                             pvacs=seq(0,1,length=41),
                             Rnots=seq(0.75,1.5,length=5))
mild <- load.numpy.megaplot(timestamp=HC.ids[2],
                            alpha.string = "0.01",
                            ptiles=seq(0,1,length=41),
                            pvacs=seq(0,1,length=41),
                            Rnots=seq(0.75,1.5,length=5))
high <- load.numpy.megaplot(timestamp=HC.ids[3],
                            alpha.string = "0.2",
                            ptiles=seq(0,1,length=41),
                            pvacs=seq(0,1,length=41),
                            Rnots=seq(0.75,1.5,length=5))

Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5))
pdf("Plots/finalsize_reduction_2pop_uncon.pdf",width=11,height=6)
plot.best.strat.pct.red(uncon)
dev.off()

pdf("Plots/finalsize_reduction_2pop_mild.pdf",width=11,height=6)
plot.best.strat.pct.red(mild)
dev.off()

pdf("Plots/finalsize_reduction_2pop_high.pdf",width=11,height=6)
plot.best.strat.pct.red(high)
dev.off()


########################################
## Making histograms of FS statistics ##
########################################

pdf("Plots/fs_stat_hists_uncon.pdf",width=11,height=6)
make.pct.better.hists(uncon$fs,uncon$timings.save,uncon$null.run)
dev.off()
pdf("Plots/fs_stat_hists_mild.pdf",width=11,height=6)
make.pct.better.hists(mild$fs,mild$timings.save,mild$null.run)
dev.off()
pdf("Plots/fs_stat_hists_high.pdf",width=11,height=6)
make.pct.better.hists(high$fs,high$timings.save,high$null.run)
dev.off()

pdf("Plots/fs_stat_hists_uncon_comb.pdf",width=11,height=6)
make.pct.better.hists.comb(uncon$fs,uncon$timings.save,uncon$null.run)
dev.off()
pdf("Plots/fs_stat_hists_mild_comb.pdf",width=11,height=6)
make.pct.better.hists.comb(mild$fs,mild$timings.save,mild$null.run)
dev.off()
pdf("Plots/fs_stat_hists_high_comb.pdf",width=11,height=6)
make.pct.better.hists.comb(high$fs,high$timings.save,high$null.run)
dev.off()

(fs_share - fs_hotspot)/(fs_share)



#################################################
## Making plots for multiyear burnin epidemics ##
#################################################

## npy_alpha0p2_140072267922
high.mult <- load.numpy.megaplot(timestamp="140072267922",
                                 alpha.string ="0.2",
                                 pvacs=seq(0,1,length=41),
                                 ptiles=seq(0,1,length=41),
                                 Rnots=seq(0.75,1.5,length=5),
                                 Rhots=seq(3,5,length=5))

pdf("Plots/subplot_high_nolab_multiyr.pdf")
make.paper.pct.better(fs=high.mult,
                      labs=FALSE,
                      ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()


## npy_alpha0p2_140072267922
high.mult.short <- load.numpy.megaplot(timestamp="140072264927",
                                 alpha.string ="0.2",
                                 pvacs=seq(0,1,length=21),
                                 ptiles=seq(0,1,length=21),
                                 Rnots=seq(0.75,1.5,length=5),
                                 Rhots=seq(3,5,length=5))
make.paper.pct.better(fs=high.mult.short,
                      labs=FALSE,
                      ptiles=seq(0,1,length=21),
                      pvacs=seq(0,1,length=21))
pdf("Plots/pctbetter_alpha0p2-mult-short.pdf",width=11,height=6)
make.pct.better.plot(fs.con=high.mult.short$fs,
                     timings.save = high.mult.short$timings.save,
                     null.run = high.mult.short$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,
                     plot.epi.curves = F,
                     pvacs=seq(0,1,length=21),
                     ptiles=seq(0,1,length=21),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()

mild.mult.short <- load.numpy.megaplot(timestamp="140073448312",
                                 alpha.string ="0.01",
                                 pvacs=seq(0,1,length=21),
                                 ptiles=seq(0,1,length=21),
                                 Rnots=seq(0.75,1.5,length=5),
                                 Rhots=seq(3,5,length=5))

pdf("Plots/pctbetter_alpha0p01-mult-short.pdf",width=11,height=6)
make.pct.better.plot(fs.con=mild.mult.short$fs,
                     timings.save = mild.mult.short$timings.save,
                     null.run = mild.mult.short$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,
                     plot.epi.curves = F,
                     pvacs=seq(0,1,length=21),
                     ptiles=seq(0,1,length=21),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()


uncon.mult.short <- load.numpy.megaplot(timestamp="140074484841",
                                 alpha.string ="0.0",
                                 pvacs=seq(0,1,length=21),
                                 ptiles=seq(0,1,length=21),
                                 Rnots=seq(0.75,1.5,length=5),
                                 Rhots=seq(3,5,length=5))

pdf("Plots/pctbetter_alpha0p0-mult-short.pdf",width=11,height=6)
make.pct.better.plot(fs.con=uncon.mult.short$fs,
                     timings.save = uncon.mult.short$timings.save,
                     null.run = uncon.mult.short$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,
                     plot.epi.curves = F,
                     pvacs=seq(0,1,length=21),
                     ptiles=seq(0,1,length=21),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()




pdf("Plots/pctbetter_alpha0p2-mult.pdf",width=11,height=6)
make.pct.better.plot(fs.con=high.mult$fs,
                     timings.save = high.mult$timings.save,
                     null.run = high.mult$null.run,
                     pct.vac.plot.lim=1,
                     plot.vac.thresh = F,
                     plot.epi.curves = F,
                     pvacs=seq(0,1,length=41),
                     ptiles=seq(0,1,length=41),
                     Rs=expand.grid(seq(1.5,2.5,length=5),
                         seq(0.75,1.5,length=5)),plot.legend=FALSE)
dev.off()


## alternate vaccine efficiacy runs
HC.VE.ids <- c(140837800038,140842873206,140847513856)
    ## load 2 pops data
uncon.VE80 <- load.numpy.megaplot(timestamp=HC.VE.ids[1],
                                  alpha.string = "0.0",
                                  ptiles=seq(0,1,length=41),
                                  pvacs=seq(0,1,length=41),
                                  Rnots=seq(0.75,1.5,length=5),
                                  VE=80)
mild.VE80 <- load.numpy.megaplot(timestamp=HC.VE.ids[2],
                                 alpha.string = "0.01",
                                 ptiles=seq(0,1,length=41),
                                 pvacs=seq(0,1,length=41),
                                 Rnots=seq(0.75,1.5,length=5),
                                 VE=80)
high.VE80 <- load.numpy.megaplot(timestamp=HC.VE.ids[3],
                            alpha.string = "0.2",
                            ptiles=seq(0,1,length=41),
                            pvacs=seq(0,1,length=41),
                            Rnots=seq(0.75,1.5,length=5),
                            VE=80)
pdf("Plots/subplot_uncon_nolab_VE80.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon.VE80,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

pdf("Plots/subplot_mild_nolab_VE80.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild.VE80,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

pdf("Plots/subplot_high_nolab_VE80.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high.VE80,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

HC.VE60.ids <- c(140893344442,140896936626,140902162129)
uncon.VE60 <- load.numpy.megaplot(timestamp=HC.VE60.ids[1],
                                  alpha.string = "0.0",
                                  ptiles=seq(0,1,length=41),
                                  pvacs=seq(0,1,length=41),
                                  Rnots=seq(0.75,1.5,length=5),
                                  VE=60)
mild.VE60 <- load.numpy.megaplot(timestamp=HC.VE60.ids[2],
                                 alpha.string = "0.01",
                                 ptiles=seq(0,1,length=41),
                                 pvacs=seq(0,1,length=41),
                                 Rnots=seq(0.75,1.5,length=5),
                                 VE=60)
high.VE60 <- load.numpy.megaplot(timestamp=HC.VE60.ids[3],
                            alpha.string = "0.2",
                            ptiles=seq(0,1,length=41),
                            pvacs=seq(0,1,length=41),
                            Rnots=seq(0.75,1.5,length=5),
                            VE=60)
pdf("Plots/subplot_uncon_nolab_VE60.pdf",height=4,width=4.5)
make.paper.pct.better(fs=uncon.VE60,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

pdf("Plots/subplot_mild_nolab_VE60.pdf",height=4,width=4.5)
make.paper.pct.better(fs=mild.VE60,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()

pdf("Plots/subplot_high_nolab_VE60.pdf",height=4,width=4.5)
make.paper.pct.better(fs=high.VE60,labs=FALSE,ptiles=seq(0,1,length=41),
                      pvacs=seq(0,1,length=41))
dev.off()


################################################################
## Plots related to time scale comments from PROC B reviewers ##
################################################################

p10 <- 5
p25 <- 11
p50 <- 21
p75 <- 31
p95 <- 39
p100 <- 41

time.cols <- AddAlpha(brewer.pal(8,"Set2"),.8)
cbind(pvacs,mild$timings.save[1,])
obj <- mild
obj <- high
make.pctile.epi <- function(obj,runnum,xlim=c(0,250),ylim,add,
                            pvacs=seq(0,1,length=41),
                            ci.cols=9:10,lty=1,...){

    tmp.timings <- round(obj$timings.save[runnum,],0)
    tmp.null <- obj$null.run[[runnum]]
    cis <- rowSums(apply(tmp.null[,ci.cols,drop=F],2,diff))
    if(add){
        lines(tmp.timings[1]:tmp.timings[p10],cis[1:tmp.timings[p10]],
             col=time.cols[1],lwd=1,,lty=lty)
    } else {
        if(missing(ylim)) ylim <- range(cis)
        plot(tmp.timings[1]:tmp.timings[p10],cis[1:tmp.timings[p10]],
             col=time.cols[1],xlim=xlim,
             ylim=ylim,type="l",lwd=1,,lty=lty,...)
    }
    lines(tmp.timings[p10]:tmp.timings[p25],
          cis[tmp.timings[p10]:tmp.timings[p25]],col=time.cols[2],lwd=1,lty=lty)
    lines(tmp.timings[p25]:tmp.timings[p50],
          cis[tmp.timings[p25]:tmp.timings[p50]],col=time.cols[3],lwd=1,lty=lty)
    lines(tmp.timings[p50]:tmp.timings[p75],
          cis[tmp.timings[p50]:tmp.timings[p75]],col=time.cols[4],lwd=1,lty=lty)
    lines(tmp.timings[p75]:tmp.timings[p95],
          cis[tmp.timings[p75]:tmp.timings[p95]],col=time.cols[5],lwd=1,lty=lty)
    lines(tmp.timings[p95]:tmp.timings[p100],
          cis[tmp.timings[p95]:tmp.timings[p100]],col=time.cols[6],lwd=1,lty=lty)
}

runs <- c(7,9,12,14,22,24)
par(mfrow=c(2,3))
make.pctile.epi(high,runnum=runs[1],add=F,ylim=c(0,60000))
for (i in runs[-1]) make.pctile.epi(high,runnum=i,add=T)
make.pctile.epi(mild,runnum=runs[1],add=F,ylim=c(0,60000))
for (i in runs[-1]) make.pctile.epi(mild,runnum=i,add=T)
make.pctile.epi(uncon,runnum=runs[1],add=F,ylim=c(0,60000))
for (i in runs[-1]) make.pctile.epi(uncon,runnum=i,add=T)

make.pctile.epi(high,runnum=runs[1],add=F,ylim=c(0,40000),ci.cols=10)
for (i in runs[-1]) make.pctile.epi(high,runnum=i,add=T,ci.cols=10)
make.pctile.epi(mild,runnum=runs[1],add=F,ylim=c(0,40000),ci.cols=10)
for (i in runs[-1]) make.pctile.epi(mild,runnum=i,add=T,ci.cols=10)
make.pctile.epi(uncon,runnum=runs[1],add=F,ylim=c(0,40000),ci.cols=10)
for (i in runs[-1]) make.pctile.epi(uncon,runnum=i,add=T,ci.cols=10)



threepop.cols <- 13:15
make.pctile.epi(high.HCC,runnum=1,add=F,ylim=c(0,50000),xlim=c(0,350),ci.cols=threepop.cols)
for (i in 2:25) make.pctile.epi(high.HCC,runnum=i,add=T,ci.cols=threepop.cols)
make.pctile.epi(mild.HCC,runnum=1,add=F,ylim=c(0,50000),xlim=c(0,350),ci.cols=threepop.cols)
for (i in 2:25) make.pctile.epi(mild.HCC,runnum=i,add=T,ci.cols=threepop.cols)
make.pctile.epi(uncon.HCC,runnum=1,add=F,ylim=c(0,50000),xlim=c(0,350),ci.cols=threepop.cols)
for (i in 2:25) make.pctile.epi(uncon.HCC,runnum=i,add=T,ci.cols=threepop.cols)


pdf("Plots/pctepidemicelapsed_high.pdf",width=9,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(high,runnum=runs[1],add=F,ylim=c(0,30000),axes=F,ci.cols=9,xlim=c(0,150))
for (i in runs[-1]) make.pctile.epi(high,runnum=i,add=T,ci.cols=9)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(high,runnum=runs[1],add=F,ylim=rev(c(0,30000)),
                axes=F,ci.cols=10,xlim=c(0,150))
for (i in runs[-1]) make.pctile.epi(high,runnum=i,add=T,ci.cols=10,
                                    ylim=rev(c(0,40000)))
axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()


## try for mild connectivitiy now
pdf("Plots/pctepidemicelapsed_mild.pdf",width=9,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(mild,runnum=runs[1],add=F,ylim=c(0,40000),axes=F,
                ci.cols=9,xlim=c(0,150),lty=1)
for (i in runs[-1]) make.pctile.epi(mild,runnum=i,add=T,ci.cols=9,lty=1)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(mild,runnum=runs[1],add=F,ylim=rev(c(0,15000)),
                axes=F,ci.cols=10,xlim=c(0,150),lty=1)
for (i in runs[-1]) make.pctile.epi(mild,runnum=i,add=T,ci.cols=10,
                                    ylim=rev(c(0,15000)),lty=1)
axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()


pdf("Plots/pctepidemicelapsed_uncon.pdf",width=9,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(uncon,runnum=runs[1],add=F,ylim=c(0,40000),axes=F,
                ci.cols=9,xlim=c(0,500),lty=1)
for (i in runs[-1]) make.pctile.epi(uncon,runnum=i,add=T,ci.cols=9,lty=2)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(uncon,runnum=runs[1],add=F,ylim=rev(c(0,13000)),
                axes=F,ci.cols=10,xlim=c(0,400),lty=1)
for (i in runs[-1]) make.pctile.epi(uncon,runnum=i,add=T,ci.cols=10,
                                    ylim=rev(c(0,13000)),lty=1)
axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()


pdf("Plots/pctepidemicelapsed_highHCC.pdf",width=7,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(high.HCC,runnum=runs[1],add=F,ylim=c(0,30000),axes=F,ci.cols=13,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(high.HCC,runnum=i,add=T,ci.cols=13)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(high.HCC,runnum=runs[1],add=F,ylim=rev(c(0,20000)),
                axes=F,ci.cols=14,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(high.HCC,runnum=i,add=T,ci.cols=14)

axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()


pdf("Plots/pctepidemicelapsed_mildHCC.pdf",width=7,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(mild.HCC,runnum=runs[1],add=F,ylim=c(0,40000),axes=F,ci.cols=13,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(mild.HCC,runnum=i,add=T,ci.cols=13)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(mild.HCC,runnum=runs[1],add=F,ylim=rev(c(0,20000)),
                axes=F,ci.cols=14,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(mild.HCC,runnum=i,add=T,ci.cols=14)

axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()



pdf("Plots/pctepidemicelapsed_highHCCCC.pdf",width=7,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(high.HCCCC,runnum=runs[1],add=F,ylim=c(0,40000),axes=F,ci.cols=21,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(high.HCCCC,runnum=i,add=T,ci.cols=21)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(high.HCCCC,runnum=runs[1],add=F,ylim=rev(c(0,20000)),
                axes=F,ci.cols=22,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(high.HCCCC,runnum=i,add=T,ci.cols=22)

axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()

pdf("Plots/pctepidemicelapsed_mildHCCCC.pdf",width=7,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(mild.HCCCC,runnum=runs[1],add=F,ylim=c(0,40000),axes=F,ci.cols=21,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(mild.HCCCC,runnum=i,add=T,ci.cols=21)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(mild.HCCCC,runnum=runs[1],add=F,ylim=rev(c(0,20000)),
                axes=F,ci.cols=22,xlim=c(0,200))
for (i in runs[-1]) make.pctile.epi(mild.HCCCC,runnum=i,add=T,ci.cols=22)

axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()

pdf("Plots/pctepidemicelapsed_unconHCCCC.pdf",width=7,height=5)
par(mfrow = c(2,1),mar=c(0,0,0,0),oma=c(4,4,4,1))
make.pctile.epi(uncon.HCCCC,runnum=runs[1],add=F,ylim=c(0,5000),axes=F,ci.cols=21,xlim=c(0,300))
for (i in runs[-1])
make.pctile.epi(uncon.HCCCC,runnum=i,add=T,ci.cols=21)
grid()
axis(3)
axis(2)
legend("topright",c("0-10%","10-25%","25-50%","50-75%","75-95%","95-100%"),
       col=time.cols[1:6],lty=1,bty="n")
make.pctile.epi(uncon.HCCCC,runnum=runs[1],add=F,ylim=rev(c(0,20000)),
                axes=F,ci.cols=22,xlim=c(0,300))
for (i in runs[-1]) make.pctile.epi(uncon.HCCCC,runnum=i,add=T,ci.cols=22)

axis(1)
grid()
axis(2)
mtext("Non-hotspot Epidemic",side=2,outer=T,at=.25,line=2.2)
mtext("Hotspot Epidemic",side=2,outer=T,at=.75,line=2.2)
mtext("Patch Specific Epidemics by Percent of Global Epidemic Elapsed",side=3,outer=T,line=2.2)
mtext("Days",side=1,outer=T,line=2.2)
dev.off()


## obj <- mild
## runnum <- 14
## tmp.timings <- round(obj$timings.save[runnum,],0)
## tmp.null <- obj$null.run[[runnum]]2## cis <- rowSums(apply(tmp.null[,ci.cols,drop=F],2,diff))
## ecdf <- cumsum(cis)/max(cumsum(cis))
