## functions to make plots from pyton runs for vaccination
## hotspots paper

##' loads data from python runs and
##' makes it into a big list with all
##' useful elemets of the simulation
##' @param timestamp unique idenifier from simluation
##' @param alpha.string connectivity as a string
##' (used in getting the right filenames)
##' @param spotconfig if more than two pactches a combination of H's
##' and C's incicating the hotspot and coldspot locs
##' @param pvacs percent vaccinated sequence used in simulations
##' @param ptiles percentiles of uncontrolled epiedmic vaccinted
##' in simulations
##' @param Rnots vector of Rhots used in sims
##' @param Rhots vector of Rnots used in sims
##' @param VE 100*VE to be used for identifying run stamp
##' @return big ol' list
##' @author andrew azman
load.numpy.megaplot <- function(timestamp=138464063999,
                                alpha.string,
                                spotconfig=NULL,
                                pvacs=seq(0,1,length=21),
                                ptiles=seq(0,1,length=21),
                                Rnots=seq(0.9,1.5,length=5),
                                Rhots=seq(1.5,2.5,length=5),
                                VE=NULL
                                ){

    require(RcppCNPy)

    times <- seq(0,300)
    Rs <- expand.grid(Rhots,Rnots)
    VEadd <- ifelse(is.null(VE),"",paste0("VE",VE,"_"))
    ## construct the base file name string
    run_stamp <- paste0(spotconfig,
                        ifelse(is.null(spotconfig) ,"","_"),
                        "alpha",gsub("\\.","p",paste(alpha.string)),
                        "_",timestamp)
    dir_name <- paste0("Generated_Data/npy_",run_stamp,"/")

    ## load final size array
    fs <- array(dim=c(nrow(Rs),length(ptiles),length(pvacs),4))
    nullrun <- vector("list",nrow(Rs))
    for (i in 0:(nrow(Rs)-1)){
        nullrun[[i+1]] <-
            npyLoad(paste0(dir_name,
                           "nullrun_r",i,"_",VEadd,run_stamp,".npy"))
        for (j in 1:ifelse(is.null(spotconfig),3,4)){
            temp_mat <- npyLoad(paste0(dir_name,
                                       "fs_r",i,"_type",j,"_",VEadd,
                                       run_stamp,".npy"))
            fs[i+1,,,j] <- temp_mat
        }}

    timings_save <- npyLoad(paste0(dir_name,
                                   "timings_save_",VEadd,run_stamp,".npy"))

    # to fix the indexing discrepancy between R and python
    # now timings are in R indexing (starting at 1)
    timings_save <- timings_save+1

    return(list(fs=fs,null.run=nullrun,
                timings.save=timings_save,Rs=Rs))
}

##' Loads data from runs with asymetric connectivity
##' @param timestamp
##' @param alpha.string
##' @param pvacs
##' @param ptiles
##' @return
##' @author andrew azman
load.numpy.asym.plots <- function(timestamp=138612379942,
                                  alpha.string="0.01",
                                  pvacs=seq(0,1,length=21),
                                  ptiles=seq(0,1,length=21)
                                  ){

    require(RcppCNPy)
    times <- seq(0,300)
    Rnots <- 1.125
    Rhots <- 2.5
    Rs <- expand.grid(Rhots,Rnots)
    travel_hs_mults <- c(.3,.5,.9,1,1/0.9,1/.5,1/.3)
    ## construct the bas file name string
    run_stamp <- paste0("asym_alpha",
                        gsub("\\.","p",
                             paste(alpha.string)),"_",timestamp)
    dir_name <- paste0("Generated_Data/npy_",run_stamp,"/")

    ## load final size array
    fs <- array(dim=c(length(travel_hs_mults),
                    length(ptiles),length(pvacs),4))
    nullrun <- vector("list",length(travel_hs_mults))
    for (i in 0:(length(travel_hs_mults)-1)){
        nullrun[[i+1]] <-
            npyLoad(paste0(dir_name,
                           "nullrun_r",i,"_",run_stamp,".npy"))
        for (j in 1:3){
            temp_mat <- npyLoad(paste0(dir_name,
                                       "fs_r",i,
                                       "_type",j,
                                       "_",run_stamp,".npy"))
            fs[i+1,,,j] <- temp_mat
        }}

    timings_save <- npyLoad(paste0(dir_name,
                                   "timings_save_",run_stamp,".npy"))

    return(list(fs=fs,null.run=nullrun,
                timings.save=timings_save,Rs=Rs,
                travel_hs_mults=travel_hs_mults))
}

##' loads objects from python 2 pop optimization runs
##' @param timestamp
##' @param alpha.string
##' @param pvacs
##' @param ptiles
##' @param Rnots
##' @param Rhots
##' @return
##' @author andrew azman
load.numpy.opt <- function(timestamp=138705875492,
                           alpha.string="0.0",
                           pvacs=seq(0,1,length=21),
                           ptiles=seq(0,1,length=21),
                           Rnots=seq(0.75,1.5,length=5),
                           Rhots=seq(1.5,2.5,length=5)
                           ){

    require(RcppCNPy)

    times <- seq(0,300)
    Rs <- expand.grid(Rhots,Rnots)

    ## construct the bas file name string
    run_stamp <- paste0("alpha",
                        gsub("\\.","p",paste(alpha.string))
                        ,"_",timestamp)
    dir_name <- paste0("Generated_Data/npy_opt_",
                       run_stamp,"/")

    ## load opt hs allocation array
    fs <- array(dim=c(nrow(Rs),length(ptiles),length(pvacs)))
    nullrun <- vector("list",nrow(Rs))
    for (i in 0:(nrow(Rs)-1)){
        nullrun[[i+1]] <-
            npyLoad(paste0(dir_name,"nullrun_r",i,"_opt_",
                           run_stamp,".npy"))
            temp_mat <- npyLoad(paste0(dir_name,"opt_r",i,
                                       "_type1_opt_",
                                       run_stamp,".npy"))
            fs[i+1,,] <- temp_mat
        }

    timings_save <- npyLoad(paste0(dir_name,
                                   "timings_save_opt_",
                                   run_stamp,".npy"))

    return(list(fs=fs,null.run=nullrun,
                timings.save=timings_save,Rs=Rs))
}



##' Gets the matrix of percent betters
##' @param fs.con - list of output with all the goodies from vac runs
##' @param num.vac.strat - number of potential vaccination strategies
##' @param run.num - run number
##' @return returns list of matrices for percent betters
##' for each vaccination strategy
##' @author andrew azman
get.stat.mats <- function(fs.con,
                          num.vac.strat,
                          run.num,
                          zlims=c(0,1)){

    # 1 = hot-spot targeted
    # 2 = non hot-spot targeted
    # 3 = pro-rata
    # 4 = non-hotspot shared targeting (when multiple non-hotspots)

    fs.con <- round(fs.con)
    best.fs <- pmin(fs.con[run.num,,,1],fs.con[run.num,,,2],
                    fs.con[run.num,,,3],fs.con[run.num,,,4],
                    na.rm=T)
    worst.fs <- pmax(fs.con[run.num,,,1],fs.con[run.num,,,2],
                     fs.con[run.num,,,3],fs.con[run.num,,,4],
                     na.rm=T)

    ## get the indices so we know which strategy it comes from
    best.fs.index <- best.fs
    best.fs.index[] <- NA

    if (num.vac.strat == 4){
        best.fs.index[fs.con[run.num,,,1]<
                      fs.con[run.num,,,2] &
                      fs.con[run.num,,,1]<
                      fs.con[run.num,,,3] &
                      fs.con[run.num,,,1]<
                      fs.con[run.num,,,4]] <- 1
        best.fs.index[fs.con[run.num,,,2]<
                      fs.con[run.num,,,3] &
                      fs.con[run.num,,,2]<
                      fs.con[run.num,,,1] &
                      fs.con[run.num,,,2]<
                      fs.con[run.num,,,4]] <- 2
        best.fs.index[fs.con[run.num,,,3]<
                      fs.con[run.num,,,1] &
                      fs.con[run.num,,,3]<
                      fs.con[run.num,,,2] &
                      fs.con[run.num,,,3]<
                      fs.con[run.num,,,4]] <- 3
        best.fs.index[fs.con[run.num,,,4]<
                      fs.con[run.num,,,1] &
                      fs.con[run.num,,,4]<
                      fs.con[run.num,,,2] &
                      fs.con[run.num,,,4]<
                      fs.con[run.num,,,3]] <- 4
    } else {

            best.fs.index[fs.con[run.num,,,1]<
                          fs.con[run.num,,,2] &
                          fs.con[run.num,,,1]<
                          fs.con[run.num,,,3]] <- 1
            best.fs.index[fs.con[run.num,,,2]<
                          fs.con[run.num,,,3] &
                          fs.con[run.num,,,2]<
                          fs.con[run.num,,,1]] <- 2
            best.fs.index[fs.con[run.num,,,3]<
                          fs.con[run.num,,,1] &
                          fs.con[run.num,,,3]<
                          fs.con[run.num,,,2]] <- 3
        }

    my_stat <- (worst.fs - best.fs)/worst.fs
    my_stat[my_stat>zlims[2]] <- zlims[2]
    my_stat[my_stat<zlims[1]] <- zlims[1]

    stat1 <- stat2 <- stat3 <- stat4 <- my_stat
    stat1[best.fs.index != 1 | is.na(best.fs.index)] <- NA
    stat2[best.fs.index !=2 | is.na(best.fs.index)] <- NA
    stat3[best.fs.index !=3 | is.na(best.fs.index)] <- NA

    if (num.vac.strat >= 4){
        stat4[best.fs.index != 4 |
              is.na(best.fs.index)] <- NA
        return(list("stat1"=stat1,"stat2"=stat2,
                    "stat3"=stat3,"stat4"=stat4,
                    "stat"=my_stat))
    }

    return(list("stat1"=stat1,"stat2"=stat2,
                "stat3"=stat3,"stat"=my_stat))
}

##' Gets the matrix of percent betters
##' this version restircted to situations
##' where we are targeting a single population
##' assumed to be strategies 1 and 2
##' @param fs.con
##' @param run.num
##' @param zlims
##' @return returns list of matrices for percent betters
##' for each vaccination strategy
get.stat.mats.restricted <- function(fs.con,
                                     run.num,
                                     zlims=c(0,1)){

    fs.con <- round(fs.con)
    best.fs <- pmin(fs.con[run.num,,,1],fs.con[run.num,,,2],
                    na.rm=T)
    worst.fs <- pmax(fs.con[run.num,,,1],fs.con[run.num,,,2],
                     na.rm=T)

    ## get the indices so we know which strategy it comes from
    best.fs.index <- best.fs
    best.fs.index[] <- NA

    best.fs.index[fs.con[run.num,,,1]<
                  fs.con[run.num,,,2]] <- 1
    best.fs.index[fs.con[run.num,,,2]<
                  fs.con[run.num,,,1]] <- 2

    my_stat <- (worst.fs - best.fs)/worst.fs
    my_stat[my_stat>zlims[2]] <- zlims[2]
    my_stat[my_stat<zlims[1]] <- zlims[1]

    stat1 <- stat2 <- my_stat
    stat1[best.fs.index != 1 | is.na(best.fs.index)] <- NA
    stat2[best.fs.index !=2 | is.na(best.fs.index)] <- NA

    return(list("stat1"=stat1,"stat2"=stat2,
                "stat"=my_stat))
}

##' Gets the matrix of percent betters
##' @title
##' @param run output from numpy run
##' @param run.num index of simultion (what Rs)
##' @return returns matrix of percent of uncontrolled epidemic averted
##' for each vaccination strategy
##' @author andrew azman
get.pct.case.reduction.mat <- function(run,
                                       run.num){
    fs.con <- round(run$fs)
    null.fs <- fs.con[run.num,1,1,1]
    best.fs <- pmin(fs.con[run.num,,,1],fs.con[run.num,,,2],
                    fs.con[run.num,,,3],fs.con[run.num,,,4],
                    na.rm=T)
    return(best.fs/null.fs)
}


##' makes paneled plots of \Theta statistic
##' @param fs.con
##' @param timings.save
##' @param null.run
##' @param zlims limits fr the statistic
##' @param ptiles percentiles of the epidemic considered
##' @param plot.epi.curves
##' @param plot.vac.thresh TRUE if we want ablines for the critical vaccination tresholds in each location
##' @param pct.vac.plot.lim
##' @param pvacs
##' @param Rs
##' @param add.local.peaks
##' @param add.global.peak
##' @param spotconfig
##' @param plot.legend
##' @param num.vac.strat
##' @param stat a statistic of the worst and best final sizes that we want to display on our heatmap
##' @return
##' @author Andrew Azman
make.pct.better.plot <- function(fs.con=fs.high,
                                 timings.save,
                                 null.run=null.run.high,
                                 zlims = c(0,1),
                                 ptiles = seq(0,1,length=21),
                                 plot.epi.curves=TRUE,
                                 plot.vac.thresh=TRUE,
                                 pct.vac.plot.lim=1,
                                 pvacs=seq(0,1,length=21),
                                 Rs=expand.grid(
                                     seq(1.5,2.5,length=5),
                                     seq(0.75,1.5,length=5)),
                                 add.local.peaks=FALSE,
                                 add.global.peak=FALSE,
                                 spotconfig=NULL,
                                 plot.legend=TRUE,
                                 num.vac.strat=2
                                 ){
    require(fields)

    ## load colors
    cols1 <- colorRampPalette(brewer.pal(9,"Reds"))(100) # hotspot targeting
    cols2 <- colorRampPalette(brewer.pal(9,"Blues"))(100) # cold spot targeting
    cols3 <- colorRampPalette(brewer.pal(9,"Greens"))(100) # sharing
    if (!is.null(spotconfig) && spotconfig == "HCH"){
        cols4 <- colorRampPalette(brewer.pal(9,"Oranges"))(100) # sharing
    } else {
        cols4 <- colorRampPalette(brewer.pal(9,"Purples"))(100) # sharing
    }

    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,6,14),
        tck=0.02,mgp=c(3,0.1,0))

    max.times <- round(apply(timings.save,1,max))
    # round to integers based final sizes
    fs.con <- round(fs.con,0)

    for (i in 1:25){

        stat.mats <- get.stat.mats(fs.con=fs.con,
                                   num.vac.strat = num.vac.strat,
                                   zlims=zlims,
                                   run.num=i)

        yaxis.max <- which.min(abs(pvacs -pct.vac.plot.lim))


        image(stat.mats$stat1[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols1)
        image(stat.mats$stat2[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols2,
                  add=TRUE)
        image(stat.mats$stat3[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols3,
              add=TRUE)
        if (!is.null(spotconfig)){
            image(stat.mats$stat4[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols4,
                  add=TRUE)
        }

        if (i %in% seq(1,26,by=5)){
            axis(2,labels=
                 round(seq(0,pvacs[yaxis.max],length=5)*500000,1)
                 ,at=seq(0,1,length=5)
                 ,col = "grey",cex.axis=0.7)} else {
                     axis(2,labels=F,tick=T,at=seq(0,1,length=5),
                          col="grey")
                 }

        if (i > 20){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }

        if (plot.vac.thresh){
            ## add the critical vaccine threshold
            abline(h=1-(1/Rs[i,]),lty=3)
        }


        ## add the epi curves
        ## need to map them to the percentiles of the global epidemic
        if (is.null(spotconfig)){
            col.nums <- 9:10
        } else {
            col.nums <- 8:10
        }

        epi.curves <- apply(null.run[[i]][,col.nums],2,diff)
        times <- pmin(round(timings.save[i,])+1,nrow(epi.curves))

        if (add.local.peaks){
            peak.indices <- apply(epi.curves,2,which.max)
            peak.pctiles <- approx(y=ptiles,x=times,xout=peak.indices)
            abline(v=peak.pctiles$y,col="grey")
        }

        if (add.global.peak){
            peak.index <- which.max(rowSums(epi.curves))
            peak.pctiles <- approx(y=ptiles,x=times,xout=peak.index)
            abline(v=peak.pctiles$y,col="grey",lty=2)
        }

        if (plot.epi.curves){
            par(new=TRUE)
                                        #recover()
            plot(ptiles,epi.curves[times,1],
                 xlim=c(0,1),
                 axes=F,
                 bty="n",
                 type="l",
                 col="grey",
                 lty=2,
                 lwd=2)
            lines(ptiles,epi.curves[times,2],col="grey",lwd=2)
        }
    }

    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
          side=1,outer=TRUE,line=1.5)
    mtext("Vaccine Courses Available",side=2,
          outer=TRUE,line=1.5)

    if (plot.legend){
        par(oma=c( 0,0,0,1))# reset margin to be much smaller.
        set.panel() # reset plotting device

        image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                   col=cols1,smallplot = c(.95,.97,.35,.6))
        image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                   col=cols2,smallplot = c(.92,.94,.35,.6),
                       axis.args=list( at=0, labels=""))
        image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                   col=cols3,smallplot = c(.89,.91,.35,.6),
                       axis.args=list( at=0, labels=""))
        if (!is.null(spotconfig) && (spotconfig == "HCH" |
                                     spotconfig == "HCC")){
            image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                       col=cols4,smallplot = c(.86,.88,.35,.6),
                           axis.args=list( at=0, labels=""))    ## load colors
    cols1 <- colorRampPalette(brewer.pal(9,"Reds"))(100) # hotspot targeting
    cols2 <- colorRampPalette(brewer.pal(9,"Blues"))(100) # cold spot targeting

    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,4,4),
        tck=0.02,mgp=c(3,0.1,0))

    max.times <- round(apply(timings.save,1,max))
    # round to integers based final sizes
    fs.con <- round(fs.con,0)

    for (i in 1:25){

        if (i == 21){

            plot(-1000,-1000,xlim=c(200,300),ylim=c(200,300),
                 bty="n",axes=FALSE)
            text(222,247,bquote(.("hotspot") ~ Theta),adj=.5,
                 srt=90,cex=.8)
            image.plot(stat.mats$stat1,zlim=zlims,legend.only=TRUE,
                       bigplot=c(-0.1,0.1,.1,.9),legend.width = 3,
                       col=cols1,labs="",
                       axis.args = list(cex.axis = .5,mgp=c(0,.3,0)))
            image.plot(stat.mats$stat2,zlim=zlims,legend.only=TRUE,
                       bigplot=c(0.3,0.5,.1,.9),legend.width = 3,
                       col=cols2,
                       axis.args = list(cex.axis = .5,mgp=c(0,.3,0)))
            text(267,247,bquote(.("non-hotspot") ~ Theta),adj=.5,
                 srt=90,cex=.8)


        } else {
        stat.mats <- get.stat.mats.restricted(fs.con=fs.con,
                                              zlims=zlims,
                                              run.num=i)

        yaxis.max <- which.min(abs(pvacs -pct.vac.plot.lim))


        image(stat.mats$stat1[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols1)
        image(stat.mats$stat2[,1:yaxis.max],zlim=zlims,axes=FALSE,col=cols2,
                  add=TRUE)

        }

        if (i %in% seq(1,26,by=5) & i != 21){
            axis(2,labels=
                 round(seq(0,pvacs[yaxis.max],length=5)*500000,1)
                 ,at=seq(0,1,length=5)
                 ,col = "grey",cex.axis=0.7)} else {
                     axis(2,labels=F,tick=T,at=seq(0,1,length=5),
                          col="grey")
                 }

        if (i > 20 & i != 21){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }


    }

    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
          side=1,outer=TRUE,line=1.5)
    mtext("Vaccine Courses Available",side=2,
          outer=TRUE,line=1.5)

    if (plot.legend){
        par(oma=c( 0,0,0,1))# reset margin to be much smaller.
        set.panel() # reset plotting device
        image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                   col=cols1,smallplot = c(.95,.97,.35,.6))
        image.plot(stat.mats$stat,zlim=zlims,legend.only=TRUE,
                   col=cols2,smallplot = c(.92,.94,.35,.6),
                       axis.args=list( at=0, labels=""))
    }
}



##' @param zlims limits fr the statistic
##' @param ptiles percentiles of the epidemic considered
##' @param plot.epi.curves
##' @param plot.vac.thresh TRUE if we want ablines for the critical vaccination tresholds in each location
##' @param pct.vac.plot.lim
##' @param pvacs
##' @param Rs
##' @param seq)
##' @param add.local.peaks
##' @param add.global.peak
##' @param python
##' @param spotconfig
##' @param plot.indices
##' @return
##' @author Andrew Azman
make.paper.pct.better <- function(subplots=c(7,9,12,14,22,24),
                                  fs,
                                  middle.text="[Unconnected Epidemics]",
                                  labs=TRUE,
                                  zlims = c(0,1),
                                  ptiles = seq(0,1,length=21),
                                  pct.vac.plot.lim=1,
                                  pvacs=seq(0,1,length=21),
                                  Rs=expand.grid(seq(1.5,2.5,length=5),
                                      seq(0.75,1.5,length=5)),
                                  num.vac.strat=3
                                 ){
    require(fields)

    ## load colors
    cols1 <- colorRampPalette(brewer.pal(9,"Reds"))(100) # hotspot targeting
    cols2 <- colorRampPalette(brewer.pal(9,"Blues"))(100) # cold spot targeting
    cols3 <- colorRampPalette(brewer.pal(9,"Greens"))(100) # sharing
    cols4 <- colorRampPalette(brewer.pal(9,"Purples"))(100) # sharing

    layout(matrix(c(1,2,1,2,3,4,3,4,5,6,5,6),ncol=2,byrow=TRUE))
    par(oma=c(4,4,2,4),tck=0.02,mgp=c(3,0.5,0),mar=c(0.2,0.2,0.2,0.2))

    fs.con <- fs$fs
    timings.save <- fs$timings.save
    null.run <- fs$null.run
    max.times <- round(apply(timings.save,1,max))
    for (i in subplots){

        stat.mats <- get.stat.mats(fs.con=fs.con,
                                   num.vac.strat = num.vac.strat,
                                   zlims=zlims,
                                   run.num=i)

        yaxis.max <-which.min(abs(pvacs -pct.vac.plot.lim))

        image(stat.mats$stat1[,1:yaxis.max],zlim=zlims,axes=FALSE,
              col=cols1)

        s1 <- stat.mats$stat1[,1:yaxis.max]
        s1 <- ifelse(is.na(s1),0,1)
        if (all(!is.na(s1))){
            contour(s1,levels=0.5,drawlabels=FALSE,add=T,col="darkgrey")
        }

        image(
            stat.mats$stat2[,1:yaxis.max],zlim=zlims,axes=FALSE,
            col=cols2,
            add=TRUE)
        image(stat.mats$stat3[,1:yaxis.max],zlim=zlims,axes=FALSE,
              col=cols3,
              add=TRUE)

        s3 <- stat.mats$stat3[,1:yaxis.max]
        s3 <- ifelse(is.na(s3),0,1)
        contour(s3,levels=0.5,drawlabels=FALSE,add=T,col="darkgrey")

        if (num.vac.strat >= 4){
            image(stat.mats$stat4[,1:yaxis.max],zlim=zlims,
                  axes=FALSE,col=cols4,
                  add=TRUE)
        }

            if (i %in% c(7,12,22)){
                axis(2,labels=
                     c("",sprintf("%.2f",seq(0,pvacs[yaxis.max],length=5)/2)[-1]),
                     ##                     round(seq(0,pvacs[yaxis.max],length=5),1), #*500000
                     at=seq(0,1,length=5)
                     ,cex.axis=1.25,las=2)#0.75)
            } else {
                axis(2,labels=F,tick=T,at=seq(0,1,length=5))
            }

        if (i > 20){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),cex.axis=1.25)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2))
        }

    }

    if (labs){
        mtext(bquote(R[0] ~ Hotspot),side=3,line=-12.5,outer=TRUE,cex=.7)
        mtext(paste0(Rs[1:5,1][c(2,4)]),side=3,at=c(.25,.75),
              outer=TRUE,cex=.7,line=-13)
        mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,adj=.3,
              outer=TRUE,cex=0.7)
        mtext(paste0(rev(round(Rs[seq(1,25,by=5),2][c(1,3,5)],2))),side=4,
              at=c(0.1,.33,.55),outer=TRUE,cex=0.7,line=0)
        mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
              side=1,outer=TRUE,line=1.5,cex=0.9)
        mtext("Vaccine Available (% of population)",side=2,adj=0.1,
              outer=TRUE,line=1.5,cex=0.9)
    }
}

##' Makes contour plots for transitions from one vaccintion
##' stratefy to another for a single set of Rs
##' Assumes output is from python runs
##' @param fs.con
##' @param timings.save
##' @param null.run
##' @param stat
##' @param zlims
##' @param ptiles
##' @param plot.epi.curves
##' @param plot.vac.thresh
##' @param pct.vac.plot.lim
##' @param pvacs
##' @param Rs
##' @param add.local.peaks
##' @param add.global.peak
##' @param python
##' @param spotconfig
##' @return
##' @author andrew azman
make.threshold.contours <- function(dat,
                                    stat=(worst.fs - best.fs)/worst.fs,
                                    zlims = c(0,1),
                                    ptiles = seq(0,1,length=21),
                                    pct.vac.plot.lim=1,
                                    pvacs=seq(0,1,length=21)
                                    ){

    fs.con <- dat$fs
    timings.save <- dat$timings.save
    null.run <- dat$null.run
   # par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,6,14),tck=0.02,mgp=c(3,0.1,0))
    max.times <- round(apply(timings.save,1,max))


    # round to integers based final sizes
    fs.con <- round(fs.con,0)
    cols1 <- colorRampPalette(brewer.pal(7,"Reds"))(7) # hotspot targeting
    cols2 <- colorRampPalette(brewer.pal(7,"Blues"))(7) # hotspot targeting

    for (i in 1:length(dat$travel_hs_mults)){
        ## get the best and worst final size values for each
        best.fs <- pmin(fs.con[i,,,2],fs.con[i,,,1],fs.con[i,,,3])
        worst.fs <- pmax(fs.con[i,,,2],fs.con[i,,,1],fs.con[i,,,3])

        ## get the indices so we know which strategy it comes from
        best.fs.index <- best.fs
        best.fs.index[] <- NA

        best.fs.index[fs.con[i,,,1]<fs.con[i,,,2] &
                      fs.con[i,,,1]<fs.con[i,,,3]] <- 1
        best.fs.index[fs.con[i,,,2]<fs.con[i,,,3] &
                      fs.con[i,,,2]<fs.con[i,,,1]] <- 2
        best.fs.index[fs.con[i,,,3]<fs.con[i,,,1] &
                      fs.con[i,,,3]<fs.con[i,,,2]] <- 3

        stat[stat>zlims[2]] <- zlims[2]
        stat[stat<zlims[1]] <- zlims[1]

        stat1 <- stat2 <- stat3 <- stat
        stat1[best.fs.index != 1 | is.na(best.fs.index)] <- NA
        stat2[best.fs.index !=2 | is.na(best.fs.index)] <- NA
        stat3[best.fs.index !=3 | is.na(best.fs.index)] <- NA

        stat1.rev <- stat1
        stat1.rev[!is.na(stat1.rev)] <- 1
        stat1.rev[is.na(stat1.rev)] <- 0

        stat2.rev <- stat2
        stat2.rev[!is.na(stat2.rev)] <- 1
        stat2.rev[is.na(stat2.rev)] <- 0

        yaxis.max <- which.min(abs(pvacs -pct.vac.plot.lim))

        if (i == 1){
            contour(stat1.rev,nlevels=1,levels=1,
                    labels=round(dat$travel_hs_mults[i],1),
                    col=cols1[i],axes=FALSE,method="simple")
            axis(1)
            axis(2,at=seq(0,1,by=0.2),labels=seq(0,.5,by=.1))
        } else {
            contour(stat1.rev,nlevels=1,levels=1,
                    labels=round(dat$travel_hs_mults[i],1),
                    col=cols1[i],lty=ifelse(i==4,2,1),
                    lwd=ifelse(i==4,2,1),add=TRUE,method="simple")
        }

        contour(stat2.rev,nlevels=1,levels=1,lty=ifelse(i==4,2,1),
                lwd=ifelse(i==4,2,1),
                labels=round(dat$travel_hs_mults[i],1),
                col=cols2[i],add=TRUE,method="simple")

    }

}
##' Makes the plot for optimal hs allocation
##' @param opt
##' @return
##' @author andrew azman
plot.2pop.opt <- function(opt,pvacs=seq(0,1,length=21)){
    cols1 <- colorRampPalette(brewer.pal(9,"Reds"))(100) # hotspot targeting
    n.vac <- pvacs*500000
    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,6,4),
        tck=0.02,mgp=c(3,0.1,0))
    zlims <- c(0,1)#500000)
    for (i in 1:25){
        per.alloc.mat <- sapply(1:length(pvacs),function(x) opt$fs[i,,x] / n.vac[x])
        per.alloc.mat <- ifelse(is.nan(per.alloc.mat),0,per.alloc.mat)
        if (i == 21){
            plot(-100,-100,xlim=c(0,100),ylim=c(0,100),
                 bty="n",axes=FALSE)
            text(15,47,"Optimal \n Vaccine to Hotspot",adj=.5,srt=90,cex=.9)
            image.plot(per.alloc.mat,#opt$fs[i,,],
                       zlim=c(0,1),#500000),
                       legend.only=TRUE,
                       bigplot=c(0.1,0.3,.1,.9),legend.width = 3,
                       col=cols1)
        } else {
            image(per.alloc.mat,#opt$fs[i,,],
                  zlim=c(0,1),#500000),
                  axes=FALSE,col=cols1)
            contour(z=per.alloc.mat,
                    #round(opt$fs[i,,]),
                    levels=seq(0,1,#500000,
                        length=5),add=T,
                    labcex = 0.3,lwd=.7)
        }

    if (i %in% c(seq(1,26,by=5)[-5],22)){
        axis(2,labels=
             sprintf("%.3f",seq(0,1,length=5)),
             ##             round(seq(0,1,length=5)*500000,1),
             at=seq(0,1,length=5)
             ,col = "grey")#,cex.axis=0.6)
    } else {
        axis(2,labels=F,tick=T,at=seq(0,1,length=5),
             col="grey")
    }

        if (i > 20 && i != 21){
            axis(1,labels=sprintf("%.3f",seq(0,1,by=0.2)),
                 at=seq(0,1,by=0.2),col="grey")#,cex.axis=0.7)
        } else if (i != 21){
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }
    }
    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(opt$Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(opt$Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
          side=1,outer=TRUE,line=1.5)
    mtext("Vaccine Courses Available",side=2,
          outer=TRUE,line=1.5)
}

##' Makes the plot for optimal hs allocation
##' Subplots
##' @param opt
##' @return
##' @author andrew azman
plot.2pop.opt.subplots <- function(opt,
                                   pvacs=seq(0,1,length=21),
                                   subplots=c(7,9,12,14,22,24),
                                   labs=FALSE
                                   ){

    cols1 <- colorRampPalette(brewer.pal(9,"Reds"))(100) # hotspot targeting
    n.vac <- pvacs*500000

    layout(matrix(c(1,2,1,2,3,4,3,4,5,6,5,6),ncol=2,byrow=TRUE))
    par(oma=c(4,4,2,4),tck=0.02,mgp=c(3,0.1,0),mar=c(0.2,0.2,0.2,0.2))

    ## par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,6,4),
    ##     tck=0.02,mgp=c(3,0.1,0))
    zlims <- c(0,1)#500000)
    for (i in subplots){
        per.alloc.mat <- sapply(1:length(pvacs),function(x) opt$fs[i,,x] / n.vac[x])
        per.alloc.mat <- ifelse(is.nan(per.alloc.mat),0,per.alloc.mat)



        ## if (i == 21){
        ##     plot(-100,-100,xlim=c(0,100),ylim=c(0,100),
        ##          bty="n",axes=FALSE)
        ##     text(15,47,"Optimal \n Vaccine to Hotspot",adj=.5,srt=90,cex=.9)
        ##     image.plot(per.alloc.mat,#opt$fs[i,,],
        ##                zlim=c(0,1),#500000),
        ##                legend.only=TRUE,
        ##                bigplot=c(0.1,0.3,.1,.9),legend.width = 3,
        ##                col=cols1)
        ## } else {
        image(per.alloc.mat,#opt$fs[i,,],
              zlim=c(0,1),#500000),
              axes=FALSE,col=cols1)
        contour(z=per.alloc.mat,
                                        #round(opt$fs[i,,]),
                levels=seq(0,1,#500000,
                    length=5),add=T,
                labcex = 0.3,lwd=.7)

                    if (i %in% c(7,12,22)){
                        axis(2,labels=
                             round(seq(0,1,length=5),1)*500000
                             ,at=seq(0,1,length=5)
                             ,col = "grey",cex.axis=0.75)
                    } else {
                        axis(2,labels=F,tick=T,at=seq(0,1,length=5),col="grey")
                    }

        if (i > 20){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }


    }


    if (labs){
        mtext(bquote(R[0] ~ Hotspot),side=3,line=-12.5,outer=TRUE,cex=.7)
        mtext(paste0(Rs[1:5,1][c(2,4)]),side=3,at=c(.25,.75),
              outer=TRUE,cex=.7,line=-13)
        mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,adj=.3,
              outer=TRUE,cex=0.7)
        mtext(paste0(rev(round(Rs[seq(1,25,by=5),2][c(1,3,5)],2))),side=4,
              at=c(0.1,.33,.55),outer=TRUE,cex=0.7,line=0)
        mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
              side=1,outer=TRUE,line=1.5,cex=0.9)
        mtext("Vaccine Courses Available",side=2,adj=0.1,
              outer=TRUE,line=1.5,cex=0.9)
    }


}


##' Makes the plot for pct reduction in final size
##' compared to controlled epidemic for all simulations
##' @param opt
##' @return
##' @author andrew azman
plot.best.strat.pct.red <- function(run){
    cols1 <- colorRampPalette(brewer.pal(9,"Greys"))(100)
    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,4,4),
        tck=0.02,mgp=c(3,0.1,0))
    zlims <- c(0,1)
    for (i in 1:25){
        ## if (i == 21){
        ##     plot(-100,-100,xlim=c(0,100),ylim=c(0,100),
        ##          bty="n",axes=FALSE)
        ##     text(15,47,"Optimal \n Vaccine to Hotspot",adj=.5,srt=90)
        ##     image.plot(fs[i,,],zlim=c(0,500000),legend.only=TRUE,
        ##                bigplot=c(0.1,0.3,.1,.9),legend.width = 3,
        ##                col=cols1)
        ## } else {
        temp.dat <- 1- get.pct.case.reduction.mat(run,i)
        image(temp.dat,zlim=c(0,1),axes=FALSE,col=cols1)
        contour(z=temp.dat,levels=seq(0,1,length=11),add=T,
                labcex = 0.4,lwd=.4,labels=round(seq(0,1,length=11),1))
                                        #        }
        text(.7,.9,paste0("Final Size = ",
                          formatC(round(run$fs[i,1,1,1]),
                                  big.mark=',',format='d')),
             bg="white")
        if (i %in% seq(1,26,by=5)[-5]){
            axis(2,labels=
                 round(seq(0,1,length=5)*500000,1)
                 ,at=seq(0,1,length=5)
                 ,col = "grey",cex.axis=0.7)} else {
                     axis(2,labels=F,tick=T,at=seq(0,1,length=5),
                          col="grey")
                 }

        if (i > 20){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }
    }
    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Vaccination Timing (as % of final size of uncontrolled epidemic)",
          side=1,outer=TRUE,line=1.5)
    mtext("Vaccine Courses Available",side=2,
          outer=TRUE,line=1.5)
}

make.pct.better.hists <- function(fs.con=fs.high,
                                  timings.save,
                                  null.run=null.run.high,
                                  zlims = c(0,1),
                                  ptiles = seq(0,1,length=41),
                                  pct.vac.plot.lim=1,
                                  pvacs=seq(0,1,length=41),
                                  Rs=expand.grid(
                                      seq(1.5,2.5,length=5),
                                      seq(0.75,1.5,length=5)),
                                  spotconfig=NULL,
                                  num.vac.strat=3
                                  ){
    require(fields)

    ## load colors
    col1 <- AddAlpha(brewer.pal(3,"Reds")[3],.4) # hotspot targeting
    col2 <- AddAlpha(brewer.pal(3,"Blues")[3],.4) # cold spot targeting
    col3 <- AddAlpha(brewer.pal(3,"Greens")[3],.4) # sharing
    col4 <- AddAlpha(brewer.pal(3,"Purples")[1],.4) # sharing


    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,4,4),
        tck=0.02,mgp=c(3,0.1,0))

    max.times <- round(apply(timings.save,1,max))
                                        # round to integers based final sizes
    fs.con <- round(fs.con,0)

    for (i in 1:25){
        stat.mats <- get.stat.mats(fs.con=fs.con,
                                   num.vac.strat = num.vac.strat,
                                   zlims=zlims,
                                   run.num=i)
        xlim <- c(0,1)
        ymax <- 0
        already.plotted <- FALSE # flag for plotting

        if (i == 21){
            plot(-100,-100,xlim=c(0,100),ylim=c(0,100),
                 bty="n",axes=FALSE)
        } else {

            ## get densities
            d1.valid <- sum(!is.na(stat.mats$stat1)) > 1
            d2.valid <- sum(!is.na(stat.mats$stat2)) > 1
            d3.valid <- sum(!is.na(stat.mats$stat3)) > 1

            if (d1.valid){
                hist1 <- hist(stat.mats$stat1,breaks=10,plot=F)
                #dens1 <- density(stat.mats$stat1,na.rm=T)
                ymax <- max(c(hist1$counts,ymax),na.rm=T)
            }

            if (d2.valid){
                hist2 <- hist(stat.mats$stat2,breaks=10,plot=F)
                #dens2 <- density(stat.mats$stat2,na.rm=T)
                ymax <- max(c(hist2$counts,ymax),na.rm=T)
            }

            if (d3.valid){
                hist3 <- hist(stat.mats$stat3,breaks=10,plot=F)
#                dens3 <- density(stat.mats$stat3,na.rm=T)
                ymax <- max(c(hist3$counts,ymax),na.rm=T)
            }

            ylim <- c(0,ymax)
            print(ylim)

            border <- "white"

            if (d1.valid){
                hist(stat.mats$stat1,xlim=xlim,border=border,
                     col=col1,axes=F,xlab="",ylab="",main="",
                     breaks=10,ylim=ylim)
                #plot(dens1, xlim = xlim, ylim = ylim, axes=FALSE,
                #     main = '',col="white",xlab="",ylab="",
                #     panel.first = grid(col=AddAlpha("grey",.2)),lwd=.5)
                #polygon(dens1, density = -1, col = col1,border=FALSE)
                already.plotted <- TRUE
            }

            if (d2.valid){
                hist(stat.mats$stat2,xlim=xlim,border=border,
                     col=col2,axes=F,xlab="",ylab="",main="",
                     add=ifelse(already.plotted,TRUE,FALSE),
                     ylim=ylim,breaks=10)
                if (!already.plotted) already.plotted <- TRUE
                    ## plot(dens2, xlim = xlim, ylim = ylim, axes=FALSE,
                    ##      main = '',col="white",xlab="",ylab="",
                    ##      panel.first = grid(col=AddAlpha("grey",.2)),lwd=.5)
                }

            #polygon(dens2, density = -1, col = col2,border=FALSE)

            if (d3.valid){
                hist(stat.mats$stat3,xlim=xlim,border=border,
                     col=col3,axes=F,xlab="",ylab="",main="",ylim=ylim,
                     add=ifelse(already.plotted,TRUE,FALSE),breaks=10)
                if (!already.plotted) already.plotted <- TRUE

            }
            ## if (!is.null(spotconfig)){
            ##     dens4 <- density(stat.mats$stat4,na.rm=T)
            ##     if (!already.plotted){
            ##         plot(dens3, xlim = xlim, ylim = ylim, axes=FALSE,
            ##              main = '',col="white",xlab="",ylab="",
            ##              panel.first = grid(),lwd=.5)
            ##         already.plotted <- TRUE
            ##     }
            ##     polygon(dens4, density = -1, col = col4,border=FALSE)

            ## }

        }
#        if (i %in% seq(1,26,by=5)){
#            axis(2,col = "grey",cex.axis=0.7)
#        } else {

        if (i == 5){
            legend('topright',c('hotspot',
                               'non-hotspot','pro-rata'),
                   fill = c(col1, col2,col3), bty = 'n',
                   border = NA,cex=.75)
        }
        if (i != 21)
            axis(2,labels=F,tick=T,col="grey")
#        }

        if ((i > 20 & i != 21) | i == 16){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }

    }

    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Relative Performace of Best Vaccination Strategy to the Worst",
          side=1,outer=TRUE,line=1.5)
    mtext("Density",side=2,
          outer=TRUE,line=1.5)

}


##'
##' @param fs.con
##' @param timings.save
##' @param null.run
##' @param zlims
##' @param ptiles
##' @param pct.vac.plot.lim
##' @param pvacs
##' @param Rs
##' @param seq)
##' @param spotconfig
##' @param num.vac.strat
##' @return
##' @author andrew azman
make.pct.better.hists.comb <- function(fs.con=fs.high,
                                      timings.save,
                                      null.run=null.run.high,
                                      zlims = c(0,1),
                                      ptiles = seq(0,1,length=41),
                                      pct.vac.plot.lim=1,
                                      pvacs=seq(0,1,length=41),
                                      Rs=expand.grid(
                                          seq(1.5,2.5,length=5),
                                          seq(0.75,1.5,length=5)),
                                      spotconfig=NULL,
                                      num.vac.strat=3
                                      ){
    require(fields)

    ## load colors
    col <- "grey"

    par(mfrow=c(5,5),mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,4,4),
        tck=0.02,mgp=c(3,0.1,0))

    max.times <- round(apply(timings.save,1,max))
    ## round to integers based final sizes
    fs.con <- round(fs.con,0)

    for (i in 1:25){
        stat.mats <- get.stat.mats(fs.con=fs.con,
                                   num.vac.strat = num.vac.strat,
                                   zlims=zlims,
                                   run.num=i)
        xlim <- c(0,1)
        ymax <- 0
        already.plotted <- FALSE # flag for plotting

        if (i == 21){
            plot(-100,-100,xlim=c(0,100),ylim=c(0,100),
                 bty="n",axes=FALSE)
        } else {


            border <- "white"

            hist(stat.mats$stat,xlim=xlim,border=border,
                 col=col,axes=F,xlab="",ylab="",main=""
                 )
        }

        if (i != 21)
            axis(2,labels=F,tick=T,col="grey")
#        }

        if ((i > 20 & i != 21) | i == 16){
            axis(1,labels=seq(0,1,by=0.2),
                 at=seq(0,1,by=0.2),col="grey",cex.axis=0.7)
        } else {
            axis(1,labels=F,tick=T,at=seq(0,1,by=0.2),col="grey")
        }

    }

    mtext(bquote(R[0] ~ Hotspot),side=3,line=1,outer=TRUE)
    mtext(paste0(Rs[1:5,1]),side=3,at=seq(0.1,.90,length=5),
          outer=TRUE,cex=.9,line=0)
    mtext(bquote(R[0] ~ Non-hotspot),side=4,line=2,
          outer=TRUE)
    mtext(paste0(rev(round(Rs[seq(1,25,by=5),2],2))),side=4,
          at=seq(0.1,.9,length=5),outer=TRUE,cex=0.9,line=0)
    mtext("Relative Performace of Best Vaccination Strategy to the Worst",
          side=1,outer=TRUE,line=1.5)
    mtext("Density",side=2,
          outer=TRUE,line=1.5)

}


##' Loads all main current runs in the global env
##' all based on ids which are in the .org file in GeneratedData/
##' @return
##' @author andrew azman
load.all.current.runs <- function(){
    HC.ids <- c(138841321629,138841339645,138841336629)
    ## load 2 pops data
    uncon <<- load.numpy.megaplot(timestamp=HC.ids[1],
                                 alpha.string = "0.0",
                                 ptiles=seq(0,1,length=41),
                                 pvacs=seq(0,1,length=41),
                                 Rnots=seq(0.75,1.5,length=5))
    mild <<- load.numpy.megaplot(timestamp=HC.ids[2],
                                alpha.string = "0.01",
                                ptiles=seq(0,1,length=41),
                                pvacs=seq(0,1,length=41),
                                Rnots=seq(0.75,1.5,length=5))
    high <<- load.numpy.megaplot(timestamp=HC.ids[3],
                                alpha.string = "0.2",
                                ptiles=seq(0,1,length=41),
                                pvacs=seq(0,1,length=41),
                                Rnots=seq(0.75,1.5,length=5))


    hcc.ids <- c(138841351653,138841353163,138841354645)

    ## now load the 3 population sims
    high.HCC <<- load.numpy.megaplot(timestamp=hcc.ids[3],
                                    alpha.string ="0.2",
                                    spotconfig="HCC",
                                    pvacs=seq(0,1,length=41),
                                    ptiles=seq(0,1,length=41))
    mild.HCC <<- load.numpy.megaplot(timestamp=hcc.ids[2],
                                    alpha.string ="0.01",
                                    spotconfig="HCC",
                                    pvacs=seq(0,1,length=41),
                                    ptiles=seq(0,1,length=41))
    uncon.HCC <<- load.numpy.megaplot(timestamp=hcc.ids[1],
                                     alpha.string ="0.0",
                                     spotconfig="HCC",
                                     pvacs=seq(0,1,length=41),
                                     ptiles=seq(0,1,length=41))

    # 4 patch models
    hccc.ids <- c(138928741228,138932831936,138937450263)

    ## load data
    high.HCCC <<-
        load.numpy.megaplot(timestamp=hccc.ids[3],
                            alpha.string ="0.2",
                            spotconfig="HCCC",
                            pvacs=seq(0,1,length=41),
                            ptiles=seq(0,1,length=41))

    mild.HCCC <<-
        load.numpy.megaplot(timestamp=hccc.ids[2],
                            alpha.string ="0.01",
                            spotconfig="HCCC",
                            pvacs=seq(0,1,length=41),
                            ptiles=seq(0,1,length=41))

    uncon.HCCC <<- load.numpy.megaplot(timestamp=hccc.ids[1],
                                      alpha.string ="0.0",
                                      spotconfig="HCCC",
                                      pvacs=seq(0,1,length=41),
                                      ptiles=seq(0,1,length=41))


    hcccc.ids <- c(138928741236,138933585811,13893903802)

    high.HCCCC <<- load.numpy.megaplot(timestamp=hcccc.ids[3],
                                      alpha.string ="0.2",
                                      spotconfig="HCCCC",
                                      pvacs=seq(0,1,length=41),
                                      ptiles=seq(0,1,length=41))

    mild.HCCCC <<- load.numpy.megaplot(timestamp=hcccc.ids[2],
                                      alpha.string ="0.01",
                                      spotconfig="HCCCC",
                                      pvacs=seq(0,1,length=41),
                                      ptiles=seq(0,1,length=41))

    uncon.HCCCC <<- load.numpy.megaplot(timestamp=hcccc.ids[1],
                                       alpha.string ="0.0",
                                       spotconfig="HCCCC",
                                       pvacs=seq(0,1,length=41),
                                       ptiles=seq(0,1,length=41))
}

## get maximum time (both in pct, and real time and in realtion to
## the global and local peaks and vaccination level where a specific
##  strategy dominates
get.max.time.and.vac <- function(run.obj,run.num,vac.strat){

    ## find the peak
    pvacs <- seq(0,499999/500000,length=41)
    real.times <- run.obj$timings.save
    peak.time <- which.max(
        diff(run.obj[[2]][[run.num]][,10])
        #diff(rowSums(run.obj[[2]][[run.num]][,9:10]))
        )


    stat.mats <- get.stat.mats(fs.con=run.obj$fs,
                               num.vac.strat = 3,
                               zlims=c(0,1),
                               run.num=run.num)

    ## assuming this will only be used for situations with three vaccination
    ## scenarios. If not will need to change.
    maxes <- sapply(vac.strat,function(y)
                    apply(which(!is.na(stat.mats[[y]]),arr.ind = T),2,max))
    mins <- sapply(vac.strat,function(y)
                    apply(which(!is.na(stat.mats[[y]]),arr.ind = T),2,min))

    max.time = real.times[run.num,maxes[1,1]]
    min.time = real.times[run.num,mins[1,1]]

    pctiles <- seq(0,1,length=41)
    max.time.pct <- pctiles[maxes[1,1]]
    min.time.pct <- pctiles[mins[1,1]]

    max.pct.vac <- pvacs[maxes[2,1]]
    min.pct.vac <- pvacs[mins[2,1]]

    list(max.time = max.time, min.time = min.time,
         max.pct.vac = max.pct.vac, min.pct.vac = min.pct.vac,
         max.time.pct = max.time.pct,
         min.time.pct = min.time.pct,
         peak.time=peak.time
         )
}
