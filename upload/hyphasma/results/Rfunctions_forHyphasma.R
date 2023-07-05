## get single coloumn #colnumber from file 'filename' from multiple folders (000...)
##do mean and sd
##bind the info of the simulation (name of param file)
##if (do_perc) --> *100
get_results <- function(filename,colnumber,headerTF,deriv,do_perc) {
    final_set<-NULL ##to return
    ##fp set in the main file; vector containing all param files to analyse
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE) ##000...
        #test
        j<-1
        set_per_run<-NULL ##for each simul
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            
            if (deriv>0) {
                ##assumes an accumulation of 'deriv' hours
                time<-xx[which(xx$V1%%deriv<1e-15),1]
                xx<-xx[which(xx$V1%%deriv<1e-15),colnumber]
                xx<-c(xx[1],diff(xx))
            } else {
                time<-xx[which(xx$V1%%1<1e-15),1]
                xx<-xx[which(xx$V1%%1<1e-15),colnumber] ##get the coloumn and sort by hour
            }
            if(do_perc) {xx<-xx*100}
#             if (deriv) {xx<-c(xx[1],diff(xx))}
            set_per_run <- cbind(set_per_run,xx) ##bind 000...
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_')) ##split by the param characteristic --> note that this is not yet general        
        
        mean<-rowMeans(set_per_run)
        sds<-rowSds(as.matrix(set_per_run))
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$time<-time
        set_per_run$mean<-mean
        set_per_run$std<-sds
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        tmp<-set_per_run[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}

##returns: (cn1+cn2)/(cn3+cn4); cn2 and cn4 can be set to 0
## if cn3=0 --> cn3=1
## 'what' can be 'percentage' or 'fraction'
get_fraction_results <- function(filename,cn1,cn2,cn3,cn4,headerTF,deriv,what) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            xx<-xx[which(xx$V1%%1<1e-15),]
            time<-xx[which(xx$V1%%1<1e-15),1]
            
            if (max(xx$V1)>=1000) {
                xx<-xx[unique(xx$V1),]
                print('THERE IS A MISTAKE IN THE OUTPUT FILE !!!!!!!!')
            }
            xx1<-xx[,cn1]
            xx2<-0
            xx3<-1
            xx4<-0
            if (cn4>0) {xx4<-xx[,cn4]}
            if (cn2>0) {xx2<-xx[,cn2]}
            if (cn3>0) {xx3<-xx[,cn3]}
            
            xx<-(xx1+xx2)/(xx3+xx4)
            if (what=='percentage') {xx<-xx*100}
            if (deriv) {xx<-c(xx[1],diff(xx))}
            
            set_per_run <- cbind(set_per_run,xx)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        mean<-rowMeans(set_per_run)
        sds<-rowSds(as.matrix(set_per_run))
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$time<-time
        set_per_run$mean<-mean
        set_per_run$std<-sds
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        tmp<-set_per_run[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}

##returns: (cn1+cn2)/(cn3+cn4+cn5); cn2 and cn4 can be set to 0
## if cn3=0 --> cn3=1
## 'what' can be 'percentage' or 'fraction'
get_fraction2_results <- function(filename,cn1,cn2,cn3,cn4,cn5,headerTF,deriv,what) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            xx<-xx[which(xx$V1%%1<1e-15),]
            time<-xx[which(xx$V1%%1<1e-15),1]
            
            if (max(xx$V1)>=1000) {
                xx<-xx[unique(xx$V1),]
                print('THERE IS A MISTAKE IN THE OUTPUT FILE !!!!!!!!')
            }
            xx1<-xx[,cn1]
            xx2<-0
            xx3<-1
            xx4<-0
            xx5<-0
            if (cn4>0) {xx4<-xx[,cn4]}
            if (cn2>0) {xx2<-xx[,cn2]}
            if (cn3>0) {xx3<-xx[,cn3]}
            if (cn5>0) {xx5<-xx[,cn5]}
            
            xx<-(xx1+xx2)/(xx3+xx4+xx5)
            if (what=='percentage') {xx<-xx*100}
            if (deriv) {xx<-c(xx[1],diff(xx))}
            
            set_per_run <- cbind(set_per_run,xx)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        mean<-rowMeans(set_per_run)
        sds<-rowSds(as.matrix(set_per_run))
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$time<-time
        set_per_run$mean<-mean
        set_per_run$std<-sds
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        tmp<-set_per_run[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}


##returns: (cn1-cn2)
get_diff_results <- function(filename,cn1,cn2,headerTF,deriv) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            xx<-xx[which(xx$V1%%1<1e-15),]
            time<-xx[which(xx$V1%%1<1e-15),1]
            
            if (max(xx$V1)>=1000) {
                xx<-xx[unique(xx$V1),]
                print('THERE IS A MISTAKE IN THE OUTPUT FILE !!!!!!!!')
            }
            xx1<-xx[,cn1]
            xx2<-xx[,cn2]
            
            xx<-xx1-xx2
            
            if (deriv) {xx<-c(xx[1],diff(xx))}
            
            set_per_run <- cbind(set_per_run,xx)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        mean<-rowMeans(set_per_run)
        sds<-rowSds(as.matrix(set_per_run))
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$time<-time
        set_per_run$mean<-mean
        set_per_run$std<-sds
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        tmp<-set_per_run[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}


##returns the coloumns in 'cn' with names 'cnnames' r-binding all simulations
get_ALLDiff_results <- function(filename,cn1,cn2,headerTF,cnnames) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            
            xx1<-xx[,cn1]
            xx2<-xx[,cn2]
            
            xxA<-xx1-xx2
            
            xx<-cbind(xx[which(xx$V1%%1<1e-15),1],xxA[which(xx$V1%%1<1e-15)])
            
            colnames(xx) <- c('time',cnnames)
            
            yy<-data.frame(xx)
            
            yy$SimN<-j ##simulation number (0...)
            
            set_per_run <- rbind(set_per_run,yy)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        final_set<-rbind(final_set,set_per_run)
        setwd(results_folder)
    }
    return(final_set)
}




##same as before but storing values for each simulation separately
get_ALL_frac_results <- function(filename,cn1,cn2,cn3,cn4,headerTF,deriv,what) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            
            xx<-xx[which(xx$V1%%1<1e-15),]
            time<-xx[which(xx$V1%%1<1e-15),1]
            
            xx1<-xx[,cn1]
            xx2<-0
            xx3<-1
            xx4<-0
            if (cn4>0) {xx4<-xx[,cn4]}
            if (cn2>0) {xx2<-xx[,cn2]}
            if (cn3>0) {xx3<-xx[,cn3]}
            
            xx<-(xx1+xx2)/(xx3+xx4)
            if (what=='percentage') {xx<-xx*100}
            if (deriv) {xx<-c(xx[1],diff(xx))}
            
            
            yy<-data.frame(xx)
            
            yy$SimN<-j ##simulation number (0...)
            
            set_per_run <- rbind(set_per_run,yy)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$time<-time
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        final_set<-rbind(final_set,set_per_run)
        setwd(results_folder)
    }
    return(final_set)
}

##returns the coloumns in 'cn' with names 'cnnames' r-binding all simulations
get_ALL_results <- function(filename,cn,headerTF,deriv,cnnames,do_perc) {
    final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            
            if (do_perc) {xx[,cn]<-100*xx[,cn]}
            
            xx<-xx[which(xx$V1%%1<1e-15),c(1,cn)]
            
            colnames(xx) <- c('time',cnnames)
            
            yy<-data.frame(xx)
            
            yy$SimN<-j ##simulation number (0...)
            
            set_per_run <- rbind(set_per_run,yy)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
        
        set_per_run<-data.frame(set_per_run)
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        final_set<-rbind(final_set,set_per_run)
        setwd(results_folder)
    }
    return(final_set)
}

##for files containing frequencies by days (histo monitor type)
##returns mean and sd by day and Bin-values
## data structure:
##Bin -- value -- total -- d0 ... last day
mean_sd_dayHisto <- function(filename,headerTF,deriv) {
    final_set<-NULL
    i<-1
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(fp[i])
        sim_folders<-list.dirs(path = "." , recursive = FALSE)
        ##test
        j<-1
        set_per_run<-NULL
        for (j in 1:1) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            SAVETHIS<-xx[,2] ##value .. to save for the plot but not needed for statistics
            xx<-xx[,-(1:2)]
            namecolo<-c('total',paste0('d',0:(ncol(xx)-2)))
            colnames(xx)<-namecolo
            set_per_run <- xx
        }
        for (j in 2:length(sim_folders)) {
            vol<-paste(sim_folders[j], filename , sep="")
            xx<-data.frame(fread(vol, header=headerTF))
            xx<-xx[,-(1:2)]
            
            colnames(xx)<-namecolo
            set_per_run <- cbind(set_per_run,xx)
        }
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_'))
#         
#         kkk<-sapply( split.default(set_per_run, names(set_per_run)), function(x) {x1<-unlist(x); data.frame(mean=apply(x1,1,rowMeans),sd=apply(x1,1,rowSds) )  }  )
# t(apply( set_per_run,1,function(x) tapply(set_per_run,colnames(set_per_run),colMeans )  ))
        complete_mean<-SAVETHIS##bin-value
        complete_sd<-SAVETHIS##bin-value
        
        for (day in unique(colnames(set_per_run))) {##mean and sd by day
            tmp <- set_per_run[,which(colnames(set_per_run)==day)]
            mm <- rowMeans(tmp)
            ss<-rowSds(as.matrix(tmp))
            complete_mean<-cbind(complete_mean,mm)
            complete_sd<-cbind(complete_sd,ss)
#             kkk<-sapply( split.default(tmp, names(tmp)), function(x) {x1<-unlist(x); data.frame(mean=rowMeans(tmp),sd=rowSds(as.matrix(tmp)) )  }  )
        }
        colnames(complete_mean)<-c('BIN',namecolo)
        complete_mean<-melt( setDT(data.frame(complete_mean)),id.vars='BIN')
        names(complete_mean)[names(complete_mean)=='value']<-'mean'
        
        colnames(complete_sd)<-c('BIN',namecolo)
        complete_sd<-melt( setDT(data.frame(complete_sd)),id.vars=c("BIN"))
        names(complete_sd)[names(complete_sd)=='value']<-'sd'
        
        ##merge mean and ssd dataframes matching bin-values and variable=tot...d21
        set_per_run<-merge(complete_mean,complete_sd,by=c('BIN','variable'))
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        
        tmp<-set_per_run#[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}

###
insphere<-function(x,y,z, cx,cy,cz, r1,r2) {
    isin<-NULL
    radmin<-r1#^2
    radmax<-r2#^2
    sph <- sqrt( (x-cx)^2+(y-cy)^2+(z-cz)^2)
#     print(sph)
    isin1<-sph-radmin ##>=0 ->1
    isin2<-sph-radmax ##>=0 ->0
#     print(summary(isin))
    isin1[which(isin1>=0)]<-1
    isin2[which(isin2<0)] <-1
    isin<-data.frame(cbind(isin1,isin2))
   isin<- isin %>% 
  mutate(C = if_else(isin1 == isin2, isin1, 0))
    
#     isin[which(isin<=0)]<-paste0(color)
#     isin[which(isin>0)]<-paste0('dark ',color)#paste0('light ',color)
#     if (sph<rad) {isin<-1}
    return(sum(isin$C))
}
##returns the frequency of a type of cell at each gc-centre distance -- for now take only one sim as example
depth_in_gc <- function(filename,headerTF,celltype,cum) {
final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(paste0(fp[i],'/000'))
        sim_folders<-list.dirs(path = "." , recursive = FALSE) ##000...
        #test
        j<-1
#         zpos<-paste(fp, filename , sep="/")
        xxa<-data.frame(fread(filename, header=F))
        colnames(xxa)<-c('time','x','y','z')
        xxa$type<-celltype
        maxspace=33
        freq_depth<-NULL
        if (cum!='cumulative') {
            for (tt in unique(xxa$time)) {
                num_per_radius<-integer(maxspace)
                data_per_time<-xxa[which(xxa$time==tt),]
                data_per_time <-data_per_time %>% distinct()
                for (rad in 0:32) {
                    num_per_radius[(rad+1)]<-insphere(data_per_time$x,data_per_time$y,data_per_time$z,32,32,32,rad,(rad+1))
                }
            nn<-sum(num_per_radius)
            if (nn>0){ num_per_radius<-num_per_radius/sum(num_per_radius)}#else {num_per_radius<-0}
            freq_depth<-rbind(freq_depth,c(round(tt),num_per_radius))
            }
        } else {
            num_per_radius<-integer(maxspace)
            data_per_time<-xxa
            data_per_time <-data_per_time %>% distinct()
            for (rad in 0:32) {
                num_per_radius[(rad+1)]<-insphere(data_per_time$x,data_per_time$y,data_per_time$z,32,32,32,rad,(rad+1))
                print(num_per_radius)
            }
            nn<-sum(num_per_radius)
            if (nn>0){ num_per_radius<-num_per_radius/sum(num_per_radius)}#else {num_per_radius<-0}
#             print(num_per_radius)
            freq_depth<-num_per_radius
        }
        
        set_per_run<-freq_depth
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_')) ##split by the param characteristic --> note that this is not yet general        
        
        set_per_run<-data.frame(set_per_run)
        if (cum!='cumulative') {
        colnames(set_per_run)<-c('time',paste0('d=',0:(length(num_per_radius)-1)))
        set_per_run<-melt( setDT(data.frame(set_per_run)),id.vars='time')
        } else {set_per_run$variable<-paste0('d=',0:(length(num_per_radius)-1))}
        set_per_run$CellType<-celltype
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        tmp<-set_per_run#[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}



##returns the frequency of a type of cell at each gc-centre distance -- for now take only one sim as example
freq_in_zone <- function(filename,headerTF,celltype) {
final_set<-NULL
    for (i in 1:length(fp)) {
        print(fp[i])
        setwd(paste0(fp[i],'/000'))
        sim_folders<-list.dirs(path = "." , recursive = FALSE) ##000...
        #test
        j<-1
#         zpos<-paste(fp, filename , sep="/")
        xxa<-data.frame(fread(filename, header=F))
        colnames(xxa)<-c('time','x','y','z')
        xxa$type<-celltype
        
        tmp_time<-NULL
        for (tt in unique(xxa$time)) {
            xxT <- xxa[which(xxa$time==tt),]
            lz <- sum(xxT$z<=32)/nrow(xxT)
            dz <- sum(xxT$z>32)/nrow(xxT)
            lzdz<-lz/dz
            tmp_time<-rbind(tmp_time, c(tt,lz,dz,lzdz))
        }
        colnames(tmp_time)<-c('time','inLZ','inDZ','fr_LZDZ')
        set_per_run<-tmp_time
        
        listpar <- unlist(strsplit(fp[i],'-all'))[1]
        
        Ref <- listpar
        listpar<-unlist(strsplit(listpar,'_')) ##split by the param characteristic --> note that this is not yet general        
        
        set_per_run<-data.frame(set_per_run)
#         if (cum!='cumulative') {
#         colnames(set_per_run)<-c('time',paste0('d=',0:(length(num_per_radius)-1)))
#         set_per_run<-melt(data.frame(set_per_run),id.vars='time')
#         } else {set_per_run$variable<-paste0('d=',0:(length(num_per_radius)-1))}
        set_per_run$CellType<-celltype
        
        set_per_run$Ref <- Ref
        set_per_run$TfrModel <- listpar[1]
        set_per_run$DivMode <- listpar[2]
        set_per_run$TfrNum <- listpar[3]
        
        #set_per_run$SelfPerc <- listpar[4]
        FIRST <- (unlist(strsplit(listpar[4],'\\.'))[1])
        if (FIRST=='Self') {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[3])
            set_per_run$SelfPerc <- paste0('p_Self=0.',SelfPerc)
        } else {
            SelfPerc<-(unlist(strsplit(listpar[4],'\\.'))[4])
            set_per_run$SelfPerc <- paste0('NO Inh; p_Self=0.',SelfPerc)
        }
        
        set_per_run$IntTime <- listpar[5]
        set_per_run$Tmove <- listpar[6]
        set_per_run$HypModel <- listpar[7]
        set_per_run$TfhNum <- listpar[8]
        set_per_run$Redeem <- listpar[9]
        set_per_run$SelfMode <- listpar[10]
        
        tfr<-as.numeric(unlist(strsplit(listpar[3],'\\.'))[2])
        tfh<-as.numeric(unlist(strsplit(listpar[8],'\\.'))[2])
        tfh_by_tfr <- round(tfh/tfr)
        if (tfh_by_tfr==Inf) {set_per_run$Tfh_by_Tfr <- paste0(tfh,':0')} else {
            set_per_run$Tfh_by_Tfr <- paste0(tfh_by_tfr,':1')
        }
        tmp<-set_per_run#[,-c(1:length(sim_folders))]
        final_set<-rbind(final_set,tmp)
        setwd(results_folder)
    }
    return(final_set)
}





