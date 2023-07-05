
library(data.table)
# library(plot3D)
# library(scatterplot3d)
# library(ggplot2)
library(dplyr)
library(matrixStats)
# library(readxl)

# library(foreach)
# library(doParallel)

HYPMODEL_folder<-getwd()#'/home/msi18/Desktop/diversity_plus_DZapop/hyphasma/results/'

setwd(HYPMODEL_folder)

folder_param = list.dirs(path = ".", recursive = FALSE)


##load functions
source('Rfunctions_forHyphasma.R')


##in 'hyphasma_model' folders
##there is one folder for each TfrModel of name type 'TfrModel'_model
##this script returns csv placed in 'TfrModel'_model/'TfrModel'_folder

allmodels<-dir(path="./",pattern='_model$',full.names=F,recursive=F)



# 
# n.cores <- min( length(allmodels), (parallel::detectCores()-1)  )
# 
# my.cluster <- parallel::makeCluster( n.cores, type='PSOCK' )
# 
# doParallel::registerDoParallel(cl=my.cluster)
curr_mod<-allmodels[1]


for (curr_mod in allmodels) {

setwd(curr_mod)
results_folder<-getwd()


model <- unlist(strsplit(curr_mod,'_model'))[1]


wheretosave<-paste0(model,'_folder')


fp<-dir(path="./",pattern='-all',full.names=F,recursive=F) ##different param files with same TfrModel


if (!dir.exists(wheretosave)) {
    dir.create(wheretosave)
    CN_file<-T
} else {
    existing <- fread(paste0(wheretosave,"/gcbc_num.csv"))
    existRef<-paste0(unique(existing$Ref),'-all')
    fp <- fp[-which(fp%in%(existRef))]
    CN_file<-F
}
# if (dir.exists(wheretosave)) {unlink(wheretosave, recursive = TRUE)}
# dir.create(wheretosave)

if (length(fp)>0) {
file_volume<-'/BC_analysis_TFR.out'
# //    time(1): live_Bcells(2): CB_DZ(3): liveCC_DZ(4): apopt_DZ(5):
# //      CB_LZ(6): liveCC_LZ(7): apopt_LZ(8):
# //      (self_CB+self_liveCC)(9): self_apopt(10):
# //      out_inGC(11): self_out_inGC(12): out(13): self_out(14):
# //      live_cell_ratio(15): live_zone_ratio(16): cell_ratio(17): zone_ratio(18):
# //      propGC_selfOut(19): prop_selfOut_inGC(20): prop_selfOut(21): CB_redeemed(22):
# //      cd138/liveBC(23) : TfrSplTfh(/apop) (24) : total_out_inGC(25) :cd138_n(26)

file_aff <- '/BCaff_TFR.out'
# //0 : 1 : 2 --> sout (self : total : nonSelf)
#     //3 --> total bcs
# //    time(1): sout--Self[m:sd](2:3): total[m:sd](4:5): nonself[m:sd](6:7)
# //      liveBCs[m:sd](8:9)

file_tfr<-'/tfr.out'

file_newly <- '/newly_gen.out'

file_division <- '/ndivtime.out'

file_ag <- '/antigen.out'

file_ccl3 <- '/ccl3_file.out'

file_frac<-'/FRAC_SEL.out'
# FRAC_SEL << time << " " << Aselect_frac_ns
#                 << " " << Aselect_frac_self << "\n";

file_apo<-'/CAUSE_OF_DEATH.out'


##selfGCBC (% GCBC)
##--mean/sd
perc_self_gcbc <- get_fraction_results(file_volume, 9,0,2,0, F,0, 'percentage')
write.table(perc_self_gcbc, paste0(wheretosave,"/perc_self_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_perc_self_gcbc <- get_ALL_frac_results(file_volume, 9,0,2,0, F,0, 'percentage')
write.table(S_perc_self_gcbc, paste0(wheretosave,"/S_perc_self_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)

##redeemedGCBC (% GCBC)
##--mean/sd
perc_redeemed_gcbc <- get_fraction_results(file_volume, 22,0,2,0, F,0, 'percentage')
write.table(perc_redeemed_gcbc, paste0(wheretosave,"/perc_redeemed_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)

##--single sim values
S_perc_redeemed_gcbc <- get_ALL_frac_results(file_volume, 22,0,2,0, F,0, 'percentage')
write.table(S_perc_redeemed_gcbc, paste0(wheretosave,"/S_perc_redeemed_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)


##selfApo (% Apo)
perc_self_apo <- get_fraction_results(file_volume, 10,0,5,8, F,0, 'percentage')
write.table(perc_self_apo, paste0(wheretosave,"/perc_self_apo.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_perc_self_apo <- get_ALL_frac_results(file_volume, 10,0,5,8, F,0, 'percentage')
write.table(S_perc_self_apo, paste0(wheretosave,"/S_perc_self_apo.csv"), row.names=FALSE, append=T,col.names=CN_file)



##GC file_volume
gcbc_num <- get_results(file_volume, 2, F, 0, 0)
write.table(gcbc_num, paste0(wheretosave,"/gcbc_num.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_gcbc_num <- get_ALL_results(file_volume, 2, F,0, 'mean', 0)
write.table(S_gcbc_num, paste0(wheretosave,"/S_gcbc_num.csv"), row.names=FALSE, append=T,col.names=CN_file)


##number of Out cells
##--mean/sd
cumul_out <- get_results(file_volume, 13, F,0, 0)
write.table(cumul_out, paste0(wheretosave,"/cumul_out.csv"), row.names=FALSE, append=T,col.names=CN_file)

##--single sim values
S_cumul_out <- get_ALL_results(file_volume, 13, F,0, 'mean', 0)
write.table(S_cumul_out, paste0(wheretosave,"/S_cumul_out.csv"), row.names=FALSE, append=T,col.names=CN_file)

daily_out <- get_results(file_volume, 13, F,72, 0)
write.table(daily_out, paste0(wheretosave,"/daily_out.csv"), row.names=FALSE, append=T,col.names=CN_file)

cumul_NSout <- get_diff_results(file_volume,13,14,F,0)
write.table(cumul_NSout, paste0(wheretosave,"/cumul_NSout.csv"), row.names=FALSE, append=T,col.names=CN_file)

##--single sim values
S_cumul_NSout <- get_ALLDiff_results(file_volume, 13,14, F, 'mean')
write.table(S_cumul_NSout, paste0(wheretosave,"/S_cumul_NSout.csv"), row.names=FALSE, append=T,col.names=CN_file)

##selfOut (% Out)
perc_self_out <- get_fraction_results(file_volume, 14,0,13,0, F,0, 'percentage')
write.table(perc_self_out, paste0(wheretosave,"/perc_self_out.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_perc_self_out <- get_ALL_frac_results(file_volume, 14,0,13,0, F,0, 'percentage')
write.table(S_perc_self_out, paste0(wheretosave,"/S_perc_self_out.csv"), row.names=FALSE, append=T,col.names=CN_file)


##affinity of output cells
##--mean and sd of means
aff_out_self <- get_results(file_aff, 2, F, 0, 0)
aff_out_self$Specificity<-'Self'
aff_out_tot <- get_results(file_aff, 4, F, 0, 0)
aff_out_tot$Specificity<-'Total'
aff_out_nonself <- get_results(file_aff, 6, F, 0, 0)
aff_out_nonself$Specificity<-'Foreign'
aff_out<- rbind(aff_out_self,aff_out_nonself,aff_out_tot)
write.table(aff_out, paste0(wheretosave,"/aff_out.csv"), row.names=FALSE, append=T,col.names=CN_file)

##--single sim values
aff_out_mean <- get_ALL_results(file_aff, c(2,4,6), F, 0, c('Self','Total','Foreign'), 0)
aff_out_sd <- get_ALL_results(file_aff, c(3,5,7), F, 0, c('Self','Total','Foreign'), 0)

meltby<-colnames(setDT(aff_out_mean))[-(2:4)]
affOutMean<-melt(setDT(aff_out_mean),id.vars=meltby)

affOutSd<-melt(setDT(aff_out_sd),id.vars=meltby)

mergeby<-colnames(affOutMean)[-which(colnames(affOutMean)=='value')]

S_aff_out<-merge(affOutMean,affOutSd,by=mergeby,all=T)
names(S_aff_out)[names(S_aff_out)%in%c('variable','value.x','value.y')]<-c('Specificity','mean','sd')

write.table(S_aff_out, paste0(wheretosave,"/S_aff_out.csv"), row.names=FALSE, append=T,col.names=CN_file)


##CD138
##--mean/sd
perc_cd138 <- get_results(file_volume, 23, F,0, 1)
write.table(perc_cd138, paste0(wheretosave,"/perc_cd138.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_perc_cd138 <- get_ALL_results(file_volume, 23, F,0, 'mean', 1)
write.table(S_perc_cd138, paste0(wheretosave,"/S_perc_cd138.csv"), row.names=FALSE, append=T,col.names=CN_file)


##GC file_volume
gcbc_num <- get_results(file_volume, 2, F, 0, 0)
write.table(gcbc_num, paste0(wheretosave,"/gcbc_num.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_gcbc_num <- get_ALL_results(file_volume, 2, F,0, 'mean', 0)
write.table(S_gcbc_num, paste0(wheretosave,"/S_gcbc_num.csv"), row.names=FALSE, append=T,col.names=CN_file)


##tfr
tfr_num<-get_results(file_tfr, 2, F, 0, 0)
write.table(tfr_num, paste0(wheretosave,"/tfr_num.csv"), row.names=FALSE, append=T,col.names=CN_file)




##--single sim values
aff_out_mean <- get_ALL_results(file_aff, c(2,4,6), F, 0, c('Self','Total','Foreign'), 0)
aff_out_sd <- get_ALL_results(file_aff, c(3,5,7), F, 0, c('Self','Total','Foreign'), 0)

meltby<-colnames(setDT(aff_out_mean))[-(2:4)]
affOutMean<-melt(setDT(aff_out_mean),id.vars=meltby)

affOutSd<-melt(setDT(aff_out_sd),id.vars=meltby)

mergeby<-colnames(affOutMean)[-which(colnames(affOutMean)=='value')]

S_aff_out<-merge(affOutMean,affOutSd,by=mergeby,all=T)
names(S_aff_out)[names(S_aff_out)%in%c('variable','value.x','value.y')]<-c('Specificity','mean','sd')

write.table(S_aff_out, paste0(wheretosave,"/S_aff_out.csv"), row.names=FALSE, append=T,col.names=CN_file)


##n TFR:BC contacts
file1<-'/nTFRcontacts_SelfCCselected.out'
Sccsel<-mean_sd_dayHisto(file1,F,0)
Sccsel$Mut<-'Self'
Sccsel$State<-'Selected'

file2<-'/nTFRcontacts_CCselected.out'
ccsel<-mean_sd_dayHisto(file2,F,0)
ccsel$Mut<-'Foreign'
ccsel$State<-'Selected'

file3<-'/nTFRcontacts_SelfCCdeleted.out'
Sccdel<-mean_sd_dayHisto(file3,F,0)
Sccdel$Mut<-'Self'
Sccdel$State<-'Apoptotic'

file4<-'/nTFRcontacts_CCdeleted.out'
ccdel<-mean_sd_dayHisto(file4,F,0)
ccdel$Mut<-'Foreign'
ccdel$State<-'Apoptotic'

tfr_b_cont <- rbind(Sccsel,Sccdel,ccsel,ccdel)

write.table(tfr_b_cont, paste0(wheretosave,"/tfr_b_cont.csv"), row.names=FALSE, append=T,col.names=CN_file)


##affinity of live GCBCs
##--mean/sd
aff_gcbc <- get_results(file_aff, 8, F, 0, 0)
write.table(aff_gcbc, paste0(wheretosave,"/aff_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)
##--single sim values
S_aff_gcbc <- get_ALL_results(file_aff, c(8,9), F, 0, c('mean','sd'), 0)
write.table(S_aff_gcbc, paste0(wheretosave,"/S_aff_gcbc.csv"), row.names=FALSE, append=T,col.names=CN_file)


###depth from GC centre
cd138<-depth_in_gc('cd138_xyz.out',F,'CD138','tt')
asc<-depth_in_gc('ASC_xyz.out',F,'ASC','tt')
if (curr_mod!='reference_model'&&curr_mod!='refnored_model'&&curr_mod!='reflong_model'&&curr_mod!='longout_model') {


##tfr
tfr_num<-get_results(file_tfr, 2, F, 0, 0)
write.table(tfr_num, paste0(wheretosave,"/tfr_num.csv"), row.names=FALSE, append=T,col.names=CN_file)



##free Tfr (% all-tfr)
##--mean/sd
perc_free_tfr <- get_fraction_results(file_tfr, 3,0,2,0, F,0, 'fract')
write.table(perc_free_tfr, paste0(wheretosave,"/perc_free_tfr.csv"), row.names=FALSE, append=T,col.names=CN_file)

##positions
tfr<-depth_in_gc('tfr_xyz.out',F,'Tfr','tt')
depth<-rbind(tfr,cd138,asc)
freq_zone_tfr<- freq_in_zone('tfr_xyz.out',F,'Tfr')

write.table(freq_zone_tfr, paste0(wheretosave,"/freq_zone_tfr.csv"), row.names=FALSE, append=T,col.names=CN_file)

} else {depth<-rbind(cd138,asc)}

write.table(depth, paste0(wheretosave,"/depth.csv"), row.names=FALSE, append=T,col.names=CN_file)

###exitpoint of ASC
exit<-depth_in_gc('exitpoint_xyz.out',F,'ASC','cumulative')

write.table(exit, paste0(wheretosave,"/exitpoint.csv"), row.names=FALSE, append=T,col.names=CN_file)


##frequency of all cells type
zpos_freq <- get_ALL_results('/Zpos.out',2:8, F, 0, c('Zpos', 'CB', 'CC', 'Tfh', 'Tfr', 'CD138','ASC'), 0)

write.table(zpos_freq, paste0(wheretosave,"/zpos_freq.csv"), row.names=FALSE, append=T,col.names=CN_file)

##mutation session
##mutation probability
mutfile<-'/mutation_time.out'
pmut <- get_ALL_results(mutfile, c(4,5), F,0, c('mean','sd'),0)

pmut$Mut<-'Foreign'

mutfileself<-'/Selfmutation_time.out'
pmut_self <- get_ALL_results(mutfileself, c(4,5), F,0, c('mean','sd'),0)

pmut_self$Mut<-'Self'
prob_mut<-rbind(pmut,pmut_self)
prob_mut<-prob_mut[which(prob_mut$time%%24<1e-15),]

write.table(prob_mut, paste0(wheretosave,"/prob_mut.csv"), row.names=FALSE, append=T,col.names=CN_file)


################################################
####
FS_nonself<-get_results(file_frac, 2, F, 0, 0)
write.table(FS_nonself, paste0(wheretosave,"/FS_nonself.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('1 of 8')

####
FS_self<-get_results(file_frac, 3, F, 0, 0)
write.table(FS_self, paste0(wheretosave,"/FS_self.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('2 of 8')

####
numFS_nonself<-get_results(file_frac, 4, F, 0, 0)
write.table(numFS_nonself, paste0(wheretosave,"/numFS_nonself.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('3 of 8')

####
numFS_self<-get_results(file_frac, 5, F, 0, 0)
write.table(numFS_self, paste0(wheretosave,"/numFS_self.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('4 of 8')

###################
####
fdcdeath<-get_results(file_apo, 2, F, 0, 0)
write.table(fdcdeath, paste0(wheretosave,"/fdcdeath.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('5 of 8')

####
fdcdeath_frac<-get_results(file_apo, 3, F, 0, 0)
write.table(fdcdeath_frac, paste0(wheretosave,"/fdcdeath_frac.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('6 of 8')

####
tfhdeath<-get_results(file_apo, 4, F, 0, 0)
write.table(tfhdeath, paste0(wheretosave,"/tfhdeath.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('7 of 8')

####
tfhdeath_frac<-get_results(file_apo, 5, F, 0, 0)
write.table(tfhdeath_frac, paste0(wheretosave,"/tfhdeath_frac.csv"), row.names=FALSE, append=T,col.names=CN_file)

print('8 of 8')


print('DONE')

setwd(HYPMODEL_folder) ##return to parental folder and proceed with next 'TfrModel'
} else {
    print("ALL DONE")

    setwd(HYPMODEL_folder) ##return to parental folder and proceed with next 'TfrModel'

    }
}

