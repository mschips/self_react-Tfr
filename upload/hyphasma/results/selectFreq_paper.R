
library(data.table)
# library(plot3D)
# library(scatterplot3d)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(readxl)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggbreak)
library(patchwork)
library(purrr)
library(readr)
library(ggprism)

###these are the plots for modelS comparison
##i.e. Tfr number -- Tfh Number -- p_Self and other parameters that were varied are fixed

HYPMODEL_folder<-getwd()#'/home/msi18/Desktop/diversity_plus_DZapop/hyphasma/results/'

setwd(HYPMODEL_folder)

##collect the relevant csv from the 'TfrModel'_model/'TfrModel'_folder

RESULTSET <- '^aponc'
LEGENDUSE <- 'Apoptosis'

######
allmodels<-dir(path="./",pattern=RESULTSET,full.names=F,recursive=F)

if (RESULTSET=='^semigate'||RESULTSET=='^aposg'||RESULTSET=='^sgmut') {
allrm<-dir(path="./",pattern=RESRM,full.names=F,recursive=F)
allmodels<-allmodels[-which(allmodels%in%allrm)]
if (RESULTSET=='^aposg') {
#allrm<-dir(path="./",pattern=RESRM2,full.names=F,recursive=F)
#allmodels<-allmodels[-which(allmodels%in%allrm)]
}
}


model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


allmodels<-c(allmodels,'reference_model') ##add results from '250:0' group


TAKE <- c('FS_nonself.csv','fdcdeath_frac.csv','tfhdeath.csv','fdcdeath.csv','tfhdeath_frac.csv','FS_self.csv')

filestobind<-filestobind[which(filestobind%in%TAKE)]


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'
select_pself<-c('p_Self=0.06','p_Self=0.02','p_Self=0.04')
select_pself<-c('p_Self=0.04')#,'p_Self=0.00')
select_tfhtfr<-c('250:0','19:1','10:1','5:1','2:1','1:1')

i<-1

for (i in 1:length(filestobind)) {
print(filestobind[i])

##list all files with path that are to bind, i.e. in different models folder collect all files with the same name
files <- list.files(path = allmodels,
  pattern = paste0("^",filestobind[i],"$"),
                 recursive = TRUE, full.names = TRUE)
#small check:
# if (length(files)!=length(allmodels)) {
#     print(paste("TAKING MORE FILES!!stop! nfiles, ",length(files)," allmod:",length(allmodels)));
#     stop;}
##the following will be the name of the data frame containing the info
cc_nn <- unlist(strsplit(filestobind[i],'.csv'))
##this is the data frame collecting all results
cc_file <- rbindlist(lapply(files, fread), fill = TRUE)

##select parameters
cc_file <- cc_file[which(cc_file$TfhNum%in%select_tfh & cc_file$SelfMode%in%select_selfmode & cc_file$SelfPerc%in%select_pself & cc_file$Tfh_by_Tfr%in%select_tfhtfr),]

## when Tfr are absent there is no difference between movements -- need to double for plotting purpouse
ref_mod <- cc_file[cc_file$TfrModel=='Reference',]
ref_mod$Tmove <- 'R5'
cc_file <- rbind(cc_file,ref_mod)
#

cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)
cc_file$model_move[which(cc_file$TfrModel=='Reference')] <- 'Reference'

if (LEGENDUSE!='Apoptosis') {
    cc_file <-cc_file[cc_file$Tmove!='R5',]
}

cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))

cc_file$Tfh_by_Tfr <- factor(cc_file$Tfh_by_Tfr, select_tfhtfr)

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind) ##these are the names of final dataframes (without extension)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



##############################################CAUSE OF DEATH

selecttime<-c(3,5,7,9,11,13,15,17,19,21)*24

fdcdeath_frac$Type<-'FDC'

tfhdeath_frac$Type<-'Tfh'


FS <- rbind(fdcdeath_frac,tfhdeath_frac)

death_plot <- ggplot(FS,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=Type, linetype=Type)) +
geom_line(aes(y=mean))+
# geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.1,col=NA)+
geom_errorbar(FS[FS$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~TfhNum+Tmove,scale='free')+
#             theme_bw()+
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',minor_breaks=seq(0,1,0.1))+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Freq. of apopt by cause') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='CAUSE')

###########################################################################

fdcdeath$Type<-'FDC'

tfhdeath$Type<-'Tfh'

FS <- rbind(fdcdeath,tfhdeath)

ag_plot <- ggplot(FS,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=Type, linetype=Type)) +
geom_line(aes(y=mean))+
# geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.1,col=NA)+
geom_errorbar(FS[FS$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~TfhNum+Tmove,scale='free')+
#             theme_bw()+
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',minor_breaks=seq(0,1,0.1))+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Num. of apo by cause') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='CAUSE')



#############################

FS_nonself$Type<-'ns--GCBC'

FS_self$Type<-'s--GCBC'

FS <- rbind(FS_nonself,FS_self)

freqsel_plot <- ggplot(FS,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=Type, linetype=Type)) +
geom_line(aes(y=mean))+
# geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.1,col=NA)+
geom_errorbar(FS[FS$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~TfhNum+Tmove,scale='free')+
#             theme_bw()+
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',minor_breaks=seq(0,1,0.1))+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Freq. of selected cells in GC') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Specificity')

###
AAA<-ggarrange(#bcaff_plot+theme(axis.title.x=element_blank()),
freqsel_plot,#+theme(axis.title.x=element_blank()), 
death_plot,#+theme(axis.title.x=element_blank()),
# ascaff_plot+theme(axis.title.x=element_blank()), 
labels='AUTO', common.legend=F,font.label = list(size = 20), ncol=1)


# ggsave(paste0('selectFreq',RESULTSET,select_pself,'.png'), h=12,w=14,, dpi=300)

ggsave(('FigS3.png'), h=12,w=14,, dpi=300)

