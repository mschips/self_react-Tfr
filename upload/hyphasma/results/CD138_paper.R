
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
# # 

# 
RESULTSET <- '^semigate'
RESRM <- '^semigate.3'
LEGENDUSE <- 'SemiGate'
# #
######
allmodels<-dir(path="./",pattern=RESULTSET,full.names=F,recursive=F)

# if (RESULTSET=='^semigate'||RESULTSET=='^aposg'||RESULTSET=='^sgmut') {
# allrm<-dir(path="./",pattern=RESRM,full.names=F,recursive=F)
# allmodels<-allmodels[-which(allmodels%in%allrm)]
# if (RESULTSET=='^aposg') {
# #allrm<-dir(path="./",pattern=RESRM2,full.names=F,recursive=F)
# #allmodels<-allmodels[-which(allmodels%in%allrm)]
# }
# }


model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


allmodels<-c(allmodels,'reference_model') ##add results from '250:0' group


TAKE <- c('gcbc_num.csv','cumul_out.csv','daily_out.csv','perc_cd138.csv','S_gcbc_num.csv','S_cumul_out.csv','S_perc_cd138.csv')

filestobind<-filestobind[which(filestobind%in%TAKE)]


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'

select_pself<-c('p_Self=0.04')
full_tfhtfr<-c('250:0','19:1','10:1','5:1','2:1','1:1')
select_tfhtfr<-c('250:0','5:1','2:1')

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

# if (LEGENDUSE!='Apoptosis') {
#     cc_file <-cc_file[cc_file$Tmove!='R5',]
# }
cc_file <-cc_file[cc_file$Tmove!='R5',]

cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))

cc_file$Tfh_by_Tfr <- factor(cc_file$Tfh_by_Tfr, select_tfhtfr)

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind) ##these are the names of final dataframes (without extension)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


COLORIND <- which(full_tfhtfr%in%select_tfhtfr)

COLORSELECTION <- cbbPalette[COLORIND]

selecttime <- c(4,14)
S_perc_cd138$time <- S_perc_cd138$time/24
sp38 <- S_perc_cd138[which(S_perc_cd138$time%in%selecttime),]



if (LEGENDUSE=='Apoptosis') {
    cd138_box_plot<- ggboxplot(data=sp38,x='Tfh_by_Tfr',y='mean',
    col='black', shape='Tfh_by_Tfr',
    facet.by=list('Tmove','time'), add='dotplot', width=.3
    ,add.params=list(size=.4,shape=21,fill='Tfh_by_Tfr'),
    )+
    stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+

    facet_wrap(~Tmove+time,ncol=2,scale='free')+
                ylab('% GC-PC') +
                
                ylim(0,15)+
                theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
                theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'),legend.position = "none")+
                scale_y_continuous(guide='prism_minor')+
                scale_fill_manual(values=COLORSELECTION) +
                scale_colour_manual(values=COLORSELECTION)+
                ggtitle(LEGENDUSE)+
                xlab('Tfh:Tfr')

    
    ggsave('FigS8.png',h=5,w=12, dpi=300)

} else {
    SG <- sp38[sp38$TfrModel%in%c('SemiGate','Reference'),]
    SG38 <- sp38[sp38$TfrModel%in%c('SemiGate.38','Reference'),]
    
    SGplot <- ggboxplot(data=SG,x='Tfh_by_Tfr',y='mean',
    col='black', shape='Tfh_by_Tfr',
    facet.by=list('Tmove','time'), add='dotplot', width=.3
    ,add.params=list(size=.4,shape=21,fill='Tfh_by_Tfr'),
    )+
    stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+

    facet_wrap(~Tmove+time,ncol=2,scale='free')+
                ylab('% GC-PC') +
                
                ylim(0,15)+
                theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
                theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'),legend.position = "none")+
                scale_y_continuous(guide='prism_minor')+
                scale_fill_manual(values=COLORSELECTION) +
                scale_colour_manual(values=COLORSELECTION)+
#                 ggtitle(unique(SG$TfrModel)[2])+
                xlab('Tfh:Tfr')
    
    SG38plot <- ggboxplot(data=SG38,x='Tfh_by_Tfr',y='mean',
    col='black', shape='Tfh_by_Tfr',
    facet.by=list('Tmove','time'), add='dotplot', width=.3
    ,add.params=list(size=.4,shape=21,fill='Tfh_by_Tfr'),
    )+
    stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+

    facet_wrap(~Tmove+time,ncol=2,scale='free')+
                ylab('% GC-PC') +
                
                ylim(0,15)+
                theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
                theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'),legend.position = "none")+
                scale_y_continuous(guide='prism_minor')+
                scale_fill_manual(values=COLORSELECTION) +
                scale_colour_manual(values=COLORSELECTION)+
                ggtitle(unique(SG38$TfrModel)[2])+
                xlab('Tfh:Tfr')
                
    AAA<-ggarrange(
        SGplot+ggtitle(unique(SG$TfrModel)[2])+theme(axis.title.x=element_blank()), 
        SG38plot+ggtitle(unique(SG38$TfrModel)[2]),#+theme(axis.title.x=element_blank()),
        labels='AUTO', common.legend=F,font.label = list(size = 20), ncol=1)
        
        
 ggsave('Figure8.png',h=10,w=12, dpi=300)

    
}


