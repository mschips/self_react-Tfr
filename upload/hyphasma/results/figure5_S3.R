
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
##i.e. Tfr # -- Tfh Number -- p_Self and other parameters that were varied are fixed

HYPMODEL_folder<-getwd()#'/home/msi18/Desktop/diversity_plus_DZapop/hyphasma/results/'

setwd(HYPMODEL_folder)

##collect the relevant csv from the 'TfrModel'_model/'TfrModel'_folder

saveA<-'^aponc'
saveB<-'ApoNC'

######
allmodels<-dir(path="./",pattern=saveA,full.names=F,recursive=F)


if (saveB=='SemiGate'||saveB=='ApoSG') {
    allmodels2<-dir(path="./",pattern=remove,full.names=F,recursive=F)
    allmodels<-allmodels[-which(allmodels%in%allmodels2)]
}


model <- unlist(strsplit(allmodels[1],'_model'))[1]



##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


allmodels<-c(allmodels,'DTRexp_model','reference_model')


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'
select_pself<-c('p_Self=0.06','p_Self=0.02','p_Self=0.04')
select_pself<-c('p_Self=0.04')

select_tfhtfr<-c('2:1','250:0')

dtrorder <- c('NONE','0','7','9')

select_model <- c(saveB,'Reference')

i<-1


aa<-c('S_gcbc_num.csv', 'S_perc_self_gcbc.csv', 'S_cumul_out.csv', 'cumul_out.csv', 'S_perc_self_out.csv', 'outself.csv', 'gcbc_num.csv', 'perc_self_gcbc.csv', 'asc_cumul.csv', 'perc_self_out.csv','tfr_num.csv','S_aff_out.csv')

filestobind <- filestobind[which(filestobind%in%aa)]


for (i in 1:length(filestobind)) {
print(filestobind[i])

##list all files with path that are to bind, i.e. in different models folder collect all files with the same name
files <- list.files(path = allmodels,
  pattern = paste0("^",filestobind[i],"$"),
                 recursive = TRUE, full.names = TRUE)
#small check:
if (length(files)!=length(allmodels)) {
    print(paste("TAKING MORE FILES!!stop! nfiles, ",length(files)," allmod:",length(allmodels)));
    stop;}
##the following will be the name the data frame containing the info
cc_nn <- unlist(strsplit(filestobind[i],'.csv'))
##this is the data frame collecting all results
cc_file <- rbindlist(lapply(files, fread), fill = TRUE)

cc_file <- cc_file[which(cc_file$TfhNum%in%select_tfh & cc_file$SelfMode%in%select_selfmode & cc_file$SelfPerc%in%select_pself & cc_file$Tfh_by_Tfr%in%select_tfhtfr),]


## when Tfr are absent there is no difference between movements -- need to double for plotting purpouse
if (saveB=='ApoNC'){
ref_mod <- cc_file[cc_file$TfrModel=='Reference',]
ref_mod$Tmove <- 'R5'
cc_file <- rbind(cc_file,ref_mod)
} else {cc_file<-cc_file[cc_file$Tmove=='NW',]}
###

# cc_file$model_move[which(cc_file$TfrModel=='Reference')] <- 'Reference'

cc_file$DTR_day <- 'NONE'

##taking depletion day for each entry
ttt <- seq_along(unlist(strsplit(cc_file$TfrModel[-which(cc_file$TfrModel %in% select_model)],'.Depl')))%%2==0

##and corresponding TfrModel
mmm <- seq_along(unlist(strsplit(cc_file$TfrModel[-which(cc_file$TfrModel %in% select_model)],'.Depl')))%%2>0

##attribute DTR=depl time
cc_file[-which(cc_file$TfrModel %in% select_model),]$DTR_day <- as.numeric(unlist( strsplit( cc_file$TfrModel[-which(cc_file$TfrModel %in% select_model)],'.Depl') )[ttt])/24
##eliminate other models
cc_file[-which(cc_file$TfrModel %in% select_model),]$TfrModel <- (unlist( strsplit( cc_file$TfrModel[-which(cc_file$TfrModel %in% select_model)],'.Depl') )[mmm])

cc_file<-cc_file[which(cc_file$TfrModel%in%select_model),]


##for reference it is depl at d0
cc_file[which(cc_file$TfrModel %in% select_model[2]),]$DTR_day <- 0



cc_file$TfrModel<-select_model[1] ##dtr at day 0 is the same for all models


cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))


cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)



cc_file$Tfh_by_Tfr <- factor(cc_file$Tfh_by_Tfr, select_tfhtfr)

cc_file$DTR_day<- factor(cc_file$DTR_day, dtrorder)

cc_file$time <- as.numeric(as.character(cc_file$time))
cc_file$mean <- as.numeric(as.character(cc_file$mean))
cc_file$std <- as.numeric(as.character(cc_file$std))
##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind) ##these are the names final dataframes (without extension)


dayofanalysis <- 18

##Tfr #
tfr_num_plot <- ggplot(tfr_num[tfr_num$Tmove=='Tfh-like',],
aes(x=time/24, col=DTR_day, fill=DTR_day, shape=DTR_day)) +
geom_line(aes(y=mean))+
geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.5)+
#                     geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean))+
xlab('Days p. GC onset') + ylab('# Tfr') +
theme_prism(base_fontface='plain',base_size=15,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=15))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
#             guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            labs(fill='Depl. day',color='Depl. day')

#tfr_num_plot

#ggsave('FigS4.png',h=4.5,w=6, dpi=300)

#####################

selecttime<-c(dayofanalysis)*24
saffout <- S_aff_out[which(S_aff_out$time%in%selecttime & S_aff_out$Specificity=='Foreign'),]

saffout <- saffout[saffout$SimN<20,]


outaff_box_plot<- ggboxplot(data=saffout,x='DTR_day',y='mean', 
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.4, color='DTR_day',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='NONE',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')
            
            facet_wrap(~Tmove,scale='free_y')+ylab(paste0('Aff. ns-ASC (day ',selecttime/24,')')) +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(legend.position="none") +
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            #scale_y_continuous(guide='prism_minor')+
            guides(colour=guide_legend(override.aes=list(linewidth=1,size=3.5)))+
            xlab('Tfr depl. condition')
            
###
sout<-S_cumul_out[which(S_cumul_out$time%in%selecttime),]

sout <- sout[sout$SimN<20,]

outnum_box_plot<- ggboxplot(data=sout,x='DTR_day',y='mean',
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.4, color='DTR_day',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='NONE',method='wilcox.test',hide.ns=T)+
            facet_wrap(~Tmove,scale='free_y')+ylab(paste0('# ASC (day ',selecttime/24,')')) +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(legend.position="none") +
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            guides(colour=guide_legend(override.aes=list(linewidth=1,size=3.5)))+
            
            xlab('Tfr depl. condition')

###
spercout<-S_perc_self_out[which(S_perc_self_out$time%in%selecttime),]

spercout <- spercout[spercout$SimN<20,]

outself_box_plot<- ggboxplot(data=spercout,x='DTR_day',y='xx', 
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.4, color='DTR_day',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='NONE',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')
            
            facet_wrap(~Tmove,scale='free_y')+ylab(paste0('%s-ASC (day ',selecttime/24,')')) +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(legend.position="none") +
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            guides(colour=guide_legend(override.aes=list(linewidth=1,size=3.5)))+
            xlab('Tfr depl. condition')
            
###
sgcbcnum<-S_gcbc_num[which(S_gcbc_num$time%in%selecttime),]

sgcbcnum <- sgcbcnum[sgcbcnum$SimN<20,]

gcbc_num_box_plot<- ggboxplot(data=sgcbcnum,x='DTR_day',y='mean', 
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.4, color='DTR_day',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='NONE',method='wilcox.test',hide.ns=T)+
            facet_wrap(~Tmove,scale='free_y')+ylab(paste0('# GCBC (day ',selecttime/24,')')) +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(legend.position="none") +
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            guides(colour=guide_legend(override.aes=list(linewidth=1,size=3.5)))+
            xlab('Tfr depl. condition')

###

gcbc_num_plot <- ggplot(gcbc_num,
aes(x=time/24, col=DTR_day, fill=DTR_day, shape=DTR_day)) +
geom_line(aes(y=mean))+
geom_errorbar(gcbc_num[gcbc_num$time%%72<1e-15,],mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+
xlab('Days p. GC onset') + ylab('# Live GCBC') +
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            xlab('Tfr depl. condition')
            
            
###
perc_self_gcbc_plot <- ggplot(perc_self_gcbc,
aes(x=time/24, col=DTR_day, fill=DTR_day))+
geom_line(aes(y=mean))+
geom_errorbar(perc_self_gcbc[perc_self_gcbc$time%%72<1e-15,], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
xlab('Days p. GC onset') + ylab('s-GCBC (% GCBC)')


###

asc_cumul_plot <- ggplot(cumul_out,
aes(x=time/24, col=DTR_day, fill=DTR_day)) +
geom_line(aes(y=mean))+
geom_errorbar(cumul_out[cumul_out$time%%72<1e-15],mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+

xlab('Days p. GC onset') + ylab('Cumulative # ASC') +

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            xlab('Tfr depl. condition')
            

perc_self_out_plot <- ggplot(perc_self_out,
aes(x=time/24, col=DTR_day, fill=DTR_day)) +
geom_line(aes(y=mean))+
geom_errorbar(perc_self_out[perc_self_out$time%%72<1e-15,],mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+
xlab('Days p. GC onset') + ylab('s-ASC (% ASC)') +
            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            
theme(prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            xlab('Days p. GC onset')
            


{
    H_PLOT <- 15
    W_PLOT <- 12
    MAIN_FIG<- ggarrange(
 outnum_box_plot+theme(axis.title.x=element_blank()), 
 outaff_box_plot+theme(axis.title.x=element_blank()), 
 outself_box_plot,
 labels='AUTO', 
 common.legend=F,font.label = list(size = 20), ncol=1 )
# KKK<-grid.arrange(grobs=list(AAA))

    ggsave('Figure5.png',h=H_PLOT,w=W_PLOT, dpi=300)
    ###
    SUPP_FIG<-ggarrange(tfr_num_plot,
 perc_self_gcbc_plot,#+theme(axis.title.x=element_blank()),
 perc_self_out_plot,
 labels='AUTO', common.legend=T,font.label = list(size = 20),ncol=1)

# KKK<-grid.arrange(grobs=list(AAA))

#     annotate_figure(KKK, bottom=text_grob('Days p. GC onset',size=20))

    ggsave('FigS3.png',h=H_PLOT,w=W_PLOT, dpi=300)
    
} 
