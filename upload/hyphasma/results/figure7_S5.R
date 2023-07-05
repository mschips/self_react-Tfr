
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


RESULTSET <- '^semigate.38'
LEGENDUSE <- 'SemiGate.38'
# #


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


TAKE <- c('aff_gcbc.csv','perc_self_gcbc.csv','cumul_out.csv','cumul_NSout.csv','perc_self_out.csv','aff_out.csv','gcbc_num.csv','S_cumul_out.csv','S_cumul_NSout.csv','S_perc_self_out.csv','S_aff_gcbc.csv','S_aff_out.csv','S_gcbc_num.csv')

filestobind<-filestobind[which(filestobind%in%TAKE)]


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'

##change the following line to generate results for different pSelf probability
select_pself<-c('p_Self=0.04')
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



selecttime<-seq(3,21,3)*24 #where to plot errorbars

jw<-.3
dw<-1.3

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)

###SUPPLEMENTARY FIGURES


perc_self_gcbc_plot <- ggplot(perc_self_gcbc,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=LEGENDUSE, linetype=LEGENDUSE)) +
geom_line(aes(y=mean))+
geom_errorbar(perc_self_gcbc[perc_self_gcbc$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('s-BCs (% BCs)') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Tfr Model')

# perc_self_gcbc_plot

####


perc_self_out_plot <- ggplot(perc_self_out,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=LEGENDUSE, linetype=LEGENDUSE)) +
geom_line(aes(y=mean))+
geom_errorbar(perc_self_out[perc_self_out$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
xlab('Days p. GC onset') + ylab('s-ASCs (% ASCs)') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Tfr Model')

# perc_self_out_plot

####

asc_cumul_plot <- ggplot(cumul_out,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=LEGENDUSE, linetype=LEGENDUSE)) +
geom_line(aes(y=mean))+
geom_errorbar(cumul_out[cumul_out$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Cumulative # of ASCs') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Tfr Model')


####

NSasc_cumul_plot <- ggplot(cumul_NSout,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=LEGENDUSE, linetype=LEGENDUSE)) +
geom_line(aes(y=mean))+
geom_errorbar(cumul_NSout[cumul_NSout$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor',minor_breaks=seq(0,21,2.5))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Cumulative # of ns-ASCs') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Tfr Model')


####################MAIN FIGURE

gcbc_num_plot <- ggplot(gcbc_num,
aes(x=time/24, col=Tfh_by_Tfr, fill=Tfh_by_Tfr, shape=LEGENDUSE, linetype=LEGENDUSE)) +
geom_line(aes(y=mean))+
# geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.1,col=NA)+
geom_errorbar(gcbc_num[gcbc_num$time%in%selecttime], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~Tmove,scale='free')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor')+#,limits=c(scale_min,scale_max))+
            scale_x_continuous(guide='prism_minor')+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
xlab('Days p. GC onset') + ylab('# Live BCs') +
            labs(fill='Tfh:Tfr',color='Tfh:Tfr',linetype='Tfr Model')
            
####

selectday <- 504

scNSout <- S_cumul_NSout[which(S_cumul_NSout$time%in%selectday),]

scale_min <- round(min(scNSout$mean))-1
scale_max <- round(max(scNSout$mean))+1


NSoutnum_box_plot<- ggboxplot(data=scNSout,x='Tfh_by_Tfr',y='mean',
# col='black', shape=21, 
facet.by=list('Tmove','time'), add='dotplot', width=.3,
add.params=list(size=.2,color='Tfh_by_Tfr',alpha=.5))+
stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')

facet_wrap(~Tmove,scale='free')+
            ylab(paste0('# ns-ASCs (day ',selectday/24,')')) +
            xlab('Tfh:Tfr')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',limits=c(scale_min,scale_max))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
            labs(fill='Tfh:Tfr')

NSoutnum_box_plot<- ggboxplot(data=scNSout,x='Tfh_by_Tfr',y='mean',
# col='black', shape=21, 
facet.by=list('Tmove','time'), add='dotplot', width=.3,
add.params=list(size=.2,color='Tfh_by_Tfr',alpha=.5))+
stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')

facet_wrap(~Tmove,scale='free')+
            ylab(paste0('# ns-ASCs (day ',selectday/24,')')) +
            xlab('Tfh:Tfr')+
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',limits=c(scale_min,scale_max))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
            labs(fill='Tfh:Tfr')

####

saffNSout <- S_aff_out[which(S_aff_out$time%in%selectday & S_aff_out$Specificity=='Foreign'),]

scale_min <- (min(saffNSout[saffNSout$time==selectday,]$mean))#-.01
scale_max <- (max(saffNSout[saffNSout$time==selectday,]$mean))#+.01

outaff_box_plot<- ggboxplot(data=saffNSout,x='Tfh_by_Tfr',y='mean', 
# col='black', shape=21, 
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.2,color='Tfh_by_Tfr',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')
            
            facet_wrap(~Tmove,scale='free')+
            ylab(paste0('Aff. ns-ASCs (day ',selectday/24,')')) + xlab('Tfh:Tfr') +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_y_continuous(guide='prism_minor',limits=c(scale_min,scale_max))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
#             ylim(scale_min,scale_max)+
            labs(fill='Tfh:Tfr')

####
spercSout <- S_perc_self_out[which(S_perc_self_out$time%in%selectday),]
scale_min <- round(min(spercSout$xx))-1
scale_max <- round(max(spercSout$xx))+1


outself_box_plot<- ggboxplot(data=spercSout,x='Tfh_by_Tfr',y='xx', 
# col='black', shape=21, 
facet.by='Tmove', add='dotplot', width=.3,
add.params=list(size=.2,color='Tfh_by_Tfr',alpha=.5))+
            stat_compare_means(label=c('p.signif'),ref.group='250:0',method='wilcox.test',hide.ns=T)+#group.by=c('TfhNum')
            
            facet_wrap(~Tmove,scale='free')+ylab(paste0('%s-ASCs (day ',selectday/24,')')) +xlab('Tfh:Tfr') +
            theme_prism(base_fontface='plain',base_size=20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), legend.text = element_text(size=20))+
            scale_fill_manual(values=cbbPalette) +
            scale_colour_manual(values=cbbPalette)+
#             ylim(0.55,0.7)+
            scale_y_continuous(guide='prism_minor',limits=c(scale_min,scale_max))+
            labs(fill='Tfh:Tfr')
#########################################################

if (LEGENDUSE!='Apoptosis') {
    H_PLOT <- 10
    W_PLOT <- 12
    MAIN_FIG<-ggarrange(gcbc_num_plot,
        NSoutnum_box_plot,#+theme(axis.title.x=element_blank()),
        outaff_box_plot,#+theme(axis.title.x=element_blank()),
        outself_box_plot,
        labels='AUTO', common.legend=T,font.label = list(size = 20))

        if (LEGENDUSE=='SemiGate') {figtitle<-'Figure6.png'; supptitle<-'FigS4.png'}
        else if (LEGENDUSE=='SemiGate.38') {figtitle<-'Figure7.png'; supptitle<-'FigS5.png'}
#     ggsave(paste0('P_',RESULTSET,select_pself,'.png'),h=H_PLOT,w=W_PLOT, dpi=300)
    ggsave(figtitle,h=H_PLOT,w=W_PLOT, dpi=300)
    ###
    SUPP_FIG<-ggarrange(
        perc_self_gcbc_plot+theme(axis.title.x=element_blank()), 
        perc_self_out_plot+theme(axis.title.x=element_blank()),
        asc_cumul_plot+theme(axis.title.x=element_blank()),
        NSasc_cumul_plot+theme(axis.title.x=element_blank()),
        labels='AUTO', common.legend=T,font.label = list(size = 20))

    KKK<-grid.arrange(grobs=list(SUPP_FIG))
    
    annotate_figure(KKK, bottom=text_grob('Days p. GC onset',size=20))

    ggsave(supptitle,h=H_PLOT,w=W_PLOT, dpi=300)
}
###

