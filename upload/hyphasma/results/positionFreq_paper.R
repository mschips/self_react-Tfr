

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

gett<-c('semigatenw_model','semigater5_model')
# gett<-c('play_model','aponcext_model','aponcnw_model','aponcr5_model')


# gett<-c('aponocontR5_model','aponocontEXT_model')
# 
# gett<-c('specialgate')#,'aponocontEXT_model')

allmodels<-gett
# dir(path="./",pattern=gett,full.names=F,recursive=F)
model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


# allmodels<-c(allmodels,'reference_model')


###multiple models to compare --> select some parameters, e.g.: TfhNum=300; p_Self=0.03; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'
select_pself<-'p_Self=0.04'
select_tfhtfr<-c('1:1')
i<-1

take <- which(filestobind%in%c('freq_zone_tfr.csv','depth.csv','exitpoint.csv','zpos_freq.csv','perc_free_tfr.csv'))

# take<-c(5,7,28)

for (i in take) {
print(filestobind[i])

##list all files with path that are to bind, i.e. in different models folder collect all files with the same name
files <- list.files(path = allmodels,
  pattern = paste0("^",filestobind[i],"$"),
                 recursive = TRUE, full.names = TRUE)
#small check:
if (length(files)!=length(allmodels)) {
    print(paste("TAKING MORE FILES!!stop! nfiles, ",length(files)," allmod:",length(allmodels)));
    stop;}
##the following will be the name of the data frame containing the info
cc_nn <- unlist(strsplit(filestobind[i],'.csv'))
##this is the data frame collecting all results
cc_file <- rbindlist(lapply(files, fread), fill = TRUE)

# print(cc_file$Tmove)

cc_file <- cc_file[which(cc_file$TfhNum%in%select_tfh & cc_file$SelfMode%in%select_selfmode & cc_file$SelfPerc%in%select_pself & cc_file$Tfh_by_Tfr%in%select_tfhtfr),]

cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)
cc_file$model_move[which(cc_file$TfrModel=='Reference')] <- 'Reference'

cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind[take]) ##these are the names of final dataframes (without extension)

pdf('playpos_plot.pdf',h=10,w=32)




##depth

ss<-as.numeric(unlist(strsplit((depth$variable),'d.') ))
ss<-ss[!is.na(ss)]
depth$variable<-ss
depth$CellType<-factor(depth$CellType, c('ASC','CD138','Tfr'))

depth_plot <- ggplot(depth[which(depth$time==240),], 
                    aes(x=variable, col=CellType, fill=CellType, shape=CellType)) +
                    geom_col(aes(y=value),alpha=.3)+
#                     geom_bar(aes(y=value),stat='identity',position=position_dodge(),alpha=.5)+
#                     geom_point(aes(y=value))+
#                     geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.3)+
#                     geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean))+
            xlab('Distance from the GC-centre') + ylab('Freq.') +
            facet_wrap(~time+Tmove,scale='free_y',ncol=2)+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')

print(depth_plot)

zpos_freq<-zpos_freq[which(zpos_freq$SimN==1)]
idtomeltby<-colnames(zpos_freq)[-which(colnames(zpos_freq)%in% c('CB', 'CC', 'Tfh', 'Tfr', 'CD138','ASC'))]

zzpos<-melt(zpos_freq, id.vars=idtomeltby)

names(zzpos)[names(zzpos)=='variable']<-'CellType'

zzpos_plot <- ggplot(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),],
                    aes(x=Zpos, col=CellType, fill=CellType, shape=CellType)) +
                    geom_col(aes(y=value),alpha=.3)+
#                     geom_bar(aes(y=value),stat='identity',position=position_dodge(),alpha=.5)+
#                     geom_point(aes(y=value))+
#                     geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.3)+
#                     geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean))+
            geom_vline(aes(xintercept=32),col='black',linetype=2)+
            xlab('Z-position in the GC') + ylab('Freq.') +
#             annotate( geom='text', x=1:nrow(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),]), y=0, label='DZ', vjust=3.5 )+
            facet_wrap(~time+Tmove,scale='free_y',ncol=2)+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')

print(zzpos_plot)


aaa <- zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CD138','ASC','Tfr') ),]
aaa$CellType<-factor(aaa$CellType, c('ASC','CD138','Tfr'))

zzpos_plot1 <- ggplot(aaa[which(aaa$time==240),], 
                    aes(x=Zpos, col=CellType, fill=CellType, shape=CellType)) +
                    geom_col(aes(y=value),alpha=.3)+
            geom_vline(aes(xintercept=32),col='black',linetype=2)+
            xlab('Z-position in the GC') + ylab('Freq.') +
#             annotate( geom='text', x=1:nrow(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),]), y=0, label='DZ', vjust=3.5 )+
            facet_wrap(~time+Tmove,scale='free_y',ncol=2)+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')

print(zzpos_plot1)

jw<-.1
dw<-1

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)

### 
##inLZ = n.LZTfr/n.Tfr
##inDZ = n.DZTfr/n.Tfr
##fr_LZDZ = n.LZTfr/n.DZTfr


tomm<-colnames(freq_zone_tfr)[-c(2:4)]
freq_zone_tfr <- melt(freq_zone_tfr,id.vars=tomm)

freqzone_plot <- ggplot(freq_zone_tfr[freq_zone_tfr$time%in%c(24,72,240),], 
                    aes(x=TfrModel, col=variable, fill=variable, shape=as.factor(round(time)))) +
                    geom_point(aes(y=value),position=posn.jd,size=3)+
                    
            xlab('TfrModel') + ylab('Freq. in zone') +
#             annotate( geom='text', x=1:nrow(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),]), y=0, label='DZ', vjust=3.5 )+
            facet_wrap(~Tmove,ncol=2,scale='free_y')+
            theme_prism(base_fontface='plain',base_size = 20
            ,base_line_size=.5
            )+
            theme( prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')

print(freqzone_plot)

AAA<-ggarrange( depth_plot,zzpos_plot1,zzpos_plot,freqzone_plot, labels='AUTO',font.label = list(size = 20))
# # %%Fig.S2
ggsave('FigS2.png', h=9,w=18, dpi=300)


freqzone_plot <- ggplot(freq_zone_tfr[freq_zone_tfr$variable=='fr_LZDZ',], 
                    aes(x=round(time/24), col=Tmove, fill=Tmove,linetype=TfrNum)) +
                    geom_point(aes(y=value),size=3)+
                    geom_line(aes(y=value))+
            xlab('TfrModel') + ylab('Freq. in zone') +
#             annotate( geom='text', x=1:nrow(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),]), y=0, label='DZ', vjust=3.5 )+
#             facet_wrap(~Tmove,ncol=2,scale='free_y')+
theme_bw()#+
#             theme_prism(base_fontface='plain',base_size = 20
#             ,base_line_size=.5
#             )+
#             theme( prism.ticks.length.y=unit(-3,'pt'))+
#             scale_y_continuous(guide='prism_minor')

print(freqzone_plot)

dev.off()



