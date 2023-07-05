

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
cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))
##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind[take]) ##these are the names of final dataframes (without extension)


zpos_freq<-zpos_freq[which(zpos_freq$SimN==1)]
idtomeltby<-colnames(zpos_freq)[-which(colnames(zpos_freq)%in% c('CB', 'CC', 'Tfh', 'Tfr', 'CD138','ASC'))]

zzpos<-melt(zpos_freq, id.vars=idtomeltby)

names(zzpos)[names(zzpos)=='variable']<-'CellType'

levels(zzpos$CellType)[levels(zzpos$CellType)=='CD138'] <- 'GCPC'

COLORS <- c('purple','blue','red','olivedrab','cyan')

zzpos_plot <- ggplot(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('CB','CC','Tfh','Tfr') ),],
                    aes(x=Zpos, col=CellType, fill=CellType)) +
                    geom_col(aes(y=value),alpha=.3)+
#                     geom_bar(aes(y=value),stat='identity',position=position_dodge(),alpha=.5)+
#                     geom_point(aes(y=value))+
#                     geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.3)+
#                     geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean))+
            geom_vline(aes(xintercept=32),col='black',linetype=2)+
            xlab('Z-position in the GC') + ylab('Freq.') +
            scale_fill_manual(values=COLORS) +
            scale_colour_manual(values=COLORS)+
            facet_wrap(~Tmove,scale='free_y',ncol=2)+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme(prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')

zzpos_plot2 <- ggplot(zzpos[which(zzpos$time==240 & zzpos$CellType%in%c('GCPC','Tfh','Tfr') ),],
                    aes(x=Zpos, col=CellType, fill=CellType)) +
                    geom_col(aes(y=value),alpha=.3)+
#                     geom_bar(aes(y=value),stat='identity',position=position_dodge(),alpha=.5)+
#                     geom_point(aes(y=value))+
#                     geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.3)+
#                     geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean))+
            geom_vline(aes(xintercept=32),col='black',linetype=2)+
            xlab('Z-position in the GC') + ylab('Freq.') +
            scale_fill_manual(values=COLORS[3:5]) +
            scale_colour_manual(values=COLORS[3:5])+
            facet_wrap(~Tmove,scale='free_y',ncol=2)+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme(prism.ticks.length.y=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')



            
AAA<-ggarrange( zzpos_plot,zzpos_plot2, labels=c('B','C'),font.label = list(size = 20),ncol=1)
# # %%Fig.S2
ggsave('FigS1B.png', h=10,w=12, dpi=300)



