
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

RESULTSET1 <- '^aponc'
LEGENDUSE1 <- 'Apoptosis'
# # 

RESULTSET2 <- '^semigate'
LEGENDUSE2 <- 'SemiGate.38'
# #

# RESULTSET <- '^semigate'
# RESRM <- '^semigate.3'
# LEGENDUSE <- 'SemiGate'
# #

# RESULTSET <- '^sgmutsel'
# LEGENDUSE <- '^SG+Mut'

######
allmodels1<-dir(path="./",pattern=RESULTSET1,full.names=F,recursive=F)
allmodels <- c(allmodels1,dir(path="./",pattern=RESULTSET2,full.names=F,recursive=F))


model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


allmodels<-c(allmodels,'reference_model') ##add results from '250:0' group


TAKE <- c('S_perc_self_gcbc.csv','S_perc_self_apo.csv','perc_self_gcbc.csv','perc_self_apo.csv')

filestobind<-filestobind[which(filestobind%in%TAKE)]


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'

select_pself<-c('p_Self=0.04')
select_tfhtfr<-c('2:1')

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


# ref_mod <- cc_file[cc_file$TfrModel=='Reference',]
# ref_mod$Tmove <- 'R5'
# cc_file <- rbind(cc_file,ref_mod)
#

cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)
cc_file$model_move[which(cc_file$TfrModel=='Reference')] <- 'Reference'

# if (LEGENDUSE!='Apoptosis') {
    cc_file <-cc_file[cc_file$Tmove!='R5',]
# }

cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))

cc_file$Tfh_by_Tfr <- factor(cc_file$Tfh_by_Tfr, select_tfhtfr)

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind) ##these are the names of final dataframes (without extension)



#######################data

S_perc_self_gcbc$Type<-'Live BCs'
S_perc_self_apo$Type <- 'Apoptotis'


perc_self_gcbc$Type<-'Live BCs'
perc_self_apo$Type <- 'Apoptotis'

Atog <- rbind(perc_self_gcbc,perc_self_apo)


{


Atog$TfrModel <- vapply(Atog$TfrModel, function(x) if (x=='ApoNC') x='Apoptosis' else x=x,'a')


ascaff_plot <- ggplot(Atog,
            aes(x=time/24,y=mean, col=Type, fill=Type)) +
		geom_line()+
            geom_errorbar(Atog[which(Atog$time%%72<1e-15),],mapping=aes(ymax=mean+std,ymin=mean-std),width=.3) +
           
            
            xlab('Days p. GC onset') + ylab('% s-BCRs') +
#             stat_compare_means(label=c('p.signif'),ref.group='Live BCs',method='wilcox.test',hide.ns=T,facet.by=list('time'))+
            
            facet_wrap(~TfrModel,scale='free')+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor')
            
ggsave('FigS8.png', h=4,w=15, dpi=300)

}
