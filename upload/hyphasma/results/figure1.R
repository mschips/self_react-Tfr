
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

RESULTSET <- '^refere'
LEGENDUSE <- ''
# # 

######
allmodels<-dir(path="./",pattern=RESULTSET,full.names=F,recursive=F)

model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)

TAKE <- c('prob_mut.csv','perc_redeemed_gcbc.csv','perc_self_out.csv','perc_self_gcbc.csv','S_aff_gcbc.csv')

filestobind<-filestobind[which(filestobind%in%TAKE)]


###multiple models to compare --> select some parameters, e.g.: TfhNum=250; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'
select_pself<-c('p_Self=0.06','p_Self=0.02','p_Self=0.04')
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
cc_file <- cc_file[which(cc_file$TfhNum%in%select_tfh & cc_file$SelfMode%in%select_selfmode & cc_file$Tfh_by_Tfr%in%select_tfhtfr),]
# 
# ## when Tfr are absent there is no difference between movements -- need to double for plotting purpouse
# ref_mod <- cc_file[cc_file$TfrModel=='Reference',]
# ref_mod$Tmove <- 'R5'
# cc_file <- rbind(cc_file,ref_mod)
# #
# 
# cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)
# cc_file$model_move[which(cc_file$TfrModel=='Reference')] <- 'Reference'
# 
# if (LEGENDUSE!='Apoptosis') {
#     cc_file <-cc_file[cc_file$Tmove!='R5',]
# }

cc_file$SelfPerc <- vapply(cc_file$SelfPerc, function(x) unlist(strsplit( x,'=' ) )[2], 'a'  )


cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))

cc_file$Tfh_by_Tfr <- factor(cc_file$Tfh_by_Tfr, select_tfhtfr)

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind) ##these are the names of final dataframes (without extension)



###mutation prob
jw<-.07
dw<-.5

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)
posn.d <- position_dodge(width = dw)

ppmut<-prob_mut[which( prob_mut$Mut=='Foreign' & prob_mut$SelfPerc=='0.04' & prob_mut$SimN==1),-17]

ppmut$what<-'Mut. Pr.'
aaff<-S_aff_gcbc[which(S_aff_gcbc$SimN==1 & prob_mut$SelfPerc=='0.04'),]
aaff$what<-'Aff.'
full<-rbind(ppmut,aaff)

###mutation prob
prob_mutc_plot1 <- ggplot(full,
            aes(x=time/24,y=mean,col=what,shape=SelfPerc)) +
            geom_line()+
                    geom_point(full[which(full$time%%24<1e-15)],mapping=aes(y=mean,col=what,fill=what),shape=21,col='black')+
                    geom_errorbar(full[which(full$time%%24<1e-15)],mapping=aes(ymin=mean-sd, ymax=mean+sd,col=what,fill=what),width=.1)+
            xlab('Days p. GC onset') + ylab('values') +
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'),prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor')
            
###% selfGCBC
perc_self_gcbc_plot <- ggplot(perc_self_gcbc, 
            aes(x=time/24, y=mean, col=SelfPerc, fill=SelfPerc
            )) +
            
            geom_line()+
            geom_errorbar(perc_self_gcbc[which(perc_self_gcbc$time%%72<1e-15),], mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+
            geom_point(perc_self_gcbc[which(perc_self_gcbc$time%%72<1e-15),], mapping=aes(y=mean),shape=21,alpha=.8,col='black')+
            
            xlab('Days p. GC onset') + ylab('s-GCBC (% of GCBC)') +
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme( prism.ticks.length.y=unit(-3,'pt'),prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
            scale_y_continuous(guide='prism_minor')+ 
            scale_x_continuous(guide='prism_minor')

            

###% selfASC
perc_self_out_plot <- ggplot(perc_self_out,
aes(x=time/24, y=mean, col=SelfPerc, fill=SelfPerc)) +
geom_line()+
geom_errorbar(perc_self_out[which(perc_self_out$time%%72<1e-15),], mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+
geom_point(perc_self_out[which(perc_self_out$time%%72<1e-15),], mapping=aes(y=mean),shape=21,alpha=.8,col='black')+

xlab('Days p. GC onset') + ylab('s-ASC (% of ASC)') +
theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
theme( prism.ticks.length.y=unit(-3,'pt'),prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+#legend.title=element_text(),
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
scale_y_continuous(guide='prism_minor')+
scale_x_continuous(guide='prism_minor')

###

perc_redeemed_gcbc_plot <- ggplot(perc_redeemed_gcbc,
aes(x=time/24, y=mean, col=SelfPerc, fill=SelfPerc)) +
geom_line()+
geom_errorbar(perc_redeemed_gcbc[which(perc_redeemed_gcbc$time%%72<1e-15),], mapping=aes(ymin=mean-std, ymax=mean+std),width=.3)+
geom_point(perc_redeemed_gcbc[which(perc_redeemed_gcbc$time%%72<1e-15),], mapping=aes(y=mean),shape=21,alpha=.8,col='black')+

xlab('Days p. GC onset') + ylab('Red-GCBC (% of GCBC)') +
theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
theme( prism.ticks.length.y=unit(-3,'pt'),prism.ticks.length.x=unit(-3,'pt'), legend.text = element_text(size=20))+
scale_y_continuous(guide='prism_minor')+
            guides(colour=guide_legend(override.aes=list(linewidth=1.5,size=3.5)))+
scale_x_continuous(guide='prism_minor')


AAA<-ggarrange(prob_mutc_plot1, perc_self_gcbc_plot, perc_self_out_plot, perc_redeemed_gcbc_plot, labels='AUTO',font.label = list(size = 20))

ggsave('Figure1.png', h=10,w=15, dpi=300)
