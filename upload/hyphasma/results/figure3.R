
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

gett<-c('aponcnw_model','aponcr5_model','semigate.38nw_model','semigate.38r5_model','semigater5_model','semigatenw_model')
# 
# gett<-c('specialgate')#,'aponocontEXT_model')

allmodels<-gett
# dir(path="./",pattern=gett,full.names=F,recursive=F)
model <- unlist(strsplit(allmodels[1],'_model'))[1]

##get all the files to use for the plot -- each folder contains the same files so that only one can be used to store the names
pp<-paste0(HYPMODEL_folder,'/',allmodels[1],'/',model,'_folder')
filestobind <- list.files(path=pp, pattern='.csv',full.names=F,recursive=T)


# allmodels<-c(allmodels,'reference_model')


###multiple models to compare --> select some parameters, e.g.: TfhNum=300; p_Self=0.02; SelfMode=0
select_tfh<-c('Tfh.250')
select_selfmode<-'SelfMode.0'
select_pself<-'p_Self=0.04'
select_tfhtfr<-c('2:1')
i<-1

take <- which(filestobind%in%c('perc_free_tfr.csv'))

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


cc_file$Tmove<-vapply(cc_file$Tmove, function(x) if (x=='NW') x='Tfh-like' else x='CC-like','a')

cc_file$Tmove <- factor(cc_file$Tmove,c('Tfh-like','CC-like'))
cc_file$TfrModel<-vapply(cc_file$TfrModel, function(x) if (x=='ApoNC') x='Apoptosis' else x=x,'a')

cc_file$model_move <- paste0(cc_file$TfrModel,'-',cc_file$Tmove)

##assigning the right name to the dataframe
assign(cc_nn, cc_file)

}
print(filestobind[take]) ##these are the names of final dataframes (without extension)


free_tfr_plot <- ggplot(perc_free_tfr,
aes(x=time/24, col=Tmove, fill=Tmove)) +
geom_line(aes(y=mean))+
#geom_ribbon(aes(ymin=mean-std, ymax=mean+std),alpha=.1,col=NA)+
geom_errorbar(perc_free_tfr[perc_free_tfr$time%%24<1e-15,], mapping=aes(y=mean,ymin=mean-std, ymax=mean+std),width=.3)+
# geom_point(S_cumul_out[which(S_cumul_out$time%in%selecttime),],mapping=aes(y=mean),position=posn.jd)+

            facet_wrap(~TfrModel,scale='free_y')+
            theme_prism(base_fontface='plain',base_size = 20,base_line_size=.5)+
            theme(legend.title=element_text(), prism.ticks.length.y=unit(-3,'pt'), prism.ticks.length.x=unit(-3,'pt'))+
            scale_y_continuous(guide='prism_minor')+
            scale_x_continuous(guide='prism_minor')+#,minor_breaks=seq(0,21,2.5)
            #scale_fill_manual(values=cbbPalette) +
            #scale_colour_manual(values=cbbPalette)+
xlab('Days p. GC onset') + ylab('Fraction of free Tfr') +
            labs(fill='Tfr motility',color='Tfr motility')

#free_tfr_plot



tag_facet2 <-  function(p, open="", close = "",
         tag_pool = toupper(letters),
         x = 0, y = 0.5,
         hjust = 0, vjust = 0.5, 
         fontface = 2, ...){
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  nm <- names(gb$layout$facet$params$rows)
  
  tags <- paste0(open,tag_pool[unique(lay$COL)],close)
  
  tl <- lapply(tags, grid::textGrob, x=x, y=y,
               hjust=hjust, vjust=vjust, gp=grid::gpar(fontface=fontface))
  
  g <- ggplot_gtable(gb)
  g <- gtable::gtable_add_rows(g, grid::unit(1,"line"), pos = 0)
  lm <- unique(g$layout[grepl("panel",g$layout$name), "l"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t=1, l=lm)
  grid::grid.newpage()
  grid::grid.draw(g)
}


tag_facet2(free_tfr_plot)

  ggsave('Figure3.png', h=4,w=15, dpi=300)
