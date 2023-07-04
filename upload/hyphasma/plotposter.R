# ggbreak v0.1.1
# 
# If you use ggbreak in published research, please cite the following
# paper:
# 
# S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively
# utilize plotting space to deal with large datasets and outliers.
# Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846


# find . -name 'trackfate.out' -exec rm -r {} \;

# ggsave("inGCSelf_compare.png", width = 40, height = 20, units = "cm", bg="transparent", dpi = 600)

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

results_folder<-getwd()#'/home/msi18/Desktop/diversity_plus_DZapop/hyphasma/results/'

setwd(results_folder)

REFFOLD<-'../REFERENCE'


folder_param = list.dirs(path = ".", recursive = FALSE)

###plot function

# https://community.rstudio.com/t/how-to-automatically-add-text-annotations-or-tags-outside-of-faceted-plots/13700/3
tag_facet2 <-  function(p, open="", close = "",
         tag_pool = letters,
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



folders<-c('Apoptosis','Trogocytosis_NW','il10sig','noplasmaonly','il10general')

######
outnumall<-outnum<-GC_refout<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    outnumall<-rbind(outnumall,data.frame(read.csv(file="outnumall.csv", header=T, sep=",")))
    outnum <-rbind(outnum,data.frame(read.csv(file='outnum.csv',header=T,sep=',')))
        if (count==1){
        GC_refout<-data.frame(read.csv(file="GC_refout.csv", header=T, sep=","))
        }
    setwd(results_folder)
}
selday<-21*24

mode.labs<-c("Flat", "Confined", "Ebb")
names(mode.labs)<-unique(outnum$SelfMode)




jw<-.05
dw<-.3

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)
posn.d <- position_dodge(width = dw)

selfout_plot <- ggplot(outnum[which(outnum$time%in%selday & outnum$Mutation=='Total'),],
            aes(x=(SelfMode),y=mean, linetype=TfhNum, fill=TfrNum
#             , linetype=Type,shape=Type
            )) +
            geom_point(outnumall[which(outnumall$time%in%selday & outnumall$Mutation=='Total'),],mapping=aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.jd,alpha=.5) +
            geom_errorbar(aes(ymin=mean-std,ymax=mean+std, shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),width=.1,position = posn.d,col='black') +
            geom_point(aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.d,size=2) +
            geom_errorbar(GC_refout[which(GC_refout$time%in%selday & GC_refout$Mutation=='Total'),],mapping=aes(ymin=mean-std,ymax=mean+std),width=.1,position = posn.d,shape=1,col='black') +
            geom_point(GC_refout[which(GC_refout$time%in%selday & GC_refout$Mutation=='Total'),], mapping=aes(x=SelfMode,y=mean),shape=23,col='black',position = posn.d)+
            
#             geom_boxplot(position=posn.d,fill=NA)+
#             facet_wrap(~Mutation)+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
#                 scale_y_cut(breaks=c(300,500,3000),which=c(1,2,3,4),scales=c(1,1,0,1)) +
                
                
                scale_x_discrete(labels = c("Flat", "Confined", "Ebb"))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
                
                xlab('Mode of generating Self-spec. clone')+
                ylab('# ASCs (d. 21)')

selfmode.labs<-c("Flat", "Confined", "Ebb")
names(selfmode.labs)<-unique(outnum$SelfMode)

outnumref<-rbind(outnum,GC_refout)

selfout_plot <- ggplot(outnum[which(outnum$time%in%selday & outnum$Mutation=='Total'),],
            aes(x=(TfrModel),y=mean, linetype=TfhNum, fill=TfrNum
#             , linetype=Type,shape=Type
            )) +
            geom_point(outnumall[which(outnumall$time%in%selday & outnumall$Mutation=='Total'),],mapping=aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.jd,alpha=.5) +
            geom_errorbar(aes(ymin=mean-std,ymax=mean+std, shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),width=.1,position = posn.d,col='black') +
            geom_point(aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.d,size=2) +
            geom_errorbar(GC_refout[which(GC_refout$time%in%selday & GC_refout$Mutation=='Total'),],mapping=aes(ymin=mean-std,ymax=mean+std),width=.1,position = posn.d,shape=1,col='black') +
            geom_point(GC_refout[which(GC_refout$time%in%selday & GC_refout$Mutation=='Total'),], mapping=aes(x=TfrModel,y=mean),shape=23,col='black',position = posn.d)+
            
#             geom_boxplot(position=posn.d,fill=NA)+
            facet_wrap(~SelfMode,labeller=labeller(SelfMode=selfmode.labs))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
#                 scale_y_cut(breaks=c(300,500,3000),which=c(1,2,3,4),scales=c(1,1,0,1)) +
                
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
                
                xlab('Tfr Models')+
                ylab('# ASCs (d. 21)')

ggsave("selfout_plot.png", width = 35, height = 12, units = "cm", bg="transparent", dpi = 600)


####

######
p_selfout<-ALLp_selfout<-p_refout<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    p_selfout<-rbind(p_selfout,data.frame(read.csv(file="p_selfout.csv", header=T, sep=",")))
    ALLp_selfout <-rbind(ALLp_selfout,data.frame(read.csv(file="ALLp_selfout.csv", header=T, sep=",")))
        if (count==1){
        p_refout<-data.frame(read.csv(file="p_refout.csv", header=T, sep=","))
        }
    setwd(results_folder)
}

selday<-max(p_selfout$time)

dayp_selfout<-p_selfout[which(p_selfout$time%in%selday),]
ALLdayp_selfout<-ALLp_selfout[which(ALLp_selfout$time%in%selday),]

p_refout<-p_refout[which(p_refout$time%in%selday),]

jw<-.05
dw<-.5

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)
posn.d <- position_dodge(width = dw)

perc_selfout_plot <- ggplot(dayp_selfout,#[which(dayp_selfout$TfhNum=='Tfh.200'),],
            aes(x=(SelfMode),y=mean, linetype=TfhNum, fill=TfrNum
#             , linetype=Type,shape=Type
            )) +
            geom_point(ALLdayp_selfout,mapping=aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.jd,alpha=.5) +
            geom_errorbar(aes(ymin=mean-std,ymax=mean+std, shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),width=.1,position = posn.d,col='black') +
            geom_point(aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.d,size=2) +
            geom_errorbar(p_refout,mapping=aes(ymin=mean-std,ymax=mean+std),width=.1,position = posn.d,shape=1) +
            geom_point(p_refout, mapping=aes(x=SelfMode,y=mean),shape=23,col='black',position = posn.d)+
            
#             geom_boxplot(position=posn.d,fill=NA)+
#             facet_wrap(~SelfPerc+TfhNum,ncol=length(spsel))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
                scale_y_cut(breaks=c(0.05,0.135)) +
                
                
                scale_x_discrete(labels = c("Flat", "Confined", "Ebb"))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
                
#                 scale_fill_discrete(name='N. T cells', limits = c("Tfh.200","Tfh.400","Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('200 Tfh; 0 Tfr','400 Tfh, 0 Tfr','40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. T cells', limits = c("Tfh.200","Tfh.400","Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('200 Tfh; 0 Tfr','400 Tfh, 0 Tfr','40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='Tfh:Tfr N.')+
                
                
                xlab('Mode of generating Self-spec. clones')+
                ylab('Self ASCs (% ASCs) (d. 21)')


selfmode.labs<-c("Flat", "Confined", "Ebb")
names(selfmode.labs)<-unique(dayp_selfout$SelfMode)

perc_selfout_plot <- ggplot(dayp_selfout,#[which(dayp_selfout$TfhNum=='Tfh.200'),],
            aes(x=(TfrModel),y=mean, linetype=TfhNum, fill=TfrNum
#             , linetype=Type,shape=Type
            )) +
            geom_point(ALLdayp_selfout,mapping=aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.jd,alpha=.5) +
            geom_errorbar(aes(ymin=mean-std,ymax=mean+std, shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),width=.1,position = posn.d,col='black') +
            geom_point(aes(shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum),position = posn.d,size=2) +
            geom_errorbar(p_refout,mapping=aes(ymin=mean-std,ymax=mean+std),width=.1,position = posn.d,shape=1) +
            geom_point(p_refout, mapping=aes(x=TfrModel,y=mean),shape=23,col='black',position = posn.d)+
            
#             geom_boxplot(position=posn.d,fill=NA)+
            facet_wrap(~SelfMode,labeller=labeller(SelfMode=selfmode.labs))+
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
                scale_y_cut(breaks=c(0.05,0.2),which=c(3,2,1),scales=c(1.5,1,0.5),space=0) +
                
                
#                 scale_x_discrete(labels = c("Flat", "Confined", "Ebb"))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
                
#                 scale_fill_discrete(name='N. T cells', limits = c("Tfh.200","Tfh.400","Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('200 Tfh; 0 Tfr','400 Tfh, 0 Tfr','40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. T cells', limits = c("Tfh.200","Tfh.400","Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('200 Tfh; 0 Tfr','400 Tfh, 0 Tfr','40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='Tfh:Tfr N.')+
                
                
                xlab('Tfr Models')+
                ylab('Self ASCs (% ASCs) (d. 21)')

ggsave("percselfout_plot.png", width = 35, height = 12, units = "cm", bg="transparent", dpi = 600)



sopna<-perc_selfout_plot+theme(legend.position='none',axis.title.x = element_blank())#+xlab('')
nout<-selfout_plot+theme(legend.position='right',legend.box='vertical',axis.title.x = element_blank())#+xlab('')

arrout <- ggarrange(sopna,nout, labels=c("a","b"),legend.grob=get_legend(nout),legend='right')

comb_fig<-annotate_figure(arrout, bottom=textGrob('Mode of generating Self-spec. clone'))

ggsave("combout_fig.png", width = 50, height = 15, units = "cm", bg="transparent", dpi = 600)


###########

aff_all<-aff_refout<-aff_mean<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    aff_all<-rbind(aff_all,data.frame(read.csv(file="aff_all.csv", header=T, sep=",")))
    aff_mean <-rbind(aff_mean,data.frame(read.csv(file="aff_mean.csv", header=T, sep=",")))
        if (count==1){
        aff_refout<-data.frame(read.csv(file="aff_refout.csv", header=T, sep=","))
        }
    setwd(results_folder)
}


aff_refout<-aff_refout[which(aff_refout$time%in%selday),]



jw<-.05
dw<-.7

posn.jd <- position_jitterdodge(jitter.width =jw, dodge.width = dw)
posn.d <- position_dodge(width = dw)


# dayaff<-aff_all[which(aff_all$time%%selday<1e-15),]
dayaff<-aff_all[which(aff_all$time%in%selday),]

# aff_all<-aff_all[which(aff_all$SimN==1),]

# 'soutself'//'soutot'//'soutNONself'//'liveBC'
onlysome<-c('Self','Foreign')
spsel<-c('p_Self=0.1','No Inh; p_Self=0.1')

selfmode.labs<-c("Flat", "Confined", "Ebb")
names(selfmode.labs)<-unique(dayaff$SelfMode)

tfhnum.labs<-c("200 Tfh", "400 Tfh")
names(tfhnum.labs)<-unique(dayaff$TfhNum)

aff_refout<-aff_refout[which(aff_refout$Type%in%onlysome),]


aff_all_plot <- ggplot(dayaff[which(dayaff$Type%in%onlysome),],
            aes(x=(TfrModel),y=mean,linetype=Type,shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum)) +
            
            geom_boxplot(aff_refout,mapping= aes(x=(TfrModel),y=mean),col='black',outlier.shape=NA,width=.1,alpha=.7,position=position_dodge(width = .5),shape=21)+
            
            
            geom_point(position = posn.jd,size=1.5) +
            geom_boxplot(position=posn.d,fill=NA,outlier.shape=NA,width=.3)+
            facet_wrap(~TfhNum+SelfMode,ncol=3,labeller=labeller(SelfMode=selfmode.labs, TfhNum=tfhnum.labs))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
#                 ggtitle('Affinity of Autor. Abs-prod. cells; mean and sd per run')+
#                 xlab('Prob. of Self-spec. mutation')+
#                 guides(colour='none')+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='Specificity')+
                
#                 scale_shape_discrete(name='Tfh:Tfr N.')+
#                 scale_fill_discrete(name='N. T cells', limits = c("Tfh.200","Tfh.400","Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('200 Tfh; 0 Tfr','400 Tfh, 0 Tfr','40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
                
                
#                 scale_colour_discrete(name='N. Tfr', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40','80','200','400'))+
#                 scale_x_discrete(labels = c("Flat", "Confined", "Ebb"))+
                xlab('Tfr Models')+
                ylab('Affinity')


ggsave("aff_all_plot.png", width = 40, height = 18, units = "cm", bg="transparent", dpi = 600)



# //    time(1): sout--Self[m:sd](2:3): totalout[m:sd](4:5): nonselfout[m:sd](6:7)
# //      liveBCs[m:sd](8:9)


aff_time_out <- ggplot(aff_all[which(aff_all$Type=='liveBC' & aff_all$time%%24<1e-15),],#[which(aff_mean$TfhNum=='Tfh.200'),],
            aes(x=time/24,y=mean,linetype=TfhNum,shape=Tfh_by_Tfr, fill=TfrNum, col=TfrNum)) +
            geom_line()+
#             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1) +
            geom_point(position = posn.jd,size=1.5,alpha=.3) +
#             geom_boxplot(aes(x=as.factor(time/24)),position=posn.d,fill=NA,outlier.shape=NA,width=.3)+
            facet_wrap(~SelfMode+TfrModel,ncol=2)+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
#                 ggtitle('Affinity of Autor. Abs-prod. cells; mean and sd per run')+
#                 xlab('Prob. of Self-spec. mutation')+
#                 scale_fill_discrete(name='pSelf', labels=c('0.1','0.3','0.4','0.5'))+scale_color_discrete(name='pSelf', labels=c('0.1','0.3','0.4','0.5'))+
#                 scale_line_discrete(name='Specificity')+
#                 scale_x_discrete(labels = c("Flat", "Confined", "Ebb"))+
                xlab('Days post GC onset')+
                ylab('Affinity')

                

sopna<-perc_selfout_plot+theme(legend.position='none',axis.title.x = element_blank())#+xlab('')
affpna<-aff_all_plot+theme(legend.position='right',legend.box='vertical',axis.title.x = element_blank())#+xlab('')

arrout <- ggarrange(sopna,affpna, labels=c("a","b"),legend.grob=get_legend(affpna),legend='right')

comb_fig<-annotate_figure(arrout, bottom=textGrob('Mode of generating Self-spec. clone'))

ggsave("comb_fig.png", width = 40, height = 15, units = "cm", bg="transparent", dpi = 600)



###############


pselfGC<-pselfGC_refout<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    pselfGC<-rbind(pselfGC,data.frame(read.csv(file="pselfGC.csv", header=T, sep=",")))
#     aff_mean <-rbind(aff_mean,data.frame(read.csv(file="aff_mean.csv", header=T, sep=",")))
        if (count==1){
        pselfGC_refout<-data.frame(read.csv(file="pselfGC_refout.csv", header=T, sep=","))
        }
    setwd(results_folder)
}



selfmode.labs<-c("Flat", "Confined", "Ebb")
names(selfmode.labs)<-unique(pselfGC$SelfMode)

selfgc_plot1 <- ggplot(pselfGC,#[which(pselfGC$TfhNum=='Tfh.200'),],
            aes(x=(time/24),y=mean, fill=TfrNum, col=TfrNum,linetype=Tfh_by_Tfr, shape=TfrModel
#             , linetype=Type,shape=Type
            )) +
            geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,col=NA) +
            geom_line()+
            geom_point(pselfGC[which(pselfGC$time%%24<1e-15),],mapping=aes(x=time/24,y=mean, shape=TfrModel),size=1.5) +
            
            geom_line(pselfGC_refout, mapping=aes(x=time/24,y=mean),linetype=1)+
            geom_ribbon(pselfGC_refout, mapping=aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1,linetype=1) +
#             geom_boxplot(position=posn.d,fill=NA)+
            facet_wrap(~TfhNum+SelfMode,labeller=labeller(SelfMode=selfmode.labs))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_linetype_discrete(name='Tfh:Tfr N.')+
#                 scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
                
                xlab('Days p. GC onset')+
                ylab('Self GCBCs (% GCBCs)')


selfgc_plot<-tag_facet2(selfgc_plot1)


ggsave("selfgc_plot.png", width = 35, height = 18, units = "cm", bg="transparent", dpi = 600)


##############

GCnum<-GC_refnum<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    GCnum<-rbind(GCnum,data.frame(read.csv(file="GCnum.csv", header=T, sep=",")))
#     aff_mean <-rbind(aff_mean,data.frame(read.csv(file="aff_mean.csv", header=T, sep=",")))
        if (count==1){
        GC_refnum<-data.frame(read.csv(file="GC_refnum.csv", header=T, sep=","))
        }
    setwd(results_folder)
}


mode.labs<-c("Flat", "Confined", "Ebb")
names(mode.labs)<-unique(GCnum$SelfMode)

gc_plot <- ggplot(GCnum,#[which(GCnum$SelfMode=='SelfMode.0'),],
            aes(x=(time/24),y=mean, fill=TfrNum, col=TfrNum,linetype=Tfh_by_Tfr, shape=TfrModel
#             , linetype=Type,shape=Type
            )) +
            geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,col=NA) +
            geom_line()+
            geom_point(GCnum[which(GCnum$time%%24<1e-15),],
            mapping=aes(x=time/24,y=mean, shape=TfrModel),size=1.5) +
            
            geom_line(GC_refnum, mapping=aes(x=(time/24),y=mean),linetype=1) +
            geom_ribbon(GC_refnum, mapping=aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1,linetype=1) +
#             geom_point(GC_refout[which(GC_refout$time%%24<1e-15),], mapping=aes(x=(time/24),y=mean)) +
#             geom_boxplot(position=posn.d,fill=NA)+
            facet_wrap(~TfhNum+SelfMode,labeller=labeller(SelfMode=selfmode.labs))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='Tfh:Tfr N.')+
                scale_linetype_discrete(name='Tfh:Tfr N.')+
                
                xlab('Days p. GC onset')+
                ylab('# GCBCs')


psgc <- selfgc_plot1+theme(legend.position='none',axis.title.x = element_blank())#+xlab('')
ngc <- gc_plot+theme(legend.position='right',legend.box='vertical',axis.title.x = element_blank())#+xlab('')

arrgc <- ggarrange(ngc,psgc, labels=c("a","b"),legend.grob=get_legend(ngc),legend='right',ncol=1)

combgc_fig<-annotate_figure(arrgc, bottom=textGrob('Days p. GC onset'))

ggsave("combgc_fig.png", width = 35, height =30, units = "cm", bg="transparent", dpi = 600)

##########

apo_self<-apo_self_ref<-NULL
count<-0
for (ff in folders) {
    count<-count+1
    setwd(ff)
    setwd('results_summary')
    print(ff)
    apo_self<-rbind(apo_self,data.frame(read.csv(file="apo_self.csv", header=T, sep=",")))
#     aff_mean <-rbind(aff_mean,data.frame(read.csv(file="aff_mean.csv", header=T, sep=",")))
        if (count==1){
        apo_self_ref<-data.frame(read.csv(file="apo_self_ref.csv", header=T, sep=","))
        }
    setwd(results_folder)
}

selfmode.labs<-c("Flat", "Confined", "Ebb")
names(selfmode.labs)<-unique(apo_self$SelfMode)

apogc_plot <- ggplot(apo_self,#[which(pselfGC$TfhNum=='Tfh.200'),],
            aes(x=(time/24),y=mean, fill=TfrNum, col=TfrNum,linetype=Tfh_by_Tfr, shape=TfrModel
#             , linetype=Type,shape=Type
            )) +
            geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,col=NA) +
            geom_line()+
            geom_point(apo_self[which(apo_self$time%%24<1e-15),],mapping=aes(x=time/24,y=mean, shape=TfrModel),size=1.5) +
            
            geom_line(apo_self_ref, mapping=aes(x=time/24,y=mean))+
            geom_ribbon(apo_self_ref, mapping=aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1) +
#             geom_boxplot(position=posn.d,fill=NA)+
            facet_wrap(~TfhNum+SelfMode,labeller=labeller(SelfMode=selfmode.labs))+#,scale='free_x') +
            theme_bw()+
                theme(strip.background=element_rect(fill='white'))+
                
                scale_fill_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_colour_discrete(name='N.Tfr', limits = c("Tfr.0", "Tfr.40", "Tfr.80", "Tfr.100", "Tfr.200",'Tfr.400'), labels=c('0 Tfr','40 Tfr','80 Tfr','100 Tfr' ,'200 Tfr','400 Tfr'))+
                scale_linetype_discrete(name='Tfh:Tfr N.')+
#                 scale_linetype_discrete(name='N.Tfh', limits = c("Tfh.200",'Tfh.400'), labels=c('200 Tfh','400 Tfh'))+
#             geom_line()+
#             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.5,size=.1) +
#             geom_point(apo_self[which(apo_self$time%%24<1e-15),],mapping=aes(x=time/24,y=mean)) +
#             geom_ribbon(apo_self_ref, mapping=aes(ymin=mean-std,ymax=mean+std),alpha=.5,size=.1) +
# #             geom_boxplot(position=posn.d,fill=NA)+
#             facet_wrap(~SelfMode,labeller=labeller(SelfMode=selfmode.labs))+#,scale='free_x') +
#             theme_bw()+
#                 theme(strip.background=element_rect(fill='white'))+
#                 scale_fill_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='Tfh:Tfr N.', limits=c('1:1','5:1'))+
#                 scale_linetype_discrete(name='Tfh:Tfr N.', limits=c('1:1','5:1'))+
                xlab('Days p. GC onset')+
                ylab('Self Apop (% Apo)')



ggsave("apogc_plot.png", width = 30, height = 18, units = "cm", bg="transparent", dpi = 600)

###
# 
# pdf('TCoccupancy_and_ratio.pdf',h=4,w=13)
# 
# 
# tb_occup<-data.frame(read.csv(file="tb_occup.csv", header=T, sep=","))
# 
# 
# mode.labs<-c("Flat", "Confined", "Ebb")
# names(mode.labs)<-unique(tb_occup$SelfMode)
# 
# 
# tbocc_plot <- ggplot(tb_occup,#[which(GCnum$SelfMode=='SelfMode.0'),],
#             aes(x=(time/24), y=mean, fill=TfrNum, col=TfrNum, shape=Tfh_by_Tfr,linetype=Tfh_by_Tfr
# #             , linetype=Type,shape=Type
#             )) +
#             geom_line()+
#             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1,col=NA) +
#             geom_point(tb_occup[which(tb_occup$time%%24<1e-15),],
#             mapping=aes(x=time/24,y=mean),size=1.5) +
#             
#             facet_wrap(~Mutation+SelfMode,ncol=3,scale='free_y',labeller=labeller(SelfMode=mode.labs))+#,scale='free_x') +
#             theme_bw()+
#                 theme(strip.background=element_rect(fill='white'))+
#                 scale_fill_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 scale_linetype_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 
#                 xlab('Days p. GC onset')+
#                 ylab('Fraction of bound Tfr')
# 
# 
# tbocc_plot
# 
# 
# tfr_free<-data.frame(read.csv(file="tfr_free.csv", header=T, sep=","))
# 
# 
# 
# mode.labs<-c("Flat", "Confined", "Ebb")
# names(mode.labs)<-unique(tfr_free$SelfMode)
# 
# 
# tbfree_plot <- ggplot(tfr_free[which(tfr_free$time%%24<1e-15),],#[which(GCnum$SelfMode=='SelfMode.0'),],
#             aes(x=(time/24), y=mean, fill=TfrNum, col=TfrNum, shape=Tfh_by_Tfr,linetype=Tfh_by_Tfr
# #             , linetype=Type,shape=Type
#             )) +
#             geom_line()+
# #             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1,col=NA) +
#             geom_point(tfr_free[which(tfr_free$time%%24<1e-15),],
#             mapping=aes(x=time/24,y=mean),size=1.5) +
#             geom_errorbar(tfr_free[which(tfr_free$time%%24<1e-15),],
#             mapping=aes(ymin=mean-std,ymax=mean+std),width=.1)+
#             
#             facet_wrap(~SelfMode,ncol=3,scale='free_y',labeller=labeller(SelfMode=mode.labs))+#,scale='free_x') +
#             theme_bw()+
#                 theme(strip.background=element_rect(fill='white'))+
#                 scale_fill_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 scale_linetype_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 
#                 xlab('Days p. GC onset')+
#                 ylab('Fraction of free Tfr')
# 
# 
# tbfree_plot
# 
# 
# 
# tfh_free<-data.frame(read.csv(file="tfh_free.csv", header=T, sep=","))
# 
# 
# mode.labs<-c("Flat", "Confined", "Ebb")
# names(mode.labs)<-unique(tfh_free$SelfMode)
# 
# 
# tbfree_plot <- ggplot(tfh_free,#[which(GCnum$SelfMode=='SelfMode.0'),],
#             aes(x=(time/24), y=mean, fill=TfrNum, col=TfrNum, shape=Tfh_by_Tfr,linetype=Tfh_by_Tfr
# #             , linetype=Type,shape=Type
#             )) +
#             geom_line()+
# #             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.3,size=.1,col=NA) +
#             geom_point(tfh_free[which(tfh_free$time%%24<1e-15),],
#             mapping=aes(x=time/24,y=mean),size=1.5) +
#             geom_errorbar(tfh_free[which(tfh_free$time%%24<1e-15),],
#             mapping=aes(ymin=mean-std,ymax=mean+std),width=.1)+
#             
#             facet_wrap(~SelfMode,ncol=3,scale='free_y',labeller=labeller(SelfMode=mode.labs))+#,scale='free_x') +
#             theme_bw()+
#                 theme(strip.background=element_rect(fill='white'))+
#                 scale_fill_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_colour_discrete(name='N. Tfr cells', limits = c("Tfr.40", "Tfr.80", "Tfr.200",'Tfr.400'), labels=c('40 Tfr','80 Tfr','200 Tfr','400 Tfr'))+
#                 scale_shape_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 scale_linetype_discrete(name='N. Tfh:Tfr')+#, labels=c('200 Tfh','400 Tfh'))+
#                 
#                 xlab('Days p. GC onset')+
#                 ylab('Fraction of free Tfh')
# 
# 
# tbfree_plot
# 
# cellratios<-data.frame(read.csv(file="cellratios.csv", header=T, sep=","))
# 
# 
# mode.labs<-c("Flat", "Confined", "Ebb")
# names(mode.labs)<-unique(cellratios$SelfMode)
# 
# ratio_plot1 <- ggplot(cellratios,#[which(cellratios$time%%selday<1e-15),], 
# aes(x=time/24, fill=TfrNum, col=TfrNum, shape=TfhNum,linetype=Tfh_by_Tfr)) +
#             geom_line(aes(y=mean)) +
#             geom_ribbon(aes(ymin=mean-std,ymax=mean+std),alpha=.35,col=NA) +
#             geom_point(cellratios[which(cellratios$time%%24<1e-15),], mapping=aes(x=time/24,y=mean)) +
#             geom_errorbar(cellratios[which(cellratios$time%%24<1e-15),], mapping=aes(x=time/24,ymin=mean-std,ymax=mean+std),width=.05) +
#             facet_wrap(~SelfMode,ncol=3,scale='free_y',labeller=labeller(SelfMode=mode.labs))+
#             theme_bw() +
#             ylim(0,5.5) +
#             theme(strip.background=element_rect(fill='white'))+ggtitle('Ratios') + xlab('Days post GC onset') + ylab('cell ratio')
# 
# 
# ratio_plot1
# 
# dev.off()
# 
# 
