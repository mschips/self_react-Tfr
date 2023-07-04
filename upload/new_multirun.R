parent_folder<-getwd() ##where to run multihyphasma
param_folder<- paste0(parent_folder,'/hyphasma/parameter_files/to_run') ##path to param files

setwd(param_folder)

# setwd('/home/msi18/Desktop/self_redeem/hyphasma/parameter_files/to_run')
files <- list.files(pattern = "\\.par$")
setwd(parent_folder)

common<-readLines("runmulticore_settings")

  ncores<- as.numeric(as.character(unlist(strsplit(common[8],'='))[2]))
  nruns<- as.numeric(as.character(unlist(strsplit(common[10],'='))[2]))
  totalruns<-ncores*nruns
  
Nfolders <- 0:(ncores-1)
multiF<-NULL

for (nn in Nfolders) {
    if (nn<10) {
        multiF<-c(multiF,paste0('hyphasma_tmp0',nn))
    } else {multiF<- c(multiF, paste0('hyphasma_tmp',nn))}
}

i<-1
for (i in 1:length(files)){
  tx  <- readLines("runmulticore_settings")
  
  tx[12]=gsub(pattern='.par',replace='',paste0('parfile=',files[i]))
  writeLines(tx,"runmulticore_settings")
  system('bash runmulticore')
  Sys.sleep(180)
  
  curr_parfile<-readLines("runmulticore_settings")
  to_remove<-paste0('rm ',param_folder,'/',as.character(unlist(strsplit(curr_parfile[12],'='))[2]),'.par')
  curr_parfile<-paste0(as.character(unlist(strsplit(curr_parfile[12],'='))[2]),'-all')
  END_SIM <- 0 ##counter for total n of sim
  while(END_SIM<totalruns) {
    Sys.sleep(180)
    END_SIM<-0
    for (ff in multiF) { ##count how many sim completed in each folder
        pp<-paste0(ff,'/results/',curr_parfile)
        TMP <- length(list.dirs(path = pp , recursive = FALSE))
        END_SIM<-END_SIM+TMP
    }
    print(paste('Not Done: ',END_SIM,' completed of ',totalruns))
  }
  print(paste('End, nsim: ',END_SIM))
  Sys.sleep(5)
  print('Collecting')
  system('bash runmulticore_collect')
  Sys.sleep(10)
  system('rm -rf hyphasma_tmp*/')
  system(to_remove)
}
