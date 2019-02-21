
# Compute loss aversion 
# (1)
# (2) 

rm(list = ls())
options(repos = c(CRAN = "http://cran.rstudio.com"))

## library ----------------------------------------------------------
name_pkg <- c(
  # Bayesian hierarchical model
  "hBayesDM",
  
  # Visualization
  "ggplot2",
  "gridExtra",
  
  # bind data files
  "readr"
)
name_pkg <- unique(name_pkg)

bool_nopkg <- !name_pkg %in% rownames(installed.packages())
if (any(bool_nopkg)) {
  install.packages(name_pkg[bool_nopkg])
}
invisible(lapply(name_pkg, library, character.only = T)) # load multiple packages

if (any(!name_pkg %in% rownames(installed.packages()))){
  stop(sprintf("Package(s) %s is(are) not available!\n", 
               paste0(name_pkg[!name_pkg %in% rownames(installed.packages())], collapse = ", ")
  )
  )
}

getVgamble <- function(gain, loss, rho, lambda){
  (gain^rho - lambda * loss^rho) / 2
}


## Working directory for the project ---------------------------------------
if(Sys.info()[1] == "Darwin"){
  # Mac
  rt_dir <- "/Users/PSO/Dropbox/Seminar/18_04.DANS/DANS_project/"
}else if(Sys.info()[1] == "Windows" & Sys.info()[7] == "PSO"){
  # Windows PC
  rt_dir <- "C:/Users/User/Dropbox/Seminar/18_04.DANS/DANS_project/"
}
setwd(rt_dir)

load(file = "./parameter_loss_aversion.RData")


## Save the behavioral data in the txt file --------------------------------
### make a list of files in behav_data_path
behav_data_path = "/Users/PSO/Desktop/dans2019-team3-master/behav_data/"  # path to behavioral data
file_list = dir(path = behav_data_path)
for(i in seq_along(file_list)){
  tmp_tsv <- read_tsv(file = paste0(behav_data_path, file_list[i]))
  tmp_tsv <-  data.frame(tmp_tsv,
                         rho = df_param[df_param$subj == as.numeric(substr(file_list[i], 5, 6)), "rho"],
                         lambda = df_param[df_param$subj == as.numeric(substr(file_list[i], 5, 6)), "lambda"]
  )
  tmp_tsv$V_gamble <- getVgamble(gain = tmp_tsv$gain, 
                                 loss = tmp_tsv$loss,
                                 rho = tmp_tsv$rho, 
                                 lambda = tmp_tsv$lambda)
  write_tsv(x = tmp_tsv, 
            path = paste0("./behav_data_with_param/", file_list[i]))
}
rm("tmp_tsv")
