## -------------------------------------------------------------------------
## Compute loss of aversion
##   (1) Toms et al.2007 : logistic model (use glm function)
##   (2) Prospect model : hierarchical Bayesian model (use hBayesDM package)
## -------------------------------------------------------------------------
rm(list = ls())
options(repos = c(CRAN = "http://cran.rstudio.com"))

## library ----------------------------------------------------------
name_pkg <- c(
  # Bayesian hierarchical model
  "hBayesDM",
  
  # Visualization
  "ggplot2",
  "gridExtra"
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

## Working directory for the project ---------------------------------------
if(Sys.info()[1] == "Darwin"){
  # Mac
  rt_dir <- "/Users/PSO/Dropbox/Seminar/18_04.DANS/DANS_project/"
}else if(Sys.info()[1] == "Windows" & Sys.info()[7] == "PSO"){
  # Windows PC
  rt_dir <- "C:/Users/User/Dropbox/Seminar/18_04.DANS/DANS_project/"
}
setwd(rt_dir)


## Save the behavioral data in txt format ----------------------------------
### make a list of files in behav_data_path
behav_data_path = "/Users/PSO/Dropbox/Seminar/18_04.DANS/dans2019-team3/behav_data/"  # path to behavioral data
file_list = dir(path = behav_data_path)

all_data = NULL 
### read each file and combine them into a single file
for (i in 1:length(file_list)) {
  ### read i_th data file
  tmpData = read.table(file.path(behav_data_path, file_list[i]), header=T, sep="\t")
  
  ### find its subject ID 
  ### From its file list, subject ID can be located as 5th~6th chracterer. Then convert it into an integer
  tmp_subjID = as.integer(substr(file_list[i], 5, 6))
  
  ### add subject ID to tmpData 
  tmpData$subjID = tmp_subjID
  
  all_data = rbind(all_data, tmpData)  
}
rm(ls()[grep(pattern = "tmp", x = ls())]) # remove temporary variables

### add 'cert' (certain outcome if subj don't gamble) column
all_data$cert = 0

### remove 'no response' trials (respcat == -1)
all_data = subset(all_data, respcat >= 0)  # select trials only if 'respcat' >= 0

### Copy 'respcat' to 'gamble' for hBayesDM. See ?hBayesDM::ra_noRA
all_data$gamble = all_data$respcat  # gamble=1 --> choose gamble, gamble=0 --> don't choose gamble

### check all_data
dim(all_data)  # it should be 3942x14
table(all_data$subjID)  # each subject should have 174-256 trials

# ### save all_data as a text file so that it can be loaded from hBayesDM
# write.table(all_data, "tom2007_behav.txt",
#             col.names = T, row.names = F, sep = "\t")


## (1) Toms et al.2007 -----------------------------------------------------
ids <- unique(all_data$subjID)
loss_av <- vector(mode = "numeric", length = length(ids))
for(i in seq_along(ids)){
  cat(sprintf("Subject ID=%d\nTable for the response variables is", ids[i]))
  print(table(all_data[all_data$subjID == ids[i],]$gamble))
  cat("\n")
  
  ### This replicates the results of the paper
  glmfit <- glm(gamble ~ loss + gain, 
                data = subset(all_data, subjID == ids[i]),
                family = "binomial")
  
  ### This gives different results from the paper
  # glmfit <- glm(gamble ~ parametric.loss + parametric.gain, 
  #               data = subset(all_data, subjID == ids[i]),
  #               family = "binomial")
  
  loss_av[i] <- -glmfit$coefficients[grepl("loss", names(glmfit$coefficients))] / glmfit$coefficients[grepl("gain", names(glmfit$coefficients))]
}
## warning is fine
range(loss_av) # same as in paper, range: 0.99 to 6.75


## (2) Prospect model ------------------------------------------------------
set.seed(6) # fix seed to make results reproducible

### fit "ra_prospect" model using hBayesDM packag
fit_prosp <- ra_prospect("tom2007_behav.txt", niter = 2000, nwarmup = 1000,
                         nchain = 4, ncore = 2, inits = "random")
# save(list = "fit_prosp", file = "./fit_prosp_hBayesDM.RData") # save the fitted object due to time 
load(file = "./fit_prosp_hBayesDM.RData")

### check if rhat is less than 1.1
rhat(fit_prosp, less = 1.1)

### extract parameters of interest
par_name <- c("rho", "lambda", "tau")
id_var <- lapply(par_name, function(x) grep(pattern = paste0(x, "\\["), 
                                            x = rownames(summary(fit_prosp$fit)[[1]])))
names(id_var) <- par_name

df_param <- data.frame(
  rho = summary(fit_prosp$fit)$summary[id_var$rho, "mean"],
  lambda = summary(fit_prosp$fit)$summary[id_var$lambda, "mean"],
  tau = summary(fit_prosp$fit)$summary[id_var$tau, "mean"]
)
df_param$subj <- 1:nrow(df_param)
df_param <- rbind(df_param, 
                  c(summary(fit_prosp$fit)$summary[paste0("mu_", par_name), "mean"], 0)
                  ) ### add group parameters
rownames(df_param) <- 1:nrow(df_param)

# save(list = c("df_param"), 
#      file = "./parameter_loss_aversion.RData")
write_tsv(x = df_param,
          path = "./parameter_loss_aversion.tsv"
)


## Draw figures ------------------------------------------------------------
### Figure options
list_pdf_option <- list()
list_pdf_option$width <- 6
list_pdf_option$height <- 4

user_theme <- theme(text = element_text(size=25),
                    axis.text = element_text(size = 20),
                    axis.title = element_text(size = 20),
                    legend.position = "top"
)

loss_av_bayes <- df_param$lambda

### boxplot
gg_fig1 <- ggplot(data = rbind(data.frame(type = "Toms", loss_av = loss_av),
                               data.frame(type = "hbayes", loss_av = loss_av_bayes)),
                  mapping = aes(x = type, y = loss_av, color = type)) + 
  geom_boxplot() +
  guides(color = FALSE) + 
  ggtitle("Boxplot of loss aversion") + 
  labs(y = "Loss aversion")

### scatterplot
gg_fig2 <- ggplot(data = data.frame(Toms = loss_av,
                                    hbayes = loss_av_bayes),
                  mapping = aes(x = hbayes, y = Toms)) + 
  geom_point() +
  geom_smooth(data = data.frame(Toms = loss_av[-which.max(loss_av)],
                                hbayes = loss_av_bayes[-which.max(loss_av)]),
              method = 'lm',formula = y ~ x, se = FALSE) + 
  ggtitle("Pairs of loss aversion")

### Save figures in pdf files
pdf(file = "prob3_fig.pdf", onefile = T, width = list_pdf_option$width, height = list_pdf_option$height)
grid.arrange(gg_fig1, gg_fig2, nrow = 1)
dev.off()
