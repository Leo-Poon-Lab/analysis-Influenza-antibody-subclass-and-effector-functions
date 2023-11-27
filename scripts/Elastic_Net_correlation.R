library(readxl)
library(tidyverse)

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)
library("pheatmap")
library(scico)
library(parallel)
library(FactoMineR)
library("factoextra")
library(mdatools) # https://mdatools.com/docs/plsda.html
library(writexl) 
library(corrr)

# read data
source(here::here("scripts/helper/deal_with_colnams.R"))
source(here::here("scripts/helper/prepare_multiplex_data.R"))
data_all <- prepare_data()
data_multiplex_aim2 <- data_all[[1]]
data_clinical_aim2 <- data_all[[2]]
data_multiplex_aim3 <- data_all[[3]]
data_clinical_aim3 <- data_all[[4]]

all_aims <- c("aim2", "aim3")
dir_rst <- "elastic_net_correlation_results"

for (this_aim in all_aims) { # hclust for log-transformed and scaled responses
	# this_aim=all_aims[2]
	if(this_aim=="aim2"){
		data_multiplex <- data_multiplex_aim2
		data_clinical <- data_clinical_aim2
	} else {
		data_multiplex <- data_multiplex_aim3
		data_clinical <- data_clinical_aim3
	}
	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)

	uniq_groups <- list(unique(data_multiplex$group)) # all groups

	for (this_group in uniq_groups) { # only one option here, i.e. for all groups
		# this_group <- uniq_groups[[1]]
		data_multiplex_subgroup <- data_multiplex %>% filter(group %in% this_group)
		data_meta <- data_clinical %>% filter(group%in%this_group)
		data_meta$sample <- as.character(data_meta$sample)

		uniq_timepoints <- as.list(unique(data_multiplex_subgroup$timepoint))
		if(this_aim=="aim3" & length(this_group)>1){
			uniq_timepoints <- list(2)
		}
		
		for (this_timepoint in uniq_timepoints){
			# this_timepoint=uniq_timepoints[1]
			cur_path <- paste0(here::here("results/"), dir_rst, "/", this_aim, "/", paste0(this_group, collapse = "-"), "/timepoint_", paste0(this_timepoint, collapse = "-"), "/")
			dir.create(cur_path, showWarnings=F, recursive=T)
			print(cur_path)

			########### EN #####################################

			data_raw <- data_multiplex_subgroup %>% filter(timepoint%in%this_timepoint)
			data_meta_reorder <- left_join(data_raw %>% dplyr::select(sample, timepoint, infection), data_meta, "sample")
			stopifnot(all(data_raw$sample==data_meta_reorder$sample))

			data_m <- data_raw %>% dplyr::select(-group, -sample, -timepoint, -infection) 
			stopifnot(length(apply(data_m, 2, function(x){which(is.na(x))}))==0) # make sure no missing values
			data_m <- data_m %>% mutate_all(as.numeric)
			data_m[data_m<0] <- 0
			data_m <- log10(data_m+1) # log transformation
			data_m <- apply(data_m, 2, pcaMethods::prep, scale="uv") # scale the data, with unit variance

			data_combined <- bind_cols(data_m,response=as.numeric(data_meta_reorder$infection))
			sampling_factor <- 0.8
			n_repeats <- 2000
			n_obs <- nrow(data_combined)
			set.seed(2023)

			variable_counts <- mclapply(seq_len(n_repeats), function(i) {
				print(i)
				idx_choosen <- sample(n_obs, round(n_obs*sampling_factor), replace=FALSE)
				Y <- data_meta_reorder$infection[idx_choosen]
				X <- makeX(data.frame(data_m[idx_choosen,]), na.impute = T)
				colnames(X) <-  colnames(data_m[idx_choosen,])
				# 10-fold CV to find the optimal lambda 
				
				enet.cv=cv.glmnet(X,Y,alpha=0.5, type="deviance", family="binomial", nfolds=10)
				## Fit lasso model with 100 values for lambda
				enet_mdl = glmnet(X,Y,alpha=0.5,nlambda=100)
				## Extract coefficients at optimal lambda
				out <- coef(enet_mdl,s=enet.cv$lambda.min)
				out <- as.matrix(out)[,1]
				return(names(out)[out!=0])
			}, mc.cores=16)
			
			count_data <- sort(table(unlist(variable_counts)), decreasing=T)
			count_data <- count_data[-which(names(count_data)=="(Intercept)")] # remove intercerpt
			out_file <- paste0(cur_path, "ranking_from_elastic_correlation.xlsx")
			write_xlsx(tibble(varibale=names(count_data), times=count_data), out_file)
			
			threshold_for_important_var <- 0.7
			imp_var <- count_data[count_data>n_repeats*threshold_for_important_var]
			imp_var <- gsub("`", "", names(imp_var))
			
			# correlation network
			min_cor_value <- 0.8
			rst_corr <- data_combined %>% select(-response) %>% correlate()

			# any_over_min <- function(x) any(x > min_cor_value, na.rm = TRUE)
			# rst_corr_filter <- rst_corr
			# rst_corr_filter[rst_corr_filter<min_cor_value] <- NA
			# check_all_na <- apply(rst_corr_filter, 2, function(x){all(is.na(x))})
			# rst_corr_filter <- rst_corr_filter[,!check_all_na]
			
			check <- rst_corr %>% select(all_of(imp_var)) %>% apply(1, function(x) {any(abs(x)>=min_cor_value)})
			vars_focus <- rst_corr$term[check]
			vars_focus <- rst_corr$term[rst_corr$term %in% c(vars_focus, imp_var)]
			rst_corr_filter_i <- rst_corr %>% focus(vars_focus, mirror = TRUE)
			
			p_out <- rst_corr_filter_i %>% network_plot()
			out_file <- paste0(cur_path, "cor_network.pdf")
			ggsave(out_file, width=10, height=10, plot=p_out, device=cairo_pdf)

			writeLines(imp_var, paste0(cur_path, "imp_vars.txt"))

			p_out <- rst_corr_filter_i  %>% rplot() + theme(axis.text.x = element_text(angle = 30, hjust=1))
			out_file <- paste0(cur_path, "cor_plot.pdf")
			ggsave(out_file, width=10, height=10, plot=p_out, device=cairo_pdf)

		}

	}
	
}


