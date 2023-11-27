library(tidyverse)
library(FactoMineR)
library("factoextra")
library(readxl)

# read data
source(here::here("scripts/helper/deal_with_colnams.R"))
source(here::here("scripts/helper/prepare_multiplex_data.R"))

data_all <- prepare_data()
data_multiplex_aim2 <- data_all[[1]]
data_clinical_aim2 <- data_all[[2]]
data_multiplex_aim3 <- data_all[[3]]
data_clinical_aim3 <- data_all[[4]]

all_aims <- c("aim2", "aim3")
all_plots <- c()
all_plots_names <- c()

for (this_aim in all_aims) { # PCA for raw responses
	# this_aim=all_aims[1]
	if(this_aim=="aim2"){
		data_multiplex <- data_multiplex_aim2
		data_clinical <- data_clinical_aim2
	} else {
		data_multiplex <- data_multiplex_aim3
		data_clinical <- data_clinical_aim3
	}

	uniq_groups <- as.list(unique(data_multiplex$group))
	uniq_groups <- c(list(unique(data_multiplex$group))) # one combinational datasets

	for (this_group in uniq_groups) {
		# this_group=uniq_groups[1]
		# this_group <- this_group[[1]]
		data_multiplex_subgroup <- data_multiplex %>% filter(group %in% this_group)
		data_meta <- data_clinical %>% filter(group%in%this_group)
		data_meta$sample <- as.character(data_meta$sample)

		uniq_timepoints <- as.list(unique(data_multiplex_subgroup$timepoint))
		if(this_aim=="aim3" & length(this_group)>1){
			uniq_timepoints <- list(2)
		}
		if(this_aim=="aim3" & length(this_group)==1 && this_group=="V1S1"){
			uniq_timepoints <- c(list(1:3), uniq_timepoints)
		}
		for (this_timepoint in uniq_timepoints){
			# this_timepoint=uniq_timepoints[1]
			# this_timepoint <- this_timepoint[[1]]
			cur_path <- paste0(here::here("results/MFA_results/"), this_aim, "/", paste0(this_group, collapse = "-"), "/timepoint_", paste0(this_timepoint, collapse = "-"), "/")
			dir.create(cur_path, showWarnings=F, recursive=T)
			print(cur_path)

			########### MFA #####################################

			data_raw <- data_multiplex_subgroup %>% filter(timepoint%in%this_timepoint)
			data_meta_reorder <- left_join(data_raw %>% select(sample, timepoint), data_meta, "sample")

      data_raw <- data_raw %>% mutate_at(vars(!any_of(c("group","sample","timepoint"))), as.numeric)
      func_log <- function(x){ # log transformation
        x[x<0] <- 0
        log10(x+1)
      }
      data_raw <- data_raw %>% mutate_at(vars(!any_of(c("group","sample","timepoint"))), func_log)
      id <- names(data_raw) %in% c("group","sample","timepoint")
      data_raw[,!id] <- pcaMethods::prep(data_raw[,!id], scale= "uv") # scale the data, with unit variance
      
      # including categorical variables
      data_raw <- left_join(data_raw, data_meta_reorder %>% select(sample, contains("pre infection HAI"), contains("post vaxx"), any_of(c("age", "male"))), "sample")
      
      ## GMT categories
			uniq_gmts <- names(data_raw)[grepl("GMT", toupper(names(data_raw)))]
			for(this_gmt in uniq_gmts){
				# this_gmt = uniq_gmts[1]
				data_raw[[this_gmt]] <- ifelse(data_raw[[this_gmt]]<=4, "<= 4", "> 4")
        data_raw[[this_gmt]] <-  factor(data_raw[[this_gmt]], levels=c("<= 4", "> 4"))
      }

      ## HAI categories
      uniq_hais <- names(data_raw)[grepl("HAI", toupper(names(data_raw)))]
      for(this_hai in uniq_hais){
        # this_hai = uniq_hais[1]
        data_raw[[this_hai]] <- ifelse(data_raw[[this_hai]]<=40, "<= 40", "> 40")
        data_raw[[this_hai]] <-  factor(data_raw[[this_hai]], levels=c("<= 40", "> 40"))
      }

      ## age categories
      data_raw$age <- cut(data_raw$age, breaks=c(0,6,8,10,12,20))
      ## gender
      data_raw$gender <- factor(data_raw$male, levels=c("0", "1"), label=c("Female", "Male"))
      data_raw <- data_raw %>% select(-male)
      ## vaccine
      data_meta_reorder$vaccine <- factor(grepl("V1", data_meta_reorder$group), levels=c(TRUE, FALSE), label=c("Vaccinated", "Unvaccinated"))

      stopifnot(all(data_raw$sample==data_meta_reorder$sample))
			data_m <- data_raw %>% select(-group, -sample, -timepoint) 
      id <- apply(data_m, 1, function(x){sum(is.na(x))})>0
      data_m <- data_m[!id,] # remove NA rows
      data_meta_reorder <- data_meta_reorder[!id,]
			stopifnot(length(apply(data_m, 2, function(x){which(is.na(x))}))==0) # make sure no missing values

      data_m <- data_m %>% select(`age`, `gender`, contains("pre infection HAI"), contains("GMT"), everything())
      length_antibody <- ncol(data_m)-1-1-length(uniq_gmts)-length(uniq_hais)

			## MFA
			res.mfa <- MFA(data_m, 
        group = c(1, 1, length(uniq_hais), length(uniq_gmts), length_antibody),
        type = c("n", "n", "n", "n", "s"),
        name.group = c(
          paste0("Age (N=1)"),
          paste0("Gender (N=1)"),
          paste0("HAI (N=", length(uniq_hais), ")"),
          paste0("GMT (N=", length(uniq_gmts), ")"),
          paste0("Antibody (N=", length_antibody, ")")),
        graph = FALSE
        )
      
      df <- as_tibble(res.mfa$eig)
			df <- bind_cols(PC=rownames(res.mfa$eig), df)
			p0 <- fviz_screeplot(res.mfa, ncp=14, addlabels = TRUE)
			file_out <- paste0(cur_path,"screet_plot.pdf")
			ggsave(file_out, width = 8, height = 6)

      group <- get_mfa_var(res.mfa, "group")
      df_contrib <- as_tibble(group$contrib)
      df_contrib$group <- rownames(group$contrib)
      df_contrib <- df_contrib %>% select(group, everything())
      write_csv(df_contrib, paste0(cur_path, "contributions.csv"))

			p1 <- fviz_ellipses(res.mfa, c("gender", "age"), repel = TRUE)
			file_out <- paste0(cur_path, "individual_plot_male_age.pdf")
			ggsave(file_out, width = 9.5, height = 8, plot=p1)

			p2 <- fviz_ellipses(res.mfa, uniq_hais, repel = TRUE)
      file_out <- paste0(cur_path, "individual_plot_HAI.pdf")
      ggsave(file_out, width = 9.5, height = 8, plot=p2)

			p3 <- fviz_ellipses(res.mfa, uniq_gmts, repel = TRUE)
      file_out <- paste0(cur_path, "individual_plot_GMT.pdf")
      ggsave(file_out, width = 9.5, height = 8, plot=p3)

      p_vaxx <- fviz_mfa_ind(res.mfa, geom = "point", col.ind = data_meta_reorder$vaccine, palette=scico::scico(length(unique(data_meta_reorder$vaccine)), end=0.8, palette = "batlow"), legend.title = "Value", title="Vaccination status", addEllipses=T)
      file_out <- paste0(cur_path, "individual_plot_vaxx.pdf")
      ggsave(file_out, width = 6, height = 5, plot=p_vaxx)

      p_groups <- fviz_mfa_ind(res.mfa, geom = "point", col.ind = data_meta_reorder$group, palette=scico::scico(length(unique(data_meta_reorder$group)), end=0.8, palette = "batlow"), legend.title = "Value", title="Groups", addEllipses=T)
      file_out <- paste0(cur_path, "individual_plot_groups.pdf")
      ggsave(file_out, width = 6, height = 5, plot=p_groups)
		}

	}
	
}
