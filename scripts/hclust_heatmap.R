library(readxl)
library(tidyverse)

# library(naturalsort)
library("pheatmap")
library(scico)
library(pcaMethods)
library(FactoMineR)
library("factoextra")

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

dir_rst <- "hclust_results"

for (this_aim in all_aims) { # hclust for log-transformed and scaled responses
	# this_aim=all_aims[1]
	if(this_aim=="aim2"){
		data_multiplex <- data_multiplex_aim2
		data_clinical <- data_clinical_aim2
	} else {
		data_multiplex <- data_multiplex_aim3
		data_clinical <- data_clinical_aim3
	}

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
			# this_timepoint <- this_timepoint[[1]]
			cur_path <- paste0(here::here("results/"), dir_rst, "/", this_aim, "/", paste0(this_group, collapse = "-"), "/timepoint_", paste0(this_timepoint, collapse = "-"), "/")
			dir.create(cur_path, showWarnings=F, recursive=T)
			print(cur_path)

			########### hclust #####################################

			data_raw <- data_multiplex_subgroup %>% filter(timepoint%in%this_timepoint)
			data_meta_reorder <- left_join(data_raw %>% dplyr::select(sample, timepoint), data_meta, "sample")
			stopifnot(all(data_raw$sample==data_meta_reorder$sample))

			data_m <- data_raw %>% dplyr::select(-group, -sample, -timepoint) 
			stopifnot(length(apply(data_m, 2, function(x){which(is.na(x))}))==0) # make sure no missing values
			data_m <- data_m %>% mutate_all(as.numeric)
			data_m[data_m<0] <- 0
			data_m <- log10(data_m+1) # log transformation
			data_completed <- prep(data_m, scale= "uv") # scale the data, with unit variance

			data_completed <- t(data_completed)
			colnames(data_completed) <- data_meta_reorder$sample

			## pheatmap
			row_group <- data.frame(antibody=sapply(strsplit(rownames(data_completed), " "), function(x) {x[1]}), subtype=sapply(strsplit(rownames(data_completed), " "), function(x) {paste(x[-1], collapse=" ")}))
			rownames(row_group) <- rownames(data_completed)
			
			names(data_meta_reorder)[grepl("HAI", toupper(names(data_meta_reorder)))]
			col_group <- data.frame(group=data_meta_reorder$group, pre_infection_HAI_pH1=ifelse(data_meta_reorder[["pre infection HAI pH1"]]<=40, "<= 40", "> 40"), pre_infection_HAI_sH1=ifelse(data_meta_reorder[["pre infection HAI sH1"]]<=40, "<= 40", "> 40"))
			rownames(col_group) <- data_meta_reorder$sample
			
			colors_ab <- scico(length(unique(row_group$antibody)), end=0.8, palette="batlow")
			names(colors_ab) <- unique(row_group$antibody)
			colors_subtype <- scico(length(unique(row_group$subtype)), end=0.8, palette="romaO")
			names(colors_subtype) <- unique(row_group$subtype)
			ann_colors = list(
				group = c("V0S0"="black", "V0S1"="red", "V1S0"="blue", "V1S1"="grey40"),
				pre_infection_HAI_pH1 = c("<= 40"="grey", "> 40"="dark red"),
				pre_infection_HAI_sH1 = c("<= 40"="grey", "> 40"="dark red"),
				antibody=colors_ab,
				subtype=colors_subtype
				)
			
			
			p0 <- pheatmap(data_completed, annotation_row = row_group, annotation_col = col_group, fontsize_row=5, color=colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize_col=5, annotation_colors = ann_colors)
			file_out <- paste0(cur_path,"hclust_", this_aim, "_", paste0(this_group, collapse = "-"), "_timpoint", this_timepoint, ".pdf")
			ggsave(file_out, width = 10, height = 10, plot=p0)
			all_plots <- c(all_plots, p0)
			all_plots_names <- c(all_plots_names, file_out)

			# (p_tmp <- fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, select.var = list

		}

	}
	
}


