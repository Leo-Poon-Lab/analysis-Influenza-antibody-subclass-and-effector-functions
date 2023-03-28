library(readxl)
library(tidyverse)

# library(naturalsort)
# library("corrplot")
library(scico)
library(pcaMethods)
library(FactoMineR)
library("factoextra")
# library("PerformanceAnalytics")
# library(ggpubr)
library(ggrepel)
library(ggforce)
# library(sp)
library(fasano.franceschini.test)
# library(patchwork)

# self defined functions for significance testing
ks_test_2d <- function(data, col_name, res_pca){
	groups <- as.character(sort(unique(data[[col_name]])))
	if(length(groups)<2){return(NA)}
	pair_groups <- combn(groups, 2)
	ks_test_2d <- apply(pair_groups, 2, function(pair_t) {
		# pair_t=pair_groups[,1]
		dim1_t <- res_pca$ind$coord[,1]
		dim2_t <- res_pca$ind$coord[,2]
		
		dim1_t1 <- dim1_t[which(data[[col_name]] == pair_t[1])]
		dim1_t2 <- dim1_t[which(data[[col_name]] == pair_t[2])]
		dim2_t1 <- dim2_t[which(data[[col_name]] == pair_t[1])]
		dim2_t2 <- dim2_t[which(data[[col_name]] == pair_t[2])]

		g1Data <- data.frame(x = dim1_t1, y = dim2_t1)
		g1Data <- g1Data[!is.na(g1Data[,1]),]
		g2Data <- data.frame(x = dim1_t2, y = dim2_t2)
		g2Data <- g2Data[!is.na(g2Data[,1]),]

		if(nrow(g1Data)>1 && nrow(g2Data)>1){
			ks_test_2d <- fasano.franceschini.test(g1Data, g2Data)
			ks_test_2d$p.value
		} else {
			NA
		}
		
	})
	df_ks_test <- as_tibble(as.data.frame(t(pair_groups)))
	df_ks_test$p_value <- ks_test_2d
	df_ks_test$`<0.05?` <- df_ks_test$p_value < 0.05
	df_ks_test$`<0.01?` <- df_ks_test$p_value < 0.01
	names(df_ks_test)[1:2] <- c("Group_1", "Group_2")
	df_ks_test$variable <- col_name
	return(df_ks_test)
}

# wilcox_test_d1 <- function(data, col_name, res_pca){
# 	groups <- as.character(sort(unique(data[[col_name]])))
# 	pair_groups <- combn(groups, 2)
# 	wilc_test_2d <- apply(pair_groups, 2, function(pair_t) {
# 		dim1_t <- res_pca$ind$coord[,1]
		
# 		dim1_t1 <- dim1_t[data[[col_name]] == pair_t[1]]
# 		dim1_t2 <- dim1_t[data[[col_name]] == pair_t[2]]
		
# 		rst_test <- wilcox.test(dim1_t1, dim1_t2)
# 		# str(ks_test_2d)
# 		rst_test$p.value
# 	})
# 	df_wilc_test <- as_tibble(as.data.frame(t(pair_groups)))
# 	df_wilc_test$p_value <- wilc_test_2d
# 	df_wilc_test$`<0.05?` <- df_wilc_test$p_value < 0.05
# 	df_wilc_test$`<0.01?` <- df_wilc_test$p_value < 0.01
# 	names(df_wilc_test)[1:2] <- c("Group_1", "Group_2")
# 	return(df_wilc_test)
# }

# read data
source("./helper/deal_with_colnams.R")

data_multiplex_aim2 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 2 multiplex data")
data_multiplex_aim2 <- deal_with_colnames(data_multiplex_aim2)
unique(sapply(strsplit(names(data_multiplex_aim2), " "), function(x){x[3]}))
data_multiplex_aim2$timepoint <- 2
names(data_multiplex_aim2) <- gsub("seas. H5", "av. H5", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("seas. H7", "av. H7", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("seas. H9", "av. H9", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("av. HA", "av.", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("H1F", "H1-stem", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("H3F", "H3-stem", names(data_multiplex_aim2), fixed=T)

sort(names(data_multiplex_aim2))

data_multiplex_aim3 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 3 multiplex data")
data_multiplex_aim3 <- deal_with_colnames(data_multiplex_aim3)
unique(sapply(strsplit(names(data_multiplex_aim3), " "), function(x){paste(x[-c(1,2)], collapse = " ")}))
names(data_multiplex_aim3) <- gsub("H3N2/Victoria/2011", "H3-2011", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("Malaysia/2004", "HA-20041", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H5", "av. H5", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H7", "av. H7", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H9", "av. H9", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("av. HA", "av.", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("H1F", "H1-stem", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("H3F", "H3-stem", names(data_multiplex_aim3), fixed=T)

split_tmp <- strsplit(data_multiplex_aim3$sample, "d", fixed=T)
data_multiplex_aim3$sample <- sapply(split_tmp, function(x){x[1]})
data_multiplex_aim3$sample <- gsub(" ", "", data_multiplex_aim3$sample, fixed=T)
data_multiplex_aim3$group <- gsub(" ", "", data_multiplex_aim3$group, fixed=T)
data_multiplex_aim3$timepoint <- sapply(split_tmp, function(x){as.numeric(x[2])})
data_multiplex_aim3$timepoint[data_multiplex_aim3$group=="V0S0"] <- 2

data_clinical_aim2 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 2 clinical data")
stopifnot(all(data_multiplex_aim2$sample %in% data_clinical_aim2$sample))
data_clinical_aim3 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 3 clinical data")
stopifnot(all(data_multiplex_aim3$sample %in% data_clinical_aim3$sample))

all_aims <- c("aim2", "aim3")
all_plots <- c()
all_plots_names <- c()

for (this_aim in all_aims) { # PCA for raw responses
	# this_aim=all_aims[2]
	if(this_aim=="aim2"){
		data_multiplex <- data_multiplex_aim2
		data_clinical <- data_clinical_aim2
	} else {
		data_multiplex <- data_multiplex_aim3
		data_clinical <- data_clinical_aim3
	}

	uniq_groups <- as.list(unique(data_multiplex$group))
	uniq_groups <- c(list(unique(data_multiplex$group)), uniq_groups) # additional combinational datasets

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
			cur_path <- paste0("../results/PCA_results/", this_aim, "/", paste0(this_group, collapse = "-"), "/timepoint_", paste0(this_timepoint, collapse = "-"), "/")
			dir.create(cur_path, showWarnings=F, recursive=T)
			print(cur_path)

			########### PCA #####################################

			data_raw <- data_multiplex_subgroup %>% filter(timepoint%in%this_timepoint)
			data_meta_reorder <- left_join(data_raw %>% select(sample, timepoint), data_meta, "sample")
			stopifnot(all(data_raw$sample==data_meta_reorder$sample))

			data_m <- data_raw %>% select(-group, -sample, -timepoint) 
			stopifnot(length(apply(data_m, 2, function(x){which(is.na(x))}))==0) # make sure no missing values
			data_m <- data_m %>% mutate_all(as.numeric)
			data_completed <- prep(data_m, scale= "uv") # scale the data, with unit variance

			## PCA
			res.pca <- PCA(data_completed, scale.unit = TRUE, ncp = 6, graph = FALSE)
			df <- as_tibble(res.pca$eig)
			df <- bind_cols(PC=rownames(res.pca$eig), df)
			p0 <- fviz_screeplot(res.pca, ncp=14, addlabels = TRUE)
			file_out <- paste0(cur_path,"screet_plot.pdf")
			ggsave(file_out, width = 8, height = 6)
			all_plots <- c(all_plots, p0)
			all_plots_names <- c(all_plots_names, file_out)

			# (p_tmp <- fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, select.var = list(contrib=5)))

			p <- fviz_contrib(res.pca, "var", axes = c(1,2))
			data_high_con <- p$data %>% filter(contrib>=1/nrow(p$data)*100) %>% arrange(contrib)
			data_scree <- p0$data
			data_scree$label <- paste0(rownames(data_scree), " (", round(data_scree$eig,2), "%)")
			write_csv(p$data %>% arrange(desc(contrib)), paste0(cur_path,"contributions.csv"))

			p <- fviz_contrib(res.pca, "var", axes = c(1,2), top=nrow(data_high_con), xtickslab.rt=90)
			file_out <- paste0(cur_path, "contributions.pdf")
			ggsave(file_out, width = 8, height = 5)
			all_plots <- c(all_plots, p)
			all_plots_names <- c(all_plots_names, file_out)

			df_plot <- as.data.frame(res.pca$var$coord)
			df_plot$label <- rownames(df_plot)
			df_plot$label_parsed <- df_plot$label
			df_plot$label_parsed <- gsub("Fcr", "Fc\u03B3", df_plot$label_parsed, fixed=T)
			
			df_plot$antibody <- factor(sapply(strsplit(df_plot$label_parsed, " "), function(x) {x[1]}))
			colors <- scico(length(unique(df_plot$antibody)), end=0.8, palette="batlow")
			df_plot$color <- colors[df_plot$antibody]
			 
			df_plot$color[!df_plot$label %in% as.character(data_high_con$name)] <- "#bababa"
			df_plot$color_seg <- ifelse(df_plot$color=="#bababa", "#bababa", "#000000")
			
			p1 <- ggplot(df_plot)+
				geom_circle(aes(x0=0, y0=0, r=1))+
				theme_minimal()+
				geom_hline(yintercept=0, linetype="dashed")+
				geom_vline(xintercept=0, linetype="dashed")+
				geom_segment(aes_string(x = 0, y = 0, xend = "Dim.1", yend = "Dim.2", color="color"),alpha=0.9, arrow = grid::arrow(length = grid::unit(0.2,"cm")))+
				geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label_parsed), color="#bababa", size=1, bg.color = "white", segment.size=0.5, data=. %>% filter(color=="#bababa"), max.overlaps=100)+
				geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label_parsed, color=color), size=2.1, bg.color = "white", segment.size=0., data=. %>% filter(color!="#bababa"), max.overlaps=100)+
				scale_color_identity("Antibody", labels=levels(df_plot$antibody), breaks=colors, guide="legend")+
				xlab(data_scree$label[1])+
				ylab(data_scree$label[2])+	
				# coord_fixed()+
				NULL
			file_out <- paste0(cur_path, "variables_plot.svg")
			ggsave(file_out, width = 9.5, height = 8)
			file_out <- paste0(cur_path, "variables_plot.jpg")
			ggsave(file_out, width = 9.5, height = 8)
			all_plots <- c(all_plots, p1)
			all_plots_names <- c(all_plots_names, file_out)

			## GMT rises
			uniq_gmts <- names(data_meta_reorder)[grepl("GMT", toupper(names(data_meta_reorder)))]
			for(this_gmt in uniq_gmts){
				# this_gmt = uniq_gmts[1]
				data_meta_reorder$gmt_group <- ifelse(data_meta_reorder[[this_gmt]]<=4, "<= 4", "> 4")
				if(length(unique(data_meta_reorder$gmt_group))>1){
					(p_gmt <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$gmt_group, palette=scico(length(unique(data_meta_reorder$gmt_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_gmt, addEllipses=T))
				} else {
					(p_gmt <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$gmt_group, palette=scico(length(unique(data_meta_reorder$gmt_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_gmt))
				}

				file_out <- paste0(cur_path, "individuals_plot_by_", this_gmt, ".pdf")
				ggsave(file_out, width = 8, height = 6)
				all_plots <- c(all_plots, p_gmt)
				all_plots_names <- c(all_plots_names, file_out)

				(df_ks_test_2d <- ks_test_2d(data=data_meta_reorder, col_name="gmt_group", res_pca=res.pca))
				if(is_tibble(df_ks_test_2d)){write_csv(df_ks_test_2d, paste0(cur_path, "individuals_plot_by_", this_gmt, "_2d_test.csv"))}
				# wilcox_test_d1(data_meta_reorder, this_gmt, res.pca)
			}
			
			## HAI
			uniq_hai <- names(data_meta_reorder)[grepl("HAI", toupper(names(data_meta_reorder)))]
			for(this_hai in uniq_hai){
				# this_hai = uniq_hai[1]
				data_meta_reorder$hai_group <- ifelse(data_meta_reorder[[this_hai]]<=40, "<= 40", "> 40")
				if(length(unique(data_meta_reorder$hai_group))>1){
					(p_hai <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$hai_group, palette=scico(length(unique(data_meta_reorder$hai_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_hai, addEllipses=T))
				} else {
					(p_hai <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$hai_group, palette=scico(length(unique(data_meta_reorder$hai_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_hai))
				}				

				file_out <- paste0(cur_path, "individuals_plot_by_", this_hai, ".pdf")
				ggsave(file_out, width = 8, height = 6)
				all_plots <- c(all_plots, p_hai)
				all_plots_names <- c(all_plots_names, file_out)

				(df_ks_test_2d <- ks_test_2d(data_meta_reorder, "hai_group", res.pca))
				if(is_tibble(df_ks_test_2d)){write_csv(df_ks_test_2d, paste0(cur_path, "individuals_plot_by_", this_hai, "_2d_test.csv"))}
				# wilcox_test_d1(data_meta_reorder, this_hai, res.pca)
			}

			## Age
			data_meta_reorder$Age_group <- cut(data_meta_reorder$age, breaks=c(0,6,8,10,12,20))
			(p_age <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$Age_group, palette=scico(length(unique(data_meta_reorder[["Age_group"]])), palette = "lapaz"), legend.title = "Age", title="Age"))
			file_out <- paste0(cur_path, "individuals_plot_by_age.pdf")
			all_plots <- c(all_plots, p_age)
			all_plots_names <- c(all_plots_names, file_out)
			ggsave(file_out, width = 8, height = 6)
			(df_ks_test_2d_agegroup <- ks_test_2d(data_meta_reorder, "Age_group", res.pca))
			if(!is.na(df_ks_test_2d_agegroup)){write_csv(df_ks_test_2d_agegroup, paste0(cur_path, "individuals_plot_by_age_2d_test.csv"))}
			# wilcox_test_d1(data_meta, "Age_group", res.pca)

			## Plots for combined data
			### if groups more than one, for aim2&aim3
			if(length(this_group)>1){
				(p_groups <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$group, palette=scico(length(unique(data_meta_reorder[["group"]])), palette = "lapaz", end=0.8), legend.title = "Group", title="Group", addEllipses=T))
				file_out <- paste0(cur_path, "individuals_plot_by_group.pdf")
				all_plots <- c(all_plots, p_groups)
				all_plots_names <- c(all_plots_names, file_out)
				ggsave(file_out, width = 8, height = 6)
				(df_ks_test_2d_group <- ks_test_2d(data_meta_reorder, "group", res.pca))
				if(!is.na(df_ks_test_2d_group)){write_csv(df_ks_test_2d_group, paste0(cur_path, "individuals_plot_by_group_2d_test.csv"))}
			}
			
			### if timepoint more than one, for aim3 v1s1 group
			if(length(this_timepoint)>1){
				data_meta_reorder$timepoint <- as.character(data_meta_reorder$timepoint)
				(p_timepoint <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$timepoint, palette=scico(length(unique(data_meta_reorder[["timepoint"]])), palette = "tokyo", end=0.8), legend.title = "Timepoint", title="Timepoint", addEllipses=T))
				file_out <- paste0(cur_path, "individuals_plot_by_timepoint.pdf")
				all_plots <- c(all_plots, p_timepoint)
				all_plots_names <- c(all_plots_names, file_out)
				ggsave(file_out, width = 8, height = 6)
				(df_ks_test_2d_timepoint <- ks_test_2d(data_meta_reorder, "timepoint", res.pca))
				if(!is.na(df_ks_test_2d_timepoint)){write_csv(df_ks_test_2d_timepoint, paste0(cur_path, "individuals_plot_by_timepoint_2d_test.csv"))}
			}

			########### PCA #####################################

		}

	}
	
}



# for (this_aim in all_aims) { # PCA for fold-change diff responses
	
# 	# as requested, we need to analyze:
# 	# (1) In the aim 2 data, is it possible to look at fold change difference in the variables plots? I.e. V1S0 vs V1S1 and V0S1 vs V0S0? To make 2 new plots? I am trying to look at them side by side to see any patterns but a plot of their differences would be more informative?
# 	# (2) Similarly could this be do with aim 3 V1S1 3x timepoint data, i.e. to make 2 plots for pre-vaccine versus post-vaccine, and post-vaccine versus post-infection?

# 	# this_aim="aim2"
# 	if(this_aim=="aim2"){
# 		data_multiplex <- data_multiplex_aim2
# 		data_clinical <- data_clinical_aim2
# 		pair_type <- "group"
# 		pairs <- list(c("V1S0", "V1S1"), c("V0S1", "V0S0"))
# 	} else {
# 		data_multiplex <- data_multiplex_aim3
# 		data_clinical <- data_clinical_aim3
# 		pair_type <- "timepoint"
# 		pairs <- list(c(1, 2), c(2, 3))
# 	}

# 	for (this_pair in pairs){
# 		# this_pair=c("V1S0", "V1S1")
# 		if(this_aim=="aim2"){
# 			this_group <- rev(this_pair)
# 			this_timepoint <- 2
			
# 			data_multiplex %>% filter(group %in% this_group) %>% group_by(sample) %>% summarise(N=n()) %>% filter(N>1)
# 			data_multiplex_subgroup1 <- data_multiplex %>% filter(group %in% this_group[1])
# 			data_multiplex_subgroup2 <- data_multiplex %>% filter(group %in% this_group[2])
# 			data_multiplex_subgroup <- 

# 		} else {
# 			this_group <- "V1S1"
# 			this_timepoint <- rev(this_pair)
# 		}

		
# 		data_multiplex_subgroup

# 		data_meta <- data_clinical %>% filter(group%in%this_group)
# 		data_meta$sample <- as.character(data_meta$sample)
		
# 		cur_path <- paste0("../results/PCA_results/", this_aim, "/", paste0(this_group, collapse = "_fold_change_"), "/timepoint_", paste0(this_timepoint, collapse = "_fold_change_"), "/")
# 		dir.create(cur_path, showWarnings=F, recursive=T)
# 		print(cur_path)

# 		########### PCA #####################################

# 		data_raw <- data_multiplex_subgroup %>% filter(timepoint%in%this_timepoint)
# 		data_meta_reorder <- left_join(data_raw %>% select(sample, timepoint), data_meta, "sample")
# 		stopifnot(all(data_raw$sample==data_meta_reorder$sample))

# 		data_m <- data_raw %>% select(-group, -sample, -timepoint) 
# 		stopifnot(length(apply(data_m, 2, function(x){which(is.na(x))}))==0) # make sure no missing values
# 		data_m <- data_m %>% mutate_all(as.numeric)
# 		data_completed <- prep(data_m, scale= "uv") # scale the data, with unit variance

# 		## PCA
# 		res.pca <- PCA(data_completed, scale.unit = TRUE, ncp = 6, graph = FALSE)
# 		df <- as_tibble(res.pca$eig)
# 		df <- bind_cols(PC=rownames(res.pca$eig), df)
# 		p0 <- fviz_screeplot(res.pca, ncp=14, addlabels = TRUE)
# 		file_out <- paste0(cur_path,"screet_plot.pdf")
# 		ggsave(file_out, width = 8, height = 6)
# 		all_plots <- c(all_plots, p0)
# 		all_plots_names <- c(all_plots_names, file_out)

# 		# (p_tmp <- fviz_pca_var(res.pca, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, select.var = list(contrib=5)))

# 		p <- fviz_contrib(res.pca, "var", axes = c(1,2))
# 		data_high_con <- p$data %>% filter(contrib>=1/nrow(p$data)*100) %>% arrange(contrib)
# 		data_scree <- p0$data
# 		data_scree$label <- paste0(rownames(data_scree), " (", round(data_scree$eig,2), "%)")
# 		write_csv(p$data %>% arrange(desc(contrib)), paste0(cur_path,"contributions.csv"))

# 		p <- fviz_contrib(res.pca, "var", axes = c(1,2), top=nrow(data_high_con), xtickslab.rt=90)
# 		file_out <- paste0(cur_path, "contributions.pdf")
# 		ggsave(file_out, width = 8, height = 5)
# 		all_plots <- c(all_plots, p)
# 		all_plots_names <- c(all_plots_names, file_out)

# 		df_plot <- as.data.frame(res.pca$var$coord)
# 		df_plot$label <- rownames(df_plot)

# 		df_plot$antibody <- factor(sapply(strsplit(df_plot$label, " "), function(x) {x[1]}))
# 		colors <- scico(length(unique(df_plot$antibody)), end=0.8, palette="batlow")
# 		df_plot$color <- colors[df_plot$antibody]
			
# 		df_plot$color[!df_plot$label %in% as.character(data_high_con$name)] <- "#bababa"
# 		df_plot$color_seg <- ifelse(df_plot$color=="#bababa", "#bababa", "#000000")

# 		p1 <- ggplot(df_plot)+
# 			geom_circle(aes(x0=0, y0=0, r=1))+
# 			theme_minimal()+
# 			geom_hline(yintercept=0, linetype="dashed")+
# 			geom_vline(xintercept=0, linetype="dashed")+
# 			geom_segment(aes_string(x = 0, y = 0, xend = "Dim.1", yend = "Dim.2", color="color"),alpha=0.9, arrow = grid::arrow(length = grid::unit(0.2,"cm")))+
# 			geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label), color="#bababa", size=1, bg.color = "white", segment.size=0.5, data=. %>% filter(color=="#bababa"), max.overlaps=100)+
# 			geom_text_repel(aes(x=Dim.1, y=Dim.2, label=label, color=color), size=2, bg.color = "white", segment.size=0., data=. %>% filter(color!="#bababa"), max.overlaps=100)+
# 			scale_color_identity("Antibody", labels=levels(df_plot$antibody), breaks=colors, guide="legend")+
# 			xlab(data_scree$label[1])+
# 			ylab(data_scree$label[2])+	
# 			# coord_fixed()+
# 			NULL
# 		file_out <- paste0(cur_path, "variables_plot.pdf")
# 		ggsave(file_out, width = 10, height = 8)
# 		all_plots <- c(all_plots, p1)
# 		all_plots_names <- c(all_plots_names, file_out)

# 		## GMT rises
# 		uniq_gmts <- names(data_meta_reorder)[grepl("GMT", toupper(names(data_meta_reorder)))]
# 		for(this_gmt in uniq_gmts){
# 			# this_gmt = uniq_gmts[1]
# 			data_meta_reorder$gmt_group <- ifelse(data_meta_reorder[[this_gmt]]<=4, "<= 4", "> 4")
# 			if(length(unique(data_meta_reorder$gmt_group))>1){
# 				(p_gmt <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$gmt_group, palette=scico(length(unique(data_meta_reorder$gmt_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_gmt, addEllipses=T))
# 			} else {
# 				(p_gmt <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$gmt_group, palette=scico(length(unique(data_meta_reorder$gmt_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_gmt))
# 			}

# 			file_out <- paste0(cur_path, "individuals_plot_by_", this_gmt, ".pdf")
# 			ggsave(file_out, width = 8, height = 6)
# 			all_plots <- c(all_plots, p_gmt)
# 			all_plots_names <- c(all_plots_names, file_out)

# 			(df_ks_test_2d <- ks_test_2d(data=data_meta_reorder, col_name="gmt_group", res_pca=res.pca))
# 			if(is_tibble(df_ks_test_2d)){write_csv(df_ks_test_2d, paste0(cur_path, "individuals_plot_by_", this_gmt, "_2d_test.csv"))}
# 			# wilcox_test_d1(data_meta_reorder, this_gmt, res.pca)
# 		}
		
# 		## HAI
# 		uniq_hai <- names(data_meta_reorder)[grepl("HAI", toupper(names(data_meta_reorder)))]
# 		for(this_hai in uniq_hai){
# 			# this_hai = uniq_hai[1]
# 			data_meta_reorder$hai_group <- ifelse(data_meta_reorder[[this_hai]]<=40, "<= 40", "> 40")
# 			if(length(unique(data_meta_reorder$hai_group))>1){
# 				(p_hai <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$hai_group, palette=scico(length(unique(data_meta_reorder$hai_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_hai, addEllipses=T))
# 			} else {
# 				(p_hai <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$hai_group, palette=scico(length(unique(data_meta_reorder$hai_group)), end=0.8, palette = "batlow"), legend.title = "Value", title=this_hai))
# 			}				

# 			file_out <- paste0(cur_path, "individuals_plot_by_", this_hai, ".pdf")
# 			ggsave(file_out, width = 8, height = 6)
# 			all_plots <- c(all_plots, p_hai)
# 			all_plots_names <- c(all_plots_names, file_out)

# 			(df_ks_test_2d <- ks_test_2d(data_meta_reorder, "hai_group", res.pca))
# 			if(is_tibble(df_ks_test_2d)){write_csv(df_ks_test_2d, paste0(cur_path, "individuals_plot_by_", this_hai, "_2d_test.csv"))}
# 			# wilcox_test_d1(data_meta_reorder, this_hai, res.pca)
# 		}

# 		## Age
# 		data_meta_reorder$Age_group <- cut(data_meta_reorder$age, breaks=c(0,6,8,10,12,20))
# 		(p_age <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$Age_group, palette=scico(length(unique(data_meta_reorder[["Age_group"]])), palette = "lapaz"), legend.title = "Age", title="Age"))
# 		file_out <- paste0(cur_path, "individuals_plot_by_age.pdf")
# 		all_plots <- c(all_plots, p_age)
# 		all_plots_names <- c(all_plots_names, file_out)
# 		ggsave(file_out, width = 8, height = 6)
# 		(df_ks_test_2d_agegroup <- ks_test_2d(data_meta_reorder, "Age_group", res.pca))
# 		if(!is.na(df_ks_test_2d_agegroup)){write_csv(df_ks_test_2d_agegroup, paste0(cur_path, "individuals_plot_by_age_2d_test.csv"))}
# 		# wilcox_test_d1(data_meta, "Age_group", res.pca)

# 		## Plots for combined data
# 		### if groups more than one, for aim2&aim3
# 		if(length(this_group)>1){
# 			(p_groups <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$group, palette=scico(length(unique(data_meta_reorder[["group"]])), palette = "lapaz", end=0.8), legend.title = "Group", title="Group", addEllipses=T))
# 			file_out <- paste0(cur_path, "individuals_plot_by_group.pdf")
# 			all_plots <- c(all_plots, p_groups)
# 			all_plots_names <- c(all_plots_names, file_out)
# 			ggsave(file_out, width = 8, height = 6)
# 			(df_ks_test_2d_group <- ks_test_2d(data_meta_reorder, "group", res.pca))
# 			if(!is.na(df_ks_test_2d_group)){write_csv(df_ks_test_2d_group, paste0(cur_path, "individuals_plot_by_group_2d_test.csv"))}
# 		}
		
# 		### if timepoint more than one, for aim3 v1s1 group
# 		if(length(this_timepoint)>1){
# 			data_meta_reorder$timepoint <- as.character(data_meta_reorder$timepoint)
# 			(p_timepoint <- fviz_pca_ind(res.pca, geom.ind = "point", col.ind = data_meta_reorder$timepoint, palette=scico(length(unique(data_meta_reorder[["timepoint"]])), palette = "tokyo", end=0.8), legend.title = "Timepoint", title="Timepoint", addEllipses=T))
# 			file_out <- paste0(cur_path, "individuals_plot_by_timepoint.pdf")
# 			all_plots <- c(all_plots, p_timepoint)
# 			all_plots_names <- c(all_plots_names, file_out)
# 			ggsave(file_out, width = 8, height = 6)
# 			(df_ks_test_2d_timepoint <- ks_test_2d(data_meta_reorder, "timepoint", res.pca))
# 			if(!is.na(df_ks_test_2d_timepoint)){write_csv(df_ks_test_2d_timepoint, paste0(cur_path, "individuals_plot_by_timepoint_2d_test.csv"))}
# 		}	
# 	}

	

	
# }