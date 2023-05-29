pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)
 
library(caret)
library(tidyverse)
library(readxl)
library(pcaMethods)
library(patchwork)
library(MASS)
library(mdatools) # https://mdatools.com/docs/plsda.html
library(writexl) 
library(scales)
library(ggallin)
# library(boot)

# read data
source("./helper/deal_with_colnams.R")

data_multiplex_aim2 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 2 multiplex data")
data_multiplex_aim2 <- deal_with_colnames(data_multiplex_aim2)
unique(sapply(strsplit(names(data_multiplex_aim2), " "), function(x){x[3]}))
data_multiplex_aim2$timepoint <- 2

data_multiplex_aim3 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 3 multiplex data")
data_multiplex_aim3 <- deal_with_colnames(data_multiplex_aim3)
unique(sapply(strsplit(names(data_multiplex_aim3), " "), function(x){paste(x[-c(1,2)], collapse = " ")}))
names(data_multiplex_aim3) <- gsub("H3N2/Victoria/2011", "H3-2011", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("Malaysia/2004", "HA-20041", names(data_multiplex_aim3), fixed=T)

split_tmp <- strsplit(data_multiplex_aim3$sample, "d", fixed=T)
data_multiplex_aim3$sample <- sapply(split_tmp, function(x){x[1]})
data_multiplex_aim3$sample <- gsub(" ", "", data_multiplex_aim3$sample, fixed=T)
data_multiplex_aim3$group <- gsub(" ", "", data_multiplex_aim3$group, fixed=T)
data_multiplex_aim3$timepoint <- sapply(split_tmp, function(x){as.numeric(x[2])})
data_multiplex_aim3$timepoint[data_multiplex_aim3$group=="V0S0"] <- 2

data_clinical_aim2 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 2 clinical data")
data_clinical_aim2$sample <- as.character(data_clinical_aim2$sample)
stopifnot(all(data_multiplex_aim2$sample %in% data_clinical_aim2$sample))
data_clinical_aim3 <- read_excel("../data/PCA_data.xlsx", sheet="Aim 3 clinical data")
data_clinical_aim3$sample <- as.character(data_clinical_aim3$sample)
stopifnot(all(data_multiplex_aim3$sample %in% data_clinical_aim3$sample))

names(data_multiplex_aim2) <- gsub("pdm H1F-2009", "H1-stem", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("seas. H3F", "H3-stem", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("seas. H1-2007", "vaxx. H1-2007", names(data_multiplex_aim2), fixed=T)

names(data_multiplex_aim3) <- gsub("pdm H1F-2009", "H1-stem", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H3F", "H3-stem", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H1-2007", "vaxx. H1-2007", names(data_multiplex_aim3), fixed=T)

# build model for aim2 and aim3
all_aims <- c("aim2", "aim3")
data_multiplex_aim3 <- data_multiplex_aim3 %>% filter(timepoint==2) # keep only the d2 data for comparison between the aim2 and aim3 data

dir_rst <- "elastic_net_plsda_results"

for (this_aim in all_aims) { # antibody + HAI
	# this_aim="aim3"
	if(this_aim=="aim2"){
		data_multiplex <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_3_en_plsda"
	} else {
		data_multiplex <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_4_en_plsda"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	cur_path <- paste0("../results/", dir_rst, "/", this_aim, "/", this_model, "/")
	dir.create(cur_path, showWarnings=F, recursive=T)
	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection, -contains("pre_infection_HAI")) %>% mutate_all(as.numeric)	
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	data_m[data_m<0] <- 0
	data_m <- log10(data_m+1) # log transformation
	data_m <- apply(data_m, 2, pcaMethods::prep, scale="uv") # scale the data, with unit variance
	data_m <- cbind(data_m, data_multiplex %>% dplyr::select(contains("pre_infection_HAI")))

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	sampling_factor <- 0.8
	n_repeats <- 2000
	n_obs <- nrow(data_combined)
	set.seed(2023)

	variable_counts <- lapply(seq_len(n_repeats), function(i) {
		print(i)
		idx_choosen <- sample(n_obs, round(n_obs*sampling_factor), replace=FALSE)
		Y <- data_multiplex$infection[idx_choosen]
		X <- model.matrix(response~.,data_combined[idx_choosen,])[,-1]
		# 10-fold CV to find the optimal lambda 
		enet.cv=cv.glmnet(X,Y,alpha=0.5, type="deviance", family="binomial", nfolds=10)
		## Fit lasso model with 100 values for lambda
		enet_mdl = glmnet(X,Y,alpha=0.5,nlambda=100)
		## Extract coefficients at optimal lambda
		out <- coef(enet_mdl,s=enet.cv$lambda.min)
		out <- as.matrix(out)[,1]
		names(out)[out!=0]
	})
	
	count_data <- sort(table(unlist(variable_counts)), decreasing=T)
	count_data <- count_data[-which(names(count_data)=="(Intercept)")] # remove intercerpt
	out_file <- paste0(cur_path, "ranking_from_elastic_net_prediction.xlsx")
	write_xlsx(tibble(varibale=names(count_data), times=count_data), out_file)
	
	## a sequential step-forward algorithm that itera- tively added a single feature into the PLSR (numerical outcome) or PLSDA (categorical outcome) model starting with the feature that had the highest fre- quency of selection, to the lowest frequency of selection. 
	summary_all_models <- lapply(seq_len(length(count_data)), function (j){
		# j=5
		variables <- names(count_data)[seq(j)]
		variables <- gsub("`","",variables)
	
		data_train <- data_combined %>% dplyr::select(all_of(variables))
		set.seed(2023)
		model <- plsda(x=data_train, c=factor(data_combined$response, levels=c(0,1), labels=c("non-infection", "infection")), cv=10)
		tmp <- capture.output(summary(model))
		idx <- grep("(infection)", tmp, fixed=T)
		tmp <- tmp[(idx+2):(idx+3)]
		tmp <- strsplit(tmp, " +")
		names(tmp) <- c("Cal", "cv")
		out <- t(data.frame(tmp))
		colnames(out) <- c("Type", "X cumexpvar", "Y cumexpvar", "TP", "FP", "TN", "FN", "Spec.", "Sens.", "Accuracy")
		out <- as_tibble(out)
		out$num_varibales <- j
		return(out)
	})
	summary_all_models <- bind_rows(summary_all_models)

	out_file <- paste0(cur_path, "plsda_accuracy.xlsx")
	best_model_summary <- summary_all_models  %>% group_by(Type) %>% arrange(Type, desc(Accuracy)) 
	write_xlsx(best_model_summary, out_file)

	best_j <- best_model_summary$num_varibales[1]
	variables <- names(count_data)[seq(best_j)]
	variables <- gsub("`","",variables)
	data_train <- data_combined %>% dplyr::select(all_of(variables))

	data_test <- data_multiplex_test %>% dplyr::select(all_of(variables)) %>% mutate_all(as.numeric)
	check <- apply(data_test, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_test <- data_multiplex_test %>% dplyr::select(all_of(variables)) %>% mutate_all(as.numeric)

	set.seed(2023)
	model <- plsda(x=data_train, c=factor(data_combined$response, levels=c(0,1), labels=c("non-infection", "infection")), cv=10)
	summary(model)
	str(model,1)
	# plotXScores(model)
	# plotPredictions(model)

	## score plot
	model$res$cal$xdecomp$expvar
	round(model$res$cal$xdecomp$expvar)
	
	df_plot <- tibble(Group=ifelse(data_combined$response==1, "infection", "non-infection"), lv1=model$res$cal$xdecomp$scores[,1], lv2=model$res$cal$xdecomp$scores[,2])
	p_out <- ggplot(df_plot)+
		geom_point(aes(x=lv1, y=lv2, color=Group, shape=Group), size=3, alpha=0.8)+
		xlab(paste0("LV1 (", round(model$res$cal$xdecomp$expvar[1],2), "%)"))+
		ylab(paste0("LV2 (", round(model$res$cal$xdecomp$expvar[2],2), "%)"))+
		scale_color_manual(values=c("dark red", "dark blue"))+
		theme_minimal()+
		theme(legend.position="bottom")+
		NULL
	out_file <- paste0(cur_path, "plot_score.pdf")
	ggsave(out_file, width=6, height=6, plot=p_out)

	## loading polot
	df_plot <- tibble(variables=names(model$xloadings[,1]), loadings=model$xloadings[,1], Group=ifelse(loadings>0, "1", "0"))
	df_plot <- df_plot %>% arrange(loadings)
	df_plot$variables <- factor(df_plot$variables, levels=df_plot$variables)
	
	p_out <- ggplot(df_plot)+
		geom_col(aes(y=loadings,x=variables,fill=Group))+
		scale_y_continuous(trans = pseudolog10_trans, breaks = c(-100, -10, 0, 10, 100 ,1000),labels = c(expression(-10^2, -10^1, 0, 10^1, 10^2, 10^3)))+
		ylab("Loadings on component 1")+
		coord_flip()+
		scale_fill_manual(values=c("dark blue", "dark red"))+
		theme_minimal()+
		guides(fill="none")
	
	out_file <- paste0(cur_path, "plot_loading.pdf")
	ggsave(out_file, width=6, height=6, plot=p_out)


	res = predict(model, data_test, factor(data_multiplex_test$infection, levels=c(FALSE,TRUE), labels=c("non-infection", "infection")))
	out_prediction <- capture.output(summary(res))
	out_file <- paste0(cur_path, "plsda_prediction.txt")
	writeLines(out_prediction, out_file)

}
