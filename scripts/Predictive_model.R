# In this script, we try to build prediction models seperately for aim2 and aim3 data,
# and we test the model with the data from the other aim.

library(caret)
library(tidyverse)
library(readxl)
library(pcaMethods)
library(patchwork)
library(MASS)
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
names(data_multiplex_aim2) <- gsub("Bris H1-2007", "vaxx H1-2007", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub("Bris H3-2007", "vaxx H3-2007", names(data_multiplex_aim2), fixed=T)

names(data_multiplex_aim3) <- gsub("pdm H1F-2009", "H1-stem", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H3F", "H3-stem", names(data_multiplex_aim3), fixed=T)
names(data_multiplex_aim3) <- gsub("seas. H1-2007", "vaxx. H1-2007", names(data_multiplex_aim3), fixed=T)

# build model for aim2 and aim3
all_aims <- c("aim2", "aim3")
data_multiplex_aim3 <- data_multiplex_aim3 %>% filter(timepoint==2) # keep only the d2 data for comparison between the aim2 and aim3 data

for (this_aim in all_aims) { # antibody alone
	# this_aim="aim2"
	# this_aim="aim3"
	if(this_aim=="aim2"){ 
		data_multiplex <- data_multiplex_aim2
		data_multiplex_test <- data_multiplex_aim3
		this_model <- "model_1"
	} else {
		data_multiplex <- data_multiplex_aim3
		data_multiplex_test <- data_multiplex_aim2
		this_model <- "model_2"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	# data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

	if(this_aim=="aim2"){
		write_csv(as_tibble(data_m), "../results/training_data_aim2.csv")
		writeLines(as.character(as.numeric(data_multiplex$infection)), "../results/response_aim2.txt")
	}

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	data_combined$response
	model <- glm(response~., data=data_combined, family="binomial")
	model %>% stepAIC(trace = FALSE) %>% summary()
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-response), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0("../results/logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0("../results/logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	if(this_aim=="aim2"){
		write_csv(as_tibble(data_new), "../results/training_data_aim3.csv")
		writeLines(as.character(as.numeric(data_multiplex_test$infection)), "../results/response_aim3.txt")
	}

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0("../results/logistic_prediction_", this_model, ".txt"))
}

for (this_aim in all_aims) { # antibody + HAI
	# this_aim="aim2"
	# this_aim="aim3"
	if(this_aim=="aim2"){
		data_multiplex <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_3"
	} else {
		data_multiplex <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_4"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	# data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

	if(this_aim=="aim2"){
		write_csv(as_tibble(data_m), "../results/training_data_hai_antibody_aim2.csv")
		writeLines(as.character(as.numeric(data_multiplex$infection)), "../results/response_aim2.txt")
	}

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-`response`), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0("../results/logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0("../results/logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	if(this_aim=="aim2"){
		write_csv(as_tibble(data_new), "../results/training_data_hai_antibody_aim3.csv")
		writeLines(as.character(as.numeric(data_multiplex_test$infection)), "../results/response_aim3.txt")
	}

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0("../results/logistic_prediction_", this_model, ".txt"))

}

for (this_aim in all_aims) { # HAI alone
	# this_aim="aim2"
	# this_aim="aim3"
	if(this_aim=="aim2"){
		data_multiplex <- left_join(data_multiplex_aim2 %>% dplyr::select(sample, group), data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim3 %>% dplyr::select(sample, group), data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_5"
	} else {
		data_multiplex <- left_join(data_multiplex_aim3 %>% dplyr::select(sample, group), data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim2 %>% dplyr::select(sample, group), data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_6"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	# data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

	if(this_aim=="aim2"){
		write_csv(as_tibble(data_m), "../results/training_data_hai_only_aim2.csv")
		writeLines(as.character(as.numeric(data_multiplex$infection)), "../results/response_aim2.txt")
	}

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	model <- glm(response~., data=data_combined, family="binomial")
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-response), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0("../results/logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0("../results/logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	if(this_aim=="aim2"){
		write_csv(as_tibble(data_new), "../results/training_hai_only_data_aim3.csv")
		writeLines(as.character(as.numeric(data_multiplex_test$infection)), "../results/response_aim3.txt")
	}

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0("../results/logistic_prediction_", this_model, ".txt"))
}

# A feature selected model
## variable selected model 1
for (this_aim in "aim2") { # antibody alone
	# this_aim="aim2"
	# this_aim="aim3"
	if(this_aim=="aim2"){ 
		data_multiplex <- data_multiplex_aim2
		data_multiplex_test <- data_multiplex_aim3
		this_model <- "model_1_select"
	} else {
		data_multiplex <- data_multiplex_aim3
		data_multiplex_test <- data_multiplex_aim2
		this_model <- "model_2_select"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	
	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	model_full <- glm(response~., data=data_combined, family="binomial")
	model <- model_full %>% stepAIC(trace = FALSE)
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-response), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0("../results/logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0("../results/logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	
	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0("../results/logistic_prediction_", this_model, ".txt"))
}

## variable selected model 3
for (this_aim in "aim2") { # antibody + HAI
	# this_aim="aim2"
	# this_aim="aim3"
	if(this_aim=="aim2"){
		data_multiplex <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_3_select"
	} else {
		data_multiplex <- left_join(data_multiplex_aim3, data_clinical_aim3 %>% dplyr::select(sample, contains("pre infection HAI")))
		data_multiplex_test <- left_join(data_multiplex_aim2, data_clinical_aim2 %>% dplyr::select(sample, contains("pre infection HAI")))
		this_model <- "model_4_select"
	}
	names(data_multiplex) <- gsub(" ", "_", names(data_multiplex))
	names(data_multiplex_test) <- gsub(" ", "_", names(data_multiplex_test))

	data_multiplex$infection <- grepl("S1", data_multiplex$group, fixed=T)
	data_multiplex_test$infection <- grepl("S1", data_multiplex_test$group, fixed=T)

	# sort(names(data_multiplex))
	# sort(names(data_multiplex_test))
	intersect(names(data_multiplex), names(data_multiplex_test))
	setdiff(names(data_multiplex), names(data_multiplex_test))

	shared_columns <- intersect(names(data_multiplex), names(data_multiplex_test))
	data_multiplex <- data_multiplex %>% dplyr::select(all_of(shared_columns))
	data_multiplex_test <- data_multiplex_test %>% dplyr::select(all_of(shared_columns))
	
	stopifnot(all(names(data_multiplex_test)==names(data_multiplex)))

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	# data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))	
	model_full <- glm(response~., data=data_combined, family="binomial")
	model <- model_full %>% stepAIC(trace = FALSE)
	
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-`response`), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0("../results/logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0("../results/logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	
	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0("../results/logistic_prediction_", this_model, ".txt"))

}

# plot importance
files_all <- list.files("../results/", "_imp_", full.names=T)
df_imps <- lapply(files_all, function(x){
	tmp <- read_csv(x)
	tmp$file <- x
	tmp
})
df_imps <- bind_rows(df_imps)
df_imps$model <- gsub("../results//logistic_imp_", "", df_imps$file, fixed=T)
df_imps$model <- gsub(".csv", "", df_imps$model, fixed=T)
df_imps$model <- gsub("_", " ", df_imps$model, fixed=T)
df_imps$positive_correlation <- ifelse(df_imps$coefficients>0, "positive", "negative")
df_imps$p_value_text <- paste0("p=", round(df_imps$p_value,3))

df_imps$variable <- gsub("`", "", df_imps$variable, fixed=T)

plots <- c()
for(this_model in unique(df_imps$model)){
	# this_model <- df_imps$model[1]
	df_tmp <- df_imps %>% filter(model %in% this_model)
	df_tmp <- df_tmp %>% arrange(desc(importance)) %>% filter(row_number()<=10)
	df_tmp$variable <- factor(df_tmp$variable, levels=df_tmp$variable)
	p <- ggplot(df_tmp, aes(y=variable))+ geom_col(aes(x=importance, fill=positive_correlation), color="black") + shadowtext::geom_shadowtext(aes(x=importance+0.1, label=p_value_text), bg.color="white", color="black", bg.r=0.1, just="left") +scale_fill_manual("Correlation", values=c("red","blue")) + xlim(c(0, 1.3)) + ggtitle(this_model) + xlab("Importance")+ylab("Variables")+theme_minimal()+theme(legend.position="bottom")
	plots <- c(plots, list(p))
	ggsave(plot=p, paste0("../results/plot_importance_top10_", this_model, ".pdf"), width=5, height=7)
}

files_all
p_out <- plots[[2]] + plots[[5]] + plots[[7]] + plot_layout(ncol=1, height=c(8,8,2), guides='collect') &  theme(legend.position='bottom')
ggsave(plot=p_out, paste0("../results/fig7abc.pdf"), width=8, height=10)