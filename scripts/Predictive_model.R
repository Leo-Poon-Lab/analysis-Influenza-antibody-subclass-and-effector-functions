# In this script, we try to build prediction models seperately for aim2 and aim3 data,
# and we test the model with the data from the other aim.

library(caret)
library(tidyverse)
library(readxl)
library(pcaMethods)
library(patchwork)
library(MASS)
library(shadowtext)

# read data
source(here::here("scripts/helper/deal_with_colnams.R"))
source(here::here("scripts/helper/prepare_multiplex_data.R"))
data_all <- prepare_data()
data_multiplex_aim2 <- data_all[[1]]
data_clinical_aim2 <- data_all[[2]]
data_multiplex_aim3 <- data_all[[3]]
data_clinical_aim3 <- data_all[[4]]

elastic_net_ranking_aim2 <- readxl::read_excel(here::here("results/elastic_net_plsda_results/aim2/model_3_en_plsda/ranking_from_elastic_net_prediction.xlsx"))
elastic_net_ranking_aim3 <- readxl::read_excel(here::here("results/elastic_net_plsda_results/aim3/model_4_en_plsda/ranking_from_elastic_net_prediction.xlsx"))

get_shared_top_N_variables <- function(data_1, data_2, N=10){
  # data_1 = elastic_net_ranking_aim2
  # data_2 = elastic_net_ranking_aim3
  # N = 10
  stopifnot(all(data_1[[1]] %in% data_2[[1]]))

  # we choose the top 10 shared variables between two models by comparing ranked variables from both models' lists. Each list contains 24 variables. In each round (we do it sequentially from rank 1 to rank 24), I will compare the same-ranked variables from both lists. If they match, the variable is selected. If not, the variables are added to a waiting list. In subsequent rounds, new variables are compared to each other and those in the waiting list. Once a match is found, the variable is selected. This process continues until 10 variables are selected.
  waiting_list <- c()
  selected_list <- c()
  for (i in 1:24){
    if (length(selected_list)<N){
      if (data_1[[1]][i]==data_2[[1]][i]){
        selected_list <- c(selected_list, data_1[[1]][i])
      } else {
        if (data_1[[1]][i] %in% waiting_list){
          selected_list <- c(selected_list, data_1[[1]][i])
          waiting_list <- waiting_list[!waiting_list %in% data_1[[1]][i]]
        } else if (data_2[[1]][i] %in% waiting_list){
          selected_list <- c(selected_list, data_2[[1]][i])
          waiting_list <- waiting_list[!waiting_list %in% data_2[[1]][i]]
        } else {
          waiting_list <- c(waiting_list, data_1[[1]][i], data_2[[1]][i])
        }
      }
    } else {
      break()
    }
  }
  selected_list
}

selected_shared_variables_15 <- get_shared_top_N_variables(elastic_net_ranking_aim2, elastic_net_ranking_aim3, N=15)
selected_shared_variables_15 <- gsub("`", "", selected_shared_variables_15)
selected_shared_variables <- selected_shared_variables_15[1:10]
writeLines(selected_shared_variables, here::here("results/logistic_reg/top10_selected_shared_variables.txt"))
selected_shared_variables_non_HAI <- selected_shared_variables_15[!grepl("pre_infection_HAI", selected_shared_variables_15)][1:10]

data_multiplex_aim2 <- data_multiplex_aim2 %>% mutate_at(vars(!one_of("group", "sample")), as.numeric)
data_multiplex_aim3 <- data_multiplex_aim3 %>% mutate_at(vars(!one_of("group", "sample")), as.numeric)
data_clinical_aim2$sample <- as.character(data_clinical_aim2$sample)
data_clinical_aim3$sample <- as.character(data_clinical_aim3$sample)

# convert the column names to be the same between aim2 and aim3 data in this analysis
names(data_multiplex_aim2) <- gsub(" sH1", " vaxx sH1-2007", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub(" pH1", " pH1-2009", names(data_multiplex_aim2), fixed=T)
names(data_multiplex_aim2) <- gsub(" H3-2007", " sH3-2007", names(data_multiplex_aim2), fixed=T)

# build model for aim2 and aim3
all_aims <- c("aim2", "aim3")
data_multiplex_aim3 <- data_multiplex_aim3 %>% filter(timepoint==2) # keep only the d2 data for comparison between the aim2 and aim3 data

dir.create(here::here("results/logistic_reg"))

for (this_aim in all_aims) { # antibody alone
	# this_aim="aim2"
	# this_aim="aim3"
  out_path <- here::here("results/logistic_reg/")
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
	data_m[data_m<0] <- 0
	data_m <- log10(data_m+1) # log transformation
	data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

  data_m <- as_tibble(data_m)
  stopifnot(all(selected_shared_variables_non_HAI %in% names(as_tibble(data_m))))
  data_m <- data_m %>% dplyr::select(all_of(selected_shared_variables_non_HAI))

	if(this_aim=="aim2"){
		write_csv(data_m, here::here("results/logistic_reg/training_data_aim2.csv"))
		writeLines(as.character(as.numeric(data_multiplex$infection)), paste0(out_path, "response_aim2.txt"))
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
	writeLines(tmp, paste0(out_path, "logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0(out_path, "logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
  data_new <- data_new %>% dplyr::select(all_of(selected_shared_variables_non_HAI))

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0(out_path, "logistic_prediction_", this_model, ".txt"))
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

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection, -contains("pre_infection_HAI")) %>% mutate_all(as.numeric)
	
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	data_m[data_m<0] <- 0
	data_m <- log10(data_m+1) # log transformation
	data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance
	data_m <- cbind(data_m, data_multiplex %>% dplyr::select(contains("pre_infection_HAI")))

  data_m <- as_tibble(data_m)
  stopifnot(all(selected_shared_variables %in% names(as_tibble(data_m))))
  data_m <- data_m %>% dplyr::select(all_of(selected_shared_variables))

	data_combined <- bind_cols(data_m,response=as.numeric(data_multiplex$infection))
	model <- glm(response~., data=data_combined, family="binomial")
	# summary(model)
	df_imp <- varImp(model)/max(varImp(model))
	df_imp <- tibble(variable=rownames(df_imp), importance=df_imp$Overall) %>% arrange(desc(importance))
	tmp <- summary(model)
	df_imp <- left_join(df_imp, tibble(variable=rownames(tmp$coefficients), coefficients=tmp$coefficients[,1], p_value=tmp$coefficients[,4]))
	probabilities <- model %>% predict(data_combined %>% dplyr::select(-`response`), type = "response")
	test_self <- (probabilities>0.5) == data_combined$response
	mean(test_self)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_combined$response==1)))
	writeLines(tmp, paste0(out_path, "logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0(out_path, "logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
  data_new <- data_new %>% dplyr::select(all_of(selected_shared_variables))

	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	if(this_aim=="aim2"){
		write_csv(as_tibble(data_new), paste0(out_path, "training_data_hai_antibody_aim3.csv"))
		writeLines(as.character(as.numeric(data_multiplex_test$infection)), paste0(out_path, "response_aim3.txt"))
	}

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0(out_path, "logistic_prediction_", this_model, ".txt"))

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

	if(this_aim=="aim2"){
		write_csv(as_tibble(data_m), paste0(out_path, "training_data_hai_only_aim2.csv"))
		writeLines(as.character(as.numeric(data_multiplex$infection)), paste0(out_path, "response_aim2.txt"))
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
	writeLines(tmp, paste0(out_path, "logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0(out_path, "logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	if(this_aim=="aim2"){
		write_csv(as_tibble(data_new), paste0(out_path, "training_hai_only_data_aim3.csv"))
		writeLines(as.character(as.numeric(data_multiplex_test$infection)), paste0(out_path, "response_aim3.txt"))
	}

	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0(out_path, "logistic_prediction_", this_model, ".txt"))
}

# A feature selected model using stepAIC
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
	data_m[data_m<0] <- 0
	data_m <- log10(data_m+1) # log transformation
	data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance

  data_m <- as_tibble(data_m)
  stopifnot(all(selected_shared_variables_non_HAI %in% names(as_tibble(data_m))))
  data_m <- data_m %>% dplyr::select(all_of(selected_shared_variables_non_HAI))

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
	writeLines(tmp, paste0(out_path, "logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0(out_path, "logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
  data_new <- data_new %>% dplyr::select(all_of(selected_shared_variables_non_HAI))
	
	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0(out_path, "logistic_prediction_", this_model, ".txt"))
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

	data_m <- data_multiplex %>% dplyr::select(-group,-sample,-timepoint,-infection, -contains("pre_infection_HAI")) %>% mutate_all(as.numeric)
	check <- apply(data_m, 1, function(x){any(is.na(x))})
	data_multiplex <- data_multiplex[!check,]
	data_m <- data_m[!check,]
	data_m[data_m<0] <- 0
	data_m <- log10(data_m+1) # log transformation
	data_m <- apply(data_m, 2, prep, scale="uv") # scale the data, with unit variance
	data_m <- cbind(data_m, data_multiplex %>% dplyr::select(contains("pre_infection_HAI")))

  data_m <- as_tibble(data_m)
  stopifnot(all(selected_shared_variables %in% names(as_tibble(data_m))))
  data_m <- data_m %>% dplyr::select(all_of(selected_shared_variables))

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
	writeLines(tmp, paste0(out_path, "logistic_accuracy_", this_model, ".txt"))
	write_csv(df_imp, paste0(out_path, "logistic_imp_", this_model, ".csv"))

	data_new <- data_multiplex_test %>% dplyr::select(-group,-sample,-timepoint,-infection) %>% mutate_all(as.numeric)
	check <- apply(data_new, 1, function(x){any(is.na(x))})
	data_multiplex_test <- data_multiplex_test[!check,]
	data_new <- data_new[!check,]
  data_new <- data_new %>% dplyr::select(all_of(selected_shared_variables))
	# data_new <- apply(data_new, 2, prep, scale="uv") # scale the data, with unit variance
	
	probabilities <- model %>% predict(as_tibble(data_new), type = "response")
	test_new <- (probabilities>0.5) == data_multiplex_test$infection
	mean(test_new)
	tmp <- capture.output(confusionMatrix(data=factor(probabilities>0.5), factor(data_multiplex_test$infection)))
	writeLines(tmp, paste0(out_path, "logistic_prediction_", this_model, ".txt"))

}

# plot importance
files_all <- list.files(out_path, "_imp_", full.names=T)
df_imps <- lapply(files_all, function(x){
	tmp <- read_csv(x)
	tmp$file <- x
	tmp
})
df_imps <- bind_rows(df_imps)
df_imps$model <- gsub(".+/logistic_imp_", "", df_imps$file)
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
	# p <- ggplot(df_tmp, aes(y=variable))+ geom_col(aes(x=importance, fill=positive_correlation), color="black") + shadowtext::geom_shadowtext(aes(x=importance+0.1, label=p_value_text), bg.color="white", color="black", bg.r=0.1, just="left", size=2) +scale_fill_manual("Correlation", values=c("red","blue")) + xlim(c(0, 1.3)) + ggtitle(this_model) + xlab("Importance")+ylab("Variables")+theme_minimal()+theme(legend.position="bottom")
	p <- ggplot(df_tmp, aes(y=variable))+ geom_col(aes(x=importance), fill="grey", color="black") + shadowtext::geom_shadowtext(aes(x=importance+0.1, label=p_value_text), bg.color="white", color="black", bg.r=0.1, just="left", size=3) + xlim(c(0, 1.3)) + ggtitle(this_model) + xlab("Importance")+ylab("Variables")+theme_minimal()+theme(legend.position="bottom")
	plots <- c(plots, list(p))
	ggsave(plot=p, paste0(out_path, "plot_importance_top10_", this_model, ".pdf"), width=5, height=7, device = cairo_pdf)
}

p_out <- plots[[2]] + plots[[1]] + plots[[5]] + plots[[4]] + plots[[7]] + plot_spacer() + plot_layout(ncol=2, height=c(8,8,2), guides='collect') &  theme(legend.position='bottom')
ggsave(plot=p_out, paste0(out_path, "fig7abcde.pdf"), width=8, height=8, device = cairo_pdf)
