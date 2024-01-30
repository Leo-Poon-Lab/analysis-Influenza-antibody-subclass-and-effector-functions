deal_with_colnames <- function(tibble_i){
	# tibble_i <- data_multiplex_aim2
	# tibble_i <- data_multiplex_aim3
	
	colnames_i <- unlist(tibble_i[1,])
	renamed_t <- gsub(" \\(\\d+\\)$", "", unlist(tibble_i[1,]))
	names(renamed_t) <- colnames_i
	renamed_t <- gsub("A/", "", renamed_t)
	renamed_t <- gsub("A/", "", renamed_t)
	renamed_t <- gsub("B/", "", renamed_t)
	renamed_t <- gsub(" \\D+/\\d+/", "-", renamed_t)
	renamed_t <- gsub(" \\D+/", "-", renamed_t)
	stopifnot(length(unique(colnames_i))==length(unique(renamed_t)))

	check_1 <- grep("A/Brisbane/59/2007", names(renamed_t), fixed=T)
	renamed_t[check_1] <- paste0("vaxx ", renamed_t[check_1])
	check_2 <- grep("A/California/04/2009", names(renamed_t), fixed=T)
	renamed_t[check_2] <- paste0("pdm ", renamed_t[check_2])
	check_2_5 <- grep("A/California/07/2009", names(renamed_t), fixed=T)
	renamed_t[check_2_5] <- paste0("pdm ", renamed_t[check_2_5])
	check_3 <- grep("chicken", tolower(names(renamed_t)), fixed=T)
	renamed_t[check_3] <- paste0("av. ", renamed_t[check_3])
	check_4_1 <- grep(" b/", tolower(names(renamed_t)), fixed=T)
	renamed_t[check_4_1] <- paste0("seas. B ", renamed_t[check_4_1])
	check_4_2 <- grep("^b/", tolower(names(renamed_t)))
	renamed_t[check_4_2] <- paste0("seas. B ", renamed_t[check_4_2])
	check_5 <- grep(" a/", tolower(names(renamed_t)), fixed=T)
	check_5 <- check_5[!check_5 %in% c(check_1, check_2, check_2_5, check_3, check_4_1, check_4_2)]
	renamed_t[check_5] <- paste0("seas. ", renamed_t[check_5])
	# unique(unlist(tibble_i[1,]))
	# unique(renamed_t)

	idx_torm <- 0
	colnames_i <- tolower(colnames_i)
	if("type" %in% colnames_i){idx_torm <- c(idx_torm, which(colnames_i=="type"))}
	if("well" %in% colnames_i){idx_torm <- c(idx_torm, which(colnames_i=="well"))}
	if(any(grepl("anti-igg-fab", colnames_i))){idx_torm <- c(idx_torm, grep("anti-igg-fab", colnames_i))}
	if(any(grepl("tetanus", colnames_i))){idx_torm <- c(idx_torm, grep("tetanus", colnames_i))}
	if(any(grepl("b/hk/330/2001", colnames_i))){idx_torm <- c(idx_torm, grep("b/hk/330/2001", colnames_i))}
	if(any(grepl("b/brisbane/2008", colnames_i))){idx_torm <- c(idx_torm, grep("b/brisbane/2008", colnames_i))}
	if(sum(idx_torm>0)>0){
		tibble_i <- tibble_i[,-idx_torm[idx_torm!=0]]
		renamed_t <- renamed_t[-idx_torm[idx_torm!=0]]
	}
	
	antibody_i <- names(tibble_i)
	antibody_i[grepl("^\\.", antibody_i)] <- ""
	antibody_i <- gsub("...\\d+$", " ", antibody_i)
  antibody_i <- gsub("FcrR", "Fc\u03B3R", antibody_i)
	real_names <- paste0(antibody_i, renamed_t)
	names(tibble_i) <- real_names
	tibble_i[-1,]
}