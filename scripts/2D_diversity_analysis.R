# Quantify antigenic diversity on different RBDs for different groups,
# and comparison between groups.

library(readxl)
library(tidyverse)
library(boot)
library(parallel)
library(scico)
library(ggrepel)
library(patchwork)
library(ggplot2)
library(ggtext)

# read data
files_input <- list.files(here::here("data/data_split/"), full.names = T)
files_input <- files_input[!grepl("~$", files_input, fixed=T)]

df_response <- lapply(files_input, function(x){
	antibody <- gsub(".+/", "", x)
	antibody <- gsub(".xlsx", "", antibody)
	tmp <- read_xlsx(x, skip=1)
	tmp$antibody <- antibody
	tmp
})

df_response <- bind_rows(df_response)
df_response <- df_response %>% filter(!is.na(sample))
df_response <- df_response %>% select(-"Tetanus (11)") %>% select(-"anti-IgG-Fab Ab (44)") %>% select(-"HA B/Brisbane/2008 (30)")

antibodys <- unique(df_response$antibody)
(antigens <- names(df_response)[5:22])

samples <- unique(df_response$sample)
types <- sapply(samples, function(x){
	if(grepl("d1$", x)){return("day 0")}
	if(grepl("d2$", x)){return("day 30")}
	if(grepl("d3$", x)){return("1 year")} else {return(NA)}
})
df_response$sample_type <- sapply(df_response$sample, function(x){
	types[names(types)==x]
}, USE.NAMES=F)

df_response_long <- df_response %>% pivot_longer(all_of(antigens), names_to="antigen", values_to="response")
df_response_long <- df_response_long %>% select(-Type, -Well)

# log transformation
do_log_transfromation <- TRUE
if(do_log_transfromation){
df_response_long$response[df_response_long$response<0] <- 0
df_response_long$response <- log10(df_response_long$response+1)
# scaled within antibody group
# df_response_long <- df_response_long %>% group_by(antibody) %>% mutate(response=scale(response))
# scaled within protein group
unique(df_response_long$antigen)
df_response_long <- df_response_long %>% group_by(antigen) %>% mutate(response=scale(response))
# df_response_long$response <- pcaMethods::prep(df_response_long$response, scale= "uv")
}

df_response_long$group_full <- paste0(df_response_long$group, " - ", df_response_long$sample_type)
df_response_long$group_full <- gsub(" - NA", "", df_response_long$group_full)
groups_all <- unique(df_response_long$group_full)

if(do_log_transfromation){
  df_response_long <- df_response_long %>% group_by(antibody, antigen) %>% mutate(neg_average_response=mean(response[group_full=="V0 S1 - day 30"]))
  df_response_long$diff <- (df_response_long$response - df_response_long$neg_average_response)
  df_response_long$Response_type <- ifelse(df_response_long$diff<log(1.2), "negative", "positive") # we manually choose the 20% threshold
} else {
  df_response_long <- df_response_long %>% group_by(antibody, antigen) %>% mutate(neg_average_response=mean(response[group_full=="V0 S1 - day 30"]))
  df_response_long$diff <- (df_response_long$response - df_response_long$neg_average_response)/df_response_long$neg_average_response
  df_response_long$Response_type <- ifelse(df_response_long$diff<0.2, "negative", "positive") # we manually choose the 20% threshold
}
# table(df_response_long$Response_type)

df_response_long <- df_response_long %>% ungroup()

# When analyzing antigenic diversity, we have two questions to be answered:
# 1. Which groups have higher magnitude responses?
# 2. Which groups have boarder antibody responses?

# before starting the analysis, we have to introduce the definition of
# "positive response", which many of the following analyses are based on.
# Only responses >= 95 percentile response of the negative group (V0 S0)
# are treated as "positive responses",
# otherwise the responses are treated as "negative responses".
# There are 21 different types (on 21 different antigens BDs) of "positive responses"
# in this study, but all negative responses (even on different RBDs)
# are classified in one "negative responses" group.


## 1. Which groups have higher magnitude responses?
## 2. Which groups have boarder antibody response?
func_boot <- function(data, indices){
	# df_tmp <- df_inhibition_i
	df_tmp <- data[indices,]

	### 1. Which groups have higher magnitude responses?
	#### we estimated the mean and 95CI of %inhibition of the "all responses"
	mean_res <- mean(df_tmp$diff)
	#### we estimated the mean and 95CI of %inhibition of the "positive responses"
	mean_pos_res <- mean(df_tmp$diff[df_tmp$Response_type=="positive"])
	
	### 2. Which groups have boarder antibody response?
	### https://academic.oup.com/ve/article/5/1/vey041/5304643
	df_tmp$RBD2 <- df_tmp$antigen
	df_tmp$RBD2[df_tmp$Response_type=="negative"] <- "Negatives"
	tmp_freq <- table(df_tmp$RBD2)
	N_total <- nrow(df_tmp)
	tmp_p <- tmp_freq/N_total

	#### we use shannon entropy to estimate the diversity
	hs <- sum(-tmp_p*log(tmp_p))

	#### we also use "antigenic diveristy (pi)" which is simialr to nucleotide diversity pi, to estimate the diversity
	hpi <- 1-(sum(tmp_freq*(tmp_freq-1))/(N_total*(N_total-1)))

	return(c(mean_res, mean_pos_res, hs, hpi))
}

### full spectrum/panel, including Sarbecovirus
num_repeat <- 1000
num_cpus <- 8
set.seed(2022)

df_boot_rst_all <- lapply(c(antibodys, "combined"), function(antibody_i){
	# antibody_i = antibodys[1]
	if(antibody_i == "combined"){
		data_i <- df_response_long
	} else {
		data_i <- df_response_long %>% filter(antibody == antibody_i)
	}

	list_boot_rst_full <- mclapply(groups_all, function(group_i){
		# print(group_i)
		# group_i <- groups_all[1]
		df_inhibition_i <- data_i %>% filter(group_full == group_i)
		rst_boot <- boot(df_inhibition_i, func_boot, R=num_repeat, parallel = "no")
		ci_mean_res <- boot.ci(rst_boot, type=c("perc"), index = 1)
		ci_mean_pos_res <- boot.ci(rst_boot, type=c("perc"), index = 2)
		ci_hs <- boot.ci(rst_boot, type=c("perc"), index = 3)
		ci_hpi <- boot.ci(rst_boot, type=c("perc"), index = 4)
		# str(rst_boot)
		# str(ci_mean_pos_res)
		rst <- c(group_i, apply(rst_boot$t, 2, mean, na.rm=T), apply(rst_boot$t, 2, sd, na.rm=T), ci_mean_res$percent[4:5], ci_mean_pos_res$percent[4:5], ci_hs$percent[4:5], ci_hpi$percent[4:5])
		print(rst)
		rst
	}, mc.cores = num_cpus)

	df_boot_rst <- do.call(rbind, list_boot_rst_full)
	colnames(df_boot_rst) <- c("group", "e_mean_res", "e_mean_pos_res", "e_hs", "e_hpi", "sd_mean_res", "sd_mean_pos_res", "sd_hs", "sd_hpi", "ci_mean_res_low", "ci_mean_res_high", "ci_mean_pos_res_low", "ci_mean_pos_res_high", "ci_hs_low", "ci_hs_high", "ci_hpi_low", "ci_hpi_high")
	df_boot_rst <- as_tibble(df_boot_rst)
	df_boot_rst$antibody <- antibody_i
	write_tsv(df_boot_rst, paste0(here::here("results/df_boot_rst_full_"), antibody_i, ".tsv"))
	df_boot_rst
})


# Draw
### full spectrum/panel, including Sarbecovirus
df_boot_rst <- bind_rows(df_boot_rst_all)
df_boot_rst <- df_boot_rst %>% mutate_at(vars(-group, -antibody), as.numeric)
df_boot_rst$group_sim <- gsub(" - .+$", "", df_boot_rst$group)
df_boot_rst$d_sim <- gsub("V. S. - ", "", df_boot_rst$group)
df_boot_rst$d_sim[df_boot_rst$d_sim=="V0 S0"] <- "day 0"
# colors_t <- scico(length(unique(df_boot_rst$group)), palette = 'batlow')

colors_t <- c("black", "red", "blue", "grey40")[factor(df_boot_rst$group_sim)]
names(colors_t) <- df_boot_rst$group_sim
colors_t <- colors_t[order(names(colors_t))]
colors_t <- colors_t[!duplicated(colors_t)]

timepoints <- c("day 0 = Pre vaccination", "day 30 = Post vaccination/pre infection", "1 year = Post infection")
df_boot_rst$d_sim <- factor(df_boot_rst$d_sim, levels=c("day 0", "day 30", "1 year"))
df_boot_rst$type_sim <- timepoints[df_boot_rst$d_sim]
shape_t <- c(1,16,15)
names(shape_t) <- timepoints

# Two dimensional illustration of the "Response fitness"
sort(unique(df_boot_rst$antibody))
df_boot_rst_for_plot <- df_boot_rst %>% filter(antibody != "combined") 
df_boot_rst_for_plot$antibody_parsed <- factor(df_boot_rst_for_plot$antibody, levels = c("IgG", "IgG1", "IgG3", "IgA1", "IgM", "FcrR2a", "FcrR3a"), labels=c("IgG", "IgG1", "IgG3", "IgA1", "IgM", expression(paste("Fc", gamma, "R2a")), expression(paste("Fc", gamma, "R3a"))))

plot_2d <- function(df, x_var, y_var, color_var, shape_var, annt_var) {
	# df = df_boot_rst_for_plot
	# x_var = "mean_res"
	# y_var = "hpi"
	# color_var = "group_sim"
	# shape_var = "type_sim"
	# annt_var = "group"
  x_lab <- ifelse(grepl("pos", x_var), "Average difference (positive responses)", "Average difference")
	
  y_lab <- ifelse(grepl("hs", y_var), "Shannon entropy", expression("HA cross-reactivity ("~pi~")"))

	x_e <- paste0("e_", x_var)
	y_e <- paste0("e_", y_var)
	x_ci_l <- paste0("ci_", x_var, "_low")
	x_ci_h <- paste0("ci_", x_var, "_high")
	y_ci_l <- paste0("ci_", y_var, "_low")
	y_ci_h <- paste0("ci_", y_var, "_high")
	
  if(do_log_transfromation){
    df[[x_e]] <- sapply(df[[x_e]], function(x){10^x-1})
    df[[x_ci_l]] <- sapply(df[[x_ci_l]], function(x){10^x-1})
    df[[x_ci_h]] <- sapply(df[[x_ci_h]], function(x){10^x-1})
  } 
	ggplot(df)+
		geom_segment(aes_string(x=x_e, xend=x_e, y=y_ci_l, yend=y_ci_h, color=color_var), alpha=0.6, size=0.5)+
		geom_segment(aes_string(x=x_ci_l, xend=x_ci_h, y=y_e, yend=y_e, color=color_var), alpha=0.6, size=0.5)+
		# geom_point(aes_string(x=x_e, y=y_e, color=color_var), alpha=0.8, shape=1, size=2, data=. %>% filter(open_circle))+
		# geom_point(aes_string(x=x_e, y=y_e, color=color_var, fill=color_var), alpha=0.8, shape=16, size=2, data=. %>% filter(!open_circle), show.legend=FALSE)+
		geom_point(aes_string(x=x_e, y=y_e, color=color_var, fill=color_var, shape=shape_var), alpha=0.8, size=2)+
		geom_text_repel(aes_string(x=x_e, y=y_e, color=color_var, label=annt_var), alpha=0.9, bg.color = "white", bg.r = 0.05, size=2, segment.size=0.3, segment.alpha=0.9, segment.color="black", arrow = arrow(length = unit(0.01, "npc")), point.padding=0.28, box.padding=0.45, show.legend=FALSE, max.overlaps = Inf, force=0.5)+
		# scale_fill_manual(name="Group", values=colors_t)+
		scale_color_manual(name="Group", values=colors_t)+
		scale_shape_manual(name="Timepoint", values=shape_t)+
		xlab(x_lab)+
		ylab(y_lab)+
		theme_bw()+
		theme(legend.position = c(0.65, 0.15), legend.direction="horizontal", legend.text = element_markdown())+
		guides(fill="none",
			color=guide_legend(override.aes = list(shape =".")),
			shape=guide_legend(nrow=length(shape_t))
			)+
		facet_wrap(vars(antibody_parsed), scales="free", labeller = label_parsed, )+
		ylim(0,1)+
		NULL
}

p0 <- plot_2d(df = df_boot_rst_for_plot,
	x_var = "mean_res",
	y_var = "hpi",
	color_var = "group_sim",
	shape_var = "type_sim",
	annt_var = "group"
	)
ggsave(here::here("results/2D_pos_res_hpi.pdf"), width=8, height=8, plot=p0)

p1 <- plot_2d(df = df_boot_rst_for_plot,
	x_var = "mean_res",
	y_var = "hs",
	color_var = "group_sim",
	shape_var = "type_sim",
	annt_var = "group"
	)
ggsave(here::here("results/2D_pos_res_hs.pdf"), width=8, height=8)



