library(tidyverse)
library(Rphenograph)
library(Rtsne)
library(ggpubr)

# this is toy script for the CRUK 2019 cancer immunology conference
# it is fairly simplified, and skips a lot of good practice for
# the purpose of illustration
# 
# Duncan Murray

patient <- read_csv("../data/CTCL_patient_data.csv")
biopsy <- read_csv("../data/CTCL_biopsy_data.csv")
flow <- read_csv("../data/CTCL_flow_data.csv")

# run tsne
# currently 4m30s on 5000 cells max per experiment
# 1000 cells: 1m
flow_tsne <-
flow %>%
	select(CD2, CD3, CD4, CD5, CD7, CD8, CD56, TCRAB, TCRGD) %>%
	Rtsne(num_threads = 4, verbose = TRUE) %>%
	pluck("Y") %>%
	cbind(flow, TSNE = .) 

# plot tsne
flow_tsne %>%
	ggplot(aes(x = TSNE.1, y = TSNE.2))+
	geom_point()

# remove some outliers to improve plotting
flow_filtered <-
flow_tsne %>%
	gather(key = "marker", value = "expression",
				 -lesionnumber, -disease, -tissue, -cell,
				 -TSNE.1, -TSNE.2) %>%
	group_by(marker) %>%
	mutate(outlier = expression < quantile(expression, 0.01) |
				 expression > quantile(expression, 0.99)) %>%
	ungroup %>% group_by(cell) %>%
	mutate(outlier_present = sum(outlier)) %>%
	filter(outlier_present == 0) %>% #cell?
	select(-outlier, -outlier_present) %>% ungroup

# plot tsnes
flow_filtered %>%
	split(.$marker) %>%
	map2(names(.), ~ggplot(data = .x,
												 aes(x = TSNE.1, y = TSNE.2, colour = expression))+
			geom_point(size = 0.5, alpha = 0.6)+
			labs(colour = .y)+
			scale_colour_gradient2()) %>%
	ggarrange(plotlist = .)

# cluster cells
flow_wide <-
flow_filtered %>%
	spread(key = marker, value = expression) 

flow_clustered <-
flow_wide %>%
	select(CD2, CD3, CD4, CD5, CD7, CD8, CD56, TCRAB, TCRGD) %>%
	Rphenograph %>%
	pluck(2, membership) %>%
	as.vector %>% as.factor %>%
	cbind(flow_wide, louvain = .)

flow_clustered <-
flow_clustered %>%
	select(CD2, CD3, CD4, CD5, CD7, CD8, CD56, TCRAB, TCRGD) %>%
	kmeans(centers = 10) %>%
	pluck("cluster") %>%
	as.factor %>%
	cbind(flow_clustered, kmean = . ) 

flow_clustered %>%
	select(lesionnumber, cell, louvain, kmean) %>%
	head

# plot clusters
flow_clustered %>%
	ggplot(aes(x = TSNE.1, y = TSNE.2, colour = louvain))+
	geom_point()

# plot clusters
flow_clustered %>%
	ggplot(aes(x = TSNE.1, y = TSNE.2, colour = kmean))+
	geom_point()

# plot clusters by sample
flow_clustered %>%
	group_by(lesionnumber, kmean) %>%
	tally() %>%
	group_by(lesionnumber) %>%
	mutate(percent = 100 * (n / sum(n))) %>%
	ggplot(aes(x = kmean))+
	geom_col(aes(y = percent, fill = kmean))+
	facet_wrap(lesionnumber ~ ., ncol = 5)

# get cluster percentages for samples
biopsy_cluster_results <-
flow_clustered %>%
	group_by(lesionnumber, kmean) %>%
	tally() %>%
	group_by(lesionnumber) %>%
	mutate(percent = 100 * (n / sum(n))) %>%
	select(-n) %>%
	pivot_wider(names_from = kmean, values_from = percent,
							names_prefix = "cluster_")

# join data with biopsy and survival data
biopsy_data <-
left_join(biopsy_cluster_results, biopsy) %>%
	left_join(patient) %>%
	select(-samplenumber, -biopsydate, -sampletype) 

save(biopsy_data, file = "../data/biopsy_data.Rdata")
save(flow_clustered, file = "../data/flow_data_clustered.Rdata")
