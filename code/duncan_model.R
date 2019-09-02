library(tidyverse)
library(ggpubr)
library(glmnet)
library(survival)

load("../data/biopsy_data.Rdata")
load("../data/flow_data_clustered.Rdata")

# reaquaint yourself with the clusters
plots <-
flow_clustered %>% as.tibble %>%
	gather(key = marker, value = expression, -lesionnumber, -disease, -tissue, -cell, -TSNE.1, -TSNE.2, -louvain, -kmean) %>%
	split(.$marker) %>%
	map2(names(.), ~ggplot(data = .x,
												 aes(x = TSNE.1, y = TSNE.2, colour = expression))+
			geom_point(size = 0.5, alpha = 0.6)+
			labs(colour = .y)) 

kmean_plot <-
flow_clustered %>% as.tibble %>%
	ggplot(aes(x = TSNE.1, y = TSNE.2, colour = kmean))+
	geom_point()

# Save image
png("../figures/tsne.png", width = 8*2*1.2, height = 8*2, units = "in", res = 100)
ggarrange(plotlist = c(list(kmean_plot), plots))
dev.off()

# Filter down to disease subtype and remove missing survival data
z <-
biopsy_data %>%
	filter(!is.na(died)) %>%
	filter(diagnosis == "MF") %>%
	mutate(fu = lastfu - diagdate) 

# Get X values
x <-
z %>%
	select_at(vars(starts_with("cluster"))) 

# Fill NAs
x[is.na(x)] <- 0

# Get Y
y <- Surv(z$fu, z$died)

# Run LASSO
cv.fit <- cv.glmnet(as.matrix(x), y, family = "cox")
fit <- glmnet(as.matrix(x), y, family = "cox")

plot(cv.fit)
plot(fit, label = TRUE)

cv.fit$lambda.min
cv.fit$lambda.1se

coefficients <- coef(fit, s = cv.fit$lambda.min)
coefficients <- coef(fit, s = cv.fit$lambda.1se)
active.index <- which(coefficients != 0)
active.coefficients <- coefficients[active.index]

predict(fit, type = 'coefficients')

coxph(Surv(z$fu, z$died) ~ x$cluster_2 + x$cluster_4 + x$cluster_10)

coxph(Surv(z$fu, z$died) ~ x$cluster_7)

x$cluster_15

cluster10_fit <- survfit(Surv(z$fu, z$died) ~ x$cluster_10 < median(x$cluster_10)) 
survdiff(Surv(z$fu, z$died) ~ x$cluster_10 < median(x$cluster_10)) 

plot(cluster10_fit,
		 col = c("blue", "red"))
legend(6000, 0.8, c("a", "b"), col = c("blue", "red"),
			 lty = 1)

cluster3_fit <- survfit(Surv(z$fu, z$died) ~ x$cluster_3 < 1.50) 
survdiff(Surv(z$fu, z$died) ~ x$cluster_3 < median(x$cluster_3)) 

plot(cluster3_fit,
		 col = c("blue", "red"))
legend(6000, 0.9, c("high CD7-", "low CD7-"), col = c("blue", "red"),
			 lty = 1)
