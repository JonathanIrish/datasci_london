library(tidyverse)
library(Rtsne)
library(Rphenograph)
library(kohonen)
library(ggpubr)
library(scales)

# look at the data
iris %>%
	as_tibble() %>%
	print()

# single dimension - look at sepal width and length
iris %>%
	mutate(flower = 1:dim(iris)[1]) %>%
	select(-Petal.Length, -Petal.Width) %>%
	gather(key = "measurement", value = "value", -Species, -flower) %>%
	ggplot(aes(x = measurement, y = value))+
	geom_dotplot(binaxis = 'y',
							 method = 'histodot',
							 stackdir = 'center',
							 position = position_jitter(0.01),
							 alpha = 0.7,
							 dotsize = 0.4)

# single dimensions - histogram
iris %>%
	mutate(flower = 1:dim(iris)[1]) %>%
	select(-Petal.Length, -Petal.Width) %>%
	gather(key = "measurement", value = "value", -Species, -flower) %>%
	ggplot(aes(value, fill = measurement, colour = measurement))+
	geom_density(alpha = 0.4)


# single dimensions and paired stripplots
iris %>%
	mutate(flower = 1:dim(iris)[1]) %>%
	select(-Petal.Length, -Petal.Width) %>%
	gather(key = "key", value = "value", -Species, -flower) %>%
	ggplot(aes(x = key, y = value))+
	geom_dotplot(binaxis = 'y',
							 method = 'histodot',
							 stackdir = 'center',
							 position = position_jitter(0.01),
							 alpha = 0.7,
							 dotsize = 0.4)+
	geom_line(aes(group = flower),
						alpha = 0.5)

# two dimensions
iris %>%
	mutate(flower = 1:dim(iris)[1]) %>%
	ggplot(aes(x = Sepal.Length, y = Sepal.Width))+
	geom_point(alpha = 0.7,
						 size = 5.5)

# removing identical results to simplify dimension reduction/clustering
iris2 <- 
iris %>%
	unique() 

# collect dimension reduction and cluster into iris3
iris3 <- iris2

tsn <- 
iris2[,1:4] %>%
	as.matrix() %>%
	Rtsne(max_iter = 5000,
				perplexity = 10)

iris3$x_tsne <- tsn$Y[,1]
iris3$y_tsne <- tsn$Y[,2]
	
iris3 %>%
	ggplot(aes(x = x_tsne, y = y_tsne))+
	geom_point(alpha = 0.7,
						 size = 5.5)

pheno <-
iris2 %>%
	select(-Species) %>%
	Rphenograph()

som <-
iris2 %>%
	select(-Species) %>%
	as.matrix %>%
	som(grid = somgrid(xdim = 3,
										 ydim = 2))

k <-
iris2 %>%
	select(-Species) %>%
	kmeans(3)

h <-
	hclust(dist(as.matrix(iris2[,1:4])))

hc <- cutree(h, 3)

p <- svd(t(iris2[1:4]))

iris3$x_pc <- p$d[1] * p$v[,1]
iris3$y_pc <- p$d[2] * p$v[,2]

iris3 %>%
	ggplot(aes(x = x_pc, y = y_pc))+
	geom_point(alpha = 0.7,
						 size = 5.5)

iris3$som <- as.factor(som$unit.classif)
iris3$pheno <- as.factor(pheno[[2]]$membership)
iris3$kmeans <- as.factor(k$cluster)
iris3$hclust <- as.factor(hc)

iris3 %>%
	head

iris3 %>%
	pivot_longer(cols = c(som, pheno, kmeans, hclust, Species),
							 names_to = "cluster_type",
							 values_to = "cluster") %>%
	rename(x_Sepal = Sepal.Length,
				 y_Sepal = Sepal.Width,
				 x_Petal = Petal.Length,
				 y_Petal = Petal.Width) %>%
	pivot_longer(cols = matches("([x,y])_"),
							 names_to = c(".value", "axis_type"),
							 names_sep = "_"
							 ) %>%
	group_by(axis_type) %>%
	mutate(x = rescale(x)) %>%
	mutate(y = rescale(y)) %>%
	ggplot(aes(x = x, y = y, colour = cluster))+
	geom_point(alpha = 0.7,
						 size = 2)+
	facet_grid(axis_type ~ cluster_type, scales = "free")



iris3 %>%
	ggplot(aes(x = Sepal.Width, y = Sepal.Length, colour = pheno))+
	geom_point(aes(colour = pheno),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = Sepal.Width, y = Sepal.Length, colour = kmeans))+
	geom_point(aes(colour = kmeans),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = Sepal.Width, y = Sepal.Length, colour = Species))+
	geom_point(aes(colour = Species),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = Sepal.Width, y = Sepal.Length, colour = hclust))+
	geom_point(aes(colour = hclust),
						 alpha = 0.7,
						 size = 5.5,
						 )

iris3 %>%
	ggplot(aes(x = tsne1, y = tsne2, colour = pheno))+
	geom_point(aes(colour = pheno),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = tsne1, y = tsne2, colour = kmeans))+
	geom_point(aes(colour = kmeans),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = tsne1, y = tsne2, colour = Species))+
	geom_point(aes(colour = Species),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = tsne1, y = tsne2, colour = hclust))+
	geom_point(aes(colour = hclust),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = pc1, y = pc2, colour = pheno))+
	geom_point(aes(colour = pheno),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = pc1, y = pc2, colour = kmeans))+
	geom_point(aes(colour = kmeans),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = pc1, y = pc2, colour = Species))+
	geom_point(aes(colour = Species),
						 alpha = 0.7,
						 size = 5.5,
						 )
iris3 %>%
	ggplot(aes(x = pc1, y = pc2, colour = hclust))+
	geom_point(aes(colour = hclust),
						 alpha = 0.7,
						 size = 5.5,
						 )

iris3 %>%
	ggplot(aes(x = pc1, y = -pc2, colour = som))+
	geom_point(aes(colour = som),
						 alpha = 0.7,
						 size = 5.5,
						 )

