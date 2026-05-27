library(ggplot2)

pca <- read.table("PCA.eigenvec.txt", header=FALSE)

colnames(pca)[1:2] <- c("FID", "IID")

colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_minimal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA - Numenius arquata")


meta <- read.table("meta_2.txt", header = TRUE)

pca <- merge(pca, meta, by = "IID")

ggplot(pca, aes(x = PC1, y = PC2, color = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA - Numenius arquata and Numenius phaeopus",
    x = "PC1",
    y = "PC2",
    color = "Population"
  )

ggplot(pca, aes(x = PC1, y = PC2, color = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA - Numenius arquata",
    x = "PC1",
    y = "PC2",
    color = "Population"
  )

ggplot(pca, aes(PC1, PC2, color = Pop)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

