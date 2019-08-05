#script to plot correlations

corrx <- cor(matrix)
corrplot(corrx, method="number", type = "upper",  order="hclust", col=col, tl.col="black", tl.srt=45,)
