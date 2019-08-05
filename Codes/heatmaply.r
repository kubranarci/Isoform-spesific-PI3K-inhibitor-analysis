#heatmaply function

heatmaply(matrix,
          margins = c(1, 400, 1, 1), 
          dedogram= "raw",
          scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", high = "red", midpoint = 0),
          grid_gap = 0.1,
          fontsize_row = 6,dendrogram = TRUE,
          fontsize_col = 8,Colv=FALSE,
          showticklabels = c(TRUE, TRUE),
          title="Heatmap",
          file= "Heatmap_for_integrative_analysis.html")
# In our study, input matrix arranged such: the sample logFC values were merged per HCC line, row names of it was ensemble gene ids and if necessary, filtrations were performed after samples were merged and before heatmap analysis. 
