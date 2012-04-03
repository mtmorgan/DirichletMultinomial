##  utility

.gradient <-                   # RColorBrewer::brewer.pal(9, "YlOrRd")
    c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
      "#FC4E2A", "#E31A1C", "#BD0026", "#800026")

.divergent <-                  # RColorBrewer::brewer.pal(9, "RdYlBu")
    c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
      "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")

.qualitative <-                # RColorBrewer::brewer.pal(10, "Paired")
    c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

csubset <- function(val, x, pheno, cidx=TRUE)
{
    ridx <- pheno %in% val
    if (!cidx)
        cidx <- colSums(x[ridx,]) != 0
    x[ridx, cidx]
}
