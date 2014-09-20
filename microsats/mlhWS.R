mlhWS <- function (input, na.string, n.digits) 
# different from mlh() since it loads from workspace, doesn´t save as text file and 
# just calculates SH
        
{
        g <- input
        chkdata(g[, 2:ncol(g)])
        h <- as.data.frame(array(dim = c(nrow(g), 2)), stringsAsFactors = FALSE)
        h[, 1] <- g[, 1]
        h[, 2] <- sh(g[, 2:ncol(g)])
        colnames(h) <- c("ID", "SH")
        h
}