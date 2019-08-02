####Sig digits function ######
sig <- function(x,dig=3,maxex=5) {
  namez <- names(x)
  x <- formatC(signif(x,digits=dig), digits=dig, format='g', flag='#')
  
  if(dig!=maxex) {
    ex <- "([]*[0-9]\\.[0-9]+)e([+][0-9]{2})"
    subit <- grepl(ex,x,perl=TRUE)
    a <- as.numeric(gsub(ex, "\\1", x))
    b <- as.numeric(gsub(ex, "\\2", x))
    subit <- subit & abs(b) < maxex
    x <- ifelse(subit,formatC(signif(as.numeric(x),digits=dig),digits=dig, format="fg",flag="#"),x)
  }
  x <- gsub("\\.$", "", x, perl=TRUE)
  names(x) <- namez
  return(x)
}