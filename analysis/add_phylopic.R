# modifed from Scott Chamberlain's rphylopic package
# https://github.com/sckott/rphylopic

add_phylopic <- function(img, alpha = 0.2, x = NULL, y = NULL, ysize = NULL,
  color = NULL, xy_ratio = 1){

  if (is.null(color)) {
    mat <- matrix(rgb(img[,,1], img[,,2], img[,,3], img[,,4] * alpha), nrow = dim(img)[1])
  } else {
    cols <- col2rgb(color)
    imglen <- length(img[,,1])
    mat <- matrix(ifelse(img[,,4] > 0, rgb(rep(cols[1,1], imglen),
      rep(cols[2,1], imglen),
      rep(cols[3,1], imglen),
      img[,,4]*255*alpha, maxColorValue = 255),
      rgb(rep(1, imglen),
        rep(1, imglen),
        rep(1, imglen),
        img[,,4]*alpha)), ## make background white for devices
      nrow = dim(img)[1])       ## that do not support alpha channel
  }
  if (!is.null(x) && !is.null(y) && !is.null(ysize)){
    aspratio <- nrow(mat) / ncol(mat) ## get aspect ratio of original image
    ymin <- y - ysize / 2
    ymax <- y + ysize / 2
    xmin <- x - (ysize / aspratio / 2 / xy_ratio)
    xmax <- x + (ysize / aspratio / 2 / xy_ratio)
  } else {
    ymin <- -Inf ## fill whole plot...
    ymax <- Inf
    xmin <- -Inf
    xmax <- Inf
  }
  rasterImage(mat, xmin, ymin, xmax, ymax)
}
