#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R.
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @param a the parameter for gibbsc
#' @param b the parameter for gibbsc
#' @param n the parameter for gibbsc
#' @return a random sample from gibbs \code{mat}
#' @importFrom stats rbinom rbeta
#' @import stats
#' @import boot
#' @import MASS
#' @import moments
#' @import RANN
#' @import energy
#' @import Ball
#' @import microbenchmark
#' @import bootstrap
#' @import graphics
#' @examples
#' \dontrun{
#'     gr <- gibbsR(100, 10, 2, 4, 16)
#'     print(gr)
#' }
#' @export
gibbsR <- function(N, thin, a, b, n){
  mat <- matrix(nrow = N, ncol = 2)
  x <- 0
  y <- 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, x + a, n - x + b)
    }
    mat[i, ] <- c(x, y)
  }
  return(mat)
}

#' @title Stochastic dominance tests for risk averters
#' @description Stochastic dominance tests for risk averters.
#' @param data1 the distribution you want to compare with data2
#' @param data2 the distribution you want to compare with data1
#' @param B the computing times for bootstrap
#' @return a list that contains tests result for FSD, SSD, TSD and some statistics \code{output}
#' @importFrom stats quantile
#' @examples
#' \dontrun{
#'     set.seed(666)
#'     data11 <- rnorm(200, 1, 1)
#'     data12 <- rnorm(100, 0, 1)
#'     test1 <- sdfra(data11, data12)
#'     print(test1$sd)
#'
#'     data21 <- rnorm(100, 0, 2)
#'     data22 <- rnorm(150, 0, 1)
#'     test2 <- sdfra(data21, data22)
#'     print(test2$sd)
#' }
#' @export
sdfra <- function(data1, data2, B=2000){

  stock1 <- subset(data1, data1 != "NA")
  stock2 <- subset(data2, data2 != "NA")
  allstock <- c(stock1,stock2)
  n1 <- length(stock1)
  n2 <- length(stock2)

  shrink <- 20/(max(allstock)-min(allstock))
  move <- -10*(max(allstock)+min(allstock))/(max(allstock)-min(allstock))
  stock1 <- stock1*shrink+move
  stock2 <- stock2*shrink+move
  allstock <- c(stock1,stock2)

  HA1 <- function(stock, x) {
    stockpos <- as.numeric(x-stock>=0)
    ha1 <- mean(stockpos)
    return(ha1)
  }

  HA2 <- function(stock, x) {
    stockpos <- (x-stock)*(x-stock>=0)
    ha2 <- mean(stockpos)
    return(ha2)
  }

  HA3 <- function(stock, x) {
    stockpos <- (x-stock)^2*(x-stock>=0)
    ha3 <- mean(stockpos)/2
    return(ha3)
  }

  VHA1 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- as.numeric(x-stock>=0)
    pos <- mean(stockpos)
    vha1 <- (pos - HA1(stock,x)^2)/ns
    return(vha1)
  }

  VHA2 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- (x-stock)^2*(x-stock>=0)
    pos <- mean(stockpos)
    vha2 <- (pos - HA2(stock,x)^2)/ns
    return(vha2)
  }

  VHA3 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- (x-stock)^4*(x-stock>=0)
    pos <- mean(stockpos)/4
    vha3 <- (pos - HA3(stock,x)^2)/ns
    return(vha3)
  }

  VA1 <- function(stocka, stockb, x){
    va1 <- VHA1(stocka, x) + VHA1(stockb, x)
    return(va1)
  }

  TA1 <- function(stocka, stockb, x){
    cc = HA1(stocka,x)-HA1(stockb,x)
    if(VA1(stocka, stockb, x)==0){
      ta1 <- 0}
    else{
      ta1 <- cc/(VA1(stocka, stockb, x))^0.5
    }
    return(ta1)
  }

  VA2 <- function(stocka, stockb, x){
    va2 <- VHA2(stocka, x) + VHA2(stockb, x)
    return(va2)
  }

  TA2 <- function(stocka, stockb, x){
    cc = HA2(stocka,x)-HA2(stockb,x)
    if(VA2(stocka, stockb, x)==0) {
      ta2 <- 0}
    else {
      ta2 <- cc/(VA2(stocka, stockb, x))^0.5
    }
    return(ta2)
  }

  VA3 <- function(stocka, stockb, x){
    va3 <- VHA3(stocka, x) + VHA3(stockb, x)
    return(va3)
  }

  TA3 <- function(stocka, stockb, x){
    cc = HA3(stocka,x)-HA3(stockb,x)
    if(VA3(stocka, stockb, x)==0) {
      ta3 <- 0}
    else {
      ta3 <- cc/(VA3(stocka, stockb, x))^0.5
    }
    return(ta3)
  }

  select_point <- seq(-9.9,10,0.1)
  TAab1 <- function(x) TA1(stock1, stock2, x)
  TAab2 <- function(x) TA2(stock1, stock2, x)
  TAab3 <- function(x) TA3(stock1, stock2, x)
  TA1re <- sapply(select_point, TAab1)
  TA2re <- sapply(select_point, TAab2)
  TA3re <- sapply(select_point, TAab3)
  l <- u <- lm <- um <- numeric(3)
  l[1] <- as.numeric(quantile(TA1re,0.05))
  l[2] <- as.numeric(quantile(TA2re,0.05))
  l[3] <- as.numeric(quantile(TA3re,0.05))
  u[1] <- as.numeric(quantile(TA1re,0.95))
  u[2] <- as.numeric(quantile(TA2re,0.95))
  u[3] <- as.numeric(quantile(TA3re,0.95))
  lm[1] <- min(TA1re)
  lm[2] <- min(TA2re)
  lm[3] <- min(TA3re)
  um[1] <- max(TA1re)
  um[2] <- max(TA2re)
  um[3] <- max(TA3re)

  A1 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TA1(sto1, sto2, x)
    res <- sapply(select_point, sab)
    a1 <- max(abs(res))
    return(a1)
  }

  A2 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TA2(sto1, sto2, x)
    res <- sapply(select_point, sab)
    a2 <- max(abs(res))
    return(a2)
  }

  A3 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TA1(sto1, sto2, x)
    res <- sapply(select_point, sab)
    a3 <- max(abs(res))
    return(a3)
  }

  methoda <- function(q){
    ff <- function(i){
      if(i>=1&i<=B){r <- A1(i)}
      else if(i>=(B+1)&i<=(2*B)){r <- A2(i)}
      else {r <- A3(i)}
      return(r)
    }
    res <- numeric(3*B)
    sam2 <- c(1:(3*B))
    res <- sapply(sam2,ff)
    return(res)
  }

  AA <- methoda(1)
  AA1 <- AA[1:B]
  AA2 <- AA[(B+1):(2*B)]
  AA3 <- AA[(2*B+1):(3*B)]
  limit <- numeric(3)
  limit[1] <- as.numeric(quantile(AA1,0.95))
  limit[2] <- as.numeric(quantile(AA2,0.95))
  limit[3] <- as.numeric(quantile(AA3,0.95))

  sd <- numeric(3)
  for(w in 1:3){
    if(u[w] > limit[w] & lm[w] > -limit[w]){
      sd[w] <- -1}
    else if(l[w] < -limit[w] & um[w] < limit[w]){
      sd[w] <- 1}
    else {
      sd[w] <- 0
    }
  }

  TA <- rbind(TA1re, TA2re, TA3re)
  bound <- rbind(um, u, l, lm, limit)
  output <- list(sd=sd,
                 TA=TA,
                 bound=bound)
  return(output)
}

#' @title Stochastic dominance tests for risk seekers
#' @description Stochastic dominance tests for risk seekers.
#' @param data1 the distribution you want to compare with data2
#' @param data2 the distribution you want to compare with data1
#' @param B the computing times for bootstrap
#' @return a list that contains tests result for FSD, SSD, TSD and some statistics \code{output}
#' @importFrom stats quantile
#' @examples
#' \dontrun{
#'     set.seed(666)
#'     data31 <- rnorm(200, 2, 1)
#'     data32 <- rnorm(100, 0, 1)
#'     test3 <- sdfrs(data31, data32)
#'     print(test3$sd)
#'
#'     data41 <- rnorm(100, 0, 2)
#'     data42 <- rnorm(150, 0, 1)
#'     test4 <- sdfrs(data41, data42)
#'     print(test4$sd)
#' }
#' @export
sdfrs <- function(data1, data2, B=2000){

  stock1 <- subset(data1, data1 != "NA")
  stock2 <- subset(data2, data2 != "NA")
  allstock <- c(stock1,stock2)
  n1 <- length(stock1)
  n2 <- length(stock2)
  shrink <- 20/(max(allstock)-min(allstock))
  move <- -10*(max(allstock)+min(allstock))/(max(allstock)-min(allstock))

  stock1 <- stock1*shrink+move
  stock2 <- stock2*shrink+move
  allstock <- c(stock1,stock2)

  HD1 <- function(stock, x) {
    stockpos <- as.numeric(x-stock<=0)
    hd1 <- mean(stockpos)
    return(hd1)
  }

  VHD1 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- as.numeric(x-stock<=0)
    pos <- mean(stockpos)
    vhd1 <- (pos - HD1(stock,x)^2)/ns
    return(vhd1)
  }

  VD1 <- function(stocka, stockb, x){
    vd1 <- VHD1(stocka, x) + VHD1(stockb, x)
    return(vd1)
  }

  TD1 <- function(stocka, stockb, x){
    cc <- HD1(stocka,x)-HD1(stockb,x)
    if(VD1(stocka, stockb, x)==0){
      td1 <- 0}
    else{
      td1 <- cc/(VD1(stocka, stockb, x))^0.5
    }
    return(td1)
  }

  HD2 <- function(stock, x) {
    stockpos <- (x-stock)*(x-stock<=0)
    hd2 <- mean(stockpos)
    return(hd2)
  }

  VHD2 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- (x-stock)^2*(x-stock<=0)
    pos <- mean(stockpos)
    vhd2 <- (pos - HD2(stock,x)^2)/ns
    return(vhd2)
  }

  VD2 <- function(stocka, stockb, x){
    vd2 <- VHD2(stocka, x) + VHD2(stockb, x)
    return(vd2)
  }

  TD2 <- function(stocka, stockb, x){
    cc <- HD2(stocka,x)-HD2(stockb,x)
    if(VD2(stocka, stockb, x)==0) {
      td2 <- 0}
    else {
      td2 <- cc/(VD2(stocka, stockb, x))^0.5
    }
    return(td2)
  }

  HD3 <- function(stock, x) {
    stockpos <- (x-stock)^2*(x-stock<=0)
    hd3 <- mean(stockpos)/2
    return(hd3)
  }

  VHD3 <- function(stock, x) {
    ns <- length(stock)
    stockpos <- (x-stock)^4*(x-stock<=0)
    pos <- mean(stockpos)/4
    vhd3 <- (pos - HD3(stock,x)^2)/ns
    return(vhd3)
  }

  VD3 <- function(stocka, stockb, x){
    vd3 <- VHD3(stocka, x) + VHD3(stockb, x)
    return(vd3)
  }

  TD3 <- function(stocka, stockb, x){
    cc <- HD3(stocka,x)-HD3(stockb,x)
    if(VD3(stocka, stockb, x)==0) {
      td3 <- 0}
    else {
      td3 <- cc/(VD3(stocka, stockb, x))^0.5
    }
    return(td3)
  }

  select_point <- seq(-10,9.9,0.1)
  TDab1 <- function(x) TD1(stock1, stock2, x)
  TDab2 <- function(x) TD2(stock1, stock2, x)
  TDab3 <- function(x) TD3(stock1, stock2, x)
  TD1re <- sapply(select_point, TDab1)
  TD2re <- sapply(select_point, TDab2)
  TD3re <- sapply(select_point, TDab3)
  l <- u <- lm <- um <- numeric(3)
  l[1] <- as.numeric(quantile(TD1re,0.05))
  l[2] <- as.numeric(quantile(TD2re,0.05))
  l[3] <- as.numeric(quantile(TD3re,0.05))
  u[1] <- as.numeric(quantile(TD1re,0.95))
  u[2] <- as.numeric(quantile(TD2re,0.95))
  u[3] <- as.numeric(quantile(TD3re,0.95))
  lm[1] <- min(TD1re)
  lm[2] <- min(TD2re)
  lm[3] <- min(TD3re)
  um[1] <- max(TD1re)
  um[2] <- max(TD2re)
  um[3] <- max(TD3re)

  D1 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TD1(sto1, sto2, x)
    res <- sapply(select_point, sab)
    d1 <- max(abs(res))
    return(d1)
  }

  D2 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TD2(sto1, sto2, x)
    res <- sapply(select_point, sab)
    d2 <- max(abs(res))
    return(d2)
  }

  D3 <- function(i){
    sto1 <- sample(allstock, n1, replace = TRUE)
    sto2 <- sample(allstock, n2, replace = TRUE)
    sab <- function(x) TD1(sto1, sto2, x)
    res <- sapply(select_point, sab)
    d3 <- max(abs(res))
    return(d3)
  }

  methodd <- function(q){
    ff <- function(i){
      if(i>=1&i<=B){r <- D1(i)}
      else if(i>=(B+1)&i<=(2*B)){r <- D2(i)}
      else {r <- D3(i)}
      return(r)
    }
    res <- numeric(3*B)
    sam2 <- c(1:(3*B))
    res <- sapply(sam2,ff)
    return(res)
  }

  DD <- methodd(1)
  DD1 <- DD[1:B]
  DD2 <- DD[(B+1):(2*B)]
  DD3 <- DD[(2*B+1):(3*B)]
  limit <- numeric(3)
  limit[1] <- as.numeric(quantile(DD1,0.95))
  limit[2] <- as.numeric(quantile(DD2,0.95))
  limit[3] <- as.numeric(quantile(DD3,0.95))

  sd <- numeric(3)
  for(w in 1:3){
    if(u[w] > limit[w] & lm[w] > -limit[w]){
      sd[w] <- (-1)^(w+1)}
    else if(l[w] < -limit[w] & um[w] < limit[w]){
      sd[w] <- (-1)^(w)}
    else {
      sd[w] <- 0
    }
  }

  TD <- rbind(TD1re, TD2re, TD3re)
  bound <- rbind(um, u, l, lm, limit)
  output <- list(sd=sd,
                 TD=TD,
                 bound=bound)
  return(output)
}
