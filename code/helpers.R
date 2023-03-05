## Source .Rmd files 

source_rmd = function(file, ...) {
  tmp_file = tempfile(fileext=".R")
  on.exit(unlink(tmp_file), add = TRUE)
  knitr::purl(file, output=tmp_file, quiet = T)
  source(file = tmp_file, ...)
}


## function to recode the IPUMS education code to years of educaiton 
recode_education <- function(df) {
  df <- df  %>%
    mutate(educ_yrs = case_when(
      EDUCD == 2 ~ 0,
      EDUCD == 14 ~ 1,
      EDUCD == 15~  2,
      EDUCD == 16 ~ 3,
      EDUCD == 17 ~ 4,
      EDUCD == 22 ~ 5,
      EDUCD == 23 ~ 6,
      EDUCD == 25 ~ 7,
      EDUCD == 26 ~ 8,
      EDUCD == 30 ~ 9,
      EDUCD == 40 ~ 10,
      EDUCD == 50 ~ 11,
      EDUCD == 60 ~ 12,
      EDUCD == 70 ~ 13,
      EDUCD == 80 ~ 14,
      EDUCD == 90 ~ 15,
      EDUCD == 100 ~ 16,
      EDUCD == 111 ~ 17,
      EDUCD == 112 ~ 17,
      EDUCD == 113 ~ 17
    ))
  return(df)
}



## verify our gomp truncated MLE estimators with simulated data using
## an input format of a vector of deaths with age labels

hgompertz.M <- function(x, b, M)
{
    b * exp(b * (x-M))
}

pgompertz.M.new <- function(x, b, M)
{
    Hx <- exp(-b*M) * ( exp(b*x) - 1 )
    lx <- exp(-Hx)
    1-lx
}

pgompertz.M.old <- function(x, b, M)
{
    Hx <- exp(b*x - b*M)
    lx <- exp(-Hx) ## used to be exp(-Hx)
    1-lx
}

p.old <- pgompertz.M.old(0:110, .1, 71)
p.new <- pgompertz.M.new(0:110, .1, 71)

pgompertz.M <- pgompertz.M.new

dgompertz.M <- function(x, b, M)
{
    hx <- hgompertz.M(x,b,M)
    lx <- 1-pgompertz.M(x,b,M)
    dx <- hx*lx
    dx
}

bM2a <- function(b, M)
{
    a = b* exp(-b * M)
    a
}

ab2M <- function(a, b)
{
    M = -log(a/b)/b
    M
}



inverse.distribution.function.gomp <- function(u, b, a)
{
    -(  log(a) - log(b - log(1-u)*b) ) / b
}
igomp <- inverse.distribution.function.gomp

idfgomp <- function(u, b, a)
{
    x = (1/b) * log(-log(1-u)*b/a + 1)
    x
}
idfgomp(u = .5, b = .1, a = .001)

rgompertz2 <- function(N, b, a)
{
    u <- runif(N)
    x <- idfgomp(u, b, a)
    x
}

rgompertz.M <- function(N, b, M)
{
    a = bM2a(b,M)
    x = rgompertz2(N, b, a)
    x
}

d.gomp.negLL <- function(par, x, x.left, x.right)
{
    ## discrete gompertz negative log likelihood
    ## for ages that are measured in discrete single years
    b = exp(par[1]) ## we search on log.b (so b > 0)
    M = exp(par[2]) ## we also search on log.M
    num.vec <- pgompertz.M(x+1, b, M) - pgompertz.M(x, b, M) ## discrete
    denom.vec <- pgompertz.M(x.right, b, M) - pgompertz.M(x.left, b, M)
    LL <- sum(log(num.vec) - log(denom.vec))
    negLL = -LL
    return(negLL)

}

d.counts.gomp.negLL <- function(par, Dx, x.left, x.right)
{
    ## discrte gomp neg log likelihood for
    ## counts of deaths by single years of age (Dx)
    if(is.null(names(Dx)))
        print("names(Dx) is NULL; please define with ages")
    b = exp(par[1])
    M = exp(par[2])
    age.x = as.numeric(names(Dx))
    num.vec.unwt <- pgompertz.M(age.x+1, b, M) - pgompertz.M(age.x, b, M)
    denom.vec.unwt <- pgompertz.M(x.right, b, M) - pgompertz.M(x.left, b, M)
    LL <- sum( Dx*( log(num.vec.unwt) - log(denom.vec.unwt)) )
    negLL = -LL
    return(negLL)
}
## d.counts.gomp.negLL(par = log(c(1/9, 80)), Dx, x.left, x.right)

counts.trunc.gomp.est <- function(Dx, x.left, x.right, b.start, M.start)
{
    par.start = c(log(b.start), log(M.start))
    my.control = list(trace = 0,
                      parscale = c(par.scale = par.start))

    fit <- optim(par = par.start,
                 fn = d.counts.gomp.negLL,
                 hessian = TRUE,
                 Dx = Dx,
                 x.left = x.left,
                 x.right = x.right,
                 control = my.control)
    return(fit)

}

## my.N = 100000
## x.full <- rgompertz.M(my.N, b = 1/10, M = 70)
## x.left = 65
## x.right = 75
## x <- x.full[x.full > x.left & x.full < x.right]
## Dx <- table(floor(x))
## ## debug(d.counts.gomp.negLL)
## out <- counts.trunc.gomp.est(Dx, x.left, x.right, b.start = 1/9, M = 80)
## print(out$par)


trunc.gomp.est <- function(x, x.left, x.right, b.start, M.start,
                           type = "discrete")
{
    par.start = c(log(b.start), log(M.start))
    my.control = list(trace = 0,
                      parscale = c(par.scale = par.start))

    ## print(par.start)
    if (type == "discrete")
        fit <- optim(par = par.start,
                     fn = d.gomp.negLL,
                     hessian = TRUE,
                     x = x,
                     x.left = x.left,
                     x.right = x.right,
                     control = my.control)

    if (type == "continuous")
        fit <- optim(par = par.start,
                     fn = c.gomp.negLL,
                     hessian = TRUE,
                     x = x,
                     x.left = x.left,
                     x.right = x.right,
                     control = my.control)
    return(fit)
}

## check to see if it's the same as with floor(x)

check = FALSE
if (check) {
    x <- floor(x)
    out.check <- trunc.gomp.est(x, x.left, x.right,
                                b.start = 1/9, M = 80, type = "discrete")
    out$value - out.check$value
    out$par - out.check$par
}

get.se <- function(est.vec, hess, log.vec = rep(TRUE, length(est.vec)))
{
    fisher_info = solve(hess)
    est.sd <- sqrt(diag(fisher_info))
    ## if estimates are not log scale then we have the SE
    ## if they are log scale than we can estimate SE
    ## as exp(beta)*SD
    est.se <- exp(est.vec)^(log.vec) * est.sd
    return(est.se)
}
## get.se(out$par, out$hess, log.vec = c(TRUE, TRUE))

get.se.from.out <- function(out) {
    get.se(out$par, out$hess, log.vec = rep(TRUE, 2))
}

## get.se.from.out(out)

summarize.out <- function(out, log.vec = rep(TRUE, length(out$par)),
                          var.names = c("b", "M"))
{
    ## put estimates, se, and conf interval in nice format
    est <- out$par
    est <- exp(out$par)[log.vec == TRUE]
    se <- get.se.from.out(out)
    ci <- paste("(", round(est - 2*se, 2), "-", round(est + 2* se, 2), ")")

    summ <- cbind("est" = round(est, 4), "se" = round(se, 4), ci)
    rownames(summ) <- var.names
    cat("   ", "est", "   se", "    CI", "\n")
    cat("b ", "", summ[1,], "\n")
    cat("M ", summ[2,], "\n")
}

if(check) {
    my.N = 100000
    x.full <- rgompertz.M(my.N, b = 1/10, M = 70)
    x.left = 65
    x.right = 75
    x <- x.full[x.full > x.left & x.full < x.right]
    Dx <- table(floor(x))
    out <- counts.trunc.gomp.est(Dx, x.left, x.right, b.start = 1/9, M = 80)
    print(out$par)
    summarize.out(out)
    ##     est    se     CI
    ## b   0.1048 0.0222 ( 0.06 - 0.15 )
    ## M  69.5572 0.6433 ( 68.27 - 70.84 )

    ## seems to work. would also like fitted counts

    plot(Dx)

    ## we can use gomp parameters to estimate truncated distribution in
    ## proportions and then we can apply these to our N = sum(Dx)
    b.hat = exp(out$par)[1]
    M.hat = exp(out$par)[2]
}
get.Dx.hat <- function(Dx, b.hat, M.hat)
{
    xx = 0:200 ## high ages so it adds up to 1
    dx.hat = pgompertz.M(xx + 1, b = b.hat, M = M.hat) -
        pgompertz.M(xx, b = b.hat, M = M.hat)
    if(sum(dx.hat) != 1)
        print(paste("warning: sum(dx.hat) = ", sum(dx.hat)))
    dx.hat.trunc <- prop.table(dx.hat[xx %in% names(Dx)])
    N <- sum(Dx)
    Dx.hat.trunc <- dx.hat.trunc * N
    names(Dx.hat.trunc) <- names(Dx)
    return(Dx.hat.trunc)
}


if (check) {
    Dx.hat.trunc <- get.Dx.hat(Dx, b.hat, M.hat)
    lines(65:74, Dx.hat.trunc, col = "red")


    ## compactly
    my.N = 10000
    x.full <- rgompertz.M(my.N, b = 1/10, M = 80)
    x.left = 65
    x.right = 75
    x <- x.full[x.full > x.left & x.full < x.right]
    Dx <- table(floor(x))
    out <- counts.trunc.gomp.est(Dx, x.left, x.right, b.start = 1/9, M = 80)
    b.hat = exp(out$par)[1]
    M.hat = exp(out$par)[2]
    Dx.hat.trunc <- get.Dx.hat(Dx, b.hat, M.hat)
    plot(Dx)
    lines(65:74, Dx.hat.trunc, col = "red")
    summarize.out(out)
}

