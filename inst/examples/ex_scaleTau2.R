# 2017-01-25 sokol@insa-toulouse.fr
# seen on http://stats.stackexchange.com/questions/4017/translate-r-to-c-eventually-with-rcpp
# convert this R function to C++, execute and compare
scaleTau2<-function (x, c1 = 4.5, c2 = 3, consistency = "TRUE", mu.too = FALSE){
    n <- length(x)
    medx <- median(x)
    x1 <- abs(x - medx)
    sigma0 <- median(x1)
    mu <- if (c1 > 0) {
        x1 <- x1/(sigma0 * c1)
        w <- 1 - x1 * x1
        w <- ((abs(w) + w)/2)^2
        sum(x * w)/sum(w)
    } else medx
    x <- (x - mu)/sigma0
    rho <- x^2
    rho[rho > c2^2] <- c2^2
    if (!identical(consistency, "FALSE")) {
        Erho <- function(b) 2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1
        Es2 <- function(c2) Erho(c2*qnorm(3/4))
        nEs2 <-ifelse(consistency == "finiteSample",n-2,n)*Es2(c2)
    } else nEs2 <- n
    #c(if (mu.too) mu, sigma0 * sqrt(sum(rho)/nEs2))
    if (mu.too)
        c(mu, sigma0 * sqrt(sum(rho)/nEs2))
    else
        c(sigma0 * sqrt(sum(rho)/nEs2))
}
# create variables for prototyping function arguments 
set.seed(7); x=rnorm(3); c1 = 4.5; c2 = 3; b=c2*qnorm(3/4); consistency = "TRUE"; mu.too = TRUE; b=c2*qnorm(3/4); mu_too=mu.too
code=rex2arma(scaleTau2, exec=1); # create and compile c++ function scaleTau2_arma_()
print(scaleTau2_arma_(x, c1, c2, consistency, mu_too)) # execute it
# compare to original R function
print(scaleTau2(x, c1, c2, consistency, mu_too))
