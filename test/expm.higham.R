expm.higham=function (A, balancing = TRUE) 
{
    d <- dim(A)
    if (length(d) != 2 || d[1] != d[2]) 
        stop("'A' must be a square matrix")
    n <- d[1]
    if (n <= 1) 
        return(exp(A))
    if (balancing) {
        baP <- balance(A, "P")
        baS <- balance(baP$z, "S")
        A <- baS$z
    }
    nA <- Matrix::norm(A, "1")
#print(nA);
    I <- diag(n)
    if (nA <= 2.1) {
        t <- c(0.015, 0.25, 0.95, 2.1)
        l <- which.max(nA <= t)
        C <- rbind(c(120, 60, 12, 1, 0, 0, 0, 0, 0, 0), c(30240, 
            15120, 3360, 420, 30, 1, 0, 0, 0, 0), c(17297280, 
            8648640, 1995840, 277200, 25200, 1512, 56, 1, 0, 
            0), c(17643225600, 8821612800, 2075673600, 302702400, 
            30270240, 2162160, 110880, 3960, 90, 1))
#print(C);
        A2 <- A %*% A
        P <- I
        U <- C[l, 2] * I
        V <- C[l, 1] * I
        for (k in 1:l) {
            P <- P %*% A2
            U <- U + C[l, (2 * k) + 2] * P
            V <- V + C[l, (2 * k) + 1] * P
        }
        U <- A %*% U
        X <- solve(V - U, V + U)
    }
    else {
        s <- log2(nA/5.4)
        B <- A
        if (s > 0) {
            s <- ceiling(s)
            B <- B/(2^s)
        }
#print(B)
        c. <- c(64764752532480000, 32382376266240000, 7771770303897600, 
            1187353796428800, 129060195264000, 10559470521600, 
            670442572800, 33522128640, 1323241920, 40840800, 
            960960, 16380, 182, 1)
        B2 <- B %*% B
        B4 <- B2 %*% B2
        B6 <- B2 %*% B4
        U <- B %*% (B6 %*% (c.[14] * B6 + c.[12] * B4 + c.[10] * 
            B2) + c.[8] * B6 + c.[6] * B4 + c.[4] * B2 + c.[2] * 
            I)
        V <- B6 %*% (c.[13] * B6 + c.[11] * B4 + c.[9] * B2) + 
            c.[7] * B6 + c.[5] * B4 + c.[3] * B2 + c.[1] * I
        X <- solve(V - U, V + U)
#print(U)
#print(V)
#print(X)
#print(s)
        if (s > 0) 
            for (t in 1:s) X <- X %*% X
#print(X)
    }
    if (balancing) {
        dd <- baS$scale
        X <- X * (rep(dd, n) * rep(1/dd, each = n))
        pp <- as.integer(baP$scale)
        if (baP$i1 > 1) {
            for (i in (baP$i1 - 1):1) {
                tt <- X[, i]
                X[, i] <- X[, pp[i]]
                X[, pp[i]] <- tt
                tt <- X[i, ]
                X[i, ] <- X[pp[i], ]
                X[pp[i], ] <- tt
            }
        }
        if (baP$i2 < n) {
            for (i in (baP$i2 + 1):n) {
                tt <- X[, i]
                X[, i] <- X[, pp[i]]
                X[, pp[i]] <- tt
                tt <- X[i, ]
                X[i, ] <- X[pp[i], ]
                X[pp[i], ] <- tt
            }
        }
    }
    X
}
