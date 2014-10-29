fun2postfix=function(f, p, pada) {
   # some calls must be translated to postfixes
   # e.g. nrow(a) -> (a).n_rows
   # In this example f is "nrow", p is ".n_rows"
   # pada is parsed data which has to be modified and returned.
   
   i=which(pada$token=="SYMBOL_FUNCTION_CALL" & pada$text==f)
   if (length(i) > 0) {
      # get opening paranthesis
      i=1L+sapply(i, function(ii) min(which(pada$id > pada$id[ii] &
         pada$token=="'('")))
      # get closing paranthesis
      i=sapply(i, function(ii) min(which(pada$line1 >= pada$line2[ii] &
         pada$col1 >= pada$col2[ii] & pada$token=="')'")))
      # append the code
      pada$text[i]=paste(pada$text[i], p, sep="")
   }
   return(pada)
}
symdim=function(st, env=parent.frame()) {
   # Get symbolic dimension of the result of an assignement statement in st
   # Statements are items in a list obtained as result of
   # e.g. quote(x=a+b)
   # return a character vector:
   # - of legth 2 for a matrix, e.g. c("nrow(mat)", "ncol(mat)")
   # - of length 1 for a vector, e.g. "length(vec)"
   # - a string "1" for a scalar
   # e.g. for a matrix a, symdim(a) should return c("nrow(a)", "ncol(a)")
   # If some error occurs, a NULL is returned
   # The argument 'env' is where objects from the statement are searched for
   #
   # Usage:
   # > lapply(parse(t="x=a%*%(b+c); y=a+b"), symdim)
   # or simply
   # > symdim(parse(t="x=a%*%(b+c)")[[1]])
   if (is.symbol(st)) {
      # we are ready to return the dim vector
      s1=as.character(st)
      obj=get(s1, env=env)
      cl=class(obj)
      len=length(obj)
      if (cl == "matrix") {
         return(c(sprintf("nrow(%s)", s1), sprintf("ncol(%s)", s1)))
      } else if (len == 1) {
         return("1")
      } else {
         return(sprintf("length(%s)", s1))
      }
   }
   s1=as.character(st[[1]])
   t1=typeof(st[[1]])
   m1=mode(st[[1]])
   if (any(s1 == c("=", "<-"))) {
      # return the dims of the RHS
      return(symdim(st[[3]]))
   } else if (any(s1 == c("+", "-", "*", "/"))) {
      # return the longest dim of all arguments
      d1=symdim(st[[2]])
      l1=length(d1)
      d2=symdim(st[[3]])
      l2=length(d2)
      if (l2 > l1) {
         return(d2)
      } else {
         return(d1)
      }
   } else if (any(s1 == c("t", "ginv"))) {
      d1=symdim(st[[2]])
      l1=length(d1)
      if (l1 == 1L) {
         # transpose of a vector
         return(c("1", d1))
      } else {
         # return the reversed dims of the first argument
         return(rev(d1))
      }
   } else if (any(s1 == c("%*%", "crossprod", "tcrossprod"))) {
      # return the dims by combyning two arguments
      d1=symdim(st[[2]])
      l1=length(d1)
      if (length(st) > 2L) {
        d2=symdim(st[[3]])
        l2=length(d2)
      } else {
         # crossprod or tcrossprod with only one argument
         d2=d1?
         l2=l1?
      }
      if (l1 == 1L && l2 == 1L) {
         # dot product of two vectors
         return("1")
      } else if (l1==2L && l2==1L) {
         # dot product mat-vec
         if (d1[1L]=="1") {
            # mat is just a t(vec)
            return("1")
         } else {
            return(d1[1L])
         }
      } else if (l1==1L && l2==2L) {
         # dot product vec-mat (equivalent to t(vec)%*%mat)
         return(c("1", d2[2L]))
      } else {
         # mat-mat dot product
         return(c(d1[1L], d2[2L]))
      }
      return(symdim(st[[2]]))
   } else {
      # we are not ready yet for all other cases
      return(NULL)
   }
}
# fastLm example
require(MASS)
y <- log(trees$Volume)
X <- cbind(1, log(trees$Girth))
fastLm_r=function(y, X) {
   df=nrow(X)-ncol(X)
   coef=qr.solve(X, y)
   res=y-X%*%coef
   s2=t(res)%*%res/df
   std_err=sqrt(s2*diag(ginv(t(X)%*%X)))
}
