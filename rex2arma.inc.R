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
symdim=function(st, dim_tab=NULL, env=parent.frame()) {
   # Get symbolic dimension of the result of an assignement statement in st
   # Statements are items in a list obtained as result of
   # e.g. quote(x=a+b)
   # return a character vector:
   # - of legth 2 for a matrix, e.g. c("nrow(mat)", "ncol(mat)")
   # - of length 1 for a vector, e.g. "length(vec)"
   # - a string "1" for a scalar
   # e.g. for a matrix a, symdim(a) should return c("nrow(a)", "ncol(a)")
   # If some error occurs, a NULL is returned
   # The argument 'dim_tab' is a named list or env where object
   # symbolic dimensions are retrieved from. 
   # The argument 'env' is where objects from the statement are searched for
   # if dim_tab is NULL or an object is not there.
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
   } else if (s1 == "%o%") {
      # tensor product of two (normally) vectors
      return(c(symdim(st[[2]], symdim(st[[3]]))))
   } else if (any(s1 == c("%*%", "crossprod", "tcrossprod"))) {
      # return the dims by combyning nrow(arg1) and ncol(arg2)
      d1=symdim(st[[2]])
      l1=length(d1)
      if (length(st) > 2L) {
         d2=symdim(st[[3]])
         l2=length(d2)
      } else {
         # crossprod or tcrossprod with only one argument
         d2=d1
         l2=l1
      }
      if (s1=="crossprod") {
         d1=rev(d1) # first argument is transposed
      } else if (s1=="tcrossprod") {
         # second argument is transposed
         if (l2 == 1L) {
            d2=c("1", d2)
            l2=2L
         } else {
            d2=rev(d2)
         }
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
         # dot product vec-mat (equivalent to t(vec)%*%mat except when nrow(mat)==1)
         if (d2[1L]=="1") {
            # mat is just a t(vec)
            return(c(d1, d2[2L]))
         } else {
            return(c("1", d2[2L]))
         }
      } else {
         # mat-mat dot product
         return(c(d1[1L], d2[2L]))
      }
      return(symdim(st[[2]]))
   } else if (any(s1 == c("sum", "prod"))) {
      # reducing to scalar functions
      return("1")
   } else {
      # we are not ready yet for all other cases
      return(NULL)
   }
}
st2arma=function(st, call2arma=c("qr.solve"="solve", "ginv"="pinv", "diag"="diagvec"), ...) {
   # Translate an R statement st (from a parsed expression) to RcppArmadillo
   # code which is returned as a string
   # Optional params '...' are passed to symdim.
   if (is.symbol(st)) {
      # just pass through
      return(as.character(st))
   }
   s1=as.character(st[[1]])
   if (any(s1==c("*", "/", "+", "-")) && length(st)==3L) {
      # plain binary operations term by term
      # mat,vec operations '*'->'%'
      if (s1 == "*" && symdim(st[[2L]], ...) != "1" && symdim(st[[3L]], ...) != "1") {
         s1="%"
      }
      return(sprintf("%s%s%s", st2arma(st[[2L]]), s1, st2arma(st[[3L]])))
   } else if (any(s1==c("%*%", crossprod, tcrossprod))) {
      # dot products
      a1=st2arma(st[[2L]])
      d1=symdim(st[[2L]], ...)
      l1=length(d1)
      if (length(st) == 3L) {
         a2=st2arma(st[[3L]])
         d2=symdim(st[[3L]], ...)
         l2=length(d2)
      }
      if (length(st) == 2L && s1=="crossprod") {
         # just one arg for crossprod
         a2=a1
         a1=sprintf("(%s).t()", a1)
         if (l1==1L) {
            d2="1"
            l2=1L
            d1=c("1", d1)
            l1=2L
         } else {
            d2=d1
            l2=l1
            d1=rev(d1)
         }
      } else if (length(st) == 2L && s1=="tcrossprod") {
         # just one arg for tcrossprod
         a2=sprintf("(%s).t()", a1)
         if (l1==1L) {
            d2=c("1", d1)
            l2=2L
         } else {
            d2=rev(d1)
            l2=l1
         }
      } else if (s1=="%*%" && l1==1L) {
         # vec%*%smth, so decide vec.t() or not
         if (l2==1L) {
            # vec%*%vec
            a1=sprintf("t(%s)", a1)
         }
      }
      res=sprintf("%s*%s", a1, a2)
      if (d1[1L] == "1" && (l2 == 1L || d2[2L]=="1")) {
         # the result is a scalar
         res=sprintf("as_scalar(%s)", res)
      }
      return(res)
   } else if (any(s1 == names(call2arma))) {
      # translate function names
      args=sapply(st[2:length(st)], st2arma, call2arma, ...)
      res=sprintf("%s(%s)", call2arma[s1], paste(args, collapse=", "))
      return()
   } else {
      # by default return st "as is"
      return(format(st))
   }
}
