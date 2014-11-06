# 2014-11-03
# Author: Serguei Sokol, sokol@insa-toulouse.fr
# Licence of RcppArmadillo applies here
# Copyright 2014, INRA, France

symdim=function(st, dim_tab=NULL, env=parent.frame()) {
   # Get symbolic dimension (as an R code) of a statement in st
   # Statements are items in a list obtained as result of
   # e.g. quote(a+b)
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
   # > symdim(parse(t="x=a%*%(b+c)")[[1L]])
#browser()
   if (is.symbol(st)) {
      # we are ready to return the dim vector
      s1=as.character(st)
      if (any(s1 == names(dim_tab))) {
         # already known
         return(dim_tab[[s1]])
      } else {
         obj=get(s1, mode="numeric", env=env)
         cl=class(obj)
         len=length(obj)
         if (cl == "matrix") {
            return(c(sprintf("nrow(%s)", s1), sprintf("ncol(%s)", s1)))
         } else if (len == 1) {
            # to be cast a scalar, it is a vector who must be of length 1, not a matrix
            return("1")
         } else {
            return(sprintf("length(%s)", s1))
         }
      }
   } else if (is.numeric(st)) {
      # a plain number
      return("1")
   }
   s1=as.character(st[[1L]])
   t1=typeof(st[[1L]])
   m1=mode(st[[1L]])
   lenst=length(st)
   if (any(s1 == c("=", "<-"))) {
      # return the dims of the RHS
      return(symdim(st[[3L]], dim_tab, env))
   }
   # dims af arguments
   dims=lapply(st[-1], symdim, dim_tab, env)
   if (any(s1 == c("+", "-", "*", "/")) && lenst == 3L) {
#browser()
      # return the longest dim of two arguments
      d1=dims[[1L]]
      l1=length(d1)
      d2=dims[[2L]]
      l2=length(d2)
      if (l2 > l1) {
         return(d2)
      } else if (l1 > l2) {
         return(d1)
      } else if (l1 == 1L) {
         # return the longest string
         return(c(d2, d1)[(nchar(d1) >= nchar(d2))+1L])
      } else {
         # two matrices => by default the first one
         return(d1)
      }
   } else if (any(s1 == c("+", "-")) && lenst == 2L) {
      # unary signs
      return(dims[[1L]])
   } else if (any(s1 == c("t", "ginv"))) {
      d1=dims[[1L]]
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
      return(c(dims[[1L]], dims[[2L]]))
   } else if (any(s1 == c("%*%", "crossprod", "tcrossprod"))) {
      # return the dims by combyning nrow(arg1) and ncol(arg2)
      d1=dims[[1L]]
      l1=length(d1)
      if (lenst > 2L) {
         d2=dims[[2L]]
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
      return(dims[[1L]])
   } else if (any(s1 == c("sum", "prod", "nrow", "ncol"))) {
      # reducing to scalar functions
      return("1")
   } else if (s1 == "c") {
      # sum of argument sizes
      return(sprintf("%s", paste(sapply(dims, function(d) if (length(d) > 1L) paste(d, collapse="*") else d), collapse="+")))
   } else if (any(s1 == c(":", "seq"))) {
#browser()
      # range size
      if (s1 == ":" && format(st[[2L]] == format(st[[3L]]))) {
         # identical begin and end => a scalar
         return("1")
      }
      # otherwise a vector
      return(sprintf("length(%s)", format(st)))
   } else if (any(s1 == c("solve", "qr.solve"))) {
      # dim of the solution of a linear system
      dims[[lenst-1L]]
   } else if (s1 == "diag") {
      # the result may be vector, may be matrix
      if (lenst == 2L && lenst[1L] == 2L) {
         # vector
         return(dims[[1L]][2L])
      } else {
         # matrix
         alast=format(st[[lenst]])
         return(c(alast, alast))
      }
   } else {
      # by default we suppose that the function operates on each term of its first argument
      return(dims[[1L]])
   }
}

st2arma=function(
st,
call2arma=c("qr.solve"="solve", "ginv"="pinv", "^"="pow"),
indent="",
known=NULL,
...) {
   # Translate an R statement st (from a parsed expression) to RcppArmadillo
   # code which is returned as a string.
   # For nested blocks in curved brackets "{...}", indent is incremented.
   # NB. newly declarated variables in these block are not visible outside them.
   # Parameters in ... are passed to symdim().
   if (is.symbol(st) || is.numeric(st)) {
      s1=as.character(st)
      if (s1 == "T") {
         return("true")
      } else if (s1 == "F") {
         return("false")
      }
      # just pass through
      return(s1)
   }
   if (isTRUE(st)) {
      return("true")
   }
   if (is.logical(st) && isTRUE(!st)) {
      return("false")
   }
   s1=as.character(st[[1L]])
   args=sapply(st[-1L], st2arma, call2arma, indent, known, ...)
   if (s1 == "=" || s1 == "<-") {
      # prepending by 'mat ', 'vec' or 'double '
      if (! args[1L] %in% known) {
         d=symdim(st[[3]], ...)
         len=length(d)
         if (len > 1L) {
            prep="mat "
         } else if (d == "1") {
            prep="double "
         } else {
            prep="vec "
         }
      }
      return(sprintf("%s%s%s=%s;\n", indent, prep, args[1L], args[2L]))
   }
   if (s1 == "if") {
      res=sprintf("%sif (%s) %s", indent, args[1L], args[2L])
      if (length(args) == 3L) {
         res=sprintf("%s else %s", res, args[3L])
      }
      return(sprintf("%s\n", res))
   }
   if (s1 == "{") {
      return(sprintf("{\n%s%s}", 
         paste(sprintf("%s   ", indent), args, collapse=""), indent))
   }
   dims=lapply(st[-1L], symdim, ...)
   lens=sapply(dims, length)
   lenst=length(st)
#browser()
   if (any(s1==c("*", "/", "+", "-", ">", "<", ">=", "<=", "==", "!=", "&&", "||")) && lenst == 3L) {
      # plain binary operations term by term
      # mat,vec operations '*'->'%'
      if (s1 == "*" && dims[[1L]][1L] != "1" && dims[[2L]][1L] != "1") {
         s1="%"
      }
      return(sprintf("%s%s%s", args[1L], s1, args[2L]))
   }
   if (s1 == "+" || s1 == "-" && lenst == 2L) {
      # unary operations
      return(sprintf("%s%s", s1, args[1L]))
   }
   if (s1 == "(") {
      # parenthesis operations
      return(sprintf("(%s)", args[1L]))
   }
   if (any(s1==c("%*%", crossprod, tcrossprod))) {
      # dot products
      a1=args[1L]
      d1=dims[[1L]]
      l1=lens[1L]
      if (lenst == 3L) {
         a2=args[2L]
         d2=dims[[2L]]
         l2=lens[2L]
      }
      if (lenst == 2L && s1=="crossprod") {
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
      } else if (lenst == 2L && s1=="tcrossprod") {
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
         if (l2==1L && d1[1L] != "1" && d2[1L] != "1") {
            # vec%*%vec
            return(sprintf("dot(%s, %s)", a1, a2))
         }
      }
      res=sprintf("%s*%s", a1, a2)
      if (d1[1L] == "1" && (l2 == 1L || d2[2L]=="1")) {
         # the result is a scalar
         res=sprintf("as_scalar(%s)", res)
      }
      return(res)
   } else if (s1 == "diag") {
      if (lenst == 2L) {
         # only one argument
         a1=args[1L]
         l1=lens[1L]
         if (dims[[1L]][1L] == "1") {
            # a scalar argument =>
            # create an identity matrix of size arg1
            return(sprintf("eye<mat>(%s, %s)", a1, a1))
         } else if (l1 == 1L) {
            # a vector argument
            # create a diagonal matrix with the vector on the main diagonal
            return(sprintf("diagmat(%s)", a1))
         } else if (l1 == 2L) {
            # a matrix argument
            # extract the main diagonal to a vector
            return(sprintf("diagvec(%s)", a1))
         }
      } else if (lenst == 3L) {
         # two arguments
         a1=args[1L]
         l1=lens[1L]
         a2=args[2L]
         l2=lens[2L]
         if (dims[[1L]][1L] == "1") {
            # first scalar argument =>
            # create a diagonal matrix of size arg2 and filled with arg1
            return(sprintf("diagmat(vec(as_scalar(%s)).fill(as_scalar(%s)))", a2, a1))
         } else if (lens[1L] == 1L) {
            # first vector argument
            # create a diagonal matrix with the vector on the main diagonal
            return(sprintf("diagmat(%s)", a1))
         } else if (lens[1L] == 2L) {
            # first matrix argument =>
            # error, it cannot be
            stop("If a matrix is the first argument to diag(), the second arguments is meaningless.")
         }
      }
   } else if (any(s1 == names(call2arma))) {
      # simply translate function names
      res=sprintf("%s(%s)", call2arma[s1], paste(args, collapse=", "))
      return(res)
   } else if (s1 == "nrow") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.n_rows", args[1L]))
      }
      return(sprintf("(%s).n_rows", args[1L]))
   } else if (s1 == "ncol") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.n_cols", args[1L]))
      }
      return(sprintf("(%s).n_cols", args[1L]))
   } else if (s1 == "t") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.t()", args[1L]))
      }
      return(sprintf("(%s).t()", args[1L]))
   } else if (s1 == "%o%") {
#browser()
      if (lens[[1L]] == 1 && dims[[1L]] == "1") {
         # first argument is a scalar
         if (!(lens[[2L]] == 1 && dims[[2L]] == "1")) {
            return(sprintf("%s*%s.t()", args[1L], args[2L]))
         } else {
            return(sprintf("%s*%s", args[1L], args[2L]))
         }
      } else if (lens[[2L]] == 1 && dims[[2L]] == "1") {
         # the second argument is a scalar
         return(sprintf("%s*%s", args[1L], args[2L]))
      }
      # all other cases
      return(sprintf("kron(vectorise(%s), vectorise(%s).t())", args[1L], args[2L]))
   } else if (any(s1 == c("c", "as.numeric", "as.double"))) {
      return(sprintf("vectorise(%s)", args[1L]))
   } else if (s1 == "list") {
#browser()
      nms=names(st[-1L])
      nms=ifelse(nchar(nms), sprintf('Named("%s")=', nms), nms)
      # convert vec in args to NumericVector
      argconv=sapply(seq_along(args), function(i) if (lens[i] > 1L || dims[[i]][1L]=="1") args[i] else sprintf("NumericVector(%s.begin(), %s.end())", args[i], args[i]))
      return(sprintf("List::create(%s)", paste(nms, argconv, sep="", collapse=", ")))
   } else {
      # by default return st as a function call
      return(sprintf("%s(%s)", s1, paste(args, collapse=", ")))
   }
}
