#' Get symbolic dimension of a statement (as an R code)
#' 
#' @param st A statements, i.e. an item in a list obtained as result of
#' e.g. quote(a+b) or expression(a+b)
#' @param env An environement where a statement can be executed to get
#'  its type and structure.
#' @param dim_tab A named list caching symbolic dimensions for R symbols.
#'
#' @return a character vector:\enumerate{
#' \item of legth 2 for a matrix \code{mat}, it returns
#'  \code{c("nrow(mat)", "ncol(mat)")}
#' \item of length 1 for a vector, e.g. "length(vec)"
#' \item a string "1" for a scalar}
#' @details
#' If some error occurs, a NA is returned
#' The objects from the statement \code{st} are searched for in \code{env}
#' if dim_tab is NULL or an object is not there.
#
#' @examples
#' a=b=c=1:2
#' lapply(parse(t="x <- a%*%(b+c); y=a+b"), symdim)
#' # or
#' symdim(expression(x <- a%*%(b+c))[[1L]])

symdim=function(st, env=parent.frame(), dim_tab=NULL) {
#browser()

   if (!is.call(st) || as.character(st[[1L]]) == "$") {
      # we are ready to return the dim vector
      s1=gsub("\\.", "_", if (is.call(st)) format(st) else as.character(st))
      if (s1 %in% names(dim_tab)) {
         # already known
         return(dim_tab[[s1]])
      } else if (s1 == "") {
         return("")
      } else {
         obj=eval(st, env=env)
         len=length(obj)
         if (len == 1) {
             # to be cast a scalar, it is a vector who must be of length 1, not a matrix
             return("1")
         }
         if (is.matrix(obj) || is(obj, "Matrix")) {
            return(c(sprintf("nrow(%s)", s1), sprintf("ncol(%s)", s1)))
         }
         if (is.vector(obj)) {
             return(sprintf("length(%s)", s1))
         }
         # for not scalar, vector, matrix
         return(NA)
      }
   }
   if (is.numeric(st) || is.logical(st) || is.character(st)) {
      # a plain number
      return("1")
   }
   s1=as.character(st[[1L]])
   t1=typeof(st[[1L]])
   m1=mode(st[[1L]])
   lenst=length(st)
   if (s1 %in% c("=", "<-")) {
      # return the dims of the RHS
      return(symdim(st[[3L]], env, dim_tab))
   }
   # dims of arguments
   dims=lapply(st[-1], symdim, env, dim_tab)
   lens=sapply(dims, length)
   if (any(s1 == c("+", "-", "*", "/", "^", ">", "<", "<=", ">=", "!=", "==")) && lenst == 3L) {
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
   }
   if (any(s1 == c("+", "-")) && lenst == 2L) {
      # unary signs
      return(dims[[1L]])
   }
   if (any(s1 == c("t", "ginv"))) {
      d1=dims[[1L]]
      l1=length(d1)
      if (l1 == 1L) {
         # transpose of a vector
         return(c("1", d1))
      } else {
         # return the reversed dims of the first argument
         return(rev(d1))
      }
   }
   if (s1 == "%o%") {
      # tensor product of two (normally) vectors
      return(c(dims[[1L]], dims[[2L]]))
   }
   if (any(s1 == c("%*%", "crossprod", "tcrossprod"))) {
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
            return(c(d1[1L], "1"))
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
   }
   if (s1 %in% c("sum", "prod", "nrow", "ncol", "length", "mean", "min", "max", "sd", "norm")) {
      # reducing to scalar functions
      return("1")
   }
   if (s1 == "var") {
      if (lens[1L] == 2L) {
         return(c(dims[[1L]][2L], dims[[1L]][2L]))
      } else {
         return("1")
      }
   }
   if (s1 == "dim") {
      return(lens[1L])
   }
   if (s1 == "c") {
      # sum of argument sizes
      return(sprintf("%s", paste(sapply(dims, function(d) if (length(d) > 1L) paste(d, collapse="*") else d), collapse="+")))
   }
   if (s1 %in% c("solve", "qr.solve")) {
      # dim of the solution of a linear system
      return(dims[[lenst-1L]])
   }
   if (s1 == "diag" || s1 == "Diagonal") {
      # the result may be vector, may be matrix
      if (lenst == 2L && lenst[1L] == 2L) {
         # vector
         return(dims[[1L]][2L])
      } else {
         # matrix
         alast=format(st[[lenst]])
         return(c(alast, alast))
      }
   }
   
   # data/code control
   if (s1 %in% c(":", "seq")) {
#browser()
      # range size
      if (s1 == ":" && st[[2L]] == st[[3L]]) {
         # identical begin and end => a scalar
         return("1")
      }
      # otherwise a vector
      return(sprintf("length(%s)", format(st)))
   }
   if (s1 == "rep") {
      # repeated vector
      args=sapply(st[-1L], function(it) if (is.call(it)) format(it) else as.character(it))
      # match argument names times and each
      times=which(!is.na(pmatch(names(st[-1L]), "times")))
      if (length(times) > 1L) {
         stop(sprintf("Argument 'times' was found more than once in '%s'",
            format(st)))
      }
      each=which(!is.na(pmatch(names(st[-1L]), "each")))
      if (length(each) > 1L) {
         stop(sprintf("Argument 'each' was found more than once in '%s'",
            format(st)))
      }
      x=which(!is.na(pmatch(names(st[-1L]), "x")))
      if (length(x) > 1L) {
         stop(sprintf("Argument 'x' was found more than once in '%s'",
            format(st)))
      }
      if (length(x) == 0L) {
         # x will be the first unnamed argument
         x=which.max(nchar(args)==0L)
      }
      if (length(times) == 0L) {
         # times will be the first unnamed argument after x removing
         times=which(nchar(args[-x])==0L)
         if (length(times)) {
            times=times[1L]
            times=times+(times>=x)
         }
      }
      # 'each' and 'times' may be empty separatly, not both
      if (length(times) && lens[times] == 1L && dims[[times]] == "1") {
         # the whole vector is repeated 'times' times
         res=sprintf("%s*(%s)", paste("(", dims[[x]], ")", sep="", collapse="*"),
            args[times])
      } else if (length(times) && lens[times] == 1L) {
         # 'times' is a vector
         res=sprintf("sum(%s)", args[times])
      } else if (length(times) == 0L) {
         res=sprintf("%s", paste("(", dims[[x]], ")", sep="", collapse="*"))
      }
      if (length(each)) {
         # repeat eache element
         res=sprintf("(%s)*(%s)", res, args[each])
      }
      return(res)
   }
   
   if (s1 == "[") {
#browser()
      # indexing operator
      if (max(lens[-1L]) > 1L) {
         # index is matrix => result of indexing is a vector
         stopifnot(length(args) != 2L) # only one matrix in index
         return(sprintf("nrow(%s)", args[2L]))
      }
      # vector index(es) of vector or matrix?
      if (lens[1L] > 1L) {
         # indexes of a matrix
         if (dims[[2L]] == "1" && dims[[3L]] == "1") {
            return("1")
         }
         return(c(dims[[2L]], dims[[3L]]))
      }
      # getting subvector
      return(dims[[2L]])
   }
   if (s1 == "(") {
      # pass through the dims of the first argument
      return(dims[[1L]])
   }
   
   # by default
   warning(sprintf("Couldn't retrive dimension for '%s'", format(st)))
   return(NA)
}

#' Translate an elementary R statement to RcppArmadillo
#'
#' @param st A statement to transalte
#' @param call2arma A named character vector giving Armadillo equivalence
#'  for some R functions. It is indexed by R function names.
#' @param indent A character string used as indentation. It is incremented
#'  by "   " (3 spaces) inside a \code{{...}} block.
#' @param iftern A logical indicating whether \code{if/else} must be
#'  considered as a ternary operator (TRUE) or a classical \code{if/else}
#'  operator (FALSE).
#' @param env An environment where statement may be executed.
#' @param ... Parameters passed through to \code{symdim()} calls.
#' @return A string with RcppArmadillo code.
#'
#' @details
#' Parameter \code{env} is also passed to \code{symdim()} with \code{...}.

st2arma=function(
st,
call2arma=c("qr.solve"="solve", "ginv"="pinv", "^"="pow", "stop"="stop",
   "ceiling"="ceil", "which.min"="which_min", "which.max"="which_max",
   "Re"="real", "Im"="imag"),
indent="",
iftern=FALSE,
env=parent.frame(),
...) {
   if (!is.call(st)) {
      s1=as.character(st)
      if (is.character(st)) {
         return(sprintf('"%s"', gsub('(["\\])', '\\\\\\1', s1)))
      }
      if (s1 == "T" || s1 == "TRUE") {
         return("true")
      } else if (s1 == "F" || s1 == "FALSE") {
         return("false")
      }
      # replace "." by "_" in symbols
      if (is.symbol(st)) {
         s1=gsub("\\.", "_", s1)
      }
      return(s1)
   }
   
#browser()
   if (is.call(st[[1L]]) && as.character(st[[1L]][[1L]]) == "::") {
      s1=sprintf("%s_%s_r_", as.character(st[[1L]][[2L]]), as.character(st[[1L]][[3L]]))
   } else {
#browser()
      s1=as.character(st[[1L]])
      # check if it is a function absent in base environment
      if (!is.function(mget(s1, mode="function", env=baseenv(), ifnot=NA)[[1L]])) {
#cat("find ", s1, "\n", sep="")
         pkg=find(s1, mode="function")
         if (length(pkg)) {
            if (! "package:base" %in% pkg) {
               pkg=pkg[1L]
               if (s1 != "ginv" && pkg != "package:base") {
                  # prepend package name for external R function
                  s1=sprintf("%s_%s_r_", substring(pkg, 9L), s1)
               }
            }
         }
      }
   }
   if (s1 == "{") {
      args=sapply(st[-1L], st2arma, call2arma, sprintf("%s   ", indent), iftern, env, ...)
   } else if (s1 == "=" || s1 == "<-") {
      args=c(st2arma(st[[2L]], call2arma, indent, iftern, env, ...),
         st2arma(st[[3L]], call2arma, indent, iftern=TRUE, env, ...))
   } else {
      args=sapply(st[-1L], st2arma, call2arma, indent, iftern, env, ...)
   }
   if (s1 == "=" || s1 == "<-") {
      return(sprintf("%s%s=%s;\n", indent, args[1L], args[2L]))
   }
   nms=names(args)
   if (!is.null(nms)) {
      i=nchar(nms) > 0L
      nms[i]=sprintf('_["%s"]=', gsub('(["\\])', '\\\\\\1', nms[i]))
   }
   # data/code control
   if (s1 == "(") {
      # parenthesis operations
      return(sprintf("(%s)", args[1L]))
   }
   if (s1 == "[") {
#browser()
      # indexing operations, decrement by 1
      inu=which(sapply(st[-(1:2)], is.numeric)) # first argument is alway false, it is a symbol
      isy=which(sapply(st[-(1:2)], is.symbol)) # first argument is alway false, it is a symbol
      dims=lapply(st[-(1:2)], symdim, env, ...)
      isca=which(sapply(dims, function(d) length(d)==1L && d=="1"))
      s1=args[1L]
      args=args[-1L]
      if (all(sort(c(inu, isy)) == isca)) {
         # individual element access
         args=sprintf("(%s)-1", args)
         args[inu]=sapply(st[inu+2], function(s) s-1L)
         return(sprintf("%s.at(%s)", s1, paste(args, collapse=", ")))
      } else if (length(dims[[1L]]) > 1) {
         stop(sprintf("Index by matrix in '%s' is not yet implemented.", format(st)))
      } else {
         args[isca]=sprintf("(%s)-1", args[isca])
         args[inu]=st[inu+2]-1
         args[args==""]="span()"
         irange=which(sapply(st[-(1:2)], function(s) is.call(s) && as.character(s[[1L]]) == ":"))
         for (i in irange) {
            b=st[[i+2L]][[2L]]
            e=st[[i+2L]][[3L]]
            args[i]=sprintf("span(%s,%s)",
               if (is.numeric(b)) b-1L else paste(b,1,sep="-"),
               if (is.numeric(e)) e-1L else paste(e,1,sep="-")
            )
         }
         return(sprintf("%s(%s)", s1, paste(args, collapse=", ")))
      }
   }
   if (s1 == ":") {
      # sequence operations
      if (is.numeric(st[[2L]]) && is.numeric(st[[3L]])) {
         len=as.integer(abs(st[[3L]]-st[[2L]])+1)
      } else {
         len=sprintf("abs(%s-%s)+1", args[2L], args[1L])
      }
      return(sprintf("linspace<ivec>(%s, %s, %s)", args[1L], args[2L], len))
   }
   if (s1 == "for") {
      if (is.call(st[[3L]]) && as.character(st[[3L]][[1L]]) == ":") {
         # Here only integer counter on integer range "begin:end" will work
         begend=sapply(st[[3L]][-1L], st2arma, call2arma, indent, iftern=TRUE, env, ...)
         counter=args[1L]
         return(sprintf("%sfor (int incr_arma_=(%s <= %s ? 1 : -1), %s=%s; %s != %s+incr_arma_; %s+=incr_arma_) %s",
            indent, begend[1L], begend[2L], counter, begend[1L],
            counter, begend[2L], counter, args[3L]))
      } else {
         # use c++11 construct
         return(sprintf("%sfor (%s : %s) %s", indent, args[1L], args[2L], args[3L]))
      }
   }
   if (s1 == "while") {
      return(sprintf("%swhile (%s) %s", indent, args[1L], args[2L]))
   }
   if (s1 == "return") {
      return(sprintf("return wrap(%s)", args[1L]))
   }
   if (s1 == "if") {
      if (iftern) {
         # ternary operator a?b:c
         return(sprintf("(%s ? %s : %s)", args[1L], args[2L], args[3L]))
      } else {
#browser()
         res=sprintf("%sif (%s) %s", indent, args[1L], args[2L])
         if (length(args) == 3L) {
            res1=substring(res, 1, nchar(res)-1)
            res=sprintf("%s%s else %s%s",
               res1, if (substring(res1, nchar(res1)) != "}") ";" else "",
               args[3L],
               if (substring(args[3L], nchar(args[3L])) != "\n") "\n" else "")
         } else {
            # end with ;\n if needed
            res=sprintf("%s%s", res,
               if (substring(res, nchar(res)) != "\n") ";\n" else "")
         }
      }
      return(res)
   }
   if (s1 == "{") {
      return(sprintf("{\n%s%s}\n", 
         paste(args, collapse=""), indent))
   }
   if (s1 == "$") {
      # list element by name
      res=sprintf('%s["%s"]', args[1L], args[2L])
      obj=eval(st, env=env)
      t=rtype2rcpparma(typeof(obj))
      d=svm(obj)
      res=sprintf("as<%s>(%s)", get_decl(t, d)["arma"], res)
      return(res)
   }
   
   dims=lapply(st[-1L], symdim, env, ...)
   lens=sapply(dims, length)
   lenst=length(st)
#browser()
   
   # binary (bilateral) operations
   if (s1 %in% c("*", "/", "+", "-", ">", "<", ">=", "<=", "==", "!=", "&&", "||") && lenst == 3L) {
      # plain binary operations term by term
      # mat,vec operations '*'->'%'
      if (s1 == "*" && dims[[1L]][1L] != "1" && dims[[2L]][1L] != "1") {
         s1="%"
      }
      imin=imax=0
      if (lens[1L] > lens[2L]) {
         imin=2
         imax=1
      } else if (lens[1L] < lens[2L]) {
         imin=1
         imax=2
      }
      if (imin != 0 && dims[imin] != "1") {
         args[imin]=sprintf("mat(vec(%s).begin(), %s.n_rows, %s.n_cols, false)",
            args[imin], args[imax], args[imax])
      }
      return(sprintf("%s %s %s", args[1L], s1, args[2L]))
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
   }
   if (s1 == "%o%") {
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
   }
   
   # unary operations
   if (s1 == "+" || s1 == "-" && lenst == 2L) {
      return(sprintf("%s%s", s1, args[1L]))
   }
   
   # data creation
   if (s1 == "c" && all(lens == 1L) && all(unlist(dims)=="1")) {
      # scalars are put in a vector
      #return(sprintf("vec(NumericVector::create(%s))", paste(nms, args, sep="", collapse=", ")))
      return(sprintf("vec({%s})", paste(nms, args, sep="", collapse=", ")))
   }
   if (s1 == "c" && lenst > 2) {
      # concat arguments in a vector
      obj=eval(st, env=env)
      return(sprintf("as<%svec>(c_r_(%s))", rtype2rcpparma(typeof(obj))["arma"],
         paste(nms, args, sep="", collapse=", ")))
   }
   
   if (s1 %in% c("c", "as.numeric", "as.double", "as.integer") && lenst == 2L) {
#browser()
      # vectorization of a single argument
      if (lens[1L] > 1L) {
         args=sprintf("vectorise(%s)", args[1L])
      } else if (dims[[1L]] == "1") {
         args=sprintf("as_scalar(%s)", args[1L])
      } else {
         args=args[1L]
      }
      switch(s1,
         "as.integer"=return(sprintf("conv_to<ivec>::from(%s)", args)),
         return(sprintf("conv_to<vec>::from(%s)", args))
      )
   }
   if (s1 == "rbind") {
      rtype=rtype2rcpparma(typeof(eval(st[[2L]], env=env)))["arma"]
      # prepare the first matrix
      if (length(dims[[1L]]) == 1L && dims[[1L]] == "1") {
         res=sprintf("%smat(1,1).fill(%s)", rtype, args[1L])
      } else if (length(dims[[1L]]) == 1L) {
         res=sprintf("(%s).st()", args[1L])
      } else {
         # return the only matrix as is
         res=args[1L]
      }
      if (length(args) == 1L) {
         return(res)
      }
      for (i in 2L:length(args)) {
         a=args[i]
         if (length(dims[[i]]) == 1L) {
            a=sprintf("(%s).st()", a)
         }
         atype=rtype2rcpparma(typeof(eval(st[[i+1L]], env=env)))["arma"]
         if (atype != rtype) {
            if (atype == "cx_" || atype == "") {
               # convert res
               rtype=atype
               res=sprintf("join_vert(conv_to<%smat>::from(%s), %s)", atype, res, a)
            } else {
               # convert a
               res=sprintf("join_vert(%s, conv_to<%smat>::from(%s))", rtype, res, a)
            }
         } else {
            # simply join
            res=sprintf("join_vert(%s, %s)", res, a)
         }
      }
      return(res)
   }
   if (s1 == "cbind") {
#browser()
      rtype=rtype2rcpparma(typeof(eval(st[[2L]], env=env)))["arma"]
      # prepare the first matrix
      if (length(dims[[1L]]) == 1L && dims[[1L]] == "1") {
         res=sprintf("%smat(1,1).fill(%s)", rtype, args[1L])
      } else if (length(dims[[1L]]) == 1L) {
         res=sprintf("%smat(%s)", rtype, args[1L])
      } else {
         # return the first matrix as is
         res=args[1]
      }
      if (length(args) == 1L) {
         return(res)
      }
      for (i in 2L:length(args)) {
         a=args[i]
         atype=rtype2rcpparma(typeof(eval(st[[i+1L]], env=env)))["arma"]
         if (atype != rtype) {
            if (atype == "cx_" || atype == "") {
               # convert res
               rtype=atype
               res=sprintf("join_horiz(conv_to<%smat>::from(%s), %s)", atype, res, a)
            } else {
               # convert a
               res=sprintf("join_horiz(%s, conv_to<%smat>::from(%s))", rtype, res, a)
            }
         } else {
            # simply join
            res=sprintf("join_horiz(%s, %s)", res, a)
         }
      }
      return(res)
   }
   if (s1 == "list") {
      # convert vec in args to NumericVector
      argconv=sapply(seq_along(args), function(i) if (lens[i] > 1L || dims[[i]][1L]=="1") args[i] else sprintf("NumericVector(%s.begin(), %s.end())", args[i], args[i]))
      return(sprintf("List::create(%s)", paste(nms, argconv, sep="", collapse=", ")))
   }
   if (s1 == "rep") {
      obj=eval(st, env=env)
      t=typeof(obj)
      return(sprintf("as<%svec>(rep_r_(%s))", rtype2rcpparma(t)["arma"],
         paste(nms, args, sep="", collapse=", ")))
   }
   
   # mono argument functions
   if (s1 == "diag") {
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
            return(sprintf("diagmat(vec(%s).fill(%s))", a2, a1))
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
   }
   
   if (s1 == "nrow") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.n_rows", args[1L]))
      }
      return(sprintf("(%s).n_rows", args[1L]))
   }
   if (s1 == "ncol") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.n_cols", args[1L]))
      }
      return(sprintf("(%s).n_cols", args[1L]))
   }
   if (s1 == "dim") {
      return(sprintf("ivec(IntegerVector::create(%s.n_rows, %s.n_cols))", args[1L], args[1L]))
   }
   if (s1 == "length") {
      return(sprintf("%s.size()", args[1L]))
   }
   if (s1 == "t") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.st()", args[1L]))
      }
      return(sprintf("(%s).t()", args[1L]))
   }
   
   # needing argument conversion functions
   if (s1 == "norm") {
      if (lenst == 2L) {
         return(sprintf("norm(%s)", args))
      } else if (lenst == 3L) {
         if (is.character(st[[3L]])) {
            ntype=c("1"=1, "O"=1, "2"=2, "I"='"inf"', "F"='"fro"')[toupper(st[[3L]])]
            if (is.na(ntype) && toupper(st[[3L]]) == "M") {
               return(sprintf("max(abs(vectorise(%s)))", args))
            } else if (is.na(ntype)) {
               stop(sprintf("unknown norm type in '%s'", format(st)))
            }
            return(sprintf("norm(%s, %s)", args[1L], ntype))
         } else {
            stop(sprintf('The norm type in "%s" must be one of litteral "O", "I", "F", "M", "2"', format(st)))
         }
      }
   }
   if (s1 == "print") {
      if (lens[1L] == 1 && dims[[1]] == "1") {
         return(sprintf('%sRcout << "%s=" << %s << endl;\n', indent, args[1L], args[1L]))
      }
      return(sprintf('%s%s.print(Rcout, "%s=");\n', indent, args[1L], args[1L]))
   }
   if (s1 %in% names(call2arma)) {
      # simply translate function names
      res=sprintf("%s(%s)", call2arma[s1], paste(nms, args, sep="", collapse=", "))
      return(res)
   }
   # by default return st as a function call
   res=sprintf("%s(%s)", s1, paste(nms, args, sep="", collapse=", "))
   n=nchar(s1)
   if (n > 3 && substring(s1, n-2) == "_r_") {
      # add as<>() converter for R call
      # get arma type of result
#browser()
      obj=eval(st, env=env)
      t=rtype2rcpparma(typeof(obj))
      d=svm(obj)
      res=sprintf("as<%s>(%s)", get_decl(t, d)["arma"], res)
   }
   return(res)
}

#' Determine if an R object is a scalar, vector or matrix
#'
#' @param obj An R object. Must be a vector or matrix.
#' @return "s" for a scalar, "v" for a vector and "m" for matrix
#'
#' @details
#' An object is considered as a scalar if it has a length 1.
svm=function(obj) {
   len=length(obj)
   if (len == 1) {
       # to be cast a scalar, it is a vector who must be of length 1, not a matrix
       return("s")
   }
   if (is.matrix(obj) || is(obj, "Matrix")) {
      return("m")
   }
   if (is.vector(obj)) {
       return("v")
   }
   # for not scalar, vector, matrix
   return(NA)
}

#' Converts R typeof() to Rcpp and Armadillo types
#'
#' @param r A string resulting from \code{typeof()} call
#' @return A named character vector of length three who's components
#' are named "r", "rcpp" and "arma".
rtype2rcpparma=function(r) {
   return(switch(r,
      "function"=c(r=r, rcpp="Function", arma=NA),
      "list"=c(r=r, rcpp="List", arma=NA),
      "character"=c(r=r, rcpp="Character", arma=NA),
      "numeric"=c(r=r, rcpp="Numeric", arma=""),
      "double"=c(r=r, rcpp="Numeric", arma=""),
      "integer"=c(r=r, rcpp="Integer", arma="i"),
      "logical"=c(r=r, rcpp="Logical", arma="i"),
      "S4"=c(r=r, rcpp="S4", arma="sp_"),
      c(r=r, rcpp="SEXP", arma=NA)
   ))
}

#' Get a string for variable declaration in RcppArmadillo code.
#' @param var An R object
#' @param env An R environment
#' @return the same as \code{\link{rtype2rcpparma}}
#' @details
#' This function is just a wrapper for \code{\link{rtype2rcpparma}}
get_vartype=function(var, env=parent.frame()) {
   # return a named character vector of r, rcpp and arma types that can be used in variable declarations
   # together with struct (one of scalar, vector, matrix if applicable)
   if (is.character(var)) var=as.symbol(var)
   obj=eval(var, env=env)
   return(rtype2rcpparma(typeof(obj)))
}

#' Get a string with variable declaration in RcppArmadillo code
#'
#' @param t A character vector as returned by \code{\link{rtype2rcpparma}}
#' @param d A character vector as returned by \code{\link{symdim}} or
#'  \code{\link{svm}}
#' @return a named vector with rcpp and arma part of declarations
get_decl=function(t, d) {
   # as scalar, vector or matrix.
   # t is a three element vector, d is a vector of symbolic dimensions
   # or a string "s", "v" or "m"
   if (t["r"] %in% c("character", "numeric", "double", "integer", "logical", "S4")) {
      if (length(d) > 1L || d == "m") {
         # matrix
         return(c(
            rcpp=sprintf("%sMatrix", t["rcpp"]),
            arma=if (!is.na(t["arma"])) sprintf("%smat", t["arma"]) else NA
         ))
      }
      if (d == "1" || d =="s") {
         # scalar
         return(c(
            rcpp=switch(t["r"],
               "character"="std:string",
               "complex"="complex",
               "numeric"="double",
               "double"="double",
               "integer"="int",
               "logical"="bool"
            ),
            arma=switch(t["r"],
               "character"="std:string",
               "complex"="complex",
               "numeric"="double",
               "double"="double",
               "integer"="int",
               "logical"="bool"
            )
         ))
      }
      # vector
      return(c(
         rcpp=sprintf("%sVector", t["rcpp"]),
         arma=if (!is.na(t["arma"])) sprintf("%svec", t["arma"]) else NA
      ))
   }
   if (is.na(t["arma"])) {
      res=c(t["rcpp"], t["rcpp"])
      names(res)=c("rcpp", "arma")
      return(res)
   }
   # by default, give just typeof
   return(t[c("rcpp", "arma")])
}

#' Get assignment operators
#'
#' @param st An R statement
#' @return A list with of assignment statements inside \code{st}
get_assign=function(st) {
   # go through the tree of statments to gather and return assignements and
   # for loop begining
   if (!is.language(st) || is.symbol(st)) {
      return(list())
   }
   s=as.character(st[[1L]])
   if (s == "=" || s == "<-") {
      return(list(st))
   }
   if (s == "for") {
      res=append(list(st[1:3]), get_assign(st[[4L]]))
      return(res)
   }
   if (s == "if") {
      res=append(list(st[1:2]), get_assign(st[[3L]]))
      if (length(st) == 4L) {
         # else block is present
         res=append(res, get_assign(st[[4L]]))
      }
      return(res)
   }
   if (s == "{") {
      res=list()
      for (i in 2:length(st)) {
         res=append(res, get_assign(st[[i]]))
      }
      return(res)
   }
   # by default, empty list
   return(list())
}
