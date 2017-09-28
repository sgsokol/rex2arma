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
   if (is.null(st)) {
      return("0")
   }
   if (!is.call(st) || as.character(st[[1L]]) == "$") {
      # we are ready to return the dim vector
      s1=gsub("\\.", "_", format1(st))
      if (s1 %in% names(dim_tab)) {
         # already known
         return(dim_tab[[s1]])
      } else if (s1 == "") {
         return("")
      } else {
         obj=eval(st, envir=env)
         len=length(obj)
         if (len == 1) {
             # to be cast a scalar, it must be a vector of length 1, not a matrix
             return("1")
         }
         if (is.matrix(obj) || methods::is(obj, "Matrix")) {
            return(c(sprintf("nrow(%s)", s1), sprintf("ncol(%s)", s1)))
         }
         if (is.vector(obj)) {
             return(sprintf("length(%s)", s1))
         }
         # for not scalar, vector, matrix
         return(NA)
      }
   }
   if (is.numeric(st) || is.logical(st) || is.character(st) || is.complex(st)) {
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
   if (s1 == "ifelse" || (s1 == "seq_along" && lenst == 2)) {
      # return the dims of the condition
      return(symdim(st[[2L]], env, dim_tab))
   }
   # dims of arguments
   dims=lapply(st[-1], symdim, env, dim_tab)
   lens=sapply(dims, length)
   if (any(s1 == c("(", "abs", "sqrt", "sin", "cos", "tan", "sinpi", "cospi", "tanpi", "atan", "exp", "sinh", "cosh", "tanh", "asin", "acos", "atan", "is.finite", "is.infinite", "is.na", "is.nan", "is.symbol", "ceiling", "floor", "!"))) {
      # pass through the dims of the first argument
      # this must be before startsWith(s1, "is.") which returns 1.
      return(dims[[1L]])
   }
   if (any(s1 == c("+", "-", "*", "/", "^", ">", "<", "<=", ">=", "!=", "==", "&", "|")) && lenst == 3L) {
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
   if (s1 == "if") {
      # return the longest dim of two parts
      if (lenst < 4)
         return(dims[[2]])
      d1=dims[[2L]]
      l1=length(d1)
      d2=dims[[3L]]
      l2=length(d2)
      return(if (l2 > l1) d2 else d1)
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
   if (s1 %in% c("sum", "prod", "nrow", "ncol", "length", "mean", "min", "max", "sd", "norm", "identical", "&&", "||") || startsWith(s1, "is.")) {
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
         alast=format1(st[[lenst]])
         return(c(alast, alast))
      }
   }
   
   # data/code control
   if (s1 == "seq" || s1 == "seq.int") {
      if (lenst == 1)
         return(format1(st[[2]]))
   }
   if (s1 == ":") {
#browser()
      # range size
      if (s1 == ":" && identical(st[[2L]], st[[3L]])) {
         # identical begin and end => a scalar
         return("1")
      }
      # otherwise a vector
      return(sprintf("length(%s)", format1(st)))
   }
   if (s1 == "rep" || s1 == "rep.int") {
      # repeated vector
      argv=sapply(st[-1L], format1)
      if (lenst == 2)
         return(argv[1])
      # match arguments
      marg=as.list(match.call(args(get(s1)), st))[-1]
      times=marg$times
      each=marg$each
      x=marg$x
      if (is.null(times)) {
         # times will be the first unnamed argument after x removing
         itimes=which(nchar(names(marg))==0L)[1+!("x" %in% names(marg))]
         times=marg[[itimes]]
         if (is.null(times))
            times=1
      }
      if (is.null(each))
         each=1
      ix=which("x" == names(marg))
      res=if (dims[[ix]] == "1" && times == 1 && each == 1) "1" else sprintf("(%s)*(%s)*(%s)", dims[[ix]], format1(times), format1(each))
      return(res)
   }
   
   if (s1 == "[") {
#browser()
      argv=sapply(st[-1L], format1)
      # indexing operator
      if (max(lens[-1L]) > 1L) {
         # index is matrix => result of indexing is a vector
         stopifnot(length(argv) == 2L) # only one matrix in index
         return(sprintf("nrow(%s)", argv[2L]))
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
   if (s1 == "head") {
      if (lenst == 2L) {
         # only one argument
         res=sprintf("min(%s, 6)", dims[[1L]][1L])
      } else {
         # the length of head can be positive or negative
         a1=format1(st[[2L]])
         a2=format1(st[[3L]])
         res=sprintf("if ((%s) > 0) min(%s, NROW(%s)) else max(0, (%s)+(%s))", a2, a2, a1, dims[[1]][1], a2)
      }
      if (lens[1L] == 2L) {
         # the first argument is matrix => head() operates on rows
         res=c(res, dims[[1L]][2L])
      }
      return(res)
   }
   if (nchar(s1) > 1 && any(substring(s1, 1, 1) == c("d", "p", "q", "r")) && any(substring(s1, 2) == c("norm", "gamma", "beta", "cauchy", "chisq", "exp", "geom", "f", "hyper", "lnorm", "multinom", "nbinom", "pois", "t", "unif", "weibull", "signrank", "tukey", "wilcox"))) {
      if (lens[1L] == 1L && dims[1] == "1") {
         return("1")
      } else {
         return(sprintf("length(%s)", format1(st[[2]])))
      }
   }
   
   # by default
   # eval and return the dim of result
   res=symdim(eval(st, envir=env), env, dim_tab)
   warning(sprintf("Couldn't figure out dimensions of '%s' so just eval-ed", format1(st)))
   return(res)
}

# translate R scalar types to C scalar types
rsca2csca=c(integer="int", double="double", logical="bool")
#' Translate an elementary R statement to RcppArmadillo
#'
#' @param st A statement to transalte
#' @param call2a A named character vector giving Armadillo equivalence
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

# function that translates just by name, i.e. the arguments are passed trough as they are
call2a=as.environment(list("qr.solve"="solve", "ginv"="pinv", "^"="pow", "stop"="stop",
   "ceiling"="arma::ceil", "floor"="arma::floor", "round"="arma::round", "trunc"="arma::trunc",
   "which.min"="which_min", "which.max"="which_max",
   "Re"="real", "Im"="imag", "integer"="ivec", "double"="vec"))
st2arma=function(
   st,
   call2arma=call2a,
   indent="",
   iftern=FALSE,
   env=parent.frame(),
   ...
) {
   if (is.null(st)) {
      return("R_NilValue")
   }
   if (typeof(st) == "logical" && is.na(st)) {
      return("datum::nan")
   }
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
      if (!(s1 %in% names(call2arma)) && !is.function(mget(s1, mode="function", envir=baseenv(), ifnotfound=NA)[[1L]])) {
#cat("find ", s1, "\n", sep="")
         pkg=utils::find(s1, mode="function")
         if (length(pkg)) {
            if (! ("package:base" %in% pkg || "package:stats" %in% pkg)) {
               pkg=pkg[1L]
               if (s1 != "ginv" && s1 != "head" && pkg != "package:base") {
                  # prepend package name for external R function
                  s1=sprintf("%s_%s_r_", gsub("\\.", "_", sub("package:", "", pkg)), s1)
               }
            }
         }
      }
   }
   # jump control
   if (s1 == "break")
      return(s1)
   if (s1 == "next")
      return("continue")
#browser()
   lenst=length(st)
   if (s1 == "{") {
      argv=sapply(st[-1L], st2arma, call2arma, indent, iftern, env, ...)
      argv=sub("^([^ ]+)", sprintf("%s\\1", indent), argv) # indent the code which is not yet
      argv=sub("^$", indent, argv)
   } else if (s1 == "=" || s1 == "<-") {
      argv=c(st2arma(st[[2L]], call2arma, indent, iftern, env, ...),
         st2arma(st[[3L]], call2arma, indent, iftern=TRUE, env, ...))
   } else {
      # put arguments in the order of the function definition
      argv=sapply(as.list(match.call(args(s1), st))[-1L], st2arma, call2arma, sprintf("%s   ", indent), iftern, env, ...)
   }
   # trim spaces for arguments of "if" and "="
   if (s1 %in% c("=", "<-", "if", "for", "while")) {
      argv=sapply(argv, sub, pattern="^\\s*", replacement="")
   }
   if (s1 == "=" || s1 == "<-") {
      # populate env with results of this operator
      eval(st, envir=env)
      rhs=eval(st[[3]], envir=env)
      if (is.function(rhs) && is.call(st[[3]]) && deparse(st[[3]][[1]]) == "function") {
         # jusr skip it. Local functions must be declared and defined out of the function body
         return("")
      }
      svm_r=svm(rhs)
      if (is.call(st[[3]]) && st[[3]][[1]] != as.symbol("<-") && st[[3]][[1]] != as.symbol("=") && svm_r == "s") {
         # rhs is a scalar
         lhs=try(eval(st[[2]], envir=env))
         use_fill = !(inherits(lhs, "try-error") || length(lhs) == 1)
#if (use_fill) browser()
         return(sprintf("%s%s%s(%s)", indent, argv[1L], if (use_fill) ".fill" else "=", argv[2L]))
      } else {
         return(sprintf("%s%s=%s", indent, argv[1L], argv[2L]))
      }
   }
   nms=names(argv)
   if (!is.null(nms)) {
      i=nchar(nms) > 0L
      nms[i]=sprintf('_["%s"]=', gsub('(["\\])', '\\\\\\1', nms[i]))
   }
   # data/code control
   if (s1 == "(") {
      # parenthesis operations
      return(sprintf("(%s)", argv[1L]))
   }
   if (s1 == "[") {
#browser()
      # indexing operations, decrement by 1
      inu=sapply(st[-(1:2)], is.numeric)
      isy=sapply(st[-(1:2)], is.symbol)
      dims=lapply(st[-(1:2)], symdim, env, ...)
      isca=sapply(dims, function(d) length(d)==1L && d=="1")
      s1=argv[1L]
      argv=argv[-1L]
      if (length(dims[[1L]]) > 1) {
         if (length(dims[[1L]]) == 2 && length(argv)==1) {
            # elements are indexed by a two column matrix (rows,cols)
            return(sprintf("at_rc(%s, (%s)-1)", format(s1), argv[1]))
         } else {
            stop(sprintf("Index by matrix in an array '%s' is not yet implemented.", format1(st)))
         }
      } else {
         # [row,col] access
         argv[!inu]=ifelse(isca[!inu], sprintf("(%s)-1", argv[!inu]), sprintf("conv_to<uvec>::from(%s)-1", argv[!inu]))
         if (any(inu))
            argv[inu]=sapply(which(inu)+2, getElement, object=st)-1
         argv[argv==""]="span()"
         irange=which(sapply(st[-(1:2)], function(s) is.call(s) && as.character(s[[1L]]) == ":"))
         for (i in irange) {
            b=st[[i+2L]][[2L]]
            e=st[[i+2L]][[3L]]
            argv[i]=sprintf("span(%s,%s)",
               if (is.numeric(b)) b-1L else paste(b,1,sep="-"),
               if (is.numeric(e)) e-1L else paste(e,1,sep="-")
            )
         }
         #return(sprintf("vectorise(%s(%s))", s1, paste(argv, collapse=", ")))
         if (is.call(st[[3]]) && deparse(st[[3]][[1]]) %in% c(">", "<", ">=", "<=", "==", "!=", "!"))
            argv=sprintf("find(%s)", paste(argv, collapse=", "))
         return(sprintf("%s(%s)", s1, paste(argv, collapse=", ")))
      }
   }
   if (s1 == "head") {
      dims=lapply(st[-1L], symdim, env, ...)
      lens=sapply(dims, length)
      if (lens[1] == 1L) {
         obj=eval(st[[2L]], envir=env)
         t=rtype2rcpparma(typeof(obj))
         d=svm(obj)
         # vector argument
         if (lenst == 2L) {
            # default length == 6 for the only vector argument
            return(sprintf("%s(%s.begin(), std::min(6, (int) (%s).n_elem))", get_decl(t, d)["arma"], argv[1L], argv[1L]))
         } else {
            return(sprintf("%s(%s.begin(), std::max(0, std::min((int) (%s >= 0 ? %s : (%s).n_elem+(%s)), (int) (%s).n_elem)))", get_decl(t, d)["arma"], argv[1L], argv[2L], argv[2L], argv[1L], argv[2L], argv[1L]))
         }
      } else {
         # matrix argument
         if (lenst == 2L) {
            # default length == 6 for the only matrix argument
            return(sprintf("(%s)(0,0,size(6, (%s).n_cols))", argv[1L], argv[1L]))
         } else {
            return(sprintf("(%s)(0,0,size((%s)>=0 ? std::min(%s, (int) (%s).n_rows) : std::max(0, (int) (%s).n_rows+(%s)), (%s).n_cols))", argv[1L], argv[2L], argv[2L], argv[1L], argv[1L], argv[2L], argv[1L]))
         }
      }
   }
   if (s1 == "seq_len") {
      return(sprintf("linspace<ivec>(1, %s, %s)", argv[1L], argv[1L]))
   }
   if (s1 == "seq_along") {
      return(sprintf("linspace<ivec>(1, (%s).n_elem, (%s).n_elem)", argv[1L], argv[1L]))
   }
   if (s1 == ":" || s1 == "seq" || s1 == "seq.int") {
#browser()
      # sequence operations
      if (lenst == 2) {
         return(sprintf("regspace<ivec>(1, %s)", argv[1L]))
      } else if (lenst == 3) {
         return(sprintf("regspace<ivec>(%s, %s)", argv[1L], argv[2L]))
      } else {
         stdarg=formals(args(s1))
         marg=as.list(match.call(get(s1), st))[-1]
         allarg=utils::modifyList(stdarg, marg)
         if (!is.null(allarg$by))
            return(sprintf("regspace<vec>(%s, %s, %s)", argv[1L], allarg$by, argv[2L]))
         else if (!is.null(allarg[["length.out"]]))
            return(sprintf("linspace<vec>(%s, %s, %s)", argv[1L], argv[2L], allarg[["length.out"]]))
         else
            stop(sprintf("Cannot yet manage this call '%s'", format1(st)))
      }
   }
   if (s1 == "for") {
#browser()
      if (substring(argv[3L], 1, 1) == "{") {
         spacer=" "
      } else {
         spacer=sprintf("\n%s   ", indent)
      }
      if (is.call(st[[3L]]) && as.character(st[[3L]][[1L]]) == ":") {
         # Here only integer counter on integer range "begin:end" will work
         begend=sapply(st[[3L]][-1L], st2arma, call2arma, indent, iftern=TRUE, env, ...)
         counter=argv[1L]
         return(sprintf("%sfor (int incr_arma_=(%s <= %s ? 1 : -1), %s=%s; %s != %s+incr_arma_; %s+=incr_arma_)%s%s",
            indent, begend[1L], begend[2L], counter, begend[1L],
            counter, begend[2L], counter, spacer, argv[3L]))
      } else {
         # use c++11 construct
         return(sprintf("%sfor (%s : %s)%s%s", indent, argv[1L], argv[2L], spacer, argv[3L]))
      }
   }
   #if (s1 == "while") {
   #   return(sprintf("%swhile (%s) %s", indent, argv[1L], argv[2L]))
   #}
   if (s1 == "return") {
      return(sprintf("return (%s)", argv[1L]))
   }
   if (s1 %in% c("if", "while")) {
#browser()
      if (s1 == "if" && iftern) {
         # ternary operator a?b:c
         return(sprintf("(%s ? (%s) : (%s))", argv[1L], argv[2L], if (length(argv) < 3) "R_NilValue" else argv[3L]))
      } else {
         f1=substring(argv[2L], 1, 1)
         if (f1 == "{")
            res=sprintf("%s%s (%s) %s", indent, s1, argv[1L], argv[2L])
         else
            res=sprintf("%s%s (%s)\n%s   %s", indent, s1, argv[1L], indent, argv[2L])
         if (s1 == "if" && length(argv) == 3L) {
            f2=substring(argv[3L], 1, 1)
            sc=if (f1 == "{") "" else ";"
            if (f2 == "{")
               res=sprintf("%s%s\n%selse %s", res, sc, indent, argv[3L])
            else
               res=sprintf("%s%s\n%selse\n%s   %s", res, sc, indent, indent, argv[3L])
         }
         # else {
         #   # end with ;\n if needed
         #   res=sprintf("%s%s", res,
         #      if (substring(res, nchar(res)) != "\n") ";\n" else "")
         #}
      }
      return(res)
   }
   if (s1 == "ifelse") {
      return(sprintf("%s ? %s : %s", argv[1L], argv[2L], if (length(argv) < 3) "R_NilValue" else argv[3L]))
   }
   if (s1 == "{") {
#browser()
      return(sprintf("{\n%s;\n%s}", 
         paste(argv, collapse=";\n"), sub("   ", "", indent, fixed=TRUE)))
   }
   if (s1 == "$") {
      # list element picked by its name
      res=sprintf('%s["%s"]', argv[1L], argv[2L])
      obj=eval(st, envir=env)
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
   if (s1 %in% c("*", "/", "+", "-", ">", "<", ">=", "<=", "==", "!=", "&&", "||", "&", "|") && lenst == 3L) {
      # plain binary operations term by term
      # mat,vec operations '*'->'%'
      if (s1 == "*" && dims[[1L]][1L] != "1" && dims[[2L]][1L] != "1") {
         s1="%"
      }
      if (s1 == "/") {
         # care off integer division
         pref=if (dims[[1L]][1L] == "1") "(double) " else "conv_to<mat>::from"
         argv[1]=sprintf("(%s(%s))", pref, argv[1])
      }
      #if (s1 == "&")
      #   s1="&&"
      #if (s1 == "|")
      #   s1="||"
      imin=imax=0
      if (lens[1L] > lens[2L]) {
         imin=2
         imax=1
      } else if (lens[1L] < lens[2L]) {
         imin=1
         imax=2
      }
      if (imin != 0 && dims[imin] != "1") {
         argv[imin]=sprintf("mat(vec(%s).begin(), %s.n_rows, %s.n_cols, false)",
            argv[imin], argv[imax], argv[imax])
      }
      return(sprintf("%s %s %s", argv[1L], s1, argv[2L]))
   }
   if (s1 == "identical") {
      # cannot figure out smth better than plain '=='
      return(sprintf("(%s) == (%s)", argv[1L], argv[2L]))
   }
   if (any(s1==c("%*%", crossprod, tcrossprod))) {
      # dot products
      a1=argv[1L]
      d1=dims[[1L]]
      l1=lens[1L]
      if (lenst == 3L) {
         a2=argv[2L]
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
            return(sprintf("%s*%s.t()", argv[1L], argv[2L]))
         } else {
            return(sprintf("%s*%s", argv[1L], argv[2L]))
         }
      } else if (lens[[2L]] == 1 && dims[[2L]] == "1") {
         # the second argument is a scalar
         return(sprintf("%s*%s", argv[1L], argv[2L]))
      }
      # all other cases
      return(sprintf("kron(vectorise(%s), vectorise(%s).t())", argv[1L], argv[2L]))
   }
   
   # unary operations
   if (s1 == "+" || s1 == "-" || s1 == "!" && lenst == 2L) {
      return(sprintf("%s(%s)", s1, argv[1L]))
   }
   
   # data creation
   if (s1 == "c" && all(lens == 1L) && all(unlist(dims)=="1")) {
      # scalars are put in a vector
      #return(sprintf("vec(NumericVector::create(%s))", paste(nms, argv, sep="", collapse=", ")))
      return(sprintf("vec({%s})", paste(nms, argv, sep="", collapse=", ")))
   }
   if (s1 == "c" && lenst > 2) {
      # concat arguments in a vector
      obj=eval(st, envir=env)
      if (is.list(obj))
         return(sprintf("c_r_(%s)", paste(nms, argv, sep="", collapse=", ")))
      else
         return(sprintf("as<%svec>(c_r_(%s))", rtype2rcpparma(typeof(obj))["arma"], paste(nms, argv, sep="", collapse=", ")))
   }
   
   if (s1 %in% c("c", "as.numeric", "as.double", "as.integer") && lenst == 2L) {
#browser()
      # vectorization of a single argument
      if (lens[1L] > 1L) {
         argv=sprintf("vectorise(%s)", argv[1L])
      } else if (dims[[1L]] == "1") {
         return(sprintf(if (s1 == "as.integer") "((int) as_scalar(%s))" else "as_scalar(%s)", argv[1L]))
      } else {
         argv=argv[1L]
      }
      if (s1 == "as.integer")
         return(sprintf("conv_to<ivec>::from(%s)", argv))
      else
         return(sprintf("conv_to<vec>::from(%s)", argv))
   }
   if (s1 == "rbind") {
      rtype=rtype2rcpparma(typeof(eval(st[[2L]], envir=env)))["arma"]
      # prepare the first matrix
      if (length(dims[[1L]]) == 1L && dims[[1L]] == "1") {
         res=sprintf("%smat(1,1).fill(%s)", rtype, argv[1L])
      } else if (length(dims[[1L]]) == 1L) {
         res=sprintf("(%s).st()", argv[1L])
      } else {
         # return the only matrix as is
         res=argv[1L]
      }
      if (length(argv) == 1L) {
         return(res)
      }
      for (i in 2L:length(argv)) {
         a=argv[i]
         if (length(dims[[i]]) == 1L) {
            a=sprintf("(%s).st()", a)
         }
         atype=rtype2rcpparma(typeof(eval(st[[i+1L]], envir=env)))["arma"]
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
      rtype=rtype2rcpparma(typeof(eval(st[[2L]], envir=env)))["arma"]
      # prepare the first matrix
      if (length(dims[[1L]]) == 1L && dims[[1L]] == "1") {
         res=sprintf("%smat(1,1).fill(%s)", rtype, argv[1L])
      } else if (length(dims[[1L]]) == 1L) {
         res=sprintf("%smat(%s)", rtype, argv[1L])
      } else {
         # return the first matrix as is
         res=argv[1]
      }
      if (length(argv) == 1L) {
         return(res)
      }
      for (i in 2L:length(argv)) {
         a=argv[i]
         atype=rtype2rcpparma(typeof(eval(st[[i+1L]], envir=env)))["arma"]
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
      # convert vec in argv to NumericVector
      argconv=sapply(seq_along(argv), function(i) if (lens[i] > 1L || dims[[i]][1L]=="1") argv[i] else sprintf("NumericVector(%s.begin(), %s.end())", argv[i], argv[i]))
      return(sprintf("List::create(%s)", paste(nms, argconv, sep="", collapse=", ")))
   }
   if (s1 == "rep" || s1 == "rep.int") {
      if (lenst == 2)
         return(sprintf(if (dims[[1]] == "1") "%s" else "vectorise(%s)", argv))
#browser()
      marg=as.list(match.call(args(get(s1)), st))[-1]
      obj=eval(st, envir=env)
      times=argv["times"]
      if (is.na(times)) {
         itimes=which(nchar(names(marg)) == 0)[1+!("x" %in% names(marg))]
         times=if (is.na(itimes)) "1" else argv[itimes]
      }
      each=if ("each" %in% names(marg)) argv["each"] else "1"
      if (lenst <= 4 && length(dims[[1]]) == 1) {
         return(sprintf("_rex_arma_rep(%s, %s, %s)", argv["x"], times, each))
      }
      t=typeof(obj)
      return(sprintf("as<%svec>(rep_r_(%s))", rtype2rcpparma(t)["arma"],
         paste(nms, argv, sep="", collapse=", ")))
   }
   if (s1 == "matrix") {
      # normalize arguments
      marg=as.list(match.call(matrix, st))[-1]
      svm_data=svm(eval(marg[["data"]], envir=env))
      argv=sapply(marg, st2arma, call2arma, indent, iftern, env, ...)
      if (is.null(marg[["data"]])) {
         stop(sprintf("no 'data' argument in the call '%s'", format1(st)))
      }
      rtype=rtype2rcpparma(typeof(eval(marg[["data"]], envir=env)))["arma"]
      if (is.null(marg[["ncol"]])) {
         # need to calculate ncol (default value=1 may be not good))
         argv["ncol"]=sprintf("(%s).n_elem/(%s)", argv[1], argv["nrow"])
      }
      if (!is.null(marg[["byrow"]])) {
         # inverse nrow and ncol, then transpose the matrix
         tmp=argv["nrow"]
         argv["nrow"]=argv["ncol"]
         argv["ncol"]=tmp
         transp=sprintf(".%st()", if (rtype == "complex") "s" else "")
      } else {
         transp=""
      }
      if (svm_data == "s") {
         # data is a scalar
         res=sprintf("%smat(%s, %s).fill(%s)%s", rtype, argv["nrow"], argv["ncol"], argv["data"], transp)
      } else {
         # data is a vector or matrix
         res=sprintf("%smat(%s).resize(%s, %s)%s", rtype, argv["data"], argv["nrow"], argv["ncol"], transp)
      }
      return(res)
   }
   
   # mono argument functions
   if (s1 == "diag") {
      if (lenst == 2L) {
         # only one argument
         a1=argv[1L]
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
         a1=argv[1L]
         l1=lens[1L]
         a2=argv[2L]
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
         return(sprintf("%s.n_rows", argv[1L]))
      }
      return(sprintf("(%s).n_rows", argv[1L]))
   }
   if (s1 == "ncol") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.n_cols", argv[1L]))
      }
      return(sprintf("(%s).n_cols", argv[1L]))
   }
   if (s1 == "dim") {
      return(sprintf("{(int) %s.n_rows, (int) %s.n_cols}", argv[1L], argv[1L]))
   }
   if (s1 == "length") {
      return(sprintf("%s.n_elem", argv[1L]))
   }
   if (s1 == "t") {
      if (is.symbol(st[[2]])) {
         return(sprintf("%s.st()", argv[1L]))
      }
      return(sprintf("(%s).t()", argv[1L]))
   }
   
   # needing argument conversion functions
   if (s1 == "norm") {
      if (lenst == 2L) {
         return(sprintf("norm(%s)", argv))
      } else if (lenst == 3L) {
         if (is.character(st[[3L]])) {
            ntype=c("1"=1, "O"=1, "2"=2, "I"='"inf"', "F"='"fro"')[toupper(st[[3L]])]
            if (is.na(ntype) && toupper(st[[3L]]) == "M") {
               return(sprintf("max(abs(vectorise(%s)))", argv))
            } else if (is.na(ntype)) {
               stop(sprintf("unknown norm type in '%s'", format1(st)))
            }
            return(sprintf("norm(%s, %s)", argv[1L], ntype))
         } else {
            stop(sprintf('The norm type in "%s" must be one of litteral "O", "I", "F", "M", "2"', format1(st)))
         }
      }
   }
   if (s1 == "solve") {
      return(sprintf("solve(%s, solve_opts::fast)", paste0(argv, collapse=", ")))
   }
   if (s1 == "print") {
      if (lens[1L] == 1 && dims[[1]] == "1") {
         return(sprintf('%sRcout << "%s=" << %s << endl;\n', indent, argv[1L], argv[1L]))
      }
      return(sprintf('%s%s.print(Rcout, "%s=");\n', indent, argv[1L], argv[1L]))
   }
   # probability functions
   if (nchar(s1) > 1 && any(substr(s1, 1, 1) == c("d", "p", "q", "r")) && any(substring(s1, 2) == c("gamma", "exp", "norm", "beta", "cauchy", "chisq", "geom", "f", "hyper", "lnorm", "multinom", "nbinom", "pois", "t", "unif", "weibull", "signrank", "tukey", "wilcox"))) {
#browser()
      # matched call
      mc=as.list(match.call(get(s1, mode="function"), st))[-1]
      mc[]=argv
      if (!is.null(mc[["rate"]])) {
         # inverse the rate
         mc[["rate"]]=sprintf("1./(%s)", mc[["rate"]])
      } else if (!is.null(mc[["scale"]])) {
         mc[["rate"]]=mc[["scale"]]
         mc[["scale"]]=NULL
      }
      vnm=c("p", "q", "x")
      for (xnm in vnm) {
         if (is.null(mc[[xnm]]))
            next
         #mc[[xnm]]=sprintf("as<NumericVector>(wrap(%s))", mc[[xnm]])
         break
      }
      # formal argv
      fa=utils::modifyList(lapply(formals(s1), format1), mc)
      fa[fa == "FALSE"]="false"
      fa[fa == "TRUE"]="true"
      fa["scale"]=NULL
      if (dims[[1]] != 1) {
         fa[[1]]=sprintf("conv_to<mat>::from(%s)", fa[[1]])
      }
      #r_pref=if (substr(s1, 1, 1) == "r") "" else "Rcpp::" # generator is not from R::
      r_pref="(vec) "
      postf="" #if (dims[[1]][[1]] == "1") "[0]" else ""
      res=sprintf("%s%s(%s)%s", r_pref, s1, paste(fa, sep="", collapse=", "), postf)
      return(res)
   }
#browser()
   if (s1 %in% names(call2arma)) {
      # simply translate function names
      res=sprintf("%s(%s)", call2arma[[s1]], paste(nms, argv, sep="", collapse=", "))
      return(res)
   }
   # by default return st as a function call
   res=sprintf("%s(%s)", s1, paste(nms, argv, sep="", collapse=", "))
   n=nchar(s1)
   if (n > 3 && substring(s1, n-2) == "_r_") {
      # add as<>() converter for R call
      # get arma type of result
#browser()
      obj=eval(st, envir=env)
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
   if (is.matrix(obj) || methods::is(obj, "Matrix")) {
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
      "environment"=c(r=r, rcpp="Environment", arma=NA),
      "function"=c(r=r, rcpp="Function", arma=NA),
      "list"=c(r=r, rcpp="List", arma=NA),
      "complex"=c(r=r, rcpp="Complex", arma="cx_"),
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
get_vartype=function(obj) {
   # return a named character vector of r, rcpp and arma types that can be used in variable declarations
   # together with struct (one of scalar, vector, matrix if applicable)
   #if (is.character(var)) var=parse(t=var)
   #obj=eval(var, envir=env)
   if (all(is.na(obj))) {
      return(rtype2rcpparma("double"))
   } else {
      return(rtype2rcpparma(typeof(obj)))
   }
}

#' Get a string with variable declaration in RcppArmadillo code
#'
#' @param t A character vector as returned by \code{\link{rtype2rcpparma}}
#' @param d A character vector as returned by \code{\link{symdim}} or
#'  \code{\link{svm}}
#' @return a named vector with rcpp and arma part of declarations
get_decl=function(t, d) {
   # as scalar, vector or matrix.
   # t is a three element vector as returned by get_vartype(), d is a vector of symbolic dimensions
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
               "character"="std::string",
               "complex"="complex",
               "numeric"="double",
               "double"="double",
               "integer"="int",
               "logical"="bool",
               "environment"="Environment"
            ),
            arma=switch(t["r"],
               "character"="std::string",
               "complex"="complex",
               "numeric"="double",
               "double"="double",
               "integer"="int",
               "logical"="bool",
               "environment"=NA
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
   # for/while loop begining
   if (!is.call(st)) {
      return(list())
   }
   s=as.character(st[[1L]])
   if (s == "=" || s == "<-") {
      return(c(get_assign(st[[3]]), list(st)))
   }
   if (s == "for") {
      res=append(list(st[1L:3L]), get_assign(st[[4L]]))
      return(res)
   }
   if (s == "while") {
      res=append(list(st[1L:2L]), get_assign(st[[3L]]))
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

#' format an R object in one string
#'
#' @param obj An R object
#' @return A character string of length 1 representing the object \code{obj}
format1=function(obj) {
   res=if (is.call(obj)) format(obj) else if (is.character(obj)) sprintf("\"%s\"", obj) else as.character(obj)
   if (length(res) > 1L) {
      res=paste(res, collapse="\n")
   }
   return(res)
}

#' recursivelly get right hand side of assignement "<-" or "=" which may be chained
#'
#' @param st An R statement
#' @return An R statement
rhs_eq=function(st) {
   if (is.call(st) && (st[[1L]] == as.symbol("<-") || st[[1L]] == as.symbol("="))) {
      return(rhs_eq(st[[3L]]))
   } else {
      return(st)
   }
}

#' translate C name nm_c to its original R counterpart (if found in dictionary dc2r, elsewhere leave it as is)
#'
#' @param nm_c (string) C variable name
#' @param dc2r (string) named vector of R names (named by C counterparts)
#' @return A string, R variable name

c2r=function(nm_c, dc2r=dc2r) {nm_r=dc2r[nm_c]; ifelse (is.na(nm_r), nm_c, nm_r)}
