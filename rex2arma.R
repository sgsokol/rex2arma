# 2014-11-03
# Author: Serguei Sokol, sokol@insa-toulouse.fr
# Licence of RcppArmadillo applies here
# Copyright 2014, INRA, France

source("rex2arma.inc.R")
rex2arma=function(text, fname="rex_arma_", exec=TRUE, copy=TRUE, rebuild=FALSE, inpvars=NULL, outvars=NULL) {
   # translate a (simple) R code (or a string having a code or an expression
   # or a function) to RcppArmadillo inline code using cppFunction().
   # It optionally executes the generated RcppArmadillo function.
   # If the parameter text is a string it is parsed. If not, it must
   # be a valid R expression or a function or a plain code.
   # if exec==FALSE, the RcppArmadillo code is not executed
   # it is just returned as a string.
   # If exec==1, the function is created but not executed.
   # If exec==2 or TRUE, than the function is called too and its output
   # is returned as the result of rex2arma() call.
   # If copy==TRUE, objects inside the cpp code are created with memory copying.
   # If copy==FLASE, the calculations are made "in place"
   # (be carefull with that! The side effects can be very surprising).
   # The argument 'rebuild' is passed through to cppFunction()
   # 'inpvars' is a vector of input variables. If NULL, it is
   # calculated from the variables appearing on the right had sides (RHS)
   # of assignements. If parameter 'text' is a function, the argument list is taken
   # from that function and inpvars is not consulted.
   # 'outvars' is a vector of variable names that must be returned in the
   # result list. If NULL, all variables that appear on the left hand side (LHS)
   # will be returned. If the vector outvars is named, than names will
   # be used to name list items. If parameter 'text' is a function,
   # the output is taken from it and outvars is not consulted.
   #
   # Usage:
   # > a=1:3; b=a+3; # NB. the input parameters must be defined before a call to rex2arma()
   # as a text:
   # > code=rex2arma("a+b", exec=F)
   # > cat(code);
   # as a function
   # > f=function(a, b) a+b
   # > code=rex2arma(f, exec=F)
   # as an R code
   # > code=rex2arma(a+b, exec=F) # NB. a+b is not executed in R,
   # # (even is passed as an argument to rex2arma()) but only in the cpp compiled code
   # or
   # > code=rex2arma({x=a+b; y=a-b}, exec=F) # NB. a+b etc. is not executed in R
   # as an expression
   # > e=parse(text="{x=a+b; y=a-b}")
   # > code=rex2arma(e, exec=F)
   # to execute the produced code:
   # > (result=eval(parse(text=code)))
   # or simply
   # > (result=rex2arma("a+b"))
   
   # Limitations:
   # - expression must include only numeric matrices and vectors;
   # - no subscripting
   # - no implicit vector recycling in terme-by-term operations
   # - no loops for, while etc.
   # - symbols "T" and "F" are converted to "true" and "false"
   # - no global assignement '<<-'; (but take care of operations that can be done in-place)
   # - and last but not least, no garanty that produced code works as
   #   expected even if it compiles without error
   # Allowed operators and calls are:
   #   binary: '+', '-', '*', '/', '%*%', '%o%'
   #   calls: t(), [qr.]solve(), ginv(), diag() (which extracts its diagonal from a matrix)
   #          nrow(), ncol()
   #   element-wise mathematical functions having the
   #     same syntaxe in R and Armadillo: sqrt(), abs() etc.;
   #
   # Code conventions:
   # R variables are considered as one of the following type (typeof(var) -> Rcpp; arma)
   # (-"- means that the type has no its own equivalent in arma and kept as in Rcpp):
   # - list -> List; -"-
   # - character -> std::string; -"-
   # - numeric -> double; -"-
   # - integer -> int; sword
   # - function -> Function; -"-
   # - logical -> int; bool
   # Depending on dimension of {numeric, integer, character, logical} variable
   # it can be one the following structures in Rcpp/arma:
   # - c++ scalar/c++ scalar
   # - {Numeric,Integer,Complex,Character}Vector/{vec,ivec,cx_vec,-"-}
   # - {Numeric,Integer,Complex,Character}Matrix/{mat,imat,cx_mat,-"-}
#browser()
   pfenv=parent.frame()
   e=substitute(text)
   if (is.symbol(e)) {
      e=eval(e, env=pfenv)
   }
   if (is.character(e)) {
      e=parse(text=text)
   }
   if (is.function(e)) {
      inpvars=names(formals(e))
      e=body(e)
      if (class(e) == "{") {
         e=e[-1L]
      } else {
         e=list(e)
      }
   }
   if (is.name(e[[1L]]) && e[[1L]]!=as.symbol("{")) {
      e=list(e)
   }
   # else if (is.expression(text) || is.call(text) || is.call(text[[1L]])) {
   #   e=text
   #} else {
   #   stop("Argument 'text' must be a string with R code or an expression")
   #}
   if (e[[1L]]==as.symbol("{")) {
      if (is.expression(e)) {
         e=e[[1L]][-1L]
      } else {
         e=e[-1L]
      }
   }
#browser()
   # make outs (to return to the caller) and inps (to be received from the caller)
   # the first occurens of an out is prepended by 'mat ' or 'vec ' or 'double ' in st2arma()
   # declaration depending on the dim() of out
   outv=c()
   inps=c()
   code=""
   ret=""
   var_dims=list()
   indent="   "
   for (i in 1:length(e)) {
      out=c() # output var in this statement
      ins=c() # input vars -"-
      s1="" # first item in st
      st=e[[i]] # current statement
      if (is.symbol(st) || is.numeric(st) ||
         as.character(s1 <- st[[1]])=="return" ||
         (s1 != "=" && s1 != "<-" && s1 != "if")) {
         # it must be the last statement before return
         if (i != length(e)) {
            stop(sprintf("Statement '%s' is not assignement. It mus be the last one in the list", st))
         }
         # prepare return
         ret=st2arma(st, indent=indent, env=pfenv, dim_tab=var_dims)
         if (s1 != "return") {
            # return just the last expression (which is not "=" neither "<-")
            ret=sprintf("%sreturn wrap(%s);\n", indent, ret)
         } else {
            ret=sprintf("%s%s;\n", indent, ret)
         }
      }
#browser()
      pada=getParseData(parse(t=format(st)))
      ieq=which(pada$token=="EQ_ASSIGN" | pada$token=="LEFT_ASSIGN")
      # gather out vars (because of "if" and "for" blocks, it may be many
      for (i in ieq) {
         out=pada$text[pada$id == pada$id[i]-1L]
         irhs=pada$line1 >= pada$line1[i+1L] &
            pada$col1 >= pada$col1[i+1L] &
            pada$line1 <= pada$line2[i+1L] &
            pada$col1 <= pada$col2[i+1L]
         # gather ins
         ins=sort(unique(pada$text[pada$token=="SYMBOL" & irhs]))
         # gather ins
         for (it in setdiff(ins, names(var_dims))) {
            # update var_dims
            var_dims[[it]]=symdim(parse(t=it)[[1L]], var_dims, pfenv)
         }
         known=names(var_dims)
         outv=c(outv, out)
         if (! out %in% known) {
            rhs=paste(pada$text[irhs], collapse="")
            d1=symdim(parse(t=rhs)[[1L]], var_dims)
            var_dims[[out]]=d1
            known=c(known, out)
            code=sprintf("%s%s%s %s;\n", code, indent,
               if (length(d1) > 1L) "mat" else if (d1 == "1") "double" else "vec", out)
         }
         # store ins that are not in previous outs neither in inps
         di=setdiff(setdiff(ins, outv), inps)
         inps=c(inps, di)
      }
#browser()
      # hart part: add a line of cpp code preceded by declarations
      if (ret == "") {
         code=sprintf("%s%s", code, st2arma(st, indent=indent, env=pfenv, dim_tab=var_dims))
      }
   }
   if (ret == "") {
      # prepare return if not yet done
      # find output vars
      outs=sort(unique(outv)) # output variables
      if (is.null(outvars)) {
         outvars=outs
      } else {
         # check that requested outvars are in outs
         di=setdiff(outvars, outs)
         if (length(di)) {
            stop(sprintf("Requested output variables '%s' are not assigned in the ''", paste(di, collapse=", "), text))
         }
      }
      # give names if needed
      if (is.null(names(outvars))) {
         names(outvars)=outvars
      } else {
         # if some names are empty, keep the value as name
         i=nchar(names(outvars))==0
         names(outvars)[i]=outvars[i]
      }
      outtype=sapply(outvars, function(v) {
         d=var_dims[[v]];
         if (length(d)==2L) {
            return(c("mat", "NumericMatrix"))
         } else if (d == "1") {
            return(c("double", "double"))
         } else {
            return(c("vec", "NumericVector"))
         }
      })
      dim(outtype)=c(2L, length(outvars))
      rownames(outtype)=c("arma", "rcpp")
      # colvec to vector conversion
#browser()
      outconv=ifelse(outtype=="vec", sprintf("NumericVector(%s.begin(), %s.end())", outvars, outvars), outvars)
      ret=sprintf("%sreturn List::create(\n%s\n   );\n", indent,
         paste('      Named("', names(outvars), '")=',
         outconv, sep='', collapse=",\n"))
   }
   if (!is.null(inpvars)) {
      # all inps must be in inpvars
      di=setdiff(inps, inpvars)
      if (length(di)) {
         stop(sprintf("The variables '%s' has to be in the inpvars",
            paste(di, collapse=", ")))
      }
      di=setdiff(inpvars, inps)
      if (length(di)) {
         stop(sprintf("The variables '%s' in the inpvars are superfluous", paste(di, sep=", ", collapse="")))
      }
   } else {
      inpvars=inps
   }
   
   # gather rcpp code
   inptype=sapply(inpvars, function(v) {
      d=var_dims[[v]];
      if (length(d)==2L) {
         return(c("mat", "NumericMatrix"))
      } else if (d == "1") {
         return(c("double", "double"))
      } else {
         return(c("vec", "NumericVector"))
      }
   })
   dim(inptype)=c(2L, length(inpvars))
   rownames(inptype)=c("arma", "rcpp")
   do_copy= if (copy) "true" else "false"
   if (length(inpvars)) {
      decl=ifelse(inptype["arma",]=="mat",
         paste("   mat ", inpvars, "(", inpvars, "_in_.begin(), ", inpvars, "_in_.nrow(), ", inpvars, "_in_.ncol(), ", do_copy, ");", sep=""),
         ifelse (inptype["arma",]=="vec", paste("   vec ", inpvars, "(", inpvars, "_in_.begin(), ",inpvars, "_in_.size(), ", do_copy, ");", sep=""),
         paste("   double ",inpvars, "=", inpvars, "_in_;", sep=""))
      )

      sig=paste(inptype["rcpp",], " ", inpvars, "_in_", sep="", collapse=",\n")
   } else {
      decl=sig=""
   }
   body=sprintf("
   using namespace arma;
   using namespace Rcpp;
   // auxiliary variables
   uword isca1, isca2;
   uvec  ivec1, ivec2;
   // Variable declarations
%s
   // Translated code starts here
%s
%s", paste(decl, collapse="\n") , code, ret)
   
   code=sprintf("
cppFunction(depends='RcppArmadillo', rebuild=%s,\n'SEXP %s(\n%s) {\n%s\n}'\n)\n",
      rebuild, fname, sig, body)
   if (isTRUE(exec) || exec==2L) {
      # create function in the parent frame
      # call it with params from the parent frame
      eval(parse(text=code), env=pfenv)
      res=do.call(fname, lapply(inpvars, get, env=pfenv), env=pfenv)
      return(res)
   } else if (exec == 1L) {
      eval(parse(text=code), env=pfenv)
      return(code)
   } else {
      return(code)
   }
}
