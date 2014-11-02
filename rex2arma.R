# 2014-10-27 sokol@insa-toulouse.fr
# Licence of RcppArmadillo applies here
# Copyright 2014, INRA, France

source("rex2arma.inc.R")
rex2arma=function(text, fname="rex_arma_", exec=TRUE, copy=FALSE, rebuild=FALSE, inpvars=NULL, outvars=NULL) {
   # parse a (simple) R expression from text argument then produce and
   # optionally execute corresponding RcppArmadillo inline code.
   # If the parameter text is a string it is parsed. If not, it must
   # be a valid R expression.
   # Produced RcppArmadillo code is returned as attr 'rcpp_code' of the result expression
   # if exec==FALSE, the RcppArmadillo code is not executed
   # it is just returned as a vector of two strings. The first one
   # is for creation of the matrix, the second one for its execution
   # and setting variables in the calling environement.
   # If exec==1, than only function creation is evaluated
   # If exec==2, than the second code string is evaluated too.
   # If copy==TRUE then objects are created with memory copying.
   # The calculations are made "in place" otherwise (be carefull with that!).
   # The argument 'rebuild' is passed through to cppFunction()
   # 'inpvars' is a vector of input variables. If NULL, it is
   # calculated from the variables appearing on the right had sides (RHS)
   # of assignements.
   # 'outvars' is a vector of variable names that must be returned in the
   # result list. If NULL, all variables that appear on the left hand side (LHS)
   # will be returned. If the vector outvars is named, than names will
   # be used to name list items.
   #
   # Usage:
   # > a=diag(3); b=a; b[]=1:9; code=rex2arma("a=a+b", exec=F); cat(code);
   # to execute the produced code:
   # > eval(parse(text=code))
   # or simply
   # > rex2arma("a=a+b") # the matrix a should be modified
   
   # Limitations:
   # - expression must include only numeric matrices and vectors;
   # - no subscripting
   # - no implicit vector recycling in terme-by-term operations
   # - no global assignement '<<-'; (but take care of operations that can be done in-place)
   # - and last but not least, no garanty that produced code works as
   #   expected even if it compiles without error
   # Allowed operators and calls are:
   #   binary: '+', '-', '*', '/', '%*%'
   #   calls: t(), [qr.]solve(), ginv(), diag() (which extracts its diagonal from a matrix)
   #   element-wise mathematical functions having the
   #     same syntaxe in R and Armadillo: sqrt(), abs() etc.;
   
   pfenv=parent.frame()
   if (is.character(text)) {
      e=parse(text=text)
   } else if (is.expression(text)) {
      e=text
      text=paste(sprintf("%s", e), collapse="; ")
   } else {
      stop("Argument 'text' must be a string with R code or an expression")
   }
#browser()
   # make outs (to return to the caller) and inps (to be received from the caller)
   # the first occurens of an out is prepended by 'mat ' or 'vec ' or 'double '
   # declaration depending on the dim() of out
   outv=c()
   inps=c()
   code=""
   ret=""
   var_dims=list()
   for (i in 1:length(e)) {
      out="" # output var in this statement
      s1="" # first item in st
      prep="" # prepend by type mat, vec or double for first occurence in out
      st=e[[i]] # current statement
      if (is.symbol(st) || is.numeric(st) ||
         as.character(s1 <- st[[1]])=="return" ||
         (s1 != "=" && s1 != "<-")) {
         # it must be the last statement before return
         if (i != length(e)) {
            stop(sprintf("Statement '%s' is not assignement and is not\nthe last one in the list"))
         }
         # prepare return
         ret=st2arma(st, env=pfenv, dim_tab=var_dims)
         if (s1 != "return") {
            # return just the last expression (which is not "=" neither "<-")
            ret=sprintf("return wrap(%s);\n", ret)
         } else {
            ret=sprintf("%s;\n", ret)
         }
      }
      if (s1 == "=" || s1 == "<-") {
         # gather out var
         out=as.character(st[[2L]])
         outv=c(outv, out)
         pada=getParseData(parse(t=format(st[-(1L:2L)])))
      } else {
         pada=getParseData(parse(t=format(st)))
      }
      # gather ins
      ins=sort(unique(pada$text[pada$token=="SYMBOL"]))
      for (it in setdiff(ins, names(var_dims))) {
         # update var_dims
         var_dims[[it]]=symdim(parse(t=it)[[1L]], var_dims, pfenv)
      }
      # prepend by 'mat ', 'vec ' or 'double '
      if (out != "" && all(out != names(var_dims))) {
         d1=symdim(st[[3]], var_dims)
         var_dims[[out]]=d1
         if (d1[1L] == "1") {
            prep="double "
         } else if (length(d1) == 1L) {
            prep="vec "
         } else {
            prep="mat "
         }
      }
#browser()
      # hart part: add a line of cpp code
      if (ret == "") {
         code=sprintf("%s\n   %s%s;", code, prep, st2arma(st, env=pfenv, dim_tab=var_dims))
      }
      # store ins that are not in previous outs neither in inps
      di=setdiff(setdiff(ins, outv), inps)
      inps=c(inps, di)
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
      ret=sprintf("return List::create(\n%s\n   );\n",
         paste('      Named("', names(outvars), '")=',
         outvars, sep='', collapse=",\n"))
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
   decl=ifelse(inptype["arma",]=="mat",
      paste("   mat ", inpvars, "(", inpvars, "_in_.begin(), ", inpvars, "_in_.nrow(), ", inpvars, "_in_.ncol(), ", do_copy, ");", sep=""),
      ifelse (inptype["arma",]=="vec", paste("   vec ", inpvars, "(", inpvars, "_in_.begin(), ",inpvars, "_in_.size(), ", do_copy, ");", sep=""),
      paste("   double ",inpvars, "=", inpvars, "_in_;", sep=""))
   )
   
   sig=paste(inptype["rcpp",], " ", inpvars, "_in_", sep="", collapse=",\n")
   body=sprintf("   using namespace arma;
   using namespace Rcpp;
   // Variable declarations
%s
   // Translated code starts here
%s
   
   %s", paste(decl, collapse="\n") , code, ret)
   
   code=sprintf("
cppFunction(depends='RcppArmadillo', rebuild=%s,\n'SEXP %s(\n%s) {\n%s\n}'\n)\n",
      rebuild, fname, sig, body)
   if (isTRUE(exec)) {
      # create function in the parent frame
      # call it with params from the parent frame
      eval(parse(text=code), env=pfenv)
      res=do.call(fname, lapply(inpvars, get, env=pfenv), env=pfenv)
      return(res)
   } else {
      return(code)
   }
}
