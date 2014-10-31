# 2014-10-27 sokol@insa-toulouse.fr
# Licence of RcppArmadillo applies here
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
   # > a=diag(3); b=a; b[]=1:9; code=expr2arma("a=a+b", exec=F); cat(code);
   # to execute the produced code:
   # > eval(parse(text=code))
   # or simply
   # > expr2arma("a=a+b") # the matrix a should be modified
   
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
   
   # class translator
   class2arma=c("matrix"="mat", "numeric"="colvec", "integer"="colvec")
   class2rcpp=c("matrix"="NumericMatrix", "numeric"="NumericVector", "integer"="NumericVector")
   
   if (is.character(text)) {
      e=parse(text=text)
   } else if (is.expression(text)) {
      e=text
      text=paste(sprintf("%s", e), collapse="; ")
   } else {
      stop("Argument 'text' must be a string with R code or an expression")
   }
   
   vars=sort(unique(pada$text[pada$token=="SYMBOL"]))
   
   # translate comments
   i=pada$token=="COMMENT"
   pada$text[i]=paste("//", substring(pada$text[i], 2), sep="")
   
#browser()
   # special calls translations (nrow(), ncol())
   pada=fun2postfix("nrow", ".n_rows", pada)
   pada=fun2postfix("ncol", ".n_cols", pada)
   pada=fun2postfix("t", ".t()", pada)
   
   # plainly translate R calls to arma calls ('<-' included)
   i=pada$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL", "'*'", "LEFT_ASSIGN")
   pada$text[i]=sapply(pada$text[i], function(item) {res=call2arma[item]; if (is.na(res)) item else res})
   pada$token[pada$token=="LEFT_ASSIGN"]="EQ_ASSIGN"
   
   ieq=which(pada$token=="EQ_ASSIGN")
   # make outs (to return to the caller) and inps (to be received from the caller)
   # the first occurens of an out is prepended by 'mat ' declaration
   outv=c()
   inps=c()
   for (i in 1:length(ieq)) {
      it=ieq[i]
      # find result var (LHS) for the current assignment
      iout=max(which(pada$id < pada$id[it] & pada$token=="SYMBOL"))
      out=pada$text[iout]
      stopifnot(length(out) == 1)
      # prepend by 'mat '
      if (! out %in% outv && ! out %in% inps) {
         pada$text[iout]=paste("mat ", pada$text[iout], sep="")
      }
      # find inputs vars (RHS) for the current assignment
      iend=max(which(pada$line1<=pada$line2[it+1L] & pada$col1<=pada$col2[it+1L]))
      ins=sort(unique(pada$text[which(pada$id <= pada$id[iend] & pada$id >= pada$id[it] & pada$token=="SYMBOL")]))
      stopifnot(length(ins) > 0)
      # find ins that are not in previous outs neither in inps
      di=setdiff(setdiff(ins, outv), inps)
      inps=c(inps, di)
      outv=c(outv, out)
   }
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
   if (!is.null(inpvars)) {
      # all inps must be in inpvars
      di=setdiff(inps, inpvars)
      if (length(di)) {
         stop(sprintf("The variables '%s' has to be in the inpvars", paste(di, sep=", ", collapse="")))
      }
      di=setdiff(inpvars, inps)
      if (length(di)) {
         stop(sprintf("The variables '%s' in the inpvars are superfluous", paste(di, sep=", ", collapse="")))
      }
   } else {
      inpvars=inps
   }
   
   # append ';' to the end of expression to the right from EQ_ASSIGN
   iend=sapply(ieq+1L, function(i) max(which(pada$line2==pada$line2[i] & pada$col2==pada$col2[i])))
   pada$text[iend]=paste(pada$text[iend], ";", sep="")
   
   # append '\n   ' to the end of rows
   # (must be the last operation in pada$text modifications)
   irows=sort(unique(pada$line1))
   iend=sapply(irows, function(i) max(which(pada$line1==i)))
   pada$text[iend]=paste(pada$text[iend], "\n   ", sep="")
   
   # prepare rcpp code
   objs=lapply(inpvars, get, env=parent.frame(), mode="numeric")
   class_r=sapply(objs, class)
   class_arma=class2arma[class_r]
   class_rcpp=class2rcpp[class_r]
   size_arma=ifelse(class_r=="matrix", paste(inpvars, "_r_.nrow(), ",
      inpvars, "_r_.ncol()", sep=""),  paste(inpvars, "_r_.size()", sep=""))
   
   sig=paste(class_rcpp, " ", inpvars, "_r_", sep="", collapse=",\n")
   decl=paste("   ", class_arma, " ", inpvars, "(", inpvars, "_r_.begin(), ",
      size_arma, ", ", if (copy) "true" else "false", ");",
      sep="", collapse="\n")
   body=sprintf("   using namespace arma;
   // Variable declarations
%s
   // Translated code starts here
   %s
   return Rcpp::List::create(\n%s\n   );\n",
      decl, paste(pada$text, sep="", collapse=""),
      paste('      Rcpp::Named("', names(outvars), '")=',
      outvars, sep='', collapse=",\n"))
   
   code1=sprintf("
cppFunction(depends='RcppArmadillo', rebuild=%s,\n'List %s(\n%s) {\n%s\n}'\n)\n",
      rebuild, fname, sig, body)
   # create function in the parent frame
   # call it in the parent frame
   # add code attribute
   code2=sprintf('%s.res<-%s(%s)\nattr(res, "rcpp_code")<-"%s"',
      fname, fname, paste(inpvars, sep="", collapse=", "),
      gsub('"', '\\\\"', code1))
   if (isTRUE(exec)) {
      eval(parse(text=code1), env=parent.frame())
      eval(parse(text=code2))
      res=get(sprintf("%s.res", fname))
      # assign result vars
      for (it in names(res)) {
         assign(it, res[[it]], parent.frame())
      }
      return(res)
   } else if (exec==1L) {
      eval(parse(text=code1), env=parent.frame())
      return(c(code1, code2))
   } else {
      return(c(code1, code2))
   }
}
