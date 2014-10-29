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
   # it is just returned as the result string.
   # If copy==TRUE then objects are created with memory copying.
   # The calculations are made "in place" otherwise.
   # The argument 'rebuild' is passed as is to cppFunction()
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
   # - no implicit vector recycling
   # - no global assignement '<<-'; (but take care of operations that can be done in-place)
   # - allowed operators and calls are: '+', '-', '*', '/', '%*%', t(),
   #   [qr.]solve(), ginv(), diag() (which extracts its diagonal from a matrix)
   #   and element-wise mathematical functions having the same syntaxe in R and Armadillo
   #   like sqrt(), abs() etc.;
   # - dot product of two vectors x%*%x maut be written as t(x)%*%x
   # - there must be at least one assignement operator in the expression;
   
   # class translator
   class2arma=c("matrix"="mat", "numeric"="colvec", "integer"="colvec")
   class2rcpp=c("matrix"="NumericMatrix", "numeric"="NumericVector", "integer"="NumericVector")
   # call translator
   call2arma=c("%*%"="*", "*"="%", "qr.solve"="solve", "ginv"="pinv",
      "t"="", "<-"="=", "nrow"="", "ncol"="", "diag"="diagvec"
   )
   
   if (is.character(text)) {
      e=parse(text=text)
   } else if (is.expression(text)) {
      e=text
      text=paste(sprintf("%s", e), collapse="; ")
   } else {
      stop("Argument 'text' must be a string with R code or an expression")
   }
   pada=getParseData(e) # parsed data
   # check that a left assign is in the expression
   stopifnot("LEFT_ASSIGN" %in% pada$token || "EQ_ASSIGN" %in% pada$token)
   # check that no global assignement
   i=pada$token=="LEFT_ASSIGN"
   if (any(i) && pada$text[i]=="<<-") {
      stop("A global assignement '<<-' is not admitted in expressions.")
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
      iend=if (i==length(ieq)) nrow(pada) else ieq[i+1]
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
   
   code=sprintf("require(Rcpp)
cppFunction(depends='RcppArmadillo', rebuild=%s,\n   'List %s(\n%s) {\n%s\n   }'\n)\n",
      rebuild, fname, sig, body)
   # create function in the parent frame
   # call it in the parent frame
   # add code attribute
   call_code=sprintf('%s\n%s.res<-%s(%s)\nattr(res, "rcpp_code")<-"%s"',
      code, fname, fname, paste(inpvars, sep="", collapse=", "), gsub('"', '\\\\"', code))
   if (exec) {
      eval(parse(text=call_code))
      res=get(sprintf("%s.res", fname))
      # assign result vars
      for (it in names(res)) {
         assign(it, res[[it]], parent.frame())
      }
      return(res)
   } else {
      return(call_code)
   }
}
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
