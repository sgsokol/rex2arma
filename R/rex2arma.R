#' Translate a (simple) R code to RcppArmadillo code.
#'
#' @param text A string with a code or an R expression
#'  or an R function or a plain code.
#' @param fname A string with a name for Rcpp function (if NULL,
#'  default to rex_arma_)
#' @param exec A logical or integer from {1, 2} (default TRUE).
#' @param copy A logical saying to make or not a local copy of input
#'  vectors/matrices (default TRUE, i.e. make a copy).
#' @param rebuild A logical saying to rebuild or not previous Rcpp code
#'  (default FALSE, i.e. no rebuild if \code{text} does not change)
#' @param inpvars A character vector with names of input variables
#'  (default NULL, i.e. this list is automatically constructed)
#' @param outvars A character vector with names of output variables.
#'  If there are many, they are put in a list. If \code{outvars} is named,
#'  The names are used for naming returned list items.
#' @return a string with generated code or the result of calculation
#'  depending on \code{exec} (cf. Details)
#'
#' @details
#' If the \code{text} is a string it is parsed. If not, it must
#' be a valid R expression or a function or a plain code (cf. Examples).
#'
#' If \code{exec==FALSE}, the RcppArmadillo code is not executed
#' it is just returned as a string.
#' If \code{exec==1}, the Rcpp function is created in the calling
#'  environment but not executed.
#' If \code{exec==2} or \code{TRUE}, than the c++ function is called too
#'  and its output is returned as the result of rex2arma() call.
#'
#' If \code{copy==TRUE}, objects inside the cpp code are created
#' with memory copying.
#' If \code{copy==FLASE}, the calculations are made "in place"
#' (be carefull with that! The side effects can be very surprising).
#' The argument \code{rebuild} is passed through to cppFunction()
#'
#' If \code{text} is a function, the argument list is taken
#' from that function and \code{inpvars} is not consulted.
#'
#' If \code{outvars} is  NULL, all variables that appear on the left hand side
#' will be returned. If \code{text} is a function,
#' the output is taken from it and \code{outvars} is not consulted.
#'
#' @section Warning:
#' The converted R code is executed in a dedicated environement.
#'  So it is better to call rex2arma when input variables are small
#'  vectors/matrices.
#'
#' Input variables must be of most generic type during \code{rex2arma()} call.
#' For example, if \code{a} is supposed to be float, don't set just
#'  \code{a <- 1:2} which will be of type \code{integer}. Instead use
#'  \code{a <- 1:2+0.1} or something alike.
#
#' @examples
#' a=1:3; b=a+3; # NB. Inputs a and b are defined before a call to rex2arma()
#' # \code{text} is a string:
#' code=rex2arma("a+b", exec=F)
#' cat(code);
#'
#' # \code{text} is a function
#' f=function(a, b) a+b
#' code=rex2arma(f, exec=F)
#'
#' # \code{text} is a plain R code
#' code=rex2arma(a+b, exec=F)
#' code=rex2arma({inner=a%*%b; outer=a%o%b}, exec=F)
#'
#' # \code{text} is an expression
#' e=parse(text="{s=a+b; d=a-b}")
#' code=rex2arma(e, exec=F)
#'
#' # to execute the produced code:
#' (result=eval(parse(text=code)))
#' # or simply
#' (result=rex2arma("a+b"))

#' @section Limitations:
#' \enumerate{
#'  \item no implicit vector recycling in term-by-term operations
#'  \item symbols "T" and "F" are converted to "true" and "false"
#'  \item no global assignement \code{<<-}
#'  \item and last but not least, no guaranty that produced code works as
#'   expected even if it compiles without error
#' }
#' Allowed operators and calls are:
#'   binary: '+', '-', '*', '/', '%*%', '%o%', logical operators
#'   calls: t(), [qr.]solve(), ginv(), diag()
#'          nrow(), ncol(), norm()
#'   element-wise mathematical functions having the
#'     same syntaxe in R and Armadillo: sqrt(), abs() etc.;
#'
#' @section Code conventions:
#' R variables are considered as one of the following type (typeof(var) -> Rcpp; arma)
#' (-"- means that the type has no its own equivalent in arma and kept as in Rcpp):
#'\enumerate{
#' \item list -> List; -"-
#' \item character -> Character; -"-
#' \item numeric -> double; -"-
#' \item double -> double; -"-"
#' \item integer -> int; sword
#' \item function -> Function; -"-
#' \item logical -> Logical; -"-
#'}
#'
#' Depending on dimension of {numeric, integer, character, logical} variable
#' it can be one the following structures in Rcpp/arma:
#' Rcpp / Armadillo types:
#' - {Integer,Numeric,Complex,Character}Vector/{ivec,vec,cx_vec,-"-}
#' - {Integer,Numeric,Complex,Character}Matrix/{imat,mat,cx_mat,-"-}

rex2arma=function(text, fname=if (is.function(text)) sprintf("%s_arma_", as.character(substitute(text))) else "rex_arma_", exec=TRUE, copy=TRUE, rebuild=FALSE, inpvars=NULL, outvars=NULL) {
#browser()
   pfenv=parent.frame()
   probenv=new.env() # execute statements here to probe typeof(out)
   if (!is.null(inpvars)) {
      # populate probenv with input variables
      for (it in inpvars) {
         assign(it, get(it, env=pfenv), env=probenv)
      }
   }
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
   var_typeof=matrix(NA, nrow=0, ncol=3)
   colnames(var_typeof)=c("r", "rcpp", "arma")
   var_decl=matrix(NA, nrow=0, ncol=2)
   fun_decl=c()
   colnames(var_decl)=c("rcpp", "arma")
   indent="   "
   for (i in 1:length(e)) {
      out=c() # output var in this statement
      ins=c() # input vars -"-
      s1="" # first item in st
      st=e[[i]] # current statement
      # declare R functin calls of type package::func()
      pada=getParseData(parse(t=if (is.call(st)) format1(st) else as.character(st)))
      # function calls
      ifu=which(pada$token == "SYMBOL_FUNCTION_CALL")
      # namespace "::" operator
      ins_get=pada$id[ifu] > 1 & pada$token[ifu-1] == "NS_GET"
      # package_function names
      pkg=pada$text[ifu[ins_get]-2]
      fun=pada$text[ifu[ins_get]]
      if (length(fun)) {
         # declare the functions if not yet done
         pfun=sprintf("%s_%s_r_", pkg, fun)
         fun_decl=c(fun_decl, sprintf('%sFunction %s=Environment("package:%s")["%s"];\n',
            indent, pfun, pkg, fun))
      }
      # alone functions other than from package:base
      fun=pada$text[ifu[!ins_get]]
      pkg=sapply(fun, function(f) { p=find(f, mode="function"); if ("package:base" %in% p) "package:base" else p[1]})
      inb=pkg != "package:base"
      fun=fun[inb]
      pkg=substring(pkg[inb], 9L)
      if (length(fun)) {
         # declare the functions if not yet done
         pfun=sprintf("%s_%s_r_", pkg, fun)
         fun_decl=c(fun_decl, sprintf('%sFunction %s=Environment("package:%s")["%s"];\n',
            indent, pfun, pkg, fun))
      }
      s1=if (is.call(st)) as.character(st[[1L]]) else as.character(st)
      if (is.symbol(st) || ! is.language(st) || s1 == "return" ||
         (s1 != "=" && s1 != "<-" && s1 != "if" && s1 != "print"
         && s1 != "for" && s1 != "if" && s1 != "while")) {
         # it must be the last statement before return
         if (i != length(e)) {
            stop(sprintf("Statement '%s' is not assignement. It mus be the last one in the list", format1(st)))
         }
         # prepare return
         ret=st2arma(st, indent=indent, iftern=TRUE, env=probenv, dim_tab=var_dims)
         if (s1 != "return") {
            # return just the last expression (which is not "=" neither "<-")
            ret=sprintf("%sreturn wrap(%s);\n", indent, ret)
         } else {
            ret=sprintf("%s%s;\n", indent, ret)
         }
         # exclude from ins list member with "$"
         isy=which(pada$token=="SYMBOL")
         idol=which(pada$token=="'$'")
         if (length(idol)) {
            idol=apply(outer(isy, idol, `-`), 2L, function(v) which.max(v > 0L))
            isy=isy[-idol]
         }
         ins=gsub("\\.", "_", sort(unique(pada$text[isy])))
         # add ins to probenv
         for (it in setdiff(ins, ls(probenv))) {
            assign(it, get(it, env=pfenv), env=probenv)
         }
         # get typeof and declaration ins
         for (it in setdiff(ins, rownames(var_decl))) {
            # update var_typeof and dims
            var_typeof=rbind(var_typeof, get_vartype(it, probenv))
            var_dims[[it]]=symdim(as.symbol(it), probenv, var_dims)
            rownames(var_typeof)[nrow(var_typeof)]=it
            # update var_decl
            var_decl=rbind(var_decl, get_decl(var_typeof[it,], var_dims[[it]]))
            rownames(var_decl)[nrow(var_decl)]=it
         }
         known=rownames(var_decl)
         # store ins that are not in previous outs neither in inps
         di=setdiff(setdiff(ins, outv), inps)
         inps=c(inps, di)
      }
browser()
      #pada=getParseData(parse(t=format1(st)))
      #ieq=which(pada$token=="EQ_ASSIGN" | pada$token=="LEFT_ASSIGN" |
      #   pada$token=="IN")
      # gather out vars in assignements and for loops (because of "if" and "for" blocks, it may be many)
      leq=get_assign(st) # list of equalities
      for (eq in leq) {
         if (length(eq) == 0L) {
            next
         }
         # get ins var
         s1=as.character(eq[[1L]])
         rhs=format1(if (s1 == "if") rhs_eq(eq[[2L]]) else rhs_eq(eq[[3L]]))
         pada=getParseData(parse(t=rhs))
         # exclude from ins list member with "$"
         isy=which(pada$token=="SYMBOL")
         idol=which(pada$token=="'$'")
         if (length(idol)) {
            idol=apply(outer(isy, idol, `-`), 2L, function(v) which.max(v > 0L))
            isy=isy[-idol]
         }
         ins=gsub("\\.", "_", sort(unique(pada$text[isy])))
         # add ins to probenv
         for (it in setdiff(ins, ls(probenv))) {
            assign(it, get(it, env=pfenv), env=probenv)
         }
         # get typeof and declaration ins
         for (it in setdiff(ins, rownames(var_decl))) {
            # update var_typeof and dims
            var_typeof=rbind(var_typeof, get_vartype(it, probenv))
            var_dims[[it]]=symdim(as.symbol(it), probenv, var_dims)
            rownames(var_typeof)[nrow(var_typeof)]=it
            # update var_decl
            var_decl=rbind(var_decl, get_decl(var_typeof[it,], var_dims[[it]]))
            rownames(var_decl)[nrow(var_decl)]=it
         }
         known=rownames(var_decl)
         # store ins that are not in previous outs neither in inps
         di=setdiff(setdiff(ins, outv), inps)
         inps=c(inps, di)
         # out update
         if (s1 == "if") {
            # no out in this operator
            next
         }
         out=lhs=c()
         eq_chain=TRUE
         while (eq_chain) {
            if (!is.symbol(eq[[2L]]) && as.character(eq[[2L]][[1L]]) == "[") {
               tmp=format1(eq[[2L]][[2L]])
               attr(tmp, "lhs")=format1(eq[[2L]])
               out=c(format1(eq[[2L]][[2L]]), out)
               lhs=c(format1(eq[[2L]]), lhs)
            } else {
               tmp=format1(eq[[2L]])
               out=c(tmp, out)
               lhs=c(tmp, lhs)
            }
            eq_chain=is.call(eq[[3L]]) && (eq[[3]][[1]] == as.symbol("<-") || eq[[3]][[1]] == as.symbol("="))
            if (eq_chain)
               eq <- eq[[3L]]
         }
         out=gsub("\\.", "_", out)
         if (! all(out %in% known)) {
            outv=c(outv, out)
         }
         # rhs for probe execution
         if (s1 =="for") {
            rhs=sprintf("head(%s, 1L)", rhs)
         }
#browser()
         for (i in seq_along(out)) {
            if (out[i] %in% known)
               next
            eval(parse(t=sprintf("%s=%s", lhs[i], rhs)), env=probenv)
            var_typeof=rbind(var_typeof, get_vartype(out[i], probenv))
            rownames(var_typeof)[nrow(var_typeof)]=out[i]
            var_dims[[out[i]]]=symdim(parse(t=out[i]), probenv, var_dims)
            var_decl=rbind(var_decl, get_decl(var_typeof[out[i],], var_dims[[out[i]]]))
            rownames(var_decl)[nrow(var_decl)]=out[i]
            known=c(known, out[i])
         }
      }
#browser()
      # hart part: add a line of cpp code (outs are already declared)
      if (ret == "") {
         code=sprintf("%s%s", code, st2arma(st, indent=indent, iftern=FALSE, env=probenv, dim_tab=var_dims))
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
         # by default name all list items
         names(outvars)=outvars
      } # if some names are empty, corresponding items will be unnamed
      # colvec to vector conversion
#browser()
      outtype=var_decl[outvars,,drop=FALSE]
      outconv=outvars
      ivec=grep("vec$", outtype[, "arma"])
      outconv[ivec]=sprintf("%s(%s.begin(), %s.end())",
         outtype[ivec, "rcpp"], outvars[ivec], outvars[ivec])
      ret=sprintf("%sreturn List::create(\n%s\n   );\n", indent,
         paste('      _("', names(outvars), '")=',
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
   inptype=var_decl[inpvars,,drop=FALSE]
   do_copy= if (copy) "true" else "false"
   if (length(inpvars)) {
      iconv=grep("(vec|mat)$", inptype[,"arma"])
      inparma=inptype[iconv,"arma"]
      imat=grep("mat$", inparma)
      ivec=grep("vec$", inparma)
      decl=inpvars[iconv]
      decl[imat]=sprintf("%s%s %s(%s_in_.begin(), %s_in_.nrow(), %s_in_.ncol(), %s);\n",
         indent, inparma[imat], decl[imat], decl[imat],
         decl[imat], decl[imat], do_copy)
      decl[ivec]=sprintf("%s%s %s(%s_in_.begin(), %s_in_.size(), %s);\n",
         indent, inparma[ivec], decl[ivec], decl[ivec],
         decl[ivec], do_copy)
#browser()
      # add "_in_" to vars that will be converted to vec mat
      ivm=grep("(vec|mat)$", inptype[,"arma"])
      inp_in=inpvars
      inp_in[ivm]=sprintf("%s_in_", inpvars[ivm])
      sig=paste(inptype[,"rcpp"], inp_in, sep=" ", collapse=",\n")
   } else {
      decl=sig=""
   }
   out_decl=""
   if (length(outv)) {
      outv=sort(unique(outv))
   }
   for (out in outv) {
      out_decl=sprintf("%s%s%s %s;\n", out_decl, indent,
         var_decl[out, "arma"], out)
   }
   
   if (length(fun_decl)) {
      fun_decl=sort(unique(fun_decl))
   }
   body=sprintf('
   // auxiliary functions
   Environment base_env_r_=Environment::base_env();
   Function rep_r_=base_env_r_["rep"];
   Function c_r_=base_env_r_["c"];
   // External R function declarations
%s
   // Input variable declarations and conversion
%s
   // Output and intermediate variable declarations
%s
   // Translated code starts here
%s
%s', paste(fun_decl, collapse=""), paste(decl, collapse=""), out_decl, code, ret)
   
   code=sprintf("
cppFunction(depends='RcppArmadillo', rebuild=%s, includes='
using namespace arma;

template <typename T>
inline unsigned which_max(T v) {
   unsigned i;
   v.max(i);
   return i+1;
}

template<typename T>
inline unsigned which_min(T v) {
   unsigned i;
   v.min(i);
   return i+1;
}
', \n\"

SEXP %s(\n%s) {\n%s\n}\"\n)\n",
      rebuild, fname, sig, gsub('"', '\\\\"', body))
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
