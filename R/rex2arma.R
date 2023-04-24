#' R Expression to RcppArmadillo
#' 
#' Translate a (simple) R code to RcppArmadillo code.
#'
#' @param text A string with a code or an R expression
#'  or an R function or a plain code.
#' @param fname A string with a name for Rcpp function (if NULL,
#'  default to rex_arma_)
#' @param exec An integer from {-1, 0, 1, 2} (default 0).
#' @param copy A logical saying to make or not a local copy of input
#'  vectors/matrices (default TRUE, i.e. make a copy).
#' @param rebuild A logical saying to rebuild or not previous Rcpp code
#'  (default FALSE, i.e. no rebuild if \code{text} does not change)
#' @param inpvars A character vector with names of input variables
#'  (default NULL, i.e. this list is automatically constructed)
#' @param outvars A character vector with names of output variables.
#'  If there are many, they are put in a list. If \code{outvars} is named,
#'  The names are used for naming returned list items.
#' @param verbose Logical scalar. It also  passed through to the Rcpp:sourceCppp() call.
#' @param includes Character vector containing code included in generated
#'  C++ code before the targeted function.
#' @param call2a An environment defining translation of R functions to C++ ones.
#' @param envir An environment where variable for parameter prototyping lives.
#' @return a string with generated code or the result of calculation
#'  depending on \code{exec} (cf. Details)
#'
#' @details
#' If the \code{text} is a string it is parsed. If not, it must
#' be a valid R expression or a function or a plain code (cf. Examples).
#'
#' If \enumerate{
#'  \item \code{exec==0}: (default) the full RcppArmadillo code is not executed
#'    it is just returned as a string;
#'  \item \code{exec==-1}: the same as 0 but only translated code is returned,
#'    i.e. without preamble, includes and Rcpp call;
#'  \item \code{exec==1}: the Rcpp function is created in the calling
#'    environment but not executed.
#'  \item \code{exec==2}: the c++ function is created and called too.
#'    Its output is returned as the result of rex2arma() call.
#' }
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
#' code=rex2arma("a+b", exec=0)
#' cat(code);
#'
#' # \code{text} is a function
#' f=function(a, b) a+b
#' code=rex2arma(f, exec=0)
#'
#' # \code{text} is a plain R code
#' code=rex2arma(a+b, exec=0)
#' code=rex2arma({inner=a%*%b; outer=a%o%b}, exec=0)
#'
#' # \code{text} is an expression
#' e=parse(text="{s=a+b; d=a-b}")
#' code=rex2arma(e, exec=0)
#'
#' # to execute the produced code:
#' #(result=eval(parse(text=code)))
#' # or simply
#' #(result=rex2arma("a+b", exec=2))

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
#' @export
rex2arma=function(text, fname=if (is.function(text) && is.symbol(substitute(text))) sprintf("%s_arma_", deparse(substitute(text))) else "rex_arma_", exec=0, copy=TRUE, rebuild=FALSE, inpvars=NULL, outvars=NULL, verbose=FALSE, includes=character(), call2a=rex2arma::call2a, envir=parent.frame()) {
#browser()
   mcall=match.call()
   if (verbose)
      cat("matched call=", format(mcall), "\n", sep="")
   fname=gsub(".", "_", fname, fixed=TRUE)
   verbose=as.logical(verbose)
   pfenv=list2env(as.list(envir, all.names=TRUE), parent=parent.frame())
   #cat("pfenv:", ls(pfenv), "\n")
   probenv=new.env() # execute statements here to probe typeof(out)
   dc2r=c() # dictionary of variables having different names in C++ and R, e.g. x.a(R)->x_a(C++)
   # populate probenv with input variables
   probenv=list2env(as.list(pfenv, all.names=TRUE)[inpvars], parent=parent.frame())
   ftype="SEXP"
   e=if (is.function(text)) text else substitute(text)
   retobj=NULL
   if (is.symbol(e)) {
      e=eval(e, envir=pfenv)
      retobj=eval(e, envir=pfenv)
   }
   if (is.character(e)) {
      e=parse(text=text)
      retobj=eval(e, envir=pfenv)
   }
   if (is.function(e)) {
      inpvars=sapply(names(formals(e)), gsub, pattern=".", replacement="_", fixed=TRUE)
#print(call2a)
      if (is.function(text) && is.symbol(substitute(text)))
         call2a[[deparse(substitute(text))]] <- fname
#print(call2a)
#stop()
      #anames=names(utils::head(as.list(args(e)), -1))
      #retobj=do.call(e, structure(lapply(anames, function(a) get(a, envir=pfenv)), names=anames))
      e=body(e)
      retobj=eval(e, envir=pfenv)
      if (inherits(e, "{")) {
         e=e[-1L]
      } else {
         e=list(e)
      }
   }
   if (is.null(retobj)) {
      retobj=eval(e, envir=pfenv)
   }
   # determine the return type
   t=rtype2rcpparma(typeof(retobj))
#browser()
   if (is.list(retobj)) {
      # list()
      ftype="List"
   } else if (is.vector(retobj)) {
      # e.g. "double" (for scalar) or "vec")
      ftype=if (length(retobj) == 1) rsca2csca[t["r"]] else sprintf("%svec", t["arma"])
   } else if (is.matrix(retobj) || is.array(retobj)) {
      # e.g. "double" (for scalar) or "mat")
      ftype=if (length(retobj) == 1) rsca2csca[t["r"]] else sprintf("%smat", t["arma"])
   } # in all other cases the return type is SEXP
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
   loc_inps=c()
   code=""
   ret=""
   var_dims=list()
   var_typeof=matrix(NA, nrow=0, ncol=3)
   colnames(var_typeof)=c("r", "rcpp", "arma")
   var_decl=matrix(NA, nrow=0, ncol=2)
   fun_decl=c()
   lfun_decl=c()
   colnames(var_decl)=c("rcpp", "arma")
   indent="   "
   fcode=""
   for (i in 1:length(e)) {
      out=c() # c output var in this statement
      ins=c() # c input vars -"-
      s1="" # first item in st
      st=e[[i]] # current statement
      sts=if (is.call(st)) format1(st) else as.character(st)
#cat("st=", format(st), "\n", sep="")
#cat("sts=", sts, "\n", sep="")
      # declare R function calls of type package::func()
#cat("parse=", format(parse(text=sts)), "\n", sep="")
      pada=utils::getParseData(parse(text=sts, keep.source=TRUE))
#cat("pada=", format(pada), "\n", sep="")
      # function calls
      ifu=which(pada$token == "SYMBOL_FUNCTION_CALL")
      # namespace "::" operator
      ins_get=pada$id[ifu] > 1 & pada$token[ifu-1] == "NS_GET"
      # package_function names
      pkg=pada$text[ifu[ins_get]-2]
      fun=pada$text[ifu[ins_get]]
#browser()
      if (length(fun)) {
         # declare the functions if not yet done
         pfun=sprintf("%s_%s_r_", pkg, fun)
         fun_decl=c(fun_decl, sprintf('%sFunction %s=Environment("package:%s")["%s"];\n',
            indent, pfun, pkg, fun))
      }
      # declare only functions other than from package:base and package:stats
      # for these two we rely on rcpp sugar
      fun=pada$text[ifu[!ins_get]]
      pkg=sapply(fun, function(f) { p=utils::find(f, mode="function"); if ("package:base" %in% p) "package:base" else p[1]})
      inb=pkg != "package:base" & pkg != "package:stats" & !(fun %in% names(call2a)) & !is.na(pkg)
      fun=fun[inb]
      pkg=pkg[inb]
      if (length(fun)) {
#browser()
         # declare the functions if not yet done
         pfun=sprintf("%s_%s_r_", sub("\\.", "_", sub("package:", "", pkg)), fun)
         fun_decl=c(fun_decl, sprintf('%sFunction %s=Environment("%s")["%s"];\n',
            indent, pfun, pkg, fun))
      }
      s1=if (is.call(st)) format1(st[[1]]) else format1(st)
      if (s1 == "<-")
         s1="="
      if (s1 == "return" || i == length(e)) {
#browser()
         # exclude from ins list member with "$"
         isy=which(pada$token=="SYMBOL") # indexes of symbols (i.e. variable names)
         idol=which(pada$token=="'$'") # indexes of dollar sign
         if (length(idol)) {
            idol=apply(outer(isy, idol, `-`), 2L, function(v) which.max(v > 0L))
            isy=isy[-idol]
         }
         insr=unique(pada$text[isy])
#cat("insr=", insr, "\n", sep="")
         ins=gsub(".", "_", insr, fixed=TRUE)
         if (!is.null(insr))
            names(insr)=ins
         dc2r=c(dc2r, insr[ins!=insr])
         # add ins to probenv
         for (it in setdiff(insr, ls(probenv))) {
            assign(it, get(it, envir=pfenv), envir=probenv)
         }
         # get typeof and declaration ins
         for (it in setdiff(ins, rownames(var_decl))) {
            obj=probenv[[c2r(it, dc2r)]]
            # update var_typeof and dims
            var_typeof=rbind(var_typeof, get_vartype(obj))
            var_dims[[it]]=symdim(as.symbol(c2r(it, dc2r)), probenv, var_dims)
            rownames(var_typeof)[nrow(var_typeof)]=it
            # update var_decl
            var_decl=rbind(var_decl, get_decl(var_typeof[it,], var_dims[[it]]))
            rownames(var_decl)[nrow(var_decl)]=it
         }
         known=rownames(var_decl)
         # store ins that are not in previous outs neither in inps
         di=setdiff(setdiff(ins, outv), inps)
         inps=c(inps, di)
         # prepare return
#browser()
         ret=st2arma(st, call2arma=call2a, indent=indent, iftern=TRUE, env=probenv, dim_tab=var_dims)
         if (s1 != "return") {
            # return just the last expression
            if (s1 == "=" || s1 == "<-") {
               ret=sprintf("%s;\n%sreturn %s;\n", ret, indent, st2arma(st[[2]], call2arma=call2a, indent="", iftern=TRUE, env=probenv, dim_tab=var_dims))
            } else if (s1 == "assign") {
               ret=sprintf("%s;\n%sreturn %s;\n", ret, indent, st2arma(eval(st[[2]], envir=probenv), call2arma=call2a, indent="", iftern=TRUE, env=probenv, dim_tab=var_dims))
            } else {
               ret=sprintf("%sreturn %s;\n", indent, ret)
            }
         } else {
            ret=sprintf("%s%s;\n", indent, ret)
         }
      }
#browser()
      #pada=utils::getParseData(parse(text=format1(st)))
      #ieq=which(pada$token=="EQ_ASSIGN" | pada$token=="LEFT_ASSIGN" |
      #   pada$token=="IN")
      
      # gather out vars in assignements and for loops (because of "if", "for" and "while" blocks, it may be many)
      leq=get_assign(st) # list of all equalities in the st
      for (eq in leq) {
         if (length(eq) == 0L) {
            next
         }
         # get ins var
         eq1=as.character(eq[[1L]])
         rhs_is_fun=length(eq) > 2 && is.call(eq[[3]]) && as.character(eq[[3]][[1]]) == "function" # fixme: function() can be at the end of assignment chain, not just one "f <- function() ..."
         loc_inps=c(loc_inps, if (rhs_is_fun) names(eq[[3]][[2]]) else c())
         rhsc=if (eq1 == "if" || eq1 == "while") rhs_eq(eq[[2L]]) else rhs_eq(eq[[3L]])
         rhs=format1(rhsc)
         # get variables which can be used in indexing on the lsh
         lhsi=if (is.call(st[[2]]) && as.character(st[[2]][[1]]) %in% c("[", "[[")) format1(st[[2]][[3]]) else character(0)
         pada=utils::getParseData(parse(text=paste(rhs, lhsi, sep=";"), keep.source=TRUE))
         # exclude from ins list member with "$"
         isy=which(pada$token == "SYMBOL")
         idol=which(pada$token == "'$'")
         if (length(idol)) {
            idol=apply(outer(isy, idol, `-`), 2L, function(v) which.max(v > 0L))
            isy=isy[-idol]
         }
         insr=unique(pada$text[isy])
         ins=gsub(".", "_", insr, fixed=TRUE)
         names(insr)=ins
         dc2r=c(dc2r, insr[ins!=insr])
         # add ins to probenv
         for (it in setdiff(insr, ls(probenv))) {
            assign(it, get(it, envir=pfenv), envir=probenv)
         }
         # get typeof and declaration ins
         for (it in setdiff(ins, rownames(var_decl))) {
            # update var_typeof and dims
            var_typeof=rbind(var_typeof, get_vartype(probenv[[c2r(it, dc2r)]]))
            var_dims[[it]]=symdim(as.symbol(c2r(it, dc2r)), probenv, var_dims)
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
         if (eq1 == "if" || eq1 == "while") {
            # no out in this operator
            next
         }
         out=lhs=c()
         eq_chain=TRUE
         while (eq_chain) {
            if (!is.symbol(eq[[2L]]) && as.character(eq[[2L]][[1L]]) == "[") {
               tmp=format1(eq[[2L]][[2L]])
               tmpr=format1(eq[[2L]])
               tmpc=gsub(".", "_", tmpr, fixed=TRUE)
               names(tmpr)=tmpc
               dc2r=c(dc2r, tmpr[tmpc!=tmpr])
               attr(tmp, "lhs")=tmpc
               out=c(tmp, out)
               lhs=c(tmpc, lhs)
            } else {
               tmp=format1(eq[[2L]])
               out=c(tmp, out)
               lhs=c(tmp, lhs)
            }
            eq_chain=is.call(eq[[3L]]) && (eq[[3]][[1]] == as.symbol("<-") || eq[[3]][[1]] == as.symbol("="))
            if (eq_chain)
               eq <- eq[[3L]]
         }
         outr=out
         out=gsub(".", "_", outr, fixed=TRUE)
         names(outr)=out
         dc2r=c(dc2r, outr[out!=outr])
         # rhs for probe execution
         if (eq1 == "for") {
            rhs=sprintf("utils::head(%s, 1L)", rhs)
         }
#browser()
         for (i in seq_along(out)) {
            if (out[i] %in% known)
               next
            lhsi=c2r(lhs[i], dc2r)
            probenv[[lhsi]]=eval(rhsc, envir=probenv)
            if (eq1 == "for")
               probenv[[lhsi]]=utils::head(probenv[[lhsi]], 1)
            #eval(parse(text=sprintf("%s=%s", lhsi, if (eq1 == "for") sprintf("(%s)[[1]]", rhs) else rhs)), envir=probenv)
            obj=probenv[[lhsi]]
            is_fun=FALSE
            if (is.function(obj)) {
               is_fun=TRUE
               obj=do.call(obj, lapply(names(utils::head(as.list(args(obj)), -1)), get, envir=probenv))
            }
            if (is_fun) {
#browser()
               ret_rex=do.call("rex2arma", list(text=as.symbol(lhsi), exec=-1, fname=out[i], copy=copy, includes=lfun_decl, rebuild=rebuild, verbose=verbose, call2a=call2a, envir=probenv), envir=probenv);
               # really create function for future prob calls
               #do.call("rex2arma", c(list(eq[[3]]), exec=1, fname=out[i], copy=copy, includes=lfun_decl, rebuild=rebuild, verbose=verbose, call2a=call2a, envir=envir))
               lfun_decl=c(lfun_decl, sprintf("%s %s(%s);\n%s",
                  ret_rex["ftype"], ret_rex["fname"], ret_rex["fsig"],
                  ret_rex["fcode"]))
            } else {
               var_typeof=rbind(var_typeof, get_vartype(obj))
               rownames(var_typeof)[nrow(var_typeof)]=out[i]
               var_dims[[out[i]]]=symdim(parse(text=c2r(out[i], dc2r)), probenv, var_dims)
               var_decl=rbind(var_decl, get_decl(var_typeof[out[i],], var_dims[[out[i]]]))
               rownames(var_decl)[nrow(var_decl)]=out[i]
               known=c(known, out[i])
               outv=c(outv, out[i])
            }
         }
      }
#browser()
      # heart part: add a line of cpp code (outs are already declared)
      if (ret == "") {
#browser()
         code=sprintf("%s%s;\n", code, st2arma(st, call2arma=call2a, indent=indent, iftern=FALSE, env=probenv, dim_tab=var_dims))
      }
   }
   if (ret == "") {
      # prepare return if not yet done
      # find output vars
      outs=unique(outv) # output variables
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
#browser()
   if (!is.null(inpvars)) {
      # all inps must be in inpvars
      di=setdiff(setdiff(inps, inpvars), loc_inps)
      if (length(di)) {
         stop(sprintf("The variables '%s' has to be in the inpvars",
            paste(sapply(di, c2r, dc2r), collapse=", ")))
      }
      di=setdiff(inpvars, inps)
      if (length(di)) {
         warning(sprintf("The variables '%s' in the inpvars are superfluous", paste0(di, collapse=", ")))
         inpvars=inpvars[!(inpvars %in% di)]
      }
   } else {
      inpvars=inps
   }
#cat("inps=", inpvars, "\n", sep="")
   
   # gather rcpp code
   do_copy= if (copy) "true" else "false"
   if (length(inpvars)) {
      inptype=var_decl[inpvars,,drop=FALSE]
      iconv=grep("(vec|mat)$", inptype[,"arma"])
      inparma=inptype[iconv,"arma"]
      imat=grep("mat$", inparma)
      ivec=grep("vec$", inparma)
      decl=inpvars[iconv]
      decl[imat]=sprintf("%s%s %s(%s_in_.begin(), %s_in_.n_rows, %s_in_.n_cols, %s);\n",
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
      sig=paste(inptype[,"arma"], inp_in, sep=" ", collapse=",\n ")
   } else {
      decl=sig=""
   }
   out_decl=""
   outv=setdiff(outv, inpvars)
   for (out in outv) {
      out_decl=sprintf("%s%s%s %s;\n", out_decl, indent,
         var_decl[out, "arma"], out)
   }
   
   if (length(fun_decl)) {
      fun_decl=unique(fun_decl)
   }
   fun_body=sprintf('
   // External R function declarations
%s
   // Input variable declarations and conversion
%s
   // Output and intermediate variable declarations
%s
   // Translated code starts here
%s
%s', paste0(fun_decl, collapse=""), paste0(decl, collapse=""), out_decl, code, ret)
   fun_code=sprintf("// [[Rcpp::export]]\n%s %s(\n %s) {\n%s\n}", ftype, fname, sig, fun_body)
   code_preamble=sprintf('
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(rex2arma)]]

// make vec to look like a vector not like a one column matrix
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR 
#include <RcppArmadillo.h>

#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// includes auxiliary definitions
#include <rex2arma.h>

// Supplied includes
%s
// Local function declarations
%s
', paste0(includes, collapse="\n"), paste0(lfun_decl, collapse="\n"))
   code=sprintf("
Rcpp::sourceCpp(rebuild=%s, \nverbose=%s,\ncode=\"%s\n\")\n",
      rebuild, format1(verbose), code=gsub('"', '\\\\"', sprintf("%s\n%s", code_preamble, fun_code))
   )
#browser()
   if (exec == 2L) {
      # create function in the parent frame
      # call it with params from the parent frame
      if (verbose)
         cat("rex2arma produced code=", code, "\n", sep="")
      eval(parse(text=code), envir=pfenv)
#if (inherits(.Last.value, "try-error"))
#   browser()
      res <- try(do.call(fname, lapply(c2r(inpvars, dc2r), get, envir=pfenv), envir=pfenv), silent=TRUE)
#if (inherits(.Last.value, "try-error"))
#   browser()
      return(res)
   } else if (exec == 1L) {
      if (verbose)
         cat("rex2arma produced code=", code, "\n", sep="")
      eval(parse(text=code), envir=pfenv)
      return(code)
   } else if (exec == 0) {
      return(code)
   } else if (exec == -1) {
#browser()
      return(c(ftype=as.vector(ftype), fname=fname, fsig=sig, fcode=fun_code))
   } else {
      stop(sprintf("Unknown value for exec '%s'", format1(exec)))
   }
}
