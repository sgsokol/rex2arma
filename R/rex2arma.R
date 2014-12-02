# 2014-11-03
# Author: Serguei Sokol, sokol@insa-toulouse.fr
# Licence of RcppArmadillo applies here
# Copyright 2014, INRA, France

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
   # - double -> double; -"-"
   # - integer -> int; sword
   # - function -> Function; -"-
   # - logical -> Logical; -"-
   # Depending on dimension of {numeric, integer, character, logical} variable
   # it can be one the following structures in Rcpp/arma:
   # - c++ scalar/c++ scalar
   # - {Numeric,Integer,Complex,Character}Vector/{vec,ivec,cx_vec,-"-}
   # - {Numeric,Integer,Complex,Character}Matrix/{mat,imat,cx_mat,-"-}
#browser()
   pfenv=parent.frame()
   probenv=new.env() # execute statements here to probe typeof(out)
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
      pada=getParseData(parse(t=if (is.call(st)) format(st) else as.character(st)))
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
         (s1 != "=" && s1 != "<-" && s1 != "if" && s1 != "print")) {
         # it must be the last statement before return
         if (i != length(e)) {
            stop(sprintf("Statement '%s' is not assignement. It mus be the last one in the list", format(st)))
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
#browser()
      #pada=getParseData(parse(t=format(st)))
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
         if (s1 == "if") {
            rhs=if (is.call(eq[[2L]])) format(eq[[2L]]) else as.character(eq[[2L]])
         } else {
            rhs=if (!is.call(eq[[3L]])) as.character(eq[[3L]]) else
               format(eq[[3L]])
         }
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
         if (!is.symbol(eq[[2L]]) && as.character(eq[[2L]][[1L]]) == "[") {
            out=as.character(eq[[2L]][[2L]])
         } else {
            out=as.character(eq[[2L]])
         }
         out=gsub("\\.", "_", out)
         if (! out %in% known) {
            outv=c(outv, out)
         }
         if (! out %in% known) {
#browser()
            # rhs for probe execution
            if (s1 =="for") {
               rhs=sprintf("head(%s, 1L)", rhs)
            }
            eval(parse(t=sprintf("%s=%s", out, rhs)), env=probenv)
            var_typeof=rbind(var_typeof, get_vartype(out, probenv))
            rownames(var_typeof)[nrow(var_typeof)]=out
            var_dims[[out]]=symdim(as.symbol(out), probenv, var_dims)
            var_decl=rbind(var_decl, get_decl(var_typeof[out,], var_dims[[out]]))
            rownames(var_decl)[nrow(var_decl)]=out
            known=c(known, out)
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
