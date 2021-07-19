#' glinvci: Confidence intervals and hypothesis testing for GLInv model
#'
#' The glinvci package provides a framework for computing the maximum-likelihood estimates
#' and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes,
#' including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic
#' comparative methods. The framework is designed to be flexible enough that the user can
#' easily specify their own parameterisation and obtain the maximum-likelihood estimates and
#' confidence intervals of their own parameters.
#' 
#' @author Hao Chi Kiang, \email{hello@hckiang.com}
#' @docType package
#' @name glinvci
#' @useDynLib glinvci, .registration=TRUE
#' @import plyr
#' @import Rcgmin
#' @import lbfgsb3c
#' @import BB
#' @import ape
#' @import stats
NULL


ndesc = function (x, ...) UseMethod('ndesc')
ndesc.glinv_gauss = function (self) .Call(Rndesc, self$ctree)

nparams = function (x, ...) UseMethod('nparams')
nparams.glinv_gauss = function (self) .Call(Rnparams, self$ctree)

#' Check if a \code{glinv_gauss} model contains trait values at their tips.
#'
#' Returns true if and only if the \code{glinv_gauss} model were initialised with \code{X=NULL}
#' and the user had never called \code{set_tips} on it.
#'
#' @param mod A \code{glinv_gauss} object.
#' @export
has_tipvals = function (mod) UseMethod('has_tipvals')

#' @rdname has_tipvals
#' @method has_tipvals glinv_gauss
#' @export
has_tipvals.glinv_gauss = function (mod) .Call(Rxavail, mod$ctree)


#' Set trait values at the tip for a \code{glinv_gauss} model.
#'
#' If a \code{glinv_gauss} or \code{glinv} object were initalised with \code{X=NULL}, methods like
#' \code{lik} will not work because it lacks actual data. In this case, the user
#' should set the trait values using this method. Any existing trait values
#' were replaced after calling this.
#'
#' This method alters an underlying C structure, therefore has a mutable-object
#' semantic (See example).
#'
#' @param mod A \code{glinv_gauss} or \code{glinv} object.
#' @param X   A matrix of trait values, in which \code{X[p,n]} stores the p-th dimension
#'            of the multivariate trait of the n-th tip of the phylogeny. 
#' @examples
#' tr = ape::rtree(10)
#' model = glinv_gauss(tr, x0=c(0,0))  # The \code{X} argument is implicitly NULL
#' traits = matrix(rnorm(20), 2, 10)
#' set_tips(model, traits)
#' # Now lik(model, ...) should work
#' 
#' @export
set_tips = function (mod, X) UseMethod('set_tips')

#' @rdname set_tips
#' @method set_tips glinv_gauss
#' @export
set_tips.glinv_gauss = function (mod, X)  {
  .Call(Rsettip, mod$ctree, clean_x(X, mod$dimtab, mod$apetree))
  invisible()
}

#' @rdname set_tips
#' @method set_tips glinv
#' @export
set_tips.glinv       = function (mod, X)    set_tips(mod$rawmod, tip_purge(X))


#' Simulate random trait values from models.
#'
#' Simulate random trait values from the Gaussian branching process specified by \code{mod}. 
#' 
#'
#' @param mod    Either a \code{glinv_gauss} or \code{glinv} object.
#' @param par    Parameters underlying the simulation, in the same format as \code{lik.glinv_gauss} or \code{lik.glinv}.
#' @param Nsamp  Number of sample point to simulate.
#' @return A list containing Nsamp elements, each of which is a length-n list, where \eqn{n} is the number of tips
#'         containing the simulated trait values. 
#' @export
rglinv = function (mod, par, Nsamp) UseMethod('rglinv')

#' @rdname rglinv
#' @export
rglinv.glinv         = function (mod, par, Nsamp) {
  ## TODO: either warn the user when they have missing values, or recover the original shape.
  rglinv.glinv_gauss(mod$rawmod, mod$gaussparams_fn(par), Nsamp)
}

#' @rdname rglinv
#' @export
rglinv.glinv_gauss = function (mod, par, Nsamp)
  .Call(Rvwphi_simul,
        mod$ctree,
        as.integer(length(mod$apetree$tip.label)),
        as.integer(mod$dimtab),
        as.double(par),
        as.integer(Nsamp),
        as.double(mod$x0))

clone_topology = function (ctree) .Call(R_clone_tree, ctree)


#' Construct an object representing a GLInv model with respect to the underlying Gaussian process parameters.
#'
#' The \code{glinv_gauss} function constructs an object of class \code{glinv_gauss}, which represents a lower-level 
#' GLInv model with respect to the underlying Gaussian process space. The likelihood Hessian of, for example, Brownian motion
#' and Ornstein-Uhlenbeck models can be computed by applying the calculus chain rule to the output of Jacobians and Hessians
#' from this class.
#' 
#' The \code{glinv_gauss} class does not include any information for dealing with evolutionary regimes, lost traits, and
#' missing data, nor does it facilitate reparametrisation. These are all functionality of the \code{glinv} class instead.
#' The member variables of the objects of the \code{glinv_gauss} class only are for the users' convenience to \emph{read}
#' the information about the model, and the user \emph{should not} modify its member variables directly.
#' 
#' For each non-root node \eqn{i} in the phylogeny, the multivariate trait vector \eqn{x_i} follows
#' a Gaussian distribution with mean \eqn{\Phi_i x_j + w_i} and variance \eqn{V_i} when conditioned on
#' the mother's trait vector \eqn{x_j}. The `parameters' of this model is, therefore, the joint of all
#' \eqn{(\Phi_i, w_i V'_i)} for all nodes \eqn{i}. The root does not have any associated parameters.
#'
#'
#' @param tree   A tree of class \code{ape::phylo}.
#' @param x0     A vector representing the root's trait vector.
#' @param dimtab An integer, a vector of integers, or NULL, specifying the number of dimensions of each nodes of the tree.
#'               If it is a vector, \code{dimtab[n]} is the trait vector dimension of node \code{n}. If it is only a single
#'               integer than all nodes are assumed to have the same amount of dimensions. If it is NULL then all
#'               nodes are asummed to have the same amount of dimensions as \code{x0}.
#' @param X      Trait values, either a matrix in which \code{X[p,n]} stores the \eqn{p}-th dimension
#'               of the multivariate trait of the \eqn{n}-th tip of the phylogeny, or a list in which
#'               \code{X[[n]]} is a numeric vector representing the multivariate trait of the \eqn{n}-th tip.
#'               The latter form is required if not all the tips has the same number of dimensions.
#' @return       An object of S3 class \code{glinv_gauss} with the following components
#'               \describe{
#'                 \item{ctree}{A pointer to an internal C structure.}
#'                 \item{apetree}{Identical to the \code{tree} argument.}
#'                 \item{x0}{The trait vector at the root of the tree.}
#'                 \item{dimtab}{Identical to the \code{dimtab} argument.}
#'                 \item{gaussdim}{The number of dimension of the parameter space of this model.}
#'               }
#' @export
glinv_gauss = function (tree, x0, dimtab=NULL, X=NULL) {
  if (!'phylo' %in% class(tree))                stop('The tree must be an ape tree')
  if (!ape::is.rooted(tree))                    stop('The input phylogenetic tree must be rooted')
  if (!is.numeric(x0))                          stop("x0 must be numeric")
  if (! (all(!is.na(x0)) && all(!is.nan(x0))))  stop("x0 must not contain NA or NaN")
  if (is.null(dimtab)) {
    dimtab = as.integer(rep(length(x0), length(tree$edge)/2+1))
  } else {
    if (!(length(dimtab) == length(tree$edge)/2+1 ||
          length(dimtab) == 1))                 stop("dimtab must be an integer vector of length N_NODES or 1")
    if (length(dimtab) == 1)                    dimtab = as.integer(rep(dimtab, length(tree$edge)/2+1))
  }
  if (length(x0) != dimtab[tree$edge[1,1]])     stop("x0's dimension doesn't match dimtab's")
  mode(x0)     = 'double'
  mode(dimtab) = 'integer'
  ctree = .Call(Rnewnode,
                { ed = t(tree$edge); mode(ed) = 'integer'; ed },
                clean_x(X, dimtab, tree),
                dimtab)
  o = list(ctree = ctree, apetree = tree, x0 = x0, dimtab = dimtab,
           gaussdim = as.integer(sum(dimtab[1:(tree$edge[1,1]-1)])))
  class(o) = 'glinv_gauss'
  o$nparams = nparams(o)
  o
}

clean_x = function (X, dimtab, tree) {
  X_clean = if      (is.list(X))      lapply(X, function (x) {mode(x)='double'; x})
            else if (is.numeric(X))   plyr::alply(matrix(data=X, ncol=length(tree$tip.label)),
                                                  2, function (y) matrix(as.double(y)))
            else if (is.matrix(X))    plyr::alply(X, 2, function (y) matrix(as.double(y)))
            else if (is.null(X))      return(NULL)
  if (length(X_clean) != length(tree$tip.label))
    stop(sprintf("X contains %d tips but the tree has %d", length(X_clean), length(tree$tip.label)))
  for (i in seq_along(X_clean))
    if (!is.numeric(X_clean[[i]]))  stop(sprintf("The %d-th tip (name: %s) isn't numeric", i, tree$tip.label[i]))
    if (!(all(!is.na(X_clean[[i]])) && all(!is.nan(X_clean[[i]]))))
        stop(sprintf('The %d-th tip (name: %s) contains NA/NaN', i, tree$tip.label[i]))
    if (length(X_clean[[i]]) != dimtab[i])
      stop(sprintf("The %d-th tip (name: %s) is %d-dimensional in `X` but `dimtab` says it's %d-dimensional",
                   i, tree$tip.label[i], length(X_clean), dimtab[i]))
  X_clean
}

#' Compute the likelihood of a GLInv model
#'
#' This is a S3 generic method. For the \code{\link{glinv}} class, which is a high-level user interface, please 
#' see \code{\link{lik.glinv}}; and for \code{\link{glinv_gauss}}, which is a lower-level facility, please see
#' \code{\link{lik.glinv_gauss}}.
#'
#' @param    mod    An object of either \code{\link{glinv}} or \code{\link{glinv_gauss}} class.
#' @param    ...    Further arguments to be passed to the S3 methods.
#' @export
lik = function (mod, ...) UseMethod('lik')

#' Compute the likelihood, gradient, and Hessian of a full GLInv model in the underlying Gaussian parameters
#'
#' The \code{lik.glinv_gauss} function computes the likelihood of a full \code{glinv_gauss} model.
#'
#' The parameter vector \code{par} should be the concatenation of all \eqn{(\Phi_i, w_i, V'_i)} in accending
#' order sorted by \eqn{i}, the node number (which is the same node numbers as in \code{tree$edge}). The matrix
#' \eqn{\Phi_i} is flattened in column-major order and \eqn{V'_i} is the lower-triangular part of V_i,
#' column-major-flattened. Since the root does not have parameters, its entry is simply skipped.
#' For example, if a binary tree has 10 non-root nodes in total and each of them are 3 dimensional, then
#' each \eqn{(\Phi_i, w_i, V'_i)} should have \eqn{9+3+6=18} elements; thus after concatenation \code{par} should
#' be a 180 elements.
#' 
#' @rdname glinv_gauss
#' @param mod A model object of class \code{glinv_gauss}.
#' @param par A vector, containing the parameters at which the likelihood should be computed.
#' @param ... Not used.
#' @examples
#' tr = ape::rtree(3)
#' model = glinv_gauss(tr, x0=c(0,0), X=matrix(rnorm(6),2,3))
#' par = unlist(
#'  list(
#'    list('Phi' = c(1,0,0,1), # Parameters for node #1, a tip
#'         'w'   = c(-1,1),
#'         'V'   = c(1,0,1)),  # Lower triangular part of a 2D identity matrix
#'    list('Phi' = c(2,0,0,2), # For node #2, a tip
#'         'w'   = c(-2,2),
#'         'V'   = c(2,0,2)),
#'    list('Phi' = c(3,0,0,3), # For node #3, a tip
#'         'w'   = c(-3,3),
#'         'V'   = c(3,0,3)),
#'    list('Phi' = c(4,0,0,4), # For node #5. Node #4 skipped as it is a root
#'         'w'   = c(-4,4),
#'         'V'   = c(4,0,4))
#'    ))
#' print(par)
#' lik(model, par)
#' grad(model, par)
#' hess(model, par)
#' @export
lik.glinv_gauss = function (mod, par=NULL, ...) {
  if (is.null(par)) {
    mod
    function (par) {
      .Call(Rndphylik, mod$ctree, par, mod$x0, mod$gaussdim)
    }
  } else {
    .Call(Rndphylik, mod$ctree, par, mod$x0, mod$gaussdim)
  }
}

#' Compute the log-likelihood gradients of GLInv models
#'
#' For the \code{\link{glinv}} class, which is a high-level user interface, please see \code{\link{grad.glinv}}; and
#' for \code{\link{glinv_gauss}}, which is a lower-level facility, please see \code{\link{grad.glinv_gauss}}.
#'
#' @param    mod    An object of either \code{\link{glinv}} or \code{\link{glinv_gauss}} class.
#' @param    ...    Further arguments to be passed to the S3 methods.
#' @export
grad = function (mod, ...) UseMethod('grad')

#' Compute the log-likelihood gradients of a full GLInv models in the underlying Gaussian parameters
#'
#' The \code{grad.glinv_gauss} function computes the log-likelihood gradients of a \code{glinv_gauss} models.
#' If \code{par} is NULL, it returns a function that, when called, returns the same thing as if \code{grad.glinv_gauss}
#' were called with \code{par} argument.
#' 
#' @rdname glinv_gauss
#' @param lik    If \code{TRUE}, \code{grad.glinv_gauss} and \code{hess.glinv_gauss} returns also the log-likelihood.
#' @param ...    Not used.
#' @export
grad.glinv_gauss = function (mod, par=NULL, lik=F, ...) {
  if (is.null(par)) {
    function (par) grad_glinv_gauss_(mod=mod, par=par, lik=lik)
  } else {
    grad_glinv_gauss_(mod=mod, par=par, lik=lik)
  }
}

grad_glinv_gauss_ = function (mod, par, lik=F, ...) {
  l = .Call(Rdphylik, mod$ctree, par, mod$x0, mod$gaussdim)
  g = c(.Call(Rextractderivvec, mod$ctree))
  if (!lik) g
  else list(lik = l, grad = g)
}

#' Compute the log-likelihood Hessian of GLInv models
#'
#' For the \code{\link{glinv}} class, which is a high-level user interface, please see \code{\link{hess.glinv}}; and
#' for \code{\link{glinv_gauss}}, which is a lower-level facility, please see \code{\link{hess.glinv_gauss}}.
#' 
#' @param    mod    An object of either \code{\link{glinv}} or \code{\link{glinv_gauss}} class.
#' @param    ...    Further arguments to be passed to the S3 methods.
#' @export
hess = function (mod, ...) UseMethod('hess')

#' Compute the log-likelihood Hessian of full  GLInv models in the underlying Gaussian parameters
#'
#' The \code{hess.glinv_gauss} function computes the log-likelihood Hessian of a \code{glinv_gauss} models.
#' 
#' @rdname glinv_gauss
#' @param grad        If \code{TRUE}, \code{hess.glinv_gauss} returns also the gradient.
#' @param directions  Either \code{NULL} or a matrix with \code{mod$nparams} many rows and arbitrarily many columns.
#'                    If \code{NULL}, `hess.glinv_gauss` returns the Hessian matrix itself, which is typically
#'                    a huge matrix; otherwise, the funciton returns a square matrix \eqn{M} such that \eqn{M_ij}
#'                    contains \eqn{d_i^T H d_j}, where \eqn{d_i} is the \eqn{i}-th column of \code{directions} and
#'                    \eqn{H} is the huge Hessian matrix, without storing \eqn{H} itself in memory.
#' @param ...         Not used.
#' @export
hess.glinv_gauss = function (mod, par=NULL, lik=F, grad=F, directions=NULL, ...) {
  mod
  if (is.null(par)) {
    function (par) {
      hess_glinv_gauss_(mod, par, lik, grad, directions, ...)
    }
  } else {
    hess_glinv_gauss_(mod, par, lik, grad, directions, ...)
  }
}

hess_glinv_gauss_ = function (mod, par, lik=F, grad=F, directions=NULL, ...) {
  if (is.null(directions)) {
    l = .Call(Rhphylik, mod$ctree, par, mod$x0, mod$gaussdim)
    h = .Call(Rextracthessall, mod$ctree)
  } else {
    h = .Call(Rhphylik_dir, mod$ctree, par, mod$x0, mod$gaussdim, directions)
    l = attr(h, 'likelihood')
    attr(h, 'likelihood') = NULL
  }
  if (!lik && !grad) h
  else {
    ans = list()
    if (lik)  ans[['lik']]  = l
    if (grad) ans[['grad']] = c(.Call(Rextractderivvec, mod$ctree))
    ans[['hess']]           = h
    ans
  }
}

#' @rdname glinv_gauss
#' @param  x    An object of class \code{glinv_gauss}.
#' @param  ...  Not used.
#' @export
print.glinv_gauss = function (x, ...) {
  y = x$apetree
  cat(paste("A `raw' GLInv model, with respect to the underlying Gaussian parameters. It has", length(y$tip.label), "tips,",
            y$Nnode, "internal nodes, and", nparams(x), "parameters. Tip values are",
            if (has_tipvals(x)) "already set." else "empty (meaning `lik()` etc. won't work).","\n"))
}


tip_purge = function (x, ...) UseMethod('tip_purge')
tip_purge.list = function (X) {
  lapply(seq_len(length(X)), function(i) {
    y = X[[i]]
    y[(!is.na(y))&(!is.nan(y))]})
}
tip_purge.matrix = function (X)
  lapply(seq_len(ncol(X)), function(i) {
    y = X[,i]
    y[(!is.na(y))&(!is.nan(y))]})


#' Construct an GLInv model with respect to user-specified parametrisation
#'
#' The \code{glinv} function construct an object of class \code{glinv}, which represents a GLInv model with respect
#' to a user-specified parametrisation.
#'
#' Models for \code{glinv} assumes one or more evolutionary regimes exists in the phylogeny. The \code{regimes} parameters defines
#' how many regimes there are, where do the regimes start, and what parameterisation function it has. If \code{regimes} were
#' NULL then a single regime starting from the root node is assumed. Multiple regimes could share the same parametrisation
#' function (and thus the same parameters) by specifying the same index; therefore the number of regimes may differs from 
#' the number of parametrisation functions. One and only one regime must start from the root of the phylogeny.
#'
#' If \code{X} contains \code{NA} in the \eqn{p}-th dimension of the \eqn{i}-th tip (whose node ID is also \eqn{i}) then \eqn{X_pi} is
#' tagged \code{MISSING}. No other tags of any other nodes are changed. The \eqn{p}-th dimension of any node \eqn{j}, regardless of
#' whether or not it is an internal node or a tips, is tagged \code{LOST} if and only if the \eqn{p}-th dimension of all tips inside
#' the clade started at \eqn{j} are \code{NaN}. Any entry that is neither \code{LOST} nor \code{MISSING} are tagged \code{OK}. These
#' tags are then passed into the user-defined functions \code{parfns} etc. as arguments; therefore the user is free to specify how
#' these tags are handled. \code{x0} cannot contain missing values, and the vectors of missingness tags passed to \code{parfns}, for
#' any nodes, are always of the same length as \code{x0}.
#'
#' Before this package calls the functions in \code{parhess}, it adds, into the function's environment, a variable named \code{INFO__} 
#' which contains some extra information.
#' 
#' Passing a single function to \code{parfns} is equivalent to passing a singleton list; and the same is true for \code{parjacs},
#' \code{parhess}, and \code{pardims}.
#'
#' @param tree    A tree of class \code{\link[ape]{phylo}}.
#' @param x0      A vector representing the root's trait vector. Must not contain \code{NA} and \code{NaN}.
#' @param X       A matrix of trait values, in which \code{X[p,n]} stores the p-th dimension
#'                of the multivariate trait of the n-th tip of the phylogeny. \code{NA} and \code{NaN}
#'                has special meanings (See Details).
#' @param parfns  A list of functions that maps from the user-parametrisation to the underlying Gaussian parameters.
#'                Each of them returns a vector of concatenated \eqn{(\Phi, w, V')}, where \eqn{V'} is the lower triangular
#'                part of \eqn{V}, and accepts four arguments: a vector of parameters whose length is specified
#'                by the \code{pardims} argument to the \code{glinv_gauss} function, the branch length leading to the currently processing node, 
#'                a vector of factors with three levels indicating which dimensions are missing or lost in the mother of
#'                the current node, and a vector of factors with the same three levels indicating missingness of the current
#'                node.
#' @param pardims A vector of integers, which has the same amount elements as the length of parfns.
#'                \code{pardims[i]} indicates the length of the parameter vector that \code{parfns[i]} accepts.
#' @param regimes A list of length-two integer vectors. Each of these length-two vectors specifies an evolutionary regime
#'                and consists of a named element \code{start}, which specifies the node ID at which an evolutionary regime 
#'                starts, and another named element \code{fn}, which is an index of \code{parfns}, indicating which parametrisation
#'                function this evolutionary regime should use.
#' @param parjacs A list of functions, which has the same amount elements as that of \code{parfns}.
#'                \code{parjacs[i]} accepts the same arguments as \code{parfns[i]} and returns the Jacobian of \code{parfns[i]}.
#' @param parhess A list of functions, which has the same amount elements as that of \code{parfn[i]}.
#'                \code{parhess[i]} accepts the same arguments as \code{parfns[i]} and returns a list of three 3D arrays,
#'                named \code{Phi}, \code{w}, \code{V} respectively inside the list. \code{((parhess[[i]])(...))$Phi[m,i,j]} contains
#'                the cross second-order partial derivative of \eqn{\Phi_m} (here we treat the matrix \eqn{\Phi} as a
#'                column-major-flattened vector) with respect to the \eqn{i}-th and\eqn{j}-th user parameters;
#'                while \code{((parhess[[i]])(...))$w[m,i,j]} and \code{((parhess[[i]])(...))$V[m,i,j]}
#'                analogously contains second-order derivative of \eqn{w_m} and \eqn{V'_m}.
#' @param x       A \code{glinv} object, for the \code{print.glinv} S3 method. 
#' @return        The \code{glinv} function returns a model object of S3 class \code{glinv}. Elements are:
#'                
#'                  \item{rawmod}{An object of class \code{glinv_gauss}.}
#'                  \item{regimes}{Identical to the \code{regimes} argument.}
#'                  \item{regtags}{An integer vector of the same length as the number of nodes. The \eqn{i}-th element is
#'                                 the regime ID (corresponding to the index in the \code{regimes} argument to the \code{glinv_gauss} function) of
#'                                 node \eqn{i}. \code{NA} at the root.}
#'                  \item{misstags}{A factor matrix with three ordered levels, \code{LOST}, \code{OK}, and \code{MISSING}. Each column
#'                                  corresponds to a node and row to a trait dimension.}
#'                  \item{nparams}{The sum of the \code{pardims} argument, an integer.}
#'                  \item{pardims}{Identical to the \code{pardims} arguemnt.}
#'                  \item{parfntags}{An integer vector of the same length as the number of nodes. The \eqn{i}-th element is
#'                                   the index of \code{parfns} that corresponds to node \eqn{i}. \code{NA} at the root.}
#'                  \item{parfns}{Identical to the \code{parfns} argument.}
#'                  \item{parjacs}{Identical to the \code{parjacs} argument.}
#'                  \item{parhess}{Identical to the \code{parhess} argument.}
#'                  \item{parsegments}{A \eqn{K}-by-2 matrix of integer indicies, where \eqn{K} is the length of \code{parfns}.
#'                                     If \code{v} is a vector that \code{\link{lik.glinv}} accepts, then
#'                                     \code{v[parsegments[k,1]:parsegments[k,2]]} is the parameter vector should \code{parfns[[k]]}
#'                                     accept.}
#'                  \item{gausssegments}{A \eqn{N}-by-2 matrix of integer indicies, where \eqn{N} is the number of nodes.
#'                                       If \code{w} is a vector that \code{\link{lik.glinv_gauss}} accepts, then
#'                                       \code{w[gausssegments[i,1]:gausssegments[i,2]]} is the concatenated \eqn{(\Phi, w, V')}
#'                                       corresponding to node \eqn{i}.}
#'                  \item{gaussparams_fn}{A function that accepts a parameter vector of length \code{nparams} and returns a
#'                                        parameter vector of length \code{rawmod$nparams}. When called, this function
#'                                        traverses the tree, calls the functions in parfns on each node, and assemble 
#'                                        the results into a format that \code{\link{lik.glinv_gauss}} accepts.}
#'                  \item{gaussparams_jac}{A function that accepts a parameter vector of length \code{nparams} and returns a
#'                                        \eqn{p}-by-\eqn{q} Jacobian matrix, where \eqn{p} is \code{rawmod$nparams} and \eqn{q}
#'                                        is \code{nparams} in this object. When called, this function traverses
#'                                        the tree, calls the functions in \code{parjacs} on each node, and row-concatenates the
#'                                        result in an order consistent with what \code{\link{lik.glinv_gauss}} accepts.}
#' @examples
#' \dontrun{
#' ### --- STEP 1: Making an example tree and trait data
#' ntips = 100
#' k     = 2                 # No. of trait dimensions
#' tr    = ape::rtree(ntips) 
#' X     = matrix(rnorm(k*ntips), k, ntips)
#' x0    = rnorm(k)
#' 
#' ### --- STEP 2: Making a model object. We use OU as an example.
#' mod = glinv(tr, x0, X,
#'             parfns = list(ou_haltlost(oupar)),
#'             pardims = list(nparams_ou(k)),
#'             parjacs = list(dou_haltlost(oujac)),
#'             parhess = list(hou_haltlost(ouhess)))
#'
#' ### --- STEP 3: Try getting the likelihood, gradient etc.
#' H     = matrix(c(1,0,0,-1), k)
#' theta = c(0,0)
#' sig   = matrix(c(0.5,0,0,0.5), k)
#' sig_x = t(chol(sig))
#' par_init = c(H=H,theta=theta,sig_x=sig_x[lower.tri(sig_x,diag=T)])
#' print(par_init)
#' print(lik(mod)(par_init))
#' print(grad(mod)(par_init))
#' print(hess(mod)(par_init))
#'
#' ### --- STEP 4: Fitting a model
#' fitted = fit(mod, par_init, method='L-BFGS-B')
#' print(fitted)
#' 
#' ### --- STEP 5: Estimating variance-covariance of the MLE
#' v_estimate = varest(mod, fitted)
#'
#' ### --- STEP 6: Get marginal confidence intervals
#' print(marginal_ci(v_estimate, lvl=0.95)) 
#' }
#' @export
glinv = function (tree, x0, X, parfns, pardims, regimes=NULL, parjacs=NULL, parhess=NULL) {
  if (!'phylo' %in% class(tree)) stop("`tree` must be of class ape::phylo")
  if (!ape::is.rooted(tree))     stop("Only rooted trees are supported.")
  ## In case there is only a single regime fns, jacs and hess is allowed to be a single function instead
  ## of a list of functions. In this case the only evolutionary regime is set to start at root automatically.
  if (is.function(parfns) || (is.list(parfns) && length(parfns) == 1 && is.function(parfns[[1]]))) {
    parfns   = c(parfns)
    parjacs  = c(parjacs)
    parhess  = c(parhess)
    if (!is.null(regimes))
      stop("`regimes` cannot and need not be specified when `parfns` contains only one function")
    regimes = list(c(start=tree$edge[1,1], fn=1))
  }
  pardims  = as.integer(unlist(pardims))
  modtpl   = glinv_gauss(tree, 0, 1, X=NULL)
  if (! is.null(X)) {
    if (is.list(X)) {
      X = unlist(X)
      dim(X) = c(length(X)/length(tree$tip.label),length(tree$tip.label))
      traitnames = names(x0)
    } else if (is.double(X)) {
      if (ncol(X) != length(tree$tip.label))
        stop('If `X` is a matrix, it must have the same number of *columns* as the number of tips.')
      traitnames = if      (!is.null(rownames(X))) rownames(X)
                   else if (!is.null(names(x0)))   names(x0)
                   else NULL
    }
    misstags = tag_missing(modtpl, X)
    dimtab   = {
      OKlvl = which(levels(misstags)=='OK')
      apply(misstags,2,function(y) sum(y == OKlvl))
    }
    if (any(dimtab == 0))
      stop(sprintf("All dimensions at Node #%d are either lost or missing",
                   which(dimtab == 0)[1]))
  } else {
    nnodes = length(unique(c(tree$edge)))
    k      = length(x0)
    misstags = factor(rep(1,k*nnodes), levels=c(0,1,2), labels=c('LOST','OK','MISSING'))
    dim(misstags) = c(k, nnodes)
    dimtab = rep(k, nnodes)
  }
  rownames(misstags) = traitnames
  if (!length(x0) == dimtab[tree$edge[1,1]]) 
    stop(sprintf("The dimension of `x0` must be %d but I got %d",
                 dimtab[tree$edge[1,1]], length(x0)))
  rawmod     = glinv_gauss(tree, x0, dimtab, X = if (is.null(X)) NULL else tip_purge(X))
  regtags    = tag_regimes(modtpl, unlist(Map(function(x) x['start'], regimes)))
  parfntags  = tag_parfns(regtags, regimes)
  parsegments= matrix(c(1L,1L+cumsum(pardims[-length(pardims)]),cumsum(pardims)),ncol=2)
  colnames(parsegments) = c('start','end')
  gausssegments = .Call(Rvwphi_paradr, rawmod$ctree)
  colnames(gausssegments) = c('start','end')
  obj    = list(rawmod        =rawmod,
                regimes       =regimes,
                regtags       =regtags,
                misstags      =misstags,
                dimtab        =dimtab,
                nparams       =sum(pardims),
                pardims       =pardims,
                parfntags     =parfntags,
                parfns        =parfns,
                parjacs       =parjacs,
                parhess       =parhess,
                parsegments   =parsegments,
                gausssegments =gausssegments)
  obj[['gaussparams_fn']]     =gaussparams(obj)
  obj[['gaussparams_jac']]    =if (is.null(parjacs)) NULL
                               else                  gaussparams_grad(obj)
  ## We don't need gaussparams_hess but perhaps add it for consistency...?
  class(obj) = 'glinv'
  obj
}

#' @rdname glinv
#' @param ...  Not used.
#' @export
print.glinv = function (x, ...) {
  mod = x
  x = mod$apetree
  cat(sprintf(paste0(
    'A GLInv model with %d regimes and %d parameters in total, %s.\n',
    'The phylogeny has %d tips and %d internal nodes.\n'),
    length(unique(mod$regtags)), mod$nparams,
    {
      if (length(mod$parsegments) > 2) {
        i = 1 # currently processing parfn ID.
        paste0(
          c(paste0(c('among which',
                     apply(mod$parsegments, 1, function (x) {
                       s = sprintf('    the %d~%d-th parameters are asociated with regime no. {%s}',
                                   x['start'], x['end'],
                                   paste0(
                                   {
                                     associated_regimes = c()
                                     for (j in seq_along(mod$regimes)) {
                                       if (mod$regimes[[j]]['fn'] == i)
                                         associated_regimes = c(associated_regimes, j)
                                     }
                                     associated_regimes
                                   }, collapse=','))
                       i <<- i+1
                       s
                     })), collapse=';\n'),
            ',\nwhere \n',
            paste0({
              j = 1
              lapply(mod$regimes, function (r) {
                s=sprintf('    regime #%d starts from node #%d%s',
                        j,
                        r['start'],
                        if (r['start'] == mod$rawmod$apetree$edge[1,1]) ', which is the root' else '')
                j <<- j+1
                s
              })},
              collapse=';\n')
            ), collapse='')
      } else {
        "all of which are associated to the only one existing regime, which starts from the root"
      }
    },
    length(mod$rawmod$apetree$tip.label), mod$rawmod$apetree$Nnode))
}



#' Compute the log-likelihood of a GLInv model
#'
#' The \code{lik.glinv} function returns a function which accepts a parameter vector, which is of length \code{mod$nparams},
#' and returns the log-likelihood.
#'
#' @rdname glinv
#' @param mod   An object of class \code{glinv}.
#' @param ...   Not used.
#' @export
lik.glinv = function (mod, ...) {
  function (x) {
    if (length(x) != mod$nparams)
      stop(sprintf("lik.glinv: your model should have %d parameters but I got %d",
                   mod$nparams, length(x)))
    lik(mod$rawmod, mod$gaussparams_fn(x))
  }
}

#' Compute the gradient of log-likelihood of a GLInv model
#'
#' The \code{grad.glinv} function returns a function which accepts a parameter vector, which is of length \code{mod$nparams},
#' and returns the gradient of log-likelihood with respect to this parametrisation.
#'
#' @rdname glinv
#' @param numDerivArgs    Arguments to pass to \code{numDeriv::\link[numDeriv]{jacobian}}. Only used the user did not specify the
#'                        \code{parjacs} arguments when creating \code{mod}.
#' @param ...             Not used.
#' @export
grad.glinv = function (mod, numDerivArgs = list(method='Richardson', method.args=list(d=.5,r=3)), ...) {
  if (is.null(mod$gaussparams_jac))
    function(x) c(do.call(numDeriv::jacobian, c(list(lik(mod), x), numDerivArgs)))
  else function (x)
    if (length(x) != mod$nparams)
      stop(sprintf("lik.glinv: your model should have %d parameters but I got %d",
                   mod$nparams, length(x)))
    else {
      c(grad(mod$rawmod, mod$gaussparams_fn(x)) %*% mod$gaussparams_jac(x))
    }
}

#' Compute the Hessian of log-likelihood of a GLInv model
#'
#' The \code{hess.glinv} function returns a function which accepts a parameter vector, which is of length \code{mod$nparams},
#' and returns the Hessian matrix of log-likelihood with respect to this parametrisation.
#'
#' @rdname glinv
#' @param numDerivArgs    Arguments to pass to \code{numDeriv::\link[numDeriv]{jacobian}}. Only used the user did not specify the
#'                        \code{parjacs} arguments when creating \code{mod}.
#' @param store_gaussian_hessian If \code{TRUE} and \code{method} is not \code{mc}, the returned list will contain
#'                             a (usually huge) Hessian matrix \code{gaussian_hessian} with respect to the Gaussian
#'                             parameters \eqn{\Phi, w, V'}. This option significantly increases the amount of memory
#'                             the function uses, in order to store the matrix.
#' @param ...             Not used.
#' @export
hess.glinv = function (mod,
                       numDerivArgs = list(method='Richardson', method.args=list(d=.5,r=3)),
                       store_gaussian_hessian = F, ...) {
  mod
  function (par) {
    if (length(par) != mod$nparams)
      stop(sprintf("lik.glinv: your model should have %d parameters but I got %d",
                   mod$nparams, length(par)))
    jac = if (is.null(mod$parjacs))
            function (x) c(do.call(numDeriv::jacobian, c(list(mod$gaussparams_fn, x), numDerivArgs)))
          else
            function (x) mod$gaussparams_jac(x)

    if (is.null(mod$parhess))
      stop("Second-order-derivatives of the user's parameterisation is not supplied")
    ## TODO: Check the type of parhess and warn the user when constructing the object!
    if (store_gaussian_hessian) {
      Hvwphi = hess(mod$rawmod, as.double(mod$gaussparams_fn(par)))
      J = jac(par)
      left = crossprod(J, Hvwphi)
      lin = left %*% J
    } else {
      J = jac(par)
      lin = hess(mod$rawmod, as.double(mod$gaussparams_fn(par)), directions = J)
    }
    ## This will be what the C function `Rcurvifyhess' will be calling
    curvifier = function (nid, par) {
      eid    = which((tr <- mod$rawmod$apetree)$edge[,2] == nid)
      if (length(eid) != 1L)
        stop(sprintf(paste0("varest(): node ID not found or it has multiple mothers. Was the ape ",
                            "tree changed after initialising ",
                            "the `glinv` object? NodeID = %d"), nid))
      pid    = tr[['edge']][eid,1]
      fid    = mod$parfntags[nid]
      hessfn = mod$parhess[[fid]]
      kp     = mod$misstags[,pid]
      kn     = mod$misstags[,nid]
      d      = mod$pardims[fid]
      pstart = mod$parsegments[fid,'start']
      pend   = mod$parsegments[fid,'end']
      el = mod[['rawmod']][['apetree']][['edge.length']][eid]
      environment(hessfn)[['INFO__']] = list(reparametrisation_jacobian = J,
                                             mod                        = mod,
                                             node_id                    = nid,
                                             parent_id                  = pid,
                                             parfn_id                   = fid)
      in_regime_res = hessfn(par[pstart:pend], el, kp, kn)
      .Call(Rchkusrhess, in_regime_res, mod$nparams, mod$pardims[fid],
            nid, pid, mod$dimtab[nid], mod$dimtab[pid])
      ## Hessfn is regime specific so we need to augment the hessians by zeros for other parameters.
      r = lapply(in_regime_res, function (r) {
        thisdim = c(dim(r))
        thisdim[c(2,3)] = mod$nparams
        fullarr = array(0., thisdim)
        fullarr[,pstart:pend, pstart:pend] = c(r)
        fullarr
      })
      r
    }
    ## lin is directly modified in the C code after this.
    .Call(Rcurvifyhess, lin, par, mod$rawmod$ctree, curvifier, environment())
    for (i in seq_along(mod$parhess))
      rm('INFO__', pos=environment(mod$parhess[[i]]))
    if (store_gaussian_hessian) {
      attr(lin, 'gaussian_hessian') = Hvwphi
    }
    lin
  }
}

gaussparams = function (mod) {
  tr    = mod$rawmod$apetree
  p     = nparams(mod$rawmod)
  function (par) {
    res   = matrix(0, nrow=p, ncol=1)
    for (j in seq_len(nrow(tr$edge))) {
      nid    = tr$edge[j,2]
      pid    = tr$edge[j,1]
      fid    = mod$parfntags[nid]
      el     = tr$edge.length[j]
      fn     = mod$parfns[[fid]]
      kp     = mod$misstags[,pid]
      kn     = mod$misstags[,nid]
      pstart = mod$parsegments[fid,'start']
      pend   = mod$parsegments[fid,'end']
      gstart = mod$gausssegments[nid,'start']
      gend   = mod$gausssegments[nid,'end']
      res[gstart:gend,] = fn(par[pstart:pend], el, kp, kn)
    }
    res
  }
}
gaussparams_grad = function (mod) {
  tr    = mod$rawmod$apetree
  p     = mod$rawmod$nparams
  blksiz= 
  function (par) {
    res   = matrix(0, nrow=p, ncol=mod$nparams)
    for (j in seq_len(nrow(tr$edge))) {
      nid    = tr$edge[j,2]
      pid    = tr$edge[j,1]
      fid    = mod$parfntags[nid]
      el     = tr$edge.length[j]
      jacfn  = mod$parjacs[[fid]]
      kp     = mod$misstags[,pid]
      kn     = mod$misstags[,nid]
      pstart = mod$parsegments[fid,'start']
      pend   = mod$parsegments[fid,'end']
      gstart = mod$gausssegments[nid,'start']
      gend   = mod$gausssegments[nid,'end']
      res[gstart:gend,pstart:pend] = jacfn(par[pstart:pend], el, kp, kn)
    }
    res
  }
}


default_ifwrong = function(f,n,default)
  function(x) tryCatch(f(x), error=function(...) rep(default,n))

## For non-optim routines
fit_fnwrap = function (f, fnscale) function (...) f(...)/fnscale
fit_grwrap = function (g, fnscale) function (...) g(...)/fnscale

#' @rdname fit.glinv
#' @export
fit = function (mod, ...) UseMethod('fit')

#' Fitting a GLInv model via numerical optimisation
#' 
#' \code{fit.glinv} finds the maximum likelihood estimate of a \code{glinv} model by solving a numerical
#' optimisation problem.
#' 
#' If \code{method} is \code{L-BFGS-B}, then \code{\link[lbfgsb3c]{lbfgsb3c}} is used for optimisation;
#' if it is \code{CG} then \code{\link[Rcgmin]{Rcgmin}} is used; if it is \code{BB} then
#' \code{\link[BB]{BBoptim}} is used, otherwise the method argument is passed to \code{\link{optim}}.
#'
#' By default, \code{L-BFGS-B} declares convergence when the change of function value is small, \code{CG}
#' tests stops when change of gradient squared-Euclidean-norm is small, \code{BB} stops when either the
#' change of function values, or the infinity norm of a project gradient, is small. These can be changed
#' through the \code{control} argument and the user should refer to the optimisation packages' respective
#' documentation for details.
#'
#' The user can opt for using \code{\link{optim}}'s version of \code{CG} and \code{L-BFGS-B}. The
#' implementation in \code{\link{optim}} of the methods does not incorporate improvements of the
#' methods in the recent decades, but they have stood the test of time.
#'
#' If \code{parinit} were not supplied and the distance between \code{lower} and \code{upper} is infinite,
#' the initialisation point of the optimisation is drawn from a uniform distribution ranging [-1,1]
#' distribution. If initalisation were not supplied, but the distance between \code{lower} and \code{upper}
#' is finite, then the initialisation is drawn from a uniform distribution ranging
#' [\code{lower}, \code{upper}].
#' 
#' 
#' @param mod          An object of class \code{\link{glinv}}.
#' @param parinit      A vector, parameter for initialisation of the optimisation routine.
#' @param method       One of \code{L-BFGS-B}, \code{CG}, \code{BB}, or any other methods which is accepted by optim.
#' @param lower        A vector of lower bounds on the parameters.
#' @param upper        A vector of upper bounds on the parameters.
#' @param use_optim    If true, use optim's version of \code{L-BFGS-B} and \code{CG}.
#' @param control      Options to be passed into each the underlying optimisation routine's \code{control}
#'                     argument.
#' @param ...          Not used.
#' @return             \code{fit.glinv} returns a list containing at least the following elements:
#'                     \item{mlepar}{The maximum likelihood estimate.}
#'                     \item{loglik}{The log-likelihood at the maximum likelihood estimate.}
#'                     \item{score}{The gradient of log-likelihood at the maximum likelihood estimate.}
#'                     \item{convergence}{Zero if the optimisation routine has converged successfully.}
#'                     \item{message}{A message from the optimisation routine.}
#' @export
fit.glinv = function (mod, parinit=NULL, method='CG', lower=-Inf, upper=Inf, use_optim=F, control=list(), ...) {
  if (is.null(parinit))
    parinit = if (is.finite(lower) && is.finite(upper))
                stats::runif(mod$nparams, lower, upper)
              else
                stats::runif(mod$nparams, -1, 1)
  result = if (!use_optim && method=='L-BFGS-B') {
             r = lbfgsb3c::lbfgsb3c(par   = parinit,
                                    fn    = default_ifwrong(fit_fnwrap(lik(mod), -1),  1,          -Inf),
                                    gr    = default_ifwrong(fit_grwrap(grad(mod),-1), mod$nparams,   NA),
                                    lower = lower,
                                    upper = upper,
                                    control= list_set_default(control, list(trace = T, maxit = 500000)))
             names(r)[which(names(r) == 'par')]   = 'mlepar'
             names(r)[which(names(r) == 'value')] = 'loglik'
             names(r)[which(names(r) == 'grad')]  = 'score'
             r
           } else if (!use_optim && method=='CG') {
             r = Rcgmin::Rcgmin(par    = parinit,
                            fn     = default_ifwrong(fit_fnwrap(lik(mod), -1),  1,          -Inf),
                            gr     = default_ifwrong(fit_grwrap(grad(mod),-1), mod$nparams,   NA),
                            lower = lower,
                            upper = upper,
                            control= list_set_default(
                              control,
                              list(trace = 1,
                                   maxit = 500000,
                                   tol=mod$nparams*mod$nparams*5e-10))) # Default is about (npar^2)*1e-16.
             names(r)[which(names(r) == 'par')]   = 'mlepar'
             names(r)[which(names(r) == 'value')] = 'loglik'
             r = c(r, list(score = grad(mod)(r[['mlepar']])))
             r
           } else if (!use_optim && method=='BB') {
             r = BB::BBoptim(par    = parinit,
                         fn     = default_ifwrong(fit_fnwrap(lik(mod), -1),  1,          -Inf),
                         gr     = default_ifwrong(fit_grwrap(grad(mod),-1), mod$nparams,   NA),
                         lower = lower,
                         upper = upper,
                         control= list_set_default(control, list(maxit = 500000, gtol=1e-3, ftol=1e-8)))
             names(r)[which(names(r) == 'par')]      = 'mlepar'
             names(r)[which(names(r) == 'value')]    = 'loglik'
             r = r[-which(names(r) == "gradient")]
             r = c(r, list(score = grad(mod)(r[['mlepar']])))
             r
           } else {
             r = stats::optim(par    = parinit,
                   fn     = default_ifwrong(fit_fnwrap(lik(mod), -1),  1,          -Inf),
                   gr     = default_ifwrong(fit_grwrap(grad(mod),-1), mod$nparams,   NA),
                   lower  = lower,
                   upper  = upper,
                   control= list_set_default(control, list(maxit=500000,trace=1)),
                   method = method)
                          names(r)[which(names(r) == 'par')]   = 'mlepar'
             names(r)[which(names(r) == 'value')] = 'loglik'
             r = c(r, list(score = grad(mod)(r[['mlepar']])))
             r
           }
  names(result$mlepar) = names(parinit)
  result
}


#' Estimate the variance-covariance matrix of the maximum likelihood estimator.
#'
#' \code{varest} estimates the uncertainty of an already-computed maximum likelihood estimate.
#' 
#' If \code{method} is \code{analytical} then the covariance matrix is estimated by inverting the
#' negative analytically-computed Hessian at the maximum likelihood estimate; if it is 
#' \code{mc} then the estimation is done by using Spall's Monte Carlo simultaneous perturbation method;
#' if it is \code{linear} then it is done by the "delta method", which approximates the user
#' parameterisation with its first-order Taylor expansion.
#'
#' The \code{analytical} method requires that \code{parhess} was specified when `mod` was created.
#' The \code{linear} method does not use the curvature of the reparameterisation and its result is
#' sometimes unreliable; but it does not require the use of \code{parhess}. The \code{mc} method also
#' does not need \code{parjacs}, but the it introduces an additional source complexity and random noise
#' into the estimation; and a large number of sample may be needed.
#'
#' The \code{control.mc} can have the following elements:
#' \describe{
#'     \item{Nsamp}{Integer. Number of Monte Carlo iteration to run. Default is 10000.}
#'     \item{c}{Numeric. Size of perturbation to the parameters. Default is 0.005.}
#'     \item{quiet}{Boolean. Whether to print progress and other information or not. Default is \code{TRUE}.}
#' }
#'
#' 
#' @param mod                  An object of class \code{glinv}
#' @param fitted               Either an object returned by \code{fit.glinv} or a vector of length \code{mod$nparams}
#'                             that contains the maximum likelihood estimate.
#' @param method               Either `analytical', `linear' or `mc'. It specifies how the covariance matrix
#'                             is computed.
#' @param numDerivArgs         Arguments to pass to \code{numDeriv::\link[numDeriv]{jacobian}}. Only used if the
#'                             user did not supply \code{parjacs} when constructing \code{mod}.
#' @param store_gaussian_hessian If \code{TRUE} and \code{method} is not \code{mc}, the returned list will contain
#'                             a (usually huge) Hessian matrix \code{gaussian_hessian} with respect to the Gaussian
#'                             parameters \eqn{\Phi, w, V'}. This option significantly increases the amount of memory
#'                             the function uses, in order to store the matrix.
#' @param control.mc           A list of additional arguments to pass to the \code{mc} method.
#' @param ...                  Not used.
#' @return                     A list containing
#'                             \item{vcov}{The estimated variance-covariance matrix of the maximum likelihood estimator.}
#'                             \item{mlepar}{The maximum likelihood estimator passed in by the user.}
#'                             \item{hessian}{The Hessian of the log-likelihood at the maximum likelihood estimate. Only exists when \code{method} is not \code{mc}}
#'                             \item{gaussian_hessian}{Optional, only exists when `store_gaussian_hessian` is TRUE.}
#' @references                 Spall JC. Monte Carlo computation of the Fisher information matrix in nonstandard settings. Journal of Computational and Graphical Statistics. 2005 Dec 1;14(4):889-909.
#' @export
varest = function (mod, ...) UseMethod('varest')
#' @rdname varest
#' @export
varest.glinv = function (mod,
                         fitted,
                         method                 = 'analytical',
                         numDerivArgs           = list(method='Richardson', method.args=list(d=.5,r=3)),
                         store_gaussian_hessian = F,
                         control.mc             = list(),
                         ...) {
  if (is.list(fitted)) {
    if (is.double(fitted$mlepar))
      mlepar = fitted$mlepar
    else
      stop('Invalid argument: `fitted`')
  } else
    mlepar = fitted
  if (length(mlepar) != mod$nparams)
    stop(sprintf("Your model should have %d parameters but I got %d", mod$nparams, length(mlepar)))
  jac = if (is.null(mod$parjacs))
          function (x) c(do.call(numDeriv::jacobian, c(list(mod$gaussparams_fn, x), numDerivArgs)))
        else
          function (x) mod$gaussparams_jac(x)
  r = switch(method,
         'linear'={
           if (store_gaussian_hessian) {
             Hvwphi = hess(mod$rawmod, as.double(mod$gaussparams_fn(mlepar)))
             J = jac(mlepar)
             H = crossprod(J, Hvwphi) %*% J
             invH = hessian_inv_warn(-H, 'The linearly approximated negative Hessian')
             rownames(invH) = colnames(invH) = rownames(H) = colnames(H) = names(mlepar)
             list(vcov = invH, hessian = H, gaussian_hessian=Hvwphi, mlepar = mlepar)
           } else {
             H = hess(mod$rawmod, as.double(mod$gaussparams_fn(mlepar)), directions = jac(mlepar))
             rownames(invH) = colnames(invH) = rownames(H) = colnames(H) = names(mlepar)
             list(vcov = hessian_inv_warn(-H), hessian = H, mlepar = mlepar)
           }
         },
         'analytical'={
           H = hess(mod, numDerivArgs = numDerivArgs, store_gaussian_hessian = store_gaussian_hessian)(mlepar)
           invH = hessian_inv_warn(-H, 'The negative Hessian')
           rownames(invH) = colnames(invH) = rownames(H) = colnames(H) = names(mlepar)
           r = list(vcov = invH, hessian = H, mlepar = mlepar)
           if (store_gaussian_hessian)
             r[['gaussian_hessian']] = attr(H, 'gaussian_hessian')
           r
         },
         'mc'={
           control.mc = list_set_default(control.mc, list(Nsamp = 2000L, c = 0.005, quiet = F))
           c = control.mc$c
           Nsamp = control.mc$Nsamp
           quiet = control.mc$quiet
           p = length(mlepar)
           g = grad.glinv(mod)
           H = matrix(0.0,p,p)
           del        = double(p)
           theta_plus = double(p)
           theta_minus= double(p)
           DG         = double(p)
           A          = double(p^2)
           dim(A)     = c(p,p)
           gpar = mod$gaussparams_fn(mlepar)
           simdat = rglinv.glinv_gauss(mod$rawmod, gpar, Nsamp)
           for (j in 1L:Nsamp) {
             if (!(quiet) && (j %% 50)==0) cat(sprintf('%d/%d\n', j, Nsamp))
             set_tips.glinv_gauss(mod$raw, simdat[[j]])
             del[] = (rbinom(p,1,0.5)-0.5)*2*c
             theta_plus[]  = mlepar + del
             theta_minus[] = mlepar - del
             DG[] = g(theta_plus) - g(theta_minus)
             A[,] = tcrossprod(matrix(DG/2),1/del)
             H[,] = H[,] - (A + t(A))/2
           }
           invH  = hessian_inv_warn(H/Nsamp, 'The Monte Carlo estimated variance-covariance matrix')
           rownames(invH) = colnames(invH) = names(mlepar)
           list(vcov = invH, mlepar = mlepar)
         })
  r
}

#' Getting marginal confidence interval for GLInv model
#'
#' \code{marginal_ci} computes the marginal confidence interval for each parameters
#' using the variance-covariance matrix output by `varest.glinv`.
#'
#' @param varest_result    The output from `varest.glinv`.
#' @param lvl              Confidence level. Default to 95 percent.
#' @return                 A matrix \eqn{p}-by-2 matrix where \eqn{p} is the number of parameters.
#'                         The first column is the lower limits and second column is the upper limits.
#' @export
marginal_ci = function (varest_result, lvl = 0.95) {
  c_alpha = qnorm(1-(1-lvl/2))
  D = diag(varest_result$vcov)
  if (any(D <= 0)) stop('Cannot compute marginal confidence interval, because some diagonal elements of the variance-covariance is non-positive.')
  upper = varest_result$mlepar + sqrt(diag(varest_result$vcov)) * c_alpha
  lower = varest_result$mlepar - sqrt(diag(varest_result$vcov)) * c_alpha
  r = cbind(Lower=lower, Upper=upper)
  rownames(r) = names(varest_result$mlepar)
  r
}
