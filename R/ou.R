#' Get the number of parameters of some pre-defined models
#' 
#' \code{nparams_ou_*} returns the number of parameters of their respective model.
#'
#' @param k    An Integer. The total number of dimensions of the multivariate traits.
#' @export
nparams_ou = function (k) k*k+k+(k*(k+1L))%/%2L
#' @rdname nparams_ou
#' @export
nparams_ou_diagH = function (k) 2L*k+(k*(k+1L))%/%2L
#' @rdname nparams_ou
#' @export
nparams_ou_logdiagH = function (k) 2L*k+(k*(k+1L))%/%2L
#' @rdname nparams_ou
#' @export
nparams_ou_symH  = function (k) k+2L*(k*(k+1L))%/%2L
#' @rdname nparams_ou
#' @export
nparams_ou_spdH  = function (k) k+2L*(k*(k+1L))%/%2L

## ou_par2list = function (par) {
##   k     = sqrt(9+24*length(par))/6-1/2
##   H     = par[1:(k*k)]
##   dim(H)= c(k,k)
##   theta = par[(k*k+1L):(k*k+k)]
##   sig_x = par[(k*k+k+1L):length(par)]
##   sig_x_M = matrix(0,k,k)
##   sig_x_M[lower.tri(sig_x_M, diag=T)] = sig_x
##   list(H     = H,
##        theta = theta,
##        sig_x = sig_x_M)
## }
## ou_list2par = function (parlist) {
##   c(parlist$H,
##     parlist$theta,
##     parlist$sig_x[lower.tri(parlist$sig_x, diag=T)])
## }


#' Parameterisation functions of Ornstein-Uhlenbeck model
#' 
#' \code{oupar} is a function that maps from the Ornstein-Uhlenbeck model
#' parameters to the Gaussian parametersation.
#'
#' By multivariate Ornstein-Uhlenbeck process, we mean
#'     \deqn{dx(t) = -H(x(t) - \theta)dt + \Sigma_x dW(t)}
#' where \eqn{H} is a \eqn{k}-by-\eqn{k} matrix with real entries,
#' \eqn{\theta} is any real \eqn{k}-vector, \eqn{\Sigma_x} is a
#' lower-triangular matrix, \eqn{W(t)} is the Brownian motion process.
#' The parameters of this model is \eqn{(H,\theta,\Sigma_x)},
#' therefore \eqn{k^2+k+k(k+1)/2} dimensional.
#' 
#' The \code{par} arguemnt is the concatenation of major-column-flattened
#' \eqn{H}, \eqn{\theta}, and the major-column-flattened lower-triangular part
#' of \eqn{\Sigma_x}.
#'
#' @param par     A numeric vector containing the joint vector of the
#'                Ornstein-Uhlenbeck drift matrix, long-term mean,
#'                and volitality matrix, which is a lower-triangular
#'                Cholesky factor.
#' @param t       Branch length of the currently processing node.
#' @param ...     Unused in these functions. Their existence is needed because
#'                \code{\link{lik.glinv}} etc. always pass us four arguments.
#'                See \code{\link{lik.glinv}} for details.
#' @return        \code{oupar} returns the a vector of concatenated \eqn{(\Phi, w, V')},
#'                where \eqn{V'} is the lower triangular part of \eqn{V}. \code{oujac}
#'                returns the Jacobian matrix of \code{oupar}. \code{ouhess} returns
#'                a list of three 3D arrays, named \code{Phi}, \code{w}, \code{V} respectively inside the list, in which
#'                \code{ouhess(...)$Phi[m,i,j]} contains
#'                the cross second-order partial derivative of \eqn{\Phi_m} (here we treat the matrix \eqn{\Phi} as a
#'                column-major-flattened vector) with respect to the \eqn{i}-th and\eqn{j}-th user parameters;
#'                and \code{ouhess(...)$w[m,i,j]} and \code{((parhess[[i]])(...))$V[m,i,j]}
#'                analogously contains second-order derivative of \eqn{w_m} and \eqn{V'_m}.
#'                
#' @export
oupar = function (par, t, ...) {
  ## L == k*k+k+k*(k+1)/2 => 3(k*k) + 3k - 2*L == 0
  ## outfmt="vector"
  opts = list(...)
  if (is.null(opts[['outfmt']])) outfmt = 'vector'
  if (is.numeric(par)) {
    k     = sqrt(9+24*length(par))/6-1/2
    H     = par[1L:(k*k)]
    theta = par[(k*k+1L):(k*k+k)]
    sig_x = par[(k*k+k+1L):length(par)]
  } else if (is.list(par)) {
    k     = length(theta)
    H     = c(par[['H']])
    theta = c(par[['theta']])
    sig_x = c(par[['sig_x']])
    if (! ((is.numeric(H)     && length(H) == k*k) ||
            is.numeric(theta) && length(theta) == k||
            is.numeric(sig_x) && length(sig_x) == (k*(k+1L))%/%2L))
      stop("Invalid argument: `par` is malformed")
  } else  stop('`par` must be either a list or a numeric vector')
  if      (identical(outfmt, 'vector')) vecfmt=T
  else if (identical(outfmt, 'list'))   vecfmt=F
  else    stop('outfmt must be one of "vector" or "list"')
  ;                      mode(H)       = 'double'
  ;                      mode(t)       = 'double'
  ;                      mode(theta)   = 'double'
  ;                      mode(sig_x)   = 'double'
  ;                      mode(k)       = 'integer'
  V = numeric((k*(k+1L))%/%2L);mode(V) = 'double'
  w = matrix(0,k);       mode(w)       = 'double'
  Phi = matrix(0,k,k);   mode(Phi)     = 'double'
  P = matrix(0,k,k);     mode(P)       = 'complex'
  invP = matrix(0,k,k);  mode(invP)    = 'complex'
  Lambda = complex(k);   mode(Lambda)  = 'complex'
  info = integer(1);     mode(info)    = 'integer'
  lwsp = 18L*k*k+19L+1L; mode(lwsp)    = 'integer'
  lzwsp = 10L*k*k;       mode(lzwsp)   = 'integer'
  wsp = double(lwsp);    mode(wsp)     = 'double'
  zwsp = complex(lzwsp); mode(zwsp)    = 'complex'
  eigavail = 0;          mode(eigavail)= 'integer'
  res = .C('d0geouvwphi_', H,k,t,theta,sig_x,V,w,Phi,P,invP,Lambda,wsp,lwsp,zwsp,lzwsp,eigavail,info)
  if (res[[17L]] != 0)    stop('Cannot eigen-decompose `H`')
  if (vecfmt)            return( unlist(res[c(8L,7L,6L)]) )
  else                   return( { R=res[c(8L,7L,6L)]; names(R)=c('Phi','w','V'); R } )
}

#' Jacobian of matrix of the Ornstein-Uhlenbeck model
#' 
#' \code{oujac} accepts the same arguments as \code{oupar} and returns the
#' Jacobian matrix of \code{oupar}.
#'
#' @rdname oupar
#' @export
oujac = function (par, t, ...) {
  if (is.numeric(par)) {
    k     = as.integer(sqrt(9+24*length(par))/6-1/2)
    H     = par[1L:(k*k)]
    theta = par[(k*k+1L):(k*k+k)]
    sig_x = par[(k*k+k+1L):length(par)]
  } else if (is.list(par)) {
    k     = length(theta)
    H     = c(par[['H']])
    theta = c(par[['theta']])
    sig_x = c(par[['sig_x']])
    if (! ((is.numeric(H)     && length(H) == k*k) ||
            is.numeric(theta) && length(theta) == k||
            is.numeric(sig_x) && length(sig_x) == (k*(k+1L))%/%2L))
        stop("Invalid argument: `par` is malformed")
  } else  stop('`par` must be either a list or a numeric vector')
  ;                                       mode(t)       = 'double'
  ;                                       mode(k)       = 'integer'
  hts = c(H,theta,sig_x);                 mode(hts)     = 'double'
  P = matrix(0,k,k);                      mode(P)       = 'complex'
  invP = matrix(0,k,k);                   mode(invP)    = 'complex'
  Lambda = complex(k);                    mode(Lambda)  = 'complex'
  info = integer(1);                      mode(info)    = 'integer'
  lzwsp = max(8L*k*k, 3L*k*k+4L*k^4L);    mode(lzwsp)   = 'integer'
  zwsp = complex(lzwsp);                  mode(zwsp)    = 'complex'
  lwsp = 18L*k*k+as.integer(3*k^4)+19L+1L;mode(lwsp)    = 'integer'
  wsp = double(lwsp);                     mode(wsp)     = 'double'
  eigavail = 0;                           mode(eigavail)= 'integer'
  jac   = matrix(0,k*k+k+(k*(k+1L))%/%2L,k*k+k+(k*(k+1L))%/%2L); mode(jac)   = 'double'
  ougejacres = .C('ougejac_', t,k,hts,P,invP,Lambda,wsp,lwsp,zwsp,lzwsp,eigavail,jac,info)
  if (ougejacres[[13]] != 0)   stop('Cannot eigen-decompose `H`')
  ougejacres[[12]]
}


#' Jacobian of matrix of the Ornstein-Uhlenbeck model
#'
#' \code{ouhess} returns accepts the same arguments as \code{oupar}
#' and returns all the second derivatives \code{oupar}. The returned
#' values are consistent with the format required by \code{link{glinv}}.
#' 
#' @rdname oupar
#' @export
ouhess = function (par, t, ...) {
  k = as.integer(sqrt(9+24*length(par))/6-1/2)
  npar = length(par)
  mode(par) = 'double'
  mode(t)   = 'double'
  sig = .C('lnunchol_', as.double(par[((k*k)+k+1L):npar]), k,
           double(k*k), k*k, matrix(0.,k,k), 0L)[[5]]
  eig = eigen(matrix(par[1L:(k*k)],k,k))
  P = as.complex(eig[['vectors']])
  invP = as.complex(solve(eig[['vectors']]))
  Lambda = as.complex(eig[['values']])
  res = .C('hphiha_', t, par[1L:(k*k)],
           k, P,invP,Lambda,
           array(0.0, c(k*k,k*k,k*k)),
           complex(k^6+2L*(k*k)+3), integer(k^6+2L*(k*k)+3),
           0L)
  if (res[[10]] != 0) stop('Error executing hphiha_()')
  hphiha = res[[7]]
  list('V' = {
    hv = array(0., c((k*(k+1L))%/%2L, npar, npar))
    hvhares = .C('hvha_', t, sig, par[1L:(k*k)],
             k, P, invP, Lambda,
             array(0., c((k*(k+1L))%/%2L,k*k,k*k)),
             double(2L*(k*k)), 2L*(k*k),
             complex(k^6+4*(k*k)+3*k),as.integer(k^6)+4*(k*k)+3L*k,
             1L,0L)
    if (hvhares[[14]] != 0) stop('Error executing hvha_()')
    hv[1L:((k*(k+1L))%/%2L),1L:(k*k),1L:(k*k)] = hvhares[[8]]
    hvdadlres = .C('hvdadl_', t, par[1L:(k*k)],
              k, par[((k*k)+k+1L):npar],
              P,invP,Lambda, array(0., c((k*(k+1L))%/%2L,k*k,(k*(k+1L))%/%2L)),
              double(4L*(k*k)), 4L*(k*k),
              complex(k^4+k*k),as.integer(k^4)+k*k,
              0L)
    if (hvdadlres[[13]] != 0) stop('Error executing hvdadl_()')
    hv[1L:((k*(k+1L))%/%2L),(k*k+k+1L):npar,1L:(k*k)] =
        aperm(hv[1L:((k*(k+1L))%/%2L),1L:(k*k),(k*k+k+1L):npar] <- hvdadlres[[8]], c(1L,3L,2L))
    hv[1L:((k*(k+1L))%/%2L),(k*k+k+1L):npar,(k*k+k+1L):npar] =
        .C('hvhl_', t, k, par[((k*k)+k+1L):npar],
           P, invP, Lambda,
           double(4L*k*k), 4L*k*k,
           complex(2L*k*k), 2L*k*k,
           array(0.0, c((k*(k+1L))%/%2L,(k*(k+1L))%/%2L,(k*(k+1L))%/%2L)))[[11]]
    hv
  },
  'w' = {
    hw = array(0., c(k,npar,npar))
    hw[1L:k,1L:(k*k),(k*k+1L):(k*k+k)] =
        aperm(hw[1L:k,(k*k+1L):(k*k+k),1L:(k*k)]<-
                  .C('hwdthetada_', k,
                     .C('dphida_', t, k,
                        P,invP,Lambda,
                        matrix(0.,k*k,k*k),
                        complex(as.integer(k^4)+k*k+2L),as.integer(k^4)+k*k+2L)[[6]],
                     array(0., c(k,k,k*k)))[[3]],
              perm=c(1L,3L,2L))
    hw[1L:k,1L:(k*k),1L:(k*k)] =
        .C('hwha_', k, hphiha, par[((k*k)+1L):(k*k+k)],
           array(0., c(k,k*k,k*k)))[[4]]
    hw
  },
  'Phi' = {
    hphi = array(0., c(k*k,npar,npar))
    hphi[1L:(k*k),1L:(k*k),1L:(k*k)] = hphiha[,,]
    hphi
  })
}

#' Restrict the drift matrix of OU model.
#'
#' \code{ou_diagH} restricts the drift matrix of an OU model to
#' diagonal matrices.
#'
#' \subsection{How reparametrisation and restriction works}{
#' 
#' In the simplest form, without any restriction or reparametrisation, the user typically
#' needs to pass \code{oupar}, \code{oujac}, \code{ouhess}, all of which are simply
#' functions which maps from the OU parameters \eqn{(H,\theta,\Sigma_x)} to the Gaussian
#' paramters \eqn{(\Phi_i,w_i,V'_i)} for each node. For example:
#' \preformatted{
#'         mod.full = glinv(tree, x0, my_data,
#'                          parfns  = oupar,
#'                          pardims = nparams_ou(k),
#'                          parjacs = oujac,
#'                          parhess = ouhess)
#' }
#' If one would like to restrict \eqn{H} to only positively definite diagonal matrices,
#' then the call should become
#' \preformatted{
#'         mod.pddiag = glinv(tree, x0, my_data,
#'                            parfns  = ou_logdiagH(oupar),
#'                            pardims = nparams_ou_logdiagH(k),
#'                            parjacs = dou_logdiagH(oujac),
#'                            parhess = hou_logdiagH(ouhess))
#' }
#' Note that there is a naming convention that \code{ou_*} should be applied to `oupar`,
#' \code{dou_*} to `oujac`, and \code{hou_*} to `ouhess`. \code{d} stands for `derivative'
#' and \code{h} stands for `Hessian'.
#' 
#' In the above call, ou_logdiagH(oupar) accepts the \code{oupar} function as argument
#' and returns a new function. This new function behaves the same way as oupar itself,
#' except that it expects its first argument (which is the model parameters) to be of
#' lower dimension, only consisting of \eqn{(h,\theta,\Sigma_x)} where \eqn{h} is the
#' diagonal vector of \eqn{H}. The following example should be illustrative:
#' \preformatted{
#'         f = ou_logdiagH(oupar)
#'         par.full = list(H     = matrix(c(3,0,0,2),2,2),
#'                         theta = c(4,5),
#'                         sig_x = c(1,0,1))
#'         par.restricted = list(H     = log(diag(par.full$H)),
#'                               theta = par.full$theta,
#'                               sig_x = par.full$sig_x)
#'         print(all.equal(f(unlist(par.restricted),1,NULL,NULL),
#'                         oupar(unlist(par.full),1,NULL,NULL)))
#'         # [1] TRUE
#' }
#' }
#' 
#' \subsection{Details about each pre-defined restrictions}{
#' The following table summarises all the pre-defined \code{ou_*} functions.
#' \tabular{ll}{
#'   \strong{R function}   \tab \strong{Restriction}\cr
#'   \code{*_diagH}    \tab \eqn{(h,\theta,\Sigma_x)}, with \eqn{h=diag(H)}\cr
#'   \code{*_logdiagH} \tab \eqn{(log(h),\theta,\Sigma_x)}, with \eqn{h=diag(H)}\cr
#'   \code{*_symH}     \tab \eqn{(L,\theta,\Sigma_x)}, with \eqn{L} being lower-triangular part of H\cr
#'   \code{*_spdH, log=F}  \tab \eqn{(L,\theta,\Sigma_x)}, with \eqn{L} being Cholesky factor of H\cr
#'   \code{*_spdH, log=T}  \tab \eqn{(L',\theta,\Sigma_x)} where \eqn{L'} equals \eqn{L}, except that on the diagonals \eqn{L'_i} = \eqn{log L_i}
#' }
#' By Cholesky factor, we mean the only the non-zero part of the lower-triangular Cholesky factor.
#' }
#'
#' @param parfn      A function that maps from the user-parametrisation to the underlying Gaussian parameters.
#'                Each of them returns a vector of concatenated \eqn{(\Phi, w, V')}, where \eqn{V'} is the lower triangular
#'                part of \eqn{V}, and accepts four arguments: a vector of parameters whose length is specified
#'                by the \code{pardims} argument to the \code{glinv_gauss} function, the branch length leading to the currently processing node, 
#'                a vector of factors with three levels indicating which dimensions are missing or lost in the mother of
#'                the current node, and a vector of factors with the same three levels indicating missingness of the current
#'                node.
#' @param jacfn     A function that accepts the same arguments as \code{parfn} and returns the Jacobian
#'                of \code{parfn}.
#' @param hessfn A function that accepts the same arguments as \code{parfns} and returns a list of three 3D arrays,
#'                named \code{Phi}, \code{w}, \code{V} respectively inside the list. \code{((hessfn)(...))$Phi[m,i,j]}
#'                contains the cross second-order partial derivative of \eqn{\Phi_m} (here we treat the matrix
#'                \eqn{\Phi} as a column-major-flattened vector) with respect to the \eqn{i}-th and\eqn{j}-th parameters
#'                in the joint \eqn{(H,\theta,\Sigma_x)} vector, and
#'                \code{((hessfn)(...))$w[m,i,j]} and \code{((hessfn)(...))$V[m,i,j]}
#'                analogously contains second-order derivative with respect to \eqn{w_m} and \eqn{V'_m}.
#' @param log     Whether or not some elements of the parameters should be passed to
#'                their logarithm in the resulting reparametrisation. See the Details
#'                section.
#' @return        \code{ou_*}, \code{dou_*}, and \code{hou_*} returns, respectively,
#'                a function that accepts the same form of arguments, and returns the
#'                same form of output, as the \code{parfn}, \code{jacfn}, and
#'                \code{hessfn} arguments. These function should be ready to be passed
#'                to the \code{parfns}, \code{parjacs}, and \code{parhess} arguments
#'                of \code{\link{lik.glinv}}.
#' @export
ou_diagH = function (parfn) {
  parfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    parfn(c(if (k==1L) par[1L] else diag(par[1L:k]),par[(k+1L):length(par)]), ...)
  }
}

#' @rdname ou_diagH
#' @export
dou_diagH = function (jacfn) {
  jacfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    J = jacfn(c(if(k==1L) par[1L] else diag(par[1L:k]),par[(k+1L):length(par)]), ...)
    J[,c((1L:k)+((1L:k)-1L)*k, (k*k+1L):ncol(J)),drop=F]
  }
}

#' @rdname ou_diagH
#' @export
hou_diagH = function (hessfn) {
  hessfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    res = hessfn(c(if(k==1L) par[1L] else diag(par[1L:k]),par[(k+1L):length(par)]), ...)
    wh = c((1L:k)+((1L:k)-1L)*k, (k*k+1L):(dim(res[['V']])[2]))
    res[['V']]   = res[['V']][,wh,wh,drop=F]
    res[['w']]   = res[['w']][,wh,wh,drop=F]
    res[['Phi']] = res[['Phi']][,wh,wh,drop=F]
    res
  }
}

#' Restrict the drift matrix of OU model to positive diagonal matrices.
#' 
#' \code{ou_logdiagH} restricts the drift matrix of an OU model to
#' positively definite diagonal matrices. A diagonal matrix is positively
#' definite if and only if all of its diagonal entries are positive;
#' therefore it parametrises the diagonals to their logarithm to
#' facilitate unconstrained optimisation.
#'
#' @rdname ou_diagH
#' @export
ou_logdiagH = function (parfn) {
  parfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    parfn(c(if (k==1L) exp(par[1L:k]) else diag(exp(par[1L:k])),par[(k+1L):length(par)]), ...)
  }
}

#' @rdname ou_diagH
#' @export
dou_logdiagH = function (jacfn) {
  jacfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    J = jacfn(c(if(k==1L) exp(par[1L:k]) else diag(exp(par[1L:k])),par[(k+1L):length(par)]), ...)
    J = J[,c((1L:k)+((1L:k)-1L)*k, (k*k+1L):ncol(J)),drop=F]
    J[,1L:k] = sweep(J[,1L:k,drop=F], MARGIN=2, exp(par[1L:k]), `*`)
    J
  }
}

#' @rdname ou_diagH
#' @export
hou_logdiagH = function (hessfn) {
  hessfn
  function (par, ...) {
    k = as.integer(sqrt(25/4+2*length(par)) - 5/2)
    eh = exp(par[1L:k])
    npar_orig = k*k+k+(k*(k+1L))%/%2L
    res = hessfn(c(if (k==1L) eh else diag(eh),par[(k+1L):length(par)]), ...)
    wh = c((1L:k)+((1L:k)-1L)*k, (k*k+1L):npar_orig)
    res[['V']]   = res[['V']][,wh,wh,drop=F]
    res[['w']]   = res[['w']][,wh,wh,drop=F]
    res[['Phi']] = res[['Phi']][,wh,wh,drop=F]
    gstart = INFO__$mod$gausssegments[INFO__$node_id,'start']
    gend   = INFO__$mod$gausssegments[INFO__$node_id,'end']
    pstart = INFO__$mod$parsegments[INFO__$parfn_id,'start']
    pend   = INFO__$mod$parsegments[INFO__$parfn_id,'end']
    res[['V']][,1L:k,]   = .Call('Reinsumijk_j_2_ijk',res[['V']][,1L:k,,drop=F],  eh)
    res[['w']][,1L:k,]   = .Call('Reinsumijk_j_2_ijk',res[['w']][,1L:k,,drop=F],  eh)
    res[['Phi']][,1L:k,] = .Call('Reinsumijk_j_2_ijk',res[['Phi']][,1L:k,,drop=F],eh)
    res[['V']][,,1L:k]   = .Call('Reinsumijk_k_2_ijk',res[['V']][,,1L:k,drop=F],  eh)
    res[['w']][,,1L:k]   = .Call('Reinsumijk_k_2_ijk',res[['w']][,,1L:k,drop=F],  eh)
    res[['Phi']][,,1L:k] = .Call('Reinsumijk_k_2_ijk',res[['Phi']][,,1L:k,drop=F],eh)
    repa = INFO__[['reparametrisation_jacobian']]
    phid = dim(res[['Phi']])[1]
    wd   = dim(res[['w']])[1]
    for (i in seq_len(k)) {
      res[['V']][,i,i] = res[['V']][,i,i]     + repa[(gstart+phid+wd):gend, pstart+i-1L,drop=F]
      res[['w']][,i,i] = res[['w']][,i,i]     + repa[(gstart+phid):(gstart+phid+wd-1L),pstart+i-1L,drop=F]
      res[['Phi']][,i,i] = res[['Phi']][,i,i] + repa[gstart:(gstart+phid-1L),pstart+i-1L,drop=F]
    }
    res
  }
}

#' Restrict the drift matrix of OU model to symmetric matricies
#' 
#' \code{ou_symH} restricts the drift matrix of an OU model to
#' symmetric matrices. It parametrises the space of symmetric
#' matrices by the lower-triangular part.
#'
#' @rdname ou_diagH
#' @export
ou_symH = function (parfn) {
  parfn
  function (par, ...) {
    k = as.integer(sqrt(4+4*length(par))/2)-1L
    M = .Call('Rsylgecpy', as.double(par[1L:((k*(k+1L))%/%2L)]), k)
    parfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
  }
}

#' @rdname ou_diagH
#' @export
dou_symH = function (jacfn) {
  jacfn
  function (par, ...) {
    k = as.integer(sqrt(4+4*length(par))/2)-1L
    M = .Call('Rsylgecpy', par[1L:(((k*(k+1L))%/%2L))], k)
    J = jacfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
    dimJ = dim(J)
    k2half = (k*(k+1L))%/%2L
    res = matrix(0., dimJ[1], k2half + (dimJ[2] - k*k))
    dH = matrix(0., k,k)
    dH_lowt_mask = lower.tri(dH,diag=T)
    for (i in seq_len(nrow(res))) {
      dH[,] = J[i,1L:(k*k)]
      dH[,] = dH+t(dH) - if (k==1L) dH else diag(diag(dH))
      res[i,1L:k2half] = dH[dH_lowt_mask]
      res[i,(k2half+1L):ncol(res)] = J[i,((k*k)+1L):dimJ[2]]
    }
    res
  }
}

#' @rdname ou_diagH
#' @export
hou_symH = function (hessfn) {
  hessfn
  function (par, ...) {
    k = as.integer(sqrt(4+4*length(par))/2)-1L
    M = .Call('Rsylgecpy', par[1L:(((k*(k+1L))%/%2L))], k)
    res = hessfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
    res[['V']]   = .Call('Rhouchnsymh', res[['V']], k)
    res[['w']]   = .Call('Rhouchnsymh', res[['w']], k)
    res[['Phi']] = .Call('Rhouchnsymh', res[['Phi']], k)
    res
  }
}

#' Restrict the drift matrix of OU model to symmetric positively definite matrices.
#'
#' \code{ou_spdH} restricts the drift matrix of an OU model to
#' symmetric positively definite matrices, either by its Cholesky factor,
#' or by a modified Cholesky factor in which the diagonals are passed
#' to their logarithm to ensure that they are all positive. The likelihood
#' surface will be multi-modal in the non-logarithm version because any
#' arbitrary sign change in the Cholesky factor's diagonals does not change
#' the likelihood value.
#' 
#' @rdname ou_diagH
#' @export
ou_spdH = function (parfn, log=T) {
  parfn
  if (log) function (par, ...) {
    k = as.integer(sqrt(4+4*length(par))/2)-1L
    M = .Call('Rlnunchol', as.double(par[1L:(((k*(k+1L))%/%2L))]), k)
    parfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
  } else function (par, ...) {
    k = as.integer(sqrt(4+4*length(par)))%/%2L-1L
    M = .Call('Runchol', as.double(par[1L:(((k*(k+1L))%/%2L))]), k)
    parfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
  }
}


#' @rdname ou_diagH
#' @export
dou_spdH = function (jacfn, log=T) {
  jacfn
  if (log) function (par, ...) {
    k = as.integer(sqrt(4L+4L*length(par))/2)-1L
    mode(par) = 'double'
    M = .Call('Rlnunchol', par[1L:(((k*(k+1L))%/%2L))], k)
    J = jacfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
    mode(J) = 'double'
    DL = .Call('Rlnchnunchol', J, par[1L:(((k*(k+1L))%/%2L))], nrow(J), k)
    cbind(DL, J[,((k*k)+1L):ncol(J),drop=F])
  } else function (par, ...) {
    k = as.integer(sqrt(4+4*length(par))/2)-1L
    M = .Call('Runchol', as.double(par[1L:(((k*(k+1L))%/%2L))]), k)
    J = jacfn(c(M,par[((k*(k+1L))%/%2L+1L):length(par)]), ...)
    mode(J) = mode(par) = 'double'
    DL = .Call('Rchnunchol', J, par[1L:(((k*(k+1L))%/%2L))], nrow(J), k)
    cbind(DL, J[,(k*k+1L):ncol(J),drop=F])
  }
}


#' @rdname ou_diagH
#' @export
hou_spdH = function (hessfn, jacfn, log=T) {
  hessfn
  jacfn
  if (log) function (par, ...) {
    k = as.integer(sqrt(4L+4L*length(par))/2)-1L
    mode(par) = 'double'
    k2half = ((k*(k+1L))%/%2L)
    npar_orig = k*k+k+k2half
    par_orig = double(k*k+k+k2half)
    par_orig[1:(k*k)] = .Call('Rlnunchol', par[1L:k2half], k)
    par_orig[(k*k+1):(k*k+k+k2half)] = par[(k2half+1L):length(par)]
    res_orig  = hessfn(par_orig, ...)
    J = jacfn(par_orig, ...)
    ku = INFO__$mod$rawmod$dimtab[INFO__$node_id]
    kv = INFO__$mod$rawmod$dimtab[INFO__$parent_id]
    list(V   = .Call('Rhoulnspdh', res_orig[['V']],   k, par, J, ku*kv+ku),
         w   = .Call('Rhoulnspdh', res_orig[['w']],   k, par, J, ku*kv),
         Phi = .Call('Rhoulnspdh', res_orig[['Phi']], k, par, J, 0L))
  } else function (par, ...) {
    k = as.integer(sqrt(4L+4L*length(par))/2)-1L
    mode(par) = 'double'
    k2half = ((k*(k+1L))%/%2L)
    npar_orig = k*k+k+k2half
    par_orig = double(k*k+k+k2half)
    par_orig[1:(k*k)] = .Call('Runchol', par[1L:k2half], k)
    par_orig[(k*k+1):(k*k+k+k2half)] = par[(k2half+1L):length(par)]
    res_orig  = hessfn(par_orig, ...)
    J = jacfn(par_orig, ...)
    ku = INFO__$mod$rawmod$dimtab[INFO__$node_id]
    kv = INFO__$mod$rawmod$dimtab[INFO__$parent_id]
    list(V   = .Call('Rhouspdh', res_orig[['V']],   k, par, J, ku*kv+ku),
         w   = .Call('Rhouspdh', res_orig[['w']],   k, par, J, ku*kv),
         Phi = .Call('Rhouspdh', res_orig[['Phi']], k, par, J, 0L))
  }
}

