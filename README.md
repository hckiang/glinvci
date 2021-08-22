# Introduction

GLInvCI is a package that provides a framework for computing the maximum-likelihood estimates and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes, including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic comparative methods. The framework is designed to be flexible enough that the user can easily specify their own parameterisation and obtain the maximum-likelihood estimates and confidence intervals of their own parameters.

The model in concern is GLInv family, in which each species' traits evolve independently of each others after branching off from their common ancestor and for every non-root node. Let (Φₖ,wₖ,Vₖ) be a child node of (Φₖ,wₖ,Vₖ), and (Φₖ,wₖ,Vₖ), (Φₖ,wₖ,Vₖ) denotes the corresponding multivariate traits. We assume that (Φₖ,wₖ,Vₖ) is a Gaussian distribution with expected value (Φₖ,wₖ,Vₖ) and variance (Φₖ,wₖ,Vₖ), where the matrices (Φₖ,wₖ,Vₖ) are parameters independent of (Φₖ,wₖ,Vₖ) but can depend other parameters including (Φₖ,wₖ,Vₖ). The traits (Φₖ,wₖ,Vₖ) and (Φₖ,wₖ,Vₖ) can have different number of dimension.

# Installation

The following command should install the latest version of the package:

    install.packages('devtools')
    devtools::install(
      'https://git.sr.ht/~hckiang/glinvci/blob/latest-tarballs/glinvci_latest_main.tar.gz')

# High-level and low-level interface

The package contains two levels of user interfaces. The high-level interface, accessible through the `glinv` function, provides facilities for handling missing traits, lost traits, multiple evolutionary regimes, and most importantly, the calculus chain rule. The lower-level interface, accessible through the `glinv_gauss` function, allows the users to operate purely in the (Φₖ,wₖ,Vₖ) parameter space.

Most users should be satisfied with the high-level interface, even if they intend to write their own custom models.

# Using the high-level interface: Brownian Motion and OU Models

To fit a model using this package, generally you will need two main pieces of input data: a rooted phylogenetic tree and a matrix of trait values. The phylogenetic tree can be non-ultrametric and can potentially contain multifurcation. The matrix of trait values should have the same number of columns as the number of tips.

    ntips = 200
    k     = 2         # No. of trait dimensions
    tr    = ape::rtree(ntips)
    X     = matrix(rnorm(k*ntips), k, ntips) # Trait matrix
    x0    = rnorm(k)  # Root value

With the above material, we are ready to make a model object. We use OU as an example. Assume H is a positively definite matrix.

    mod = glinv(tr, x0, X,
                parfns = list(ou_logdiagH(ou_haltlost(oupar))),
                pardims = list(nparams_ou_diagH(k)),
                parjacs = list(dou_logdiagH(dou_haltlost(oujac))),
                parhess = list(hou_logdiagH(hou_haltlost(ouhess))))

The following code demostrates how to computing the model's likelihood, gradient, and Hessian at an arbitrarily specified pararmenter:

    H     = matrix(c(1,0,0,-1), k)
    theta = c(0,0)
    sig   = matrix(c(0.5,0,0,0.5), k)
    sig_x = t(chol(sig))
    par_init = c(H=diag(H),theta=theta,sig_x=sig_x[lower.tri(sig_x,diag=T)])
    print(par_init)
    print(lik(mod)(par_init))
    print(grad(mod)(par_init))
    print(hess(mod)(par_init))

The maximum likelihood estimates can be obtained by calling the `fit.glinv` method, as follows.

    fitted = fit(mod, par_init)
    print(fitted)

Once the model is fitted, one can estimate the variance-covariance matrix of the maximum-likelihood estimator using `varest`.

    v_estimate = varest(mod, fitted)

The marginal confidence interval can be obtained by calling `marginal_ci` on the object returned by `varest`.

    print(marginal_ci(v_estimate, lvl=0.95))
