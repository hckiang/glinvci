# Introduction

GLInvCI is a package that provides a framework for computing the maximum-likelihood estimates and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes, including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic comparative methods. The framework is designed to be flexible enough that the user can easily specify their own parameterisation and obtain the maximum-likelihood estimates and confidence intervals of their own parameters.

The model in concern is GLInv family, in which each species' traits evolve independently of each others after branching off from their common ancestor and for every non-root node. Let k be a child node of j, and zₖ, zⱼ denotes the corresponding multivariate traits. We assume that zₖ|zⱼ is a Gaussian distribution with expected value wₖ+Φₖzⱼ and variance Vₖ, where the matrices (Φₖ,wₖ,Vₖ) are parameters independent of zₖ but can depend other parameters including tₖ. The traits zₖ and zⱼ can have different number of dimension.

# Installation

The following command should install the latest version of the package:

    install.packages('devtools')
    devtools::install(
      'https://git.sr.ht/ hckiang/glinvci/blob/latest-tarballs/glinvciₗatestₘain.tar.gz')

# High-level and low-level interface

The package contains two levels of user interfaces. The high-level interface, accessible through the `glinv` function, provides facilities for handling missing traits, lost traits, multiple evolutionary regimes, and most importantly, the calculus chain rule. The lower-level interface, accessible through the `glinvgauss` function, allows the users to operate purely in the (Φη,wη,Vη) parameter space.

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
                parfns = list(ouₗogdiagH(ouₕaltlost(oupar))),
                pardims = list(nparamsₒudiagH(k)),
                parjacs = list(douₗogdiagH(douₕaltlost(oujac))),
                parhess = list(houₗogdiagH(houₕaltlost(ouhess))))

    ### — STEP 3: Try getting the likelihood, gradient etc.
    H     = matrix(c(1,0,0,-1), k)
    theta = c(0,0)
    sig   = matrix(c(0.5,0,0,0.5), k)
    sigₓ = t(chol(sig))
    parᵢnit = c(H=diag(H),theta=theta,sigₓ=sigₓ[lower.tri(sigₓ,diag=T)])
    print(parᵢnit)
    print(lik(mod)(parᵢnit))
    print(grad(mod)(parᵢnit))
    print(hess(mod)(parᵢnit))

    ### — STEP 4: Fitting a model
    fitted = fit(mod, parᵢnit)
    print(fitted)

    ### — STEP 5: Estimating variance-covariance of the MLE
    vₑstimate = varest(mod, fitted)

    ### — STEP 6: Get marginal confidence intervals
    print(marginalci(vₑstimate, lvl=0.95))
