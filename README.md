# Introduction

GLInvCI is a package that provides a framework for computing the maximum-likelihood estimates and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes, including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic comparative methods. The framework is designed to be flexible enough that the user can easily specify their own parameterisation and obtain the maximum-likelihood estimates and confidence intervals of their own parameters.

The model in concern is GLInv family, in which each species' traits evolve independently of each others after branching off from their common ancestor and for every non-root node. Let <img alt="$k$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/63bb9849783d01d91403bc9a5fea12a2.svg" align="middle" width="9.075367949999992pt" height="22.831056599999986pt"/> be a child node of <img alt="$j$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710416999999989pt" height="21.68300969999999pt"/>, and <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" align="middle" width="14.91068204999999pt" height="14.15524440000002pt"/>, <img alt="$z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/6fb0e03671d387b5d3116a9ee6e84ad3.svg" align="middle" width="13.74916289999999pt" height="14.15524440000002pt"/> denotes the corresponding multivariate traits. We assume that <img alt="$z_{k}|z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/689db23bebd4a7d6b48eb94ad677f530.svg" align="middle" width="34.04798429999999pt" height="24.65753399999998pt"/> is a Gaussian distribution with expected value <img alt="$w_{k}+\Phi_{k}z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/98abb1e31234a45dfd3839a728eae15b.svg" align="middle" width="73.65692894999998pt" height="22.465723500000017pt"/> and variance <img alt="$V_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/f75d10508c43d5007ef4f62dc4d36e47.svg" align="middle" width="16.855081649999992pt" height="22.465723500000017pt"/>, where the matrices <img alt="$\left(\Phi_{k},w_{k},V_{k}\right)$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/e20c8bad79d9b9bad3f20ca2b90dd947.svg" align="middle" width="84.8907675pt" height="24.65753399999998pt"/> are parameters independent of <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" align="middle" width="14.91068204999999pt" height="14.15524440000002pt"/> but can depend other parameters including <img alt="$t_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/c33bbdc9525803b09428bbe4e6223a67.svg" align="middle" width="13.20212684999999pt" height="20.221802699999984pt"/>. The traits <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" align="middle" width="14.91068204999999pt" height="14.15524440000002pt"/> and <img alt="$z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/6fb0e03671d387b5d3116a9ee6e84ad3.svg" align="middle" width="13.74916289999999pt" height="14.15524440000002pt"/> can have different number of dimension.

# Installation

The following command should install the latest version of the package:

    install.packages('devtools')
    devtools::install(
      'https://git.sr.ht/~hckiang/glinvci/blob/latest-tarballs/glinvci_latest_main.tar.gz')

# High-level and low-level interface

The package contains two levels of user interfaces. The high-level interface, accessible through the `glinv` function, provides facilities for handling missing traits, lost traits, multiple evolutionary regimes, and most importantly, the calculus chain rule. The lower-level interface, accessible through the `glinv_gauss` function, allows the users to operate purely in the <img alt="$\left(\Phi_{\eta},w_{\eta},V_{\eta}\right)$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/3b58684a4d580c22f8932bd74829bcc6.svg" align="middle" width="84.64997144999998pt" height="24.65753399999998pt"/> parameter space.

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

    ### --- STEP 3: Try getting the likelihood, gradient etc.
    H     = matrix(c(1,0,0,-1), k)
    theta = c(0,0)
    sig   = matrix(c(0.5,0,0,0.5), k)
    sig_x = t(chol(sig))
    par_init = c(H=diag(H),theta=theta,sig_x=sig_x[lower.tri(sig_x,diag=T)])
    print(par_init)
    print(lik(mod)(par_init))
    print(grad(mod)(par_init))
    print(hess(mod)(par_init))

    ### --- STEP 4: Fitting a model
    fitted = fit(mod, par_init)
    print(fitted)

    ### --- STEP 5: Estimating variance-covariance of the MLE
    v_estimate = varest(mod, fitted)

    ### --- STEP 6: Get marginal confidence intervals
    print(marginal_ci(v_estimate, lvl=0.95))
