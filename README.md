# Introduction

GLInvCI is a package that provides a framework for computing the maximum-likelihood estimates and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes, including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic comparative methods. The framework is designed to be flexible enough that the user can easily specify their own parameterisation and obtain the maximum-likelihood estimates and confidence intervals of their own parameters.

The model in concern is GLInv family, in which each species' traits evolve independently of each others after branching off from their common ancestor and for every non-root node. Let <img alt="$k$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/63bb9849783d01d91403bc9a5fea12a2.svg" width="17.2945773pt" height="11.4155283pt"/> be a child node of <img alt="$j$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/36b5afebdba34564d884d347484ac0c7.svg" width="15.929626349999998pt" height="14.0378568pt"/>, and <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" width="23.129891399999998pt" height="9.54335085pt"/>, <img alt="$z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/6fb0e03671d387b5d3116a9ee6e84ad3.svg" width="21.968372249999998pt" height="11.780795399999999pt"/> denotes the corresponding multivariate traits. We assume that <img alt="$z_{k}|z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/689db23bebd4a7d6b48eb94ad677f530.svg" width="42.267193649999996pt" height="17.031940199999998pt"/> is a Gaussian distribution with expected value <img alt="$w_{k}+\Phi_{k}z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/98abb1e31234a45dfd3839a728eae15b.svg" width="81.8761383pt" height="15.936036599999998pt"/> and variance <img alt="$V_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/f75d10508c43d5007ef4f62dc4d36e47.svg" width="25.074291pt" height="13.698590399999999pt"/>, where the matrices <img alt="$\left(\Phi_{k},w_{k},V_{k}\right)$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/e20c8bad79d9b9bad3f20ca2b90dd947.svg" width="93.10997685pt" height="16.438356pt"/> are parameters independent of <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" width="23.129891399999998pt" height="9.54335085pt"/> but can depend other parameters including <img alt="$t_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/c33bbdc9525803b09428bbe4e6223a67.svg" width="21.4213362pt" height="12.57663pt"/>. The traits <img alt="$z_{k}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/99cf206e5f6e47758950ec2c9da266d1.svg" width="23.129891399999998pt" height="9.54335085pt"/> and <img alt="$z_{j}$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/6fb0e03671d387b5d3116a9ee6e84ad3.svg" width="21.968372249999998pt" height="11.780795399999999pt"/> can have different number of dimension.

# Installation

The following command should install the latest version of the package:

    install.packages('devtools')
    devtools::install(
      'https://git.sr.ht/~hckiang/glinvci/blob/latest-tarballs/glinvci_latest_main.tar.gz')

# High-level and low-level interface

The package contains two levels of user interfaces. The high-level interface, accessible through the `glinv` function, provides facilities for handling missing traits, lost traits, multiple evolutionary regimes, and most importantly, the calculus chain rule. The lower-level interface, accessible through the `glinv_gauss` function, allows the users to operate purely in the <img alt="$\left(\Phi_{\eta},w_{\eta},V_{\eta}\right)$" src="https://git.sr.ht/~hckiang/glinvci/blob/readme_svgs/svgs/3b58684a4d580c22f8932bd74829bcc6.svg" width="92.8691808pt" height="17.031940199999998pt"/> parameter space.

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
