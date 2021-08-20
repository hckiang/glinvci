# Introduction

GLInvCI is a package that provides a framework for computing the maximum-likelihood estimates and asymptotic confidence intervals of a class of continuous-time Gaussian branching processes, including the Ornstein-Uhlenbeck branching process, which is commonly used in phylogenetic comparative methods. The framework is designed to be flexible enough that the user can easily specify their own parameterisation and obtain the maximum-likelihood estimates and confidence intervals of their own parameters.

The model in concern is ![\\mathscr{G}\_{Linv}](https://latex.codecogs.com/png.latex?%5Cmathscr%7BG%7D_%7BLinv%7D "\mathscr{G}_{Linv}") model, in which each species' traits evolve independently of each others after branching off from their common ancestor and for every non-root node. Let ![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\eta") is a child node of ![u](https://latex.codecogs.com/png.latex?u "u"), and ![z\_{\\eta}](https://latex.codecogs.com/png.latex?z_%7B%5Ceta%7D "z_{\eta}"), ![z\_{u}](https://latex.codecogs.com/png.latex?z_%7Bu%7D "z_{u}") denotes the multivariate traits. We assume that ![z\_{\\eta}\|z\_{u}](https://latex.codecogs.com/png.latex?z_%7B%5Ceta%7D%7Cz_%7Bu%7D "z_{\eta}|z_{u}") is a Gaussian distribution with expected value ![w\_{\\eta}+\\Phi\_{\\eta}z\_{u}](https://latex.codecogs.com/png.latex?w_%7B%5Ceta%7D%2B%5CPhi_%7B%5Ceta%7Dz_%7Bu%7D "w_{\eta}+\Phi_{\eta}z_{u}") and variance ![V\_{\\eta}](https://latex.codecogs.com/png.latex?V_%7B%5Ceta%7D "V_{\eta}"), where the matrices ![\\left(\\Phi\_{\\eta},w\_{\\eta},V\_{\\eta}\\right)](https://latex.codecogs.com/png.latex?%5Cleft%28%5CPhi_%7B%5Ceta%7D%2Cw_%7B%5Ceta%7D%2CV_%7B%5Ceta%7D%5Cright%29 "\left(\Phi_{\eta},w_{\eta},V_{\eta}\right)") are independent of ![z\_{\\eta}](https://latex.codecogs.com/png.latex?z_%7B%5Ceta%7D "z_{\eta}") but can depend on ![t\_{\\eta}](https://latex.codecogs.com/png.latex?t_%7B%5Ceta%7D "t_{\eta}"), ![t\_{u}](https://latex.codecogs.com/png.latex?t_%7Bu%7D "t_{u}") and other parameters. Note that ![z\_{\\eta}](https://latex.codecogs.com/png.latex?z_%7B%5Ceta%7D "z_{\eta}") and ![z\_{u}](https://latex.codecogs.com/png.latex?z_%7Bu%7D "z_{u}") can be of different dimension.

# Installation

The following command should install the latest version of the package:

    install.packages('devtools')
    devtools::install(
      'https://git.sr.ht/~hckiang/glinvci/blob/latest-tarballs/glinvci_latest_main.tar.gz')

# High-level and low-level interface

The package contains two levels of user interfaces. The high-level interface, accessible through the `glinv` function, provides facilities for handling missing traits, lost traits, multiple evolutionary regimes, and most importantly, the calculus chain rule. The lower-level interface, accessible through the `glinv_gauss` function, allows the users to operate purely in the ![\\left(\\Phi\_{\\eta},w\_{\\eta},V\_{\\eta}\\right)](https://latex.codecogs.com/png.latex?%5Cleft%28%5CPhi_%7B%5Ceta%7D%2Cw_%7B%5Ceta%7D%2CV_%7B%5Ceta%7D%5Cright%29 "\left(\Phi_{\eta},w_{\eta},V_{\eta}\right)") parameter space.

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

