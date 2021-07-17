
    ## THESE TESTS ARE GENERATED AUTOMATICALLY. DO NOT MODIFY.
    test_that("With restrictions", {


        gen_data = function (seed, treesize, k, missing_ranges, lost_ranges) 
{
    RNGkind(sample.kind = "Rejection")
    set.seed(seed)
    tr = ape::rtree(treesize)
    x = matrix(rnorm(k * treesize), k, treesize)
    x0 = rnorm(k)
    H = matrix(rnorm(k * k), k)
    H[lower.tri(H, diag = F)] = H[lower.tri(H, diag = F)]/5
    H[upper.tri(H, diag = F)] = H[upper.tri(H, diag = F)]/5
    theta = rnorm(k)
    sig_x_frs = matrix(0, k, k)
    sig_x_frs[upper.tri(sig_x_frs, diag = T)] = rnorm(k * (k + 
        1)/2)
    diag(sig_x_frs) = sqrt(abs(diag(sig_x_frs))) + 1
    sig = sig_x_frs %*% t(sig_x_frs)
    tmp = t(chol(sig))
    diag(tmp) = log(diag(tmp))
    sig_x = tmp[lower.tri(sig, diag = T)]
    if (!is.null(missing_ranges)) 
        for (rgidx in seq_along(missing_ranges)) x[rgidx, missing_ranges[[rgidx]]] = NA
    if (!is.null(lost_ranges)) 
        for (rgidx in seq_along(lost_ranges)) {
            x[rgidx, lost_ranges[[rgidx]]] = NaN
        }
    list(tr = tr, k = k, x = x, x0 = x0, H = H, theta = theta, 
        sig_x = sig_x, sig_x_frs = sig_x_frs, sig = sig)
}
        par = list(list(treesize = 2, k = 1, miss = NULL, lost = NULL), list(
    treesize = 4, k = 1, miss = NULL, lost = NULL), list(treesize = 29, 
    k = 2, miss = list(9, 10), lost = list(1:7, 19:20)))
        repar = list(list(fn = "ou_diagH(ou_zaplost(oupar))", jac = "dou_diagH(dou_zaplost(oujac))", 
    hess = "hou_diagH(hou_zaplost(ouhess))", npar = "nparams_ou_diagH(D$k)", 
    parform = "c(diag(D$H),D$theta,D$sig_x)"), list(fn = "ou_diagH(ou_haltlost(oupar))", 
    jac = "dou_diagH(dou_haltlost(oujac))", hess = "hou_diagH(hou_haltlost(ouhess))", 
    npar = "nparams_ou_diagH(D$k)", parform = "c(diag(D$H),D$theta,D$sig_x)"), 
    list(fn = "ou_logdiagH(ou_zaplost(oupar))", jac = "dou_logdiagH(dou_zaplost(oujac))", 
        hess = "hou_logdiagH(hou_zaplost(ouhess))", npar = "nparams_ou_logdiagH(D$k)", 
        parform = "c(diag(D$H),D$theta,D$sig_x)"), list(fn = "ou_logdiagH(ou_haltlost(oupar))", 
        jac = "dou_logdiagH(dou_haltlost(oujac))", hess = "hou_logdiagH(hou_haltlost(ouhess))", 
        npar = "nparams_ou_logdiagH(D$k)", parform = "c(diag(D$H),D$theta,D$sig_x)"), 
    list(fn = "ou_symH(ou_zaplost(oupar))", jac = "dou_symH(dou_zaplost(oujac))", 
        hess = "hou_symH(hou_zaplost(ouhess))", npar = "nparams_ou_symH(D$k)", 
        parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"), 
    list(fn = "ou_symH(ou_haltlost(oupar))", jac = "dou_symH(dou_haltlost(oujac))", 
        hess = "hou_symH(hou_haltlost(ouhess))", npar = "nparams_ou_symH(D$k)", 
        parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"), 
    list(fn = "ou_spdH(ou_zaplost(oupar))", jac = "dou_spdH(dou_zaplost(oujac))", 
        hess = "hou_spdH(hou_zaplost(ouhess), dou_zaplost(oujac))", 
        npar = "nparams_ou_spdH(D$k)", parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"), 
    list(fn = "ou_spdH(ou_haltlost(oupar))", jac = "dou_spdH(dou_haltlost(oujac))", 
        hess = "hou_spdH(hou_haltlost(ouhess), dou_haltlost(oujac))", 
        npar = "nparams_ou_spdH(D$k)", parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"), 
    list(fn = "ou_spdH(ou_zaplost(oupar),log=F)", jac = "dou_spdH(dou_zaplost(oujac),log=F)", 
        hess = "hou_spdH(hou_zaplost(ouhess), dou_zaplost(oujac),log=F)", 
        npar = "nparams_ou_spdH(D$k)", parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"), 
    list(fn = "ou_spdH(ou_haltlost(oupar),log=F)", jac = "dou_spdH(dou_haltlost(oujac),log=F)", 
        hess = "hou_spdH(hou_haltlost(ouhess), dou_haltlost(oujac),log=F)", 
        npar = "nparams_ou_spdH(D$k)", parform = "c(D$H[lower.tri(D$H,diag=T)],D$theta,D$sig_x)"))
        i = 1
        D = gen_data(i*314, par[[i]]$treesize, par[[i]]$k,
                     par[[i]]$miss, par[[i]]$lost)
            j = 1
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7657127690)
        
            j = 2
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7657127690)
        
            j = 3
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -4.6109789708)
        
            j = 4
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -4.6109789708)
        
            j = 5
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7657127690)
        
            j = 6
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7657127690)
        
            j = 7
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -6.1075126139)
        
            j = 8
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -6.1075126139)
        
            j = 9
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -0.8864052261)
        
            j = 10
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -0.8864052261)
        
        i = 2
        D = gen_data(i*314, par[[i]]$treesize, par[[i]]$k,
                     par[[i]]$miss, par[[i]]$lost)
            j = 1
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -5.6146244595)
        
            j = 2
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -5.6146244595)
        
            j = 3
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7723015968)
        
            j = 4
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -15.7723015968)
        
            j = 5
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -5.6146244595)
        
            j = 6
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -5.6146244595)
        
            j = 7
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -14.2432881326)
        
            j = 8
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -14.2432881326)
        
            j = 9
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -3.0712190930)
        
            j = 10
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -3.0712190930)
        
        i = 3
        D = gen_data(i*314, par[[i]]$treesize, par[[i]]$k,
                     par[[i]]$miss, par[[i]]$lost)
            j = 1
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -64.2839302123)
        
            j = 2
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -64.2839302123)
        
            j = 3
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -87.8881551963)
        
            j = 4
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -87.8881551963)
        
            j = 5
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -74.3482596935)
        
            j = 6
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -81.3334552972)
        
            j = 7
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -118.9767257510)
        
            j = 8
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -124.8499942739)
        
            j = 9
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -201.9112320177)
        
            j = 10
            mod = glinv(D$tr, D$x0, D$x,
                    parfns = list(eval(parse(text=repar[[j]]$fn))),
                    pardims = list(eval(parse(text=repar[[j]]$npar))),
                    parjacs = list(eval(parse(text=repar[[j]]$jac))),
                    parhess = list(eval(parse(text=repar[[j]]$hess))))
            varest = suppressWarnings(varest(mod, eval(parse(text=repar[[j]]$parform)),
                                      method="analytical"))
            expect_equal(sum(varest$hessian), -200.0100172472)
        
    })
