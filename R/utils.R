## Invert a Hessian matrix and warn if it is not symmetric positive definite
hessian_inv_warn = function (X, matrix_name = 'The matrix', extra_info_nonpd = '', extra_info_noninvt = '') {
  tryCatch({
    L = chol(X)
    chol2inv(L)
  }, error = function (e_) {
    tryCatch({
      Y = solve(X)
      warning(sprintf('%s is numerically non-positive-definite. %s', matrix_name, extra_info_nonpd))
      Y
    }, error = function (f_) {
      stop(sprintf('%s is numerically non-invertible. %s', matrix_name, extra_info_noninvt))
    })
  })
}


list_set_default = function (L, defaults) {
  default_names = names(defaults)
  L_names       = names(L)
  for (i in seq_along(defaults)) {
    if (!(default_names[[i]] %in% L_names)) {
      L[[default_names[[i]]]] = defaults[[i]]
    }
  }
  L
}
