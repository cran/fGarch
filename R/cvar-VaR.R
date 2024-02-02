## GNB TODO: this is to get thing going
qfun_fGarch <- function(dist, parameters){
    skew  <- parameters["skew"]
    shape <- parameters["shape"]
    switch(dist,
           norm = qnorm,
           snorm = function(p, ...) qsnorm(p, xi = skew, ...), # TODO: do we need '...'?
           ged = function(p, ...) qged(p, nu = shape, ...),
           sged = function(p, ...) qsged(p, nu = shape, xi = skew, ...),
           std = function(p, ...)  qstd(p, nu = shape, ...),
           sstd = function(p, ...) qsstd(p, nu = shape, xi = skew, ...),
           snig = function(p, ...) qsnig(p, zeta = shape, rho = skew, ...),
           ## default
           stop("distribution 'dist' not implemented here")
           )
}



VaR.fGARCH <- function(dist, p_loss = 0.05, ..., tol = .Machine$double.eps^0.5) {
    stopifnot(inherits(dist, "fGARCH"))
    
    mu_t <- dist@fitted
    sigma_t <- dist@sigma.t
    cond_dist <- dist@fit$params$cond.dist

    qf <- qfun_fGarch(cond_dist, coef(dist))

    cvar::VaR_qf(qf, p_loss, intercept = mu_t, slope = sigma_t, tol = tol)
}

ES.fGARCH <- function (dist, p_loss = 0.05, ...) {
    stopifnot(inherits(dist, "fGARCH"))
    
    mu_t <- dist@fitted
    sigma_t <- dist@sigma.t
    cond_dist <- dist@fit$params$cond.dist

    qf <- qfun_fGarch(cond_dist, coef(dist))

    cvar::ES(qf, p_loss, intercept = mu_t, slope = sigma_t)
}
