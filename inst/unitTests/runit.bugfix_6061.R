test.bugfix_6061 <- 
    function()
{       
    ## 2022-07-27 GB check fixes related to issue #6061, see NEWS and the source code
    stopifnot( all( diff(qsnorm(c(0.49, 0.5, 0.51))) > 0) )
    stopifnot( all( diff(psnorm(qsged(c(0.49, 0.5, 0.51)))) > 0) )
    
    stopifnot( all( diff(qsstd(c(0.49, 0.5, 0.51))) > 0) )
    stopifnot( all( diff(psstd(qsged(c(0.49, 0.5, 0.51)))) > 0) )
    
    stopifnot( all( diff(qsged(c(0.49, 0.5, 0.51))) > 0) )
    stopifnot( all( diff(psged(qsged(c(0.49, 0.5, 0.51)))) > 0) )

    # Return Value:
    return()    
}
