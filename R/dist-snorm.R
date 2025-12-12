
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA


################################################################################ 
# FUNCTION:              DESCRIPTION:
#  dsnorm                 Density for the skew normal Distribution
#  psnorm                 Probability function for the skew NORM
#  qsnorm                 Quantile function for the skew NORM
#  rsnorm                 Random Number Generator for the skew NORM
# FUNCTION:              DESCRIPTION:
#  .dsnorm                Internal, density for the skew normal Distribution
#  .psnorm                Internal, probability function for the skew NORM
#  .qsnorm                Internal, quantile function for the skew NORM
#  .rsnorm                Internal, random Number Generator for the skew NORM
################################################################################


dsnorm <- 
function(x, mean = 0, sd = 1, xi = 1.5, log = FALSE)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the density function of the skew normal distribution
    
    # Arguments:
    #   x - a numeric vector of quantiles.
    #   mean, sd, xi - location parameter, scale parameter, and 
    #       skewness parameter.
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 3) {
        xi = mean[3]
        sd = mean[2]
        mean = mean[1]
    } 
    
    # Shift and Scale:
    result = .dsnorm(x = (x-mean)/sd, xi = xi) / sd
    
    # Log:
    if(log) result = log(result)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


psnorm <- 
function(q, mean = 0, sd = 1, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the distribution function of the 
    #   skew normal distribution
    
    # Arguments:
    #   q - a numeric vector of quantiles.
    #   mean, sd, xi - location parameter, scale parameter, and 
    #       skewness parameter.
    
    # FUNCTION:
    
    # Shift and Scale:
    result = .psnorm(q = (q-mean)/sd, xi = xi)
          
    # Return Value:
    result
}


# ------------------------------------------------------------------------------    

   
qsnorm <- 
function(p, mean = 0, sd = 1, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the quantile function of the 
    #   skew normal distribution
    
    # Arguments:
    #   p - a numeric vector of probabilities.
    #   mean, sd, xi - location parameter, scale parameter, and 
    #       skewness parameter.
    
    # FUNCTION:
    
    # Shift and Scale:
    result = .qsnorm(p = p, xi = xi) * sd + mean
    
    # Return Value:
    result
}

    
# ------------------------------------------------------------------------------
 

rsnorm <- 
function(n, mean = 0, sd = 1, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Generate random deviates from the 
    #   skew normal distribution
    
    # Arguments:
    #   n - an integer value giving the number of observation.
    #   mean, sd, xi - location parameter, scale parameter, and 
    #       skewness parameter.
    
    # FUNCTION:
    
    # Shift and Scale:
    result = .rsnorm(n = n, xi = xi) * sd + mean
    
    # Return Value:
    result
}


################################################################################


.dsnorm <-  
function(x, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the density function of the "normalized" skew 
    #   normal distribution
    
    # FUNCTION:

    # Standardize:
    m1 = 2/sqrt(2*pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    z = x*sigma + mu  
    # Compute:
    Xi = xi^sign(z)
    g = 2 / (xi + 1/xi) 
    Density = g * dnorm(x = z/Xi)  
    # Return Value:
    Density * sigma 
}


# ------------------------------------------------------------------------------
   

.psnorm <- 
function(q, xi) 
{   
    # A function implemented by Diethelm Wuertz 
    ##
    ## fixed by GNB, see section 'CHANGES in fGarch VERSION 4021.87, 2022-08-06', subsection
    ## 'BUG fixes' in NEWS.Rd.

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Standardize:
      m1 = 2/sqrt(2*pi)
      mu = m1 * (xi - 1/xi)
      sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
      z = q*sigma + mu
    # Compute:  
      sig <- ifelse(z >= 0, 1, -1) # note: 1 for z = 0; was sign(z)
    Xi = xi^sig  # not sign(z)
    g = 2  / (xi + 1/xi)
    ## was: Probability = Heaviside(z) - sig * g * Xi * pnorm(q = -abs(z)/Xi)
      Probability = ifelse(z >= 0, 1, 0) - sig * g * Xi * pnorm(q = -abs(z)/Xi)
    # Return Value:
      Probability 
}
    
.qsnorm <- 
function(p, xi) 
{   
    ## A function implemented by Diethelm Wuertz
    ##
    ## Corrected at the centre part by Georgi N. Boshnakov on 2022-07-27 to fix
    ## bug [#6061], see section 'CHANGES in fGarch VERSION 4021.87, 2022-08-06',
    ## subsection 'BUG fixes' in NEWS.Rd.
    ##
    ## The old version of this and similar functions were temporarilly kept with
    ## suffix '_orig' but removed after the release of 4022.89, due to lack of
    ## (reported) problems with the fix. The examples below are kept for reference.
    ##
    ##
    ## Compare
    ##    plot(function(p) .qsnorm(p, xi = 1.5), from  = 0, to = 1)
    ##    plot(function(p) .qsnorm_orig(p, xi = 1.5), from  = 0, to = 1,
    ##         col = "blue", add = TRUE)
    ##
    ## Create a quantile function by numerically inverting psnorm:
    ##     f <- function(x, ...){ sapply(x, function(p)
    ##                                        gbutils::cdf2quantile(p, cdf = psnorm, ...))}
    ## It agrees with the fixed qsnorm (the 2nd is right over the first, 3rd is the orig.):
    ##   plot(function(p) .qsnorm(p, xi = 1.5), from = 0.49, to = 0.51)
    ##   plot(f, from = 0.49, to = 0.51, add = TRUE, col = "red")
    ##   plot(function(p) .qsnorm_orig(p, xi = 1.5), from = 0.49, to = 0.51, add = TRUE)
    ##   
    ##   plot(function(p) .qsnorm(p, xi = 1.5), from = 0, to = 1)
    ##   plot(f, from = 0, to = 1, add = TRUE, col = "red")
    ##   plot(function(p) .qsnorm_orig(p, xi = 1.5), from = 0, to = 1, add = TRUE)

    ## Description:
    ##   Internal Function
    
    ## FUNCTION:
    
    ## Standardize:
    m1 = 2/sqrt(2*pi)
    mu = m1 * (xi - 1/xi)
    sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    ## Compute:  
    g = 2 / (xi + 1/xi)
    pxi <- p - (1 / (1 + xi^2)) # not p - 1/2
    sig <- sign(pxi)  # not p - 1/2
    Xi = xi^sig
    p = (Heaviside(pxi) - sig * p) / (g * Xi)  # pxi, not p - 1/2
    ## GNB:
    ## Fixed: BUG [#6061] fGarch::qsnorm() incorrect around p=0.5
    ##
    ## This has, in general,  discontinuity for p = 1/2, since then sig = 0.
    ## Note that the original p = 1/2 is transformed above to 1/(2*g*Xi),
    ## so qnorm() doesn't necessarilly give 1/2 when p = 1/2.
    ##
    ## Note also that p can be a vector.
    ##
    ## Further note: the issue at p = 0.5 is a separate problem. Zooming in
    ## shows that the quantile is not continuous at p = 0.5 and to the right of
    ## 0.5 the values are smaller than just to the left of 0.5 up to around 0.51.
    ##
    ## SOLUTION(?): The error seems to be that sign() and Heaviside should compare to
    ##    1/(1+1/xi^2), not 0.5 which is correct only for xi = 1. 
    ## 
    Quantile = (-sig * qnorm(p = p, sd = Xi) - mu ) / sigma

    ## Return Value:
    Quantile 
}
    
 
# ------------------------------------------------------------------------------
  

.rsnorm <- 
function(n, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Generate Random Deviates:
      weight = xi / (xi + 1/xi)
      z = runif(n, -weight, 1-weight)
      Xi = xi^sign(z)
      Random = -abs(rnorm(n))/Xi * sign(z)  
    # Scale:
      m1 = 2/sqrt(2*pi)
      mu = m1 * (xi - 1/xi)
      sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
      Random = (Random - mu ) / sigma   
    # Return value:
      Random 
}
 

################################################################################

          
