# https://link.springer.com/article/10.3758/s13428-012-0289-7
import numpy as np
import math
import scipy
from scipy import stats


def compare_correlations(r1, r2, n1, n2):
    '''
    Finds a p-value for the significance of the difference between two correlations. Uses Fisher-Z transformation. See https://link.springer.com/article/10.3758/s13428-012-0289-7 for more details on this method.
    
    Parameters:
    r1 (float): pearson correlation for group 1
    r2 (float): pearson correlation for group 2
    n1 (int): size of group 1. Must be >= 4, otherwise will return Nan
    n2 (int): size of group 2. Must be >= 4, otherwise will return Nan
    
    Returns:
    float: p-value for difference in correlation
    
    '''
    rp1 = np.arctanh(r1)
    rp2 = np.arctanh(r2)    
    if n1 < 4 or n2 < 4:
        return(np.nan)
    Sr12 = math.sqrt((1/(n1-3))+(1/(n2-3)))
    z = (rp1-rp2) / Sr12
    p = scipy.stats.norm.sf(abs(z))*2
    return (p)
 