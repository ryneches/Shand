import pandas
import skbio
from scipy.stats import pearsonr
from itertools import combinations
import random
from numpy import array, zeros

def foundon( df ) :
    s = df.unstack()
    return s[s > 0].to_dict().keys()

def host_guest_distances( host_dmatrix, 
                          guest_dmatrix,
                          links,
                          shuffled=False ) :

    assert type( host_dmatrix ) == skbio.stats.distance._base.DistanceMatrix
    assert type( guest_dmatrix ) == skbio.stats.distance._base.DistanceMatrix
    assert type( links ) == pandas.core.frame.DataFrame
    
    assert host_dmatrix.shape[0] > 3
    assert guest_dmatrix.shape[0] > 3

    assert set( host_dmatrix.ids ) == set( links.index )
    assert set( guest_dmatrix.ids ) == set( links.columns )

    nlinks = ( links.values > 0 ).sum()
    a = zeros( nlinks * ( nlinks - 1 ) / 2 )
    b = zeros( nlinks * ( nlinks - 1 ) / 2 )
    
    if not shuffled :
        F = foundon( links )
        
        assert len(F) > 3
        
        for n,((i,j),(k,l)) in enumerate( combinations( F, 2 ) ) :
            a[n] = host_dmatrix[ j, l ]
            b[n] = guest_dmatrix[ i, k ]
        
    else :
        
        for n in range( nlinks ) :
            i = random.choice( guest_dmatrix.ids )
            j = random.choice( host_dmatrix.ids )
            k = random.choice( guest_dmatrix.ids )
            l = random.choice( host_dmatrix.ids )
            a[n] = host_dmatrix[ j, l ]
            b[n] = guest_dmatrix[ i, k ]
 
    return ( a, b )

def hommola( host_dmatrix,
             guest_dmatrix,
             links,
             permutations=100 ) :

    a,b = host_guest_distances( host_dmatrix, guest_dmatrix, links )
    pcc, p = pearsonr( a, b )
    
    prm = zeros( permutations )
    for n in range(permutations) :
        a,b = host_guest_distances( host_dmatrix,
                                    guest_dmatrix,
                                    links,
                                    shuffled=True )
        prm[n] = ( pearsonr( a, b )[0] )

    p_value = ( ( array(prm) >= pcc ).sum() + 1 ) / float( permutations + 1 )
    
    return pcc, p_value, prm
