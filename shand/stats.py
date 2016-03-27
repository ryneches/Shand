import pandas
import skbio
from scipy.stats import pearsonr
from itertools import combinations
import random
from numpy import array

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

    F = foundon( links )

    assert len(F) > 3
    
    a,b = [],[]
    for (i,j),(k,l) in combinations( F, 2 ) :
        if shuffled :
            i = random.choice( guest_dmatrix.index )
            j = random.choice( host_dmatrix.index )
            k = random.choice( guest_dmatrix.index )
            l = random.choice( host_dmatrix.index )
        a.append( host_dmatrix[ j, l ] )
        b.append( guest_dmatrix[ i, k ] )
    return ( a, b )

def hommola( host_dmatrix,
             guest_dmatrix,
             links,
             permutations=100 ) :

    a,b = host_guest_distances( host_dmatrix, guest_dmatrix, links )
    pcc, p = pearsonr( a, b )
    
    prm = []
    for n in range(len(permutations)) :
        a,b = host_guest_distances( host_dmatrix,
                                    guest_dmatrix,
                                    links,
                                    shuffle=False )
        prm.append( pearsonr( a, b )[0] )
    
    p_value = ( ( array(prm) >= pcc ).sum() + 1 ) / ( permutations + 1 )
    
    return pcc, p_value, prm
