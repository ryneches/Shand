from nose.tools import assert_equal, assert_almost_equal
from shand import quicktree
from dendropy import Tree
import numpy

test_tree = 'shand/tests/test.tree'
dpt = Tree.get( file=open(test_tree), schema='newick' )
for n,node in enumerate( dpt.inorder_node_iter() ) :
    node.label = n

def test_init() :
    T = quicktree.QuickTree( test_tree )
    assert_equal( type(T), quicktree.QuickTree )

def test_get_children() :
    T = quicktree.QuickTree( test_tree )
    for node in dpt.inorder_node_iter() :
        if not node.taxon :
            left, right = [ n.label for n in node.child_nodes() ]
        else :
            left, right = -1, -1
        L,R = T.get_children( node.label )
        assert_equal( L, left )
        assert_equal( R, right )

def test_get_distance_to_root() :
    T = quicktree.QuickTree( test_tree )
    for leaf in dpt.leaf_node_iter() :
        assert_almost_equal( T.get_distance_to_root( leaf.label ),
                             leaf.distance_from_root(),
                             places=4 )

def test_distance() :
    T = quicktree.QuickTree( test_tree )
    for leaf_a in dpt.leaf_node_iter() :
        for leaf_b in dpt.leaf_node_iter() :
            d1 = T.distance( leaf_a.label, leaf_b.label )
            d2 = abs( leaf_a.distance_from_root() 
                    - leaf_b.distance_from_root() )
            assert_almost_equal( d1, d2, places=4 )
   
def test_distances() :
    T = quicktree.QuickTree( test_tree )
    ids = []
    d1 = []
    for leaf_a in dpt.leaf_node_iter() :
        for leaf_b in dpt.leaf_node_iter() :
            d = abs( leaf_a.distance_from_root() 
                   - leaf_b.distance_from_root() )
            ids.append( [leaf_a.label, leaf_b.label] )
            d1.append( d )
    result = T.distances( numpy.array( ids, dtype=numpy.int64 ) )
    for D1,D2 in zip( d1,result ) :
        assert_almost_equal( D1, D2, places=3 )

#def test_dump_array() :
#    T = quicktree.QuickTree( test_tree )
#    assert T.dump_array()
