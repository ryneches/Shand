from shand import quicktree
import skbio

test_tree = skbio.TreeNode.read( 'tests/test.tree' )

def test_quicktree_init() :
    T = quicktree.QuickTree( test_tree )
    assert type(T) is QuickTree
