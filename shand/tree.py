import subprocess
from os.path import splitext

def fasttree( alignment_file, threads=1 ) :
    
    base = splitext(alignment_file)[0]
    tree_file = base + '_fasttree.tree'
    log_file = base + '_fasttree.log'

    args = ['FastTreeMp', '-nt', '-gtr', 
            '-log', log_file ] 

    with open( tree_file, 'w' ) as outfile:
        with open( alignment_file, 'r' ) as infile:
            proc = subprocess.Popen( args, stdout=outfile, stdin=infile )
            proc.wait()
