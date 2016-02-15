import subprocess
from os.path import splitext
from os import environ

def fasttree( alignment_file, threads=1 ) :
    
    base = splitext(alignment_file)[0]
    tree_file = base + '_fasttree.tree'
    log_file = base + '_fasttree.log'
    
    e = environ
    e['OMP_NUM_THREADS'] = str(threads)
    
    args = ['FastTreeMP', '-nt', '-gtr', 
            '-log', log_file ] 

    with open( tree_file, 'w' ) as outfile:
        with open( alignment_file, 'r' ) as infile:
            proc = subprocess.Popen( args, env=e, stdout=outfile, stdin=infile )
            proc.wait()

    return tree_file
