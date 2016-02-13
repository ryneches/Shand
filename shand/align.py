import subprocess
from os.path import splitext

def clustalo( fasta_file, threads=1 ) :
    
    base = splitext(fasta_file)
    alignment_file = base + '_clustalo.fasta'
    log_file = '_clustalo.log'

    args = ['clustalo', '-vvv',
            '-i', fasta_file,
            '--threads=' + str(threads),
            '-o', alignment_file,
            '-l', log_file ]

    return subprocess.check_output( args )
