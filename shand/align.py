import subprocess
from os.path import splitext

def clustalo( fasta_file, threads=1 ) :
    
    base = splitext(fasta_file)[0]
    alignment_file = base + '_clustalo.fasta'
    log_file = base + '_clustalo.log'

    args = ['clustalo', '-v', '--force',
            '-i', fasta_file,
            '--outfmt=fasta',
            '--threads=' + str(threads),
            '-o', alignment_file,
            '-l', log_file ]

    subprocess.call( args )
    return alignment_file
