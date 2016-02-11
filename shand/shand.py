import pandas as pd
from screed import read_fasta_sequences, ScreedDB
from hat_trie import Trie
import ProgressBar
import skbio
from os.path import splitext, exists

class Problem(object) :
    def __init__( self ) :
        pass
    def add_reads( self, reads, read_name_sep='_' ) :
        if exists( reads + '_screed' ) :
            print 'reads previously indexed.'
        else :
            print 'indexing records...'
            read_fasta_sequences(reads)
        print 'building database...'
        db = ScreedDB(reads)
        self.db = db
        self.read_name_sep = read_name_sep
        self.reads_path = reads
    def add_metadata( self, metadata_file, 
                      sample_id_col=None, host_col='Host', 
                      sep='\t', drop_cols=None ) :
        '''
        One of the columns must have host names which match the names
        found in the host tree. By default it's 'Host', but you can
        change it.
        '''
        df = pd.DataFrame.from_csv( metadata_file, sep=sep )
        if drop_cols :
            for col in drop_cols :
                del df[col]
        self.metadata = df
        self.host_col = host_col
        if sample_id_col :
            self.sample_ids = df[sample_id_col]
        else :
            self.sample_ids = df.index
    def add_host_tree( self, host_tree_file, shear=True ) :
        tree = skbio.tree.TreeNode.read(host_tree_file)
        if shear :
            tree = tree.shear( list( set( self.metadata[self.host_col] ) ) )
        self.host_tree = tree
    def run( self, cutoff=2 ) :
        print 'building trie...\n'
        self.trie = Trie()
        p = ProgressBar.ProgressBar( len( self.db ) )
        for n,name in enumerate( self.db ) :
            if n % 10000 == 0 : p.animate( n + 1 )
            demul_name, readnumber = name.split( self.read_name_sep )
            if demul_name in self.sample_ids :
                seq = unicode( self.db[name].sequence )
                if not self.trie.__contains__( seq ) : self.trie[seq] = []
                self.trie[seq].append( name ) 
        print 'writing uniqued records with at least ' + str(cutoff) + ' instances...\n'
        p = ProgressBar.ProgressBar(len(self.trie.keys()))
        basename = splitext( self.reads_path )[0] + '_unique_' + str(cutoff)
        with open( basename + '.fasta', 'w' ) as f1, \
             open( basename + '.txt',   'w' ) as f2 :
            for n,seq in enumerate( self.trie.keys() ) :
                if n % 10000 == 0 : p.animate( n + 1 )
                records = self.trie[seq]
                if len(records) >= cutoff :
                    f1.write( '>' + records[0] + '\n' + seq + '\n' )
                    f2.write( ','.join(records) + '\n' )
        print 'bulding count table...\n'
        p = ProgressBar.ProgressBar( len( self.trie.keys() ) )
        counts = {}
        for n,record in enumerate( self.trie.keys() ) :
            if n % 10000 == 0 : p.animate( n + 1 )
            OTUs = self.trie[record]
            counts[ OTUs[0] ] = map( lambda x : map( lambda x : x.split(self.read_name_sep)[0], OTUs ).count(x), self.sample_ids )    
        self.count_table = pd.DataFrame( counts, index=self.sample_ids )
        self.abundance_table = self.count_table.div( self.count_table.sum( axis=1 ), axis=0 )
        # Take the OTU counts for host taxa with more than one 
        # sample, and merge them (basically, and inner join)
        self.host_count_table = pd.merge( self.count_table,
                                          pd.DataFrame( self.metadata[self.host_col] ),
                                          right_index=True,
                                          left_index=True ).groupby( self.host_col ).sum()
        self.host_abundance_table = self.host_count_table.div( self.host_count_table.sum(axis=1), axis=0)
