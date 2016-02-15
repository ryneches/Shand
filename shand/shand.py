import pandas as pd
from screed import read_fasta_sequences, ScreedDB
from hat_trie import Trie
import ProgressBar
import skbio
from os.path import splitext, exists
from align import clustalo
from tree import fasttree

class Problem(object) :
    def __init__( self, threads=1 ) :
        self.thrads = threads
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
        
        # fail if there are NEWICK reserved characters in sample names
        newick_reserved = set( [ '[', ']', '(', ')', ',', ';', ':', ' ', '\t' ] )
        newick_clash = reduce( lambda a,b : a|b, map( set, self.sample_ids ) ) & newick_reserved
        if newick_clash :
            raise Exception('sample IDs contain reserved characters : ' + ' '.join( newick_clash ) )
        
    def add_host_tree( self, host_tree_file ) :
        tree = skbio.tree.TreeNode.read(host_tree_file)
        # fail if there are missing taxa in the host tree
        leftovers = set(self.metadata[self.host_col]) - set([ tip.name for tip in tree.tips() ])
        if not leftovers :
            tree = tree.shear( list( set( self.metadata[self.host_col] ) ) )
            self.host_tree = tree
        else :
            raise Exception('metadata contains species not found in host tree : ' + ', '.join(leftovers))
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
        
        self.unique_seq_file = basename + '.fasta'
        self.unique_seq_to_sample_file = basename + '.txt'
        with open( self.unique_seq_file,           'w' ) as f1, \
             open( self.unique_seq_to_sample_file, 'w' ) as f2 :
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
            if not len(OTUs) >= cutoff : continue
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
        
        # build alignment
        print 'building alignment...'
        self.alignment_file = clustalo( self.unique_seq_file, threads=self.threads )
        
        # build tree
        print 'bulding guest tree...'
        self.guest_tree_file = fasttree( self.alignment_file, threads=self.threads )
        
        # load guest tree
        print 'loading guest tree...'
        self.guest_tree = skbio.tree.NodeTree.read( self.guest_tree_file )
        print 'computing patristic distances...'
        self.guest_tree_dmatrix = self.guest_tree.tip_tip_distances()


