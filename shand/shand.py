import pandas as pd
from screed import read_fasta_sequences, ScreedDB
from hat_trie import Trie
import pyprind
import skbio
from os.path import splitext, exists
from align import clustalo
from tree import fasttree
import stats
from multiprocessing import Process, Queue, current_process

class Problem(object) :

    def __init__( self, name, threads=1 ) :
        self.name = name
        self.threads = threads

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
            self.sample_ids = list(df[sample_id_col])
            self.sample_id_col = sample_id_col
            self.metadata.index = self.metadata[ sample_id_col ]
        else :
            self.sample_ids = list(df.index)
            self.sample_id_col = df.index.name
        
        # fail if there are NEWICK reserved characters in sample names
        newick_reserved = set( [ '[', ']', '(', ')', ',', ';', ':', ' ', '\t' ] )
        newick_clash = reduce( lambda a,b : a|b, map( set, self.sample_ids ) ) & newick_reserved
        if newick_clash :
            raise Exception('sample IDs contain reserved characters : ' + str(newick_clash) )
        
    def add_host_tree( self, host_tree_file ) :
        tree = skbio.tree.TreeNode.read(host_tree_file)
        # fail if there are missing taxa in the host tree
        leftovers = set(self.metadata[self.host_col]) - set([ tip.name for tip in tree.tips() ])
        if not leftovers :
            tree = tree.shear( list( set( self.metadata[self.host_col] ) ) )
            self.host_tree = tree
            self.host_tree_dmatrix = tree.tip_tip_distances()
        else :
            raise Exception('metadata contains species not found in host tree : ' + ', '.join(leftovers))
        
    def find_unique_reads( self, cutoff ) :
        bar_title = 'building trie...'
        self.trie = Trie()
        p = pyprind.ProgBar( len( self.db ), monitor=True, title=bar_title )
        for n,name in enumerate( self.db ) :
            demul_name, readnumber = name.split( self.read_name_sep )
            if demul_name in self.sample_ids :
                seq = unicode( self.db[name].sequence )
                if not self.trie.__contains__( seq ) : self.trie[seq] = []
                self.trie[seq].append( name ) 
            p.update()
        print(p)
        basename = self.name + '_unique_' + str(cutoff)
        self.unique_seq_file = basename + '.fasta'
        self.unique_seq_to_sample_file = basename + '.txt'
        bar_title = 'writing uniqued records with at least ' + str(cutoff) + ' instances...'
        p = pyprind.ProgBar( len(self.trie.keys()), monitor=True, title=bar_title )
        with open( self.unique_seq_file,           'w' ) as f1, \
             open( self.unique_seq_to_sample_file, 'w' ) as f2 :
            for n,seq in enumerate( self.trie.keys() ) :
                records = self.trie[seq]
                if len(records) >= cutoff :
                    f1.write( '>' + records[0] + '\n' + seq + '\n' )
                    f2.write( ','.join(records) + '\n' )
                p.update()
        print(p)

    def build_count_tables( self, cutoff ) :
        bar_title = 'bulding count table...'
        p = pyprind.ProgBar( len( self.trie.keys() ), monitor=True, title=bar_title )
        counts = {}
        for n,record in enumerate( self.trie.keys() ) :
            OTUs = self.trie[record]
            if not len(OTUs) >= cutoff : continue
            counts[ OTUs[0] ] = map( lambda x : map( lambda x : x.split(self.read_name_sep)[0], OTUs ).count(x), self.sample_ids )    
            p.update()
        print(p)
        self.count_table = pd.DataFrame( counts, index=self.sample_ids )
        self.count_table.index.name = self.sample_id_col
        self.abundance_table = self.count_table.div( self.count_table.sum( axis=1 ), axis=0 )
        # Take the OTU counts for host taxa with more than one 
        # sample, and merge them (basically, and inner join)
        self.host_count_table = self.count_table.join( 
                                    self.metadata[ self.host_col ] ).groupby( 
                                        self.metadata[ self.host_col] ).sum()
         
        self.host_abundance_table = self.host_count_table.div( self.host_count_table.sum(axis=1), axis=0)
        
        # save tables
        self.count_table.to_csv( self.name + '_count_table.tsv', sep='\t' )
        self.abundance_table.to_csv( self.name + '_abundance_table.tsv', sep='\t' )
        self.host_count_table.to_csv( self.name + '_host_count_table.tsv', sep='\t' )
        self.host_abundance_table.to_csv( self.name + '_host_abundance_table.tsv', sep='\t')

    def build_guest_tree( self ) :
        # build alignment
        print 'building alignment...'
        self.alignment_file = clustalo( self.unique_seq_file, threads=self.threads )
        
        # build tree
        print 'bulding guest tree...'
        self.guest_tree_file = fasttree( self.alignment_file, threads=self.threads )
        
        # load guest tree
        print 'loading guest tree...'
        self.guest_tree = skbio.tree.TreeNode.read( self.guest_tree_file, 
                                                    convert_underscores=False )
        self.guest_tree.assign_ids()
        self.guest_tree.index_tree()
         
    def predict_cospeciation( self, max_tree_size ) :
                
        def worker( work_queue, done_queue ) :
            for task in iter( work_queue.get, 'STOP' ) :
                h  = task['host_dmatrix']
                ct = task['clade_tree']
                l  = task['links']
                M  = task['permutations']
                c  = ct.tip_tip_distances()
                try :
                    t = stats.all_tests( h, c, l, permutations=M )
                    t['pid'] = current_process().name
                    t['node_id'] = task['node_id']
                    t['n_links'] = task['n_links']
                    t['clade_size'] = task['clade_size']
                    done_queue.put(t)
                except AssertionError :
                    done_queue.put(False)
            return True
 
        work_queue = Queue()
        done_queue = Queue()
        processes = []
        
        internal_nodes = len( list( self.guest_tree.non_tips() ) )
        bar_title = 'building work queue...'
        progbar = pyprind.ProgBar( internal_nodes, monitor=True, title=bar_title )
        for node in self.guest_tree.non_tips() :
            progbar.update()
            clade = node.copy()
            clade.index_tree()
            clade_leafs = [ tip.name for tip in clade.tips() ]
            clade_size = len(clade_leafs)
            if clade_size <= 3 : continue
            if clade_size >= max_tree_size : continue
            links = self.host_count_table[ clade_leafs ]
            n_links = ( links.values > 0 ).sum()
            if n_links <= 3 : continue
            task = { 'host_dmatrix' : self.host_tree_dmatrix,
                     'clade_tree'   : clade,
                     'links'        : links,
                     'permutations' : self.permutations,
                     'node_id'      : node.id,
                     'n_links'      : n_links,
                     'clade_size'   : clade_size }
            work_queue.put(task)
        
        print 'creating worker threads...'
        for w in xrange( self.threads ) :
            work_queue.put( 'STOP' )
            p = Process( target = worker, args = ( work_queue, done_queue ) )
            p.start()
            processes.append( p )
        
        print 'launching ' + str(self.threads) + ' threads for ' + str(work_queue.qsize()) + ' tasks...'
        for p in processes :
            p.join()

        done_queue.put( 'STOP' )
        
        n_results = done_queue.qsize() - 1
        bar_title = 'writing results...'
        progbar = pyprind.ProgBar( n_results, monitor=True, title=bar_title )
        with open( self.name + '_cospeciation_results_table.tsv', 'w' ) as f :
            cols = [ 'node_id', 'pid', 'n_links', 'clade_size', 'r',
                     'p_r', 'roh', 'p_roh', 'tau', 'p_tau' ]
            f.write( '\t'.join( cols ) + '\n' )
            for task in iter( done_queue.get, 'STOP' ) :
                if not task : continue
                result = []
                for item in cols :
                    result.append( str(task[item]) )
                f.write( '\t'.join( result ) + '\n' )
                progbar.update()
    
    def run( self, cutoff=2, permutations=10, max_tree_scale=0.1 ) :
        self.permutations = permutations
        self.find_unique_reads( cutoff )
        self.build_count_tables( cutoff ) 
        self.build_guest_tree()
        max_tree_size = len(self.guest_tree.subset()) * max_tree_scale
        self.predict_cospeciation( max_tree_size )        
        print '\nrun complete.'
