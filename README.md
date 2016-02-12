# Shand

A pipeline for investigating cospeciation in microbiomes

### Installation

Installing the pipeline is very simple :

    pip install shand

This will also install some other key packages, if they are not
already installed on your system :

* [`pandas`](http://pandas.pydata.org/)
* [`scikit-bio`](http://scikit-bio.org/)
* [`hat-trie`](https://github.com/kmike/hat-trie)
* [`screed`](https://github.com/ctb/screed)

You will also need to install a sequence aligner and a tree builder.
For now, Shand supports [`clustal-omega`](http://www.clustal.org/omega/)
for alignment and [`fasttree`](http://www.microbesonline.org/fasttree/)
for tree building. You can install these anywhere, as long as they can
be located using what's in your `PATH` variable. Before launching
Shand, you can check to make sure they're available by typing :

    $ which clustalo
    /usr/local/bin/clustalo
    $ which fasttree
    /usr/local/bin/fasttree

Shand will check for these when you launch it and report an error if
it can't find them.

### Preparing your data

Most sequencing projects are complicated, and so your metadata
probably isn't ready to just drop into a pipeline like
[`mothur`](http://www.mothur.org/), [`QIIME`](http://qiime.org/) or Shand
without meeting some criteria. Shand requires three pices of data :

* A `fasta` file of demultiplexed reads that have been scrubbed for chimeras
* A phylogeny of the host organisms in [`NEWICK`](http://scikit-bio.org/docs/latest/generated/skbio.io.format.newick.html#module-skbio.io.format.newick) format
* A CSV or TSV file that links sample names to the taxa in the host tree

There is more than one way to achieve each of these things, and so I
will outline in more detail what to do for each one.

#### Removing chimeras

I use [`vsearch`](https://github.com/torognes/vsearch) to perform
chimera checking. This step goes a lot faster if a reference database
of full-length 16S genes is used (I suggest the [`SILVA SSU
Ref`](http://www.arb-silva.de/download/arb-files/) database, but it
probably doesn't matter very much which one is used). The rest of the
analysis is reference-free.

Here is how I do it :

    wget http://ftp.arb-silva.de/release_123/Exports/SILVA_123_SSURef_tax_silva_trunc.fasta.gz
    
    vsearch --threads [n] --uchime_ref [your_reads.fasta.gz] \
    --db SILVA_123_SSURef_tax_silva_trunc.fasta.gz           \
    --chimeras [your_reads_chimeras.fasta]                   \
    --nonchimeras [your_reads_nochimeras.fasta]
    
    gzip [your_reads_chimeras.fasta]
    gzip [your_reads_nochimeras.fasta]

This takes a couple of days on my rather elderly 16 core machine (Xeon
E5520 @ 2.27GHz), circa 2010.

#### Getting a host phylogeny

Building good phylogenies of plants and animals can be very difficult.
My personal heuristic is that if the most recent common ancestor to
the members of the clade emerged closer to the present than to the 
[K-T event](https://en.wikipedia.org/wiki/Cretaceous%E2%80%93Paleogene_extinction_event),
I will use a tree from the literature rather than building my own. If
you are confident in your ability to build a reliable phylogeny of 
the hosts, then by all means, do so. Nevertheless, be mindful that
Shand will assume that the host phylogeny provided is correct. If 
you provide a lousy tree, it will give you lousy results.

To make it easier to use phylogenies from the literature, Shand does
not require you to prune your host tree. You must, however, take note
of how the taxa names are formatted.

#### Link the microbiome metadata to the host tree

Shand can use `QIIME`-style mapping files, though it ignores most of
the information contained in them. Suppose the reads are named like
this in the `fasta` file :

    >LOBOLAB1_50760

Shand will split the read name on the underscore (`'_'`). The string
on the left (`LOBOLAB1` in this example) should be a sample name in your
mapping file. Not all the sample names in your data need to appear in
your mapping file, but all sample names in your mapping file must
appear in your data. If your data contains reads from several
projects, there is no need to split the data into multiple `fasta`
files to use with Shand. It will take longer to index, but that only
happens once, and the index is saved to disk.

The mapping file must have two columns in it that link these sample
names to the taxa names that appear in your host tree. For example :

    #SampleID   Host
    LOBOLAB1    Lobochilotes_labiatus

By default, Shand will assume the first column contains the sample
names, the column labeled 'Host' contains the host name, and the
seperator is a tab. If you have a `QIIME`-style mapping file, all you
have to do is add a column named 'Host' and populate it with the host
taxa names. If not, you can override the sample name column label, 
the host name column label and the separator.

In more formal terms, you must create a many-to-one relationship from
the set of sample names to the set of host taxa.

### Running the pipeline

You can run Shand from within Python or from the command line. I
recommend running from withing Python using a [`jupyter`] notebook.
If you have to manipulate your data to get it ready to put into Shand,
a notebook will help keep a record of that.

#### Create problem plan 

    import shand
    myproject = shand.Problem() 

##### Attach the read data

First, attach the read data :

    myproject.add_reads( 'mydata.fasta' )

If the read names are separated from the read number by something
other than an underscore (`'_'`), you can specify a different
separator :

    myproject.add_reads( 'mydata.fasta.gz', read_name_sep='|' )

This will use `screed` to index the reads. A new file with the suffix
`_screed` will be created in the same directory as the reads, so you
will need to be able to write to this directory. If your data lives in
a directory to which you do not have write access, create a symbolic
link to your data in a directory where you can write, and use the path
to the link with Shand.

##### Attach the metadata

Next, attach the metadata. This could be a `QIIME` mapping file, or
just a CSV file that maps sample names to host taxa. If you have a
`QIIME`-like mapping file with tabs as separators and the host taxa
names in a column named 'Host' :

    myproject.add_metadata( 'mapping.txt' )    

If you have a separate column for the names that appear in the `fasta`
file containing your reads :

    myproject.add_metadata( 'mapping.txt', sample_id_col='ReadNames' )

If the column containing host taxa assignments is named something
other than 'Host', you can specify that :

    myproject.add_metadata( 'mapping.txt', host_col='Species' )

If you aren't using a tab delimited table, you can specify a different
separator :

    myproject.add_metadata( 'mapping.csv', sep=',' )

You may, of course, combine any of these options.

##### Attach the host tree

Next, attach the host tree :

    myproject.add_host_tree( 'host_tree.nwk' )

The host tree must be in `NEWICK` format, and the names must exactly
match those that appear in your metadata. The host tree may contain
taxa that are not in your sample list; the tree will be sheared to
include only the taxa that appear in your metadata. However, if your
sample list contains taxa that do not appear in your host tree, Shand
will report this as an error.

#### Run the pipeline

Now that all the data is attached, run the pipeline :

    myproject.run()

The following operations will then be carried out :

##### Identify unique reads

All unique sequences that appear more than once in your data will be
identified. This discards sequences that cannot be distinguished from
sequencing errors. If you want to change this threshold, run the
pipeline like so :

    myproject.run( cutoff=4 ) 

This will discard sequences that appear four or fewer times.

After unique sequnces are identified, four new `pandas` `DataFrames`
will be attached to the problem object. 

##### Align unique reads

##### Build phylogeny of reads

##### Compute Hommola correlation for all subtrees

### Post-analysis

#### Build JPrIME-DLTRS histories for candidate clades


### Name

Shand is named for James Shand, the inventor of the
[modern chainsaw](http://www.bac-lac.gc.ca/eng/discover/patents-1869-1919/Pages/item.aspx?IdNumber=186260&).

![Chainsaw](http://i.imgur.com/I1avJN5.jpg)
