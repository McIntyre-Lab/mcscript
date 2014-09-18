import gffutils
import os.path
from collections import defaultdict

class _Anno(object):
    """ General class for handling annotation database created from GFFUtils 
    
    Arguments:
    filename (str) = filename of your gff/gtf file
    filetype {'gff', 'gtf'} = file extension
    force (bool) = if True force the recreation of the database
    merge {'merge', 'create_unique', 'error', 'warning'} = passed to the
    merge_strategey option in gffutils.create_db()

    Attributes:
    db = GFFUtils database
    _featureKey = a dictionary to convert between feature keywords in the GFFUtils database
    """

    def __init__(self, filename, filetype, force=False, merge=''):
        # Create database name
        if os.path.splitext(filename)[1] == '.' + filetype:
            dbname = filename + '.db'
        else:
            print "Please provide a {0} file.".format(filetype)
            raise ValueError

        # Create or connect to GFFUtils database
        if not force:
            try:
                self.db = gffutils.FeatureDB(dbname)
            except:
                self.db = gffutils.create_db(filename, dbname, force=force, merge_strategy=merge)
        else:
                self.db = gffutils.create_db(filename, dbname, force=force, merge_strategy=merge)

        # Create a general feature key
        self._featureKey = { 'mRNA': 'mRNA',
                            'exon': 'exon',
                            'intron': 'intron',
                            'CDS': 'CDS',
                            '3UTR': '3UTR',
                            '5UTR': '5UTR'}

    def _get_annotation(self, name, childtype=''):
        """ Method that returns a list of children given some parent """
        try:
            return [i for i in self.db.children(name, featuretype=childtype, order_by='start')]
        except:
            print "Sorry I could not find " + name + " in the database, are you sure it is there?"

    def get_genes(self):
        """ Return a generator of a list of genes """
        return self.db.features_of_type('gene')

    def get_transcripts(self, name):
        """ Method that returns a list of transcripts for a given gene. """
        return self._get_annotation(name, self._featureKey['mRNA'])

    def get_exons(self, name):
        """ Method that returns a list of exons for a given transcript ID """
        return self._get_annotation(name, self._featureKey['exon'])

    def get_introns(self, name):
        """ Method that returns a list of introns for a given transcript ID """
        return self._get_annotation(name, self._featureKey['intron'])

    def get_cds(self, name):
        """ Method that returns a list of CDS for a given transcript ID """
        return self._get_annotation(name, self._featureKey['CDS'])

    def get_utrs(self, name):
        """ Method that returns a list of three prime UTRs and a list of five
        prime UTRs for a given transcript ID 
        """
        return (self._get_annotation(name, self._featureKey['3UTR']), self._get_annotation(name, self._featureKey['5UTR']))

    def get_coords(self, annotations, order_by='start'):
        """ Method that returns a list of coords for a given feature """
        coords = []
        for row in annotations:
            coords.append((row.start, row.end))

        if order_by == 'start':
            coords.sort(key=lambda tup: tup[0])
        elif order_by == 'end':
            coords.sort(key=lambda tup: tup[1])
        return coords

class FlyGff(_Anno):
    """ Flybase specific class for handling an annotation database created from
    GFFUtils using a FlyBase GFF file. 

    Arguments:
    filename (str) = filename of your gff/gtf file
    force (bool) = if True force the recreation of the database
    merge {'merge', 'create_unique', 'error', 'warning'} = passed to the
    merge_strategy option in gffutils.create_db()

    Attributes:
    symbol2fbgn = a dictionary to convert from gene symbol to fbgn
    """

    def __init__(self, filename, force=False, merge='merge'):
        # Import initilize _Anno class
        super(FlyGff, self).__init__(filename, 'gff', force, merge)

        # Update Feature key
        self._featureKey['3UTR'] = 'three_prime_UTR'
        self._featureKey['5UTR'] = 'five_prime_UTR'

        # Create key to go from gene symbol to fbgn
        self.symbol2fbgn = {}
        for gene in self.db.features_of_type('gene'):
            self.symbol2fbgn[gene['Name'][0]] = gene.id

    def get_transcripts(self, name):
        """ Override the get_transcripts method to allow translation of gene
        symbol to FBgn. 
        """
        # If a gene symbol was given look up the FBgn
        if name in self.symbol2fbgn:
            name = self.symbol2fbgn[name]
        return self._get_annotation(name, self._featureKey['mRNA'])

class _Gene(object):
    """ A general class to build a complete gene object with transcript and
    exon information.

    Arguments:
    geneName (str) = a gene name
    gffObject = a mclib gene object

    Attributes:
    gene = geneName
    db = gffObject
    transID = List of transcript ids
    transCnt = Number of transcripts
    transcript = Dictionary of transcript IDs that contain the following:
        exons = List of all exon coords (start, end).
        cds = List of all coding sequence coords (start, end).
        utr = List of UTR region coordinates [3'UTR coords, 5'UTR coords].
    """
    def __init__(self, geneName='', gffObject=None):
        self.gene = geneName

        # Get gene information 
        if self.gene in gffObject.symbol2fbgn:
            self.id = gffObject.symbol2fbgn[self.gene]
        else:
            self.id = self.gene

        try:
            self.chrom = gffObject.db[self.id].chrom
            self.start = gffObject.db[self.id].start
            self.end = gffObject.db[self.id].end
            self.strand = gffObject.db[self.id].strand
        except:
            print """There is something wrong creating your gene. Please check
            and make sure the geneName you used is in your gffDB. If this does
            not fix the problem, check that your gffDB is formated correctly.
            """
            raise ValueError

        # Pull list of transcripts for that gene
        self._transOb = gffObject.get_transcripts(self.gene)
        self.transCnt = len(self._transOb)

        # Create dictionary where key is transcript id and values are
        # annotations for that transcript
        self.transcript = defaultdict(dict)
        for ts in self._transOb:
            self.transcript[ts.id]['exons'] = [(exon.start, exon.end) for exon in gffObject.get_exons(ts)]
            self.transcript[ts.id]['introns'] = [(intron.start, intron.end) for intron in gffObject.get_introns(ts)]
            self.transcript[ts.id]['cds'] = [(cds.start, cds.end) for cds in gffObject.get_cds(ts)]
            self.transcript[ts.id]['utr'] = [gffObject.get_coords(utr) for utr in gffObject.get_utrs(ts)]
            self.transcript[ts.id]['tsStart'] = ts.start
            self.transcript[ts.id]['tsEnd'] = ts.end

class FlyGene(_Gene):
    """ A class to build a complete gene object with transcript and exon
    information from the FlyGFF object 

    Arguments:
    geneName (str) = a gene name
    gffObject = a mclib gene object

    Attributes:
    gene = geneName
    db = gffObject
    transID = List of transcript ids
    transCnt = Number of transcripts
    transcript = Dictionary of transcript IDs that contain the following:
        exons = List of all exon coords (start, end).
        cds = List of all coding sequence coords (start, end).
        utr = List of UTR region coordinates [3'UTR coords, 5'UTR coords].
    """
    def __init__(self, geneName='', gffObject=None):
        # If a gene symbol was given look up the FBgn
        if geneName in gffObject.symbol2fbgn:
            name = gffObject.symbol2fbgn[geneName]
        else:
            name = geneName

        super(FlyGene, self).__init__(name, gffObject)
