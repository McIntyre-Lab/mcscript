import gffutils
import os.path

class FlyGff(object):
    def __init__(self, filename, force=False, merge='merge'):
        """ A class for dealing with Flybase GFF files.

        filename = filename of your gff file
        force = if True force the recreation of the database
        merge = passed to the merge_strategey option in gffutils.create_db()
        """
        # Database name
        if os.path.splitext(filename)[1] == '.gff':
            dbname = filename + '.db'
        else:
            print "Please provide a gff file."
            raise ValueError

        # Create or connect the GFF database
        if not force:
            try:
                self.db = gffutils.FeatureDB(dbname)
            except:
                self.db = gffutils.create_db(filename, dbname, force=force, merge_strategy=merge)
        else:
                self.db = gffutils.create_db(filename, dbname, force=force, merge_strategy=merge)

        # Create key to go from gene symbol to fbgn
        self.symbol2fbgn = {}
        for gene in self.db.features_of_type('gene'):
            self.symbol2fbgn[gene['Name'][0]] = gene.id

    def _get_annotation(self, name, childtype=''):
        """ Method that returns a list of children given some parent """
        try:
            return [i for i in self.db.children(name, featuretype=childtype, order_by='start')]
        except:
            print "Sorry I could not find " + name + " in the database, are you sure it is there?"

    def get_coords(self, annotations, order_by='start'):
        """ Method that returns a list of coords for a given feature """
        """ Method to return a list of tuples with exon start and stops for a given transcript """
        coords = []
        for row in annotations:
            coords.append((row.start, row.end))

        if order_by == 'start':
            coords.sort(key=lambda tup: tup[0])
        elif order_by == 'end':
            coords.sort(key=lambda tup: tup[1])
        return coords

    def get_transcripts(self, name):
        """ Method that returns a list of transcripts for a given gene. You can use gene symbol or FBgn """
        # If a gene symbol was given look up the FBgn
        if name in self.symbol2fbgn:
            name = self.symbol2fbgn[name]
        return self._get_annotation(name, 'mRNA')

    def get_exons(self, name):
        """ Method that returns a list of exons for a given transcript ID """
        return self._get_annotation(name, 'exon')

    def get_introns(self, name):
        """ Method that returns a list of introns for a given transcript ID """
        return self._get_annotation(name, 'intron')

    def get_cds(self, name):
        """ Method that returns a list of CDS for a given transcript ID """
        return self._get_annotation(name, 'CDS')

    def get_utrs(self, name):
        """ Method that returns a list of three prime UTRs and a list of five
        prime UTRs for a given transcript ID 
        """
        return (self._get_annotation(name, 'three_prime_UTR'), self._get_annotation(name, 'five_prime_UTR'))

    def get_mirs(self):
        """ Return a list of all miRNAs """
        return [mir for mir in self.db.features_of_type('miRNA')]

