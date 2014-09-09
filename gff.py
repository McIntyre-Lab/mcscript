import gffutils
import os.path

class FlyGff(object):
    def __init__(self, filename, force=False, merge='merge'):
        """ A class for dealing with Flybase GFF files.

        filename = filename of your gff file
        force = if True force the recreation of the database
        merge = passed to the merge_strategey option in gffutils.create_db()
        """

        # Create or connect the GFF database
        dbname = filename + '.db'
        self.db = gffutils.create_db(filename, dbname, force=force, merge_strategy=merge)

        # Create key to go from gene symbol to fbgn
        self.symbol2fbgn = {}
        for gene in self.db.features_of_type('gene'):
            self.symbol2fbgn[gene['Name']] = gene.id

    def _get_annotation(self, name, childtype=''):
        """ Method that returns a list of children given some parent """
        try:
            return [i for i in self.db.children(name, featuretype=childtype, order_by='start')]
        except:
            print "Sorry I could not find " + name + " in the database, are you sure it is there?"

    def _get_coords(self, name, featuretype=''):
        """ Method that returns a list of coords for a given feature """
        """ Method to return a list of tuples with exon start and stops for a given transcript """
        coords = []
        for row in self._get_annotation(name, ):
            coords.append((row.start, row.end))
        return coords

    def get_transcripts(self, name):
        """ Method that returns a list of transcripts for a given gene. You can use gene symbol or FBgn """
        # If a gene symbol was given look up the FBgn
        if name in self.symbol2fbgn:
            name = self.symbol2fbgn[name]
        return self._get_annotation(name, 'mRNA')

    def get_exons(self, name, withCoords=False):
        """ Method that returns a list of exons for a given transcript ID """
        if not withCoords:
            return self._get_annotation(name, 'exon')
        else:
            exons = self._get_annotation(name, 'exon')
            coords = self._get_coords(name, 'exon')
            return exons, coords

    def get_introns(self, name, withCoords=False):
        """ Method that returns a list of introns for a given transcript ID """
        if not withCoords:
            return self._get_annotation(name, 'intron')
        else:
            introns = self._get_annotation(name, 'intron')
            coords = self._get_coords(name, 'intron')
            return introns, coords

    def get_cds(self, name, withCoords=False):
        """ Method that returns a list of CDS for a given transcript ID """
        if not withCoords:
            return self._get_annotation(name, 'CDS')
        else:
            cds = self._get_annotation(name, 'CDS')
            coords = self._get_coords(name, 'CDS')
            return cds, coords

    def get_three_prime_utr(self, name, withCoords=False):
        """ Method that returns a list of three prime UTRs for a given transcript ID """
        if not withCoords:
            return self._get_annotation(name, 'three_prime_UTR')
        else:
            utr = self._get_annotation(name, 'three_prime_UTR')
            coords = self._get_coords(name, 'three_prime_UTR')
            return utr, coords

    def get_five_prime_utr(self, name, withCoords=False):
        """ Method that returns a list of five prime UTRs for a given transcript ID """
        if not withCoords:
            return self._get_annotation(name, 'five_prime_UTR')
        else:
            utr = self._get_annotation(name, 'five_prime_UTR')
            coords = self._get_coords(name, 'five_prime_UTR')
            return utr, coords

    def get_mirs(self):
        """ Return a list of all miRNAs """
        return [mir for mir in self.db.features_of_type('miRNA')]

