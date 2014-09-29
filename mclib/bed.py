import re

class Bed:
    """ Basic class to handle the reading and access of BED files. A BED object
    can be genereted by:

    myBedObject = mclib.Bed(filename)
    """
    def __init__(self, filename, keys=None):
        # List of possible bed columns with their index, from UCSC website
        self._keys2index = { 'chrom' : 0, 'chromStart': 1, 'chromEnd': 2, 'name': 3,
                             'score': 4, 'strand': 5, 'thickStart': 6, 'thickEnd': 7,
                             'itemRgb': 8, 'blockCount': 9, 'blockSizes': 10, 'blockStarts': 11,
                             'start': 1, 'end':2 }

        with open(filename, 'r') as FH:
            self._rows = [BedRow(line, parent=self) for line in FH.readlines() ]

    def __getitem__(self, rowIndex, colIndex=None):
        if colIndex == None:
            return self._rows[rowIndex]
        return self._rows[rowIndex][colIndex]

    def get_all_rows(self):
        """ Returns all rows of the BED file """
        return [row for row in self._rows ]

    def get_rows(self, name='', exact=True):
        """ Returns a specific row based on a name given by the user. Looks
        first in the 'name' column 4, then in the 'chrom' column 1. If no name
        is given then all rows are returned.

        myBedObject.get_rows(name="F10001_SI")
        """
        if name == '':
            self.get_all_rows()
        else:
            if exact:
                try:
                    return [row for row in self._rows if (name == row[self._keys2index['name']]) or (name == row[self._keys2index['chrom']])]
                except:
                    return [row for row in self._rows if (name == row[self._keys2index['chrom']])]
            else:
                try:
                    return [row for row in self._rows if (re.search(name, row[self._keys2index['name']])) or (re.search(name, row[self._keys2index['chrom']]))]
                except:
                    return [row for row in self._rows if (re.search(name, row[self._keys2index['chrom']]))]

class BedRow:
    """ Basic class to represnt a row in a BED file. """
    def __init__(self, line, parent=None):
        self._contents = self._split_bedRow(line) 
        self._parent = parent

    def _split_bedRow(self, line):
        # split the row by tabs
        cols = line.strip().split('\t')

        # convert numeric values to intergers
        self._make_int(cols,1)
        self._make_int(cols,2)
        self._make_int(cols,4)
        self._make_int(cols,6)
        self._make_int(cols,7)
        self._make_int(cols,9)

        # Split itemRgb, blockSizes and blockStarts into tuples
        self._split_col(cols,8)
        self._split_col(cols,10)
        self._split_col(cols,11)
        return cols
    
    def _make_int(self, cols, index):
        try:
            cols[index] = int(cols[index])
        except:
            pass
        
    def _split_col(self, cols, index):
        try:
            cols[index] = tuple(int(x) for x in cols[index].split(','))
        except:
            pass

    def __repr__(self):
        return "\t".join(str(x) for x in self._contents)

    def __getitem__(self, item):
        if type(item) in [int, long]:
            return self._contents[item]
        try:
            return self._contents[self._parent._keys2index[item]]
        except:
            print "Your bed file did not have this column."

    def get_coords(self):
        return (self['chromStart'], self['chromEnd'])

    def get_thickRange(self):
        return (self['thickStart'], self['thickEnd'])
