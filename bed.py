#!/usr/bin/env python
import re

class Bed:
    def __init__(self, filename, keys=None):
        # List of possible bed columns with their index, from UCSC website
        self._keys = { 'chrom' : 0, 'chromStart': 1, 'chromEnd': 2, 'name': 3,
                'score': 4, 'strand': 5, 'thickStart': 6, 'thickEnd': 7,
                'itemRgb': 8, 'blockCount': 9, 'blockSizes': 10, 'blockStarts': 11 }

        with open(filename, 'r') as FH:
            self._rows = [BedRow(line, parent=self) for line in FH.readlines() ]

    def __getitem__(self, rowIndex, colIndex=None):
        if colIndex == None:
            return self._rows[rowIndex]
        return self._rows[rowIndex][colIndex]

    def get_all_rows(self):
        return [row for row in self._rows ]

    def get_rows(self, name='', exact=True):
        if name == '':
            self.get_all_rows()
        if exact:
            return [row for row in self._rows if (name == row[self._keys['name']]) or (name == row[self._keys['chrom']]) ]
        else:
            return [row for row in self._rows if (re.search(name, row[self._keys['name']])) or (re.search(name, row[self._keys['chrom']]))]

class BedRow:
    def __init__(self, line, parent=None):
        self._contents = line.strip().split('\t') 
        self._parent = parent

    def __repr__(self):
        return "\t".join(self._contents)

    def __getitem__(self, item):
        if type(item) in [int, long]:
            return self._contents[item]
        return self._contents[self._parent._key2index[item]]

    def get_thickRange(self):
        return (self['thickStart'], self['thickEnd'])
