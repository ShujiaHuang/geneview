# pylint disable: W0622,C0103,R0913,R0902
"""
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""
from __future__ import division
from six import string_types

####
from ._io_util import Gzip

class Fastq(object):
    """
    A class to hold features from fastq reads.
    """
    def __init__(self, name='', seq='', strand='+', qual='', conv=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.conv = conv
        self.i = int()
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)
        assert isinstance(qual, string_types)

    def __iter__(self):
        return self

    def next(self):
        if self.i < len(self):
            value, self.i = self[self.i], self.i + 1
            return value
        else:
            raise StopIteration()

    def __getitem__(self, key):
        if self.conv:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key], self.conv[key])
        else:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key])

    def __next__(self):
        return self.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name[0] != '@':
            self.name = ''.join(['@', self.name])
        if self.conv:
            return '\n'.join(['{0}:YM:Z:{1}'.format(self.name, self.conv),
                              self.seq, self.strand, self.qual])
        else:
            return '\n'.join([self.name, self.seq, self.strand, self.qual])

    def __len__(self):
        return len(self.seq)

    def gc(self):
        """ Return the GC content of self as an int: int((g+c)/len(self)*100)
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT', 
        ...           qual=' # # ##EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.gc()
        30
        """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return int((g + c) / len(self) * 100)


class FastqReader(object):
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """
    def __init__(self, f):
        # 
        if f[-3:] == '.gz':
            self.file = Gzip(f, 'r')
        else:
            self.file = open(f)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        # Return Fastq
        try:
            name = next(self.file).strip().split()[0]  # remove whitespace
            seq = next(self.file).strip()
            strand = next(self.file).strip()
            qual = next(self.file).strip()
            if name.count(':YM:Z:') > 0:
                tag, dtype, data = name.split(':')[-3:]
                name = ':'.join(name.split(':')[:-3])
                return Fastq(name=name, seq=seq, strand=strand, qual=qual, 
                             conv=data)
            else:
                return Fastq(name=name, seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n):
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                if name.count(':YM:Z:') > 0:
                    tag, dtype, data = name.split(':')[-3:]
                    name = ':'.join(name.split(':')[:-3])
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual, conv=data)
                else:
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual)

    def fileno(self):
        return self.file.fileno()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()
