"""
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""
from __future__ import division
import re

from itertools import islice
from subprocess import Popen, PIPE
from io import TextIOWrapper
from collections import defaultdict
try:
    from collections import Counter
except ImportError:
    from ._backport import Counter

from six.moves import range, zip


class Gzip(object):
    """ 
    Call system gzip and maintain interface compatibility with python
    gzip module 
    """
    def __init__(self, filename, mode='r'):
        self.stream, self.p = self.open(filename, mode)
        self.mode = mode
        self.filename = filename

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return next(self.stream)

    def open(self, filename, mode):
        if 'r' in mode:
            self.fh = open(filename, 'rb', 0)
            p = Popen(['gzip', '-dc', filename], stdout=PIPE, stderr=PIPE)
            if 'b' in mode:
                fh = p.stdout
            else:
                try:
                    fh = TextIOWrapper(p.stdout)
                except AttributeError:  # python2.7?
                    fh = p.stdout
        elif 'w' in mode:
            self.fh = open(filename, 'wb', 0)
            p = Popen(['gzip', '-c'], stdin=PIPE, stdout=self.fh)
            fh = p.stdin
        return (fh, p)

    def write(self, s):
        self.stream.write(s.encode('utf-8'))

    def read(self, s):
        self.stream.read(s)

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.p.communicate()
        if self.fh:
            self.fh.close()


class Stats:
    """ Counter for characterization of NGS reads """
    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
        self.conv = defaultdict(lambda: defaultdict(int))

    def evaluate(self, seq, qual, conv=None):
        """ 
        Evaluate read object at each position, and fill in nuc and 
        qual dictionaries 
        """
        self.gc[gc(seq)] += 1
        if conv:
            cpgs = cpg_map(seq)
        for i in range(1, len(seq) + 1):
            self.depth[i] += 1
            self.nuc[i][seq[i-1]] += 1
            self.qual[i][qual[i-1]] += 1
            if conv:
                if cpgs[i-1] != 'N':
                    self.conv[i][conv[i-1]] += 1

    def kmercount(self, seq, k=5):
        """ Count all kmers of k length in seq and update kmer counter.
        """
        for kmer in window(seq, n=k):
            self.kmers[kmer] += 1

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def gc(seq):
    """ Return the GC content of as an int
    >>> x = tuple('TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT')
    >>> gc(x)
    30
    """
    g = seq.count('G')
    c = seq.count('C')
    return int((g + c) / len(seq) * 100)


def padbases(bases):
    """
    For each base call in dictionary D, add an entry base:0 if key base 
    does not exist.
    """
    def inner(D):
        for key in bases:
            if key not in D:
                D[key] = 0
    return inner


def percentile(D, percent):
    """
    modified from: http://stackoverflow.com/a/2753343/717419
    Find the percentile of a list of values.

    N - is a dictionary with key=numeric value and value=frequency.
    percent - a float value from 0.0 to 1.0.

    outlier removal: http://www.purplemath.com/modules/boxwhisk3.htm

    return the percentile of the values
    """
    N = sorted(D.keys())  # dict keys
    P = [D[n] for n in N]  # dict values
    if not N:
        return None
    k = (sum(P)) * percent
    l = (sum(P)) * 0.25  # lower quartile
    u = (sum(P)) * 0.75  # upper quartile
    e = int()
    for n,p in zip(N, P):  # find percentile
        e += p
        if e >= k:
            z = n  # value at percentile
            break
    e = int()
    for n,p in zip(N, P):  # find upper quartile
        e += p
        if e >= u:
            uz = n  # value at quartile
            break
    e = int()
    for n,p in zip(N, P):  # find lower quartile
        e += p
        if e >= l:
            lz = n  # value at quartile
            break
    iqd = 1.5 * (uz - lz)  # 1.5 times the inter-quartile distance
    if (z) & (z < lz - iqd):
        return int(lz - iqd)
    elif (z) & (z > uz + iqd):
        return int(uz + iqd)
    elif z:
        return int(z)
    else:
        return N[-1]


def window(seq, n=2):
    """ Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)


def cpg_map(seq):
    """ Return tuple of C/G/N.

    >>> cpg_map('CGCGTAGCCG')
    'CGCGNNNNCG'
    """
    starts = (x.start() for x in re.finditer('CG', ''.join(['N', seq, 'N'])))
    cpgs = ['N'] * (len(seq) + 2)
    for start in starts:
        cpgs[start] = 'C'
        cpgs[start+1] = 'G'
    return ''.join(cpgs[1:-1])

