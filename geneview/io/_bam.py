# pylint disable: W0622,C0103,R0913,R0902
"""
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""

from __future__ import division
import sys
import os
from six.moves import range
from itertools import groupby
from subprocess import Popen, PIPE
from io import TextIOWrapper


class Sam(object):
    """ Store fields in each line of a SAM file, provided as a tuple. """
    __slots__ = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 
                 'pnext', 'tlen', 'seq', 'qual', 'tags', '_tags', '_cigars']

    def __init__(self, fields):
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]
        self.tags = None
        self._tags = fields[11:]
        self._cigars = None

    def __gt__(self, other):
        if self.rname != other.rname:
            return self.rname > other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos > other.pos
        else:
            return str(self) > str(other)

    def __lt__(self, other):
        if self.rname != other.rname:
            return self.rname < other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos < other.pos
        else:
            return str(self) < str(other)

    def __eq__(self, other):
        if ((self.rname == other.rname) and (self.pos == other.pos) 
                and (str(self) != str(other))):
            return str(self) == str(other)
        else:
            return self.pos == other.pos

    def __str__(self):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return '\t'.join((self.qname, str(self.flag), self.rname, str(self.pos),
                          str(self.mapq), str(self.cigar), self.rnext, str(self.pnext),
                          str(self.tlen), ''.join(self.seq), ''.join(self.qual)) + \
                          tuple(':'.join((tag, self.tags[tag][0], str(self.tags[tag][1]))) for tag in sorted(self.tags.keys()))) + '\n'
    def __repr__(self):
        return "Sam({0}:{1}:{2})".format(self.rname, self.pos, self.qname)

    def __len__(self):
        return sum(c[0] for c in self.cigars if c[1] in
                   ("M", "D", "N", "EQ", "X", "P"))

    def __getitem__(self, key):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return self.tags[key][1]

    def __setitem__(self, key, value):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        self.tags[key] = value

    def cigar_split(self):
        """ CIGAR grouping function modified from:
        https://github.com/brentp/bwa-meth
        """
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    @property
    def conv(self):
        return self['YM']

    @property
    def cigars(self):
        if not self._cigars:
            self._cigars = tuple(self.cigar_split())
        return self._cigars

    @property
    def mapped(self):
        return not (self.flag & 0x4)

    @property
    def secondary(self):
        return bool(self.flag & 0x100)

    @property
    def reverse(self):
        return bool(self.flag & 0x10)

    @property
    def duplicate(self):
        return bool(self.flag & 0x400)

    def gapped(self, string, gap_char='-'):
        """ Return string with all deletions wrt reference
         represented as gaps '-' and all insertions wrt reference
         removed.
        i: sequence index
        """
        seq = []
        i = 0
        for n, t in self.cigars:
            if t in ("M", "N", "EQ", "X", "P"):
                seq.extend(string[i:i+n])
                i += n
            elif t in ("D",):
                seq.extend(('-',) * n)
            elif t in ("I",):
                i += n
        return ''.join(seq)

    @property
    def coords(self):
        return range(self.pos, self.pos + len(self))


class Reader(object):
    """ Read SAM/BAM format file using iterator. """
    def __init__(self, f):
        name, ext = os.path.splitext(f.name)
        if ext == '.bam':
            BamReaderSamtools.__init__(self, f)
        else:
            SamReader.__init__(self, f)

    def next(self):
        try:
            line = next(self.file).rstrip('\n\r')
            return Sam(tuple(line.split('\t')))
        except StopIteration:
            raise StopIteration

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def subsample(self, n):
        """ Draws every nth read from self. Returns Sam. """
        for i, line in enumerate(self.file):
            if i % n == 0:
                yield Sam(tuple(line.rstrip().split('\t')))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class SamReader(Reader):
    """ Read SAM format file using iterator. """
    def __init__(self, f):
        self.header = []
        self.file = f
        for line in self.file:
            if line[0] == '@':
                self.header.append(line.rstrip('\n\r'))
            else:
                break


class BamReaderSamtools(Reader):
    """ Read BAM format file using iterator. """
    def __init__(self, f):
        pline = ['samtools', 'view', '-H', f.name]
        try:
            p = Popen(pline, bufsize=-1, stdout=PIPE,
                      stderr=PIPE)
        except OSError:
            sys.stderr.write('Samtools must be installed for BAM file support!\n')
            sys.exit(1)
        self.header = [line.decode('utf-8').rstrip('\n\r') for line in p.stdout]
        p.wait()
        pline = ['samtools', 'view', f.name]
        self.p = Popen(pline, bufsize=-1, stdout=PIPE,
                       stderr=PIPE)
        self.file = TextIOWrapper(self.p.stdout)

    def __exit__(self, *args):
        self.p.wait()


def bam_read_count(bamfile):
    """ Return a tuple of the number of mapped and unmapped reads in a bam file """
    p = Popen(['samtools', 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)


def parse_sam_tags(tagfields):
    """ Return a dictionary containing the tags """
    return dict((tag, (dtype, data)) for tag, dtype, data in (decode_tag(x) for x in tagfields))


def encode_tag(tag, data_type, data):
    """ Write a SAM tag in the format ``TAG:TYPE:data``
    >>> encode_tag('YM', 'Z', '#""9O"1@!J')
    'YM:Z:#""9O"1@!J'
    """
    value = ':'.join(list((tag.upper(), data_type.upper(), data)))
    return value


def decode_tag(tag):
    """ Parse a SAM format tag to a (TAG, TYPE, data) tuple.

    TYPE in A, i, f, Z, H, B

    >>> decode_tag('YM:Z:#""9O"1@!J')
    ('YM', 'Z', '#""9O"1@!J')
    >>> decode_tag('XS:i:5')
    ('XS', 'i', 5)
    """
    values = tag.split(':')
    if len(values) != 3:
        values = (values[0], values[1], ':'.join(values[2:]))
    if values[1] == 'i':
        values[2] = int(values[2])
    elif values[1] == 'f':
        values[2] = float(values[2])
    elif values[1] == 'H':
        raise(NotImplementedError, "Hex array SAM tags are currently not parsed.")
    elif values[1] == 'B':
        raise(NotImplementedError, "Byte array SAM tags are currently not parsed.")
    return tuple(values)
