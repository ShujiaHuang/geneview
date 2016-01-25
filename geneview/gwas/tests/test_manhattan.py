"""
    Copytright (c) Shujia Huang
    Date: 2016-01-23

    Plot a manhattan plot of the input file(s).
    python %prog [options] files
"""

import optparse
import sys
from itertools import groupby, cycle
from operator import itemgetter
from matplotlib import pyplot as plt
import numpy as np

def _gen_data(fhs, columns, sep):
    """
    iterate over the files and yield chr, start, pvalue
    """
    for fh in fhs:
        for line in fh:
            if line[0] == '#': continue
            toks = line.strip('\n').split(sep) if sep else line.strip('\n').split()
            yield toks[columns[0]], int(toks[columns[1]]), float(toks[columns[2]])

def chr_cmp(a, b):
    a = a.lower().replace('_', '') 
    b = b.lower().replace('_', '')
    achr = a[3:] if a.startswith('chr') else a
    bchr = b[3:] if b.startswith('chr') else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        # X Y
        return cmp(achr, bchr)

def chr_loc_cmp(alocs, blocs):
    return chr_cmp(alocs[0], blocs[0]) or cmp(alocs[1], blocs[1])

def manhattan(fhs, columns, image_path, no_log, colors, sep, title, lines, ymax):

    if ',' in colors: colors = colors.split(',')
    colors = cycle(colors)

    last_x = 0
    data = sorted(_gen_data(fhs, columns, sep), cmp=chr_loc_cmp)

    xs, ys, cs = [], [], []
    xs_by_chr = {}
    for seqid, rlist in groupby(data, key=itemgetter(0)):
        color = colors.next()
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        xs.extend(region_xs)
        ys.extend([r[2] for r in rlist])
        cs.extend([color] * len(rlist))

        xs_by_chr[seqid] = (region_xs[0] + region_xs[-1]) / 2

        # keep track so that chrs don't overlap.
        last_x = xs[-1]

    cs = np.array(cs)
    xs = np.array(xs)
    ys = np.array(ys) if no_log else -np.log10(ys)
    #idx = ys>4; xs = xs[idx]; ys = ys[idx]; cs = cs[idx]

    xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), cmp=chr_cmp)]

    # Plotting the manhattan image
    #plt.style.use('ggplot')
    plt.close() # in case plot accident
    f, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 8), tight_layout=True) 

    ax.tick_params(labelsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if title is not None: plt.title(title, fontsize=18)

    if lines:
        ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
    else:
        ax.scatter(xs, ys, s=20, c=cs, alpha=0.8, edgecolors='none')

    # plot 0.05 line after multiple testing.
    #ax.axhline(y=-np.log10(0.05 / len(data)), color='0.5', linewidth=1)
    ax.axhline(y=5, color='b')
    ax.axhline(y=7, color='r')

    ax.set_xlim(0, xs[-1])
    ax.set_ylim(ymin=0)
    if ymax is not None: plt.ylim(ymax=ymax)

    xtick_label_set = set(range(15) + [16,18,20,22])
    plt.xticks([c[1] for c in xs_by_chr if int(c[0]) in xtick_label_set], 
               [c[0] for c in xs_by_chr if int(c[0]) in xtick_label_set])
    #ax.set_xticks([c[1] for c in xs_by_chr if int(c[0]) in xtick_label_set])
    #ax.set_xticklabels([c[0] for c in xs_by_chr if int(c[0]) in xtick_label_set])

    ax.set_xlabel('Chromosome', fontsize=18)
    ax.set_ylabel('-Log10 (P-value)', fontsize=18)
    print >> sys.stderr, 'saving to: %s' % image_path

    plt.show()
    plt.savefig(image_path)

def get_filehandles(args):
    return (open(a) if a != "-" else sys.stdin for a in args)

def main():
    COLORFUL = '#6DC066,#FD482F,#8A2BE2,#3399FF'
    p = optparse.OptionParser(__doc__)
    p.add_option("--no-log", dest="no_log", help="don't do -log10(p) on the value",
                 action='store_true', default=False)
    p.add_option("--cols", dest="cols", help="zero-based column indexes to get "
                 "chr, position, p-value respectively e.g. %default", 
                 default="0,1,2")
    p.add_option("--colors", dest="colors", help="cycle through these colors",
                 default="#000000,#969696")
    p.add_option("--colorful", default=False, dest="colorful", action="store_true",
                 help="if set than enhance to cycle through these colors: "
                 "%s" % COLORFUL)
    p.add_option("--image", dest="image", 
                 help="save the image_path to this file. e.g. %default",
                 default="manhattan.png")
    p.add_option("--title", help="title for the image.", default=None, 
                 dest="title")
    p.add_option("--ymax", help="max (logged) y-value for plot", dest="ymax", 
                 type='float')
    p.add_option("--sep", help="data separator, default is any space",
                 default=None, dest="sep")
    p.add_option("--lines", default=False, dest="lines", action="store_true",
                 help="plot the p-values as lines extending from the x-axis "
                 "rather than points in space. plotting will take longer "
                 "with this option.")

    opts, args = p.parse_args()
    if opts.colorful: opts.colors = COLORFUL
    if (len(args) == 0):
        sys.exit(not p.print_help())

    fhs = get_filehandles(args)
    columns = map(int, opts.cols.split(","))
    manhattan(fhs, columns, opts.image, opts.no_log, opts.colors, opts.sep,
              opts.title, opts.lines, opts.ymax)

if __name__ == "__main__":
    main()
