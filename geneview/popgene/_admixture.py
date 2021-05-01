"""Plot admixture

Author: Shujia Huang
Date: 2021-05-01
"""
from itertools import cycle

from ..algorithm import hierarchical_cluster


def admixtureplot(
        data=None,
        group_order=None,
        xticklabels=None,
        xticklabel_kws=None,
        palette=None,
        ax=None
):
    """Plot admixture.

    Parameters
    ----------
        data : dict

    """
    TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
             "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"] 
    palette = cycle(TAB10 if palette is None else palette)
        
    if group_order is None:
        group_order = set(data.keys())
    
    n = 0
    for g in group_order:
        n += len(data[g])
        hc = hierarchical_cluster(data=data[g], axis=0)
        data[g] = hc.data.iloc[hc.reordered_index]
        
    x = np.arange(n)
    K = data[group_order[0]].columns
    
    add = np.zeros(len(x))
    for k in K:
        c = next(palette)  # one color for one 'k'

        start_g_pos = 0
        tmp = []
        for g in group_order:  # make group order
            y = data[g][k]
            end_g_pos = start_g_pos + len(y)
            ax.bar(x[start_g_pos:end_g_pos],
                   y,
                   bottom=add[start_g_pos:end_g_pos],
                   color=c, 
                   width=1.0, 
                   linewidth=0, 
                   edgecolor=None)
            start_g_pos = end_g_pos
            tmp.append(y)
        add += np.array(pd.concat(tmp, axis=0))
    
    g_pos = 0
    xticks_pos = []
    for g in group_order:  # make group order
        g_size = len(data[g][k])
        xticks_pos.append(g_pos + 0.5 * g_size)
        
        g_pos += g_size
        ax.axvline(x=g_pos, color="k", lw=1.5)

    ax.spines["left"].set_linewidth(2)
    ax.spines["top"].set_linewidth(2)
    ax.spines["right"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    
    ax.set_xlim([x[0], x[-1]])
    ax.set_ylim([0, 1.0])
    
    if xticklabel_kws is None:
        xticklabel_kws = {}
        
    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(group_order if xticklabels is None else xticklabels, 
                       fontsize=params["axes.labelsize"], 
                       # params["xtick.labelsize"],  # params["axes.labelsize"], 
                       **xticklabel_kws)
    
    ax.tick_params(bottom=False,top=False,left=False,right=False)
    ax.set_yticks([])
    ax.set_xlabel(None)
    ax.set_ylabel("K=%d" % len(K), rotation=0, ha="right")
    
    return ax
