#!/usr/bin/env python

"""
Plotting functions in toytree/toyplot to generate HTML 
and PNG format figures.
"""

import numpy as np
import toytree
import toyplot
import toyplot.html, toyplot.png, toyplot.browser


def dist_RI_scatterplot(ridata, outdir="/tmp", jitter=0.01):
    """
    Scatterplot of RI values x genetic dist with the 
    genetic distance values jittered slightly to show
    more clearly despite overlap.
    """

    # get canvas and axes
    canvas = toyplot.Canvas(width=300, height=250)
    axes = canvas.cartesian(
        label="Observations",
        xlabel="Genetic dist.",
        ylabel="Reproductive isolation",
    )

    # jitter points w/ RI=0 to increase visibiliby (constrain in 0-1)
    xdata = ridata.dist[ridata.RI == 0]
    xdata += np.random.uniform(-jitter, jitter, size=xdata.size)
    xdata[xdata > 1] = 1
    xdata[xdata < 0] = 0    
    mark = axes.scatterplot(
        xdata,
        ridata.RI[ridata.RI == 0],
        size=10,
        opacity=0.1,
        marker="|",
        mstyle={
            "stroke": 'black',
            "stroke-width": 3,
        },
    )

    # jitter points w/ RI=1 to increase visibiliby (constrain in 0-1)    
    xdata = ridata.dist[ridata.RI == 1]
    xdata += np.random.uniform(-jitter, jitter, size=xdata.size)
    xdata[xdata > 1] = 1
    xdata[xdata < 0] = 0    
    mark = axes.scatterplot(
        xdata,
        ridata.RI[ridata.RI == 1],
        size=10,
        opacity=0.1,
        marker="|",
        mstyle={
            "stroke": 'red',
            "stroke-width": 3,
        },
    )
    axes.x.ticks.show = True
    axes.y.ticks.show = True
    toyplot.browser.show(canvas)
    # return canvas, axes, mark



def heatmap_tree_plot(tree, clades, ridata, **kwargs):
    """
    Sanity check: does the data make sense?
    """
    # get canvas size
    canvas = toyplot.Canvas(
        width=kwargs.get('width', 600),
        height=kwargs.get('height', 600),
    )

    # add tree to canvas
    ax0 = canvas.cartesian(
        bounds=("5%", "25%", "5%", "95%"),
        show=False
    )
    tree.draw(
        axes=ax0, 
        layout='r', 
        tip_labels=False,
        edge_colors=tree.get_edge_values_mapped({
            j: toytree.colors[i] for i,j in enumerate(clades)
        })
    )

    # add heatmap to canvas
    ax1 = canvas.table(
        rows=tree.ntips, 
        columns=tree.ntips, 
        bounds=("27%", "95%", "5%", "95%"),
        margin=20
    )
    ax1.cells.cell[:].style = {
        "fill": 'lightgrey',
        "stroke": "none"
    }

    # set values for data
    for idx in ridata.index:

        # get the tip indices
        ridx, cidx = map(int, ridata.loc[idx, ['sidx0', 'sidx1']])

        # row indexing starts from the bottom, not the top, so flip
        rridx = tree.ntips - ridx - 1
        ccidx = tree.ntips - cidx - 1

        # convert all indexers to ints
        ridx, cidx, rridx, ccidx = map(int, [ridx, cidx, rridx, ccidx])

        # get RI for this cell
        value = ridata.at[idx, 'RI']
        if value:
            color = "red"
        if not value:
            color = "black"

        # fill color for observation on both sides of diagonal
        ax1.cells.cell[rridx, cidx].style = {
            "fill": color,
            "stroke": "none"
        }
        # ax1.cells.cell[cidx, ridx].style = {
        #     "fill": color,
        #     "stroke": "none"
        # }

        # # fill hover for observation on both sides of diagonal        
        # hover = (
        #     "{} x {}; logit={:.3f}".format(
        #         tree.idx_dict[ridx].name,
        #         tree.idx_dict[ccidx].name,
        #         logit,
        #     ))
        # ax1.cells.cell[ccidx, ridx].title = hover
        # ax1.cells.cell[ridx, ccidx].title = hover        

    # dividers
    ax1.body.gaps.columns[...] = 0.
    ax1.body.gaps.rows[...] = 0.

    toyplot.browser.show(canvas)
    # return canvas
