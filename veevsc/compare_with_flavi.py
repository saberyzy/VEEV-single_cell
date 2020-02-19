# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/11/19
content:    Compare VEEV with flavis
'''
import os
import sys
import numpy as np
import pandas as pd
import loompy

import ternary
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


def get_trajectories(ds, gene, normalize=True):
    x = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])
    y = np.log10(0.1 + ds.counts.loc[gene])
    v = ds.samplesheet['virus']

    xbins = np.linspace(np.log10(2), 5, 30)
    xbinc = 0.5 * (xbins[:-1] + xbins[1:])
    xker = np.ones(5)
    xint = np.convolve(xbinc, xker, mode='valid') / xker.sum()
    yav = {}
    for vir in ['DENV', 'ZIKV', 'VEEV']:
        tmp = []
        for ib in range(len(xbins) - 1):
            yi = y[(v == vir) & (x >= xbins[ib]) & (x < xbins[ib+1])]
            tmp.append(yi.mean())
            yav[vir] = tmp
        yint = np.convolve(tmp, xker, mode='valid') / xker.sum()
        yav[vir] = yint

    if normalize:
        ysum = sum([(10**v) - 0.1 for v in yav.values()])
        for vir in yav:
            yav[vir] = (10**yav[vir] - 0.1) / ysum

    return {'x': xint, 'ys': yav}


def get_points_from_traj(data, order):
    points = []
    for i in range(len(data['x'])):
        points.append([data['ys'][key][i] for key in order])
    points = np.array(points)
    return points


if __name__ == '__main__':

    print('Load VEEV + DENV + ZIKA')
    ds0 = singlet.Dataset(dataset={
        'path': '../data/dengue_Zika/all3.loom',
        'index_samples': 'CellID',
        'index_features': 'GeneName',
        'bit_precision': 32,
        'format': 'loom',
        })
    ds = ds0.query_samples_by_metadata('virus_reads_per_million >= 2')

    print('Find genes that start at similar levels and diverge')
    dsp = ds.split('virus')
    genes = pd.Series(np.ones(ds.n_features, bool), index=ds.featurenames)
    tmps = {}
    for key1 in dsp:
        ind1 = dsp[key1].samplesheet['virus_reads_per_million'] <= 10
        cou1 = dsp[key1].counts.loc[:, ind1].mean(axis=1)
        tmps[key1] = cou1
        for key2 in dsp:
            if key1 == key2:
                continue
            if key2 in tmps:
                cou2 = tmps[key2]
            else:
                ind2 = dsp[key2].samplesheet['virus_reads_per_million'] <= 10
                cou2 = dsp[key2].counts.loc[:, ind2].mean(axis=1)
                tmps[key2] = cou2
            genes &= np.abs(np.log10(cou1 + 0.1) - np.log10(cou2 + 0.1)) < 1

    print('Find genes with the most different correlation')
    corrp = {
        vir: dsp[vir].correlation.correlate_features_phenotypes('virus_reads_per_million').fillna(0)
        for vir in dsp
        }
    genes_top = set()
    for key1 in dsp:
        for key2 in dsp:
            if key1 == key2:
                continue
            cdiff = np.abs(corrp[key1] - corrp[key2])
            gcand = cdiff.nlargest(5).index
            gcand = gcand[genes[gcand]].tolist()

            order = ['DENV', 'ZIKV', 'VEEV']
            gcand2 = []
            for gene in gcand:
                data = get_trajectories(ds, gene)
                stddevs = get_points_from_traj(data, order).std(axis=1)
                if stddevs[-1] > stddevs[0]:
                    gcand2.append(gene)
            gcand = gcand2

            genes_top |= set(gcand)
    genes_top = sorted(genes_top)
    genes_plotall = genes_top

    print('Get trajectories')
    points_dict_nn = {
        gene: get_points_from_traj(get_trajectories(ds, gene, normalize=False), order)
        for gene in genes_plotall}
    points_dict = {gene: get_points_from_traj(get_trajectories(ds, gene), order) for gene in genes_plotall}
    # Group genes for clarity
    groups = [[], [], []]
    for gene, points in points_dict.items():
        if points[-1][0] > 0.6:
            groups[0].append(gene)
        elif points[-1][1] > 0.6:
            groups[1].append(gene)
        else:
            groups[2].append(gene)

    if False:
        print('Plot 3D')
        order = ['DENV', 'ZIKV', 'VEEV']
        fig = plt.figure(figsize=(12, 6))
        axs = []
        axs.append(fig.add_subplot(131, projection='3d'))
        axs.append(fig.add_subplot(132, projection='3d'))
        axs.append(fig.add_subplot(133, projection='3d'))
        for igp, genes_plot in enumerate(groups):
            ax = axs[igp]
            colors = sns.color_palette('husl', n_colors=len(genes_plot))
            for ig, (gene, color) in enumerate(zip(genes_plot, colors)):
                ls = ['-', '--'][ig % 2]
                marker = ['*', 's'][ig % 2]
                markersize = [12, 7][ig % 2]
                points = points_dict_nn[gene]
                xs, ys, zs = points.T
                ax.plot(xs, ys, zs, lw=2, color=color, alpha=0.5)
                ax.scatter(
                    [xs[-1]], [ys[-1]], [zs[-1]],
                    color=color, s=[markersize],
                    ls=ls, zorder=10, alpha=0.75, label=gene)
            ax.set_xlabel(order[0])
            ax.set_ylabel(order[1])
            ax.set_zlabel(order[2])
            ax.set_xlim(-1, 5)
            ax.set_ylim(-1, 5)
            ax.set_zlim(-1, 5)
            ax.legend(
                    ncol=1, title='Gene:',
                    loc='upper left',
                    bbox_to_anchor=[0.01, -0.1],
                    bbox_transform=ax.transAxes,
                    )
            ax.set_xticks([-1, 0, 1, 2, 3, 4, 5])
            ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
            ax.set_zticks([-1, 0, 1, 2, 3, 4, 5])
            ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_zticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        fig.tight_layout(w_pad=0.3)

    if False:
        print('Plot ternary')
        order = ['DENV', 'ZIKV', 'VEEV']
        for igp, genes_plot in enumerate(groups):
            fig, tax = ternary.figure(scale=1.0)
            fig.set_size_inches(5, 3)
            tax.boundary()
            tax.gridlines(multiple=0.2, color="black")
            colors = sns.color_palette('husl', n_colors=len(genes_plot))
            for ig, (gene, color) in enumerate(zip(genes_plot, colors)):
                ls = ['-', '--'][ig % 2]
                marker = ['*', 's'][ig % 2]
                markersize = [12, 7][ig % 2]
                points = points_dict[gene]

                # Plot with increasing alpha
                ngroups = 5
                alphas = np.linspace(0.3, 0.6, ngroups)
                gs = len(points) // ngroups
                for igr in range(ngroups):
                    if igr == 0:
                        pointsi = points[:gs]
                    elif igr == ngroups - 1:
                        pointsi = points[(igr * gs) - 1:]
                    else:
                        pointsi = points[(igr * gs) - 1: (igr + 1) * gs]
                    tax.plot(pointsi, linewidth=2.0, alpha=alphas[igr], color=color, ls=ls)
                tax.plot([points[-1]], marker=marker, ls=ls, label=gene, color=color, markersize=markersize, zorder=10, alpha=0.75)
            tax.ticks(axis='lbr', multiple=0.2, linewidth=1, tick_formats="%.1f", offset=0.02)
            tax.get_axes().axis('off')
            tax.clear_matplotlib_ticks()

            offset = 0.14
            tax.left_axis_label("VEEV", offset=offset)
            tax.right_axis_label("ZIKV", offset=offset)
            tax.bottom_axis_label("DENV", offset=offset)
            tax.right_corner_label("DENV\nonly", offset=2.6 * offset, va='top')
            tax.top_corner_label("ZIKV\nonly", offset=2 * offset)
            tax.left_corner_label("VEEV\nonly", offset=2.2 * offset, va='top')

            tax.legend(
                    ncol=1, title='Gene:',
                    loc='upper left',
                    bbox_to_anchor=[1.17, 1.01],
                    bbox_transform=fig.get_axes()[0].transAxes,
                    )
            tax.show()
            fig.tight_layout()
            fig.savefig('../figures/ternary_plot_group_{:}.png'.format(igp+1), dpi=600)
            fig.savefig('../figures/ternary_plot_group_{:}.pdf'.format(igp+1))
            fig.savefig('../figures/ternary_plot_group_{:}.svg'.format(igp+1))

    if False:
        print('Scatter a few genes')
        #genes = ['DDIT3', 'ACTG1', 'XBP1', 'HSPA5', 'TAF7', 'DUSP14', 'DDIT4', 'DNAJB9',
        #         'EIF1', 'EIF5', 'TMED2', 'SKIL', 'SSR2', 'SEC61B', 'SOX4']
        genes = genes_top
        fig, axs = plt.subplots(3, 5, figsize=(15, 8), sharex=True)
        axs = axs.ravel()
        for i, (ax, gene) in enumerate(zip(axs, genes)):
            x = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])
            y = np.log10(0.1 + ds.counts.loc[gene])
            cmap = {'DENV': 'tomato', 'ZIKV': 'steelblue', 'VEEV': 'green'}
            v = ds.samplesheet['virus']
            c = v.map(cmap).values

            xbins = np.linspace(np.log10(2), 5, 30)
            xbinc = 0.5 * (xbins[:-1] + xbins[1:])
            yav = {}
            for vir in ['DENV', 'ZIKV', 'VEEV']:
                tmp = []
                for ib in range(len(xbins) - 1):
                    yi = y[(v == vir) & (x >= xbins[ib]) & (x < xbins[ib+1])]
                    tmp.append(yi.mean())
                    yav[vir] = tmp

                xker = np.ones(5)
                xint = np.convolve(xbinc, xker, mode='valid') / xker.sum()
                yint = np.convolve(tmp, xker, mode='valid') / xker.sum()
                ax.plot(xint, yint, lw=2, color=cmap[vir])

            ax.scatter(x, y, s=12, alpha=0.08, c=c)
            ax.grid(True)
            ax.set_xlim(np.log10(2), 5)
            ax.set_ylim(-1, 5)
            ax.set_xlabel('log10[vRNA cpm]')
            if (i % 2) == 0:
                ax.set_ylabel('log10[mRNA cpm]')
            ax.set_title(gene)
        fig.tight_layout()
        fig.savefig('../figures/scatter_plot_for_ternary_allgroups.png', dpi=600)
        fig.savefig('../figures/scatter_plot_for_ternary_allgroups.pdf')
        fig.savefig('../figures/scatter_plot_for_ternary_allgroups.svg')

    if False:
        print('Scatter fewer genes')
        genes = ['HSPA5', 'NRBF2', 'SERP1']
        fig, axs = plt.subplots(1, 3, figsize=(9, 3), sharex=True)
        axs = axs.ravel()
        for i, (ax, gene) in enumerate(zip(axs, genes)):
            x = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])
            y = np.log10(0.1 + ds.counts.loc[gene])
            cmap = {'DENV': 'tomato', 'ZIKV': 'steelblue', 'VEEV': 'green'}
            v = ds.samplesheet['virus']
            c = v.map(cmap).values

            xbins = np.linspace(np.log10(2), 5, 30)
            xbinc = 0.5 * (xbins[:-1] + xbins[1:])
            yav = {}
            for vir in ['DENV', 'ZIKV', 'VEEV']:
                tmp = []
                for ib in range(len(xbins) - 1):
                    yi = y[(v == vir) & (x >= xbins[ib]) & (x < xbins[ib+1])]
                    tmp.append(yi.mean())
                    yav[vir] = tmp

                xker = np.ones(5)
                xint = np.convolve(xbinc, xker, mode='valid') / xker.sum()
                yint = np.convolve(tmp, xker, mode='valid') / xker.sum()
                ax.plot(xint, yint, lw=2, color=cmap[vir])

            ax.scatter(x, y, s=12, alpha=0.08, c=c)
            ax.grid(True)
            ax.set_xlim(np.log10(2), 5)
            ax.set_ylim(-1, 5)
            ax.set_xticks([0, 1, 2, 3, 4, 5])
            ax.set_xticklabels(['$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
            ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            ax.set_xlabel('vRNA [cpm]')
            if i == 0:
                ax.set_ylabel('Gene expression [cpm]')
            ax.set_title(gene)
        fig.tight_layout()
        fig.savefig('../figures/scatter_plot_for_ternary_3genes.png', dpi=600)
        fig.savefig('../figures/scatter_plot_for_ternary_3genes.pdf')
        fig.savefig('../figures/scatter_plot_for_ternary_3genes.svg')

    if True:
        print('Plot ternary of 3 genes')
        order = ['DENV', 'ZIKV', 'VEEV']
        genes_plot = ['HSPA5', 'NRBF2', 'SERP1']
        colors = ['tomato', 'steelblue', 'green']
        scale = 100
        fig, tax = ternary.figure(scale=scale)
        fig.set_size_inches(5, 3)
        tax.boundary()
        for ig, (gene, color) in enumerate(zip(genes_plot, colors)):
            ls = ['-', '-'][ig % 2]
            marker = ['*', '*'][ig % 2]
            markersize = [12, 12][ig % 2]
            points = points_dict[gene]

            # Plot with increasing alpha
            ngroups = 5
            alphas = np.linspace(0.3, 0.6, ngroups)
            gs = len(points) // ngroups
            for igr in range(ngroups):
                if igr == 0:
                    pointsi = points[:gs]
                elif igr == ngroups - 1:
                    pointsi = points[(igr * gs) - 1:]
                else:
                    pointsi = points[(igr * gs) - 1: (igr + 1) * gs]
                tax.plot(scale * pointsi, linewidth=2.0, alpha=alphas[igr], color=color, ls=ls, zorder=9)
            tax.plot([scale * points[-1]], marker=marker, ls=ls, label=gene, color=color, markersize=markersize, zorder=10, alpha=0.75)

        def generate_heatmap_data(scale=100, alpha=0.03):
            from matplotlib.colors import to_rgba
            from ternary.helpers import simplex_iterator
            def color_point(i, j, k, scale):
                if i > 0.5 * scale:
                    col = to_rgba(colors[0])
                elif j > 0.5 * scale:
                    col = to_rgba(colors[1])
                elif k > 0.5 * scale:
                    col = to_rgba(colors[2])
                else:
                    col = to_rgba('grey')
                return (col[0], col[1], col[2], alpha)

            d = dict()
            for (i, j, k) in simplex_iterator(scale=scale):
                d[(i, j, k)] = color_point(i, j, k, scale)
            return d
        hdata = generate_heatmap_data(scale=scale)
        tax.heatmap(hdata, use_rgba=True, scale=scale, colorbar=False)

        tax.ticks(axis='lbr', multiple=0.5 * scale, linewidth=2, tick_formats="%.0f", offset=0.04)
        #tax.gridlines(color="blue", multiple=0.25, linewidth=0.5)
        tax.gridlines(color="black", multiple=0.5 * scale)
        tax.get_axes().axis('off')
        tax.clear_matplotlib_ticks()

        offset = 0.14
        tax.left_axis_label("% VEEV", offset=offset*2)
        tax.right_axis_label("% ZIKV", offset=offset*2)
        tax.bottom_axis_label("% DENV", offset=offset*2)
        tax.right_corner_label("\n\n\nDENV\nonly", offset=4.8 * offset, va='top')
        tax.top_corner_label("ZIKV\nonly", offset=2 * offset)
        tax.left_corner_label("\n\nVEEV\nonly", offset=4.8 * offset, va='top')

        tax.legend(
                ncol=1, title='Gene:',
                loc='upper left',
                bbox_to_anchor=[1.17, 1.01],
                bbox_transform=fig.get_axes()[0].transAxes,
                )
        tax.show()
        fig.tight_layout()
        fig.savefig('../figures/ternary_plot_group_3genes.png', dpi=600)
        fig.savefig('../figures/ternary_plot_group_3genes.pdf')
        fig.savefig('../figures/ternary_plot_group_3genes.svg')


    plt.ion()
    plt.show()
