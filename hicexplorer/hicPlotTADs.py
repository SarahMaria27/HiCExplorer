
"""

[hic]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/hic_data/s2/hic-norm/RF/merge.npz
title =
colormap = RdYlBu_r
depth = 100000
#min_value =2.8
#max_value = 3.0
transform = log1p
boundaries_file = /home/ramirez/p/HiC-ThomasLing/doc/hic_paper/revision/figure_boundaries_method/conductance_vs_hic/boundaries_all.bed
show_masked_bins = yes
x labels = yes
type = arcplot
type = interaction

[x-axis]

[spacer]

[bigwig]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/external/Graveley_mRNA-seq/GSM390060_Kc167-4_spa.bw
title = Kc RNA-seq (Cherbas et al.)
color = black
min_value = 0
#max_value = auto
width = 1.5
number of bins = 500
nans to zeros = True

[genes]
file = dm3.genes_with_symbol_chrX.bed
title = genes
color = darkblue
width = 5
type = genes

[chrom states]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/ChromatinStates/chromatinStates_kc.bed
title =
color = black
display = collapsed
width = 0.3


[bedgraph matrix]
# a bedgraph matrix file is like a bedgraph, except that per bin there
# are more than one value separated by tab: E.g.
# chrX	18279	40131	0.399113	0.364118	0.320857	0.274307
# chrX	40132	54262	0.479340	0.425471	0.366541	0.324736
file = spectra_conductance.bm ,
title = conductance spectra ,
width = 1.5,
orientation = inverted
min_value = 0.10
max_value = 0.70
"""

from __future__ import division
import sys

import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import enlarge_bins
from hicexplorer._version import __version__

from scipy import sparse
from scipy.sparse import triu, dia_matrix

import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import mpl_toolkits.axisartist as axisartist
from bx.intervals.intersection import IntervalTree, Interval

def parseArguments(args=None):
    parser = argparse.ArgumentParser(
        description='Plots the diagonal,  and some values close to '
        'the ediagonal of a  HiC matrix. The diagonal of the matrix is '
        'plotted horizontally for a region. I will not draw the diagonal '
        'for the whole chromosome')

    # define the arguments
#    parser.add_argument('--tracks',
#                        help='List of files to plot. Type of plot is decided '
#                        'based on the file ending. ',
#                        nargs='+',
#                        )
    parser.add_argument('--tracks',
                        help='File containing the instructions to plot the tracks ',
                        type=argparse.FileType('r'),
                        required=True,
                        )

    parser.add_argument('--width',
                        help='figure width in inches :p',
                        type=float,
                        default=20)

    parser.add_argument('--height',
                        help='figure height in inches :p',
                        type=float)

    parser.add_argument('--title', '-t',
                        help='Plot title',
                        required=True)

    parser.add_argument('--scoreName', '-s',
                        help='Score name',
                        required=True)

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image. ',
                        type=argparse.FileType('w'),
                        required=True)

    parser.add_argument('--region',
                        help='Plot only this region. The format is '
                        'chr:start-end ',
                        required=True
                        )

    parser.add_argument('--zMax',
                        help='zMax',
                        type=float,
                        default=None)

    parser.add_argument('--vlines',
                        help='Genomic cooordindates separated by space. E.g. '
                        '--vlines 150000 3000000 124838433 ',
                        type=int,
                        nargs='+'
                        )

    parser.add_argument('--fontSize',
                        help='Font size for the labels of the plot',
                        type=float,
                        )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def plot_boundaries(ax, file_name, region):
    """
    Plots the boundaries as triangles in the given ax.

    :param ax:
    :param file_name: boundaries file
    :param region: chromosome region to plot
    :return:
    """
    chrom_region, start_region, end_region = region
    try:
        file_h = open(file_name, 'r')
    except IOError:
        sys.stderr.write("Boundaries file not found:\n{}".format(file_name))
        return

    prev_start = -1
    prev_chrom = None
    prev_line = None
    x = []
    y = []
    for line in file_h.readlines():
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = None
        start = None
        end = None
        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(detail)
            exit(msg)

        start = int(start)
        end = int(end)
        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

        assert start > 0, \
            "negative values found in bed file in line {}".format(line)

        if chrom == chrom_region and end > start_region and \
                start < end_region:

            if prev_start is None:
                # draw only half a triangle
                length = start - prev_start
                x1 = prev_start
                y1 = length
                x2 = start
                y2 = 0
                x.extend([x1, x2])
                y.extend([y1, y2])
            else:
                x1 = prev_start
                x2 = x1 + (start - prev_start) / 2
                x3 = start
                y1 = 0
                y2 = (start - prev_start)
                y3 = 0
                x.extend([x1, x2, x3])
                y.extend([y1, y2, y3])

            prev_start = start
            prev_line = line
            prev_chrom = chrom
    file_h.close()
    ax.plot(x, y,  color='black')


def plot_matrix(ax, label_ax, cbar_ax, matrix_properties, region):

    matrix_file = matrix_properties['file']
    hic_ma = HiCMatrix.hiCMatrix(matrix_file)
    chrom, region_start, region_end = region
    hic_ma.keepOnlyTheseChr(chrom)
    if 'show_masked_bins' in matrix_properties and \
            matrix_properties['show_masked_bins'] == 'yes':
        pass
    else:
        hic_ma.maskBins(hic_ma.nan_bins)
    if 'orientation' in matrix_properties and \
            matrix_properties['orientation'] == 'inverted':
        plot_inverted = True
    else:
        plot_inverted = False

    new_intervals = enlarge_bins(hic_ma.cut_intervals)
    hic_ma.interval_trees, hic_ma.chrBinBoundaries = \
        hic_ma.intervalListToIntervalTree(new_intervals)

    hic_ma.cut_intervals = new_intervals

    # expand region to plus depth on both sides
    # to avoid a 45 degree 'cut' on the edges

    # get bin id of start and end of region in given chromosome
    chr_start_id, chr_end_id = hic_ma.getChrBinRange(chrom)
    chr_start = hic_ma.cut_intervals[chr_start_id][1]
    chr_end = hic_ma.cut_intervals[chr_end_id-1][1]
    start_bp = max(chr_start, region_start - matrix_properties['depth'])
    end_bp = min(chr_end, region_end + matrix_properties['depth'])


    idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                           enumerate(hic_ma.cut_intervals)
                           if x[0] == chrom and x[1] >= start_bp
                           and x[2] <= end_bp])

    idx = idx[0:-1]
    # select only relevant matrix part
    hic_ma.matrix = hic_ma.matrix[idx, :][:, idx]
    try:
        new_nan_bins = hic_ma.nan_bins[np.in1d(hic_ma.nan_bins, idx)]
        hic_ma.nan_bins = new_nan_bins - idx[0]
    except:
        pass
    """
    print "filling matrix gaps..."
    hic_ma.matrix, _ = fill_gaps(hic_ma, depth_bins, False)
    """

    # fill the main diagonal, otherwise it looks
    # not so good. The main diagonal is filled
    # with an array containing the max value found
    # in the matrix
    if sum(hic_ma.matrix.diagonal()) == 0:
        print "filling main diagonal because is empty and " \
            "otherwise it looks bad..."
        max_value = hic_ma.matrix.data.max()
        main_diagonal = dia_matrix(([max_value]*hic_ma.matrix.shape[0], [0]),
                                shape=hic_ma.matrix.shape)
        hic_ma.matrix = hic_ma.matrix + main_diagonal
    """
    # code to select some meaningful max and min values
    min_zscore = -5
    max_zscore = None

    mean = np.mean(hic_ma.matrix.data)
    std = np.std(hic_ma.matrix.data)
    z_score = (hic_ma.matrix.data - mean) / std
    hic_ma.matrix.data[z_score < min_zscore] = 0
    if max_zscore is not None:
        min_ = hic_ma.matrix.data[z_score >= max_zscore].min()
        hic_ma.matrix.data[z_score >= max_zscore] = min_
        print "{}, {}".format(mean, min_)
    hic_ma.matrix.eliminate_zeros()
    """

    # select only the upper triangle of the matrix
    hic_ma.matrix = triu(hic_ma.matrix, k=0, format='csr')

    matrix = np.asarray(hic_ma.matrix.todense().astype(float))

    norm = None

    if 'transform' in matrix_properties:
        if matrix_properties['transform'] == 'log1p':
            matrix += 1
            norm = LogNorm()

        elif matrix_properties['transform'] == '-log':
            mask = matrix == 0
            matrix[mask] = matrix[mask == False].min()
            matrix = -1 * np.log(matrix)

    if 'max_value' in matrix_properties and matrix_properties['max_value'] != 'auto':
        vmax = matrix_properties['max_value']

    else:
        # try to use a 'aesthetically pleasant' max value
        vmax = np.percentile(matrix.diagonal(1), 80)

    if 'min_value' in matrix_properties and matrix_properties['min_value'] != 'auto':
        vmin = matrix_properties['min_value']
    else:
        bin_size = hic_ma.getBinSize()
        depth_bins = int(matrix_properties['depth'] / bin_size)
        vmin = np.median(matrix.diagonal(depth_bins))

    sys.stderr.write("setting min, max values to: {}, {}\n".format(vmin, vmax))
    cmap = cm.get_cmap(matrix_properties['colormap'])
    cmap.set_bad('white')

    img = pcolormesh_45deg(matrix, ax, start_pos, vmax=vmax,
                     vmin=vmin, cmap=cmap, norm=norm)

    img.set_rasterized(True)
    if plot_inverted:
        ax.set_ylim(matrix_properties['depth'], 0)
    else:
        ax.set_ylim(0, matrix_properties['depth'])

    # if a boundaries file is given, plot the
    # tad boundaries as line delineating the TAD triangles
    if 'boundaries_file' in matrix_properties:
        plot_boundaries(ax, matrix_properties['boundaries_file'], region)

    ax.set_xlim(region_start, region_end)
    if 'x labels' in matrix_properties and matrix_properties['x labels'] != 'no':
        ticks = ax.get_xticks()
        labels = ["{:.2f}".format((x / 1e6))
                  for x in ticks]
        labels[-1] = labels[-1] + "Mbp"
        ax.get_xaxis().set_tick_params(
            which='both',
            bottom='on',
            top='off',
            direction='out')

        ax.set_xticklabels(labels)
    else:
        ax.get_xaxis().set_tick_params(
            which='both',
            bottom='off',
            top='off',
            direction='out')
        ax.axes.get_xaxis().set_visible(False)

#    ax.xaxis.tick_top()

    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    cbar_ax.patch.set_alpha(0.0)
    cobar = plt.colorbar(img, ax=cbar_ax, fraction=0.95)
    cobar.solids.set_edgecolor("face")
    label_ax.text(0.3, 0.0, matrix_properties['title'],
            horizontalalignment='left', size='large',
            verticalalignment='bottom', transform=label_ax.transAxes)
    # plt.subplots_adjust(wspace=0, hspace=0.1, top=0.95,
    #                    bottom=0.05, left=0.1, right=0.5)


def plot_x_axis(ax, region, properties):
    chrom_region, region_start, region_end = region
    ax.set_xlim(region_start, region_end)
    ticks = ax.get_xticks()
    if ticks[-1] - ticks[1] <= 1e5:
        labels = ["{:.3f}".format((x / 1e3))
                  for x in ticks]
        labels[-1] = labels[-1] + "Kbp"

    elif 1e5 < ticks[-1] - ticks[1] < 4e6:
        labels = ["{:.3f}".format((x / 1e3))
                  for x in ticks]
        labels[-1] = labels[-1] + "Kbp"
    else:
        labels = ["{:.0f}Mbp".format((x / 1e6))
                  for x in ticks]
        #labels[-1] = labels[-1] + "Mbp"



    ax.axis["x"] = ax.new_floating_axis(0,0.5)
    
    ax.axis["x"].axis.set_ticklabels(labels)
    ax.axis['x'].axis.set_tick_params(which='minor',bottom='on')

    if 'fontsize' in properties:
        ax.axis["x"].major_ticklabels.set(size=int(properties['fontsize']))

    if 'where' in properties and properties['where']=='top':
        ax.axis["x"].set_axis_direction("top")


def pcolormesh_45deg(C, ax, start_pos_vector, vmin=None,
                     vmax=None, cmap=None, norm=None):
    """
    Turns the matrix 45 degrees and adjusts the
    bins to match the actual start end positions.
    """
    import itertools
    # code for rotating the image 45 degrees
    n = C.shape[0]
    # create rotation/scaling matrix
    t = np.array([[1, 0.5], [-1, 0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1], i[0])
                         for i in itertools.product(start_pos_vector[::-1],
                                                    start_pos_vector)]),t)
    # this is to convert the indices into bp ranges
    X = A[:, 1].reshape(n+1, n+1)
    Y = A[:, 0].reshape(n+1, n+1)
    # plot
    im = ax.pcolormesh(X, Y, np.flipud(C),
                       vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)
    return im


def fill_gaps(hic_ma, depth, fill_contiguous=False):
    """ try to fill in the gaps in the matrix by adding the average values of
    the neighboring rows and cols.

    This method produces best results when the
    missing data is not consecutive.

    """
    from scipy import sparse

    M, N = hic_ma.matrix.shape
    # get only a hic_ma.matrixtrix that contains only
    # the [deph:-depth] diagonals
    hic_ma.matrix = sparse.triu(hic_ma.matrix, k=-depth)
    hic_ma.matrix = sparse.tril(hic_ma.matrix, k=depth, format='csr')

    fill_ma = hic_ma.matrix.copy().tolil()
    if fill_contiguous is True:
        good_nans = hic_ma.nan_bins
        consecutive_nan_idx = np.array([])
    else:
        # find stretches of consecutive nans
        consecutive_nan_idx = np.flatnonzero(np.diff(hic_ma.nan_bins) == 1)
        # the banned list of indices is equal to the actual list
        # and the list plus one, to identify both consecutive nans
        consecutive_nan_idx = np.concatenate([consecutive_nan_idx,
                                              consecutive_nan_idx+1])
        # find the nans that are not consecutive
        good_nans = [x for idx, x in enumerate(hic_ma.nan_bins)
                     if idx not in consecutive_nan_idx]

    for missing_bin in good_nans:
        if 0 < missing_bin < M - 1:
            # the new row value is the mean between the upper
            # and lower rows
            fill_ma[missing_bin, :] = (hic_ma.matrix[missing_bin - 1, :] +
                                       hic_ma.matrix[missing_bin + 1, :]) / 2

            # same for cols
            fill_ma[:, missing_bin] = (hic_ma.matrix[:, missing_bin - 1] +
                                       hic_ma.matrix[:, missing_bin + 1]) / 2

    # identify the intersection points of the failed regions because their
    # neighbors get wrong values
    for bin_a in good_nans:
        for bin_b in good_nans:
            if 0 < bin_a < M and \
                    0 < bin_b < M:
                # the fill value is the average over the
                # neighbors that do have a value

                fill_value = np.mean([
                    hic_ma.matrix[bin_a-1, bin_b-1],
                    hic_ma.matrix[bin_a-1, bin_b+1],
                    hic_ma.matrix[bin_a+1, bin_b-1],
                    hic_ma.matrix[bin_a+1, bin_b+1],
                    ])

                fill_ma[bin_a, bin_b] = fill_value

    # return the matrix and the bins that continue to be nan
    return fill_ma.tocsr(), np.sort(consecutive_nan_idx)


def plot_hic_arcs(ax, label_ax, properties, region):

    """
    Makes and arc connecting two points on a linear scale representing
    interactions between Hi-C bins.
    :param ax: matplotlib axis
    :param label_ax: matplotlib axis for labels
    :param region: tuple containing (chrom, start, end) that will limit the view to
    only the genomic positions between start and end

    :param properties: Dictionary of values that should include:
        at least the file. Optionally: lower threshold, value to filter out low scoring
        enrichments, alpha: (transparency) value for arc lines,  color: default color
        for arcs unless specified in the intervals
        intervals: list of the format chr1:1000000:200000:red
         chr3:13000000:14000000:blue. Only the arcs starting or ending on
         such regions will be plotted. The last field field of the intervals used
         to set the color of the arcs.
    """
    from matplotlib.colors import colorConverter
    from matplotlib.patches import Arc

    hic_matrix = HiCMatrix.hiCMatrix(properties['file'])
    chrom, region_start, region_end = region
    hic_matrix.keepOnlyTheseChr(chrom)
    hic_matrix.diagflat()
    hic_matrix.keepOnlyTheseChr(chrom)

    # filter low scoring enrichments
    hic_matrix.matrix.data[np.isnan(hic_matrix.matrix.data)] = 0
    if 'lower threshold' not in properties:
        sys.stderr.write("Plotting arcs without a lower threshold. This can result "
                         "in too many interactions being plotted. Consider adding a "
                         "'lower threshold' value to the configuration file.")
    else:
        try:
            lower_threshold = float(properties['lower threshold'])
        except ValueError:
            sys.exit("lower threshold value is invalid: {} for {}".format(properties['lower threshold'],
                                                                          properties['file']))
        hic_matrix.matrix.data[hic_matrix.matrix.data < lower_threshold] = 0

    hic_matrix.matrix.eliminate_zeros()
    mat = triu(hic_matrix.matrix, k=0, format='coo')
    max_radius = 0
    count = 0

    if properties['color']:
        color_main = properties['color']
    else:
        color_main = 'blue'

    if properties['alpha']:
        alpha = float(properties['alpha'])
    else:
        alpha = 0.8

    # expand view point
    vp_chrom, vp_start, vp_end = properties['view point'].split(':')
    if vp_chrom == chrom:
        vp_intval = IntervalTree()
        vp_intval.insert_interval(Interval(int(vp_start), int(vp_end)))
    else:
        vp_intval = None
    # process intervals, if any, to plot
    # only interactions from those regions
    intervals = None
    if 'highlights' in properties:
        intervals = {chrom: IntervalTree()}
        if properties['highlights'].endswith(".bed"):
            # the highlights are a bed file
            with open(properties['highlights']) as f:
                for line in f.readlines():
                    if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')

                    try:
                        chrom_vp, start, end = fields[0:3]
                    except Exception as detail:
                        msg = "Error reading line: {}\nError message: {}".format(detail)
                        exit(msg)
                    # skip regions not in the same chromosome
                    if chrom_vp != chrom:
                        continue
                    start = int(start)
                    end = int(end)
                    # insert the interval only if it does not overlap with the view point
                    if vp_intval is not None:
                        if len(vp_intval.find(start, end+1)) > 0:
                            sys.stderr.write("skipping region\n".format(start, end))
                            continue

                    intervals[chrom_vp].insert_interval(Interval(start, end, value='red'))
        else:
            for rangebp in properties['highlights']:
                fields = rangebp.split(":")
                chrom_vp = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if len(fields) == 4:
                    color = fields[3]
                else:
                    color = 'blue'
                print (chrom, start, end, color)
                intervals[chrom_vp].insert_interval(Interval(start, end, value=color))

    for idx, (row, col) in enumerate(np.vstack([mat.row, mat.col]).T):
        chrom1, start1, end1, _ = hic_matrix.cut_intervals[row]
        chrom2, start2, end2, _ = hic_matrix.cut_intervals[col]
        if intervals:
            match1 = intervals[chrom1].find(start1, end1)
            match2 = intervals[chrom2].find(start2, end2)
            if len(match1) == 0 and len(match2) == 0:
                #continue
                color = color_main
            else:
                if len(match1) > 0:
                    color = colorConverter.to_rgba(match1[0].value, alpha=alpha)
                else:
                    color = colorConverter.to_rgba(match2[0].value, alpha=alpha)

        center = start1 + float(start2 - start1) / 2
        radius = start2 - start1
        if radius > max_radius:
            max_radius = radius
        count += 1
        ax.plot([center], [radius])
        if 'line width' in properties:
            line_width = float(properties['line width'])
        else:
            line_width = 0.5*np.sqrt(mat.data[idx])
        ax.add_patch(Arc((center, 0), radius,
                         radius*2, 0, 0, 180, color=color, lw=line_width))

    print "{} arcs plotted".format(count)
    chrom, region_start, region_end = region
    if 'orientation' in properties and properties['orientation'] == 'inverted':
        ax.set_ylim(region_end, -1)
    else:
        ax.set_ylim(-1, region_end)
    ax.set_xlim(region_start, region_end)
    ax.axis('off')

    label_ax.text(0.3, 0.0, properties['title'],
            horizontalalignment='left', size='large',
            verticalalignment='bottom', transform=label_ax.transAxes)

def plot_interactions(ax, matrix_properties, region, args):
    """
    Plots arrows from a viewpoint to nearby
    positions that are in 'regions'
    """
    matrix_file = matrix_properties['file']
    hic_ma = HiCMatrix.hiCMatrix(matrix_file)
    bin_size = hic_ma.getBinSize()
    chrom, region_start, region_end = region
    hic_ma.keepOnlyTheseChr(chrom)

    # hard coded threshold
    pval_threshold = float(matrix_properties['extra'][1])
    # hard coded view point (rox2)
    viewpoint = (matrix_properties['extra'][2],
                 int(matrix_properties['extra'][3]),
                 int(matrix_properties['extra'][3])+1)
    viewpoint_bin = hic_ma.getRegionBinRange(*viewpoint)

    # extend view point to include left and right bins
    mat = hic_ma.matrix[viewpoint_bin[0]-1:viewpoint_bin[0]+2, :]
    mat.data[mat.data <= -np.log(pval_threshold)] = 0
    mat.data[np.isnan(mat.data)] = 0
    mat.eliminate_zeros()
    mat = mat.tocoo()

    for target in np.unique(mat.col):
        t_chrom, t_start, t_end, _ = hic_ma.cut_intervals[target]
        x_pos = t_end - (t_end - t_start)/2
        rad = '-0.1' if x_pos < viewpoint[1] else '0.1'
        ax.annotate(" ", xy=(x_pos, 0), xytext=(viewpoint[1], 0),
                    xycoords='data', color='black',
                    arrowprops=dict(arrowstyle="simple", lw=10,
                                    facecolor='black', edgecolor='none',
                                    connectionstyle="arc3,rad={}".format(rad))
                    )
    ax.set_ylim(-5, 0)
    ax.set_xlim(region_start, region_end)

    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)


def get_gene_parts(fields, bed_properties):
    chrom, start, end = fields[0:3]
    start = int(start)
    end = int(end)
    try:
        name = fields[3]
    except:
        name = ''
    try:
        strand = fields[5]
    except:
        strand = "."
    try:
        thick_start = int(fields[6])
    except:
        thick_start = int(fields[1])
    try:
        thick_end = int(fields[7])
    except:
        thick_end = int(fields[2])

    try:
        rgb = fields[8]
        rgb = rgb.split(',')
        if len(rgb) == 3:
            rgb = [float(x)/255 for x in rgb]
            edgecolor = 'black'
        else:
            rgb = bed_properties['color']
            edgecolor = bed_properties['color']
    except IndexError:
        rgb = bed_properties['color']
        edgecolor = bed_properties['color']
    try:
        block_count = int(fields[9])
        block_sizes = [int(x) for x in fields[10].split(',') if x.strip() != '']
        block_starts = [int(x) for x in fields[11].split(',') if x.strip() != '']
#    except IndexError, ValueError:
    except:
        block_count = 0
        block_sizes = []
        block_starts = []

    return chrom, start, end, name, strand, thick_start, thick_end, rgb, \
           block_count, block_sizes, block_starts

def draw_gene_simple(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor):
    """
    draws a gene using different styles
    """
    (chrom, start, end, name, strand, thick_start, thick_end, rgb_, block_count,
     block_sizes, block_starts) = get_gene_parts(fields, bed_properties)
    # prepare the vertices for all elements that will be drawn
    # the region length without the tip of the end arrow whose length is 'small_relative'
    if strand == '+':
        x0  = start
        x1  = end - small_relative
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        -----------------\
        ---------------- /

        """

        vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + small_relative, y0 + 50 ), (x1, y0)]

    elif strand == '-':
        x0  = start + small_relative
        x1  = end
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        /--------___---------_
        \--------   ----------
        """
        vertices = [(x0, y0), (x0 -small_relative, y0 + 50), (x0, y1), (x1, y1), (x1, y0)]

    ax.add_patch(matplotlib.patches.Polygon(vertices, closed=True, fill=True,
                                            edgecolor='black',
                                            facecolor=rgb))

    center = start + float(end - start)/2
    ax.text(center, ypos + 125, fields[3], size='small',
            horizontalalignment='center', verticalalignment='top')


def draw_gene_with_introns(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor,
                           fontproperties=None):
    """
    draws a gene using different styles
    """
    (chrom, start, end, name, strand, thick_start, thick_end, rgb_, block_count,
     block_sizes, block_starts) = get_gene_parts(fields, bed_properties)
    # draw a line from the start until the end of the gene
    if block_count == 0 and thick_start == start and thick_end == end:
        draw_gene_simple(ax, fields, ypos,  bed_properties, small_relative, rgb, edgecolor)
        return

    ax.plot([start, end], [ypos + 50, ypos + 50], 'black', linewidth=0.5, zorder = -1)
    if strand == '+':
        x1  = thick_start
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        -----------------\
        -----------------/

        """
        start_box = [(start, y0), (start, y1), (x1, y1), (x1, y0)]
        end_box = [(thick_end, y0), (thick_end, y1), (end - small_relative, y1), (end, y0 + 50), (end-small_relative, y0)]

    elif strand == '-':
        x0  = start + min(small_relative, thick_start - start)
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        /--------___---------_
        \--------   ----------
        """
        start_box = [(x0, y0), (start, y0 + 50), (x0, y1), (thick_start, y1), (thick_start, y0)]
        end_box = [(thick_end, y0), (thick_end, y1), (end, y1),  (end, y0)]


    for idx in range(0, block_count):
        x0 = start + block_starts[idx]
        x1 = x0 + block_sizes[idx]
        if x1 < thick_start or x0 > thick_end:
            y0 = ypos + 25
            y1 = ypos + 75
        else:
            y0 = ypos
            y1 = ypos + 100

        if x0 < thick_start < x1:
            vertices = ([(x0, ypos+25), (x0, ypos+75), (thick_start, ypos+75), (thick_start, ypos+100),
                         (thick_start, ypos+100), (x1, ypos+100), (x1, ypos), (thick_start, ypos), (thick_start, ypos+25)])

        elif x0 < thick_end < x1:
            vertices = ([(x0, ypos), (x0, ypos+100), (thick_end, ypos+100), (thick_end, ypos+75),
                         (x1, ypos+75), (x1, ypos+25), (thick_end, ypos+25), (thick_end, ypos)])
        else:
            vertices = ([(x0,y0), (x0, y1), (x1, y1), (x1, y0)])

        ax.add_patch(matplotlib.patches.Polygon(vertices, closed=True, fill=True,
                                                linewidth=0.1,
                                                edgecolor='none',
                                                facecolor=rgb))

        if idx < block_count - 1:
            intron_length =  block_starts[idx+1]  - (block_starts[idx] + block_sizes[idx])
            marker = 5 if strand == '+' else 4
            if intron_length > 3*small_relative:
                pos = np.arange(x1 + 1*small_relative, x1 + intron_length +small_relative, int(2*small_relative))
                ax.plot(pos, np.zeros(len(pos))+ypos+50, '.', marker=marker,
                    fillstyle='none', color='blue', markersize = 3)

            elif intron_length > small_relative:
                intron_center = x1 + int(intron_length)/2
                ax.plot([intron_center], [ypos+50], '.', marker=5,
                    fillstyle='none', color='blue', markersize = 3)


    center = start + float(end - start)/2
    ax.text(center, ypos + 125, fields[3],
            horizontalalignment='center', verticalalignment='top', fontproperties=fontproperties)


def plot_bed(ax, label_ax, bed_properties, region):
    file_h = open(bed_properties['file'], 'r')
    chrom_region, start_region, end_region = region
    counter = 0
    small_relative = 0.003 * (end_region-start_region)
    prev_start = -1
    prev_line = None
    prev_chrom = None
    max_num_row = 1
    region_intervals = IntervalTree()
    colormap = None
    # check if the color given is a color map
    from matplotlib import cm
    color_options = [m for m in cm.datad]

    ax.set_frame_on(False)
    # to improve the visualization of the genes
    # it is good to have an estimation of the label
    # length. In the following code I try to get
    # the length of a 'w'.

    if 'fontsize' in bed_properties:
        from matplotlib import font_manager
        fp = font_manager.FontProperties(size=bed_properties['fontsize'])
    else:
        fp = None

    if 'type' in bed_properties and \
        bed_properties['type'] == 'genes':
        t = matplotlib.textpath.TextPath((0,0), 'w', prop=fp)
        len_w = t.get_extents().width * 300
    else:
        len_w = 1

    if bed_properties['color'] in color_options:
        import matplotlib as mpl
        import matplotlib.cm as cm
        norm = mpl.colors.Normalize(vmin=bed_properties['min_value'],
                                    vmax=bed_properties['max_value'])
        cmap = cm.get_cmap(bed_properties['color'])
        colormap  = cm.ScalarMappable(norm=norm, cmap=cmap)

    for line in file_h.readlines():
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')

        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(detail)
            exit(msg)
        start = int(start)
        end = int(end)
        if prev_chrom is None:
            prev_chrom = chrom
        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}\n{} ".format(prev_line, line)

        prev_start = start
        prev_line = line
        assert start > 0, \
            "negative values found in bed file in line {}".format(line)

        if chrom == chrom_region and end > start_region and \
                start < end_region:
            # Rectangle(xy, width, height, angle=0.0, **kwargs)
            counter += 1
            try:
                strand = fields[5]
            except:
                strand = "."
            try:
                rgb = fields[8]
                rgb = rgb.split(',')
                rgb[2]  # check that rgb has three elements
                rgb = [float(x)/255 for x in rgb]
                edgecolor = bed_properties['color']
            except IndexError:
                if colormap:
                    rgb = colormap.to_rgba(float(fields[4]))
                    edgecolor = colormap.to_rgba(float(fields[4]))
                else:
                    rgb = bed_properties['color']
                    edgecolor = bed_properties['color']

            except Exception as e:
                exit("Error occurred: {}".format(e))

            if 'type' in bed_properties and \
                    bed_properties['type'] == 'domain':
                ypos = 100 if counter % 2 == 0 else 1
                ax.add_patch(matplotlib.patches.Rectangle(
                        (start, ypos),
                        end-start,
                        100, edgecolor='black',
                        facecolor=bed_properties['color']))

            # check for overlapping features
            match = region_intervals.find(start, end)
            if len(match) == 0:
                min_free_row = 0
            else:
                rows_used = np.zeros(max_num_row + 2)
                for x in match:
                    rows_used[x.value] = 1
                min_free_row = min(np.flatnonzero(rows_used==0))

            if 'type' in bed_properties and \
                bed_properties['type'] == 'genes' and \
                end - start < len(fields[3]) * len_w:
                region_intervals.add_interval(Interval(start, start + (len(fields[3]) * len_w), min_free_row))
            else:
                region_intervals.add_interval(Interval(start, end+small_relative, min_free_row))
            if min_free_row > max_num_row:
                max_num_row = min_free_row
            if min_free_row > 0:
                # this means we are plotting 
                # on top of previous genes
                ypos = min_free_row * 230
            else:
                ypos = 0

            if 'display' in bed_properties and \
                    bed_properties['display'] == 'collapsed':
                ypos = 0

            # give direction to genes
            if 'type' in bed_properties and \
                    bed_properties['type'] == 'genes' and \
                    strand in ['+','-']:
#                draw_gene_simple(ax, fields, ypos, bed_properties, small_relative)
                draw_gene_with_introns(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor,
                                       fontproperties=fp)

            else:
                ax.add_patch(matplotlib.patches.Rectangle(
                        (start, ypos), end-start,
                        100, edgecolor=edgecolor,
                        facecolor=rgb))

    if 'type' in bed_properties and \
            bed_properties['type'] == 'domain':
        ax.set_ylim(-5, 205)
    elif 'display' in bed_properties and \
        bed_properties['display'] == 'collapsed':
        ax.set_ylim(-5, 105)
    else:
        ax.set_ylim((max_num_row+1)*230, -25)

    ax.set_xlim(region[1], region[2])

    label_ax.text(0.15, 1.0, bed_properties['title'],
            horizontalalignment='left', size='large',
            verticalalignment='top', transform=label_ax.transAxes)


def plot_bedgraph(ax, label_ax, bedgraph_properties, region):
    file_h = open(bedgraph_properties['file'], 'r')
    chrom_region, start_region, end_region = region
    score_list = []
    pos_list = []

    for line in file_h.readlines():
        if line.startswith('browser') or line.startswith('track'):
            continue
        chrom, start, end, score = line.split('\t')
        start = int(start)
        end = int(end)
        if chrom == chrom_region and start_region -100000 <= start and \
                end_region + 100000 >= end:
            score_list.append(float(score))
            pos_list.append(start + (end - start)/2)

    if 'extra' in bedgraph_properties and \
            bedgraph_properties['extra'][0] == '4C':
        # draw a vertical line for each fragment region center
        ax.fill_between(pos_list, score_list,
                        facecolor=bedgraph_properties['color'],
                        edgecolor='none')
        ax.vlines(pos_list, [0], score_list, color='olive', linewidth=0.5)
        ax.plot(pos_list, score_list, '-', color='slateblue', linewidth=0.7)
    else:
        try:
            ax.fill_between(pos_list, score_list,
                            facecolor=bedgraph_properties['color'])
        except ValueError:
            exit("Invalid color {} for {}".format(bedgraph_properties['color'], bedgraph_properties['file']))
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim(region[1], region[2])

    ymin, ymax = ax.get_ylim()
    if 'max_value' in bedgraph_properties and bedgraph_properties['max_value'] != 'auto':
        ymax = bedgraph_properties['max_value']
    if 'min_value' in bedgraph_properties and bedgraph_properties['min_value'] != 'auto':
        ymin = bedgraph_properties['min_value']

    if float(ymax) % 1 == 0:
        ymax_print = int(ymax)
    else:
        ymax_print = "{:.1f}".format(ymax)
    ax.set_ylim(ymin, ymax)
    ydelta = ymax - ymin
    small_x = 0.01 * (end_region - start_region)

    if 'show data range' in bedgraph_properties and \
            bedgraph_properties['show data range'] == 'no':
        pass
    else:
        # by default show the data range
        ax.text(start_region-small_x, ymax - ydelta * 0.2,
                "[{}-{}]".format(ymin, ymax_print),
                horizontalalignment='left', size='small',
                verticalalignment='bottom')

    label_ax.text(0.15, 0, bedgraph_properties['title'],
                horizontalalignment='left', size='large',
                verticalalignment='bottom', transform=label_ax.transAxes)

    """
    if 'extra' in bedgraph_properties :

        ticks = ax.get_xticks()
        labels = ["{:.1f} Mbp".format((x / 1e6))
                  for x in ticks]

        ax.set_xticklabels(labels, size='large')
        ax.axes.get_xaxis().set_visible(True)
    """

def plot_bedgraph_matrix(ax, label_ax, properties, region):
    """
    Plots a bedgraph matrix file, that instead of having
    a single value per bin, it has several values.

    :param ax:
    :param label_ax:
    :param properties:
    :param region:
    :return:
    """

    fh = open(properties['file'], 'r')
    chrom_region, start_region, end_region = region
    start_pos = []
    matrix_rows = []
    for line in fh:
        line = line.strip()
        region = line.split('\t')
        chrom = region[0]
        start = int(region[1])
        end = int(region[2])
        if chrom == chrom_region and start_region -100000 <= start and \
            end_region + 100000 >= end:
            start_pos.append(start)
            matrix_rows.append(np.fromiter(region[3:], np.float))

    matrix = np.vstack(matrix_rows).T
    if 'orientation' in properties and \
        properties['orientation'] == 'inverted':
        matrix = np.flipud(matrix)

    vmin = None
    vmax = None
    if 'max_value' in properties and properties['max_value'] != 'auto':
        vmax = properties['max_value']


    if 'min_value' in properties and properties['min_value'] != 'auto':
        vmin = properties['min_value']


    X,Y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
    img = ax.pcolormesh(X, Y, matrix, vmin=vmin, vmax=vmax, shading='gouraud')
    img.set_rasterized(True)
    ax.set_xlim(start_region, end_region)
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    label_ax.text(0.15, 0, properties['title'],
            horizontalalignment='left', size='large',
            verticalalignment='bottom', transform=label_ax.transAxes)


def plot_bigwig(ax, label_ax, bigwig_properties, region):

    from bx.bbi.bigwig_file import BigWigFile
    bw = BigWigFile(open(bigwig_properties['file'], 'r'))
    chrom, region_start, region_end = region
    # compute the score in bins of 10000 SLOW
#    score = np.array([bw.query(region[0], x, x+10000,1)[0]['mean']
#                      for x in range(region[1], region[2], 10000)])

    num_bins = 700
    if 'number of bins' in bigwig_properties:
        try:
            num_bins = int(bigwig_properties['number of bins'])
        except TypeError:
            exit("'number of bins' value: {} for bigwig file {} "
                 "is not valid.".format(bigwig_properties['number of bins'],
                                        bigwig_properties['file']))

    scores = bw.get_as_array(chrom, region_start, region_end)
    if 'nans to zeros' in bigwig_properties and bigwig_properties['nans to zeros'] is True:
        scores[np.isnan(scores)] = 0

    scores = np.ma.masked_invalid(scores)

    if region_end - region_start < 2e6:
        if scores is None:
            chrom = chrom.replace('chr', '')
            scores = np.ma.masked_invalid(bw.get_as_array(chrom, region_start, region_end))
        if scores is None:
            sys.stderr.write("could not find values for region {}\n".format(region))

        else:
            lins = np.linspace(0, len(scores), num_bins).astype(int)
            scores_per_bin = [np.mean(scores[lins[x]:lins[x+1]]) for x in range(len(lins)-1)]
            _x = lins + region_start
            x_values = [float(_x[x] + _x[x+1])/2 for x in range(len(lins)-1)]
            ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                            color=bigwig_properties['color'],
                            facecolor=bigwig_properties['color'])

    else:
        # this method produces shifted regions. It is not clear to me why this happens.
        # Thus I only activate the faster but shifted method for large regions
        # when the previous method would be to slow
        score = bw.query(chrom, region_start, region_end, num_bins)
        if score is None:
            chrom = chrom.replace('chr', '')
            score = bw.query(chrom, region_start, region_end, num_bins)
        if score is None:
            sys.stderr.write("could not find values for region {}\n".format(region))
        else:
            score = [x['mean'] for x in score]
            x_values = np.linspace(region_start, region_end, num_bins)
            ax.fill_between(x_values, score, linewidth=0.1,
                            color=bigwig_properties['color'],
                            facecolor=bigwig_properties['color'])


    ax.set_xlim(region[1], region[2])
    ymin, ymax = ax.get_ylim()
    if 'max_value' in bigwig_properties and ['max_value'] != 'auto':
        ymax = bigwig_properties['max_value']
    if 'min_value' in bigwig_properties and bigwig_properties['min_value'] != 'auto':
        ymin = bigwig_properties['min_value']

    if 'orientation' in bigwig_properties and \
        bigwig_properties['orientation'] == 'inverted':

        ax.set_ylim(ymax, ymin)
    else:
        ax.set_ylim(ymin, ymax)

#    ax.set_yticks([ymax])
    ydelta = ymax - ymin

#    ax.set_yticklabels(["{}-{}".format(int(ymin), int(ymax))], size='large')
    # set min max
    if float(ymax) % 1 == 0:
        ymax_print = int(ymax)
    else:
        ymax_print = "{:.1f}".format(ymax)
    small_x = 0.01 * (region_end - region_start)
    if 'show data range' in bigwig_properties and \
        bigwig_properties['show data range'] == 'no':
        pass
    else:
        # by default show the data range
        ax.text(region_start-small_x, ymax - ydelta * 0.2,
                "[{}-{}]".format(int(ymin), ymax_print),
                horizontalalignment='left', size='small',
                verticalalignment='bottom')

    """
    ax.text(region_end, ymax - ydelta * 0.2, bigwig_properties['title'],
            horizontalalignment='right', size='large',
            verticalalignment='bottom')

    """
    label_ax.text(0.15, 0, bigwig_properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom')
                  #transform=label_ax.transAxes)


def get_region(region_string):
    """
    splits a region string into
    a chrom, start_region, end_region tuple
    The region_string format is chr:start-end
    """
    if region_string:
        region_string = region_string.translate(
            None, ",.;|!{}()").replace("-", ":")
        region = region_string.split(":")
        chrom = region[0]
        try:
            region_start = int(region[1])
        except IndexError:
            region_start = 0
        try:
            region_end = int(region[2])
        except IndexError:
            region_end = 1e15  # a huge number
        region_string = [chrom, region_start, region_end]
    return chrom, region_start, region_end

def parse_tracks(tracks_file):
    """
    Parses a configuration file

    :param tracks_file: file path containing the track configuration
    :return: array of dictionaries. Each
    """
    from ConfigParser import SafeConfigParser
    from ast import literal_eval
    parser = SafeConfigParser()
    parser.readfp(tracks_file)

    track_list = []
    for section_name in parser.sections():
        track_options = dict()
        if section_name in ['spacer', 'x-axis']:
            track_options[section_name] = True
        for name, value in parser.items(section_name):
            if name in ['max_value', 'min_value', 'depth', 'width'] and value != 'auto':
                track_options[name] =  literal_eval(value)
            else:
                track_options[name] =  value

        track_list.append(track_options)

    return track_list

def main():

    args = parseArguments().parse_args()

    region = get_region(args.region)
    chrom, region_start, region_end = region
    if region_end <= region_start:
        exit("Please check that the region end is larger than the region start.\n"
             "Values given:\nstart: {}\nend: {}\n".format(region_start, region_end))

    track_properties = parse_tracks(args.tracks)

    # prepare layout based on the tracks given.
    # The main purpose of the following loop is
    # to get the height of each of the tracks


    track_height = []
    for track_dict in track_properties:
        if 'width' in track_dict:
            track_height.append(track_dict['width'])
        elif 'depth' in track_dict:
            height = (track_dict['depth'] * args.width /
                      (1.8*(region_end - region_start)))
            track_height.append(height)
        else:
            track_height.append(0.3)

    if args.height:
        fig_height = args.height
    else:
        fig_height = sum(track_height)

    print (args.width, fig_height)
    fig = plt.figure(figsize=(args.width, fig_height))
    fig.suptitle(args.title)

    if args.fontSize:
        fontsize = args.fontSize
    else:
        fontsize = float(args.width) * 0.7

    font = {'size': fontsize}
    matplotlib.rc('font', **font)

    grids = gridspec.GridSpec(len(track_properties), 2,
                              height_ratios=track_height,
                              width_ratios=[1, 0.05])

    # iterate again to plot each track
    axis_list = []
    for idx, properties in enumerate(track_properties):
        if 'spacer' in properties:
            continue
        axis = axisartist.Subplot(fig, grids[idx, 0])
        fig.add_subplot(axis)
        axis.axis[:].set_visible(False)

        if 'x-axis' in properties:
            # ideally the axisartis would allow
            # to have a floating axis but this is not
            # working
            plot_x_axis(axis, region, properties)
            continue
        else:
            label_axis = plt.subplot(grids[idx, 1])
            label_axis.set_axis_off()

        if properties['file'].endswith('.bed'):
            plot_bed(axis, label_axis, properties, region)
        elif properties['file'].endswith('.bg'):
            plot_bedgraph(axis, label_axis, properties, region)
        elif properties['file'].endswith('.bw'):
            plot_bigwig(axis, label_axis, properties, region)
        elif properties['file'].endswith('.npz'):
            if 'type' in properties:
                if properties['type'] == 'interaction':
                    plot_interactions(axis, properties, region, args)
                elif properties['type'] == 'arcplot':
                    plot_hic_arcs(axis, label_axis, properties, region)
                else:
                    exit("Hi-C plot type invalid ({}) for {}".format(properties['type'],
                                                                     properties['file']))
            else:
                # to avoid the color bar to span all the
                # width of the axis I pass two axes
                # to plot_matrix
                cbar_axis = label_axis
                label_axis = plt.subplot(grids[idx, 1])
                label_axis.set_axis_off()
                plot_matrix(axis, label_axis, cbar_axis,
                            properties, region)
        elif properties['file'].endswith('.bm'):
            plot_bedgraph_matrix(axis, label_axis, properties, region)
        axis_list.append(axis)
    if args.vlines:
        from matplotlib.patches import ConnectionPatch
        a_ymax = axis_list[0].get_ylim()[1]
        b_ymin = axis_list[-1].get_ylim()[0]
        for start_pos in args.vlines:
            con = ConnectionPatch(xyA=(start_pos, a_ymax),
                                  xyB=(start_pos, b_ymin),
                                  coordsA="data", coordsB="data",
                                  axesA=axis_list[0],
                                  axesB=axis_list[-1],
                                  arrowstyle="-",
                                  linestyle='dashed',
                                  linewidth=0.5,
                                  zorder=10)
            axis_list[0].add_artist(con)

    plt.subplots_adjust(wspace=0, hspace=0.1, top=0.9,
                        bottom=0.12, left=0.04, right=0.92)

#    plt.tight_layout()
    plt.savefig(args.outFileName.name, dpi=300)
