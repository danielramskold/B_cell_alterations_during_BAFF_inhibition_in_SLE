from __future__ import division
import argparse, csv, os, random, numpy
from collections import defaultdict
import pylab as pyplot

def CMY_to_hexcode(C, M, Y):
	r = 1-C
	g = 1-M
	b = 1-Y
	if any(not 0<=v<=1 for v in (r, g, b)): print C, M, Y, r, g, b
	return '#' + ''.join(['%02x'%round(255*v) for v in (r, g, b)])

filename_to_col = [defaultdict(random.random) for i in range(3)]

def cell_to_value(row, cm, ci):
	if cm is None:
		return 0
	elif cm == 'file' or cm == 'population':
		if o.highlight_cluster is None:
			return filename_to_col[ci][row[cm]]
		elif row[cm] in o.highlight_cluster:
			return filename_to_col[ci][row[cm]] / 3 + 0.67
		else:
			return filename_to_col[ci][row[cm]] / 3
	elif ',' in cm:
		return max(float(row['*'+marker]) for marker in cm.split(','))
	else:
		return float(row['*'+cm])

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('accense_2DSNE', nargs='?')
	parser.add_argument('accense_csv')
	parser.add_argument('-c', '--cyan', metavar='marker', help='if "file" the channel will be random')
	parser.add_argument('-m', '--magenta', metavar='marker', help='if "file" the channel will be random')
	parser.add_argument('-y', '--yellow', metavar='marker', help='if "file" the channel will be random')
	parser.add_argument('-o', '--plot', default='tsne.pdf')
	parser.add_argument('--marker', default='.', help='dot shape, e.g. "o"')
	parser.add_argument('--minint', default=0.2, type=float)
	parser.add_argument('--zeroint', default=0.05, type=float)
	parser.add_argument('-s', '--samples', nargs='+')
	parser.add_argument('--rescale', type=float, default=1)
	parser.add_argument('--highlight_cluster', nargs='+')
	parser.add_argument('--name_clusters', nargs='?', type=int, help='takes optional font size', const=12)
	o = parser.parse_args()
	
	random.seed(100)
	
	coordinates_ri = list()
	if o.accense_2DSNE is None:
		with open(o.accense_csv, 'rU') as infh:
			for ri, row in enumerate(csv.DictReader(infh)):
				coordinates_ri.append((float(row['SNEx'])*o.rescale, float(row['SNEy'])*o.rescale))
	else:
		with open(o.accense_2DSNE, 'rU') as infh:
			for ri, row in enumerate(csv.DictReader(infh)):
				coordinates_ri.append((float(row['y1'])*o.rescale, float(row['y2'])*o.rescale))
	
	samples = None if o.samples is None else set(o.samples)
	
	skip_ri = set()
	levels_ri = list()
	cluster_ri = defaultdict(set)
	with open(o.accense_csv, 'rU') as infh:
		for ri, row in enumerate(csv.DictReader(infh)):
			#levels_thisrow = [0 if cm is None else float(row['*'+cm]) for ci, cm in enumerate((o.cyan, o.magenta, o.yellow))]
			levels_thisrow = [cell_to_value(row, cm, ci) for ci, cm in enumerate((o.cyan, o.magenta, o.yellow))]
			levels_ri.append(levels_thisrow)
			if o.name_clusters: cluster_ri[row['population']].add(ri)
			
			if o.samples is not None:
				if os.path.split(row['file'])[-1] not in samples:
					skip_ri.add(ri)
					
	
	maxlevels = [max(zip(*levels_ri)[mi])+1 for mi in range(3)]
	
	minintensity = o.minint
	for ri in range(len(coordinates_ri)):
		if ri in skip_ri: continue
		CMY = []
		for mi, marker in enumerate((o.cyan, o.magenta, o.yellow)):
			if marker is None: CMY.append(0)
			else:
				maxlvl = maxlevels[mi]
				level = levels_ri[ri][mi]
				CMY.append(o.zeroint if level == 0 else level/maxlvl * (1-minintensity) + minintensity)
		x,y = coordinates_ri[ri]
		c = CMY_to_hexcode(*CMY)
		pyplot.plot([x], [y], o.marker, color=c)
	title = ', '.join(['%s: %s'%(c,m) for m,c in ((o.cyan, 'cyan'), (o.magenta, 'magenta'), (o.yellow, 'yellow'))])
	
	if o.name_clusters is not None:
		for cluster_name, ri_arr in cluster_ri.items():
			x_mid = numpy.mean([coordinates_ri[ri][0] for ri in ri_arr])
			y_mid = numpy.mean([coordinates_ri[ri][1] for ri in ri_arr])
			if o.highlight_cluster and cluster_name not in o.highlight_cluster:
				textcolour = '#555555'
			else:
				textcolour = '#000000'
			pyplot.text(x_mid, y_mid, cluster_name, horizontalalignment='center', verticalalignment='center', fontdict={'size':o.name_clusters, 'family':'Arial', 'color':textcolour})
			
	
	pyplot.title(title)
	pyplot.xlim(-50, 50)
	pyplot.ylim(-50, 50)
	pyplot.savefig(o.plot)
