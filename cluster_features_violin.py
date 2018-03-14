from __future__ import division, print_function
import argparse, csv, dr_tools, pylab, numpy, math
from scipy import stats
from collections import defaultdict

def log(v):
	if o.cellfraction:
		return 1 if v>0 else 0
	return math.log(max(0.25, v), 2)

def up(V):
	return [v+10 for v in V]

def boxplot(data, **kwargs):
	for d,x in zip(data, kwargs['positions']):
		if numpy.percentile(d, 0.25) == -2:
			pylab.plot([x], [max(d)], '_')
	pylab.boxplot(data, **kwargs)

def violins(data, positions, colour):
	ax = pylab.gca()
	pos = positions
	 # from http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
	from matplotlib.patches import Rectangle
	from scipy.stats import gaussian_kde
	from numpy.random import normal
	from numpy import arange
	
	dist = max(pos)-min(pos)
	w = min(0.15*max(dist,1.0),0.5)
	for i, (d,p) in enumerate(zip(data,pos)):
		k = gaussian_kde(d) #calculates the kernel density
		m = k.dataset.min() #lower bound of violin
		M = k.dataset.max() #upper bound of violin
		x = arange(m,M,(M-m)/500.) # support for violin
		v = k.evaluate(x) #violin profile (density curve)
		v = v/v.max()*w #scaling the violin to the available space
		
		ax.fill_betweenx(x,p,v+p,facecolor=colour,alpha=0.3)
		ax.fill_betweenx(x,p,-v+p,facecolor=colour,alpha=0.3)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('matrix_accense')
	parser.add_argument('--cluster', required=True, nargs='+')
	parser.add_argument('--cmpcluster', nargs='+')
	parser.add_argument('--samples', nargs='+')
	parser.add_argument('--cellfraction', action='store_true')
	parser.add_argument('--pdf')
	parser.add_argument('--reqpos', nargs='+', default=[])
	parser.add_argument('--markers_plotted', type=int)
	o = parser.parse_args()
	
	cutoff = 0
	
	levels_cluster = defaultdict(list)
	levels_elsewhere = defaultdict(list)
	levels_everywhere = defaultdict(list)
	with open(o.matrix_accense, 'rb') as infh:
		for row in csv.DictReader(infh):
			sample = row['file'].rsplit('/', 1)[-1].split('.csv')[0]
			if o.samples and sample not in o.samples: continue
			cluster = row['population']
			if cluster in o.cluster and all(float(row['*'+m])>cutoff for m in o.reqpos):
				levels = levels_cluster
			elif o.cmpcluster is None or cluster in o.cmpcluster:
				levels = levels_elsewhere
			else:
				continue
			for parameter, value in row.items():
				if parameter in ('file', 'population'): continue
				levels[parameter].append(log(float(value)))
				levels_everywhere[parameter].append(log(float(value)))
			
			
	if o.samples: samples = set(o.samples)
	
	# give p-values
	output_items = []
	marker_order = []
	for parameter in levels_cluster:
		#print map(len, (levels_cluster[parameter], levels_elsewhere[parameter]))
		p_ranksums = stats.ranksums(levels_cluster[parameter], levels_elsewhere[parameter])[1]
		averages = map(numpy.mean, (levels_cluster[parameter], levels_elsewhere[parameter]))
		output_items.append((parameter, p_ranksums, averages, averages[0]-averages[1]))
	for item in sorted(output_items, key=(lambda I: abs(I[-1])), reverse=True):
		if item[1] < 0.01: print(*item)
		marker_order.append(item[0])
	
	
	if o.pdf:
		if o.markers_plotted is not None:
			marker_order = marker_order[:o.markers_plotted]
	
		#print map(len, [levels_cluster[P] for P in levels_cluster])
		#print levels_cluster['*CD27']
		#dr_tools.violin_plot(pylab.axes(), [[v+10 for v in levels_cluster['*CD27']]], [0])
		#dr_tools.violin_plot(pylab.axes(), [[v+10 for v in levels_elsewhere['*CD27']]], [1])
		#dr_tools.violin_plot(pylab.axes(), [up(levels_cluster[P]) for P in levels_cluster], [i*3 for i,P in enumerate(levels_cluster)])
		#dr_tools.violin_plot(pylab.axes(), [up(levels_elsewhere[P]) for P in levels_cluster], [i*3+1 for i,P in enumerate(levels_cluster)])
		pylab.xticks([i*3+0.5 for i,P in enumerate(marker_order)], [P.lstrip('*') for i,P in enumerate(marker_order)], rotation=90)
		violins([levels_cluster[P] for P in marker_order], positions=[i*3 for i,P in enumerate(marker_order)], colour='y')
		violins([levels_elsewhere[P] for P in marker_order], positions=[i*3+1 for i,P in enumerate(marker_order)], colour='r')
		pylab.xlim(-2, len(marker_order)*3)
		pylab.gca().xaxis.set_ticks_position('none')
		pylab.savefig(o.pdf)
