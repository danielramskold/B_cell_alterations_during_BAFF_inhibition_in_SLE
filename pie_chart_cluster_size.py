#! python3
from matplotlib import pyplot
import argparse, pandas, numpy
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('csvfile')
	parser.add_argument('-o', '--pdf_out', required=True)
	parser.add_argument('-g', '--group_pops', nargs='+', action='append', default=[], type=int)
	parser.add_argument('--pop_colname', default='population')
	parser.add_argument('--explode', action='store_true')
	parser.add_argument('--order', type=int, nargs='+')
	parser.add_argument('--req_in_file_name')
	o = parser.parse_args()
	
	df = pandas.read_csv(o.csvfile)
	
	if o.req_in_file_name is not None:
		df = df[df['file'].str.contains(o.req_in_file_name)]
	
	if not o.explode:
		for group in o.group_pops:
			name = ','.join([str(v) for v in group])
			for cluster in group:
				df.loc[df[o.pop_colname] == cluster, o.pop_colname] = name
	
	counts = df[o.pop_colname].value_counts()
	if o.order is not None:
		counts = counts.reindex(o.order)
	fractions = [c/sum(counts) for c in counts]
	names = [str(n) for n in counts.index]
	print(names, fractions)
	
	
	wedges, texts = pyplot.pie(fractions, labels=names)
	if o.explode:
		# from https://stackoverflow.com/questions/20549016/explode-multiple-slices-of-pie-together-in-matplotlib
		# untested
		groups = [list(sorted(g)) for g in o.group_pops]
		radfraction = 0.1
		patches = []
		for i in groups:
		  ang = numpy.deg2rad((wedges[i[-1]].theta2 + wedges[i[0]].theta1)/2,)
		  for j in i:
		    we = wedges[j]
		    center = (radfraction*we.r*np.cos(ang), radfraction*we.r*np.sin(ang))
		    patches.append(mpatches.Wedge(center, we.r, we.theta1, we.theta2))
		colors = numpy.linspace(0, 1, len(patches))
		collection = PatchCollection(patches, cmap=pyplot.cm.hsv)
		collection.set_array(numpy.array(colors))
		pyplot.add_collection(collection)
	
	
	pyplot.axes().set_aspect('equal', 'datalim')
	pyplot.savefig(o.pdf_out)
