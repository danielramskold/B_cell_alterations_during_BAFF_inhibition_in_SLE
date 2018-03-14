from __future__ import division
import argparse, pylab, dr_tools, numpy, random
from collections import defaultdict
from scipy import stats

class Change:
	def __init__(self):
		self.ends = []
		self.startpoint = None
		self.endpoint = None
		self.name = None
		self.pvalue = 2
		self.qvalue = 2
	
	def point(self):
		return (self.endpoint - self.startpoint)
	
	def err(self):
		return abs(self.ends[1]-self.ends[0])/2
	
	def __repr__(self):
		return 'Change<%s-%s-%s>'%(self.name, self.startpoint, self.endpoint)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile1', help='output from expression_time_dependence with setting --correlation paired_confint')
	parser.add_argument('infile2', help='output from expression_time_dependence with setting --correlation paired_confint')
	parser.add_argument('outplot')
	parser.add_argument('--lim', type=float, nargs=2)
	parser.add_argument('--input_nodivision', action='store_true', help='if --nodivision_confint was used for the input file generation')
	parser.add_argument('--shorten_names', action='store_true')
	parser.add_argument('--random_seed', type=float)
	parser.add_argument('-c', '--marker_colour', action='append', nargs=2, default=[])
	o = parser.parse_args()
	
	if o.random_seed is not None: random.seed(o.random_seed)
	
	changesD = defaultdict(lambda: [Change(), Change()])
	for ii, infile in enumerate((o.infile1, o.infile2)):
		for p in dr_tools.splitlines(infile):
			name = p[2]
			c = changesD[name][ii]
			c.name = name.upper() if 'Fox' in name else name
			if o.shorten_names:
				c.name = c.name.split('_')[-1].rstrip('+')
			bottom_ch, top_ch = map(float, p[4].strip('()').split(', '))
			change = numpy.mean([top_ch, bottom_ch]) # relative value, i.e. start*(1+change) = end
			avg_abundance = float(p[5]) # = (start+end)/2
			if o.input_nodivision:
				c.startpoint = avg_abundance + change / 2
				c.endpoint = c.startpoint + change
				c.ends.append(top_ch + c.startpoint)
				c.ends.append(bottom_ch + c.startpoint)
			else:
				c.startpoint = 2 * avg_abundance / (2 + change)
				c.endpoint = c.startpoint * (1 + change)
				c.ends.append((top_ch + 1) * c.startpoint)
				c.ends.append((bottom_ch + 1) * c.startpoint)
			c.pvalue = float(p[0])
			c.qvalue = float(p[1])
	
	xarr = []
	yarr = []
	colours = defaultdict(dr_tools.randomcolour)
	colours.update(dict(o.marker_colour))
	for c1, c2 in changesD.values():
		if c1.name is None or c2.name is None: continue
		x = c1.point()
		xerr = c1.err()
		y = c2.point()
		yerr = c2.err()
		ebar_intensity = 0.001/min(xerr,yerr)
		#ebar_intensity = 0.004/(xerr+yerr)
		rcol = colours[c1.name]
		pylab.errorbar([x], [y], [yerr], [xerr], ecolor=dr_tools.mixcolours(['#000000', '#ffffff'], [ebar_intensity, 1-ebar_intensity]), fmt='', capsize=0, zorder=0, alpha=0.8)
		pylab.plot([x], [y], '.', zorder=2000, color=dr_tools.mixcolours(['#000000', rcol], [0.5, 0.5]), alpha=0.9)
		xarr.append(x)
		yarr.append(y)
		pylab.text(x, y, c1.name, fontsize=10, color=rcol, fontname='Arial', zorder=100, alpha=0.6)
	
	pylab.ylabel(o.infile2)
	pylab.xlabel(o.infile1)
	
	r, p = stats.spearmanr(xarr, yarr)
	pylab.title('Spearman rho = %f, p = %e'%(r, p))
	
	if o.lim is not None:
		pylab.xlim(*o.lim)
		pylab.ylim(*o.lim)
		
	
	# draw diagonal line
	lowerlim = max(pylab.gca().get_xlim()[0], pylab.gca().get_ylim()[0])
	upperlim = min(pylab.gca().get_xlim()[1], pylab.gca().get_ylim()[1])
	pylab.plot([lowerlim, upperlim], [lowerlim, upperlim], 'k:')
	
	print ' '.join('-c "%s" "%s"'%I for I in colours.items())
	
	pylab.savefig(o.outplot)
