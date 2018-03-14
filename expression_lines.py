from __future__ import division
import argparse, csv, dr_tools, pylab, numpy, math, random, os
from scipy import stats, interpolate
from collections import defaultdict

def log(v): return math.log(max(0.25, float(v)), 2)

def up_pos(V): return [v+3 for v in V]
def up_neg(V): return [v+3 for v in V]

def existance(V): return [int(v>0) for v in V]
def nonexistance(V): return [int(not v) for v in existance(V)]

def product(V):
	ret = 1
	for v in V: ret *= v
	return ret

def splitlines(infile):
	with open(infile, 'rU') as infh:
		for line in infh:
			yield line.rstrip('\r\n').split('\t')

def slidingwindow(seq, n):
    # http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    from itertools import islice
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def interpolation(xarr, yarr):
	if o.interpolate_method == '1m_lin':
		return Nearest_linear_mix(xarr, yarr, o.extrapolate_area)
	if o.interpolate_method == 'lin_e1m':
		return Linear_extrapnearest(xarr, yarr, o.extrapolate_area)
	if o.interpolate_method == 'lin_erel':
		return Linear_extrapnearest_rel(xarr, yarr, o.extrapolate_area)
	try:
		return interpolate.interp1d(xarr, yarr, kind=o.interpolate_method)
	except ValueError:
		return interpolate.interp1d(xarr, yarr, kind='linear', bounds_error=False)


class Nearest_linear_mix:
	# does interpolation/extrapolation, using nearest if near enough (1 month, if o.extrapolate_area=1 as is default), otherwise linear if inbetween points, or gives NaN when more than 1 month off the min or max value
	def __init__(self, xarr, yarr, maxnear):
		self.arr = zip(xarr, yarr)
		self.linear_f = interpolate.interp1d(xarr, yarr, kind='linear')
		self.maxnear = maxnear
		self.fillvalue = None
	
	def __call__(self, X):
		Y = []
		for xp in X:
			near_x = None
			lastdist = self.maxnear
			for x,y in self.arr:
				if abs(x - xp) <= lastdist:
					lastdist = abs(x - xp)
					near_x = (x,y)
			
			if near_x is None:
				try: y = float(self.linear_f(xp))
				except ValueError: y = self.fillvalue
			else:
				y = near_x[1]
			Y.append(y)
		return Y

class Linear_extrapnearest:
	# does interpolation/extrapolation, using linear inbetween points, and nearest for up to 1 month outside the points
	def __init__(self, xarr, yarr, maxnear):
		self.arr = zip(xarr, yarr)
		self.linear_f = interpolate.interp1d(xarr, yarr, kind='linear')
		self.maxnear = maxnear
		self.fillvalue = None
	
	def __call__(self, X):
		Y = []
		for xp in X:
			try:
				y = float(self.linear_f(xp))
			except ValueError:
				near_x = None, self.fillvalue
				lastdist = self.maxnear
				for x,y in self.arr:
					if abs(x - xp) <= lastdist:
						lastdist = abs(x - xp)
						near_x = (x,y)
				y = near_x[1]
			Y.append(y)
		return Y

class Linear_extrapnearest_rel:
	# does interpolation/extrapolation, using linear inbetween points, and nearest for up to distance outside the points 
	def __init__(self, xarr, yarr, maxnear):
		self.arr = zip(xarr, yarr)
		self.linear_f = interpolate.interp1d(xarr, yarr, kind='linear')
		self.maxnear_rel = maxnear
		self.fillvalue = None
	
	def __call__(self, X):
		Y = []
		for xp in X:
			try:
				y = float(self.linear_f(xp))
			except ValueError:
				near_x = None, self.fillvalue
				lastdist = None
				for x,y in self.arr:
					if abs(x - xp) <= abs(x * self.maxnear_rel) and (lastdist is None or abs(x - xp) <= lastdist):
						near_x = (x,y)
						lastdist = abs(x - xp)
				y = near_x[1]
			Y.append(y)
		return Y


category_to_number = {'Yes':1, 'No':0, 'Male':0, 'Female':1}

def midY(levels_s):
	if o.invertY:
		if o.calc != 'coexist': raise Exception
		return 1-numpy.mean(levels_s)
	else: return numpy.mean(levels_s)

class Bin:
	def __init__(self, minx, maxx):
		self.minx = minx
		self.maxx = maxx
		self.midx = (self.minx+self.maxx)/2
		self.yarr = []
		self.perID = dict()
	
	def add1within(self, xarr, yarr, ID=None):
		within = []
		for x, y in zip(xarr, yarr):
			if self.minx <= x <= self.maxx:
				within.append((x, y))
		if len(within) == 1:
			y = within[0][1]
			if y is not None:
				if o.ylim is not None and y > o.ylim[1]*1e5: y = o.ylim[1]*1.01 # to fix an Illustrator problem
				self.yarr.append(y)
				self.perID[ID] = y
		elif len(within) > 1:
			raise Exception # multiple data points from same patient in one bin
			# add the closest value to mid
			#self.yarr.append(sorted(within, key=lambda w: abs(w[0]-self.midx))[0][1])
	
	def midy(self):
		return numpy.median(self.yarr)
	
	def uppery(self):
		return numpy.percentile(self.yarr, 75)
	
	def lowery(self):
		return numpy.percentile(self.yarr, 25)
	
	def Pvalue_text(self, compare_bin, test='signedrank'):
		if test == 'ranksums':
			Pvalue = stats.ranksums(self.yarr, compare_bin.yarr)[1]
		elif test == 'signedrank':
			commonIDs = list(set(self.perID) & set(compare_bin.perID))
			Pvalue = stats.wilcoxon([self.perID[ID] for ID in commonIDs], [compare_bin.perID[ID] for ID in commonIDs])[1]
			print [self.perID[ID] for ID in commonIDs], [compare_bin.perID[ID] for ID in commonIDs], Pvalue
		if Pvalue < 0.001: return '***'
		elif Pvalue < 0.01: return '**'
		elif Pvalue < 0.05: return '*'
		elif Pvalue >= 0.05: return 'ns'
		#elif Pvalue >= 0.05: return 'P=%.2f'%Pvalue
		else:
			# happens when the two lists are identical
			return 'ns'
	
	def saypatients(self):
		if o.quite: return
		print ''
		print self.midx, ':'
		for ID in sorted(self.perID):
			print ID

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--csv_in', nargs='+', required=True)
	parser.add_argument('-p', '--proteins', nargs='*', required=True)
	parser.add_argument('-n', '--proteins_neg' ,nargs='*', default=[])
	parser.add_argument('-o', '--pdf')
	parser.add_argument('--calc', choices=['coexist', 'sumlog', 'prodlog'], default='coexist')
	parser.add_argument('-t', '--patienttable')
	parser.add_argument('-x', '--x_column', default='Time_from_baseline_months')
	parser.add_argument('-q', '--quite', action='store_true')
	parser.add_argument('-c', '--cellnumbers')
	parser.add_argument('--invertY', action='store_true')
	parser.add_argument('--patientnorm', action='store_true')
	parser.add_argument('--patientnorm_div', action='store_true')
	parser.add_argument('--noplot', action='store_true')
	parser.add_argument('--divxbaseline', action='store_true')
	parser.add_argument('--divybaseline', action='store_true')
	parser.add_argument('--colour_patients', nargs='+')
	parser.add_argument('--xbin', nargs=2, type=float, action='append', default=[])
	parser.add_argument('--title')
	parser.add_argument('--ylim', nargs=2, type=float)
	parser.add_argument('--use_from_list', choices=['0', '1'], help='which to select if -i specific a list of patients followed by 0 or 1 for each')
	parser.add_argument('--show_splines', action='store_true', help='not tested with --xbin')
	parser.add_argument('--interpolate_method', choices=['lin_erel', 'lin_e1m', '1m_lin', 'linear', 'cubic', 'nearest', 'zero', 'quadratic'], default='cubic')
	parser.add_argument('--interpolate_points', nargs='+', type=float)
	parser.add_argument('--extrapolate_area', type=float, help='use with --interpolate_method lin_e1m or 1m_lin or lin_erel')
	parser.add_argument('--writepatient', action='store_true')
	o = parser.parse_args()
	
	if o.extrapolate_area is None:
		if o.interpolate_method in ('lin_e1m', '1m_lin'):
			o.extrapolate_area = 1   # means 1 month
		elif o.interpolate_method in ('lin_erel',):
			o.extrapolate_area = 0.07    # means 7% of the x value
	
	if o.cellnumbers is not None:
		cellnumbers = dict()
		for sample, nstr in splitlines(o.cellnumbers):
			if nstr:
				cellnumbers[sample] = float(nstr)
	
	proteins_pos = [p.strip('*') for p in o.proteins]
	proteins_neg = [p.strip('*') for p in o.proteins_neg]
	random.seed(0)
	if o.pdf is None:
		o.pdf = 'expression_' + ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg]) + '.pdf'
	elif o.title is None:
		o.title = ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
	
	samples_csv = o.csv_in
	if o.use_from_list is not None:
		samples_to_use = set()
		samples_possible = [f for f in os.listdir('.') if f.endswith('.csv')]
		for p in dr_tools.splitlines(o.csv_in[0]):
			if p[1] == o.use_from_list:
				samples_to_use.update([s for s in samples_possible if s == p[0] or s.startswith(p[0]+'_')])
		samples_csv = list(samples_to_use)
		o.csv_in = o.csv_in[1:]
	if len(o.csv_in) == 1 and o.csv_in[0].endswith('.txt'):
		with open(o.csv_in[0], 'rU') as infh:
			samples_csv = [f.strip() for f in infh.read().split()]
	
	levels = defaultdict(list)
	for csv_in in samples_csv:
		with open(csv_in, 'rb') as infh:
			sample = csv_in.rsplit('/', 1)[-1].split('.csv')[0]
			if sample.count('_') == 2: sample = sample.rsplit('_',1)[0]
			for row in csv.DictReader(infh):
				values_pos = [log(row[p]) for p in proteins_pos]
				values_neg = [log(row[p]) for p in proteins_neg]
				if o.calc == 'sumlog':
					levels[sample].append(sum(values_pos)-sum(values_neg))
				elif o.calc == 'prodlog':
					levels[sample].append(product(up_pos(values_pos))/product(up_neg(values_neg)))
				elif o.calc == 'coexist':
					levels[sample].append(product(existance(values_pos) + nonexistance(values_neg)))
				else: raise Exception
	
	time_since_baseline_per_sample = defaultdict(dict)
	if o.patienttable:
		pylab.xlabel(o.x_column)
		for li, p in enumerate(splitlines(o.patienttable)):
			if li == 0:
				time_i = p.index(o.x_column)
				patient_i = p.index('CMM_ID')
				#patient_i = p.index('Patient_ID') # doesn't have the same number of timepoints per patient as the cytof samples, must be wrong
				visit_i = p.index('Visit_ID')
			else:
				visit = int(p[visit_i])-1
				patient = p[patient_i]
				if p[time_i] in category_to_number:
					time = category_to_number[p[time_i]]
				else:
					try: time = float(p[time_i])
					except ValueError:
						if not o.quite and patient:
							print 'Missing data:', patient, visit
						continue
				if visit not in time_since_baseline_per_sample[patient]:
					time_since_baseline_per_sample[patient][visit] = time
	
	bins = [Bin(*boundaries) for boundaries in o.xbin]
	locations_x = []
	locations_y = []
	per_patient_locations = []
	#direction_continuation = []
	samples = sorted(levels.keys())
	patient_samples = defaultdict(list)
	for sample in samples:
		patient_samples[sample.split('_time')[0]].append(sample)
	for patient in sorted(patient_samples.keys()):
		if patient == 'UnknownDonor': continue
		samples = patient_samples[patient]
		
		if o.cellnumbers is None:
			yarr = [midY(levels[sample]) for sample in samples]
		else:
			samples = [sample for sample in samples if sample in cellnumbers]
			if len(samples) == 0: continue
			yarr = [midY(levels[sample])*cellnumbers[sample] for sample in samples]
		
		if o.divybaseline:
			mid_y = max(1e-16, yarr[0])
			yarr = [y / mid_y for y in yarr]
		if o.patientnorm_div:
			mid_y = numpy.mean(yarr)
			yarr = [y / mid_y for y in yarr]
		if o.patientnorm:
			mid_y = numpy.mean(yarr)
			yarr = [y - mid_y for y in yarr]
		
		
		#xarr = range(len(samples))
		xarr = [int(sample.split('_time')[1].split('_')[0].split('.')[0]) for sample in samples]
		if o.patienttable:
			try:
				earliest_timepoint = min(time_since_baseline_per_sample[patient].keys())
				
				if o.divxbaseline:
					xarr = [time_since_baseline_per_sample[patient][visit]/time_since_baseline_per_sample[patient][earliest_timepoint] for visit in xarr]
				else:
					xarr = [time_since_baseline_per_sample[patient][visit] for visit in xarr]
			except KeyError, ValueError:
				if not o.quite: print patient
				continue
		
		if o.interpolate_points is not None and not o.show_splines:
			f_splines = interpolation(xarr, yarr)
			xarr, yarr = o.interpolate_points, f_splines(o.interpolate_points)
		
		
		#for a, b, c in slidingwindow(yarr, 3):
		#	if (b-a)*(c-b) > 0: direction_continuation.append(1)
		#	elif (b-a)*(c-b) < 0: direction_continuation.append(0)
		if o.colour_patients is None:
			c = dr_tools.randomcolour()
		elif o.colour_patients[0].endswith('.txt'):
			for p in dr_tools.splitlines(o.colour_patients[0]):
				if p[0] == patient:
					c = {'1':'#3465a4', '0':'#73d216', '2':'#d9b98e'}[p[1]]
					break
			else:
				c = '#ff8888'
			#c = dr_tools.mixcolours([c, dr_tools.randomcolour()], [0.8, 0.2])
		else:
			c = '#5555ff' if patient in o.colour_patients else '#ccffcc'
			c = dr_tools.mixcolours([c, dr_tools.randomcolour()], [0.67, 0.33])
		if bins:
			for b in bins:
				b.add1within(xarr, yarr, patient)
		else:
			if o.show_splines:
				if o.interpolate_points:
					Xsplines = o.interpolate_points
				else:
					Xsplines = pylab.arange(min(xarr), max(xarr), 0.2)
				if len(xarr) >= 4:
					f_splines = interpolation(xarr, yarr)
					pylab.plot(Xsplines, f_splines(Xsplines), '--', color=c)
					pylab.plot(xarr, yarr, 'o', label=patient, color=c)
			else:
				pylab.plot(xarr, yarr, 'o-', label=patient, color=c)
			if o.writepatient:
				pylab.text(xarr[-1], yarr[-1], patient, fontsize=6, color=c) # places patient id to the right of the line
		locations_x.extend(xarr)
		locations_y.extend(yarr)
		per_patient_locations.append(zip(xarr, yarr))
	#pylab.legend()
	if bins:
		
		pylab.errorbar(x=[b.midx for b in bins], y=[b.midy() for b in bins], yerr=[[b.midy()-b.lowery() for b in bins], [b.uppery()-b.midy() for b in bins]], fmt='ko-', ecolor='k')
		for bi, b in enumerate(bins):
			if bi == 0:
				text = ''
			else:
				text = b.Pvalue_text(bins[0])
			text += ' n='+str(len(b.yarr))
			if text:
				y_txt = b.uppery()
				if o.ylim is not None: y_txt = min(y_txt, o.ylim[1])
				pylab.text(x=b.midx, y=y_txt, s=text)
		
		for b in bins: b.saypatients()
		
		bin_diff_pvalues = []
		for bin1, bin2 in zip(bins[:-1], bins[1:]):
			bin_diff_pvalues.append(stats.ranksums(bin1.yarr, bin2.yarr)[1])
		if not o.quite: print bin_diff_pvalues, 'P_adjacent_bins'
	
	#if o.calc == 'coexist': pylab.ylabel('Fraction of T cells')
	#pylab.xticks([0, 1], ['0m', '6m'])
	
	
	if o.ylim is not None: pylab.ylim(*o.ylim)
	if o.title is not None: pylab.title(o.title)
	if not o.noplot: pylab.savefig(o.pdf)
	
	if not(o.interpolate_points is not None and not o.show_splines):
		print stats.spearmanr(locations_x, locations_y), numpy.mean(locations_y)
	if not o.quite:
		if o.patienttable is None:
			print stats.spearmanr([abs(2.5-x) for x in locations_x], locations_y)
			print numpy.median([y for x,y in zip(locations_x, locations_y) if x==0]), 'median y at x=0'
			print numpy.median([y for x,y in zip(locations_x, locations_y) if x>=6]), 'median y at x>=6'
			print numpy.mean([y for x,y in zip(locations_x, locations_y) if 11<x<13]), 'avg y at 11<x<13'
			print numpy.mean([y for x,y in zip(locations_x, locations_y) if x>=23]), 'avg y at x>=23'
		#print stats.binom_test([direction_continuation.count(1), direction_continuation.count(0)]), direction_continuation.count(1), direction_continuation.count(0)  # flawed test, not stat independence between data points
		
		# a patient-wise test, to see it's not driven by 1 patient or so
		def flatten_and_correlation(locationsZ):
			xloc = [x for arr in locationsZ for x,y in arr]
			yloc = [y for arr in locationsZ for x,y in arr]
			return stats.spearmanr(xloc, yloc)[0]
		if not o.quite: print 'Patient-wise bootstrap P-value(rho<0), min_rho, max_rho: ', dr_tools.bootstrap(flatten_and_correlation, (per_patient_locations,), controls=1000, processes=4)
		
