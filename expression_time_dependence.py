from __future__ import division
import argparse, csv, dr_tools, pylab, numpy, math, itertools, os, sys, gzip
from scipy import stats
from collections import defaultdict, namedtuple
from dr_tools import MINE_MIC

def log(v): return math.log(max(0.25, float(v)), 2)

def up_pos(V): return [v+3 for v in V]
def up_neg(V): return [v+3 for v in V]

def existance(V, genes=None):
	if genes:
		return [int(v>cutoffsD.get(gene, cutoff)) for v, gene in zip(V, genes)]
	else:
		return [int(v>cutoff) for v in V]

def nonexistance(V, genes=None): return [int(not v) for v in existance(V, genes)]

def product(V):
	ret = 1
	for v in V: ret *= v
	return ret

def combinations_maxN(pick_from, max_n, exclusion=set()):
	restricted_pick_from = set(pick_from) - set(exclusion)
	for n in range(max_n+1):
		for chosen in itertools.combinations(restricted_pick_from, n):
			yield chosen

class PVal:
	def __init__(self, p, text):
		self.p = p
		self.q = None
		self.text = text
	
	def __repr__(self):
		if self.p >= 0.0002:
			return '%.4f\t%.4f\t%s'%(self.p, self.q, self.text)
		elif  self.p >= 0.000002:
			return '%.6f\t%.6f\t%s'%(self.p, self.q, self.text)
		else:
			return '%.1e\t%.1e\t%s'%(self.p, self.q, self.text)

NaN = float('nan')
distfunctions = {'across':float, 'mid':(lambda t: 2.5-abs(2.5-t)), 'start':(lambda t: min(2,t)), 'end':(lambda t: 2-min(2,5-t)), 'late':(lambda t: 3-min(3,5-t)), 'after1year':(lambda months: max(0, months - 12)), 'before6months': (lambda months: min(6, months)), '2groups':(lambda x:x), 'before6months_cut': (lambda months: NaN if months > 6 else months), 'after1year_cut':(lambda months: NaN if months < 12 else months), 'after6months_cut': (lambda months: NaN if months <= 6 else months), 'after6months':(lambda months: max(0, months - 6))}

def calc_level(calc, row, pos, neg):
	if calc in ('coexist', 'sum'):
		try:
			values_pos = [float(row[p]) for p in pos]
			values_neg = [float(row[p]) for p in neg]
		except:
			print >>sys.stderr, pos, neg, row
			raise
		if calc == 'coexist':
			return product(existance(values_pos, pos) + nonexistance(values_neg, neg))
		elif calc == 'sum':
			return sum(values_pos) - sum(values_neg)
	
	values_pos = [log(row[p]) for p in pos]
	values_neg = [log(row[p]) for p in neg]
	if calc == 'sumlog':
		return sum(values_pos)-sum(values_neg)
	elif calc == 'prodlog':
		return product(up_pos(values_pos))/product(up_neg(values_neg))
	else:
		raise Exception

def split_pos_neg(protline):
	proteins_pos = []
	proteins_neg = []
	pos_split = protline.split('+')
	for pos in pos_split[:-1]:
		proteins_pos.append(pos)
	neg_split = pos_split[-1].split('-')
	ni = 0
	while ni < len(neg_split)-1:
		if neg_split[ni] in ('HLA', 'Ki'):
			proteins_neg.append(neg_split[ni]+'-'+neg_split[ni+1])
			ni+=2
		else:
			proteins_neg.append(neg_split[ni])
			ni+=1
	return proteins_pos, proteins_neg

def splitlines(infile):
	with openfile(infile, 'rU') as infh:
		for line in infh:
			yield line.rstrip('\r\n').split('\t')
category_to_number = {'Yes':1, 'No':0, 'Male':0, 'Female':1}

def filter_listing_loop(row_keys, o, filter_req_proteins):
	if o.set_filter is not None:
		for protline in o.set_filter:
			filter_pos, filter_neg = split_pos_neg(protline)
			yield set(filter_pos)|set(filter_neg), tuple(filter_pos), tuple(filter_neg)
	else:
		for filter_pos in combinations_maxN(row_keys, o.filter_pos):
			exclude_1 = set(filter_pos)
			for filter_neg in combinations_maxN(row_keys, o.filter_neg, exclude_1):
				exclude_2 = exclude_1|set(filter_neg)
				if not exclude_2 and not skip_filter: continue
				if o.filter_req_proteins and (set(filter_pos)|set(filter_neg))&filter_req_proteins != filter_req_proteins:
					continue
				yield exclude_2, filter_pos, filter_neg

def paired_ttest(val1, val0):
	if all(v==0 for v in val1) and all(v==0 for v in val0):
		return 0, float('nan')
	try: stat, P = stats.ttest_rel(val1, val0)
	except:
		print >>sys.stderr, val1, val0
		raise
	stat = numpy.mean(val1)-numpy.mean(val0)
	return stat, P

def paired_confint(val1, val0, confint=0.95):
	if not o.nodivision_confint: diff = [(v1-v0)/v0 for v1,v0 in zip(val1, val0)] # might need a pseudocount to not crash
	#diff = [0 if v1==v0 else (v1-v0)/numpy.mean([v1,v0]) for v1,v0 in zip(val1, val0)] # will crash if --patientnorm is set, and will work poorly with negative values
	if o.nodivision_confint: diff = [v1-v0 for v1,v0 in zip(val1, val0)]
	stat, P = stats.ttest_1samp(diff, 0)
	halfwidth = numpy.std(diff) / math.sqrt(len(diff)) * stats.t.ppf(1-(0.5-abs(0.5-confint))/2, len(diff)-1) # changed from v8
	mid = numpy.mean(diff)
	return (mid-halfwidth, mid+halfwidth), P

def unpaired_confint(val1, val0, confint=0.95):
	# http://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals5.html
	n1minus1 = len(val1)-1
	n0minus1 = len(val0)-1
	s1 = numpy.std(val1)
	s0 = numpy.std(val0)
	Sp = math.sqrt((n0minus1*s0**2+ n1minus1*s1**2)/(n1minus1+n0minus1))
	halfwidth = stats.t.ppf(1-(0.5-abs(0.5-confint))/2, n1minus1+n0minus1) * Sp * math.sqrt(1/len(val0) + 1/len(val1))
	mid = numpy.mean(val1)-numpy.mean(val0)
	mid0 = numpy.mean(val0)
	stat, P = stats.ttest_ind(val1, val0)
	if o.nodivision_confint:
		return ((mid-halfwidth), (mid+halfwidth)), P
	else:
		return ((mid-halfwidth)/mid0, (mid+halfwidth)/mid0), P

def perm_mean(val1, val0):
	P = dr_tools.permutationtest(numpy.mean, val1, val0, controls=1000)
	stat = numpy.mean(val1)-numpy.mean(val0)
	return P, stat

def openfile(filepath, *moreargs, **kwargs):
	if filepath.endswith('.gz'): return gzip.open(filepath, *moreargs, **kwargs)
	else: return open(filepath, *moreargs, **kwargs)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--csv_in', nargs='+', required=True)
	parser.add_argument('-p', '--num_pos', type=int, default=1)
	parser.add_argument('-n', '--num_neg', type=int, default=0)
	parser.add_argument('--filter_pos', type=int, default=0)
	parser.add_argument('--filter_neg', type=int, default=0)
	parser.add_argument('-P', '--proteins', nargs='+')
	parser.add_argument('-R', '--req_proteins', nargs='+', default=[])
	parser.add_argument('-F', '--filter_req_proteins', nargs='+', default=[])
	parser.add_argument('--calc', choices=['coexist', 'sumlog', 'prodlog', 'sum'], default='coexist')
	parser.add_argument('--acrosscells', choices=['mean', 'sum'], default='mean')
	parser.add_argument('--peakpoint', choices=distfunctions.keys(), nargs='+')
	parser.add_argument('--cutoff', type=float, default=0)  # cutoff 1 corresponds to _v3
	parser.add_argument('-q', '--maxq', type=float, default=1)
	parser.add_argument('--correlation', choices=['pearson', 'spearman', 'MIC', 'ttest', 'ranksums', 'paired_ttest', 'perm_mean', 'signedrank', 'paired_confint', 'paired_confint996', 'paired_confint998', 'unpaired_confint', 'unpaired_confint998'], default='spearman')
	parser.add_argument('--patientnorm', action='store_true')
	parser.add_argument('--maxlines', type=int)
	parser.add_argument('-c', '--given_combo', nargs='+', help='e.g. reqCD3e+_CD27+ or IgD+IgA-')
	parser.add_argument('-t', '--patienttable')
	parser.add_argument('--simpletable', action='store_true')
	parser.add_argument('-x', '--x_column', default='Time_from_baseline_months')
	parser.add_argument('--divxbaseline', action='store_true')
	parser.add_argument('--cellnumbers')
	parser.add_argument('--use_from_list', choices=['0', '1'], help='which to select if -i specific a list of patients followed by 0 or 1 for each')
	parser.add_argument('--set_filter', nargs='+')
	parser.add_argument('--saycombinedP', action='store_true')
	parser.add_argument('--skiplineprobability', type=float, default=0)
	parser.add_argument('-m', '--marker_cutoff', default=[], action='append', nargs=2, help='argument are pair of marker and value, several instances can be used e.g. -m CD127 2 -m CD3e 0')
	parser.add_argument('--output_summary_values', action='store_true')
	parser.add_argument('--nodivision_confint', action='store_true')
	o = parser.parse_args()
	
	# load cell numbers
	if o.cellnumbers is not None:
		cellnumbers = dict()
		for sample, nstr in splitlines(o.cellnumbers):
			if nstr:
				cellnumbers[sample] = float(nstr)
	
	correlation, pairwise_test = {'spearman':(stats.spearmanr, False), 'pearson':(stats.pearsonr, False), 'MIC':(MINE_MIC, False), 'ttest':(stats.ttest_ind, True), 'ranksums':(stats.ranksums, True), 'paired_ttest':(paired_ttest, True), 'perm_mean':(perm_mean, True), 'signedrank':(stats.wilcoxon, True), 'paired_confint': (paired_confint, True), 'paired_confint998': ((lambda a,b: paired_confint(a,b, 0.998)), True), 'paired_confint996': ((lambda a,b: paired_confint(a,b, 0.996)), True), 'unpaired_confint': (unpaired_confint, True), 'unpaired_confint998': ((lambda a,b: unpaired_confint(a,b, 0.998)), True)}[o.correlation]
	
	if o.peakpoint is None: o.peakpoint = ['across'] if not pairwise_test else ['2groups']
	
	cutoff = o.cutoff
	cutoffsD = dict((k, float(v)) for k, v in o.marker_cutoff)
	
	if o.skiplineprobability:
		import random
		random.seed(0)
	
	skip_filter = not (o.filter_pos or o.filter_neg or o.set_filter)
	
	# look for pre-defined marker combinations
	levels_per_combo = dict()
	combos = []
	if o.given_combo is not None:
		for given_combo in o.given_combo:
			if given_combo.startswith('req'):
				reqline = given_combo.split('_')[0][3:]
				filter_pos, filter_neg = split_pos_neg(reqline)
				protline = given_combo.split('_')[1]
			else:
				filter_pos, filter_neg = tuple(), tuple()
				protline = given_combo
			proteins_pos, proteins_neg =  split_pos_neg(protline)
			combos.append(tuple(map(tuple, (filter_pos, filter_neg, proteins_pos, proteins_neg))))
			levels_per_combo[combos[-1]] = defaultdict(list)
	if o.num_pos + o.num_neg < 1 and not combos: raise Exception
	req_proteins = set(o.req_proteins)
	filter_req_proteins = set(o.filter_req_proteins)
	
	# select samples
	samples_csv = o.csv_in
	if o.use_from_list is not None:
		samples_to_use = set()
		samples_possible = [f for f in os.listdir('.') if f.endswith('.csv') or f.endswith('.csv.gz')]
		for p in dr_tools.splitlines(o.csv_in[0]):
			if p[1] == o.use_from_list:
				samples_to_use.update([s for s in samples_possible if s == p[0] or s.startswith(p[0]+'_')])
		samples_csv = list(samples_to_use)
		o.csv_in = o.csv_in[1:]
	if len(o.csv_in) == 1 and o.csv_in[0].endswith('.txt'):
		with openfile(o.csv_in[0], 'rU') as infh:
			if o.use_from_list is not None:
				samples_csv = [f.strip() for f in infh.read().split() if f.strip() in samples_to_use]
			else:
				samples_csv = [f.strip() for f in infh.read().split()]
	
	CD3double = False
	if not combos:
		if o.proteins:
			row_keys = set(o.proteins)
		else:
			# collect list of marker names
			row_keys = set()
			for csv_in in samples_csv:
				with openfile(csv_in, 'rb') as infh:
					for ri, row in enumerate(csv.DictReader(infh)):
						row_keys.update(set(row))
						break # added 2 Jan 2016, should make no difference
			row_keys -= set(['DNAIr', 'DNAIr.1', 'EQBeads', 'EQBeads.1', 'EQBeads.2', 'Cisplatin'])
			
			
			if 'CD3' in row_keys and 'CD3e' in row_keys:
				CD3double = True
				row_keys.remove('CD3e')
			
			# make marker combinations
			
			filter_generator = filter_listing_loop(row_keys, o, filter_req_proteins)
			
			for exclude_2, filter_pos, filter_neg in filter_generator:
				for proteins_pos in combinations_maxN(row_keys, o.num_pos, exclude_2):
					exclude_3 = exclude_2|set(proteins_pos)
					for proteins_neg in combinations_maxN(row_keys, o.num_neg, exclude_3):
						if o.req_proteins and (set(proteins_pos)|set(proteins_neg))&req_proteins != req_proteins:
							continue
						combos.append((filter_pos, filter_neg, proteins_pos, proteins_neg))
						levels_per_combo[combos[-1]] = defaultdict(list)
				
	
	
	# assign values to marker combinations, using the data
	for csv_in in samples_csv:
		with openfile(csv_in, 'rb') as infh:
			sample = csv_in.rsplit('/', 1)[-1].split('.csv')[0]
			for ri, row in enumerate(csv.DictReader(infh)):
				if o.maxlines and ri >= o.maxlines:break
				if o.skiplineprobability and o.skiplineprobability > random.random(): continue
				if CD3double and 'CD3' not in row and 'CD3e' in row:
					row['CD3'] = row['CD3e']
				for combo in combos:
					filter_pos, filter_neg, proteins_pos, proteins_neg = combo
					levels = levels_per_combo[combo]
					
					if skip_filter or calc_level('coexist', row, filter_pos, filter_neg):
						try:
							levels[sample].append(calc_level(o.calc, row, proteins_pos, proteins_neg))
						except KeyError:
							#raise
							continue
	
	
	# load time points etc
	time_since_baseline_per_sample = defaultdict(dict)
	if o.simpletable:
		for li, p in enumerate(splitlines(o.patienttable)):
			if p[1] in category_to_number:
				time_since_baseline_per_sample[p[0]][0] = category_to_number[p[1]]
			else:
				time_since_baseline_per_sample[p[0]][0] = float(p[1])
	elif o.patienttable:
		pylab.xlabel(o.x_column)
		for li, p in enumerate(splitlines(o.patienttable)):
			if li == 0:
				xtime_i = p.index(o.x_column)
				patient_i = p.index('CMM_ID')
				visit_i = p.index('Visit_ID')
				rtime_i = p.index('Time_from_baseline_months')
			else:
				visit = int(p[visit_i])-1
				patient = p[patient_i]
				rtime = float(p[rtime_i])
				
				if p[xtime_i] in category_to_number:
					xtime = category_to_number[p[xtime_i]]
				else:
					try: xtime = float(p[xtime_i])
					except ValueError:
						#if o.vocal and patient:
						#	print 'Missing data:', patient, visit
						continue
				
				
				if visit not in time_since_baseline_per_sample[patient]:
					time_since_baseline_per_sample[patient][visit] = xtime
	
	pvals = []
	for combo in combos:
		# place samples by patients and summarize across cells
		
		filter_pos, filter_neg, proteins_pos, proteins_neg = combo
		if not skip_filter:
			comboname = 'req' + ''.join([p+'+' for p in filter_pos]) + ''.join([p+'-' for p in filter_neg]) +'_' + ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
		else:
			comboname = ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
		levels = levels_per_combo[combo]
		locations_x = []
		locations_y = []
		samples = sorted(levels.keys())
		patient_samples = defaultdict(list)
		non_norm_y = []
		for sample in samples:
			patient_samples[sample.split('_time')[0]].append(sample)
		for patient in sorted(patient_samples.keys()):
			samples = patient_samples[patient]
			if o.cellnumbers is not None:
				samples = [sample for sample in samples if sample in cellnumbers]
			xarr = [int(sample.split('_time')[1].split('_')[0].split('.')[0]) for sample in samples]
			if o.patienttable:
				
				try:
					earliest_timepoint = min(time_since_baseline_per_sample[patient].keys()) # this does not take into account which samples are included with -i, so will choose baseline even if it's not an included sample
					if o.divxbaseline:
						xarr = [time_since_baseline_per_sample[patient][visit]/time_since_baseline_per_sample[patient][earliest_timepoint] for visit in xarr]
					else:
						xarr = [time_since_baseline_per_sample[patient][visit] for visit in xarr]
				except KeyError:
					print patient, time_since_baseline_per_sample[patient]
					print earliest_timepoint
					raise
				except ValueError:
					pass
			
			yarr = [{'mean':numpy.mean, 'sum':sum}[o.acrosscells](levels[sample]) for sample in samples]
			if o.cellnumbers is not None:
				if o.acrosscells != 'mean': raise Exception
				yarr = [numpy.mean(levels[sample])*cellnumbers[sample] for sample in samples]
			else:
				yarr = [{'mean':numpy.mean, 'sum':sum}[o.acrosscells](levels[sample]) for sample in samples]
			non_norm_y.extend(yarr)
			if o.patientnorm:
				mid_y = numpy.mean(yarr)
				yarr = [y - mid_y for y in yarr]
			locations_x.extend(xarr)
			locations_y.extend(yarr)
		
		
		
		# calculate P-values
		abundance = numpy.mean(non_norm_y) # different from _v3
		for distname in o.peakpoint:
			distfunc = distfunctions[distname]
			if pairwise_test:
				if distname == '2groups':
					categories = set(locations_x)
					if not len(categories) == 2:
						continue
					g1,g2 = sorted(categories)
					V1 = [y for x,y in zip(locations_x, locations_y) if x==g1 and not math.isnan(y)]
					V2 = [y for x,y in zip(locations_x, locations_y) if x==g2 and not math.isnan(y)]
					try:
						stat,p = correlation(V2, V1) # changed from v6
					except  ZeroDivisionError:
						stat = float('nan')
						p = float('nan')
					except ValueError:
						# comes when len(V1)!=len(V2), maybe because some samples have zero cells meeting the filter requirements?
						print >>sys.stderr, combo
						stat = float('nan')
						p = float('nan')
					try:stat = float(stat) # to deal with 1-element array values
					except TypeError: pass
					n = len(V1)
					if o.output_summary_values:
						n = repr(zip(V1, V2))
			else:
				transformed_x, clean_y = zip(*[(distfunc(x),y) for x,y in zip(locations_x, locations_y) if not math.isnan(x) and not math.isnan(distfunc(x))])
				stat, p = correlation(transformed_x, clean_y)
				n = len(clean_y)
			if str(p) == 'nan':
				#print combo, len(V1), len(V2), V1[:10], V2[:10]
				continue
			pvals.append(PVal(p, dr_tools.join(comboname, distname, str(stat), abundance, n))) # different from _v4: n
			if isinstance(stat, tuple):
				pvals[-1].r = numpy.mean(stat)
			else:
				pvals[-1].r = stat
	
	# false discovery rate and output
	for test_inst, q in zip(pvals, dr_tools.globalFDR([test_inst.p for test_inst in pvals])):
		test_inst.q = q
	for test_inst in sorted(pvals, key=lambda obj: (obj.p, -abs(obj.r)), reverse=False):
		if test_inst.q < o.maxq or o.maxq >= 1: print test_inst
	
	if o.saycombinedP:
		print 'combined P:', dr_tools.combinedP([test_inst.p for test_inst in pvals])
