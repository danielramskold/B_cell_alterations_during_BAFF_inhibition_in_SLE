from __future__ import division
import argparse, csv, dr_tools, numpy, math, itertools, os, random
from collections import defaultdict, namedtuple

def log(v): return math.log(max(0.25, float(v)), 2)

def up_pos(V): return [v+3 for v in V]
def up_neg(V): return [v+3 for v in V]

def existance(V): return [int(v>cutoff) for v in V]
def nonexistance(V): return [int(not v) for v in existance(V)]

def product(V):
	ret = 1
	for v in V: ret *= v
	return ret

def combinations_maxN(pick_from, max_n, exclusion=set()):
	restricted_pick_from = set(pick_from) - set(exclusion)
	for n in range(max_n+1):
		for chosen in itertools.combinations(restricted_pick_from, n):
			yield chosen

def make_comboname(combo):
	filter_pos, filter_neg, proteins_pos, proteins_neg = combo
	if filter_pos or filter_neg:
		comboname = 'req' + ''.join([p+'+' for p in filter_pos]) + ''.join([p+'-' for p in filter_neg]) +'_' + ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
	else:
		comboname = ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
	return comboname

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
	if calc == 'coexist':
		values_pos = [float(row[p]) for p in pos]
		values_neg = [float(row[p]) for p in neg]
		return product(existance(values_pos) + nonexistance(values_neg))
	
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
	with open(infile, 'rU') as infh:
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
	stat, p = stats.ttest_rel(val1, val0)
	stat = numpy.mean(val1)-numpy.mean(val0)
	return stat, p

def perm_mean(val1, val0):
	p = dr_tools.permutationtest(numpy.mean, val1, val0, controls=1000)
	stat = numpy.mean(val1)-numpy.mean(val0)
	return p, stat


set_order_patients = '''1
1
1
1
1
1
1
2
3
3
3
3
3
3
3
3
4
4
4
4
4
4
4
4
5
5
5
5
6
6
6
6
6
6
6
6
7
7
7
7
7
7
7
7
7
9
9
9
9
9
9
9
10
10
10
10
11
11
11
11
11
12
12
12
12
12
12
13
13
13
13
13
13
15
15
15
16
16
16
16
16
16
3001
3001
3001
3001
3001
19
19
19
19
19
3002
3002
3002
3002
17
17
17
17
17
17
18
18
18
18
18
20
20
20
20
20
21
21
21
21
21
3003
3003
3003
22
22
22
22
22
23
23
23
23
23'''.split()


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
	parser.add_argument('--calc', choices=['coexist', 'sumlog', 'prodlog'], default='coexist')
	parser.add_argument('--acrosscells', choices=['mean', 'sum'], default='mean')
	parser.add_argument('--cutoff', type=float, default=0)
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
	parser.add_argument('--shuffle', action='store_true')
	parser.add_argument('--shuffle_name', action='store_true')
	o = parser.parse_args()
	
	if o.cellnumbers is not None:
		cellnumbers = dict()
		for sample, nstr in splitlines(o.cellnumbers):
			if nstr:
				cellnumbers[sample] = float(nstr)
	
	
	cutoff = o.cutoff
	
	skip_filter = not (o.filter_pos or o.filter_neg or o.set_filter)
	
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
			if o.use_from_list is not None:
				samples_csv = [f.strip() for f in infh.read().split() if f.strip() in samples_to_use]
			else:
				samples_csv = [f.strip() for f in infh.read().split()]
	
	for csv_in in samples_csv:
		with open(csv_in, 'rb') as infh:
			sample = csv_in.rsplit('/', 1)[-1].split('.csv')[0].split('_B')[0].split('_T')[0]
			for ri, row in enumerate(csv.DictReader(infh)):
				if o.maxlines and ri >= o.maxlines:break
				if not combos:
					row_keys = set(row) - set(['DNAIr', 'DNAIr.1', 'EQBeads', 'EQBeads.1', 'EQBeads.2', 'Cisplatin'])
					if o.proteins:
						row_keys = set(o.proteins)
					
					filter_generator = filter_listing_loop(row_keys, o, filter_req_proteins)
					
					for exclude_2, filter_pos, filter_neg in filter_generator:
						for proteins_pos in combinations_maxN(row_keys, o.num_pos, exclude_2):
							exclude_3 = exclude_2|set(proteins_pos)
							for proteins_neg in combinations_maxN(row_keys, o.num_neg, exclude_3):
								if o.req_proteins and (set(proteins_pos)|set(proteins_neg))&req_proteins != req_proteins:
									continue
								combos.append((filter_pos, filter_neg, proteins_pos, proteins_neg))
								levels_per_combo[combos[-1]] = defaultdict(list)
				
				for combo in combos:
					filter_pos, filter_neg, proteins_pos, proteins_neg = combo
					levels = levels_per_combo[combo]
					
					if skip_filter or calc_level('coexist', row, filter_pos, filter_neg):
						try:
							levels[sample].append(calc_level(o.calc, row, proteins_pos, proteins_neg))
						except KeyError:
							#raise
							continue
	
	
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
	
	column_order = ['CMM_ID', 'time_from_baseline_months']
	table = dict()
	for combo in combos:
		filter_pos, filter_neg, proteins_pos, proteins_neg = combo
		if not skip_filter:
			comboname = 'req' + ''.join([p+'+' for p in filter_pos]) + ''.join([p+'-' for p in filter_neg]) +'_' + ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
		else:
			comboname = ''.join([p+'+' for p in proteins_pos]) + ''.join([p+'-' for p in proteins_neg])
		levels = levels_per_combo[combo]
		locations_x = []
		locations_y = []
		patient_by_sample = []
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
			patient_by_sample.extend([patient for x in xarr])
		
		if o.shuffle:
			random.shuffle(locations_y)
		
		table[comboname] = locations_y
		table['time_from_baseline_months'] = locations_x
		table['CMM_ID'] = patient_by_sample
		if o.shuffle_name:
			random.shuffle(table['CMM_ID'])
		column_order.append(comboname)
		
		
		
	print dr_tools.join(column_order)
	transposed_table = zip(*[table[c] for c in column_order])
	
	last_patient = ''
	for patient in set_order_patients:
		if patient == last_patient:
			print dr_tools.join(patient)
		else:
			for row in transposed_table:
				if row[0] == patient:
					print dr_tools.join(row)
					break
			else:
				print patient
		last_patient = patient
	
	
	
	#for row in transposed_table:
	#	print dr_tools.join(row)
