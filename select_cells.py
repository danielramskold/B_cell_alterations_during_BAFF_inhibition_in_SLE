from __future__ import division
import argparse, operator, csv

comparison_op = {'>=':operator.ge, '<=':operator.le, '>':operator.gt, '<':operator.lt}
calculation_op = {'+':operator.add, '-':operator.sub}

def match_colname(colname):
	names = colname.strip().split()
	
	possibilities = set()
	for name in names:
		for colname_csv in csvheader:
			if name in colname_csv.strip('" ').split(':'):
				possibilities.add(colname_csv)
	
	if len(possibilities) == 1:
		return list(possibilities)[0]
	else:
		print possibilities, colname
		raise Exception
	

class Criterion:
	def __init__(self, line):
		for opsym_c in ['>=', '<=', '>', '<']:
			if opsym_c in line:
				parts = line.rsplit(opsym_c,1)
				self.cutoff = float(parts[-1].strip())
				self.comparison = comparison_op[opsym_c]
				break
		else:
			raise Exception, 'no cutoff operator (>=,<=,>,<)'
			
		# could use asteval module for more complicated cases: http://newville.github.io/asteval/basics.html
		# allow e.g. "Ir193Di + Ir191Di >=9" as line
		for opsym_a in ['+', '-']:
			if opsym_a in line:
				self.calculation = calculation_op[opsym_a]
				self.columns = [match_colname(n) for n in parts[0].split(opsym_a,1)]
				self.representation = ' '.join([self.columns[0], opsym_a, self.columns[1], opsym_c, str(self.cutoff)])
				break
		else:
			self.calculation = None
			self.columns = [match_colname(parts[0])]
			self.representation = ' '.join([self.columns[0], opsym_c, str(self.cutoff)])
		self.reset_counts()
	
	def matches(self, csv_row):
		try:
			if self.calculation is None:
				value = float(csv_row[self.columns[0]])
			else:
				value = self.calculation(float(csv_row[self.columns[0]]), float(csv_row[self.columns[1]]))
			result = self.comparison(value, self.cutoff)
		except KeyError:
			if o.pass_missing_marker:
				result = True
			else:
				raise
		
		if result: self.passes += 1
		else: self.fails += 1
		return result
	
	def pass_rate(self):
		tot = self.passes+self.fails
		if tot == 0: return float('nan')
		return self.passes/tot
	
	def reset_counts(self):
		self.passes = 0
		self.fails = 0

# no normalisation done, so counts in -> counts out

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--gating_list', required=True)
	parser.add_argument('-i', '--csv_in', nargs='+', required=True)
	parser.add_argument('-m', '--output_mapping')
	parser.add_argument('-H', '--imitate_header_csv')
	parser.add_argument('-P', '--proteins_only', action='store_true')
	parser.add_argument('--pass_missing_marker', action='store_true')
	o = parser.parse_args()
	
	if o.proteins_only and o.imitate_header_csv is None:
		o.imitate_header_csv = o.csv_in[0]
	
	# prepare csv column names for match_colname() called by Criterion()
	with open(o.csv_in[0], 'rU') as infh:
		for row in csv.DictReader(infh):
			csvheader = row.keys()
			break
	
	# parse gating criteria
	criteria = []
	with open(o.gating_list, 'rU') as infh:
		for line in infh:
			if line[0] == '#': continue
			if not line.strip(): continue
			criteria.append(Criterion(line))
	
	# load output name to input name mapping
	if o.output_mapping is not None:
		with open(o.output_mapping, 'rU') as infh:
			output_mapping = dict(line.rstrip().split('\t') for line in infh)
	
	# prepare changing the header
	if o.imitate_header_csv is not None:
		with open(o.imitate_header_csv, 'rU') as infh:
			for row in csv.DictReader(infh):
				model_csvheader = row.keys()
				break
		csvheader_renaming = dict() # not done
		for target_colname in model_csvheader:
			metal = target_colname.strip('" ').split(':')[0]
			try: source_colname = [n for n in csvheader if metal == n.strip('" ')][0]
			except IndexError: source_colname = target_colname
			if o.proteins_only:
				target_colname = target_colname.split(':')[-1].strip('" ')
				if not target_colname: continue
				if target_colname == source_colname: continue
				if target_colname in ['Cisplatin','EQBeads', 'DNAIr', 'Time']: continue
			csvheader_renaming[source_colname] = target_colname
	
	for csv_in in o.csv_in:
		for criterion in criteria:
			criterion.reset_counts()
		
		# output file name
		if o.output_mapping is None:
			csv_out = csv_in.split('/')[-1]
			if csv_out == csv_in:
				csv_out = 'o_' + csv_out
		else:
			csv_out = output_mapping[csv_in.split('/')[-1]]
		
		# parse, select and output
		with open(csv_in, 'rU') as infh:
			with open(csv_out, 'w') as outfh:
				if o.imitate_header_csv is None:
					writer = csv.DictWriter(outfh, fieldnames=csvheader)
				else:
					writer = csv.DictWriter(outfh, fieldnames=[csvheader_renaming[k] for k in csvheader if k in csvheader_renaming])
				writer.writeheader()
				
				for row in csv.DictReader(infh):
					try:
						if all(criterion.matches(row) for criterion in criteria):
							if o.imitate_header_csv is None:
								writer.writerow(row)
							else:
								writer.writerow(dict((csvheader_renaming[k],v) for k,v in row.items() if k in csvheader_renaming))
					except KeyError:
						print csv_in, row
						raise
		
		# say statistics
		print csv_out
		for criterion in criteria:
			print '%s : %.1f%%'%(criterion.representation, criterion.pass_rate()*100)
		print 'total %d (%.1f%%)'%(criteria[-1].passes, 100*criteria[-1].passes/(criteria[0].passes + criteria[0].fails))
