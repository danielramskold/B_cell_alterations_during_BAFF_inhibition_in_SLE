from __future__ import division
import argparse, dr_tools

class Sample:
	def __init__(self, numname, descrname, cells_total):
		self.names = [descrname, numname]
		self.cells_total = cells_total
		self.cells_cytof_all = None
		self.cells_cytof_B = None
		if o.vocal: print self.names
	
	def has_all_info(self):
		return not any(n is None for n in(self.cells_total, self.cells_cytof_all, 
self.cells_cytof_B))
	
	def est_Bcells(self):
		return self.cells_total * self.cells_cytof_B / self.cells_cytof_all

def count_cells(csv_path):
	count = 0
	with open(csv_path, 'rU') as infh:
		for line in infh:
			count += 1
	return count-1

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allcells_csv', nargs='+', required=True)
	parser.add_argument('-b', '--Bcells_csv', nargs='+', required=True)
	parser.add_argument('-c', '--cellnumbertable', required=True)
	parser.add_argument('--vocal', action='store_true')
	o = parser.parse_args()
	
	samples = []
	for li, p in enumerate(dr_tools.splitlines(o.cellnumbertable)):
		samples.append(Sample(str(li+1), p[0], int(p[1]) if p[1] else None))
	
	for csv_path in o.allcells_csv:
		name = csv_path.split('/')[-1].split('_1.')[0].split('.txt')[0]
		matching_samples = [s for s in samples if name in s.names]
		if len(matching_samples) > 1: raise Exception
		if len(matching_samples) == 0:
			if o.vocal:
				print name, 'A'
			continue
		matching_samples[0].cells_cytof_all = count_cells(csv_path)
	
	for csv_path in o.Bcells_csv:
		name = csv_path.split('/')[-1].split('.')[0]
		matching_samples = [s for s in samples if name in s.names]
		if len(matching_samples) > 1: raise Exception
		elif len(matching_samples) == 0:
			print name, 'B'
			continue
		matching_samples[0].cells_cytof_B = count_cells(csv_path)
	
	for sample in samples:
		if sample.has_all_info():
			print dr_tools.join(sample.names[0], sample.est_Bcells())
		elif o.vocal:
			print sample.names, 'C'
