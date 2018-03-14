from __future__ import division
import argparse, numpy
from scipy import stats
from collections import defaultdict

def splitlines(infile):
	with open(infile, 'rU') as infh:
		for line in infh:
			yield line.rstrip('\r\n').split('\t')

category_to_number = {'Yes':1, 'No':0, 'Male':0, 'Female':1}

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', '--patienttable', required=True)
	parser.add_argument('-i', '--csv_list')
	parser.add_argument('-x', '--x_column', default='SLEDAI_score', metavar='SLEDAI_score')
	parser.add_argument('-y', '--y_column', default='Time_from_baseline_months', metavar='Time_from_baseline_months')
	parser.add_argument('--patientnorm', action='store_true')
	parser.add_argument('--divxbaseline', action='store_true')
	parser.add_argument('--divybaseline', action='store_true')
	o = parser.parse_args()
	o.quite = True
	
	if o.csv_list:
		allowed_patient_visits = set()
		with open(o.csv_list, 'rU') as infh:
			samples_csv = [f.strip() for f in infh.read().split()]
		for sample in samples_csv:
			patient = sample.split('/')[-1].split('_time')[0]
			visit = int(sample.split('/')[-1].split('_time')[1].split('_')[0].split('.')[0])
			allowed_patient_visits.add((patient, visit))
	
	
	
	time_since_baseline_per_sample = defaultdict(dict)
	sledai_since_baseline_per_sample = defaultdict(dict)
	if o.patienttable:
		for li, p in enumerate(splitlines(o.patienttable)):
			if li == 0:
				time_i = p.index(o.y_column)
				patient_i = p.index('CMM_ID')
				sledai_i = p.index(o.x_column)
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
				if p[sledai_i] in category_to_number:
					sledai = category_to_number[p[sledai_i]]
				else:
					try: sledai = float(p[sledai_i])
					except ValueError:
						if not o.quite and patient:
							print 'Missing data:', patient, visit
						continue
				
				
				if visit not in time_since_baseline_per_sample[patient]:
					time_since_baseline_per_sample[patient][visit] = time
					sledai_since_baseline_per_sample[patient][visit] = sledai
	
	locations_x = []
	locations_y = []
	for patient in time_since_baseline_per_sample:
		if o.csv_list:
			yarr = [time_since_baseline_per_sample[patient][visit] for visit in sorted(time_since_baseline_per_sample[patient]) if (patient, visit) in allowed_patient_visits]
			xarr = [sledai_since_baseline_per_sample[patient][visit] for visit in sorted(time_since_baseline_per_sample[patient]) if (patient, visit) in allowed_patient_visits]
			if not xarr:
				continue
		else:
			yarr = [time_since_baseline_per_sample[patient][visit] for visit in sorted(time_since_baseline_per_sample[patient])]
			xarr = [sledai_since_baseline_per_sample[patient][visit] for visit in sorted(time_since_baseline_per_sample[patient])]
		if o.divxbaseline:
			earliest_timepoint = min(sledai_since_baseline_per_sample[patient].keys())
			if not any(sledai_since_baseline_per_sample[patient].values()): # for Anti_dsDNA_IFL
				#xarr = [0 for x in xarr]
				continue
			else:
				while sledai_since_baseline_per_sample[patient][earliest_timepoint] == 0: # for Anti_dsDNA_IFL
					earliest_timepoint += 1
				xarr = [x/sledai_since_baseline_per_sample[patient][earliest_timepoint] for x in xarr]
		if o.divybaseline:
			earliest_timepoint = min(time_since_baseline_per_sample[patient].keys())
			yarr = [y/time_since_baseline_per_sample[patient][earliest_timepoint] for y in yarr]
		if o.patientnorm:
			xarr = [x-numpy.mean(xarr) for x in xarr]
		locations_x.extend(xarr)
		locations_y.extend(yarr)
	
	print stats.spearmanr(locations_x, locations_y)
	
