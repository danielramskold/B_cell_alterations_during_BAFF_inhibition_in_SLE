import phenograph # run in python3
import csv, argparse, numpy

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('accense_2DSNE', nargs='?')
	parser.add_argument('accense_csv')
	parser.add_argument('-o', '--csv_out', required=True)
	o = parser.parse_args()
	
	coordinates_ri = list()
	if o.accense_2DSNE is None:
		with open(o.accense_csv, 'rU') as infh:
			for ri, row in enumerate(csv.DictReader(infh)):
				coordinates_ri.append((float(row['SNEx']), float(row['SNEy'])))
	else:
		with open(o.accense_2DSNE, 'rU') as infh:
			for ri, row in enumerate(csv.DictReader(infh)):
				coordinates_ri.append((float(row['y1']), float(row['y2'])))
	
	data = numpy.array(coordinates_ri)
	
	communities, graph, Q = phenograph.cluster(data)
	# For a dataset of *N* rows, `communities` will be a length *N* vector of integers specifying a community assignment for each row in the data. Any rows assigned `-1` were identified as *outliers* and should not be considered as a member of any community.
	# `graph` is a *N* x *N* `scipy.sparse` matrix representing the weighted graph used for community detection. 
	# `Q` is the modularity score for `communities` as applied to `graph`.
	
	assert len(coordinates_ri) == len(communities)
	
	with open(o.accense_csv, 'rU') as infh:
		with open(o.csv_out, 'w') as outfh:
			outf_csv = csv.writer(outfh)
			for li, row in enumerate(csv.reader(infh)):
				if li == 0:
					population_i = row.index('population')
					outf_csv.writerow(row)
				else:
					row[population_i] = communities[li-1]
					outf_csv.writerow(row)
