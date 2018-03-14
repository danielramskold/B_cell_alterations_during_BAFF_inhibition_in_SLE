import argparse, random, os

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--csv_in', nargs='+', required=True)
	parser.add_argument('-o', '--folder_out', required=True)
	parser.add_argument('-n', type=int, required=True)
	o = parser.parse_args()
	
	random.seed(1000)
	
	if not os.path.exists(o.folder_out):
		os.mkdir(o.folder_out)
	
	for infile in o.csv_in:
		lines = 0
		with open(infile, 'rU') as infh:
			for li, line in enumerate(infh):
				lines += 1
		if lines <= o.n:
			chosen_lines = set(range(0, lines))
			print lines, infile
		else:
			chosen_lines = set(random.sample(range(1, lines), o.n))
			chosen_lines.add(0) # header
		
		outfile = os.path.join(o.folder_out, os.path.basename(infile))
		with open(infile, 'rU') as infh:
			with open(outfile, 'w') as outfh:
				for li, line in enumerate(infh):
					if li in chosen_lines:
						outfh.write(line)
		
