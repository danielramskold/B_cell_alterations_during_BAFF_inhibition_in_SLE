#! python3
import argparse, pandas, dr_tools, os, random
from matplotlib import pyplot

def shades_gen():
	for c in ['#ff0000', '#0000ff', '#00ff00']:
		yield c
	while True:
		yield dr_tools.randomcolour()

def style_gen():
	for s in ['-.', '--', ':', '-']: yield s
	raise Exception


if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--csv_in', nargs='*', action='append', required=True)
	parser.add_argument('-m', '--marker', '--channel', default='CD57')
	parser.add_argument('-o', '--plot_out', default='distribution.pdf')
	parser.add_argument('--base_folder_in', default='')
	parser.add_argument('--xlim', type=float, default=3)
	parser.add_argument('--ylim', type=float)
	parser.add_argument('--step', type=float, default=0.01)
	parser.add_argument('--noncumulative', action='store_true')
	o = parser.parse_args()
	
	random.seed(80)
	
	for files_in, bg_shade in zip(o.csv_in, shades_gen()):
		for csv_in in files_in:
			df = pandas.read_csv(os.path.join(o.base_folder_in, csv_in))
			values = list(df[o.marker])
			
			xarr, yarr = dr_tools.bin(values, start=o.step, end=o.xlim, step=o.step, cumulative=0 if o.noncumulative else -1, fractions=True)
			pyplot.plot(xarr, [y*100 for y in yarr], '-', color=dr_tools.mixcolours([bg_shade, dr_tools.randomcolour()], [0.5, 0.5]))
	
	pyplot.xlabel('Signal' if o.noncumulative else 'Cutoff')
	pyplot.ylabel('% of cells in bin'  if o.noncumulative else '% of cells above cutoff')
	if o.ylim is not None: pyplot.ylim(0, o.ylim)
	pyplot.title(o.marker)
	pyplot.savefig(o.plot_out)
