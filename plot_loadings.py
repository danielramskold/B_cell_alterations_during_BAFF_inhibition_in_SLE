#! python3
import argparse, pylab, itertools, dr_tools
from collections import defaultdict

def get_loadings(filepath):
	loadingsdict = dict((key, float(val)) for key,val in dr_tools.splitlines(filepath))
	name = filepath.split('/')[-1].rsplit('_loadings_',1)[-1].split('.')[0]
	return loadingsdict, name

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('loadingsfile', nargs='+')
	o = parser.parse_args()
	
	namecolourdict = defaultdict(dr_tools.randomcolour)
	
	for file1, file2 in itertools.zip_longest(o.loadingsfile[::2], o.loadingsfile[1::2]):
		loadings1, name1 = get_loadings(file1)
		
		if file2 is None:
			loadings2 = defaultdict(lambda: 0)
			name2 = 'none'
		else:
			loadings2, name2 = get_loadings(file2)
			if set(loadings1.keys()) != set(loadings2.keys()): raise Exception
		
		nameout = '%s_plotloadings_%s_%s.pdf'%(file1.rsplit('_loadings_',1)[0], name1, name2)
		
		pylab.clf()
		pylab.axhline(0, color='#eeeeee')
		pylab.axvline(0, color='#eeeeee')
		colnames = list(loadings1.keys())
		xarr = [loadings1[c] for c in colnames]
		yarr = [loadings2[c] for c in colnames]
		for x, y, c in zip(xarr, yarr, colnames):
			pylab.text(x,y,c.replace('.','-'), fontsize=8, color=namecolourdict[c])
		for x, y, c in zip(xarr, yarr, colnames):
			pylab.plot([x], [y], '.', color=namecolourdict[c])
		pylab.xlabel(name1+' loadings')
		pylab.ylabel(name2+' loadings')
		pylab.savefig(nameout)
