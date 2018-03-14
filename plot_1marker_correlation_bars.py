import argparse, pylab

def get_5FDR_markers(filename):
	markers = []
	with open(filename, 'rU') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			if float(p[1]) < 0.05 or o.allP or (o.nomP and float(p[0]) < 0.05):
				markers.append(p[2].rstrip('+'))
	return markers

def ordered_union(*lists):
	together = []
	for l in lists:
		for entry in l:
			if entry not in together:
				together.append(entry)
	return together

def plot_marker_names(names, x):
	for y, name in enumerate(names):
		pylab.text(x, y, name, fontsize=o.fontsize)

def plot_correlations(filename, names, x, as_text, colourscheme, mid0, invertsign):
	marker_to_colour = dict()
	marker_to_r = dict()
	marker_to_p = dict()
	with open(filename, 'rU') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			marker = p[2].rstrip('+')
			if marker in marker_to_r: continue # might change to raisingException
			marker_to_r[marker] = float(p[4])
			pval = float(p[0])
			qval = float(p[1])
			marker_to_p[marker] = pval
			if marker_to_r[marker] > 0:
				marker_to_colour[marker] = colourscheme[2] if pval > 0.05 else colourscheme[1] if qval > 0.05 else colourscheme[0]
			else:
				marker_to_colour[marker] = colourscheme[5] if pval > 0.05 else colourscheme[4] if qval > 0.05 else colourscheme[3]
		for y, name in enumerate(names):
			if name == 'CD3' and name not in marker_to_colour: name = 'CD3e'
			if as_text:
				pylab.text(x, y, '%.2f'%marker_to_r[name], color=marker_to_colour[name])
			elif o.mid0:
				pylab.barh(y, marker_to_r[name]/(-2.0 if invertsign else 2.0), color=marker_to_colour[name], linewidth=0, left=x+0.5)
				
			else:
				pylab.barh(y, abs(marker_to_r[name]), color=marker_to_colour[name], linewidth=0, left=x)
				if marker_to_p[name] <= 0.05:
					if marker_to_r[name] >= 0:
						pylab.text(x+ abs(marker_to_r[name]), y+0.20, '+')
					else:
						pylab.text(x+ abs(marker_to_r[name]) + 0.01, y+0.20, '-')

def get_coordinates(filename):
	marker_to_r = dict()
	with open(filename, 'rU') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			marker = p[2].rstrip('+')
			if marker in marker_to_r: continue # might change to raisingException
			marker_to_r[marker] = float(p[4])
	return marker_to_r

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--table_in', required=True, nargs='+')
	parser.add_argument('-o', '--pdf_out', default='1marker_correlations.pdf')
	parser.add_argument('-P', '--proteins', nargs='+')
	parser.add_argument('--textonly', action='store_true')
	parser.add_argument('-c', '--colourscheme', choices=['fade', 'redshift', 'faderev'], default='fade')
	parser.add_argument('--allP', action='store_true')
	parser.add_argument('--nomP', action='store_true')
	parser.add_argument('--fontsize', type=int, default=20)
	parser.add_argument('--scatterplot', action='store_true')
	parser.add_argument('--order_by_direction', action='store_true')
	parser.add_argument('-m', '--mid0', '--minustoplus', action='store_true')
	parser.add_argument('-I', '--invertsign', '--invertplusminus', action='store_true')
	o = parser.parse_args()
	
	if o.proteins is not None:
		significant_markers = list(reversed(o.proteins))
	else:
		significant_markers = list(reversed(ordered_union(*map(get_5FDR_markers, o.table_in))))
	
	if o.order_by_direction:
		marker_to_r = get_coordinates(o.table_in[0])
		significant_markers = sorted(significant_markers, key=lambda g: -marker_to_r[g])
	
	colours = {'fade':['#ff0000', '#dd8060', '#eec4a0', '#0000ff', '#6080dd', '#a0c4ee'], 'redshift':['#000077', '#9900dd', '#dd88ff', '#007700', '#77dd00', '#c8ff22'], 'faderev':['#0000e1', '#6089dd', '#b0d4ee', '#e10000', '#dd8960', '#eed4b0'],}[o.colourscheme]
	
	if o.scatterplot:
		coordsD1 = get_coordinates(o.table_in[0])
		coordsD2 = get_coordinates(o.table_in[1])
		xarr, yarr = [coordsD1[m] for m in significant_markers], [coordsD2[m] for m in significant_markers]
		pylab.plot(xarr, yarr, 'k.')
		for x, y, m in zip(xarr, yarr, significant_markers):
			pylab.text(x, y, m, fontsize=o.fontsize)
		pylab.xlim(-1, 1)
		pylab.ylim(-1, 1)
	else:
		xticks_pos = []
		xticks_txt = []
		plot_marker_names(significant_markers, 0)
		for x, filename in enumerate(o.table_in):
			plot_correlations(filename, significant_markers, (x*1.14+0.8), o.textonly, colours, o.mid0, o.invertsign)
			if o.mid0:
				xticks_pos.extend([x*1.14+0.8, x*1.14+1.05, x*1.14+1.3, x*1.14+1.55, x*1.14+1.8])
				if o.invertsign:
					xticks_txt.extend(['1', '0.5', '0', '-0.5', '-1'])
				else:
					xticks_txt.extend(['-1', '-0.5', '0', '0.5', '1'])
			else:
				xticks_pos.extend([x*1.14+0.8, x*1.14+1.3, x*1.14+1.8])
				xticks_txt.extend(['0', '0.5', '1'])
		pylab.xlim(0, len(o.table_in)+1)
		pylab.ylim(-0.5, len(significant_markers)+0.5)
		pylab.xticks(xticks_pos, xticks_txt)
		pylab.yticks([])
	pylab.savefig(o.pdf_out)
