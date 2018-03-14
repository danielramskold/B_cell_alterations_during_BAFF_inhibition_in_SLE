import dr_tools, argparse, itertools, networkx, numpy, math, random
from collections import *
from matplotlib import pyplot, colors
import scipy.cluster.hierarchy as hcluster

def accentuate(diff):
	return diff
	return math.copysign(abs(diff)**0.5, diff)

def absmax(V):
	return max(map(abs, V))

def signmax(x):
	return max(min(x), max(x), key=abs)

def list_duplicates(l):
    # http://stackoverflow.com/questions/9835762/find-and-list-duplicates-in-python-list
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in l if x in seen or seen_add(x) )
    # turn the set into a list (as requested)
    return list( seen_twice )

def sign(v): return cmp(v, 0)

def make_distmatrix(protein_order, correlations_together, correlations_opposing, addl_func=float, twosided=False, metrics='difference', diagonal=0):
	distmatrix = numpy.zeros((len(protein_order), len(protein_order)))
	for p1i, p1 in enumerate(protein_order):
		for p2i, p2 in enumerate(protein_order):
			C = frozenset([p1, p2])
			P1 = frozenset([p1])
			P2 = frozenset([p2])
			distance = None
			if len(C) == 1:
				distance = diagonal
			elif twosided and p1i<p2i: # top right
				if metrics == 'numbers':
					distance = 1 - addl_func(signmax(correlations_opposing[C]))
				elif metrics == 'difference':
					distance = 1 - addl_func(absmax(correlations_together[C]) - max(absmax(correlations_together[P]) for P in (P1,P2)))
				#distance = addl_func(max(correlations_opposing[C]))
				#distance = 1 - addl_func(max(correlations_together[C]) - max(max(correlations_together[P]) for P in (P1,P2)))
				#distance = 1 - addl_func(max(correlations_together[C]) - max(correlations_opposing[C]))
			elif twosided: # bottom left
				if metrics == 'numbers':
					distance = 1 - addl_func(signmax(correlations_together[C]))
				elif metrics == 'difference':
					distance = 1 - addl_func(absmax(correlations_together[C]) - absmax(correlations_opposing[C]))
				
			else:
				if metrics == 'numbers':
					distance = 1 - addl_func(absmax(correlations_opposing[C]))
				elif metrics == 'difference':
					distance = 1 - addl_func(absmax(correlations_together[C]) - absmax(correlations_opposing[C]))
				#distance = 1 - addl_func(max(correlations_together[C]) - max(max(correlations_together[P]) for P in (P1,P2)))
				#distance = 1 - addl_func(max(correlations_together[C]))
			distmatrix[p1i][p2i] = distance
	return distmatrix

def make_markers_matrix(protein_order, correlations_together, correlations_opposing, qvalues_together, qvalues_opposing):
	import rpy2.robjects as robjects
	R = robjects.r
	R(r"library(psych)")
	
	cutoff = 0.2
	diff_cutoff = 0.1
	distmatrix = numpy.zeros((len(protein_order), len(protein_order)))
	for p1i, p1 in enumerate(protein_order):
		for p2i, p2 in enumerate(protein_order):
			marker = 0
			C = frozenset([p1, p2])
			P1 = frozenset([p1])
			P2 = frozenset([p2])
			try:
				if len(C) == 1:
					if abs(signmax(correlations_together[C])) > 0.05 or min(qvalues_together[C]) < 0.05:
						marker = 3 if sign(signmax(correlations_together[C])) == 1 else 2
				else:
					pvalue_different_r = list(R(r"paired.r(%f, %f, n=%i)[['p']]"%(absmax(correlations_together[C]), absmax(correlations_opposing[C]), o.n)))[0] # I checked with http://vassarstats.net/rdiff.html and it gives the same p-value (twotailed)
					if pvalue_different_r < 0.05:
						marker = 1
			except:
				print C
				raise
			distmatrix[p1i][p2i] = marker
	return distmatrix

def build_edges_1(protein_order, correlations_together, correlations_opposing, G):
	#print correlations_together, correlations_opposing
	for p1i, p1 in enumerate(protein_order):
		for p2i, p2 in enumerate(protein_order):
			C = frozenset([p1, p2])
			P1 = frozenset([p1])
			P2 = frozenset([p2])
			#print C, correlations_together[C], correlations_opposing[C]
			cutoff = 0.15
			if len(C) == 1:
				if correlations_together[C] >= cutoff:
					G.add_node(p1)
			elif C in correlations_together:
				if max(correlations_together[C]) - max(correlations_opposing[C]) > cutoff:
					G.add_edge(p1, p2)
			else:
				print C

def build_edges(protein_order, correlations_together, correlations_opposing, G):
	for p1i, p1 in enumerate(protein_order):
		listing = []
		for p2i, p2 in enumerate(protein_order):
			C = frozenset([p1, p2])
			P1 = frozenset([p1])
			P2 = frozenset([p2])
			cutoff = 0.2
			if len(C) == 1:
				if correlations_together[C] >= cutoff:
					pass#G.add_node(p1)
			elif C in correlations_together:
				d = max(correlations_together[C]) - max(correlations_opposing[C])
				if d >= cutoff:
					listing.append((d, random.random(), p2))
		#for d, rand, p2 in sorted(listing)[:2]:
		for d, rand, p2 in sorted(listing):
			G.add_edge(p1, p2)


def split_pos_neg(protline):
	proteins_pos = []
	proteins_neg = []
	pos_split = protline.split('+')
	for pos in pos_split[:-1]:
		proteins_pos.append(pos)
	neg_split = pos_split[-1].split('-')
	ni = 0
	while ni < len(neg_split)-1:
		if neg_split[ni]+'-'+neg_split[ni+1] in ('HLA-DR', 'Ki-67', 'CTLA-4', 'PD-1'):
			proteins_neg.append(neg_split[ni]+'-'+neg_split[ni+1])
			ni+=2
		else:
			proteins_neg.append(neg_split[ni])
			ni+=1
	return proteins_pos, proteins_neg

def flatten(matrix):
	return [a for b in matrix for a in b]


def secondstep_even(max_pos):
	return flatten((a+0.5, a+0.5) for a in range(-1, max_pos))[:-1]

def secondstep_odd(max_pos):
	return flatten((a-0.5, a+0.5) for a in range(0, max_pos+1))[:-1]

def markmiddle():
	pyplot.autoscale(False)
	pyplot.plot(secondstep_odd(len(protein_order)), secondstep_even(len(protein_order)), 'k-', lw=3, linewidth=0.5)
	pyplot.plot(secondstep_even(len(protein_order)), secondstep_odd(len(protein_order)), 'k-', lw=3, linewidth=0.5)

if '__main__' == __name__:
	parser =argparse.ArgumentParser()
	parser.add_argument('file_in', nargs='+')
	parser.add_argument('--proteins', nargs='+')
	parser.add_argument('--no_clustering', action='store_true')
	parser.add_argument('--coexpression', action='store_true')
	parser.add_argument('--peakpoint')
	parser.add_argument('-n', type=int, default=97) # 97 for both round1 and round2
	parser.add_argument('--max_dm', type=float)
	parser.add_argument('--not_clusterabs', action='store_false', dest='clusterabs')
	parser.add_argument('--reverse_order', action='store_true')
	o = parser.parse_args()
	o.network = False
	
	proteins = set()
	correlations_together = defaultdict(lambda: defaultdict(list))
	correlations_opposing = defaultdict(lambda: defaultdict(list))
	qvalues_together = defaultdict(lambda: defaultdict(list))
	qvalues_opposing = defaultdict(lambda: defaultdict(list))
	for file_in in o.file_in:
		for p in dr_tools.splitlines(file_in):
			if p[0].startswith('(('): continue
			peakpoint = p[3]
			r = float(p[4])
			protline = p[2]
			
			if o.coexpression:
				r = float(p[5])
			
			if protline.startswith('req'):
				filter_pos, filter_neg = split_pos_neg(protline.split('_')[0][3:])
				proteins_pos, proteins_neg = split_pos_neg(protline.split('_')[1])
			else:
				proteins_pos, proteins_neg = split_pos_neg(protline)
			
			qvalue_col = 1
			proteins_l = frozenset(proteins_pos+proteins_neg)
			if len(proteins_pos) == 2 and len(proteins_neg) == 0:
				correlations_together[peakpoint][proteins_l].append(r)
				qvalues_together[peakpoint][proteins_l].append(float(p[qvalue_col]))
			elif len(proteins_pos) == 0 and len(proteins_neg) == 2:
				correlations_together[peakpoint][proteins_l].append(r)
				raise Exception
			elif len(proteins_pos) == 1 and len(proteins_neg) == 1:
				correlations_opposing[peakpoint][proteins_l].append(r)
				qvalues_opposing[peakpoint][proteins_l].append(float(p[qvalue_col]))
			elif len(proteins_pos) == 1 and len(proteins_neg) == 0:
				correlations_together[peakpoint][proteins_l].append(r)
				qvalues_together[peakpoint][proteins_l].append(float(p[qvalue_col]))
			proteins.update(proteins_l)
	
	for peakpoint in (correlations_opposing if o.peakpoint is None else [o.peakpoint]):
		proteins = sorted(proteins)
		if o.proteins:
			if set(o.proteins) - set(proteins):
				print 'Incorrect/missing:', set(o.proteins) - set(proteins)
			if len(set(o.proteins)) < len(o.proteins):
				print 'Duplicates in --proteins:', list_duplicates(o.proteins)
			proteins = o.proteins
		if o.no_clustering:
			protein_order = o.proteins if o.proteins else proteins
		else:
			distance_matrix = make_distmatrix(proteins, correlations_together[peakpoint], correlations_opposing[peakpoint], abs if o.clusterabs else float, False, 'difference')
			hclinks = hcluster.linkage(distance_matrix, method=('complete' if o.clusterabs else 'complete'))
			draw_order = hcluster.leaves_list(hclinks)
			protein_order = [proteins[i] for i in draw_order]
			
			pyplot.clf()
			hcluster.dendrogram(hclinks, labels=proteins)
			pyplot.savefig('dendrogram_%s.pdf'%peakpoint)
		
		if o.reverse_order: protein_order = list(reversed(protein_order))
		
		if o.network:
			pyplot.clf()
			G = networkx.Graph()
			build_edges(protein_order, correlations_together[peakpoint], correlations_opposing[peakpoint], G)
			graph_pos = networkx.spring_layout(G)
			networkx.draw_networkx_nodes(G, graph_pos, node_size=2, node_color='#cccccc')
			networkx.draw_networkx_edges(G, graph_pos, width=1, edge_color='#777777')
			networkx.draw_networkx_labels(G, graph_pos, font_size=7, font_family='Arial')
			pyplot.savefig('correlation_difference_network_%s.pdf'%peakpoint)
		else:
			distance_matrix = make_distmatrix(protein_order, correlations_together[peakpoint], correlations_opposing[peakpoint], accentuate if o.coexpression else accentuate, False, 'difference', diagonal=1)
			fig = pyplot.figure(figsize=(5, 5))
			#pyplot.imshow(distance_matrix, cmap="PRGn", interpolation="none", norm=colors.Normalize(0, 2))
			absmax_dm = max(abs(numpy.amin(distance_matrix)-1), abs(numpy.amax(distance_matrix)-1))
			if o.max_dm is None:
				pass
			elif o.max_dm >= absmax_dm:
				absmax_dm = o.max_dm
			else:
				raise Exception   # too low --max_dm
			pyplot.imshow(distance_matrix, cmap="PRGn", interpolation="nearest", norm=colors.Normalize(1-absmax_dm, 1+absmax_dm))
			pyplot.colorbar(fraction=.03)
			markmiddle()
			
			if not o.coexpression:
				marker_matrix = make_markers_matrix(protein_order, correlations_together[peakpoint], correlations_opposing[peakpoint], qvalues_together[peakpoint], qvalues_opposing[peakpoint])
				for (x, y), markernum in numpy.ndenumerate(marker_matrix):
					if markernum == 1:
						pyplot.plot([x], [y], 'w.', lw=3)
					elif markernum == 2:
						pyplot.plot([x], [y], 'k_', lw=3)
					elif markernum == 3:
						pyplot.plot([x], [y], 'k+', lw=3)
			
			pyplot.xticks(range(len(protein_order)), protein_order, rotation=90, fontsize=7)
			pyplot.yticks(range(len(protein_order)), protein_order, rotation=0, fontsize=7)
			pyplot.gca().tick_params(axis=u'both', which=u'both',length=0)
			if o.coexpression:
				pyplot.savefig('coexpression_difference_matrix_%s.png'%peakpoint)
			else:
				pyplot.savefig('correlation_difference_matrix_%s.pdf'%peakpoint)
				pyplot.savefig('correlation_difference_matrix_%s.png'%peakpoint)
		
		
		
			
			distance_matrix = make_distmatrix(protein_order, correlations_together[peakpoint], correlations_opposing[peakpoint], float if o.coexpression else float, True, 'numbers')
			fig = pyplot.figure(figsize=(5, 5))
			#pyplot.imshow(distance_matrix, cmap="PRGn", interpolation="none", norm=colors.Normalize(0, 2))
			#print min(distance_matrix)
			absmax_dm = max(abs(numpy.amin(distance_matrix)-1), abs(numpy.amax(distance_matrix)-1))
			pyplot.imshow(distance_matrix, cmap="PRGn", interpolation="nearest", norm=colors.Normalize(1-absmax_dm, 1+absmax_dm))
			markmiddle()
			pyplot.xticks(range(len(protein_order)), protein_order, rotation=90, fontsize=7)
			pyplot.yticks(range(len(protein_order)), protein_order, rotation=0, fontsize=7)
			if o.coexpression:
				pyplot.savefig('coexpression_numbers_matrix_%s.png'%peakpoint)
			else:
				pyplot.savefig('correlation_numbers_matrix_%s.png'%peakpoint)
		
		
		if o.peakpoint is not None and not o.no_clustering: print ' '.join(protein_order)
		
