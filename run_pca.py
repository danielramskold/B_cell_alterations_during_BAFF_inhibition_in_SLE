from __future__ import division
import pylab, argparse, numpy, dr_tools, math
from rpy2.robjects import r as R
from collections import defaultdict

def listRightIndex(alist, value):
    # http://stackoverflow.com/questions/9836425/equivelant-to-rindex-for-lists-in-python
    return len(alist) - alist[-1::-1].index(value) -1

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('merged_csv')
	parser.add_argument('-o', '--out_prefix', required=True)
	parser.add_argument('--log', const=0.1, type=float, nargs='?')
	parser.add_argument('--eachsample', action='store_true')
	parser.add_argument('--plotformat', choices=['pdf', 'png'], default='pdf')
	parser.add_argument('--marker', default='.')
	parser.add_argument('--eachpatient', action='store_true')
	parser.add_argument('--timecorr', action='store_true')
	parser.add_argument('--generaloutput', '--general', action='store_true')
	parser.add_argument('--nolegend', action='store_true')
	parser.add_argument('--eachchannel', action='store_true')
	parser.add_argument('--shown_samples', nargs='+')
	parser.add_argument('--alpha', '--opacity', type=float, default=1)
	o = parser.parse_args()
	
	R(r"dat = read.csv('%s', header=T, row.names=1)"%o.merged_csv)
	if o.log is not None: R(r"dat = log(dat+%f)"%o.log)
	R(r"pcaresult = prcomp(dat)")
	names = list(R(r"row.names(dat)"))
	#coordinates = numpy.array(R(r"pcaresult$x"))
	#print coordinates.shape # (32875, 31)
	sample_per_row = [n.split('c')[0] for n in names]
	samples = sorted(set(sample_per_row)) if o.shown_samples is None else o.shown_samples
	
	for components in [(1, 2), (3, 4), (5,6), (7,8)]:
		pylab.clf()
		xarr = list(R(r"pcaresult$x[,%i]"%components[0]))
		yarr = list(R(r"pcaresult$x[,%i]"%components[1]))
		if o.generaloutput:
			for sample in samples:
				first_i = sample_per_row.index(sample)
				last_i = listRightIndex(sample_per_row, sample)
				pylab.plot(xarr[first_i:last_i+1], yarr[first_i:last_i+1], o.marker, label=sample, alpha=o.alpha)
			pylab.xlabel('PC%i'%components[0])
			pylab.ylabel('PC%i'%components[1])
			pylab.savefig('%s_%i_%i.%s'%(o.out_prefix, components[0], components[1], o.plotformat))
		if o.eachsample:
			for sample_show in samples:
				pylab.clf()
				'''for sample in samples:
					if sample.split('t')[0] != sample_show.split('t')[0]:
						first_i = sample_per_row.index(sample)
						last_i = listRightIndex(sample_per_row, sample)
						pylab.plot(xarr[first_i:last_i+1], yarr[first_i:last_i+1], '.', label=sample, color='#cccccc')'''
				for sample in samples:
					if sample.split('t')[0] == sample_show.split('t')[0] and sample != sample_show:
						first_i = sample_per_row.index(sample)
						last_i = listRightIndex(sample_per_row, sample)
						pylab.plot(xarr[first_i:last_i+1], yarr[first_i:last_i+1], o.marker, label=sample, color='#888888', alpha=o.alpha)
				for sample in samples:
					if sample == sample_show:
						first_i = sample_per_row.index(sample)
						last_i = listRightIndex(sample_per_row, sample)
						pylab.plot(xarr[first_i:last_i+1], yarr[first_i:last_i+1], o.marker, label=sample, color='#000000', alpha=o.alpha)
				pylab.xlabel('PC%i'%components[0])
				pylab.ylabel('PC%i'%components[1])
				pylab.savefig('%s_%i_%i_%s.%s'%(o.out_prefix, components[0], components[1], sample_show, o.plotformat))
		if o.eachpatient:
			patients = sorted(set(sample.split('t')[0] for sample in samples))
			for patient in patients:
				pylab.clf()
				for sample in samples:
					if sample.split('t')[0] == patient:
						first_i = sample_per_row.index(sample)
						last_i = listRightIndex(sample_per_row, sample)
						try: timepoint = float(sample.split('t')[1].split('T')[0].split('B')[0])
						except: print sample; raise
						#colour = dr_tools.mixcolours(['#ff0000', '#0000ff'], [timepoint/6, 1-timepoint/6])
						colour = dr_tools.mixcolours(["#eeffcc", "#0000ee"], [timepoint/6, 1-timepoint/6])
						pylab.plot(xarr[first_i:last_i+1], yarr[first_i:last_i+1], o.marker, label=sample, color=colour, alpha=o.alpha) # would be better to colour e.g. from red (t0) to blue (t5)
				pylab.xlabel('PC%i'%components[0])
				pylab.ylabel('PC%i'%components[1])
				if not o.nolegend: pylab.legend()
				pylab.savefig('%s_%i_%i_%s.%s'%(o.out_prefix, components[0], components[1], patient, o.plotformat))
		if o.eachchannel:
			genes =  list(R(r"colnames(dat)"))
			for channel in genes:
				pylab.clf()
				shades = [dr_tools.mixcolours(["#eeffcc", "#0000ee"], [1-f, f]) for f in pylab.arange(0, 1, 1/10)]
				shadebinfunc = lambda v: int(min(len(shades)-1, 0 if v== 0 else max(1, math.log(v,2)+2)))
				colvalues = list(R('dat[["%s"]]'%channel))
				shadebins = defaultdict(list)
				for i, v in enumerate(colvalues):
					shadebins[shadebinfunc(v)].append(i)
				for shadenum in sorted(shadebins.keys()):
					iarr = shadebins[shadenum]
					pylab.plot([xarr[i] for i in iarr], [yarr[i] for i in iarr], o.marker, color=shades[shadenum], label='0' if shadenum == 0 else '0-1' if shadenum == 1 else  '>128' if shadenum == 9 else '%i-%i'%(2**(shadenum-2), 2**(shadenum-1)), alpha=o.alpha)
				pylab.xlabel('PC%i'%components[0])
				pylab.ylabel('PC%i'%components[1])
				if not o.nolegend: pylab.legend()
				pylab.savefig('%s_%i_%i__%s.%s'%(o.out_prefix, components[0], components[1], channel, o.plotformat))
	
	if o.timecorr:
		from scipy import stats
		distfunctions = {'across':float, 'mid':(lambda t: 2.5-abs(2.5-t)), 'start':(lambda t: min(2,t)), 'end':(lambda t: 2-min(2,5-t))}
		for component in range(1, 9):
			pos_arr = list(R(r"pcaresult$x[,%i]"%component))
			patients = sorted(set(sample.split('t')[0] for sample in samples))
			whole_posarr = []
			whole_timearr = []
			for patient in patients:
				patient_posarr = []
				patient_timearr = []
				for sample in samples:
					if sample.split('t')[0] == patient:
						first_i = sample_per_row.index(sample)
						last_i = listRightIndex(sample_per_row, sample)
						try: timepoint = float(sample.split('t')[1].split('T')[0].split('B')[0])
						except: print sample; raise
						patient_posarr.extend(pos_arr[first_i:last_i+1])
						patient_timearr.extend([timepoint for i in range(first_i, last_i+1)])
				whole_posarr.extend(patient_posarr)
				whole_timearr.extend(patient_timearr)
				
				for distname, distfunc in distfunctions.items():
					print component, patient, distname, stats.pearsonr(patient_posarr, map(distfunc, patient_timearr))
			for distname, distfunc in distfunctions.items():
				print component, 'whole', distname, stats.pearsonr(whole_posarr, map(distfunc, whole_timearr))
	
	if o.generaloutput:
		sdev = list(R(r"pcaresult$sdev"))
		with open(o.out_prefix + '_relconrib.txt', 'w') as outfh:
			for v in [v/sum(sdev) for v in sdev]:
				print >>outfh, v
		for component in ['PC%i'%(i+1) for i in range(len(sdev))]:
			loadings = list(R(r"pcaresult$rotation[,'%s']"%component))
			genes =  list(R(r"colnames(dat)"))
			with open(o.out_prefix + '_loadings_%s.txt'%component, 'w') as outfh:
				for l,g in sorted(zip(loadings, genes)):
					print >>outfh, str(g)+'\t'+str(l)
