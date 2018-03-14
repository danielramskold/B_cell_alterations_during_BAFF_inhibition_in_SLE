#! python3

grep_output = '''pcaBrun2_loadings_PC1.txt:CD57	0.00213289704675
pcaBrun2_loadings_PC10.txt:CD57	-0.00467347239758
pcaBrun2_loadings_PC11.txt:CD57	0.00938763731316
pcaBrun2_loadings_PC12.txt:CD57	0.00843216341306
pcaBrun2_loadings_PC13.txt:CD57	0.00491552467384
pcaBrun2_loadings_PC14.txt:CD57	0.0106112959522
pcaBrun2_loadings_PC15.txt:CD57	-0.00891889465053
pcaBrun2_loadings_PC16.txt:CD57	-0.0103199345092
pcaBrun2_loadings_PC17.txt:CD57	0.00473435509884
pcaBrun2_loadings_PC18.txt:CD57	0.00353951289479
pcaBrun2_loadings_PC19.txt:CD57	-0.00969625567717
pcaBrun2_loadings_PC2.txt:CD57	-0.00168987384451
pcaBrun2_loadings_PC20.txt:CD57	0.00687262941844
pcaBrun2_loadings_PC21.txt:CD57	0.0126756237955
pcaBrun2_loadings_PC22.txt:CD57	0.0169472538713
pcaBrun2_loadings_PC23.txt:CD57	-0.0305180422953
pcaBrun2_loadings_PC24.txt:CD57	-0.0215581518625
pcaBrun2_loadings_PC25.txt:CD57	0.075819834444
pcaBrun2_loadings_PC26.txt:CD57	-0.0124826379954
pcaBrun2_loadings_PC27.txt:CD57	-0.0593580004791
pcaBrun2_loadings_PC28.txt:CD57	0.000319397912724
pcaBrun2_loadings_PC29.txt:CD57	-0.00473694255768
pcaBrun2_loadings_PC3.txt:CD57	0.00228368577093
pcaBrun2_loadings_PC30.txt:CD57	0.993869869655
pcaBrun2_loadings_PC4.txt:CD57	-0.000332621659945
pcaBrun2_loadings_PC5.txt:CD57	-0.0037767560514
pcaBrun2_loadings_PC6.txt:CD57	-0.00218052018599
pcaBrun2_loadings_PC7.txt:CD57	0.00525855810522
pcaBrun2_loadings_PC8.txt:CD57	-0.00156770252581
pcaBrun2_loadings_PC9.txt:CD57	0.0136523917998'''
# from round2_cytof/PCA$ grep CD57 pcaBrun2_loadings_PC*txt

relcontrib_text = '''0.130778033803
0.0948222809751
0.0683424098086
0.0652072631226
0.0501814605215
0.0460346095117
0.0436089920985
0.0430033888921
0.0385127582143
0.0348960151937
0.0337891713553
0.0313004776818
0.0306874249908
0.0301318700594
0.0294141203858
0.0267216435451
0.0255711423014
0.0221040373167
0.0208621446639
0.0167412624558
0.0151454206116
0.0142016390834
0.0137898065158
0.013665930192
0.0130735944307
0.0116570111332
0.0108802979659
0.00943842130502
0.00856335306815
0.00687401879705'''
# from round2_cytof/PCA$ open pcaBrun2_relconrib.txt

from matplotlib import pyplot
splitlines = [line.rstrip().split('\t') for line in grep_output.split('\n')]
heightsD = {int(p[0].split('PC')[1].split('.txt')[0]) : float(p[1]) for p in splitlines}
order = list(range(1, 30+1))
relcontrib = [float(line.strip()) for line in relcontrib_text.split('\n')]
pyplot.plot([x+0.4 for x in order], relcontrib, 'r-')
pyplot.bar(order, [heightsD[v] for v in order], linewidth=0)
pyplot.savefig('CD57_loading_per_PC.pdf')
