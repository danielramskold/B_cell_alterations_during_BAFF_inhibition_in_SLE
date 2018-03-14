library(psych)
library(argparser)
p = arg_parser('')
p = add_argument(p, 'neg_file', help='')
p = add_argument(p, 'pos_file', help='')
p = add_argument(p, 'out_file', help='')
p = add_argument(p, '--all_q', flag=TRUE, help='')
p = parse_args(p)

negdata = read.table(p$neg_file, header=F, sep='\t', row.names=3)
posdata = read.table(p$pos_file, header=F, sep='\t', row.names=3)
bothdata = merge(posdata, negdata, by=0)
if(!p$all_q) {
  significant_somewhere = bothdata$V2.x < 0.05 | bothdata$V2.y < 0.05
  bothdata = bothdata[significant_somewhere,]
}
res=paired.r(bothdata$V5.x, bothdata$V5.y, n=bothdata$V7.x, n2=bothdata$V7.y, twotailed=T)
outdata = data.frame(Markers=bothdata$Row.names,
                     Ppos=bothdata$V1.x, rpos=bothdata$V5.x, 
                     Pneg=bothdata$V1.y, rneg=bothdata$V5.y,
                     Pdiff=res$p, FDRdiff=p.adjust(res$p))
outdata = outdata[order(outdata$Pdiff),]
write.table(outdata, sep='\t', quote=F, row.names=F, file=p$out_file)
paste('rpos > rneg:', sum(abs(bothdata$V5.x) > abs(bothdata$V5.y)))
paste('rpos < rneg:', sum(abs(bothdata$V5.x) < abs(bothdata$V5.y)))