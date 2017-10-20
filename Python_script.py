######################################
#########   1.对文件进行分割   #########
######################################
def filesplit(input,each):
	"""split file by the lines we want"""
	"""input:file need to split,each:lines you wanted in each file"""
	"""需要完整路径哦,each是每个文件多少行的参数"""

	name = input.split('/')[-1][:-4]  ## 具体的文件的后缀有多长可以修改(4)
	filepath = '/'.join(input.split('/')[:-1])
	each = int(each)
	c = 0
	outf = open("%s/%s_%s.sam" %(filepath,name,str(c/each)),'w')  ## sam这个后缀可以改哦
	for line in open(input,'r').xreadlines():
		outf.write(line)
		c += 1
		if c%each == 0:
			outf.close()
			outf = open("%s/%s_%s.sam" %(filepath,name,str(c/each)),'w')
	
	outf.close()
	os.system("mkdir %s/%s_split" %(filepath,name))                            ##### create a directory at samfile path
	os.system("mv %s/%s_*.sam %s/%s_split" %(filepath,name,filepath,name))    ##### move files to samfile_split



######################################################
#########   2.利用bed文件和bw文件进行信号的绘制   #########
######################################################

import sys
import math
import numpy
from optparse import OptionParser
from bx.bbi.bigwig_file import BigWigFile

def profile_on_bed(wig, bed, bin = 400, upstream = 2000, downstream = 2000):
	
	signal = BigWigFile(open(wig, 'rb'))
	ave_signal = {}
	for line in open(bed, 'r').xreadlines():
		line = line.strip().split('\t')
		sum_signal = signal.summarize(line[0], max(0,int(line[1]) - upstream), int(line[1]) + downstream, bin)
		if not sum_signal:
			continue
		if (int(line[1]) - upstream) < 0:
			continue
		else:
			val_signal = sum_signal.sum_data/sum_signal.valid_count
			for i in range(len(val_signal)):
				if math.isnan(val_signal[i]):
					val_signal[i] = 0
			if line[5] == '+':
				ave_signal['\t'.join(line)] = val_signal
			elif line[5] == '-':
				ave_signal['\t'.join(line)] = val_signal[::-1]
	
	ave_sum = [sum(t) for t in ave_signal.values()]
	upper = sorted(ave_sum)[int(len(ave_sum)*0.99)]
	lower = sorted(ave_sum)[int(len(ave_sum)*0.01)]	
	out_signal = []
	outf = open(bed.split('/')[-1]+'-'+wig.split('/')[-1]+'.txt','w')
	for k,v in ave_signal.iteritems():
		if sum(v) > upper or sum(v) < lower:
			continue
		else:
			out_signal.append(v)
			print >>outf, k+'\t'+'\t'.join([str(t) for t in v])
	outf.close()		

	return [str(sum(t)*1.0/len(t)) for t in zip(*out_signal)]

def multiple_bed(bedfiles, wig, name):
	
	outf = open(name,"w")
	for bed in bedfiles:
		signal = profile_on_bed(wig, bed)
		print >>outf, '%s<-c('%bed.split('/')[-1][:-4]+','.join(signal)+')'
	print >>outf, 'x<-seq(-1995,1995,10)'
	print >>outf, 'cr<-colorRampPalette(col=c("#C8524D", "#BDD791", "#447CBE", "#775A9C"), bias=1)'
	print >>outf, 'linecols<-cr(%d)' %len(bedfiles)
	print >>outf, 'plot(x,%s,type="l",xlab="Relative Distance from center (bp)",ylab="Mean ChIP-Seq signal",col=linecols[1],xaxt="s",yaxt="s",lwd=2)' %bedfiles[0].split('/')[-1][:-4] 
	for i in range(1,len(bedfiles)):
		print >>outf, 'lines(x,%s,type="l",col=linecols[%d],lwd=2)'%(bedfiles[i].split('/')[-1][:-4], i+1)
	print >>outf, 'legend("topleft",legend=c(%s),col=linecols,pch=15,bty="o",box.lty=0)' %','.join(['"'+t.split('/')[-1][:-4]+'"' for t in bedfiles])
	print >>outf, 'abline(v=0.000000,lty=2,col=c("black"))'
	outf.close()

def multiple_wig(bed, wigfiles, name):

	outf = open(name,"w")
	print >>outf, 'pdf("%s_sitepro_around_TSS.pdf")'%(wigfiles[0].split('/')[-1][4:-3])
	for wig in wigfiles:
		signal = profile_on_bed(wig, bed)
		print >>outf, '%s<-c('%wig.split('/')[-1][:-3]+','.join(signal)+')'
	print >>outf, 'x<-seq(-1995,1995,10)'
	print >>outf, 'cr<-colorRampPalette(col=c("#C8524D", "#BDD791", "#447CBE", "#775A9C"), bias=1)'
	print >>outf, 'linecols<-cr(%d)' %len(wigfiles)
	print >>outf, 'plot(x,%s,type="l",xlab="Relative Distance from center (bp)",ylab="Mean ChIP-Seq signal",col=linecols[1],xaxt="s",yaxt="s",lwd=2)' %(wigfiles[0].split('/')[-1][:-3])
	for i in range(1,len(wigfiles)):
		print >>outf, 'lines(x,%s,type="l",col=linecols[%d],lwd=2)'%(wigfiles[i].split('/')[-1][:-3], i+1)
	print >>outf, 'legend("topleft",legend=c(%s),col=linecols,pch=15,bty="o",box.lty=0)' %','.join(['"'+t.split('/')[-1][:-3]+'"' for t in wigfiles])
	print >>outf, 'abline(v=0.000000,lty=2,col=c("black"))'
	print >>outf, 'title("%s_sitepro_around_TSS.pdf")'%(wigfiles[0].split('/')[-1][4:-3])
	print >>outf, 'dev.off()'
	outf.close()

def main():
	
	usage = "python %prog <options> "
	description = """Multiple sitepro on bed/wig files."""

	optparser = OptionParser(version="1.0000",description=description,usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

	optparser.add_option("-b","--bed",dest="bed",type="str",help="Bed file path, multifiles should be splited by commma.")
	optparser.add_option("-w","--wig",dest="wig",type="str",help="Wig file path, multifiles should be splited by commma.")
	optparser.add_option("-n","--name",dest="name",type="str",help="Name of the output r file.")

	(options,args) = optparser.parse_args()
	bedfiles = options.bed.split(',')
	wigfiles = options.wig.split(',')
	if len(bedfiles) > 1 and len(wigfiles) > 1:
		print "Error! You can not input multiple bed and wig files at the same time."
		sys.exit(1)
	elif len(bedfiles) > 1:
		multiple_bed(bedfiles, wigfiles[0], options.name)
	else:
		multiple_wig(bedfiles[0], wigfiles, options.name)	

if __name__ == "__main__":
	main()




'''
利用bed文件和bw文件进行信号的绘制,修改版
'''
import sys
import math
import numpy
from optparse import OptionParser
from bx.bbi.bigwig_file import BigWigFile

def signal_on_bed(wig, bed, bin = 1, upstream = 2000, downstream = 2000):
	signal = BigWigFile(open(wig, 'rb'))
	ave_signal = {}
	for line in open(bed, 'r').xreadlines():
		line = line.strip().split('\t')
		sum_signal = signal.summarize(line[0], max(0,int(line[1])), int(line[2]), bin)
		if not sum_signal:
			continue
		if (int(line[1]) - upstream) < 0:
			continue
		else:
			val_signal = sum_signal.sum_data/sum_signal.valid_count
			for i in range(len(val_signal)):
				if math.isnan(val_signal[i]):
					val_signal[i] = 0
			ave_signal['\t'.join(line)] = str(val_signal[0])
	return(ave_signal)

def more_bws_signal_on_bed(wigs, bed, name):
	fout = open(name,'w')
	bws = wigs.strip().split(',')
	MII_K4me3 = signal_on_bed(bws[0],bed)
	MII_K27me3 = signal_on_bed(bws[1],bed)
	MII_K9me3 = signal_on_bed(bws[2],bed)
	_2cell_K4me3 = signal_on_bed(bws[3],bed)
	_2cell_K27me3 = signal_on_bed(bws[4],bed)
	_2cell_K9me3 = signal_on_bed(bws[5],bed)
	_4cell_K4me3 = signal_on_bed(bws[6],bed)
	_4cell_K27me3 = signal_on_bed(bws[7],bed)
	_4cell_K9me3 = signal_on_bed(bws[8],bed)
	ICM_K4me3 = signal_on_bed(bws[9],bed)
	ICM_K27me3 = signal_on_bed(bws[10],bed)
	ICM_K9me3 = signal_on_bed(bws[11],bed)
	print >>fout, '\t'.join(['chrom','tss','tts','transcript','gene','strand','MII_K4me3','MII_K27me3','MII_K9me3','2cell_K4me3','2cell_K27me3','2cell_K9me3','4cell_K4me3','4cell_K27me3','4cell_K9me3','ICM_K4me3','ICM_K27me3','ICM_K9me3'])
	for i in MII_K4me3.keys():
		print >>fout, i+'\t'+'\t'.join([MII_K4me3[i],MII_K27me3[i],MII_K9me3[i],_2cell_K4me3[i],_2cell_K27me3[i],_2cell_K9me3[i],_4cell_K4me3[i],_4cell_K27me3[i],_4cell_K9me3[i],ICM_K4me3[i],ICM_K27me3[i],ICM_K9me3[i]])
	fout.close()

more_bws_signal_on_bed('/mnt/Storage/home/xiuwenchao/MouseEmbryo/MII.K4me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/MII.K27me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/MII.K9me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/2cell.K4me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/2cell.K27me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/2cell.K9me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/4cell.K4me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/4cell.K27me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/4cell.K9me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/ICM.K4me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/ICM.K27me3.bw,/mnt/Storage/home/xiuwenchao/MouseEmbryo/ICM.K9me3.bw','/mnt/Storage/home/wangcf/annotations/mm9.promoter.bed','/mnt/Storage/home/xiuwenchao/MouseEmbryo/promoter_histone_signal.txt')



################################################
#########   3.获取nk的promoter bed文件   #########
################################################

def get_promoter_nk(promoter2k,promoternk,length):
	'''
	get nk promoter annotaion base on 2k promoter annotation
	'''
	fout = open(promoternk,'w')
	for i in open(promoter2k,'r').xreadlines():
		linei = i.strip().split('\t')
		print >>fout, '\t'.join([linei[0],str(int(linei[1])+2000-int(length)),str(int(linei[2])-2000+int(length)),linei[3],linei[4],linei[5]])

get_promoter_nk('/homea2/jxfeng/xiuwenchao/TSC_TET1_TET2/annotaion/mm9.promoter.bed','/homea2/jxfeng/xiuwenchao/TSC_TET1_TET2/annotaion/mm9.promoter3k.bed','3000')




##########################################
#########   4.利用python进行绘图   #########
##########################################
'利用python进行散点图绘制'
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0.0, 2.0, 0.01)
s = np.sin(2*np.pi*t)
plt.plot(t, s)

plt.xlabel('time (s)')
plt.ylabel('voltage (mV)')
plt.title('About as simple as it gets, folks')
plt.grid(True)
plt.savefig("test.png")
plt.show()

'利用python进行绘制density plot'
import seaborn as sns
sns.set_style('whitegrid')
sns.kdeplot(np.log2(np.array(regulate_M())))
sns.plt.show()
plt.xlabel("gene A expression(log2)")
plt.ylabel("probability")
plt.title("Regulate gene A")
plt.savefig("Regulate_gene_A.pdf")


''' python绘制二维通过颜色界定深浅的山地图 '''
library('ggplot2')
'绘制散点图,有深度那种'
setwd("~/Desktop/sc_RNA-seq/hESC_differentiation(0_12_24_36_72_96h)en_derived/data")
load("./Days.Rdata")
gene1 <- 'HAND1'
gene2 <- 'GATA6'

pdf('HAND1_GATA6.pdf')
for (i in 1:6){
par(mfrow=c(2,1))
plot(density(as.numeric(Days[[i]]['HAND1',])),ylim=c(0,1.2),main='x',xlim=c(-1,12))
plot(density(as.numeric(Days[[i]]['GATA6',])),ylim=c(0,1.2),main='y',xlim=c(-1,12))
x <- sample(as.numeric(Days[[i]]['HAND1',]),10000, replace = TRUE)
y <- sample(as.numeric(Days[[i]]['GATA6',]),10000, replace = TRUE)
df <- data.frame(x = x, y = y,
  d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
print(p)
}
dev.off()





##########################################################
#########   5.利用two big reader进行基因序列的读取   #########
##########################################################
'two big reader的使用方法'
import twobitreader
sequence  = "/mnt/Storage/home/xiuwenchao/annotation/mm9.2bit"
genome = twobitreader.TwoBitFile(sequence) 
genome[chrom][start:end]      # start is 0


#########################################################
#########   6.对基因区域分类,区分高/中/低甲基化区域   #########
#########################################################

import twobitreader
import re
def CGI_split(sequence ,CGI ,HCP ,ICP ,LCP,r1 = 0.75, r2 = 0.5):
	fHCP = open(HCP,'w')
	fICP = open(ICP,'w')
	fLCP = open(LCP,'w')
	genome = twobitreader.TwoBitFile(sequence) 
	for i in open(CGI,'r').xreadlines():
		linei = i.strip().split('\t')
		seq = genome[linei[0]][int(linei[1]),int(linei[2])]
		nC = len(re.findall('C',seq,re.IGNORECASE))
		nG = len(re.findall('G',seq,re.IGNORECASE))
		nCG = len(re.findall('CG',seq,re.IGNORECASE))
		seq_len = int(linei[2])-int(linei[1])
		ratio = (nCG*seq_len*1.0)/(nC*nG)
		if ratio >= r1:
			print >>fHCP, '\t'.join(linei)
		elif ratio <= r2:
			print >>fLCP, '\t'.join(linei)
		else:
			print >>fICP, '\t'.join(linei)
	fHCP.close()
	fLCP.close()
	fICP.close()


CGI_split('/mnt/Storage/home/xiuwenchao/annotation/mm9.2bit','/mnt/Storage/home/wangcf/annotations/mm9.CGI.bed'
			,'/mnt/Storage/home/xiuwenchao/mm9.HCP.bed','/mnt/Storage/home/xiuwenchao/mm9.ICP.bed','/mnt/Storage/home/xiuwenchao/mm9.LCP.bed')





#######################################################
#########   7.refseq ID 转换成 genesymbol ID   #########
#######################################################

'refseq 装换成 genesymbol'
def refseq2genesymbol(ref,infile,outfile):
	fout = open(outfile,'w')
	refseq = []
	genesymbol = []
	for i in open(ref,'r').xreadlines():
		linei = i.strip().split('\t')
		refseq.append(linei[0])
		genesymbol.append(linei[1])
	for i in open(infile,'r').xreadlines():
		if i.startswith('gene') or i.startswith('log'):#####文件开头名称,忽略的行
			print >>fout, i.strip()
			continue
		else:
			linei = i.strip().split('\t')
			for j in range(len(refseq)):
				if linei[0].upper() == refseq[j].upper():###linei[0]对应的是refseq的列，可以修改index对不同文件应用
					print >>fout, '\t'.join([genesymbol[j]]+linei[1:])
	fout.close()

				
refseq2genesymbol('refseq_genesymbol.txt','count/all_WT_Tet2_hetero_.txt','count/WT_Tet2_hetero.txt')
refseq2genesymbol('refseq_genesymbol.txt','count/all_WT_Tet2_homo_.txt','count/WT_Tet2_homo.txt')



####################################################
#########   8.利用Python进行常微分方程的求解   #########
#####################################################

'求解常微分方程'
from scipy.integrate import odeint 
import numpy as np 
import matplotlib.pyplot as plt
import time

#np.random.seed(1)
start1 = time.clock()
def hill_p(x1,x2,x3,x4):
	'''
	建立希尔方程,在正调节(λ为正)情况下
	'''
	return (x4+(1.0-x4)/(1.0+pow((x1/x2),x3)))/x4

def hill_m(x1,x2,x3,x4):
	'''
	建立希尔方程,在正调节(λ为负)情况下
	'''
	x4 = 1.0/x4
	return x4+(1.0-x4)/(1.0+pow((x1/x2),x3))

def hill(x1,x2,x3,x4):
	'''
	建立希尔方程
	'''
	return x4+(1.0-x4)/(1.0+pow((x1/x2),x3))

def regulate(w, t, g, k, treshold, hill_n, Lambda): 
    '''
    对于受其他基因调控的基因，
    进行微分方程的建立(a分别受b和c调节,b抑制a表达,c促进a表达)
    '''
    a,b,c = w
    ga,gb,gc = g
    ka,kb,kc = k
    ba0,ca0 = treshold
    nba,nca = hill_n
    lambda_ba,lambda_ca = Lambda
    return np.array([ga*hill_m(b,ba0,nba,lambda_ba)*hill_p(c,ca0,nca,lambda_ca)-ka*a,gb-kb*b,gc-kc*c])

def regulate_M(w=(1.0,1.0,1.0), t=np.arange(0, 300, 0.01), G=(1,100), K=(0.1,1), Treshold=(0.02,1.98), Hill_n=(1,6), LAMBDA=(1,100), nparamiter=1000):
	'''
	对于受其他基因调控的基因进行expression的分布进行微分方程求解,
	得到不同随机参数下的基因的表达值list,并求解表达值的中位数
	'''
	expression = []
	M_b,M_c = (isolate_M(),isolate_M())
	for i in range(nparamiter):
		g = np.random.uniform(G[0],G[1],size=3)
		k = np.random.uniform(K[0],K[1],size=3)
		treshold = np.random.uniform(Treshold[0],Treshold[1],size=2)*np.array([M_b,M_c])
		hill_n = np.random.randint(Hill_n[0],Hill_n[1]+1,size=2)
		Lambda = np.random.uniform(LAMBDA[0],LAMBDA[1],size=2)
		track = odeint(regulate, w, t, args=(g, k, treshold, hill_n, Lambda))
		expression.append(round(float(track[:,0][-1]),3))
	return np.median(expression)
#	return expression



##########################################################
#########   9.将所有fpkm文件进行合并，合成一个table   #########
##########################################################

import os
import commands
import numpy as np
fout = open('merge_fpkm.txt','w')
files = commands.getoutput('ls *.gene.fpkm').strip().split('\n')

values = []
for i in files:
	name = []
	value = []
	for j in open(i,'r').xreadlines():
		linej = j.strip().split('\t')
		name.append(linej[0])
		value.append(linej[1])
	values.append(value)

for k in range(len(name)):
	print >>fout, '\t'.join([name[k]]+(np.array(values)[:,k]).tolist())

fout.close()	







############################################################################
#########   10.对GSEA annotation中的数据进行展开,内部用0,1表示是否存在   #########
############################################################################

def unfold(fold_file,unfold_file):
	fout = open(unfold_file,'w')
	genelist = []
	genedic = {}
	termlist = []
	n=0
	m=-1
	for i in open(fold_file,'r').xreadlines():
		linei = i.strip().split('\t')
		genelist += linei[2:]
		n += 1
	genelist = sorted(list(set(genelist)))
	for j in genelist:
		genedic[j]=[0]*n
	for l in open(fold_file,'r').xreadlines():
		linel = l.strip().split('\t')
		termlist.append(linel[0])
		m += 1
		for p in linel[2:]:
			genedic[p][m]+=1
	print >>fout, 'gene_symol'+'\t'+'\t'.join(termlist)
	for q in sorted(genedic.keys()):
		print >>fout, q+'\t'+'\t'.join([str(o) for o in genedic[q]])
	fout.close()
		

foldfile = '/Users/xiu/Desktop/c5.bp.v5.2.symbols.gmt'
unfoldfile = '/Users/xiu/Desktop/c5.bp.v5.2.symbols.unfold.txt'

unfold(foldfile,unfoldfile)







