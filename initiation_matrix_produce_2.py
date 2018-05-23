from Bio import SeqIO
import re
import datetime
import math
def addword2dict(thedict, key_a, key_b, val):
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})
start_time=datetime.datetime.now()
inf_mat=open(r"C:\Users\wbzhe\Paramecium2\PP_RNA_v2.matrix")
ok_dict=dict()
for line in inf_mat:
	mobj=re.match("(.*ig.*?)\t",line)
	if mobj:
		ok_dict.update({mobj.group(1):1})
inf_mat.close()
output=open(r"C:\Users\wbzhe\Paramecium2\PP_MNASE.matrix","w") #temporary file
pro_inf=open(r"C:\Users\wbzhe\Paramecium2\mnase_reads_mapping_ptet_51.pro") #profile input
mpile=dict()
for line in pro_inf:
	mobj=re.match(">ptetraurelia_mac_51_scaffold51_(\d+)",line)
	if mobj:
		this_tag=mobj.group(1)
	else:
		mobj_1=re.match("(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)",line)
		if mobj_1:
			pos=mobj_1.group(1)
			#print (pos)
			A=int(mobj_1.group(2))
			C=int(mobj_1.group(3))
			G=int(mobj_1.group(4))
			T=int(mobj_1.group(5))
			sums=A+C+G+T
			addword2dict(mpile,this_tag,pos,sums)
pro_inf.close()
seqs_dict=SeqIO.to_dict(SeqIO.parse(r"C:\Users\wbzhe\Paramecium\annotations\++\ptetraurelia_++.txt","fasta")) #intergenic region input
seqs_dict=sorted(seqs_dict.items(),key=lambda x:(len(str(x[1].seq))))
for key in seqs_dict:
	this_tag=key[1].id
	mobj=re.match("scaffold51_(\d+)\|(.*?)\|(\d+)\.\.(\d+)",this_tag)
	chr=mobj.group(1)
	start=mobj.group(3)
	end=mobj.group(4)
	length_this=len(str(key[1].seq))
	tag=chr + "_" + mobj.group(2) + "_" + start + "_" + end + "_" + str(length_this);
	
	
	c_start=int(end)-600 # ++
	#c_start=int(start)-100 # --
	#c_start=int(start)-100 # +-
	#c_start=int(start)-100 # -+
	if length_this<=500:
		this_list=list()
		for i in range(700):
			try:
				depth=mpile[chr][str(c_start+i)]
			except:
				depth=0
			this_list.append(depth)
		summ=sum(this_list)
		print ("%s\t" % (tag),file=output,end="")
		for i in range(len(this_list)):
			try:
				dp=math.log(this_list[i]/summ*1000+1)
			except:
				dp=0
			print ("%.3f\t" % (dp),file=output,end="")
		print ("\n" % (),file=output,end="")

inf1=open(r"C:\Users\wbzhe\Paramecium2\PP_MNASE.matrix") #temporary file
out1=open(r"C:\Users\wbzhe\Paramecium2\PP_MNASE_v2.matrix","w") #output matrix
for i in range(700):
	print ("\t%d" % (i),file=out1,end="")
print ("",file=out1)
for line in inf1:
	mobj=re.match("(.*ig.*?)\t",line)
	if mobj:
		if mobj.group(1) in ok_dict:
			print ("%s" % (line),file=out1,end="")
end_time=datetime.datetime.now()
inf1.close()
print ((end_time-start_time).seconds)