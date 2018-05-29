import re
from Bio import SeqIO
def addword2dict(thedict, key_a, key_b, val):
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})
def find_cloest_summit(chrs,sc_pos,types):
	#print (chr,sc_pos,type)
	np_limit=10
	average_gene_length=1000
	if types=='+':
		np=0
		for i in range(1,average_gene_length):
			p_a=str(sc_pos-i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t-%d\t%s\t-%s\t%s\t%s\t%s" % (chrs,sc_pos,np,p_a,i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a]),file=output1)
				if np>=np_limit:
					break
		np=0
		for i in range(average_gene_length):
			p_a=str(sc_pos+i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t+%d\t%s\t+%s\t%s\t%s\t%s" % (chrs,sc_pos,np,p_a,i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a]),file=output1)
				if np>=np_limit:
					break
	if types=='-':
		np=0
		for i in range(1,average_gene_length):
			p_a=str(sc_pos+i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t-%d\t%s\t-%s\t%s\t%s\t%s" % (chrs,sc_pos,np,p_a,i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a]),file=output1)
				if np>=np_limit:
					break
		np=0
		for i in range(average_gene_length):
			p_a=str(sc_pos-i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t+%d\t%s\t+%s\t%s\t%s\t%s" % (chrs,sc_pos,np,p_a,i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a]),file=output1)
				if np>=np_limit:
					break
xls_inf=open(r"C:\Users\wbzhe\Paramecium2\pcaud\pooled\N_dc2_projects_ciliate_Paramecium_caud_pcaud_mnase.smooth.positions.xls")
output1=open(r"C:\Users\wbzhe\Paramecium2\pcaud\pooled\start_codon_summit.txt","w")
smt_v_dict=dict()
smt_f_dict=dict()
for line in xls_inf:
	mobj=re.match(".*sca.*?_(\d+)\t\d+\t\d+\t(\d+)\t(\d+\.\d+)\t(\d+\.\d+)",line)
	if mobj:
		chr=mobj.group(1)
		smt=mobj.group(2)
		smt_v=mobj.group(3)
		smt_f=mobj.group(4)
		addword2dict(smt_v_dict,chr,smt,smt_v)
		addword2dict(smt_f_dict,chr,smt,smt_f)
gff_inf=open(r"C:\Users\wbzhe\Paramecium\annotations\caudatum.gff3");
name_former='NO'
for line in gff_inf:
	mobj=re.match("(.*?)\t.*?\t(.*?)\t(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\n",line)
	name_this=mobj.group(1)
	tag_this=mobj.group(2)
	start_this=mobj.group(3)
	end_this=mobj.group(4)
	frame_this=mobj.group(5)
	if (name_this == name_former):
		if (tag_this == 'CDS'):
			if (tag_former == 'CDS'):
				intergene_start=end_this
				intergene_type_f=frame_this
				tag_former=tag_this;
			else:
				mobj2=re.match(".*?_(\d+)",name_this)
				chr=mobj2.group(1)
				intergene_end=start_this
				intergene_type_l=frame_this
				if intergene_type_f == '-':
					start_codon_pos=int(intergene_start)
					find_cloest_summit(chr,start_codon_pos,intergene_type_f)
				if intergene_type_l == '+':
					start_codon_pos=int(intergene_end)
					find_cloest_summit(chr,start_codon_pos,intergene_type_l)
				intergene_start=end_this
				intergene_type_f=frame_this
				tag_former=tag_this
		else:
			tag_former=tag_this
		name_former=name_this
	else:
		if (tag_this == 'CDS'):
			name_former=name_this
			intergene_start=end_this
			intergene_type_f=frame_this
			tag_former=tag_this
gff_inf.close()

