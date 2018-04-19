import re
import numpy
import math
from sklearn import preprocessing
#print (len(pos_arr))
#out=open(r"C:\Users\wbzhe\Stentor\intergenic\reads_count_all_less20_intergenic.txt","w")
pro_inf=open(r"C:\Users\wbzhe\Stentor\intergenic\PP.pro")
former_title=0
new=1
m=numpy.array([])
for line in pro_inf:
	mobj=re.match(">(NODE.*)",line)
	if mobj:
		if new==1:
			new=0
			former_title=0
			this_list=[0]*421
		else:
			if former_title!=0:
				former_title=0
				if m.any():
					m=numpy.row_stack((this_list,m))
				else:
					m=numpy.mat(this_list)
				this_list=[0]*421
			else:
				former_title=0
				this_list=[0]*421
	else:
		former_title=1
		mobj_1=re.match("(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)",line)
		pos=int(mobj_1.group(1))
		#print (pos)
		A=int(mobj_1.group(2))
		C=int(mobj_1.group(3))
		G=int(mobj_1.group(4))
		T=int(mobj_1.group(5))
		sum=A+C+G+T
		this_list[pos]+=math.log(sum)
		if not pos in range(len(this_list)):
			print (pos)
m_scaled=preprocessing.scale(m,axis=1)
min_max_scaler=preprocessing.MinMaxScaler(feature_range=(0,10))
m_minmax=min_max_scaler.fit_transform(m)
i,j=m.shape
for ni in range(i):
	for nj in range(j):
		print ("%d\t" % m_minmax[ni,nj],end="")
	print ("\n")