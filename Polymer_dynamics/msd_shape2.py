# Calculates msd of the chain orientation displacement
# Needs to run shape_code first
# python3 shape_chain2.py  
# python3  msd_shape.py #no_of bin used
# Edit: To include the polymer only after the first count.

import math
import numpy as np 
import sys

nbin=int(sys.argv[1])
num_avg=sys.argv[2].lower()=="true"
out=open("out_msd.dat",'w')

trj_value = np.genfromtxt('out_raw_ocf.dat', usecols=(0,), dtype=int) #just the first column to know the number of cluster
clus_traj=[]
trj_first=trj_value[0]
print(trj_value)
clus_size=0
for i in range(len(trj_value)):
	if(trj_first==trj_value[i]):
		clus_size+=1		
	else:
		trj_first=trj_value[i]
		clus_traj.append(clus_size)
		clus_size =1
clus_traj.append(clus_size) # for the last cluster
ntraj=len(clus_traj)
print("ntraj ", ntraj)
print("clusters in each traj ", clus_traj)
print("average no of cluster in the system", sum(clus_traj)/ntraj)
msdx       = [[0.0 for _ in range(nbin)] for _ in range(ntraj)]  #[[0.0]*ntraj]*nbin
msdy       = [[0.0 for _ in range(nbin)] for _ in range(ntraj)]  #[[0.0]*nbin]*ntraj
msdz       = [[0.0 for _ in range(nbin)] for _ in range(ntraj)]  #[[0.0]*nbin]*ntraj
norm_value = [[0.0 for _ in range(nbin)] for _ in range(ntraj)]  #[[0.0]*nbin]*ntraj
#("out:", len(norm_value))
#print("in:", len(norm_value[0]))

head_skip  = 1
#print("dimension", len(msdx))
for trj_count in range(ntraj):
	dat=open("out_raw_ocf.dat",'r')
	clus_id_ref=[]
	clus_size=[]
	xref=[]
	yref=[]
	zref=[]
	for i in range(head_skip):
		dat.readline()  # skipping the head lines
	nclus_ref=clus_traj[trj_count]
	for iclus in range(nclus_ref):
		line=dat.readline().split()
		clus_id_ref.append(int(line[1]))
		clus_size.append(float(line[2]))
		xref.append(float(line[7]))
		yref.append(float(line[8]))
		zref.append(float(line[9]))
	head_skip+=nclus_ref
	#print("clus_id_ref:", clus_id_ref)
	time_count=0
	for frame_count in range(trj_count+1,ntraj):
		time_count+=1
		for jclus in range(clus_traj[frame_count]):  ### index of clus_traj
			line=dat.readline().split()
			try:
				ii = clus_id_ref.index(int(line[1]))
			except:
				ii=-1
			if(ii==-1):
				continue
			iz= int(line[6])
			imass=float(line[2])
			if (imass > clus_size[ii]): # cluster decresing is okay
				print("cluster_mass increased", imass, clus_size[ii])
				continue
			lx=float(line[10]) ; lx_inv=1.0/lx
			ly=float(line[11]) ; ly_inv=1.0/ly
			lz=float(line[12]) ; lz_inv=1.0/lz
			dx=xref[ii]-float(line[7]) ; dx-=round(dx*lx_inv)*lx ### why MIC
			dy=yref[ii]-float(line[8]) ; dy-=round(dy*ly_inv)*ly ## particle jump box
			dz=zref[ii]-float(line[9]) ; dz-=round(dz*lz_inv)*lz
			if(num_avg):				
				norm_value[time_count][iz]+=imass
				msdx[time_count][iz]+=pow(dx,2)*imass
				msdy[time_count][iz]+=pow(dy,2)*imass
				msdz[time_count][iz]+=pow(dz,2)*imass
			else:
				norm_value[time_count][iz]+=1
				msdx[time_count][iz]+=pow(dx,2)
				msdy[time_count][iz]+=pow(dy,2)
				msdz[time_count][iz]+=pow(dz,2) 		
	dat.close()

print(norm_value)
for trj_count in range(1,ntraj):
	out.writelines(f"{trj_count}")
	for iz in range(nbin):
		value= norm_value[trj_count][iz]
		if(int(value)!=0):
			value= 1.0/value
		x=msdx[trj_count][iz]*value
		y=msdy[trj_count][iz]*value
		z=msdz[trj_count][iz]*value
		#out.writelines(f" {z} {norm_value[trj_count][iz]}\t ")
		out.writelines(f" {x+y+z} {x} {y} {z} \t ")
	out.writelines(f"\n ")
	
