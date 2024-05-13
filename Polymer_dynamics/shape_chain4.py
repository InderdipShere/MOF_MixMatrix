# To calculate the shape factor of polymer backbone
# Details of polymer, get the number of chain, coordinate of backbone
# calculate the calculate the shape factor: gyration tensor, eigen value, Rg, Orientation
# Update 17Oct23: To include Orientation mobility
# Update 18OCt23: To include for MSD
# Update 21Nov23:  identify the chain first and let them be same throught out
import math
import numpy as np

inp=open("input_shape_chain.dat",'r')
nchain=int(inp.readline().split()[0])
ndx_file=inp.readline().split()[0]
line=inp.readline()
natom=int(line.split()[0])
index_skip=int(line.split()[1])
line=inp.readline().split()
trj_file=line[0]
nhead=int(line[1])
ntail=int(line[2])
line=inp.readline().split()
ntraj=int(line[0])
nskip=int(line[1])
size_cut=int(inp.readline().split()[0]) #cut-off size
line=inp.readline().split()
nbin=int(line[0])
#with_zone=bool(line[1])
with_zone=line[1].strip().lower() == "true"
if(with_zone):
	zone=inp.readline().split(",")
num_avg=inp.readline().split()[0].strip().lower()=="true"
inp.close()

temp=open("Temp.ndx",'w')
trj=open(trj_file,'r')
log=open("log_shape_chain.dat",'w')
out=open("out_shape_chain.dat",'w')
ocf=open("out_raw_ocf.dat",'w')
out.writelines("#z iz, Rg^2, Shape_factor(3), ang(z,x)\n")
log.writelines(f"Total atoms {natom}, with index skip of {index_skip} \n")
log.writelines(f"Total number of chains {nchain}\n")
log.writelines(f"Trajectory file {trj_file} with head = {nhead} and tail={ntail}\n")
log.writelines(f"Total trajectory {ntraj}, with skip {nskip}\n")
log.writelines(f"Total number of bins {nbin} \n")
if (with_zone):
	log.writelines(f"z is sliced in {zone} \n")
log.writelines(f"\n")
ocf.writelines("#Traj Seg_chainID Len Principal_direction[3] ZBin# COM[3] Box[3]\n")


line_skip=(nskip-1)*(natom+nhead+ntail)
line_skip=max(line_skip,0)
max_chain_bin=nchain*nbin

c_index=[]
c_len=[]
c_srt=[]
c_srt.append(0)
ndx=open(ndx_file,'r')
for i in range(nchain):
	ndx.readline() #chainID
	aline=ndx.readline().split()
	c_len.append(len(aline))
	for j in aline:
		c_index.append(int(j)-1+index_skip)
	c_srt.append(len(c_index)) #for next chain
	ndx.readline() # blank
ndx.close()
print("ncahin", nchain, c_srt[0], c_index[c_srt[0]])
catom=len(c_index)
#rx_ext=[0.0]*natom 
#ry_ext=[0.0]*natom 
#rz_ext=[0.0]*natom
rx_ext,ry_ext,rz_ext=np.zeros(catom),np.zeros(catom),np.zeros(catom)
bin_index=np.zeros(natom)
z_count=np.zeros(nbin)
seg_index=np.zeros(catom)   ### should the segmet be []
seg_chain=0
seg_len=[]
seg_Zid=[]
seg_beg=[]
seg_beg.append(0)
rx,ry,rz=np.zeros(natom),np.zeros(natom),np.zeros(natom)

#Calculates the segments of the chain
print("First frame calcuation")
#First frame
trj.readline() #blank
natom=int(trj.readline()) #why are we reading natom
for i in range(natom):
	xyz=trj.readline()
	rx[i]=float(xyz[20:28].split()[0])
	ry[i]=float(xyz[28:36].split()[0])
	rz[i]=float(xyz[36:44].split()[0])
	bin_index[i]=-1
box=trj.readline().split()
lx,ly,lz=float(box[0]),float(box[1]),float(box[2])
lx_inv,ly_inv,lz_inv = 1/lx, 1/ly, 1/lz
lx2,ly2,lz2=lx*0.5,ly*0.5,lz*0.5
bin_size=lz/float(nbin)
bin_size_inv=1.0/bin_size
print("First trajectory read")
# Put all the particle inside the box
for i in range(natom):		
	x=rx[i]-lx2 ; x-=round(x*lx_inv)*lx ; rx[i]=x+lx2
	y=ry[i]-ly2 ; y-=round(y*ly_inv)*ly ; ry[i]=y+ly2
	z=rz[i]-lz2 ; z-=round(z*lz_inv)*lz ; rz[i]=z+lz2
print("PBC applied")
# get the index of each atoms
for iz in range(nbin):
	if(with_zone):
		zbeg=float(zone[iz])
		zend=float(zone[iz+1])
	else:
		zbeg=float(iz)*bin_size
		zend=zbeg+bin_size
	for i in range(nchain):
		ii = c_srt[i]
		nlen=c_len[i]
		for j in c_index[ii:ii+nlen]:
			z=rz[j]
			if(z>=zbeg and z<zend):
				bin_index[j]=iz
print("Z-Index of all the chain atoms assigned")
#Unwrap all the chains
for i in range(nchain):  ## what is the need to unwrap?
	ii=c_srt[i]
	nlen=c_len[i]
	kr=ii
	for j in range(nlen):
		ak=c_index[kr]
		bk=c_index[ii+j]
		xr,yr,zr=rx[ak], ry[ak], rz[ak]
		dx=rx[bk]-xr ; dx-=round(dx*lx_inv)*lx
		dy=ry[bk]-yr ; dy-=round(dy*ly_inv)*ly
		dz=rz[bk]-zr ; dz-=round(dz*lz_inv)*lz
		rx[bk]=xr+dx
		ry[bk]=yr+dy
		rz[bk]=zr+dz
		kr=ii+j
print("Chains are unwraped")
bond_dis2_max=0.4
#Identify the segments
seg_chain = 0
seg_index=[]
seg_beg=[]
seg_beg.append(0) # the first one
for i in range(nchain):
	ii=c_srt[i]
	nlen=c_len[i]
	kr=c_index[ii]
	cur_seg_count=1
	seg_index.append(ii)
	for j in c_index[ii+1:ii+nlen]:
		ak=bin_index[kr]
		bk=bin_index[j]
		if (ak==bk):
			cur_seg_count+=1
			seg_index.append(j)
		else:
			if(cur_seg_count<size_cut): #discard
				for k in range(cur_seg_count):
					seg_index.pop()
				seg_index.append(j)
				cur_seg_count=1
				print(cur_seg_count, j)
			else:
				seg_chain+=1
				seg_len.append(cur_seg_count)
				seg_beg.append(seg_beg[-1]+cur_seg_count)
				seg_index.append(j) # for the first of the next chain
				#print("seg_chain:", seg_chain, seg_beg[seg_chain-1], cur_seg_count, len(seg_index), seg_index)
			cur_seg_count=1
		#print(i, kr, j, int(ak),int(bk), cur_seg_count, seg_chain, len(seg_index))
		kr=j
	if(cur_seg_count<size_cut):
		pass #do nothing
	else:
		seg_chain+=1
		seg_len.append(cur_seg_count)
		seg_beg.append(seg_beg[-1]+cur_seg_count)
		#print("seg_chain::", seg_chain, seg_beg[seg_chain-1], cur_seg_count, len(seg_index), seg_index)
print("Identification of chain segment done")
print("Total number of chain segments =", seg_chain)
log.writelines("Seg_cahin#, seg_len, seg_beg, index \n")
for i in range(seg_chain):
	log.writelines(f"{i} {seg_len[i]} {seg_beg[i]} :\t {seg_index[seg_beg[i]: seg_beg[i]+seg_len[i]]}\n")
print("len seg_index", len(seg_index))
log.writelines("\n")
log.writelines("Traj#: Bin#, Chain#( Chain_len_counted) ZCOM Rg^2_chain eigenvalues(3) Ang_Z, Ang_X \n")

Rg=np.zeros((seg_chain,nbin))
shape_1,shape_2,shape_3=np.zeros((seg_chain,nbin)),np.zeros((seg_chain,nbin)),np.zeros((seg_chain,nbin))
ang_z,ang_x=np.zeros((seg_chain,nbin)),np.zeros((seg_chain,nbin))
for i in range(seg_chain): #first frame calcuation
	ii=seg_beg[i]
	nlenz=seg_len[i]
	iz=int(bin_index[seg_index[ii]])
	for j in range(nlenz):
		jj=int(seg_index[ii+j])
		rx_ext[j]=rx[jj]  # unwrapped
		ry_ext[j]=ry[jj]
		rz_ext[j]=rz[jj]
	nlen_inv=1.0/nlenz
	xc=sum(rx_ext[0:nlenz])*nlen_inv
	yc=sum(ry_ext[0:nlenz])*nlen_inv
	zc=sum(rz_ext[0:nlenz])*nlen_inv
	xx, yy, zz=0.0, 0.0, 0.0
	xy, yz, xz=0.0, 0.0, 0.0
	for k in range(nlenz):
		dx=rx_ext[k]-xc
		dy=ry_ext[k]-yc
		dz=rz_ext[k]-zc
		xx+=dx*dx
		yy+=dy*dy
		zz+=dz*dz
		xy+=dx*dy
		yz+=dy*dz
		xz+=dx*dz
	xx*=nlen_inv
	yy*=nlen_inv
	zz*=nlen_inv
	xy*=nlen_inv
	yz*=nlen_inv
	xz*=nlen_inv
	xc-=lx2 ;  xc-=round(xc*lx_inv)*lx ; xc+=lx2 #in box
	yc-=ly2 ;  yc-=round(yc*ly_inv)*ly ; yc+=ly2
	zc-=lz2 ;  zc-=round(zc*lz_inv)*lz ; zc+=lz2
	mat=[[xx,xy,xz], [xy,yy,yz], [xz,yz,zz]]
	eigenvalues, eigenvectors = np.linalg.eig(np.array(mat))
	idx_max = np.argmax(eigenvalues)
	eigenvalues=-np.sort(-eigenvalues)
	principal_dirc=eigenvectors[:, idx_max]
	angz=math.acos(np.dot([0,0,1], principal_dirc))*180/math.pi    
	angx=math.acos(np.dot([1,0,0], principal_dirc))*180/math.pi
	angz=min(angz,180-angz)
	angx=min(angx,180-angx)
	rg=sum(eigenvalues)
	eigenvalues=np.multiply(eigenvalues,1/rg)
	ocf.writelines(f"{1} {i} {nlenz} \t {principal_dirc[0]} {principal_dirc[1]} {principal_dirc[2]} {iz} {xc} {yc} {zc} {lx} {ly} {lz}\n")
	log.writelines(f"{1}: {iz+1} {i+1} ({nlenz}) {zc} {rg} \t [{eigenvalues}] \t {angz}, {angx}\n")
	if(num_avg):
		Rg[i][iz]     +=rg*nlenz
		shape_1[i][iz]+=eigenvalues[0]*nlenz
		shape_2[i][iz]+=eigenvalues[1]*nlenz
		shape_3[i][iz]+=eigenvalues[2]*nlenz
		ang_z[i][iz]  +=angz*nlenz
		ang_x[i][iz]  +=angx*nlenz
		z_count[iz]           +=nlenz
	else:
		Rg[i][iz]     +=rg
		shape_1[i][iz]+=eigenvalues[0]
		shape_2[i][iz]+=eigenvalues[1]
		shape_3[i][iz]+=eigenvalues[2]
		ang_z[i][iz]  +=angz
		ang_x[i][iz]  +=angx
		z_count[iz]           +=1

for i in range(line_skip): #End of the first frame
		trj.readline() # skip lines

print("First frame calculation done")


for trj_count in range(1,ntraj):
	line=trj.readline()
	line=trj.readline()
	natom=int(line) # why are we reading natom every frame
	for i in range(natom):
		xyz=trj.readline()
		rx[i]=float(xyz[20:28].split()[0])
		ry[i]=float(xyz[28:36].split()[0])
		rz[i]=float(xyz[36:44].split()[0])
	box=trj.readline().split()
	lx,ly,lz=float(box[0]),float(box[1]),float(box[2])
	lx_inv,ly_inv,lz_inv = 1/lx, 1/ly, 1/lz
	lx2,ly2,lz2=lx*0.5,ly*0.5,lz*0.5
	bin_size=lz/float(nbin)
	bin_size_inv=1.0/bin_size
	
	for i in range(natom): # Put all the particle inside the box
		x=rx[i]-lx2 ; x-=round(x*lx_inv)*lx ; rx[i]=x+lx2
		y=ry[i]-ly2 ; y-=round(y*ly_inv)*ly ; ry[i]=y+ly2
		z=rz[i]-lz2 ; z-=round(z*lz_inv)*lz ; rz[i]=z+lz2
	for seg_count in range (seg_chain):
		ii=seg_beg[seg_count]
		nlenz=seg_len[seg_count]
		nlen_inv=1.0/nlenz
		for j in range(nlenz):
			jj=int(seg_index[ii+j])
			rx_ext[j]=rx[jj]
			ry_ext[j]=ry[jj]
			rz_ext[j]=rz[jj]  #unwrap the chain
		iz=int(bin_index[int(seg_index[ii])])
		kr=0
		for k in range(1,nlenz):
			xr,yr,zr=rx_ext[kr],ry_ext[kr],rz_ext[kr]
			dx=rx_ext[k]-xr ; dx-=round(dx*lx_inv)*lx
			dy=ry_ext[k]-yr ; dy-=round(dy*ly_inv)*ly
			dz=rz_ext[k]-zr ; dz-=round(dz*lz_inv)*lz
			rx_ext[k]=xr+dx
			ry_ext[k]=yr+dy
			rz_ext[k]=zr+dz
			kr=k
		xc=sum(rx_ext[0:nlenz])*nlen_inv
		yc=sum(ry_ext[0:nlenz])*nlen_inv
		zc=sum(rz_ext[0:nlenz])*nlen_inv
		xx, yy, zz=0.0, 0.0, 0.0
		xy, yz, xz=0.0, 0.0, 0.0
		for k in range(nlenz):
			dx=rx_ext[k]-xc
			dy=ry_ext[k]-yc
			dz=rz_ext[k]-zc
			xx+=dx*dx
			yy+=dy*dy
			zz+=dz*dz
			xy+=dx*dy
			yz+=dy*dz
			xz+=dx*dz
		xx*=nlen_inv
		yy*=nlen_inv
		zz*=nlen_inv
		xy*=nlen_inv
		yz*=nlen_inv
		xz*=nlen_inv
		xc-=lx2 ;  xc-=round(xc*lx_inv)*lx ; xc+=lx2  # centre in box
		yc-=ly2 ;  yc-=round(yc*ly_inv)*ly ; yc+=ly2
		zc-=lz2 ;  zc-=round(zc*lz_inv)*lz ; zc+=lz2
		mat=[[xx,xy,xz], [xy,yy,yz], [xz,yz,zz]]
		eigenvalues, eigenvectors = np.linalg.eig(np.array(mat))
		idx_max = np.argmax(eigenvalues)
		eigenvalues=-np.sort(-eigenvalues)
		principal_dirc=eigenvectors[:, idx_max]
		angz=math.acos(np.dot([0,0,1], principal_dirc))*180/math.pi    
		angx=math.acos(np.dot([1,0,0], principal_dirc))*180/math.pi
		angz=min(angz,180-angz)
		angx=min(angx,180-angx)
		rg=sum(eigenvalues)
		eigenvalues=np.multiply(eigenvalues,1/rg)
		i=seg_count
		ocf.writelines(f"{trj_count+1}  {i}  {nlenz} \t {principal_dirc[0]} {principal_dirc[1]} {principal_dirc[2]} {iz} {xc} {yc} {zc} {lx} {ly} {lz}\n")
		log.writelines(f"{trj_count+1}: {iz+1} {i+1}({nlenz}) {zc} {rg} {eigenvalues} {angz}, {angx}\n")
		if(num_avg):
			z_count[iz]+=nlenz
			Rg[i][iz]+=rg*nlenz
			shape_1[i][iz]+=eigenvalues[0]*nlenz
			shape_2[i][iz]+=eigenvalues[1]*nlenz
			shape_3[i][iz]+=eigenvalues[2]*nlenz
			ang_z[i][iz]+=angz*nlenz
			ang_x[i][iz]+=angx*nlenz
		else:
			z_count[iz]+=1
			Rg[i][iz]+=rg
			shape_1[i][iz]+=eigenvalues[0]
			shape_2[i][iz]+=eigenvalues[1]
			shape_3[i][iz]+=eigenvalues[2]
			ang_z[i][iz]+=angz
			ang_x[i][iz]+=angx
		
	for i in range(line_skip):
		trj.readline() # skip lines

for iz in range(nbin):
	if(z_count[iz]==0):
		continue
	nor=1.0/z_count[iz]
	out.writelines(f"{(iz+0.5)*bin_size} {iz+1} {sum(Rg[:,iz])*nor} \
	{sum(shape_1[:,iz])*nor} {sum(shape_2[:,iz])*nor} {sum(shape_3[:,iz])*nor} \
	{sum(ang_z[:,iz])*nor} {sum(ang_x[:,iz])*nor}   {z_count[iz]}\n")
