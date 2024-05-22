#  ! converts gro file to cif file
# ! It is assumed that the all the atoms in the gro are fully occupied.

import re
import sys
import math

gro_file=sys.argv[1]
r = re.compile("([a-zA-Z]+)([0-9]+)")
gro=open(gro_file,"r")

head=gro.readline()
natom=int(gro.readline())

lab=[]
rx=[]
ry=[]
rz=[]

for i in range(natom):
	line=gro.readline()
	lab.append(line[11:16].split()[0])
	rx.append(float(line[21:29].split()[0]))
	ry.append(float(line[29:37].split()[0]))
	rz.append(float(line[37:45].split()[0]))
line=gro.readline().split()
box_vec=[[0]*3 for i in range(3)]

for i in range(3):
	
	box_vec[i][i]=float(line[i])
if (len(line)==9):
	box_vec[0][1] = float(line[3])
	box_vec[0][2] = float(line[4])
	box_vec[1][0] = float(line[5])
	box_vec[1][2] = float(line[6])
	box_vec[2][0] = float(line[7])
	box_vec[2][1] = float(line[8])

lxyz=[]
lxyz.append(box_vec[0][0])
lxyz.append(pow(pow(box_vec[1][1],2)+pow(box_vec[1][0],2), 0.5))
lxyz.append(pow(pow(box_vec[2][2],2)+pow(box_vec[2][0],2)+pow(box_vec[2][1],2), 0.5))
angle=[]
pi=3.141592653589793238462
angle.append(math.acos((box_vec[1][0]*box_vec[2][0]+box_vec[1][1]*box_vec[2][1])/lxyz[1]/lxyz[2])*180.0/pi)
angle.append(math.acos(box_vec[2][0]/lxyz[2])*180.0/pi)
angle.append(math.acos(box_vec[1][0]/lxyz[1])*180.0/pi)
gro.close()

# writing a cif file
cif_file=gro_file.split('.')[0]+".cif"
cif=open(cif_file,'w')
cif.writelines(head)
cif.writelines(f"_cell_length_a  {lxyz[0]*10.0}\n")
cif.writelines(f"_cell_length_b  {lxyz[1]*10.0}\n")
cif.writelines(f"_cell_length_c  {lxyz[2]*10.0}\n")
cif.writelines(f"_cell_angle_alpha {angle[0]}\n")
cif.writelines(f"_cell_angle_beta  {angle[1]}\n")
cif.writelines(f"_cell_angle_gamma {angle[2]}\n")
cif.writelines(f"_space_group_name_H-M_alt    'P 1'\n")
cif.writelines(f"space_group_IT_number      1  \n \n")
cif.writelines("loop_ \n")
cif.writelines("    _atom_site_label \n")
cif.writelines("    _atom_site_occupancy \n")
cif.writelines("    _atom_site_fract_x \n")
cif.writelines("    _atom_site_fract_y \n")
cif.writelines("    _atom_site_fract_z \n")
cif.writelines("    _atom_site_adp_type \n")
cif.writelines("    _atom_site_U_iso_or_equiv \n")
cif.writelines("    _atom_site_type_symbol \n")
print("Total atoms", len(lab))
'''
fractional coordinates (https://en.wikipedia.org/wiki/Fractional_coordinates)
X =  x/a -cot(gamma)*y/a
Y =  y/b
Z =  z/c
Simplified for the trigonal matrix
'''
for  i in range(natom):
	try:
		cif.writelines(f"{lab[i]} 1.0  \
		{rx[i]/lxyz[0]-ry[i]/math.tan(angle[2]*pi/180.0)/lxyz[0]} \
		{ry[i]/box_vec[1][1]}  {rz[i]/lxyz[2]} \
		Uiso   ? \
		{r.match(lab[i]).group(1)} 	\n" )
	except:
		cif.writelines(f"{lab[i]} 1.0  \
		{rx[i]/lxyz[0]-ry[i]/math.tan(angle[2]*pi/180.0)/lxyz[0]} \
		{ry[i]/box_vec[1][1]}  {rz[i]/lxyz[2]} \
		Uiso   ? \
		{lab[i]} 	\n" )
		

cif.close()

	
		
