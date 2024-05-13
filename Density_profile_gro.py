# -*- coding: utf-8 -*-
"""
Code calcuates density profile of selected atoms for specific duration
Trajectories are in gromacs formate


"""

#import numpy as np
import input_Density_profile as In
""" file_input=open("input_Density_profile.py",'r')
    Natom
    Ntrack
    Orientation
    Box
    Nbin
    traj_file
    traj_head, traj_tail
    Ntraj 
    Ntraj_beg,  Ntraj_skip    
    track_file
"""
log=open("log_Density_profile.dat",'w')
log.write(f"Total number of atoms    = {In.Natom}\n")
#log.write(f"Number of atoms to track = {In.Ntrack}\n")
log.write(f"Box dimesion = {In.Box}\n")
log.write(f"Direction of orientiation(x=1, y=2,z=3)= {In.Orientation}\n")
log.write(f"Number of bin = {In.Nbin}\n")
log.write(f"Name of the trajectory file = {In.traj_file}\n")
log.write(f"Start of the file, and file skip = {In.Ntraj_beg}, {In.Ntraj_skip}\n")
log.write(f"File containing list of track atom = {In.track_file}\n Atoms")


Atom_track=[]
track=open(In.track_file,'r')
Ntrack=int(track.readline().split()[0]) # heading
for i in range(Ntrack):
    k=track.readline().split()
    j=int(k[0])
    Atom_track.append(j)
    log.write(f"{j}\n")
track.close()
log.write(f"Total atoms read to be tracked = {len(Atom_track)}")
    
bin_thick=In.Box[In.Orientation-1]/float(In.Nbin)
bin_thick_inv=1.0/bin_thick
density_profile=[]
for i in range(In.Nbin):
    density_profile.append(0)

traj=open(In.traj_file, 'r')
One_frame=In.Natom+In.traj_head+In.traj_tail
Nhead_skip=One_frame*(In.Ntraj_beg-1)
for i in range(Nhead_skip):
    traj.readline()
Ntraj_skip=One_frame*(In.Ntraj_skip-1)


for traj_count in range(1,In.Ntraj+1):
    atom_count=0
    rx=[]
    ry=[]
    rz=[]
    for i in range(In.traj_head):
        traj_line=traj.readline()
        #print("head",i,  traj_line)
    for i in range(1,In.Natom+1):
        traj_line=traj.readline()
        #print(traj_count, i, atom_count, Atom_track[atom_count], traj_line)
        if (atom_count <Ntrack-1):
            if(i==Atom_track[atom_count]) :
                atom_count = atom_count + 1
                x=float(traj_line[20:28].split()[0])
                y=float(traj_line[28:36].split()[0])
                z=float(traj_line[36:44].split()[0])
                rxyz=[]
                rxyz.append(x)
                rxyz.append(y)
                rxyz.append(z)
                # apply PBC
                for j in range(3):
                    rxyz[j]=rxyz[j]-round(rxyz[j]/In.Box[j])*In.Box[j]                 
                ig= int(rxyz[In.Orientation-1]*bin_thick_inv)
                density_profile[ig]=density_profile[ig] + 1
    for i in range(In.traj_tail):
        traj_line=traj.readline()
        #print("tail", i, traj_line)        
    for i in range(Ntraj_skip):
        traj.readline
bin_thick_2= bin_thick /2.0
den=open("Density.dat",'w')
for i in range(In.Nbin):
	den.write(f"{float(i)*bin_thick+bin_thick_2} \
	{float(density_profile[i])/float(In.Ntraj*Ntrack) } \n")

den.close()
