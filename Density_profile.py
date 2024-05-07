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
log.write(f"Number of atoms to track = {In.Ntrack}\n")
log.write(f"Box dimesion = {In.Box}\n")
log.write(f"Direction of orientiation(x=1, y=2,z=3)= {In.Orientation}\n")
log.write(f"Number of bin = {In.Nbin}\n")
log.write(f"Name of the trajectory file = {In.traj_file}\n")
log.write(f"Start of the file, and file skip = {In.Ntraj_beg}, {In.Ntraj_skip}\n")
log.write(f"File containing list of track atom = {In.track_file}\n Atoms")


Atom_track=[]
track=open(In.track_file,'r')
for i in range(In.Ntrack):
    k=track.readline()
    j=int(k[0])
    Atom_track.append(j)
    log.write(f"{j}")
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
Ntraj_skip=One_frame*(In.Ntraj_skip)


for traj_count in range(1,In.Ntraj+1):
    atom_count=0
    rx=[]
    ry=[]
    rz=[]
    for i in range(In.traj_head):
        traj.readline()
    for i in range(1,In.Natom+1):
        traj_line=traj.readline()
        if(i==Atom_track[atom_count]) :
            atom_count = atom_count + 1
            lab,x,y,z=traj_line.split()
            rx.append(float(x))
            ry.append(float(y))
            rz.append(float(z))
            ig= int(float(z)*bin_thick_inv)
            density_profile[ig]=density_profile[ig] + 1
    for i in range(In.traj_tail):
        traj.readline
    for i in range(Ntraj_skip):
        traj.readline
        
    
    


