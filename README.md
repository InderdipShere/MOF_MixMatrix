Project Name:  MOF Polymer Mix Matrix Composite

Property Calculations-
A. Atomic density distribution
Code:  Density_profile_gro.py
Input file:  input_Density_profile.py
    Natom      = <Total number of atoms in the Trajectory files>
    Ntrack     = 0 
    Orientation= <Density distribution along which direction: 1-x, 2-y,3-z>
    Nbin       = <Nummber of bins required along the selected direction- affects the bin size>
    Ntraj      = <Total number of trajectories to be sampled>
    Box        = <[5.0, 5.0, 7.8 ], box dimension in nm>
    Ntraj_beg,Ntraj_skip = <Frames to skip in the begning>, <Interval of sampling>
    traj_head,traj_tail  = <2, Head lines>, <1, tail lines of the single frame>
    traj_file  = <Trajectory file, in gromacs formate>
    track_file = <"atom_track.dat", file of atom id>
Supporting files:  atom_track.dat
    Atom_track.dat containts <First line:  total number of tracked atoms, in the single column atom id of tracked atoms>
Trajectory file: *.gro (gromacs format)
Excecution:-
" python3 Density_profile_gro.py"

Output:  Density.dat
Two columns (distance (nm), Number density (N Å-3) )
