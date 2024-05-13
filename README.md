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


Polymer dynamics
Polymer behaves differently in MOF and outside MOF.
To understand its dynamics we calculate several properties like shape (shape factor), Rg, Orientation of the polymer. This is acheived by the code shape_chain4.py
Code:  shape_chain4.py
Input file: input_shape_chain.dat

    1                   # No of chains
    PVDF_BBchain.ndx    # chain index file
    16202 6120          # Natom  Index to skip before
    traj_DUT4_PVDF_nvt10.gro 2 1 # traj_file, head+tail lines
    10  1               # ntraj nskip
    5		            # cut_off size
    3 True              # nbins in z direction, To zone or not
    0.0, 4.272, 7.985, 8.622, #if zone-dimension
    True                # To use number average in calculation

Supplimentry files: index.ndx (It contains serial number of atoms in backbone chain )
    1
    34 45 201 ....
     #blank
    2
    35 44 ...
trajetory_file.gro (Trajectory file in gromacs formate)