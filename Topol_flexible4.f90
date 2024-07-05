!PREPARE TOPOLOGY FOR THE flexible MOF 
!! We are going to add the more MOF 
! inputs: input_MAKE_TOPOL.dat (21)
! INPUT:  PARAMETER FILE       (22)
! input:  position of fluid    (23)
! input:  position+bonding MOF (24)
! INPUT:  MOF_CHARGE           (25)
! INPUT:  FLUID_COORDINATE     (26)

!OUTPUT:  NEW_TOPOL.top        (31)
!OUTPUT:  NEW_SYSTEM.gro       (32)
!OUTPUT:  TOPOL_DIHEDRAL.dat   (33)
!!! Charge arrangement
! C1C2H1 => (3) : C1C2,C1C2, C1C2, H1,H1,H1  allocates charges like that (1)
! C1C2H1 => (3) : C1C1C1,C2C2C2, H1H1H1  allocates charges like that (2)
! Charge arrangement can be modified.
!Edit:  25July2023-
! Added Angle type 3 (harmonic)
! ADDED FILE  34,35,36 FOR New_topol2.top generation

PROGRAM MAKE_TOPOL
  IMPLICIT NONE
  REAL*8, ALLOCATABLE, DIMENSION(:)::MASS_A, SIG_A, EPS_A, RX,RY,RZ,CHARG_MOF,&
           R0_BOND, KB_BOND, &
           T0_ANG,  KA_ANG , TX_ANG, &
           T0_DHI,  KA_DHI , TX_DHI, &
           T0_IMP,  KA_IMP , TX_IMP

  REAL*8:: FUDGELJ, FUDGEQQ, REP_ORD, B0, KB, X, Y, Z, CHARG,&
           AX,BY,CZ,BX,CX,CY, DIS_VECTOR(3)
  REAL*8, PARAMETER:: GAS_CONST=8.31446261815324E-3 , CAL2JOUL=4.1858, DEG2RAD=ASIN(1.0D0)/90.0D0

  INTEGER, ALLOCATABLE:: NMOL(:), BOND_TYPES(:,:), ATOM_1MOL(:),&
            IF_THIS(:), THEN_THIS(:), ATOM2ATOMT(:),&
            TYPE_ANG(:), TYPE_DHI(:), TYPE_IMP(:),&
            ATOM_IN_MOF(:), MOF_ID(:)
  INTEGER:: I,J,K, NBFUNC, COMB_RULE, NTYPE_A, AI,AJ,AK,AL, NTYPE_MOL, &
            NREXC, FUNCT_PARA, NBOND, ATOM_COUNT, MOL_ATOM,&
            NCONSTRAIN, HEAD_MOF_SC_CAR, HEAD_MOF_UC_CHARG, HEAD_FLUID_PARA, HEAD_FLUID_COORD,&
            NBOND_TYPES_MAX, NATOM_MOF, ATOM_COUNT2, NCONDITION, NBONDT,&
            NANGLET, NDIHEDRALT, NIMPROPERT,&
            NANGLE , NDIHEDRAL , NIMPROPER,&
            NPAIRT, NLAMMPS_ATOMT, MOL_COUNT, &
            NTYPE_MOF, NTYPE_FLUID, MOF_COUNT, OVERALL_COUNT, FLUID_COUNT, ATOMT_IN_MOF ,&
            REPEAT_STYLE
  CHARACTER*50:: PARA_FILE, FLUID_MOLECULE_CHARGE, MOF_UNIT_CELL_CHARG_FILE, MOF_SUPERCELL_CAR, &
            GROFORMAT, FLUID_COORDINATES
  CHARACTER*8:: GEN_PAIR, LAB, ATLAB, lab1,lab2,lab3
  CHARACTER*8, ALLOCATABLE, DIMENSION(:):: LAB_A, LAB_AT, MOL_LAB
  LOGICAL, ALLOCATABLE::CONNECT(:,:)

  OPEN(21,FILE="input_MAKE_TOPOL.dat")
  READ(21,*) PARA_FILE 
  GROFORMAT  = '(i5,2a5,i5,3f8.3,3f8.4)'
  

  OPEN(22,FILE=TRIM(PARA_FILE), STATUS='OLD')
  OPEN(31,FILE="New_topol.top")
  OPEN(34,FILE="New_topol2.top")
  OPEN(32,FILE="NEW_SYSTEM.gro")
  WRITE(31,*) "; TOPOLOGY FILE GENERATED FROM MAKE_TOPOL PROGRAM @INDERDIP"
  WRITE(31,*) "; INPUT PARMETERS FILE : ", TRIM(PARA_FILE)
  WRITE(34,*) "; TOPOLOGY FILE GENERATED FROM MAKE_TOPOL PROGRAM @INDERDIP"
  WRITE(34,*) "; INPUT PARMETERS FILE : ", TRIM(PARA_FILE)
  READ(22,*) ! DEFAULTS PARAMENTERS
  READ(22,*) NBFUNC
  READ(22,*) COMB_RULE
  READ(22,*) GEN_PAIR
  READ(22,*) FUDGELJ
  READ(22,*) FUDGEQQ
  READ(22,*) REP_ORD
  WRITE(31,*) "[ defaults ]"
  WRITE(31,*) "; NBFUNC, COMB-RULE, GEN_PAIR, FUDGELJ, FUDGEQQ, REPULSE_N"
  WRITE(31,'(2I3,1x, A4, 3F9.5)') NBFUNC, COMB_RULE, TRIM(GEN_PAIR), FUDGELJ, FUDGEQQ, REP_ORD
  WRITE(34,*) "[ defaults ]"
  WRITE(34,*) "; NBFUNC, COMB-RULE, GEN_PAIR, FUDGELJ, FUDGEQQ, REPULSE_N"
  WRITE(34,'(2I3,1x, A4, 3F9.5)') NBFUNC, COMB_RULE, TRIM(GEN_PAIR), FUDGELJ, FUDGEQQ, REP_ORD

  READ(22,*)! ATOMTYPES
  READ(22,*) NTYPE_A  
  ALLOCATE(LAB_A(NTYPE_A),LAB_AT(NTYPE_A), MASS_A(NTYPE_A), &
  	      SIG_A(NTYPE_A), EPS_A(NTYPE_A), ATOM_1MOL(NTYPE_A), MOF_ID(NTYPE_A))  
  PRINT*, "TOTAL TYPE OF ATOMS IN SYSTEM ", NTYPE_A
  DO I = 1, NTYPE_A
    READ(22,*) LAB_AT(I), LAB_A(I), MASS_A(I), SIG_A(I), EPS_A(I), ATOM_1MOL(I), MOF_ID(I)
  END DO 
  PRINT*, "ATOMS INFORMATION READ "
  SIG_A = SIG_A /10.0D0
  EPS_A = EPS_A *GAS_CONST ! K * kJ/mol/K
  WRITE(31,*) "[ atomtypes ]"
  WRITE(31,*) "; ATOM_LAB, ELEMENT, MASS, CHARG, PTYPE, SIG(nm), EPS(kJ/mol)"
  WRITE(34,*) "[ atomtypes ]"
  WRITE(34,*) "; ATOM_LAB, ELEMENT, MASS, CHARG, PTYPE, SIG(nm), EPS(kJ/mol)"
  DO I = 1, NTYPE_A
    WRITE(31,'(2A8, 2F10.5,A2,1x, 2F12.7)') &
    LAB_AT(I),LAB_A(I), MASS_A(I), 0.0D0, "A", SIG_A(I), EPS_A(I)
    WRITE(34,'(2A8, 2F10.5,A2,1x, 2F12.7)') &
    LAB_AT(I),LAB_A(I), MASS_A(I), 0.0D0, "A", SIG_A(I), EPS_A(I)
  END DO 
  PRINT*, "ATOM INFORMATION IS WRITTEN"
  
  ! DEFINING MOLECULES
  ATOM_COUNT    = 0
  OVERALL_COUNT = 0
  ATOMT_IN_MOF  = 0
  READ(21,*) NTYPE_MOF, NTYPE_FLUID 
  READ(22,*) ! blank
  NTYPE_MOL =NTYPE_MOF+NTYPE_FLUID
  ALLOCATE(MOL_LAB(NTYPE_MOL), NMOL(NTYPE_MOL), ATOM_IN_MOF(NTYPE_MOF))
  PRINT*, "TOTAL TYPES OF MOLECULES IN THE SYSTEM ", NTYPE_MOL
  CLOSE(34) ! WILL BE UPDATED DURING THE RUN
  
  DO MOF_COUNT= 1, NTYPE_MOF
	  OPEN(35,FILE="TOP_TEMP.top")
	  OPEN(36,FILE="DIHEDRAL_TEMP.dat")
	  READ(21,*) MOL_LAB(MOF_COUNT), NMOL(MOF_COUNT),  ATOM_IN_MOF(MOF_COUNT)!MOL_ATOM     
	  READ(21,*) MOF_UNIT_CELL_CHARG_FILE, HEAD_MOF_UC_CHARG, REPEAT_STYLE
	  READ(21,*) MOF_SUPERCELL_CAR       , HEAD_MOF_SC_CAR
	  READ(22,*) MOL_LAB(MOF_COUNT), NREXC
	  PRINT*, "LABLE OF MOF ", MOL_LAB(MOF_COUNT)
	  PRINT*, "ATOMS PER MOF ", ATOM_IN_MOF(MOF_COUNT) !MOL_ATOM
	  PRINT*, "TOTAL NUMBER OF THE MOLECULE OF MOF ", NMOL(MOF_COUNT)
	  PRINT*, "REPEAT STYLE", REPEAT_STYLE
	  PRINT*, "1= C1C2C1C2H1H1"
	  PRINT*, "2= C1C1C2C2H1H1"
	  WRITE(31,*) ! 
	  WRITE(31,*) ! MOLECULE DEFINATION
	  WRITE(31,*) "[ moleculetype ]"
	  WRITE(31,*) "; MOLENAME, NEAR EXCLUSION"
	  WRITE(31,*) TRIM(MOL_LAB(MOF_COUNT)), NREXC 

	  WRITE(31,*) ! ATOMS
	  WRITE(31,*) "[ atoms ]"
	  WRITE(31,*) "; nr, at_type, res_nr, res_name, at_name, cgnr, charge MASS"
	  WRITE(35,*) ! 
	  WRITE(35,*) ! MOLECULE DEFINATION
	  WRITE(35,*) "[ moleculetype ]"
	  WRITE(35,*) "; MOLENAME, NEAR EXCLUSION"
	  WRITE(35,*) TRIM(MOL_LAB(MOF_COUNT)), NREXC 

	  WRITE(35,*) ! ATOMS
	  WRITE(35,*) "[ atoms ]"
	  WRITE(35,*) "; nr, at_type, res_nr, res_name, at_name, cgnr, charge MASS"
  
    !charges 
    NATOM_MOF = NMOL(MOF_COUNT)*ATOM_IN_MOF(MOF_COUNT)  !MOL_ATOM
    PRINT*, "TOTAL ATOMS FOR MOF ", MOL_LAB(MOF_COUNT), NATOM_MOF
	  IF (ALLOCATED(RX)) THEN
      DEALLOCATE(RX,RY,RZ,CHARG_MOF,ATOM2ATOMT)
    ENDIF
    ALLOCATE(RX(NATOM_MOF), RY(NATOM_MOF),RZ(NATOM_MOF),CHARG_MOF(NATOM_MOF),&
		       ATOM2ATOMT(NATOM_MOF))
    ATOM2ATOMT = 0
    OPEN(25, FILE=TRIM(MOF_UNIT_CELL_CHARG_FILE),STATUS='OLD')
    DO I = 1, HEAD_MOF_UC_CHARG
      READ(25,*) !
    END DO 
	  ATOM_COUNT = 0 
    IF (REPEAT_STYLE == 1 ) THEN 
      DO K = 1, NTYPE_A
        IF (MOF_ID(K)==MOF_COUNT) THEN 
          DO I = 1, ATOM_1MOL(K)
            READ(25,*) LAB, X, Y, Z, CHARG
            DO J = 1, NMOL(MOF_COUNT)    	
              CHARG_MOF(ATOM_COUNT+(J-1)*ATOM_1MOL(K)+I ) = CHARG
            END DO
          END DO 
          ATOM_COUNT = ATOM_COUNT+ATOM_1MOL(K)*NMOL(MOF_COUNT)      
        END IF 
      END DO
    ELSEIF (REPEAT_STYLE ==2) THEN  !C1C1C2C2H1H1
      DO K = 1, NTYPE_A
        IF (MOF_ID(K)==MOF_COUNT) THEN 
          DO I = 1, ATOM_1MOL(K)
            READ(25,*) LAB, X, Y, Z, CHARG
            DO J = 1, NMOL(MOF_COUNT)     
              CHARG_MOF(ATOM_COUNT+J) = CHARG
            END DO
            ATOM_COUNT = ATOM_COUNT+NMOL(MOF_COUNT) 
          END DO     
        END IF 
      END DO
    ELSE
      PRINT*, "REPEAT STYLE IS NOT CORRECT", REPEAT_STYLE
      STOP 
    END IF 
    PRINT*, "CHARGE INFORMATION READ FOR MOF ", MOF_COUNT
    PRINT*, "TOTAL ATOMS IN MOF", ATOM_COUNT
    CLOSE(25)
    ! position lables and bonding 
    READ(21,*) NCONDITION
    IF(ALLOCATED(IF_THIS)) THEN
      DEALLOCATE(IF_THIS, THEN_THIS)
    END IF 
    ALLOCATE(IF_THIS(NCONDITION),THEN_THIS(NCONDITION))
    DO I = 1, NCONDITION
      READ(21,*) IF_THIS(I), THEN_THIS(I)
    END DO 
    open(37,file="temp.sh")
      write(37,*) 'sed -i " " "s/\//_/g" ', trim(MOF_SUPERCELL_CAR)
      call system("chmod +x temp.sh ;  ./temp.sh ")
    close(37)
    OPEN(24, FILE=TRIM(MOF_SUPERCELL_CAR),STATUS='OLD')  
    DO I = 1, HEAD_MOF_SC_CAR
      READ(24,*) !HEADERS
    END DO  
    READ(24,*) J 
    READ(24,*) NBOND
    READ(24,*) NANGLE 
    READ(24,*) NDIHEDRAL
    READ(24,*) NIMPROPER
    READ(24,*) !
    READ(24,*) NLAMMPS_ATOMT 
    READ(24,*) NBONDT
    READ(24,*) NANGLET
    READ(24,*) NDIHEDRALT
    READ(24,*) NIMPROPERT
    PRINT*, "ATOMS",    J,          NLAMMPS_ATOMT
    PRINT*, "BONDS",    NBOND,      NBONDT
    print*, "ANGLES",   NANGLE,     NANGLET
    print*, "DIHEDRAL", NDIHEDRAL,  NDIHEDRALT
    print*, "IMPROPER", NIMPROPER,  NIMPROPERT
    READ(24,*)  X, AX 
    READ(24,*)  X, BY 
    READ(24,*)  X, CZ
    READ(24,*)  BX, CX, CY   
    AX = AX/10.0; BY = BY/10.0; CZ = CZ/10.0
    BX =BX/10.0;  CX = CX/10.0; CY = CY/10.0
    READ(24,*) !BOX DIMENSION
    READ(24,*) ! MASSES
    READ(24,*)
    DO I = 1, NLAMMPS_ATOMT 
  	  READ(24,*) K, MASS_A(I)
    END DO  
    READ(24,*)!
    IF(ALLOCATED(R0_BOND)) THEN
      DEALLOCATE(R0_BOND,KB_BOND)
    END IF
    ALLOCATE(R0_BOND(NBONDT), KB_BOND(NBONDT))    
    READ(24,*) !BOND PARAMETER
    READ(24,*) !BLANK
    DO I = 1, NBONDT
      READ(24,*) J, KB_BOND(I), R0_BOND(I)
      PRINT*, I, "BOND HARMONIC", KB_BOND(I), R0_BOND(I)
    END DO
    R0_BOND =  R0_BOND /10.0  ! nm
    KB_BOND =  KB_BOND*200.0*CAL2JOUL  ! kJ/MOL/nm2  ** CONFIRM FOR KCAL / KJ
    READ(24,*) ! BLANK
    READ(22,*) ! BOND TYPE
    PRINT*, "BOND PARAMETER READ OF MOF", MOF_COUNT
     
    IF(ALLOCATED(T0_ANG)) THEN 
      DEALLOCATE(T0_ANG,KA_ANG, TX_ANG, TYPE_ANG)
    END IF 
    ALLOCATE(T0_ANG(NANGLET), KA_ANG(NANGLET), TX_ANG(NANGLET), TYPE_ANG(NANGLET))
    READ(24,*) !ANGLE PARAMETER
    READ(24,*) !BLANK  
    READ(22,*) ! ANGLES 
    DO I = 1, NANGLET
  	  READ(22,*) J, K 
      TYPE_ANG(I) = K 
  	  IF (K==1) THEN  ! COSIN/PERIODIC (C, B, n)
        READ(24,*) J, LAB1, X, J, K 
        PRINT*, I,"COSINE/PERIODIC ", X, J, K  
        Z = 0.0      
        IF (-J*(-1)**K<0) Z = 180.0 
        T0_ANG(I) = 360.0/DBLE(K)
        KA_ANG(I) = X/2025.0*CAL2JOUL/DEG2RAD/DEG2RAD 
        TX_ANG(I) = K  
        PRINT*, "GRO>ANGLE(1)", T0_ANG(I), KA_ANG(I), INT(TX_ANG(I))
      ELSEIF(K==2) THEN ! FOURIER (K, C0, C1, C2)
    	  READ(24,*) J, LAB, KB, X, Y, Z 
    	  PRINT*, I, "FOURIER", KB, X,Y,Z
    	  T0_ANG(I) = ACOS(-Y/4.0/Z)/DEG2RAD
    	  KA_ANG(I) = 2.0*KB*(X-Y+Z)/(Y/4.0/Z-1.0)**2 * CAL2JOUL
    	  TX_ANG(I) = 0.0
    	  PRINT*, "GRO>ANGLE(2)", T0_ANG(I),KA_ANG(I)
      ELSEIF(K==3) THEN ! HARMONIC (K, C0)
    	  READ(24,*) J, KB, X 
    	  PRINT*, I, "HARMONIC", KB, X
    	  T0_ANG(I) = X
    	  KA_ANG(I) = 2.0*KB* CAL2JOUL
    	  PRINT*, "GRO>ANGLE(3)", T0_ANG(I),KA_ANG(I)
      END IF 
    END DO 
    READ(24,*) ! BLANK
    PRINT*, "ANGLES PARAMETERS READ FOR MOF", MOF_COUNT

    IF (ALLOCATED(T0_DHI)) THEN 
      DEALLOCATE(T0_DHI, KA_DHI, TX_DHI, TYPE_DHI)
    END IF 
    ALLOCATE(T0_DHI(NDIHEDRALT), KA_DHI(NDIHEDRALT),&
             TX_DHI(NDIHEDRALT), TYPE_DHI(NDIHEDRALT))
    READ(24,*) ! DIHEDRAL
    READ(24,*) ! BLANK
    READ(22,*) ! DIHEDRAL
    DO I = 1, NDIHEDRALT
  	  READ(22,*) J, K 
      TYPE_DHI(I) = K 
  	  IF (K==1) THEN  ! COSIN/PERIODIC (C, B, n)
        READ(24,*) J, X, J, K 
        PRINT*, I, "COSIN/PERIODIC", X, J, K 
        Z = 0.0      
        IF (-J*(-1)**K<0) Z = 180.0 
        T0_DHI(I) = Z 
        KA_DHI(I) = 2.0/DBLE(K*K)*X*CAL2JOUL 
        TX_DHI(I) = K 
        PRINT*, "GRO>Dihedral(9)", T0_DHI(I),KA_DHI(I),INT(TX_DHI(I))
      ELSEIF(K==2) THEN ! FOURIER (K, C0, C1, C2)  !!!! MULTIPLE DIHEDRAL
    	  READ(24,*) J, Z, KB, X, Y 
    	  PRINT*, I, 'FOURIER', KB, X, Y 
    	  T0_DHI(I) = Y 
    	  KA_DHI(I) = KB*CAL2JOUL
    	  TX_DHI(I) = X
    	  PRINT*, "GRO>Dihedral(9)", T0_DHI(I),KA_DHI(I),TX_DHI(I), 0.0,0.0,0.0
      ELSEIF(K==3) THEN ! HARMONIC (K, d, n)
    	  READ(24,*) J, KB, J, K 
    	  PRINT*, I, 'HARMONIC', KB, J , K
    	  Z = 0.0
    	  IF(J<0) Z = 180.0
    	  T0_DHI(I) = Z 
    	  KA_DHI(I) = KB *CAL2JOUL
    	  TX_DHI(I) = K 
    	  PRINT*, "GRO>Dihedral(9)", T0_DHI(I), KA_DHI(I), INT(TX_DHI(I))
      END IF 
    END DO 
    READ(24,*) ! BLANK
    PRINT*, "DIHEDRAL PARAMETERS READ FOR MOF", MOF_COUNT
     
    IF(ALLOCATED(T0_IMP)) THEN
      DEALLOCATE(T0_IMP, KA_IMP, TX_IMP, TYPE_IMP)
    END IF 
    ALLOCATE(T0_IMP(NIMPROPERT), KA_IMP(NIMPROPERT),&
             TX_IMP(NIMPROPERT), TYPE_IMP(NIMPROPERT))
    READ(24,*) ! IMPROPER
    READ(24,*) ! BLANK
    READ(22,*) ! IMPROPER
    DO I = 1, NIMPROPERT
  	  READ(22,*) J, K 
      TYPE_IMP(I) = K 
  	  IF (K==1) THEN  ! COSIN/PERIODIC (C, B, n)
        READ(24,*) J, X, J, K 
        PRINT*, I, "COSINE/PERIODIC", X, J, K 
        Z = 0.0      
        IF (-J*(-1)**K<0) Z = 180.0 
        T0_IMP(I) = Z 
        KA_IMP(I) = 2.0/DBLE(K*K)*X*CAL2JOUL 
        TX_IMP(I) = K
        PRINT*, "GRO>Dihedral(9)", T0_IMP(I), KA_IMP(I), INT(TX_IMP(I))
      ELSEIF(K==2) THEN ! FOURIER (K, C0, C1, C2)
    	  READ(24,*) J, KB, X, Y, Z   
    	  PRINT*, I, "FOURIER", KB, X, Y, Z 
    	  T0_IMP(I) = KB*(X - Z)*CAL2JOUL 
    	  KA_IMP(I) =-KB*Y      *CAL2JOUL
    	  TX_IMP(I) = 2.0*Z*KB  *CAL2JOUL
    	  PRINT*, "GRO>Dihedral(3)", T0_IMP(I),KA_IMP(I),TX_IMP(I), 0.0,0.0,0.0
      ELSEIF(K==4) THEN ! IMPROPER FOURIER (K, C0, C1, C2) 
    	  READ(24,*) J, KB, X, Y, Z 
    	  PRINT*, I, "IMPROPER_FOURIER", KB, X, Y, Z     	
    	  T0_IMP(I) = ACOS(-Y/4.0/Z)/DEG2RAD    	
    	  KA_IMP(I) =KB*(X-Y+Z)/(Y/4.0/Z-1.0)**2  *CAL2JOUL
    	  TX_IMP(I) = 0.0
    	  PRINT*, "GRO>2 ANGLES(2)", T0_IMP(I),KA_IMP(I)
      ELSEIF(K==5) THEN ! IMPROPER FOURIER (K, C0, C1)
    	  READ(24,*) J, KB, X, Y
    	  PRINT*, I, "IMPROPER_FOURIER ZERO", KB, X, Y
    	  T0_IMP(I) = 0.0    	
    	  KA_IMP(I) =KB*X/8100.0 *CAL2JOUL/DEG2RAD/DEG2RAD
    	  TX_IMP(I) = 0.0
    	  PRINT*, "GRO>2 ANGLES(1)", T0_IMP(I),KA_IMP(I)     
      END IF 
    END DO 
    PRINT*, "IMPROPER PARAMETERS READ FOR MOF", MOF_COUNT
    READ(24,*) ! BLANK
    READ(24,*) ! PAIR
    DO I = 1, NLAMMPS_ATOMT
      READ(24,*) J, X, Y  ! kcal/mol , sig
      PRINT*, I, "SIG (A) = ", Y, "EPS (Kcal/mol) =", X
      SIG_A(I) =  Y/10.0D0
      EPS_A(I) =  X*CAL2JOUL 
      PRINT*, "GRO>", SIG_A(I), EPS_A(I), LAB_AT(I)
    END DO 
    READ(24,*) !
    PRINT*, "PAIR PARAMETERS READ FOR MOF", MOF_COUNT
    
    ! WRITING TOPOLOGY DATA
    READ(24,*) ! ATOMS
    READ(24,*) ! BLANK
    READ(21,*) (DIS_VECTOR(I), I =1, 3)
    PRINT*, "DISPLACE VECTOR FOR MOR", MOF_COUNT, DIS_VECTOR(:) 
    DIS_VECTOR = DIS_VECTOR /10.0
	  DO I = 1,NATOM_MOF 
      READ(24,*) J, J, K, X,  RX(I),RY(I),RZ(I)    
      DO J = 1, NCONDITION
        IF(K==IF_THIS(J)) K=THEN_THIS(J)
      END DO   
      K= K+ATOMT_IN_MOF
      ATOM2ATOMT(I) = K   
	    WRITE(31,'(I6, A8, I3, 2A5,i3, F10.5)') &
	    I, TRIM(LAB_AT(K)), MOF_COUNT, TRIM(MOL_LAB(MOF_COUNT)), TRIM(LAB_A(K)) , MOF_COUNT, CHARG_MOF(I)    
	    WRITE(35,'(I6, A8, I3, 2A5,i3, F10.5)') &
	    I, TRIM(LAB_AT(K)), MOF_COUNT, TRIM(MOL_LAB(MOF_COUNT)), TRIM(LAB_A(K)) , MOF_COUNT, CHARG_MOF(I)    
	  END DO  
	  READ(24,*) ! 
	  WRITE(31,*) !
	  WRITE(35,*) !
	  RX = RX/10.0 + DIS_VECTOR(1)
	  RY = RY/10.0 + DIS_VECTOR(2)
	  RZ = RZ/10.0 + DIS_VECTOR(3)
    PRINT*, "ATOM POSITION READ AND WRITTEN FOR MOF", MOF_COUNT

	  ! GET CONNECTION	
	  READ(24,*) ! BOND
	  READ(24,*)
	  WRITE(31,*) "[ bonds ]"
    WRITE(31,*) ";ai, aj, funct, B0(nm), KB(kJ/mol/nm2)"
	  WRITE(35,*) "[ bonds ]"
    WRITE(35,*) ";ai, aj, funct, B0(nm), KB(kJ/mol/nm2)"
    DO I = 1, NBOND !! CHECH THE NUMBERS
      READ(24,*) J, K, AI, AJ
      WRITE(31,*) AI, AJ, 1, R0_BOND(K), KB_BOND(K), ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ))
      WRITE(35,*) AI, AJ, 1, R0_BOND(K), KB_BOND(K), ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ))
    END DO 
    READ(24,*) !
    WRITE(31,*)!
    WRITE(35,*)!
    PRINT*, "BOND INFORMATION WRITTEN FOR MOF", MOF_COUNT
    WRITE(31,*) ! ANGLES
    WRITE(31,*) "; ANGLES"
    WRITE(35,*) ! ANGLES
    WRITE(35,*) "; ANGLES"
    READ(24,*)  ! ANGLES
    READ(24,*)  !
    WRITE(31,*) "[ angles ]"  
    WRITE(31,*) ";a1, a2, a2, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol/deg2)"  
    WRITE(35,*) "[ angles ]"  
    WRITE(35,*) ";a1, a2, a2, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol/deg2)"  
    DO I = 1, NANGLE
      READ(24,*) J, K, AI, AJ, AK  !, FUNCT_PARA, B0, KB
      IF(TYPE_ANG(K) == 1) THEN ! COSINE > harmonic(1)
        WRITE(31,*) AI, AJ, AK, 1, T0_ANG(K), KA_ANG(K), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! ANG, KA
        WRITE(35,*) AI, AJ, AK, 1, T0_ANG(K), KA_ANG(K), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! ANG, KA
      ELSEIF(TYPE_ANG(K)==2) THEN ! FOURIER > cosine harmonic(2)
      	WRITE(31,*) AI, AJ, AK, 2, T0_ANG(K), KA_ANG(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! Ang, KA
      	WRITE(35,*) AI, AJ, AK, 2, T0_ANG(K), KA_ANG(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! Ang, KA
      ELSEIF(TYPE_ANG(K)==3) THEN ! HARMONIC >  harmonic(3)
      	WRITE(31,*) AI, AJ, AK, 1, T0_ANG(K), KA_ANG(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! Ang, KA
      	WRITE(35,*) AI, AJ, AK, 1, T0_ANG(K), KA_ANG(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), K ! Ang, KA
      END IF 
    END DO 
    READ(24,*) !BLANK
    PRINT*, "ANGLE INFORMATION WRITTEN"
    OPEN(33, FILE="TOPOL_DIHEDRAL.dat")
    WRITE(33,*) ! DIHEDRAL
    WRITE(33,*) "[ dihedrals ]"
    WRITE(36,*) ! DIHEDRAL
    WRITE(36,*) "[ dihedrals ]"
    READ(24,*)  ! DIHEDRAL
    READ(24,*)  !
    WRITE(33,*) ";a1, a2, a3, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol/deg2)"  
    WRITE(36,*) ";a1, a2, a3, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol/deg2)"  
    DO I = 1, NDIHEDRAL
      READ(24,*) J, K, AI, AJ, AK, AL  !, FUNCT_PARA, B0, KB
      IF(TYPE_DHI(K) == 1) THEN ! COSINE >  periodic dihedral
        WRITE(33,*) AI, AJ, AK, AL, 9, T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !
        WRITE(36,*) AI, AJ, AK, AL, 9, T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !
      ELSEIF(TYPE_DHI(K)==2) THEN ! FOURIER >  RB    !!!! FOURIER DIHEDRAL >  APPROX PERIODIC DIHEDRAL
    	  WRITE(33,*) AI, AJ, AK, AL, 9,  T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)),  ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    	  WRITE(36,*) AI, AJ, AK, AL, 9,  T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)),  ";",&
    	  LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      ELSEIF(TYPE_DHI(K)==3)	THEN ! HARMONIC >  periodic dihedral 
        WRITE(33,*) AI, AJ, AK, AL, 9, T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !     
        WRITE(36,*) AI, AJ, AK, AL, 9, T0_DHI(K), KA_DHI(K), INT(TX_DHI(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !     
      END IF 
    END DO 
    READ(24,*) !BLANK  
    PRINT*, "DIHEDRAL INFORMATION WRITTEN FOR MOF", MOF_COUNT

    WRITE(31,*) ! IMPROPER
    WRITE(31,*) "; IMPROPER"
    WRITE(35,*) ! IMPROPER
    WRITE(35,*) "; IMPROPER"
    READ(24,*)  ! IMPROPER
    READ(24,*)  !
    WRITE(31,*) ";a2, a1, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol)"  
    WRITE(31,*) ";a3, a1, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol)"  
    WRITE(35,*) ";a2, a1, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol)"  
    WRITE(35,*) ";a3, a1, a4, FUNCT_PARA, THETA(deg), KB(kJ/mol)"  
    DO I = 1, NIMPROPER
      READ(24,*) J, K, AI, AJ, AK, AL  !, FUNCT_PARA, B0, KB
      IF(TYPE_IMP(K) == 1) THEN ! COSINE > periodic dihedral 
        WRITE(31,*) AI, AJ, AK, AL, 9, T0_IMP(K), KA_IMP(K), INT(TX_IMP(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !
        WRITE(35,*) AI, AJ, AK, AL, 9, T0_IMP(K), KA_IMP(K), INT(TX_IMP(K)), ";",&
        LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !IJJK, ANG, K, N"    !
      ELSEIF(TYPE_IMP(K)==2) THEN ! FOURIER > RB
      	WRITE(31,*) AI, AJ, AK, AL, 3, T0_IMP(K), KA_IMP(K), TX_IMP(K), 0.0, 0.0, 0.0, ";",&
      	LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      	WRITE(35,*) AI, AJ, AK, AL, 3, T0_IMP(K), KA_IMP(K), TX_IMP(K), 0.0, 0.0, 0.0, ";",&
      	LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      ELSEIF(TYPE_IMP(K)==4) THEN ! IMPROPER_FOURIER > COS HARMONIC
      	WRITE(31,*) ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL))
      	WRITE(31,*)  AJ, AI, AL, 2, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    	  WRITE(31,*)  AK, AI, AL, 2, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      	WRITE(35,*) ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL))
      	WRITE(35,*)  AJ, AI, AL, 2, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    	  WRITE(35,*)  AK, AI, AL, 2, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    
      ELSEIF(TYPE_IMP(K)==5) THEN ! IMPROPER_FOURIER ZERO > COS HARMONIC
      	WRITE(31,*) ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL))
      	WRITE(31,*)  AJ, AI, AL, 1, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    	  WRITE(31,*)  AK, AI, AL, 1, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      	WRITE(35,*) ";", LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AL))
      	WRITE(35,*)  AJ, AI, AL, 1, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AJ)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
    	  WRITE(35,*)  AK, AI, AL, 1, T0_IMP(K), KA_IMP(K), ";",&
    	  LAB_AT(ATOM2ATOMT(AK)), LAB_AT(ATOM2ATOMT(AI)), LAB_AT(ATOM2ATOMT(AL)), K !C0, C1, C2, c3c4c5"
      END IF 
    END DO 
    !READ(24,*) !BLANK  
    PRINT*, "IMPROPER INFORMATION WRITTEN FOR MOF", MOF_COUNT      
 
    DO I = 1, NATOM_MOF
  	  J = ATOM2ATOMT(I)  

      WRITE(32,GROFORMAT) MOF_COUNT, TRIM(ADJUSTL(MOL_LAB(MOF_COUNT))), TRIM(ADJUSTL(LAB_A(J))), OVERALL_COUNT+I, RX(I), RY(I), RZ(I)
    END DO 
    OVERALL_COUNT = OVERALL_COUNT + NATOM_MOF
    ATOMT_IN_MOF = ATOMT_IN_MOF + COUNT(MOF_ID==MOF_COUNT)
    PRINT*, "COUNT: ", COUNT(MOF_ID==MOF_COUNT)
    READ(22,*) !
    ! transfer of dihedral
    CLOSE(35) ! TEMP NEW_TOPOL2 
    CLOSE(36) ! TEMP TOPOL
    !CALL SYSTEM('cat TOPOL_DIHEDRAL.dat >> New_topol.top')
    CALL SYSTEM('cat TOP_TEMP.top >> New_topol2.top')
    CALL SYSTEM('cat DIHEDRAL_TEMP.dat >> New_topol2.top')
  END DO  ! NEXT MOF 
  CLOSE(31) ! TOPOLOGY FILE WRITTEN


  ATOM_COUNT = DOT_PRODUCT(NMOL(1:NTYPE_MOF), ATOM_IN_MOF)
  PRINT*, "TOTAL MOF ATOMS", ATOM_COUNT
  OPEN(36, FILE="DIHEDRAL_TEMP.dat")
  READ(21,*) !
  READ(22,*) !
  MOL_COUNT = 1
	DO FLUID_COUNT = 1, NTYPE_FLUID
    J=NTYPE_MOF+FLUID_COUNT
		ATOM_COUNT2 = 0
	  READ(21,*) MOL_LAB(J), NMOL(J), MOL_ATOM
	  READ(21,*) FLUID_MOLECULE_CHARGE, HEAD_FLUID_PARA
	  READ(21,*) FLUID_COORDINATES,     HEAD_FLUID_COORD
	  READ(21,*) DIS_VECTOR(1),DIS_VECTOR(2),DIS_VECTOR(3)
	  READ(22,*) MOL_LAB(J), NREXC
	  PRINT*, J, "FLUID DISPLACE", DIS_VECTOR
	  WRITE(33,*) !	  
	  WRITE(33,*) "[ moleculetype ]"
	  WRITE(33,*) "; MOLENAME, NEAR EXCLUSION"
	  WRITE(33,*) TRIM(MOL_LAB(J)), NREXC 
	  WRITE(36,*) !	  
	  WRITE(36,*) "[ moleculetype ]"
	  WRITE(36,*) "; MOLENAME, NEAR EXCLUSION"
	  WRITE(36,*) TRIM(MOL_LAB(J)), NREXC 
	  OPEN(23,FILE=TRIM(FLUID_MOLECULE_CHARGE), STATUS='OLD')
	  DO I = 1, HEAD_FLUID_PARA
	    READ(23,*) !
	  END DO 
	  WRITE(33,*) ! ATOMS
	  WRITE(33,*) "[ atoms ]"
	  WRITE(33,*) "; nr, at_type, res_nr, res_name, at_name, cgnr, charge MASS"		  
	  WRITE(36,*) ! ATOMS
	  WRITE(36,*) "[ atoms ]"
	  WRITE(36,*) "; nr, at_type, res_nr, res_name, at_name, cgnr, charge MASS"		  
	  DO I = 1, MOL_ATOM
	  	ATOM_COUNT2 = ATOM_COUNT2 + 1
	    READ(23,*) LAB, ATLAB, X, Y, Z, CHARG 
	    WRITE(33,'(I6, A8, I3, 2A5,i3, F10.5)') &
	    ATOM_COUNT2, TRIM(ATLAB), J, TRIM(MOL_LAB(J)), TRIM(LAB) , 1, CHARG         
	    WRITE(36,'(I6, A8, I3, 2A5,i3, F10.5)') &
	    ATOM_COUNT2, TRIM(ATLAB), J, TRIM(MOL_LAB(J)), TRIM(LAB) , 1, CHARG         
	  END DO 
    WRITE(33,*) ! BONDS
	  WRITE(33,*) "[ bonds ]"
    WRITE(33,*) ";ai, aj, funct, B0(nm), KB(kJ/mol/nm2)"
    WRITE(36,*) ! BONDS
	  WRITE(36,*) "[ bonds ]"
    WRITE(36,*) ";ai, aj, funct, B0(nm), KB(kJ/mol/nm2)"
    READ(23,*)  NBOND
    DO I = 1, NBOND
      READ(23,*) AI, AJ, FUNCT_PARA,  B0, KB
      WRITE(33,'(2I6,I2,2F15.5)') AI, AJ, FUNCT_PARA, B0/10.0D0, KB*100.0D0
      WRITE(36,'(2I6,I2,2F15.5)') AI, AJ, FUNCT_PARA, B0/10.0D0, KB*100.0D0
    END DO 

    WRITE(33,*) ! ANGLES
    WRITE(33,*) "[ angles ]"  
    WRITE(33,*) ";a1, a2, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol)"
    WRITE(36,*) ! ANGLES
    WRITE(36,*) "[ angles ]"  
    WRITE(36,*) ";a1, a2, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol)"
    READ(23,*) NANGLE
    DO I = 1, NANGLE
      READ(23,*) AI, AJ, AK, FUNCT_PARA, B0, KB
      WRITE(33,*) AI, AJ, AK, FUNCT_PARA, B0, KB      
      WRITE(36,*) AI, AJ, AK, FUNCT_PARA, B0, KB      
    END DO 

    WRITE(33,*) ! DIHEDRALS
    WRITE(33,*) "[ dihedrals ]"  
    WRITE(33,*) ";a1, a2, a3, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol)"
    WRITE(36,*) ! DIHEDRALS
    WRITE(36,*) "[ dihedrals ]"  
    WRITE(36,*) ";a1, a2, a3, a3, FUNCT_PARA, THETA(deg), KB(kJ/mol)"
    READ(23,*) NDIHEDRAL
    DO I = 1, NDIHEDRAL
      READ(23,*)  AI, AJ, AK, AL,  FUNCT_PARA, B0, KB
      WRITE(33,*) AI, AJ, AK, AL,  FUNCT_PARA, B0, KB
      WRITE(36,*) AI, AJ, AK, AL,  FUNCT_PARA, B0, KB
    END DO 

    WRITE(33,*) ! CONSTRAINTS
    WRITE(33,*) "[ constraints ]"  
    WRITE(33,*) ";a1, a2, FUNCT_PARA, DIST(nm) "
    WRITE(36,*) ! CONSTRAINTS
    WRITE(36,*) "[ constraints ]"  
    WRITE(36,*) ";a1, a2, FUNCT_PARA, DIST(nm) "
    READ(23,*) NCONSTRAIN
    DO I = 1, NCONSTRAIN
      READ(23,*)  AI, AJ, FUNCT_PARA, B0
      WRITE(33,*) AI, AJ,  FUNCT_PARA, B0/10.0D0
      WRITE(36,*) AI, AJ,  FUNCT_PARA, B0/10.0D0
    END DO 
    ! CONFIGURATION
    OPEN(UNIT=26,FILE=TRIM(FLUID_COORDINATES), STATUS='OLD')
    DO I = 1, HEAD_FLUID_COORD
    	READ(26,*) !
    END DO 
    DO K = 1, NMOL(J)
    	MOL_COUNT = MOL_COUNT + 1
    	DO I = 1, MOL_ATOM
    		ATOM_COUNT = ATOM_COUNT + 1
        READ(26,*) LAB1, X, Y, Z, LAB2, AI, LAB3, LAB 
        WRITE(32,GROFORMAT) MOL_COUNT, TRIM(ADJUSTL(MOL_LAB(J))), TRIM(ADJUSTL(LAB)), &
        ATOM_COUNT, (X+DIS_VECTOR(1))/10.0, (Y+DIS_VECTOR(2))/10.0, (Z+DIS_VECTOR(3))/10.0 
      END DO  
    END DO 
    CLOSE(26)
	END DO
  
  WRITE(33,*) !
  WRITE(33,*) "[ system ]"
  WRITE(33,*) (TRIM(MOL_LAB(I)), I = 2, NTYPE_MOL), " IN ", TRIM(MOL_LAB(1))
  WRITE(33,*) !
  WRITE(33,*) "[ molecules ]"
  WRITE(33,*) ";MOL NAME , NUMBER"
  WRITE(36,*) !
  WRITE(36,*) "[ system ]"
  WRITE(36,*) (TRIM(MOL_LAB(I)), I = 2, NTYPE_MOL), " IN ", TRIM(MOL_LAB(1))
  WRITE(36,*) !
  WRITE(36,*) "[ molecules ]"
  WRITE(36,*) ";MOL NAME , NUMBER"
  DO I = 1, NTYPE_MOF
    WRITE(33,*) TRIM(MOL_LAB(I)), 1  ! MOF WILL BE ALWAYS BE ONE.
    WRITE(36,*) TRIM(MOL_LAB(I)), 1  ! MOF WILL BE ALWAYS BE ONE.
  END DO  
  DO J = 1, NTYPE_FLUID
     I = NTYPE_MOF + J
     WRITE(33,*) TRIM(MOL_LAB(I)), NMOL(I)
     WRITE(36,*) TRIM(MOL_LAB(I)), NMOL(I)
  END DO  

  WRITE(32,*) AX, BY, CZ, 0.0, 0.0, BX, 0.0, CX, CY 
  CLOSE(33)
  CLOSE(32)
  CALL SYSTEM('cat TOPOL_DIHEDRAL.dat >> New_topol.top')
  CALL SYSTEM('cat DIHEDRAL_TEMP.dat >> New_topol2.top')

  OPEN(32,FILE="GRO_TEMP.dat")
  WRITE(32,*) MOL_COUNT, "NEW SYSTEM"  
  WRITE(32,*) ATOM_COUNT
  CLOSE(32)
  CALL SYSTEM('cat NEW_SYSTEM.gro >> GRO_TEMP.dat ; &
  	           mv GRO_TEMP.dat NEW_SYSTEM.gro') 

END PROGRAM MAKE_TOPOL
