! CALCULATES THE RDF BETWEEN A AND B
! READ LIST OF A, AND B AND COMPAIRE
! Edit : For intra species RDF
! Edit :  using the distance trigonal dimension.
PROGRAM RDF
  IMPLICIT NONE
  REAL*8, ALLOCATABLE, DIMENSION(:):: AX,AY,AZ,&
                                      BX,BY,BZ
  REAL*8:: X, Y, Z, BOXL(3), BOXL_INV(3), DR, DR_INV,DR_2,&
           RCUT, RCUT2, RIJ(3), RIJ2, COM_FACT, XY, XZ, YZ, &
           BOX_VEC(3,3)
  INTEGER, ALLOCATABLE:: LOCAL_NUM(:)
  INTEGER:: I, J, K, NBIN, HEAD_TRAJ,TAIL_TRAJ, NTRAJ, NSKIP,&
            NATOM, A_ATOM, B_ATOM, A_COUNT, B_COUNT,&
            FRAME_SIZE, TRAJ_COUNT, IG, DISP_FREQ, ABOX(3)
  CHARACTER*50:: TRAJFILE, LIST_FILE_A, LIST_FILE_B, LAB  ,groformat
  LOGICAL, ALLOCATABLE:: A_LIST(:), B_LIST(:)
  LOGICAL:: INTRA_SPECIE
  

  groformat='(i5,2a5,i5,3f8.3)'
  OPEN(21,FILE='input_RDF.dat')
  READ(21,*) TRAJFILE, HEAD_TRAJ, TAIL_TRAJ
  READ(21,*) NTRAJ, NSKIP  
  READ(21,*) NATOM  
  READ(21,*) (BOXL(I), I = 1, 3)
  READ(21,*) XY, XZ, YZ
  READ(21,*) NBIN
  READ(21,*) LIST_FILE_A   
  READ(21,*) LIST_FILE_B
  READ(21,*) INTRA_SPECIE
  CLOSE(21)  
  BOXL_INV = 1.0D0 /BOXL
  BOX_VEC = 0.0D0
  DO I = 1, 3
    BOX_VEC(I,I) =  BOXL(I)
  END DO
  BOX_VEC(1,2) = XY
  BOX_VEC(1,3) = XZ
  BOX_VEC(2,1) = YZ
  
  ALLOCATE(A_LIST(NATOM), B_LIST(NATOM))   
  A_LIST = .FALSE.
  B_LIST = .FALSE.
  !IF(LIST_FILE_A == LIST_FILE_B) INTRA_SPECIE = .TRUE.

  OPEN(23, FILE=TRIM(LIST_FILE_A),STATUS='OLD')
  READ(23,*) A_ATOM
  DO I = 1, A_ATOM
    READ(23,*) J 
    A_LIST(J) = .TRUE.    
  END DO 
  CLOSE(23)
  IF (INTRA_SPECIE) THEN
    B_LIST = A_LIST
    B_ATOM = A_ATOM
  ELSE
    OPEN(23, FILE=TRIM(LIST_FILE_B),STATUS='OLD')
    READ(23,*) B_ATOM  
    DO I = 1, B_ATOM
      READ(23,*) J 
      B_LIST(J) = .TRUE.
    END DO 
    CLOSE(23)
  END IF
  ALLOCATE(AX(A_ATOM), AY(A_ATOM), AZ(A_ATOM), &
  	       BX(B_ATOM), BY(B_ATOM), BZ(B_ATOM), &
  	       LOCAL_NUM(0:NBIN))
  
  RCUT   = MINVAL(BOXL)/2.0D0
  RCUT2  = RCUT*RCUT
  DR     = RCUT/DBLE(NBIN)
  DR_INV = 1.0D0/DR 
  DR_2   = DR/2.0D0
  NSKIP  = NSKIP - 1

  PRINT*, "TOTAL BINS :   ", NBIN
  PRINT*, "BIN THICKNESS: ", DR
  PRINT*, "DR_INV         ", DR_INV
  PRINT*, "HALF DR        ", DR_2
  PRINT*, "TOTAL REF ATOMS", A_ATOM
  PRINT*, "TOTAL OTH ATOMS", B_ATOM
  PRINT*, "TOTAL FRAMES   ", NTRAJ
  PRINT*, "SKIP INTERVAL  ", NSKIP + 1
  PRINT*, "RCUT, RCUT2    ", RCUT, RCUT2
  PRINT*, "TRAJ FILE      ", TRIM(TRAJFILE)
  PRINT*, "REF LIST FILE  ", TRIM(LIST_FILE_A)
  PRINT*, "SAP LIST FILE  ", TRIM(LIST_FILE_B)
  PRINT*, "IS SAME SPECIES", INTRA_SPECIE
  

  OPEN(22,FILE=TRIM(TRAJFILE),STATUS='OLD')
  FRAME_SIZE = NSKIP*(HEAD_TRAJ + TAIL_TRAJ + NATOM)
  LOCAL_NUM  = 0
  DISP_FREQ  = INT(NTRAJ/100)
  DO TRAJ_COUNT = 1, NTRAJ
  	A_COUNT = 0 
  	B_COUNT = 0
    DO I = 1, HEAD_TRAJ
       READ(22,*) !SKIP
    END DO 
    DO I = 1, NATOM
      READ(22,groformat) j, LAB, LAB,J, X, Y, Z 
  	  IF (A_LIST(I)) THEN
	    A_COUNT = A_COUNT + 1
	    AX(A_COUNT) = X 
	    AY(A_COUNT) = Y 
	    AZ(A_COUNT) = Z
      END IF
      IF(.NOT. INTRA_SPECIE) THEN
        IF (B_LIST(I)) THEN
	      B_COUNT = B_COUNT + 1
	      BX(B_COUNT) = X 
	      BY(B_COUNT) = Y 
	      BZ(B_COUNT) = Z
        END IF
      END IF
    END DO 
    DO I = 1, TAIL_TRAJ
       READ(22,*) !SKIP
    END DO 
    DO J = 1, FRAME_SIZE
        READ(22,*) ! SKIP
    END DO  
    ! RDF CALCUATION
    IF (INTRA_SPECIE) THEN
      DO I = 1, A_ATOM
        X = AX(I)
        Y = AY(I)
        Z = AZ(I)
        DO  J = I+1, A_ATOM
          RIJ(1) = AX(J) - X
          RIJ(2) = AY(J) - Y
          RIJ(3) = AZ(J) - Z
          !RIJ    = RIJ - ANINT(RIJ*BOXL_INV)*BOXL
          CALL MIN_IMAGE_CONV(RIJ)
          RIJ2   = DOT_PRODUCT(RIJ, RIJ)
          IF (RIJ2 < RCUT2 ) THEN
            IG =  INT(SQRT(RIJ2)*DR_INV)
            LOCAL_NUM(IG) =  LOCAL_NUM(IG) + 2
          END IF 
        END DO  
      END DO      
    ELSE ! INTER SPEICES = NOT SAME
      DO I = 1, A_ATOM
        X = AX(I)
        Y = AY(I)
        Z = AZ(I)
        DO  J = 1, B_ATOM
          RIJ(1) = BX(J) - X
          RIJ(2) = BY(J) - Y
          RIJ(3) = BZ(J) - Z
          !RIJ    = RIJ - ANINT(RIJ*BOXL_INV)*BOXL
          CALL MIN_IMAGE_CONV(RIJ)          
          RIJ2   = DOT_PRODUCT(RIJ, RIJ)
          IF (RIJ2 < RCUT2 ) THEN
            IG =  INT(SQRT(RIJ2)*DR_INV)
            LOCAL_NUM(IG) =  LOCAL_NUM(IG) + 1
          END IF 
        END DO  
      END DO
    END IF  
  END DO  ! NTRAJ

  
  LOCAL_NUM(0) = 0 ! SELF COUNTING
  COM_FACT = BOXL(1)*BOXL(2)*BOXL(3)/&
             (4.0D0*ACOS(-1.0D0)*DR*DBLE(A_ATOM)*DBLE(B_ATOM)*DBLE(NTRAJ))
  PRINT*, "TOTAL TRAJ SAMPLED", TRAJ_COUNT, NTRAJ
  PRINT*, "COMMON FACTOR     ", COM_FACT 
  OPEN(31,FILE="RDF.dat")
  WRITE(31,*) "# R, GOR, LOCAL_NUM, CUMM "
  X=-DR_2
  K = 0
  DO I = 0, NBIN-1
  	X = X+DR        
  	J = LOCAL_NUM(I)
    K = K + J
    WRITE(31,*) X, DBLE(J)/X/X*COM_FACT, J, K
  END DO

  CONTAINS
  SUBROUTINE MIN_IMAGE_CONV(RIJ)
    REAL*8, INTENT(INOUT):: RIJ(3)
    ! MIC FOR THE TRIGONAL BOX
    
    ABOX(3) = ANINT( RIJ(3)*BOXL_INV(3))
    RIJ     = RIJ - ABOX(3)*BOX_VEC(:,3)
    ABOX(2) = ANINT( RIJ(2)*BOXL_INV(2))
    RIJ     = RIJ - ABOX(2)*BOX_VEC(:,2) 
    ABOX(1) = ANINT( RIJ(1)*BOXL_INV(1))
    RIJ     = RIJ - ABOX(1)*BOX_VEC(:,1)
    RETURN    
  END SUBROUTINE  MIN_IMAGE_CONV 
  
END PROGRAM RDF
