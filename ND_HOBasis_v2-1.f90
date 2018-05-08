!===================================
!    HO Basis for ND-CAHO v2.0
!   Created by: Tom  Halverson
!   Date: 5-29-2013
!
!      -Parallel
!===================================



    !+++++++++++++++++++++++++++++++++++++++++
    !++++++++++  Main Program Loop +++++++++++
    !+++++++++++++++++++++++++++++++++++++++++
Program Main
!----Modules----
 Use Inputs
 Use Scalapack_Params
!---------------

  Include 'mpif.h'

  Double Precision :: t1, t2, tot_time
  Character(len = 5) :: units

  !Total Time for the entire program
  t1 = MPI_Wtime()

! ------BLACS Intialization------
  Call BLACS_PINFO(IAM, NPROCS)

  If (IAM==0) Then
    Print *,'========================================================'
    Print *,'|     Harmonic Oscillator Basis Solver v2.0            |'
    Print *,'|     Created by: Thomas Halverson                     |'
    Print *,'|     Date: 5-29-13                                    |'
    Print *,'|                                                      |'
    Print *,'|  Solving N-D coupled/anharmonic oscillators using    |'
    Print *,'|  the harmonic oscillator basis                       |'
    Print *,'========================================================'
    Print *
  End If

! -----Gets Grid info from file------
  Call GetInputs() 

! -----Grid Intialization-----
  Call SL_INIT(ICTXT, NPROWS, NPCOLS)
  Call BLACS_GRIDINFO(ICTXT, NPROWS, NPCOLS, MYROW, MYCOL)

  Call GetOneDArrays()
  Call CreateMatrix() 
  Call Diagonalize()
  Call WriteVals()

  t2 = MPI_Wtime()

  If(t2-t1 > 18000.0d0) then
    tot_time = (t2-t1)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t2-t1 > 300) then
    tot_time = (t2-t1)/(60.0d0)
    units = 'mins.'
  Else 
    tot_time = (t2-t1)
    units = 'sec.'
  End If

  If(IAM == 0) Then
    Print 9999, tot_time, units
    Print *,'*****************************************************************'
    Print *
    Print *
  End If

  9999 Format ('  Total time taken: ', f10.5, 2x, a)
  
  Call BLACS_GRIDEXIT(ICTXT)
  Call BLACS_EXIT(0)
       
End Program Main
    !+++++++++++++++++++++++++++++++++++++++++
    !++++++  End of Main Program Loop ++++++++
    !+++++++++++++++++++++++++++++++++++++++++




    !************************************************************!
    !********************    SUBROUTINES   **********************!
    !************************************************************!
!==========================================================================================
  Subroutine GetInputs()  
   !----Modules-----
    Use Inputs
    Use ParameterArrays
    Use Scalapack_Params
   !----------------

    Implicit None
    integer :: i,j,k
    Character(len = 40) :: input_file
 
    Double Precision :: x1, test, h
   
    Call GETARG(1, input_file)
    input_file = trim(input_file)

    Open(UNIT = 15, FILE = input_file, STATUS='OLD',  ACTION ='Read')

    If(IAM==0) then
      Print *,'*****************************************************************'
      Print *,'STEP 1:'
      Print '(x,a,2x,a20)',"Reading data from: ", input_file
    End If

   !   ---- Read in files names (and directories) of user inputs ----  
    Open(UNIT = 15, FILE = input_file, STATUS='OLD',  ACTION ='Read')
    
      Read (15,'(A)') ODA_FN
      ODA_FN = trim(ODA_FN)
      Read (15,'(A)') PE_FN
      PE_FN = trim(PE_FN)
      Read (15,'(A)') masses_FN
      masses_FN = trim(masses_FN)
      Read (15,'(A)') SysParams_FN
      SysParams_FN = trim(SysParams_FN) 
      Read (15,'(A)') points_FN
      points_FN = trim(points_FN)   
      Read (15,'(A)') Out_FN
      Out_FN = trim(Out_FN)   

    Close(15)
!   -----------------------------------------


    If(Iam==0) then    
      Print *, ' Use requested files: '
      Print 1002, "PE:                       ", PE_FN
      Print 1002, "Masses:                   ", masses_FN
      Print 1002, "Energy truncation points: ", points_FN
      Print 1002, "Output:                   ", Out_FN
      Print *
      Print *, ' Reading data from suggested files'
    End If


!   ---- Read in user input system paramters from file ----  
    Open(UNIT = 15, FILE = SysParams_FN, STATUS='OLD',  ACTION ='Read')
    
      Read (15,*) NPROWS
      Read (15,*) NPCOLS
      Read (15,*) BLOCK_ROW
      BLOCK_COL = BLOCK_ROW
    
    Close(15)
!   -----------------------------------------

  

!   -----Read in PE info from file-----
    Open(UNIT = 15, FILE = PE_FN, STATUS='OLD',  ACTION ='Read')
    
      Read(15,*) Dmax, coeffs_len

      Allocate(Coeffs(coeffs_len))
      Allocate(Powers(coeffs_len, Dmax))
     
      Do i = 1, coeffs_len
        Read (15,*) (Powers(i,j), j=1, Dmax), Coeffs(i)
      End Do 
    
    Close(15)
!   -----------------------------------


!   ----Read masses from file-----
    Open(15, file = masses_FN, STATUS='OLD',  ACTION ='Read')
    
      Allocate(Masses(DMax))
     
      Do i = 1, Dmax
        Read(15,*) Masses(i)
        !Print *, Masses(i)
      End Do
            
    Close(15)
!   ------------------------------



!   ---- Read phase space points from file-----
    Open(15, file = points_FN, STATUS='OLD',  ACTION ='Read')

      Read(15,*) Emax, iMax
      Allocate(MandN(iMax,Dmax))
    
      Do i = 1, iMax
        Read(15,*) (MandN(i,j), j=1, Dmax)
      End Do
                
    Close(15)
!   ------------------------------

    If (IAM==0) then
      Print *, " Data read: SUCCESS"
      Print *
      Print *, "    _________________________________________________________"
      Print 1000, "Dimension:             ", Dmax
      Print 1000, "Number of parameters:  ", coeffs_len
      Print 1001, "Energy cuttoff:        ", EMax
      Print 1000, "Basis size:            ", iMax
      Print *, "    _________________________________________________________"
      Print *
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
    End If

!      Do i = 1, iMax
!        Do j = 1, iMax
!          Do k = 1, Dmax
!            test = x1(MandN(k,i),MandN(k,j),1.0d0,Masses(k),coeffs(k))
!            h = Sqrt( 1.0d0/(2.0d0*Masses(k)*coeffs(k)) )
!            Write(*,'(i1,2x,i1,2x,f10.5,5x)',Advance='NO') MandN(i,k),MandN(j,k),test
!          End Do
!          Print *
!        End Do
!      End do
      
      1000 Format(10x, a, i13)
      1001 Format(7x, a, f13.1)
      1002 Format(7x, a, a)
      1003 Format(10x, a, i13)
   
  End Subroutine GetInputs
!==========================================================================================





!==========================================================================================
Subroutine GetOneDArrays()

   !----Modules-----
    Use Inputs
    Use ParameterArrays
    Use Scalapack_Params
   !----------------

    Implicit NONE
    
    !Declare Variables
    Integer :: m,n 	!Indices for the Loop
    Integer :: index

    Double Precision :: Delta

    Character(len = 128) :: X1_FN   
    Character(len = 128) :: X2_FN    
    Character(len = 128) :: X3_FN       
    Character(len = 128) :: X4_FN    
    Character(len = 128) :: P2_FN

    If (IAM==0) then
      Print *,'*****************************************************************'
      Print *, '        STEP 2:'
      Print *,"   Initializing 1-D arrays"
      Print *
    End If
 
    Open(UNIT = 10, FILE = ODA_FN, STATUS='OLD',  ACTION ='Read')
      Read (10,'(A)') X1_FN
      X1_FN = trim(X1_FN)
      Read (10,'(A)') X2_FN
      X2_FN = trim(X2_FN)
      Read (10,'(A)') X3_FN
      X3_FN = trim(X3_FN) 
      Read (10,'(A)') X4_FN
      X4_FN = trim(X4_FN)   
      Read (10,'(A)') P2_FN
      P2_FN = trim(P2_FN)
    Close(10)   

    If (IAM==0) then
      Print *, '  ------------------------------------'
      Print *, '   Use requested files: '
      Print 1002, "X1: ", X1_FN
      Print 1002, "X2: ", X2_FN
      Print 1002, "X3: ", X3_FN
      Print 1002, "X4: ", X4_FN
      Print 1002, "P2: ", P2_FN
      Print *
      Print *, '   Reading data from suggested files'
      Print *, '  ------------------------------------'
      Print *
    End If
    
   
    Open(UNIT = 10, FILE = X1_FN ,STATUS='OLD',  ACTION ='Read')
    Open(UNIT = 11, FILE = X2_FN ,STATUS='OLD',  ACTION ='Read')
    Open(UNIT = 12, FILE = X3_FN ,STATUS='OLD',  ACTION ='Read')
    Open(UNIT = 13, FILE = X4_FN ,STATUS='OLD',  ACTION ='Read')
    Open(UNIT = 14, FILE = P2_FN ,STATUS='OLD',  ACTION ='Read')
  

    Do m = 1,100
      Do n = 1, 100
        OneDArrays_Array(m,n,1) = Delta(m,n)
      End Do
    End Do
       
   
      Do m = 1,100
        Do n = 1, 100
          Do index = 2, 6

            Read(index+8,*) OneDArrays_Array(m,n,index)
          
          End Do
        End Do
      End Do

    
    Close(UNIT = 10)
    Close(UNIT = 11)
    Close(UNIT = 13)
    Close(UNIT = 14)
   

    !Do index = 1, 6
    !  Do m = 1,7
    !    Do n = 1,7
    !     
    !      Write(*,'(f7.3)',Advance='NO') OneDArrays_Array(m,n,index)
    !      
    !    End Do
    !    Print *
    ! End Do
    !  Print *
    !  Print *
    !End Do



    If (IAM==0) Then
      Print*,'   1-D array initialization: SUCCESS'
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
    End If
    
    1002 Format(7x, a, a)

End Subroutine GetOneDArrays
!==========================================================================================




!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine CreateMatrix()
!   - CreateMatrix subroutine builds the Hamiltonian and overlap matrices
!     from the given list of truncated phase space points MandN
!
!-----------------------------------------------------------------------------
  Subroutine CreateMatrix()

   !----Modules-----
    Use Inputs
    Use ParameterArrays
    Use Outputs
    Use Scalapack_Params
   !----------------   
  
    Implicit NONE

    Integer :: i, j, power
    Integer :: eye, jay, kay
    Double Precision :: hbar
    Double Precision :: PE, KE, dummy
    Double Precision :: t1, t2, tot_time

    Character(len = 5) :: units
    
    !-------------------------BLACS Vars----------------------------
    Integer :: RSRC, CSRC, LLD, INFO
    Integer :: GLOBAL_rindex, GLOBAL_cindex
    Integer :: LOC_ROW, LOC_COL
    Integer :: LOC_roffset, LOC_coffset, LOC_rindex, LOC_cindex
    Integer :: ROW_first, COL_first, ROW_end, COL_end, ROW_start, COL_start
    !---------------------------------------------------------------
  
    !External  and BLACS functions
    Double Precision MPI_Wtime
    Integer numroc, indxg2p
    
    Integer, Dimension(Dmax) :: ems
    Integer, Dimension(Dmax) :: ens
    
    Double Precision, Dimension(iMax) :: Work !Needed for PDLAPRT
  
    hbar = 1.0d0

    t1 = MPI_Wtime()
    
    !----------------Intitialize Descriptor-------------------
    !-----------Preparing for block-cylic distribution--------
    !Proc offset
    RSRC = 0
    CSRC = 0
    
    !------------Find local array size----------------
    !-------------------------------------------------
    LOC_COL = numroc(iMax, BLOCK_COL, NPCOLS, CSRC, NPCOLS)
    LOC_COL = max(1,LOC_COL)
    LOC_ROW = numroc(iMax, BLOCK_ROW, NPROWS, RSRC, NPROWS)
    LLD = max(LOC_ROW,1)
    !-------------------------------------------------
  
    Call DESCINIT(DESCH, iMax, iMax, BLOCK_ROW, BLOCK_COL, RSRC, CSRC, ICTXT, LLD, INFO)
       
    Allocate(H(LOC_ROW, LOC_COL))
    
  
    If (IAM == 0) Then
      Print *,'*****************************************************************'
      Print *, '        STEP 3:'
      Print *, 'Preparing to create and distribute Hamiltonian matrix'
      Print *
      Print *, ' Distribution type: block-cyclic'
      Print 1000, 'Total number of processors: ', NPROCS
      Print 1001, 'Proccess grid:              ', NPROWS, ' x ', NPCOLS
      Print 1001, 'Global array size:          ', imax, ' x ', imax
      Print 1001, 'Block size:                 ', BLOCK_ROW, ' x ', BLOCK_COL
      Print 1001, 'Local array size:           ', LOC_ROW, ' x ', LOC_COL
      Print *
    End If
  
    1000 Format(2x, a, 6x, i4)
    1001 Format(2x, a, i10, 5x, a, i10)
  
    !---------------Begin Block-Cyclic Distribution---------
    !-------------------------------------------------------
  
    !Establish the starting point on each processor
    If (MYROW >= DESCH(RSRC_)) Then
      ROW_first = (MYROW - DESCH(RSRC_))*DESCH(BLOCK_ROWS_) + 1
    Else
      ROW_first = (MYROW + (NPROWS - DESCH(RSRC_))) * DESCH(BLOCK_ROWS_) + 1
    End If
  
    If (MYCOL >= DESCH(CSRC_)) Then
      COL_first = (MYCOL - DESCH(CSRC_)) * DESCH(BLOCK_COLs_) + 1
    Else
      COL_first = (MYCOL + (NPCOLS - DESCH(CSRC_))) * DESCH(BLOCK_COLS_) + 1
    End If
  
    call blacs_barrier(ICTXT, 'All')
  
    !Move through A by block
    Do COL_start = COL_first, DESCH(COLS_), NPCOLS * DESCH(BLOCK_COLS_)
      Do ROW_start = ROW_first, DESCH(ROWS_), NPROWS * DESCH(BLOCK_ROWS_)
  
        !find the last index in the block
        ROW_end = min( DESCH(ROWS_), ROW_start + DESCH(BLOCK_ROWS_)-1)
        COL_end = min( DESCH(COLS_), COL_start + DESCH(BLOCK_COLS_)-1)
  
        !find the local offset
        call infog2l(ROW_start, COL_start, DESCH, NPROWS, NPCOLS, MYROW, MYCOL, &
          &  LOC_roffset, LOC_coffset, RSRC, CSRC)
  
        !Move through each block
        Do GLOBAL_cindex = COL_start, COL_end
          Do GLOBAL_rindex = ROW_start, ROW_end
  
            LOC_rindex = LOC_roffset + (GLOBAL_rindex - ROW_start)
            LOC_cindex = LOC_coffset + (GLOBAL_cindex - COL_start)
        
            !Fills up the m and n indices arrays
            
            Do eye = 1, Dmax
              ems(eye) = MandN(GLOBAL_rindex,eye)+1
              ens(eye) = MandN(GLOBAL_cindex,eye)+1
            End Do 
            
  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Kenetic Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
        
            KE = 0.0d0
            Do jay = 1, Dmax
        
              dummy = (hbar**2)/(2.0d0*Masses(jay))*OneDArrays_array(ems(jay),ens(jay),6)
          
              Do kay = 1, jay - 1
                dummy = dummy * OneDArrays_array(ems(kay),ens(kay),1)
              End Do
              Do kay = jay + 1, Dmax            
                dummy = dummy * OneDArrays_array(ems(kay),ens(kay),1)
              End Do
          
              KE = KE + dummy

            End Do
              
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Potential Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
            PE = 0.0d0
            Do jay = 1, coeffs_len
              dummy = 1.0d0
              Do kay = 1, Dmax
            
                power = Powers(jay,kay) + 1
                dummy = dummy * OneDArrays_Array(ems(kay),ens(kay), power)
        
              End Do
              PE = PE + Coeffs(jay)*dummy
            End Do
            
            H(LOC_rindex,LOC_cindex) =  PE + KE
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
        
          End Do
        End Do
      End Do
    End Do

    !Do eye = 1, Imax
    !  Do jay= 1, iMax
    !    Write (*,'(2x,f4.2)',Advance = 'NO') H(i,j)
    !  End Do
    !  Print *
    !End Do
  
    call blacs_barrier(ICTXT, 'All')
  
    !Call pdlaprnt(Imax, Imax, H, 1, 1, DESCH, 0, 0, 'H', 6, work)

    t2 = MPI_Wtime()

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If
    
    If (IAM == 0) then
      Print *, 'Distribution: COMPLETE'
      Print *
      Print 9999, tot_time,units
      Print *,'*****************************************************************'  
      Print *
      Print *
      Print *  
    End If  

    Deallocate(Coeffs)
    Deallocate(Powers)
    Deallocate(Masses)
    Deallocate(MandN)
  
    9999 Format ('  Time taken to create and distribute: ',f10.5, 2x, a)
  
  End SubRoutine CreateMatrix
!-----------------------------------------------------------------------------
    



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine Diagonalize()
!   - Diagonalize subroutine calculates the eigenvalues of the Hamiltonian
!    using the Scalapack routine PDSYGVX
!
!------------------------------------------------------------------------------
SubRoutine Diagonalize()

   !----Modules-----
    Use Inputs
    Use Outputs
    Use Scalapack_Params
   !----------------  

	Implicit NONE
    
	Integer           :: IA, JA
    Integer           :: IZ, JZ, LWORK, INFO
    Character*1       :: JOBZ, UPLO
    Double Precision  :: t1, t2, tot_time
    Integer           :: NB, NN, NP0, NNP

    Character(len = 5) :: units
    
    Integer,                       Dimension(DLEN)            :: DESCZ
	Double Precision, Allocatable, Dimension(:,:)             :: Z
    Double Precision, Allocatable, Dimension(:)               :: WORK
    
    
	Double Precision PDLAMCH, MPI_Wtime, NUMROC
   

    JOBZ = 'N'
    UPLO = 'U'
	IA = 1
    JA = 1
    IZ = 1
    JZ = 1
   
    NB = DESCH( BLOCK_ROWS_ )
    NN = MAX( iMax, NB, 2 )
    NP0 = NUMROC( NN, NB, 0, 0, NPROWS )
    NNP = MAX( iMax, NPROWS*NPCOLS + 1, 4 )
    

    t1 = MPI_Wtime()
 
    If(IAM==0) Then
      Print *,'*****************************************************************'
      Print *, '        STEP: 4'
      Print *, 'Diagonalizing Hamiltonian'
      Print *
      Print *, ' Using direct diagonalization'
      Print *, ' Scalapack routine:   PDSYGVX'
      Print *
    End If

    Allocate(EigenVals(iMax))
    !Allocate(Z(DESCH(LLD_),DESCH(LLD_)))
    Allocate(Z(10,10))

    Call DESCINIT(DESCZ, iMax, iMax, BLOCK_ROW, BLOCK_COL, DESCH(RSRC_), DESCH(CSRC_), &
                   &  ICTXT, DESCH(LLD_), INFO)

    If(iMax < 75000) then
      
      LWORK = -1
      Allocate(WORK(1))
      
      Call PDSYEV( JOBZ, UPLO, iMax, H, IA, JA, DESCH, EigenVals, Z, IZ, JZ, DESCZ, WORK, LWORK, INFO )
                        
      LWORK = WORK(1)
      
      Deallocate(WORK)
      Allocate(WORK(LWORK))
      

    Else 
    
      LWORK = 9 * iMax + MAX( 5 * NN, NB * ( NP0 + 1 ) )

      Allocate(WORK(LWORK))
      

    End If
      
    If(IAM==0) Then
      Print *, DESCH(LLD_)
      Print *, '       Workspace Sizes:'
      Print *, '------------------------------------'
      Print *, '  Work array length: ',LWORK
      Print *, '------------------------------------'
      Print *
    End If
    

    call blacs_barrier(ICTXT, 'All')

    !Print *,'GOTHERE'

    Call PDSYEV( JOBZ, UPLO, iMax, H, IA, JA, DESCH, EigenVals, Z, IZ, JZ, DESCZ, WORK, LWORK, INFO )

    t2 = MPI_Wtime()

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If

    If (IAM == 0) Then
      Print *, 'Diagonalization: COMPLETE'
      Print *
      Print 9999, tot_time, units
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
    End If
    
    9999 Format ('  Time taken to diagonalize: ',f10.5,2x,a)

End SubRoutine
!-----------------------------------------------------------------------------



!==========================================================================================
  Subroutine WriteVals()
    
   !----Modules-----
    Use Inputs
    Use Outputs
    Use Scalapack_Params
    !----------------  

    Implicit NONE
    Integer	:: i
   
    If (IAM==0) then
      Print *,'*****************************************************************'
      Print *,'SUMMARY'
      Print *
      Print *, "   ",Dmax,'-D Coupled Anharmonic Oscillator'
      Print 1004, Dmax
      Print 1001, EMax
      Print 1002, iMax
    

      Open(UNIT = 15, FILE = Out_FN, ACTION='WRITE')    
    
      Do i = 1, iMax
        Write (15,*) EigenVals(i)
      End Do

      Close(UNIT = 15)
    
      Print 1003, Out_FN
      Print *
      Print *
    End If

    1001 Format ('  Energy truncation:           ', 3x, f7.1)
    1002 Format ('  Total number of eigenvalues: ', 3x, i7)
    1003 Format ('  Eigenvalues read to file:    ', a)
    1004 Format ('  Problem dimensionality:      ', i7)
    
    

End Subroutine WriteVals


    

!======================================
!==============Functions===============
!======================================

!===============================
  Double Precision Function Delta(m,n)

    Implicit NONE
    Integer, Intent(IN) :: m
    Integer, Intent(IN) :: n

    If (m==n) Then
      Delta = 1.0d0
    Else 
      Delta = 0.0d0
    End If

  End Function Delta
!===============================    