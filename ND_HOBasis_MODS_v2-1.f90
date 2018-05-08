Module Inputs
  Implicit NONE
  Integer :: iMax
  Integer :: DMax         !Problem Dimensionality
  Integer :: coeffs_len   !Total number of coulping/anharmonic constants
  Double Precision :: Emax
  
  Character(len = 128) :: Out_FN       !Output file: eigenvalues
  Character(len = 128) :: PE_FN        !Input file: potential energy
  Character(len = 128) :: ODA_FN       !Input file: 1-d Arrays
  Character(len = 128) :: masses_FN    !Input file: masses
  Character(len = 128) :: points_FN    !Input file: masses
  Character(len = 128) :: SysParams_FN  !Input file: Scalapack Grid Info

  !Blacs grid info 
  Integer		   :: NPROWS       !Number of processor rows
  Integer		   :: NPCOLS       !Number of processor columns
  Integer		   :: BLOCK_ROW    !Number of rows in each block
  Integer          :: BLOCK_COL    !NUmber of cloumns in each block
End Module Inputs

Module ParameterArrays
  Implicit NONE
  Integer,          Allocatable, Dimension(:,:) :: MandN
  Double Precision, Allocatable, Dimension(:)   :: Coeffs
  Integer,          Allocatable, Dimension(:,:) :: Powers
  Integer,          Allocatable, Dimension(:,:) :: Points
  Double Precision, Allocatable, Dimension(:)   :: Masses
  Double Precision, Dimension(100,100,6) :: OneDArrays_Array  
End Module ParameterArrays

  
Module Outputs
  Implicit NONE
  Double Precision, Allocatable, Dimension(:,:) :: H
  Double Precision, Allocatable, Dimension(:)   :: EigenVals
  !Double Precision, Allocatable, Dimension(:)  :: Ezero     
End Module Outputs


!-------------------------------------------------------------------------------------------
Module Scalapack_Params
  
  Implicit NONE

  !Blacs grid vars
  Integer            :: ICTXT      !BLACS context variable
  Integer 			 :: IAM, NPROCS
  Integer 			 :: MYROW, MYCOL
  
  Integer, Parameter :: DLEN = 9   !Length of descriptor DESC

  !Variables for each element of DESC_Variable_() array
  !DESC_Var_={DTYPE,ICTXT,ROWS,COLS,BLOCK_ROWS,BLOCK_COLS,RSRC,CSRC,LLD}
  Integer, Parameter :: BLOCK_CYCLIC_2D = 1
  Integer, Parameter :: DTYPE_          = 1 
  Integer, Parameter :: ICTXT_          = 2
  Integer, Parameter :: ROWS_ 			= 3     
  Integer, Parameter :: COLS_ 			= 4
  Integer, Parameter :: BLOCK_ROWS_		= 5
  Integer, Parameter :: BLOCK_COLS_		= 6
  Integer, Parameter :: RSRC_			= 7
  Integer, Parameter :: CSRC_			= 8
  Integer, Parameter :: LLD_			= 9
  
  Integer, Dimension(DLEN) :: DESCH

End Module Scalapack_Params
!-------------------------------------------------------------------------------------------