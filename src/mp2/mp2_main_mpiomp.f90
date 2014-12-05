      PROGRAM MP2_Main_MPIOMP
!
      USE MP2_Module, ONLY : InCore
      USE MP2_Basis_Module, ONLY : Spherical
      USE MPI_Module, ONLY : MyRank
!
!     o Main MP2 routine
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
      REAL(8) :: TimeBgn, TimeEnd, WTimeBgn, WTimeEnd
!
!     o Initialize MPI
!
      CALL Util_InitMPI(.TRUE.)
!
!     o Read MP2 input
!
      CALL MP2_Read_Input_MPI
!
!     o Determine parameters
!
      CALL MP2_MaxParam_MPI
      CALL RIMP2_MaxParam_MPI
!
!     o Read SCF_Info
!
      CALL MP2_Read_SCFInfo_MPI
!
!     o Memory allocation
!
      CALL MP2_Allocate
      CALL MP2_Basis_Allocate
      CALL RIMP2_Basis_Allocate
!
!     o Determine IOff
!
      CALL MP2_IOff
!
!     o Make index of basis sets
!
      CALL MP2_IndBasis
      CALL RIMP2_IndBasis
!
!     o Read geometry and basis information
!
      CALL MP2_Read_Geom_MPI
      CALL MP2_Read_Basis_MPI
      CALL RIMP2_Read_Basis_MPI
!
!     o Make density factors
!
      CALL MP2_FacDens
!
!     o Construct spherical Harmonics coefficients
!
      IF (Spherical) THEN
         CALL MP2_SphCoef
      END IF
!
!     o RI-MP2 energy calculation
!
      IF (MyRank == 0) THEN
         PRINT '(" ...   Enter    (RIMP2_Driver    )")'
      END IF
      WTimeBgn = MPI_WTIME()
      CALL CPU_TIME(TimeBgn)
      IF (InCore) THEN
!
!        o In-Core RI-MP2 algorithm
!
         CALL RIMP2_Driver_InCore_MPIOMP
!
      ELSE
!
!        o Semi-Direct RI-MP2 algorithm
!
         CALL RIMP2_Driver_SemiDirect_MPIOMP
!
      END IF
      CALL CPU_TIME(TimeEnd)
      WTimeEnd = MPI_WTIME()
      IF (MyRank == 0) THEN
         PRINT '(" ...   CPU time (RIMP2_Driver    ) :", F12.2)', TimeEnd - TimeBgn
         PRINT '(" ...   WALL time (RIMP2_Driver    ) :", F12.2)', WTimeEnd - WTimeBgn
      END IF
!
!     o Memory deallocation
!
      CALL MP2_Deallocate
      CALL MP2_Basis_Deallocate
      CALL RIMP2_Basis_Deallocate
!
!     o Finalize MPI
!
      CALL Util_FinMPI
!
      END PROGRAM
