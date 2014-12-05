      SUBROUTINE RIMP2_Inv2c_MPI
!
      USE RIMP2_Module, ONLY : NBF_RI, RI2cInv, RIInt2c
!
!     o Construct inverse of 2-center ERIs
!
      IMPLICIT NONE
!
      INTEGER :: Info   ! LAPACK
!
!     o Diagonalization of 2-center ERIs (P|Q) using Cholesky decomposition
!
!     --- Cholesky decomposition
!
      CALL Util_Tr1To2(RIInt2c, RI2cInv, NBF_RI)
      Info = 0
      ! write(*,*) "[RIMP2_Inv2c_MPI] DPOTRF, N:", NBF_RI  ! debug
      CALL DPOTRF('U', NBF_RI, RI2cInv, NBF_RI, Info)   ! LAPACK
      IF (Info /= 0) THEN
         CALL Util_AbortMPI('Error: DPOTRF')
      END IF
!
!     --- Inverse matrix
!
      Info = 0
      ! write(*,*) "[RIMP2_Inv2c_MPI] DTRTRI, N:", NBF_RI  ! debug
      CALL DTRTRI('U', 'N', NBF_RI, RI2cInv, NBF_RI, Info)   ! LAPACK
!!!      CALL DPOTRI('U', NBF_RI, RI2cInv, NBF_RI, Info)   ! LAPACK
      IF (Info /= 0) THEN
         CALL Util_AbortMPI('Error: DPOTRI')
      END IF
      CALL Util_LowTrMZero(RI2cInv, NBF_RI, NBF_RI)
!DebugWRITE(*, *) '+++++ RI2cInv +++++'
!DebugCALL Util_MatOut(RI2cInv, NBF_RI, NBF_RI)
!
      END SUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RIMP2_Inv2c_MPI_Old
!
      USE RIMP2_Module, ONLY : NBF_RI, RI2cInv, RIInt2c
      USE MP2_Constant_Module, ONLY : Zero, One
      USE MPI_Module, ONLY : MPIMain
!
!     o Construct inverse of 2-center ERIs
!
      IMPLICIT NONE
!
      REAL(8), PARAMETER :: ThrRI = 1.0D-13
      INTEGER :: I, J
      INTEGER :: LWork, LIWork, Info
      INTEGER :: NBF_RI0
      REAL(8) :: VD, VLinDepMax
      INTEGER, ALLOCATABLE :: IWork1(:)
      REAL(8), ALLOCATABLE :: VDiag(:)
      REAL(8), ALLOCATABLE :: RWork1(:)
      REAL(8), ALLOCATABLE :: RWork2A(:,:), RWork2B(:,:), RWork2C(:,:)
!
      ALLOCATE(RWork2A(NBF_RI,NBF_RI))
!
!     o Diagonalization of 2-center ERIs (P|Q)
!
      LWork = 1 + 6 * NBF_RI + NBF_RI ** 2
      LIWork = 3 + 5 * NBF_RI
      ALLOCATE(VDiag(NBF_RI))
      ALLOCATE(RWork1(LWork))
      ALLOCATE(IWork1(LIWork))
      Info = 0
      CALL DSPEVD('V', 'U', NBF_RI, RIInt2c, VDiag, RWork2A, NBF_RI, RWork1, LWork, IWork1, LIWork, Info)
      IF (Info /= 0) THEN
         CALL Util_AbortMPI('Error: Not converge in RIInt2_Inv2c')
      END IF
      DEALLOCATE(RWork1)
      DEALLOCATE(IWork1)
!
!     o Construct (P|Q)^(-1)
!
!     --- V^(-1) = C * VDiag^(-1) * CT
!
      ALLOCATE(RWork2B(NBF_RI,NBF_RI))
      ALLOCATE(RWork2C(NBF_RI,NBF_RI))
!
      NBF_RI0 = 0
      VLinDepMax = Zero
      DO J = 1, NBF_RI
         IF (VDiag(J) >= ThrRI) THEN
            NBF_RI0 = NBF_RI0 + 1
            VDiag(NBF_RI0) = One / VDiag(J)
            VD = VDiag(NBF_RI0)
            DO I = 1, NBF_RI
               RWork2B(I,NBF_RI0) = RWork2A(I,J) * VD
               RWork2C(I,NBF_RI0) = RWork2A(I,J)
            END DO
         ELSE
            VLinDepMax = MAX(VLinDepMax, ABS(VDiag(J)))
         END IF
      END DO
      IF (NBF_RI0 /= NBF_RI .AND. MPIMain) THEN
         WRITE(*, *) '--- Remove linear-dependence : ', NBF_RI, NBF_RI0, VLinDepMax
      END IF
      CALL DGEMM('N', 'T', NBF_RI, NBF_RI, NBF_RI0, One, RWork2B(1:NBF_RI,1:NBF_RI0), NBF_RI, &
     &   RWork2C(1:NBF_RI,1:NBF_RI0), NBF_RI, Zero, RI2cInv, NBF_RI)
!
      DEALLOCATE(VDiag)
      DEALLOCATE(RWork2A)
      DEALLOCATE(RWork2B)
      DEALLOCATE(RWork2C)
!
      END SUBROUTINE
