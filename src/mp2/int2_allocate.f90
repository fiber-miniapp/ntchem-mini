      SUBROUTINE Int2_Allocate
!
      USE MP2_Parameter_Module, ONLY : MaxAngl, MaxSgmt, MaxSgm2
      USE Int2_Module, ONLY : PreFactAB, PreFactCD, CCoefAB, CCoefCD, CContAB, CContCD, PX, PY, PZ, &
     &   PAX, PAY, PAZ, PBX, PBY, PBZ, QX, QY, QZ, QCX, QCY, QCZ, QDX, QDY, QDZ, &
     &   ExpntA, ExpntB, ExpntC, ExpntD, ExpntP, ExpntQ
!
      IMPLICIT NONE
!
      ALLOCATE(PreFactAB(MaxSgm2))
      ALLOCATE(PreFactCD(MaxSgm2))
      ALLOCATE(CCoefAB(MaxSgm2))
      ALLOCATE(CCoefCD(MaxSgm2))
      ALLOCATE(CContAB(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl))
      ALLOCATE(CContCD(MaxSgm2,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl,0:MaxAngl))
      ALLOCATE(PX(MaxSgm2))
      ALLOCATE(PY(MaxSgm2))
      ALLOCATE(PZ(MaxSgm2))
      ALLOCATE(PAX(MaxSgm2))
      ALLOCATE(PAY(MaxSgm2))
      ALLOCATE(PAZ(MaxSgm2))
      ALLOCATE(PBX(MaxSgm2))
      ALLOCATE(PBY(MaxSgm2))
      ALLOCATE(PBZ(MaxSgm2))
      ALLOCATE(QX(MaxSgm2))
      ALLOCATE(QY(MaxSgm2))
      ALLOCATE(QZ(MaxSgm2))
      ALLOCATE(QCX(MaxSgm2))
      ALLOCATE(QCY(MaxSgm2))
      ALLOCATE(QCZ(MaxSgm2))
      ALLOCATE(QDX(MaxSgm2))
      ALLOCATE(QDY(MaxSgm2))
      ALLOCATE(QDZ(MaxSgm2))
      ALLOCATE(ExpntA(MaxSgmt))
      ALLOCATE(ExpntB(MaxSgmt))
      ALLOCATE(ExpntC(MaxSgmt))
      ALLOCATE(ExpntD(MaxSgmt))
      ALLOCATE(ExpntP(MaxSgm2))
      ALLOCATE(ExpntQ(MaxSgm2))
!
      END SUBROUTINE
