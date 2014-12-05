      SUBROUTINE RIMP2_Get_BatchInfo_MPI
!
      USE MP2_Module, ONLY : MP2BatchLv, NOccBat, LenOccBat, IOccBat, NOccBat_per_Proc, UHF, NAlpBet, NActO, NActV
      USE MPI_Module, ONLY : NProcsMO, MyRank
!
      IMPLICIT NONE
!
      INTEGER :: MXNActO, ModNActO, LenOccBat0
      INTEGER :: IAlpBet, IiBat, IiBat_Proc, IBgn
      INTEGER :: Irank, Icnt
!
      IF ((Mod(NProcsMO, 2) == 0) .AND. (MP2BatchLv == 0)) THEN
         MP2BatchLv = 1
      END IF
!
      NOccBat_per_Proc = 2 ** MP2BatchLv
      NOccBat = NOccBat_per_Proc * NProcsMO
      MXNActO = MAX(NActV(1), NActV(2))
      LenOccBat = MXNActO / NOccBat
      ModNActO = MOD(MXNActO, NOccBat)
      IF (ModNActO > 0) THEN
         LenOccBat = LenOccBat + 1
      END IF
!
      ALLOCATE(IOccBat(2,NOccBat,2))
!
      DO IAlpBet = 1, NAlpBet
         LenOccBat0 = NActV(IAlpBet) / NOccBat
         ModNActO = MOD(NActV(IAlpBet), NOccBat)
         DO IiBat = 1, NOccBat
            IOccBat(2,IiBat,IAlpBet) = LenOccBat0
         END DO
         IF (ModNActO > 0) THEN
            Icnt = ModNActO
            DO IiBat_Proc = 1, NOccBat_per_Proc
               DO Irank = 0, (NProcsMO - 1)
                  IiBat = Irank * NOccBat_per_Proc + IiBat_Proc
                  IOccBat(2,IiBat,IAlpBet) = LenOccBat0 + 1
                  Icnt = Icnt - 1
                  IF (Icnt == 0) GO TO 100
               END DO
            END DO
  100       CONTINUE
         END IF
!
         IBgn = 0
         DO IiBat = 1, NOccBat
            IOccBat(1,IiBat,IAlpBet) = IBgn
            IBgn = IBgn + IOccBat(2,IiBat,IAlpBet)
         END DO
!
      END DO
!
      IF (MyRank == 0) THEN
         WRITE(*, *) 'MP2BatchLv       = ', MP2BatchLv
         WRITE(*, *) 'NOccBat          = ', NOccBat
         WRITE(*, *) 'LenOccBat        = ', LenOccBat
         WRITE(*, *) 'NOccBat_per_proc = ', NOccBat_per_Proc
         WRITE(*, *) 'Batch infomation'
         WRITE(*, *) 'Alpha MO case'
         WRITE(*, *) 'Batch NO., NO. of AOs'
         DO IiBat = 1, NOccBat
            WRITE(*, '(2I10)') IiBat, IOccBat(2,IiBat,1)
         END DO
         IF (UHF) THEN
            WRITE(*, *) 'Beta MO case'
            WRITE(*, *) 'Batch NO., NO. of AOs'
            DO IiBat = 1, NOccBat
               WRITE(*, '(2I10)') IiBat, IOccBat(2,IiBat,1)
            END DO
         END IF
      END IF
!
      END SUBROUTINE
