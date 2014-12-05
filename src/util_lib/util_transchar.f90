      SUBROUTINE Util_TransChar(Chara, MLeng)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=MLeng) :: Chara
      INTEGER :: MLeng
!
      INTEGER :: ILeng
!
!     o Transform lower-case letters to upper-case (capital) ones
!
      DO ILeng = 1, MLeng
         IF ((Chara(ILeng:ILeng) >= 'a') .AND. (Chara(ILeng:ILeng) <= 'z')) THEN
            Chara(ILeng:ILeng) = CHAR(ICHAR(Chara(ILeng:ILeng)) - (ICHAR('a') - ICHAR('A')))
         END IF
      END DO
!
      END SUBROUTINE
