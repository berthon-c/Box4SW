MODULE Conditions_limites

  !!$ Module qui gere les conditions limites                      !
  !---------------------------------------------------------------!
  !!$ Auteurs : C. Berthon
  !---------------------------------------------------------------!

  USE var_types
  USE Physique
 
  IMPLICIT NONE

CONTAINS

  !*******************************************************************!
  SUBROUTINE CL_evaluation_UA2 (Ua1, Ua2, icl)

    REAL(PR), DIMENSION(3),  INTENT(in)  :: Ua1
    REAL(PR), DIMENSION(3),  INTENT(out) :: Ua2

    INTEGER, INTENT(in) :: icl


    SELECT CASE (icl)

    CASE(3) 
       !CALL CL_Neumann(Ua1, Ua2)
       CALL CL_wall(Ua1, Ua2)

    CASE DEFAULT
       WRITE(6,*) " La condition aux limites  ", icl, " n est pas traitee"
       STOP
                
    END SELECT
                     
  END SUBROUTINE CL_evaluation_UA2

  !*******************************************************************!

  !************************************************************************!
  SUBROUTINE CL_Neumann(Ua1, Ua2)

    REAL(PR), DIMENSION(3),  INTENT(in)  :: Ua1
    REAL(PR), DIMENSION(3),  INTENT(out) :: Ua2
    
    Ua2 = Ua1

  END SUBROUTINE CL_Neumann
  !**************************************************************************!

  !************************************************************************!
  SUBROUTINE CL_wall(Ua1, Ua2)

    REAL(PR), DIMENSION(3),  INTENT(in)  :: Ua1
    REAL(PR), DIMENSION(3),  INTENT(out) :: Ua2
    
    Ua2(1) = Ua1(1)
    Ua2(2) = 0._PR
    Ua2(3) = 0._PR

  END SUBROUTINE CL_wall
  !**************************************************************************!

END MODULE Conditions_limites
