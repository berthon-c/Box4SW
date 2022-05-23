MODULE Physique

  !!$ Module qui gere la physique pour Shallw-water  !
 
  USE var_types

  IMPLICIT NONE

CONTAINS

!***********************************************************************!
  FUNCTION ValeurPropre(Ua, nx, ny) result (Mvp)
    
    REAL(PR), DIMENSION(:), INTENT(IN) :: Ua
    REAL(PR)              , INTENT(IN) :: nx ,ny
    REAL(PR)               :: Mvp
    REAL(PR)               :: unorm
    REAL(PR)               :: ux, uy, c, h

    h  = Ua(1)
    ux = vitesse(h,Ua(2))
    uy = vitesse(h,Ua(3))

    c  = sqrt( gravitation*h )

    unorm = ux*nx + uy*ny

    Mvp   = abs(unorm) + c


  END FUNCTION ValeurPropre
!***********************************************************************!

!**********************************************************************!
  function vitesse(h,hu)
    
    real(PR), intent(in) :: h,hu
    real(PR)             :: vitesse
    
    vitesse = 0._PR
    if (h>eps) vitesse = hu/h
    
  end function vitesse
!**********************************************************************!

END MODULE Physique
