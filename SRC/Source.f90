MODULE Source

  !!$ Module qui gere les les termes source  !
  !---------------------------------------------------------------!
  !!$ Auteurs : C. Berthon
  !---------------------------------------------------------------!

  USE var_types
  USE Physique
 
  IMPLICIT NONE

CONTAINS

  !*******************************************************************!
  SUBROUTINE TermeSource (Stopo, Sfriction, Spluie, &
                  &    nx, ny, Ua1, Ua2, Z1, Z2, dx, Source)

    INTEGER,                 INTENT(in)    :: Stopo, Sfriction, Spluie
    REAL(PR),                INTENT(in)    :: nx, ny
    REAL(PR), DIMENSION(3),  INTENT(in)    :: Ua1, Ua2
    REAL(PR),                INTENT(in)    :: Z1, Z2, dx
    REAL(PR), DIMENSION(3),  INTENT(out)   :: Source

    REAL(PR), DIMENSION(3)                 :: UL, UR

    REAL(PR)    :: pluie
    REAL(PR)    :: friction_n, friction_t
    REAL(PR)    :: source_topo

    Source      = 0._PR
    pluie       = 0._PR
    friction_n  = 0._PR
    friction_t  = 0._PR
    source_topo = 0._PR
    
    !calcul des variables primitives dans le repere de l'interface
    UL(1) = Ua1(1)
    UL(2) = nx*vitesse(Ua1(1),Ua1(2)) + ny*vitesse(Ua1(1),Ua1(3))
    UL(3) =-ny*vitesse(Ua1(1),Ua1(2)) + nx*vitesse(Ua1(1),Ua1(3))

    UR(1) = Ua2(1)
    UR(2) = nx*vitesse(Ua2(1),Ua2(2)) + ny*vitesse(Ua2(1),Ua2(3))
    UR(3) =-ny*vitesse(Ua2(1),Ua2(2)) + nx*vitesse(Ua2(1),Ua2(3))

    !----------------------------------------
    ! terme source de topographie
    ! Stopo = 0 -> fond plat
    SELECT CASE (Stopo) 

    CASE(1) ! source topo centree
       CALL Source_topo_centree(UL,UR,Z1,Z2,dx,Source_topo)

    END SELECT

    Source(2) = Source(2) + nx*source_topo
    Source(3) = Source(3) + ny*source_topo

    !----------------------------------------
    ! terme source de friction
    ! Sfriction = 0 -> pas de friction
    SELECT CASE (Sfriction) 

    CASE(1) ! friction centree sur l'interface
       CALL Friction_centree(UL,UR,friction_n,friction_t)
       
    END SELECT

    Source(2) = Source(2) + nx*friction_n - ny*friction_t
    Source(3) = Source(3) + ny*friction_n + nx*friction_t

    !----------------------------------------
    ! terme source de pluie
    ! Spluie = 0 -> pas de pluie
    SELECT CASE (Spluie) 

    CASE(1) ! pluie constante
       CALL Pluie_constante(pluie)

    END SELECT

    Source(1) = Source(1) + pluie

  END SUBROUTINE TermeSource

  !*******************************************************************!

  !************************************************************************!
  SUBROUTINE Pluie_constante(P)

    REAL(PR), INTENT(out)   :: P
    
    P = 3._PR

  END SUBROUTINE Pluie_constante
  !**************************************************************************!

  !***************************************************************************!
  SUBROUTINE Friction_centree(UL,UR,fx,fy)
    
    REAL(PR), DIMENSION(:), INTENT(in)  :: UL, UR
    REAL(PR),               INTENT(out) :: fx, fy
    
    REAL(PR) :: hL, uxL, uyL
    REAL(PR) :: hR, uxR, uyR
    REAL(PR) :: norme_uL, norme_uR, norme_u
    REAL(PR) :: ux, uy
    REAL(PR) :: coeff_friction

    coeff_friction = 1._PR

    hL  = UL(1)
    uxL = UL(2)
    uyL = UL(3)

    hR  = UR(1)
    uxR = UR(2)
    uyR = UR(3)

    norme_uL = sqrt(uxL*uxL + uyL*uyL)
    norme_uR = sqrt(uxR*uxR + uyR*uyR)
    norme_u  = 0.5_PR*(norme_uL+norme_uR)

    ux = 0.5_PR * (uxL + uxR)
    uy = 0.5_PR * (uyL + uyR)

    fx = -coeff_friction * norme_u * ux
    fy = -coeff_friction * norme_u * uy

  END SUBROUTINE Friction_centree
!**********************************************************************!

  !***************************************************************************!
  SUBROUTINE Source_topo_centree(UL,UR,ZL,ZR,dx,Source_topo)
    
    REAL(PR), DIMENSION(:), INTENT(in)  :: UL, UR
    REAL(PR),               INTENT(in)  :: ZL, ZR, dx
    REAL(PR),               INTENT(out) :: Source_topo
    
    REAL(PR) :: hL, uxL, uyL
    REAL(PR) :: hR, uxR, uyR

    hL  = UL(1)
    uxL = UL(2)
    uyL = UL(3)

    hR  = UR(1)
    uxR = UR(2)
    uyR = UR(3)

    Source_topo = -gravitation*0.5_PR*(hL+hR)*(ZR-ZL)/dx

  END SUBROUTINE Source_topo_centree
!**********************************************************************!


END MODULE Source
