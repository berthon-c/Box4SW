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
                  &    nx, ny, Ua1, Ua2, Z1, Z2, dx, Source,coor1,coor2)

    INTEGER,                 INTENT(in)    :: Stopo, Sfriction, Spluie
    REAL(PR),                INTENT(in)    :: nx, ny
    REAL(PR), DIMENSION(3),  INTENT(in)    :: Ua1, Ua2
    REAL(PR),                INTENT(in)    :: Z1, Z2, dx
    REAL(PR), DIMENSION(3),  INTENT(out)   :: Source
    REAL(PR), DIMENSION(2),  INTENT(in)    :: coor1, coor2

    REAL(PR), DIMENSION(3)                 :: UL, UR

    REAL(PR)    :: pluie, xp, yp
    REAL(PR)    :: friction_n, friction_t
    REAL(PR)    :: source_topo

    Source      = 0._PR
    pluie       = 0._PR
    friction_n  = 0._PR
    friction_t  = 0._PR
    source_topo = 0._PR

    xp =0.5_PR*(coor1(1)+coor2(1))
    yp =0.5_PR*(coor1(2)+coor2(2))
    
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
    CASE(2) ! source topo Meissa Eq.2.26
       CALL Source_topo_M(UL,UR,Z1,Z2,dx,Source_topo)
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
    CASE(2) ! pluie intereieur
       CALL Pluie_interieur(pluie,xp,yp)

    END SELECT

    Source(1) = Source(1) + pluie

  END SUBROUTINE TermeSource

  !*******************************************************************!

  !************************************************************************!
  SUBROUTINE Pluie_constante(P)

    REAL(PR), INTENT(out)   :: P
    
    P = 0.1_PR

  END SUBROUTINE Pluie_constante
  !**************************************************************************!


    !************************************************************************!
  SUBROUTINE Pluie_interieur(P,x,y)

    REAL(PR), INTENT(out)   :: P
    REAL(PR), INTENT(in)    :: x, y
    
    P = 0.0_PR
    if ( (x>0.1_PR).AND.(x<1.9_PR).AND.(y>0.1_PR).AND.(y<0.9_PR) ) THEN
       P = 0.1_PR
    END if
 
  END SUBROUTINE Pluie_interieur
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



   !***************************************************************************!
  SUBROUTINE Source_topo_M(UL,UR,ZL,ZR,dx,Source_topo)
    
    REAL(PR), DIMENSION(:), INTENT(in)  :: UL, UR
    REAL(PR),               INTENT(in)  :: ZL, ZR, dx
    REAL(PR),               INTENT(out) :: Source_topo
    
    REAL(PR) :: hL, uxL, uyL, qxL, BL
    REAL(PR) :: hR, uxR, uyR, qxR, BR
    REAL(PR) :: q_tilde, Froude, h_bar, qL, qR, k, p, norme_uL, norme_uR, norme_u
    
    hL  = UL(1)
    uxL = UL(2)
    uyL = UL(3)

    hR  = UR(1)
    uxR = UR(2)
    uyR = UR(3)

    ! variable debit 
    qxL = hL*uxL
    qxR = hR*uxR

    p = 0.5_PR
    h_bar = (hL+hR)*0.5_PR

    CALL Bernouilli(UL,UR, ZL, ZR, BL, BR)
 !   Br = compute_bernouilli(WR)  
 !   Bl = compute_bernouilli(WL) 
    
    norme_uL = sqrt(uxL*uxL + uyL*uyL)
    norme_uR = sqrt(uxR*uxR + uyR*uyR)
    norme_u  = 0.5_PR*(norme_uL+norme_uR)

    qL = 0.5_PR * (qxL + qxR)
    qR = 0.5_PR * (qxL + qxR)
    
    q_tilde =  (qL**2 + qR**2)*0.5_PR

    Froude = ( q_tilde*h_bar )/( gravitation*(hL**2)*(hR**2) )
        
    if ((hL.lt. 1.e-10_PR).or.(hR.lt. 1.e-10_PR)) then 
       k = 0._PR
    else        
       k = sqrt(abs(BR-BL)**2 + abs(qR-qL)**2)/dx
    endif
    
    !if (Froude == 1._PR) then 
     ! print*, 4._PR*(hL**2)*(hR**2)*((1._PR - Froude)**2 + (dx**p)*k)
   !endif

    
    if( (hL.lt. 1.e-10_PR).and.(hR.lt. 1.e-10_PR) ) then 
       Source_topo = 0._PR 
    elseif((hL.lt. 1.e-10_PR).or.(hR.lt. 1.e-10_PR)) then  
       Source_topo = -gravitation*h_bar*(ZR-ZL)
    elseif ( (1._PR - Froude)**2 + (dx**p)*k < 1.e-10_PR ) then
        Source_topo = 0._PR
    else  
       Source_topo = -(gravitation*h_bar + q_tilde*(hR-hL)**2/(4._PR*(hL**2)*(hR**2)*((1._PR-Froude)**2+(dx**p)*k)))*(ZR-ZL)
       Source_topo = Source_topo/dx 
    endif

  END SUBROUTINE Source_topo_M
  !**********************************************************************!

  !----------------------------------------------------------------------------------------------!
  !         Subroutine Bernoulli
  !----------------------------------------------------------------------------------------------!

  SUBROUTINE Bernouilli(UL, UR, ZL, ZR, BL, BR)

    implicit none

    REAL(PR), DIMENSION(:), INTENT(in)  :: UL, UR
    REAL(PR),               INTENT(in)  :: ZL, ZR
    REAL(PR),               INTENT(out) :: BL, BR
      
    REAL(PR) :: hL, uxL
    REAL(PR) :: hR, uxR
    
    hL  = UL(1)
    uxL = UL(2)
    !uyL = UL(3)

    hR  = UR(1)
    uxR = UR(2)
    !uyR = UR(3)
    
    BL = 0.5_PR*uxL**2+gravitation*(hL+ZL) 
    BR = 0.5_PR*uxR**2+gravitation*(hR+ZR) 

  END SUBROUTINE Bernouilli

    

END MODULE Source
