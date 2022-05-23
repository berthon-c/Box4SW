MODULE Flux_Numeriques

  !!$ Module qui gere les flux numeriques  !
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  !!$ listes des subroutines : 
  !!$ - GodunovTypeFlux
  !!$ - fluxhll

  USE var_types
  USE Physique

  IMPLICIT NONE

CONTAINS


!*******************************************************************!
  SUBROUTINE Calcul_Flux_numerique (Iflux, nx, ny, Ua1, Ua2, Flux)
    
    INTEGER,                 INTENT(in)    :: Iflux
    REAL(PR),                INTENT(in)    :: nx, ny
    REAL(PR), DIMENSION(3),  INTENT(in)    :: Ua1, Ua2
    REAL(PR), DIMENSION(3),  INTENT(out)   :: Flux

    REAL(PR), DIMENSION(3)                 :: UL, UR, Flux_interface

    !calcul des variables primitives dans le repere de l'interface
    UL(1) = Ua1(1)
    UL(2) = nx*vitesse(Ua1(1),Ua1(2)) + ny*vitesse(Ua1(1),Ua1(3))
    UL(3) =-ny*vitesse(Ua1(1),Ua1(2)) + nx*vitesse(Ua1(1),Ua1(3))

    UR(1) = Ua2(1)
    UR(2) = nx*vitesse(Ua2(1),Ua2(2)) + ny*vitesse(Ua2(1),Ua2(3))
    UR(3) =-ny*vitesse(Ua2(1),Ua2(2)) + nx*vitesse(Ua2(1),Ua2(3))
   
    ! calcul du flux avec UL = Ua1 et Ur = Ua2avec des variables primitives !
    SELECT CASE (Iflux)

    CASE(1)
       CALL flux_HLL_1d(UL,UR,Flux_interface)

    CASE(2)
       CALL flux_Relaxation_1d(UL,UR,Flux_interface)

    CASE DEFAULT 
       PRINT*,'Flux non code '
       STOP
    END SELECT

    ! rotation du flux dans le repere cartesien
    flux(1) = flux_interface(1)
    flux(2) = nx*flux_interface(2) - ny*flux_interface(3)
    flux(3) = ny*flux_interface(2) + nx*flux_interface(3)

  END SUBROUTINE Calcul_Flux_numerique
!***************************************************************************!


!***************************************************************************!
  SUBROUTINE flux_HLL_1d(UL,UR,Flux)
    
    REAL(PR), DIMENSION(:),  INTENT(in) :: UL, UR
    REAL(PR), DIMENSION(:),  INTENT(out):: Flux
    
    REAL(PR) :: hL, uxL, uyL, cL, vpL, presL
    REAL(PR) :: hR, uxR, uyR, cR, vpR, presR
    REAL(PR) :: valp

    hL  = UL(1)
    uxL = UL(2)
    uyL = UL(3)
    presL  = 0.5_PR*gravitation*hL*hL
    cL  = sqrt(gravitation*hL)

    hR  = UR(1)
    uxR = UR(2)
    uyR = UR(3)

    presR  = 0.5_PR*gravitation*hR*hR

    cR  = sqrt(gravitation*hR)

    vpL  = abs(uxL) + cL
    vpR  = abs(uxR) + cR
    valp = vpL
    if (valp<vpR) valp=vpR

    Flux(1) = 0.5_PR*( hL*uxL + hR*uxR                   - valp*(hR-hL))
    Flux(2) = 0.5_PR*( hL*uxL**2+presL + hR*uxR**2+presR - valp*(hR*uxR-hL*uxL))
    Flux(3) = 0.5_PR*( hL*uyL*uxL +hR*uyR*uxR            - valp*(hR*uyR-hL*uyL))

  END SUBROUTINE flux_HLL_1d
!**********************************************************************!

!**********************************************************************!
  SUBROUTINE flux_Relaxation_1d(UL,UR,Flux)

    REAL(PR), DIMENSION(:), INTENT(in) :: UL, UR
    REAL(PR), DIMENSION(:), INTENT(out):: Flux

    REAL(PR) :: a, cL, cR
    REAL(PR) :: hL, uxL, uyL, presL
    REAL(PR) :: hR, uxR, uyR, presR
    REAL(PR) :: us, pis, tau1,tau2
    REAL(PR) :: vp1,vp2,vp3

    REAL(PR), DIMENSION(1:4) :: W, WL, WR, W1, W2

    hL    = UL(1)
    uxL   = UL(2)
    uyL   = UL(3)
    presL = 0.5_PR*gravitation*hL*hL
    cL    = sqrt(gravitation*hL)

    hR    = UR(1)
    uxR   = UR(2)
    uyR   = UR(3)
    presR = 0.5_PR*gravitation*hR*hR
    cR    = sqrt(gravitation*hR)

    flux = 0._PR

    a = calcul_a(hL,cL,hR,cR)

    us   = 0.5_PR*(uxL+uxR)+0.5_PR/a*(presL-presR)
    pis  = 0.5_PR*(presL+presR)+0.5_PR*a*(uxL-uxR)
    tau1 = (us-uxL)/a+1._PR/hL
    tau2 = (uxR-us)/a+1._PR/hR

    WL(1) = hL    ; WR(1) = hR    ; W1(1) = 1._PR/tau1 ; W2(1) = 1._PR/tau2
    WL(2) = uxL   ; WR(2) = uxR   ; W1(2) = us         ; W2(2) = us
    WL(3) = uyL   ; WR(3) = uyR   ; W1(3) = uyL        ; W2(3) = uyR
    WL(4) = presL ; WR(4) = presR ; W1(4) = pis        ; W2(4) = pis

    vp1 = uxL - a/hL
    vp2 = us
    vp3 = uxR + a/hR

    W = WL
    if ((vp1.le.0._PR).and.(vp2.gt.0._PR)) W = W1
    if ((vp2.le.0._PR).and.(vp3.gt.0._PR)) W = W2
    if (vp3.le.0._PR)                      W = WR
    
    flux(1) = W(1)*W(2)
    flux(2) = W(1)*W(2)**2+W(4)
    flux(3) = W(1)*W(2)*W(3)

  end SUBROUTINE flux_Relaxation_1d
!**********************************************************************!

!**********************************************************************!
  function calcul_a(hL,cL,hR,cR)
    
    real(PR), intent(in) :: hL,cL
    real(PR), intent(in) :: hR,cR
    real(PR)             :: calcul_a
    
    real(PR)             :: a
    
    a = max(hL*cL, hR*cR)
    
    calcul_a=a*1.0001_PR + 1.d-6
    
  end function calcul_a
!**********************************************************************!

END MODULE Flux_Numeriques
