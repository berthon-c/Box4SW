MODULE Initialisation

!!$ Module gerant les entrees/sorties et l'initialisation

  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  ! listes des subroutines :
  ! - Tube_a_chocs_1D_En_x
  ! - InitRiemann
  ! - Solution_Constante
  ! - DonNumMeth
  ! - DonClimites
  ! - Reprise_in
  ! - Reprise_out
  ! - CellVertexVtk


  USE var_types
  USE Physique

  IMPLICIT NONE

CONTAINS

  !******************************************************!
  !*    Choix de la donnee initiale                     *!
  !******************************************************!
  SUBROUTINE Donnees_initiales(DATA, Mesh, Var)
     
    TYPE(Donnees),     INTENT(in)       :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    SELECT CASE(DATA%Icas)
    CASE(0) ! donnee initiale constante
       CALL Tube_a_chocs_1D_En_x(DATA, Mesh, Var)
       
    CASE(1) ! initialisation d'un pb de riemann
       CALL Tube_a_chocs_1D_En_x(DATA, Mesh, Var)
       
    CASE(2) ! initialisation rupture de barrage sur fond sec
       CALL Tube_a_chocs_1D_En_x(DATA, Mesh, Var)
       
    CASE DEFAULT
       print*,"initialisation demandee non codee Icas=",DATA%Icas
       STOP
    END SELECT

  END SUBROUTINE Donnees_Initiales

  !******************************************************!
  !*    Tube A Choc 1D                                  *!
  !******************************************************!
  SUBROUTINE Tube_a_chocs_1D_En_x(DATA, Mesh, Var)

    TYPE(Donnees),     INTENT(in)       :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    INTEGER            :: is
    REAL(PR)           :: r
    REAL(PR)           :: Xmin, Xmax, Xmil
    REAL(PR)           :: hL, uxL, uyL
    REAL(PR)           :: hR, uxR, uyR


    IF (Data%impre > 5 ) WRITE(6,*) ' ..... Initialisation  &
         &    Tube a chocs 1D Suivant x .... '

    SELECT CASE(DATA%Icas)
    CASE(0)
       CALL Init_Constante(hL,uxL,uyL)
       hR  = hL
       uxR = uxL
       uyR = uyL
       
    CASE(1)
       CALL Init_Riemann(hL,uxL,uyL,hR,uxR,uyR)

    CASE(2)
       CALL Init_barrage_fond_sec(hL,uxL,uyL,hR,uxR,uyR)

     CASE DEFAULT
        print*,"initialisation UL et UR non codee Icas=",DATA%Icas
        STOP
    END SELECT

    DO is = 1, Var%NCells
       Var%UA(1,is) = max(0._PR, hL-DATA%Z(is))
       Var%UA(2,is) = Var%UA(1,is)*uxL
       Var%UA(3,is) = Var%UA(1,is)*uyL
    END DO

    Xmax = 0.0_PR
    Xmin = 0.0_PR

    DO is = 1, Var%Ncells
       Xmax = Max(Xmax, Mesh%coor(1,is) )
       Xmin = Min(Xmin, Mesh%coor(1,is) )
    END DO
    Xmil = 0.5_PR*(Xmin + Xmax)
    Xmil = 0.70_PR

    DO is = 1, Var%NCells
       r = Mesh%coor(1,is)
     
       IF( r > Xmil ) THEN
          Var%UA(1,is) = max(0._PR, hR-DATA%Z(is))
          Var%UA(2,is) = Var%UA(1,is)*uxR
          Var%UA(3,is) = Var%UA(1,is)*uyR
       END IF
    END DO

    IF (Data%impre > 5 ) WRITE(6,*) ' UA  min =  ' , MINVAL( Var%UA, DIM=2 )
    IF (Data%impre > 5 ) WRITE(6,*) ' UA  max =  ' , MAXVAL( Var%UA, DIM=2 )
    IF (Data%impre > 5 ) WRITE(6,*) ' ..... Fin de l''initialisation .... '


  END SUBROUTINE Tube_a_chocs_1D_En_x
  !************************************************************************!


  !************************************************************************!
  SUBROUTINE Init_Riemann(hL,uxL,uyL,hR,uxR,uyR)

    REAL(PR), INTENT(out)   :: hL,uxL,uyL
    REAL(PR), INTENT(out)   :: hR,uxR,uyR

    hL  = 2.0_PR
    uxL = 0.0_PR
    uyL = 0.0_PR
    
    hR  = 1.0_PR
    uxR = 0.0_PR
    uyR = 0.0_PR

  END SUBROUTINE Init_Riemann
  !************************************************************************!
  
  !************************************************************************!
  SUBROUTINE Init_barrage_fond_sec(hL,uxL,uyL,hR,uxR,uyR)

    REAL(PR), INTENT(out)   :: hL,uxL,uyL
    REAL(PR), INTENT(out)   :: hR,uxR,uyR

    hL  = 2.0_PR
    uxL = 0.0_PR
    uyL = 0.0_PR
    
    hR  = 0.0_PR
    uxR = 0.0_PR
    uyR = 0.0_PR

  END SUBROUTINE Init_Barrage_Fond_Sec
  !************************************************************************!

  !************************************************************************!
  SUBROUTINE Init_Constante(h,ux,uy)

    REAL(PR), INTENT(out)   :: h,ux,uy
    
    h  = 2.0_PR
    ux = 0.0_PR
    uy = 0.0_PR

  END SUBROUTINE Init_Constante
  !**************************************************************************!

END MODULE Initialisation

