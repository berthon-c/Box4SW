MODULE Topographie

!!$ Module gerant la topographie

  USE var_types

  IMPLICIT NONE

CONTAINS


  !******************************************************!
  !*   initialisation de la topographie                 *!
  !******************************************************!
  SUBROUTINE initialisation_topo(DATA, Mesh, Var)

    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    IF (Data%impre > 5 ) WRITE(6,*) ' ..... Initialisation  de la topographie .... '

    SELECT CASE(data%Itopo)
       
       case(1)
          call topo_rectangle_bosse(DATA, Mesh, Var)

       case(3)
          call topo_parabolic(DATA, Mesh, Var)

       case(4)
          call topo_canal_bosse(DATA, Mesh, Var)

       case default
          call topo_plate(DATA)

    END SELECT

  END SUBROUTINE initialisation_topo

  !************************************************************************!

  SUBROUTINE topo_plate(DATA)

    TYPE(Donnees),     INTENT(inout)    :: DATA

    DATA%Z = 0._PR    ! Z0 de reference

  END SUBROUTINE topo_plate
  !************************************************************************!

  !************************************************************************!

  SUBROUTINE topo_rectangle_bosse(DATA, Mesh, Var)

    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    INTEGER            :: is
    REAL(PR)           :: x, y, r, r2
    REAL(PR)           :: Xmin, Xmax, Xmil
    REAL(PR)           :: Ymin, Ymax, Ymil

    DATA%Z = 0.1_PR    ! Z0 de reference

    Xmax = 0.0_PR
    Xmin = 0.0_PR

    DO is = 1, Var%Ncells
       Xmax = Max(Xmax, Mesh%coor(1,is) )
       Xmin = Min(Xmin, Mesh%coor(1,is) )
    END DO
    Xmil = 0.25_PR*Xmin + (1._PR-0.25_PR)*Xmax

    Ymax = 0.0_PR
    Ymin = 0.0_PR

    DO is = 1, Var%Ncells
       Ymax = Max(Ymax, Mesh%coor(2,is) )
       Ymin = Min(Ymin, Mesh%coor(2,is) )
    END DO
    Ymil = 0.5_PR*(Ymin + Ymax)

    DO is = 1, Var%NCells
       x = Mesh%coor(1,is)
       y = Mesh%coor(2,is)
       r = min(Xmax-Xmin, Ymax-Ymin)/5._PR
    
       r2 = (x-Xmil)**2+(y-Ymil)**2
       IF( r2 < r**2 ) THEN
          DATA%Z(is) = DATA%Z(is) + 1.5_PR*(r**2-r2)/r**2
       END IF
    END DO

    IF (Data%impre > 1 ) WRITE(6,*) ' Z  min =  ' , MINVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' Z  max =  ' , MAXVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' ..... Fin de l''initialisation topo.... '

  END SUBROUTINE topo_rectangle_bosse
  !************************************************************************!

  !************************************************************************!

  SUBROUTINE topo_canal_bosse(DATA, Mesh, Var)

    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    INTEGER            :: is
    REAL(PR)           :: x, y, r, r2, ht
    REAL(PR)           :: Xmil, Ymil

    DATA%Z = 0.1_PR    ! Z0 de reference

    Xmil = 30._PR
    Ymil = 20._PR
    r    = 5._PR
    ht   = 1.5_PR

    DO is = 1, Var%NCells
       x = Mesh%coor(1,is)
       y = Mesh%coor(2,is)
 
       r2 = (x-Xmil)**2+(y-Ymil)**2
       IF( r2 < r**2 ) THEN
          DATA%Z(is) = DATA%Z(is) + ht*(r**2-r2)/r**2
       END IF
    END DO

    IF (Data%impre > 1 ) WRITE(6,*) ' Z  min =  ' , MINVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' Z  max =  ' , MAXVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' ..... Fin de l''initialisation topo.... '

  END SUBROUTINE topo_canal_bosse
  !************************************************************************!

  !************************************************************************!

  SUBROUTINE topo_parabolic(DATA, Mesh, Var)

    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var

    INTEGER            :: is
    REAL(PR)           :: x, y
    REAL(PR)           :: Xmin, Xmax, Xmil
    REAL(PR)           :: Ymin, Ymax, Ymil

    DATA%Z = 0.1_PR    ! Z0 de reference

    Xmax = 0.0_PR
    Xmin = 0.0_PR

    DO is = 1, Var%Ncells
       Xmax = Max(Xmax, Mesh%coor(1,is) )
       Xmin = Min(Xmin, Mesh%coor(1,is) )
    END DO
    Xmil = 0.5_PR*(Xmin + Xmax)

    Ymax = 0.0_PR
    Ymin = 0.0_PR

    DO is = 1, Var%Ncells
       Ymax = Max(Ymax, Mesh%coor(2,is) )
       Ymin = Min(Ymin, Mesh%coor(2,is) )
    END DO
    Ymil = 0.5_PR*(Ymin + Ymax)

    DO is = 1, Var%NCells
       x = Mesh%coor(1,is)
       y = Mesh%coor(2,is)

       DATA%Z(is) = DATA%Z(is) + 0.1_PR*(1._PR - y)**2      
       
    END DO

    IF (Data%impre > 1 ) WRITE(6,*) ' Z  min =  ' , MINVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' Z  max =  ' , MAXVAL( DATA%Z )
    IF (Data%impre > 1 ) WRITE(6,*) ' ..... Fin de l''initialisation topo.... '

  END SUBROUTINE topo_parabolic
  !************************************************************************!


END MODULE Topographie
