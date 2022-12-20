MODULE Un_Pas_Temps

  !!$ Module qui gere les iterations en temps !
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  !!$ listes des subroutines : 
  !!$ - boucle_en_temps 
  !!$ - affiche_resu
  !!$ - pas_de_temps

  USE readmesh
  USE Var_Types
  USE IO
  USE Topographie
  USE Flux_Numeriques
  USE Source
  USE Conditions_limites

  IMPLICIT NONE

CONTAINS

!***********************************************************************!
  SUBROUTINE boucle_en_temps(Mesh, Seg, Var, DATA, vtkName)
    !---------------------------------------------
    !     role:   
    !

    TYPE(MeshDef),INTENT(inout)      :: Mesh
    TYPE(MeshSeg),INTENT(in)         :: Seg
    TYPE(Variables),INTENT(inout)    :: Var
    TYPE(Donnees),INTENT(inout)      :: DATA

    CHARACTER(LEN=70),INTENT(inout)  :: vtkName
    CHARACTER(LEN=70)                :: vtkName_topo, vtkname_aux

    INTEGER  :: Nsegmtot, Ncells, Nvar, compt_segfr
    INTEGER  :: is, is1, is2, iseg, file_count
    INTEGER  :: NbVois1, NbVois2
    INTEGER  :: Impre, icl
    INTEGER  :: kt, kt0
    REAL(PR) :: Temps, CFL, Dt, rdenom1, rdenom2, nx, ny
    REAL(PR) :: Z1, Z2, dx1, dx2

    INTEGER,  DIMENSION(:,:,:), POINTER          :: TabSeg
    REAL(PR), DIMENSION(1:DATA%Nvar)             :: Ua1, Ua2, Umin, Umax

    REAL(PR), DIMENSION(:),   ALLOCATABLE        :: Flux
    REAL(PR), DIMENSION(:),   ALLOCATABLE        :: Source12, Source21, Source11, Source22


    OPEN( UNIT=16, FILE='results/evolve.time', FORM='formatted')
    REWIND(16)

    Nsegmtot    = Seg%Nsegmtot
    Ncells      = Var%Ncells
    Nvar        = Data%Nvar
    Impre       = Data%Impre
    Temps       = DATA%T
    kt0         = DATA%KT
    CFL         = Data%Cfl
    file_count  = (kt0/Data%ifres)+1

    ALLOCATE(Flux( 1:Nvar    ))   ! flux numerique aux interfaces!
    ALLOCATE(Source12( 1:Nvar))   ! source aux interfaces entre i1 et i2
    ALLOCATE(Source21( 1:Nvar))   ! source aux interfaces entre i2 et i1
    ALLOCATE(Source11( 1:Nvar))   ! source sur la cellule i1
    ALLOCATE(Source22( 1:Nvar))   ! source sur la cellule i2

    TabSeg => Seg%TabSeg
 
    IF( Impre > 2 ) THEN
       WRITE(6,*)
       WRITE(6,*)  "::::::::::::::::::::::::::::::::::::::::::::::::::::::"
       WRITE(6,*)  "::   Debut de la boucle en temps :: Maillage Fixe   ::"
       WRITE(6,*)  "::::::::::::::::::::::::::::::::::::::::::::::::::::::"
       WRITE(6,*)
    END IF

    !----------------------------------------------------!
    !            DEBUT DE LA BOUCLE EN TEMPS             !
    !----------------------------------------------------!

    DO kt = kt0+1, DATA%Ktmax

       compt_segfr = 0
       
       ! verification de la positivite de la hauteur d'eau !

       DO is = 1,Ncells
          IF (Var%Ua(1,is) < 0._PR) THEN
             WRITE(6,*) 'Erreur verif : hauteur d''eau negative',var%Ua(1,is),' a l indice is = ',is
             WRITE(6,*) 'Arret a l iteration', kt,' au temps t:',temps
             WRITE(6,*)
             WRITE(6,*) 'ARRET DU CODE Un_pas_remps.f90/boucle_en_temps'
             STOP
          END IF

       END DO

       !----------------------------------------------------!
       !    CALCUL DU FLUX NUMERIQUE AUX INTERFACES         !
       !----------------------------------------------------!

       CALL PasDeTemps(Seg, Var, Data%Dt, CFL)
       Dt = MIN(Data%Dt,  Data%Tmax - Temps)

       DO is = 1, Ncells
          Var%Un(1,is) = Var%Ua(1,is)
          Var%Un(2,is) = Var%Ua(2,is)
          Var%Un(3,is) = Var%Ua(3,is)
       ENDDO

        ! Boucle sur les segments ! 

       DO iseg = 1,Nsegmtot

          Flux     = 0.0_PR
          Source12 = 0.0_PR
          Source21 = 0.0_PR
          source11 = 0._PR
          source22 = 0._PR
        
          is1 = TabSeg(1,iseg,1)
          is2 = TabSeg(2,iseg,1)

          rdenom1 = Seg%LenSeg(iseg) / Var%CellVol(is1)
          rdenom2 = Seg%LenSeg(iseg) / Var%CellVol(is2)

          nx = Seg%Normext(1,iseg)
          ny = Seg%Normext(2,iseg)

          NbVois1 = Var%NbVois(is1)
          NbVois2 = Var%NbVois(is2)

          Z1 = DATA%Z(is1)
          Z2 = DATA%Z(is2)

          dx1 = 1._PR/rdenom1/NbVois1
          dx2 = 1._PR/rdenom2/NbVois2
   
          ! cas ou le segment est a l'interieur !
          IF(TabSeg(1,iseg,2) .eq. 0) THEN 
            
             Ua1 = Var%Ua(:,is1) 
             Ua2 = Var%Ua(:,is2) 
             
             CALL Calcul_Flux_numerique (DATA%Iflux, nx, ny, Ua1, Ua2, Flux)

             CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &     nx,  ny, Ua1, Ua2, Z1, Z2, dx1, Source12, Mesh%coor(1:2,is1), Mesh%coor(1:2,is2))
            
              CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &     nx,  ny, Ua1, Ua1, Z1, Z1, dx1, Source11, Mesh%coor(1:2,is1), Mesh%coor(1:2,is1))

             CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &    -nx, -ny, Ua1, Ua2, Z2, Z1, dx2, Source21, Mesh%coor(1:2,is2), Mesh%coor(1:2,is1))

             CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &    -nx, -ny, Ua2, Ua2, Z2, Z2, dx2, Source22, Mesh%coor(1:2,is2), Mesh%coor(1:2,is2))

             !----------------------------------------------------!
             !      SCHEMA NUMERIQUE (interieur du domaine)       !
             !----------------------------------------------------!
             Var%Un(:,is1) = Var%Un(:,is1) - Dt*rdenom1 * Flux(:) + &
                  &          0.5_PR*Dt/NbVois1 * (Source12(:)+Source11(:))

             Var%Un(:,is2) = Var%Un(:,is2) + Dt*rdenom2 * Flux(:) + &
                  &          0.5_PR*Dt/NbVois2 * (Source21(:)+Source22(:))

          ! cas ou le segment est sur la frontiere !
          ELSE
            
             compt_segfr = compt_segfr + 1
             icl = TabSeg(2,iseg,2)

             Ua1 = Var%Ua(:,is1)
             CALL CL_evaluation_UA2 (Ua1, Ua2, icl)

             CALL Calcul_Flux_numerique (DATA%Iflux, nx, ny, Ua1, Ua2, Flux)

             CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &     nx,  ny, Ua1, Ua2, Z1, Z2, dx1, Source12, Mesh%coor(1:2,is1), Mesh%coor(1:2,is2))
             
             CALL TermeSource (Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie, &
                  &     nx,  ny, Ua1, Ua1, Z1, Z1, dx1, Source11, Mesh%coor(1:2,is1), Mesh%coor(1:2,is1))

             !----------------------------------------------------!
             !         SCHEMA NUMERIQUE (bord du domaine)         !
             !----------------------------------------------------!

             Var%Un(:,is1) = Var%Un(:,is1) - Dt*rdenom1*Flux(:) + &
                  &          0.5_PR*Dt/NbVois1 * (Source12(:)+Source11(:))
           
          ENDIF
                             
       ENDDO

       !-----------------------------------------------------------!
       !                 MISE A JOUR de U avec le terme source     !
       !-----------------------------------------------------------!
       Var%Ua = Var%Un

       DO is = 1, Ncells

          ! verification de la positivite de la hauteur d'eau
          IF (Var%Ua(1,is) < 0._PR) THEN
             WRITE(6,*) 'hauteur <0', Var%Ua(1,is)
             WRITE(6,*) 'arret dans la mise a jour : Un_Pas_Tremps.f90'
             STOP
          END IF

          ! si la hauteur d'eau est sous la precision machine on impose hauteur = 0
          IF (Var%Ua(1,is) < eps) Var%Ua(1,is) = 0._PR
       ENDDO

       !----------------------------------------------------!
       !         AFFICHAGE ET ENREGISTREMENT                !
       !----------------------------------------------------!
    
       Temps       = Temps + Dt
       DATA%T      = Temps
       DATA%kt     = kt

       IF (Impre >2  .AND. MOD(kt,Data%ifre) .EQ. 0 ) THEN
          CALL afficheResu(6)
          CALL afficheResu(16)
       END IF
   
       IF( ABS(Temps - Data%Tmax) .LT. eps) EXIT

       IF (MOD(kt,Data%ifres) .EQ. 0 ) THEN
      
          vtkName_topo="surface"
          CALL rename_Box4SW(file_count, 4, VtkName_aux, vtkname_topo)
          CALL CellVertexvtk(DATA, Mesh, Var, vtkName_aux, 0 )

          if (Data%Itopo .eq. 2) then
             vtkName_topo="topo"
             CALL rename_Box4SW(file_count, 4, VtkName_aux, vtkname_topo)
             CALL CellVertexvtk(DATA, Mesh, Var, vtkName_aux, 1 )
          end if

          file_count = file_count + 1

          CALL Reprise_out(DATA, Mesh, Var, vtkName)

       END IF

    END DO
    ! fin de la boucle en temps !
       
    IF( Impre > 2 ) WRITE(6,*) " Fin de la boucle en temps "
    CALL afficheResu(6)
    CALL afficheResu(16)

    !-----------------------------------------!
    !   ENREGISTREMENT DES RESULTATS FINAUX   !
    !-----------------------------------------!

    vtkName_topo="surface"
    CALL rename_Box4SW(file_count, 4, VtkName_aux, vtkname_topo)
    CALL  CellVertexvtk(DATA, Mesh, Var, vtkName_aux, 0 )

    vtkName_topo="topo"
    CALL CellVertexvtk(DATA, Mesh, Var, vtkname_topo, 1 )


  CONTAINS

    !********************************!
    SUBROUTINE  afficheResu(argunit)

      INTEGER, INTENT(in) :: argunit
      
199   FORMAT(":::::::::::::::::::::::::::::::::::::::::::::::::&
           &:::::::::::::::::::::::::::::::")
189   FORMAT("::  ----------------------------&
           &------------------                            ::" )
187   FORMAT("::                                               &
           &                             ::")

      WRITE(argunit,*)
      WRITE(argunit,199)
      WRITE(argunit,189)
      WRITE(argunit,'(":: ",&
           &"|  Temps = ", E10.3, " ||      ",&
           &"  kt = " ,I7,"   |                           ::"  )' ) Temps, kt
      WRITE(argunit,'(":: |  Dt =  ", E10.3,1x, &
           &   "  ||    Cfl = ", E10.3, "   |&
           &                           ::"  )' ) Dt , Dt/Data%Dt
      WRITE(argunit,189)

      WRITE(argunit,187)

      Umin = MINVAL( Var%Ua, DIM=2 )
      Umax = MAXVAL( Var%Ua, DIM=2 )
      IF(kt==1)  WRITE(6,'(":: Ua  min     =  ", 4(3x,E12.5)," Me =  ", I5)') Umin
      IF(kt==1)  WRITE(6,'(":: Ua  max     =  ", 4(3x,E12.5)," Me =  ", I5)') Umax

      WRITE(argunit,'("::  Ua  min   =", 3(3x,E12.5), &
               &"                  ::" )') Umin
      WRITE(argunit,'("::  Ua  max   =", 3(3x,E12.5), &
               &"                  ::" )') Umax
      WRITE(argunit,199) 
      WRITE(argunit,*) 


    END SUBROUTINE afficheResu
    !**********************************!
    

  END SUBROUTINE boucle_en_temps
  !***********************************************************************!


  !***********************************************************************!
  SUBROUTINE PasDeTemps(Seg, Var, Dt, CFL)
   
    TYPE(Variables)              :: Var
    TYPE(MeshSeg)   ,INTENT(IN)  :: Seg
    REAL(PR), INTENT(OUT)        :: Dt
    REAL(PR), INTENT(IN)         :: CFL

    INTEGER        :: is1, is2, iseg
    REAL(PR)       :: haut, minhaut
    REAL(PR)       :: nx, ny, maxVP

    maxVP = 0._PR
    
    do iseg = 1,Seg%Nsegmt

       nx = Seg%Normext(1,iseg)
       ny = Seg%Normext(2,iseg)

       is1 = Seg%Nubo(1,iseg)
       is2 = Seg%Nubo(2,iseg)

       maxVP = max( maxVP, ValeurPropre(Var%Ua(:,is1), nx, ny) )
       maxVP = max( maxVP, ValeurPropre(Var%Ua(:,is2), nx, ny) )
       
    end do

    ! calcul de la hauteur de triangle la plus petite
    minhaut = 1.E10_PR

    DO iseg = 1, Seg%Nsegmtot

       if(Seg%TabSeg(1,iseg,2) == 0) then
          is1    = Seg%TabSeg(1,iseg,1)
          is2    = Seg%TabSeg(2,iseg,1)
          
          haut    = Var%CellVol(is1) / Seg%LenSeg(iseg)
          minhaut = MIN(minhaut,haut)
          
          haut    = Var%CellVol(is2) / Seg%LenSeg(iseg)
          minhaut = MIN(minhaut,haut)
       endif

    END DO

    Dt = CFL*minhaut/maxVP

  END SUBROUTINE PasDeTemps
  !***********************************************************************!


END MODULE Un_Pas_Temps
