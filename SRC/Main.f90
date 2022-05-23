PROGRAM Main
 
  ! Code qui résoud le système M1 du Transfert Radiatif avec des méthodes !
  ! de volumes finis en 2D sur un maillage non structuré !

  !! ************************************************************ !!
  !! Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !!
  !! ************************************************************ !!

  USE Var_Types
  USE IO 
  USE Initialisation 
  USE Topographie 
  USE ReadMesh 
  USE SegmtsGeom
  USE ParmGeom 
  USE Physique 
  USE Flux_Numeriques
  USE Un_Pas_Temps 

  IMPLICIT NONE

  TYPE(MeshDef)          :: Mesh
  TYPE(MeshSeg)          :: Seg
  TYPE(Variables)        :: Var
  TYPE(Donnees)          :: DATA

  CHARACTER(LEN=70)      :: BaseMesh, MeshFile, vtkName, MeshType

  INTEGER                :: Nvar
  INTEGER                :: lensd
  REAL(PR), DIMENSION(2) :: Xmin, Xmax


  ! -----------------------------------------------------------------------! 
  !     CREATION DU REPERTOIRE DES RESULTATS S'IL N'EXISTE PAS ENCORE      !
  ! -----------------------------------------------------------------------!
  CALL system("mkdir results")
  
  ! -----------------------------------------! 
  !     LECTURE DES DONNEES DU PROBLEME      !
  ! -----------------------------------------!
  CALL DonNumMeth(DATA)
  MeshType     = "MESH"
  Nvar         = Data%Nvar
 
  WRITE(6,*) 
  WRITE(6,*) " ______________________________________________________________________"
  WRITE(6,*) " |                              Box 4 Shallow water                   |"
  WRITE(6,*) " ----------------------------------------------------------------------"
  WRITE(6,*)
  WRITE(6,*) "  ********************************************************************* "
  WRITE(6,*) "  **************         Cas TEST ",Data%PbName(1:20),"   ************** " 
  WRITE(6,*) "  ********************************************************************* "
  WRITE(6,*) 

 
  ! -----------------------------------------! 
  !     LECTURE ET CREATION DU MAILLAGE      !
  ! -----------------------------------------! 
  BaseMesh = Data%PbName
  MeshFile=" "
  vtkName = Data%PbName
  IF( Data%Impre > 2 ) WRITE(6, *) ' Nom local du probleme :: ',vtkName 

  SELECT CASE( TRIM(MeshType) )
  CASE("MESH")

     lensd                     = INDEX(BaseMesh,' ') - 1
     BaseMesh(lensd+1:lensd+7) = ".mesh " 
     MeshFile = BaseMesh
     IF( Data%Impre > 2 )  WRITE(6,*) "Fichier de maillage MeshFile : ",MeshFile
     CALL ReadMESH2D( Mesh, MeshFile )

  CASE DEFAULT
     WRITE(6, *) " Maillage de type ",TRIM(MeshType) ," non traite"
     STOP
  END SELECT

  Xmin = MINVAL( Mesh%coor, DIM=2)
  Xmax = MAXVAL( Mesh%coor, DIM=2)

  ! -----------------------------------------! 
  !   CREATION DE LA GEOMETRIE(SEGMENTS)     !
  ! -----------------------------------------! 
  Var%NCells  =  Mesh%Npoint

  CALL Segments2DP( Mesh, Seg, DATA, Var)

  ! ---------------------------------------------------! 
  !   LECTURE ET CREATION DES CONDITIONS AUX LIMITES   !
  ! ---------------------------------------------------!

  CALL DonClimites( DATA , Mesh, Seg)

  ! -----------------------------------------! 
  !        TOPOGRAPHIE                       !
  ! -----------------------------------------! 
  ALLOCATE( DATA%Z(1:Var%NCells) )
  CALL initialisation_topo(DATA, Mesh, Var)

  ! -----------------------------------------! 
  !        INITIALISATION OU REPRISE         !
  ! -----------------------------------------! 
  ALLOCATE(  Var%Ua(   1:Nvar , 1:Var%NCells) )     ! inconnues dans les cellules!
  ALLOCATE(  Var%Un(   1:Nvar , 1:Var%NCells) )     ! inconnues temporaires !

  IF( DATA%Reprise ) THEN
     CALL Reprise_in(DATA, Mesh, Var, vtkName)
     IF( DATA%Icas >= 30 .AND. DATA%Icas < 40 ) THEN
        DATA%KT = 0
        DATA%T  = 0.0_PR
     END IF

  ELSE

     CALL Donnees_initiales(DATA, Mesh, Var)
     
     Data%T  = 0.0_PR
     Data%KT = 0

  END IF
  
  ! -----------------------------------------! 
  !      CALCUL DES AIRES DES CELLULES,      !
  !    DES NORMALES ET LONGUEURS DES SEG     !
  ! -----------------------------------------! 

  CALL CellVertexMvvno(Mesh, Seg)
  CALL GeomAires(Mesh)

  NULLIFY( Var%CellVol )
  Var%CellVol => Mesh%Airs

  ! -----------------------------------------! 
  !   ENREGISTREMENT DES DONNEES INITIALES   !
  ! -----------------------------------------!
 
  CALL  CellVertexvtk(DATA, Mesh, Var, vtkName, 0 )

  ! -----------------------------------------! 
  !           ITERATIONS EN TEMPS            !
  !      + ENREGISTREMENT DES RESULTATS      !
  ! -----------------------------------------!

  CALL boucle_en_temps(Mesh, Seg,  Var, DATA, vtkName)
  
  ! ---------------------------------------------------!
  !       CREATION DU FICHIER DE REPRISE               !
  ! ---------------------------------------------------!

  CALL Reprise_out(DATA, Mesh, Var, vtkName)

  ! ---------------------------------------------------!
  !            LIBERATION DE LA MEMOIRE                !
  ! ---------------------------------------------------!

  DEALLOCATE ( Data%lgn, Data%lgf1, Data%lgf2, Data%lgc, Data%lgd, Data%Fr )
  DEALLOCATE ( Mesh%Ndegre,Mesh%Typoint,Mesh%Nu,Mesh%Coor,Mesh%Airt,Mesh%Airs)
  DEALLOCATE ( Seg%TabSeg, Seg%Vno, Seg%LenSeg, Seg%NbVois , Seg%Numsegfr)
  DEALLOCATE ( Seg%LogFac, Seg%Nubo, Seg%Nusv, Seg%Nutv, Seg%Normext)
  DEALLOCATE ( Var%Ua, Var%Un, Var%NbVois, Var%NuVois )

  print*, 'Fin d''execution de Box4SW'
  print*,


END PROGRAM Main
