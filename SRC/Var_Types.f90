MODULE Var_Types

  !!$ Module qui gere les variables globales et les types  !
 
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!
  

  IMPLICIT NONE

  INTEGER, PARAMETER     :: PR=KIND(1.0D0)

  REAL(PR), PARAMETER    :: gravitation = 10._PR         ! gravitation

  REAL(PR), PARAMETER    :: eps       = 1.E-13_PR        ! tolerance 
  INTEGER, PARAMETER     :: Dirichlet = 4, Neumann   = 5 ! CL

  TYPE cellule
     INTEGER                :: val
     TYPE(cellule), POINTER :: suiv
  END TYPE cellule

  TYPE rpointer
     REAL(PR), DIMENSION(:) , POINTER :: Ua
  END TYPE rpointer

  TYPE Donnees
  ! lecture dans NumMeth.data, dimensions et donnees du probleme 
  ! **************************************************************
     LOGICAL             :: Reprise
     INTEGER             :: Impre
     INTEGER             :: Icas, Itopo
     INTEGER             :: Iflux
     INTEGER             :: Isource_topo, Isource_friction, Isource_pluie
     INTEGER             :: Ifre,  Ifres
     INTEGER             :: Ndim, Nvar, KT, Ktmax 
     CHARACTER(LEN=70)   :: PbName, MeshFile
     REAL(PR)            :: DT, T, Tmax, Cfl
   
     ! lecture dans  Climites.data, CL
     ! *********************************
     INTEGER                               :: Nbord
     INTEGER       , DIMENSION(:), POINTER :: lgn, lgf1, lgf2
     INTEGER       , DIMENSION(:), POINTER :: lgc, lgd
     TYPE(rpointer), DIMENSION(:), POINTER :: Fr

     ! Topographie
     ! *********************************
     REAL(PR)        , DIMENSION(:), POINTER :: Z

  END TYPE Donnees

  TYPE MeshDef
     ! contient les caractérisiques du maillage primal
     !**************************************************
     INTEGER                              :: Npoint, Nelemt, NsegmtFr
     INTEGER, DIMENSION(:),    POINTER    :: Ndegre, Logfac, Typoint
     INTEGER, DIMENSION(:,:),  POINTER    :: Nu 
     REAL(PR), DIMENSION(:,:), POINTER    :: Coor
     REAL(PR), DIMENSION(:),   POINTER    :: AirT  , AirS
  
  END TYPE MeshDef

  TYPE MeshSeg
     ! contient les caractéristiques des segments du maillage vertex-centered
     !************************************************************************
     INTEGER                              :: Nsegmt, Nsegmtot
     INTEGER, DIMENSION(:),       POINTER :: NbVois , NbVoit, Numsegfr, LogFac
     INTEGER, DIMENSION(:,:),     POINTER :: Nubo, Nusv, Nutv
     INTEGER, DIMENSION(:,:,:),   POINTER :: TabSeg
     REAL(PR)   , DIMENSION(:),   POINTER :: LenSeg
     REAL(PR)   , DIMENSION(:,:), POINTER :: Vno, Normext
  END TYPE MeshSeg

  TYPE variables
     ! contient les inconnues du système 
     !***********************************
     INTEGER                                :: NCells
     REAL(PR),    DIMENSION(:)    , POINTER :: CellVol
     REAL(PR),    DIMENSION(:,:)  , POINTER :: Ua, Un
    
     ! pour la gestion des voisins
     INTEGER, DIMENSION(:)    , POINTER :: NbVois  ! nombre de sommets voisins 
                                                   ! de chaque sommet !
     INTEGER, DIMENSION(:,:)  , POINTER :: NuVois  ! numero des sommets voisins 
                                                   ! de chaque sommet !
  END TYPE variables


END MODULE Var_Types
