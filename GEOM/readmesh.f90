MODULE ReadMesh

  !!$ Module qui gere la lecture du maillage primal !
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  !!$ listes des subroutines : 
  !!$ - Rename
  !!$ - ReadMESH2D


CONTAINS

!------------------------------------------------------------------------------!
  SUBROUTINE Rename_Box4SW(ipd , nthd , namesd, patname)
    !     -----------------------------------------------
    !     role: concatener patname avec ipd dans namesd !
    !     namesd = 'patname''ipd'
  
    IMPLICIT NONE
 
    INTEGER,INTENT(in)              :: ipd  , nthd 
    CHARACTER*1                     :: str
    CHARACTER(LEN=70),INTENT(inout) :: namesd
    CHARACTER(LEN=70),INTENT(in)    :: patname
    !!
    !!  Local variables definition
    !!
    INTEGER   ::   i, ii , ii2, ii3, ip2, inum
    INTEGER   ::   lensd, lensd2, lensd1
 

    inum              = nthd
    lensd2            = INDEX(namesd,' ') - 1
    lensd1            = INDEX(patname,' ') - 1
  
    IF( lensd1 .NE. 0 ) THEN 
       namesd(1:lensd1) = patname(1:lensd1)
    END IF
 
    lensd             = lensd1 
    !!
    !!  Forming the name of the file associated with
    !!  the submesh definition of the current pthread
    !! 
    namesd(lensd+1:lensd+1) = '-'

    ii                = 10
    ip2               = ipd
    DO  i=1,inum
       ii2            = (ipd/ii)
       ii3            = ip2 - ii2*10
       ii3            = MAX(ii3,0)
       WRITE (str,77) ii3
       namesd(lensd+inum+2-i:lensd+inum+2-i) = str
       ii             = ii*10
       ip2            = ii2

    END DO
    !!
    lensd             = lensd + inum + 1
    namesd(lensd+1:lensd+2) = '  '

  77  FORMAT(i1.1)
    !!
  END SUBROUTINE Rename_Box4SW
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
  SUBROUTINE ReadMESH2D( Mesh, FileName)

    !     -----------------------------------------------
    !     role: lire le maillage dans le fichier Filename 
    !
    !==============================================================
    ! Coor(i,is) ::  Ieme coord du sommet is du maillage primal
    ! LogFr(is) :: contrainte (CL) du sommet is 
    ! Nu(k,jt) ::  numero global du sommet k dans le triangle jt
    ! Ndegre(jt) :: nombre de sommets de l'element jt (3, 4, ...)
    !==============================================================

    USE var_types
    IMPLICIT NONE

    ! ****************
    !   declarations
    ! ****************
    ! variables d'appel
    TYPE(MeshDef), INTENT(out)       :: Mesh
    CHARACTER(LEN=*), INTENT(in)     :: FileName

    !variables locales
    INTEGER                          :: MyUnit, j, ljt
    INTEGER                          :: NsegmtFr
    CHARACTER(LEN=30)                :: str

    !variables locales
    INTEGER, DIMENSION(:),   POINTER :: LogFr

    !***************************
    !  lecture de la geometrie
    !***************************

    MyUnit = 11

     !ouverture du fichier
    OPEN(unit=MyUnit, file=FileName)
    READ(MyUnit,*) str,j
    READ(MyUnit,'(A)') str
    READ(MyUnit,*) j
    READ(MyUnit,'(A)') str
    READ(MyUnit,*) Mesh%Npoint
    
    !allocation memoire
    ALLOCATE( Mesh%Coor(2,Mesh%Npoint) )
   
    !allocation memoire
    ALLOCATE(LogFr(Mesh%Npoint))

    NsegmtFr = 0

    !lecture des coord des sommets et du type de point (interieur ou frontiere)
    DO j = 1,Mesh%Npoint
       READ(MyUnit,*) Mesh%Coor(1,j), Mesh%Coor(2,j), LogFr(j)
       IF(LogFr(j) .NE. 0) THEN  ! si le point est la frontiere !
          NsegmtFr = NsegmtFr + 2
       ENDIF
    END DO

    Mesh%NsegmtFr = NsegmtFr
    Mesh%Typoint => LogFr

    READ(MyUnit,'(A)') str
    READ(MyUnit,*) Mesh%Nelemt
  
    ALLOCATE(  Mesh%Nu(4,Mesh%Nelemt)   )
    ALLOCATE(  Mesh%Ndegre(Mesh%Nelemt) )

    !lecture des numeros des sommets du maillage primal 	
    DO j = 1,Mesh%Nelemt
       READ(MyUnit,*) Mesh%Nu(1,j), Mesh%Nu(2,j), Mesh%Nu(3,j) , ljt
    END DO
 
    CLOSE(MyUnit)

    ! triangles !
    Mesh%Ndegre(:) = 3
    
   
  END SUBROUTINE ReadMESH2D
!------------------------------------------------------------------------------!


END MODULE ReadMesh
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
