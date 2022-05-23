MODULE IO

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


  !**************************************************************************!
  SUBROUTINE DonNumMeth(DATA)
   
    TYPE(Donnees),     INTENT(inout)     :: DATA

    DATA%Itopo = 0
    
    OPEN(66,file='NumMeth.data',FORM="formatted", STATUS="old" )

    READ(66,*) DATA%PbName
    READ(66,*) Data%Reprise
    READ(66,*) Data%Impre  , Data%Icas, Data%Itopo
    READ(66,*) Data%Isource_topo, Data%Isource_friction, Data%Isource_pluie
    READ(66,*) Data%Ndim   
    READ(66,*) Data%Iflux  , Data%Ifre  ,  Data%Ifres
    READ(66,*) Data%CFl 
    READ(66,*) Data%Tmax   , Data%Ktmax  

    CLOSE(66)

    Data%Nvar   = Data%Ndim+1
    DATA%kt     = 0
    DATA%T      = 0.0_PR

    ! verification de certaines compatibilite
    IF( (Data%Itopo.NE.0) .AND. (Data%Isource_topo .EQ. 0) ) THEN
       WRITE(6,*) 'ATTENTION probleme de compatibilite dans le cas test'
       WRITE(6,*) 'fond non plat : Data%Itopo = ',Data%Itopo
       WRITE(6,*) 'terme source de topo non active : Data%Isource_topo =',Data%Isource_topo
       STOP
    END IF
    
  END SUBROUTINE DonNumMeth
!***************************************************************************!

!***************************************************************************!
  SUBROUTINE DonClimites(DATA, Mesh, Seg)
    
    TYPE(Donnees),     INTENT(inout)     :: DATA
    TYPE(MeshDef),     INTENT(inout)     :: Mesh
    TYPE(MeshSeg),     INTENT(inout)     :: Seg

    LOGICAL   :: Found
    INTEGER   :: ibord, MaxLogfac, ilog, ifac
    REAL(PR)  :: h, u_x, u_y
  
    REAL(PR), DIMENSION(:,:), POINTER    :: UaFrl

    MaxLogfac = 0
    IF( Mesh%NsegmtFr <= 0 ) RETURN
    OPEN(unit=66,file='Climites.data',FORM="formatted", STATUS="old" )

    ! nombre de CL pour ce cas test !
    READ(66,*) Data%Nbord
    !
    IF( Data%Nbord <= 0 ) RETURN

    ALLOCATE( Data%Lgn(Data%Nbord)  )
    ALLOCATE( Data%Lgf1(Data%Nbord) , Data%Lgf2(Data%Nbord) )
    ALLOCATE( Data%Lgc(Data%Nbord)  , Data%Lgd(Data%Nbord) )
    ALLOCATE( Seg%LogFac(1:Mesh%NsegmtFr))
    ALLOCATE( UaFrl(Data%Nvar, Data%Nbord) )
    Data%Lgc = 0 ; Data%Lgd=0
    UaFrl = 0.0_PR   

    ! lgn : compteur pour nbord !
    ! Lgf1 : indice de LogFac de debut !
    ! Lgf2 : indice de LogFac de fin   ! 
    ! lgc : type de CL affecte aux LogFac( Lgf1 : Lgf2) !

    DO  ibord = 1 , Data%nbord 
       READ(66,*) Data%lgn(ibord), Data%Lgf1(ibord),&
            & Data%Lgf2(ibord), Data%lgc(ibord), Data%lgd(ibord)

       MaxLogfac = MAX( MaxLogfac, Data%lgc(ibord) )
       SELECT CASE(Data%lgc(ibord))
       CASE(4)
          READ(66,*) UaFrl(1:Data%Nvar, ibord)
       END SELECT
       
    END DO

    CLOSE(66)

    ! verification des Logiques de face (CL) !
    DO ifac = 1, Mesh%Nsegmtfr
       ilog = Seg%TabSeg(2,Seg%Numsegfr(ifac),2)
       Found = .FALSE.
       DO  ibord = 1, Data%nbord 
          IF(ilog>= Data%Lgf1(ibord) .AND. ilog <= Data%Lgf2(ibord)  ) THEN
             Seg%TabSeg(2,Seg%Numsegfr(ifac),2) = Data%lgc(ibord)
             Seg%LogFac(ifac) = ibord
             Found             = .TRUE.
             EXIT
          END IF
       END DO
       IF( .NOT. Found ) THEN
          WRITE(6,*) " Logique de face non defini " , ilog,  &
               &    " Data%nbord = ",Data%nbord 
          WRITE(6,*) " Data%lgn(ibord),         Data%Lgf1(ibord),       &
               &    Data%lgf2(ibord),         Data%lgc(ibord),       &
               &    Data%lgd(ibord) "
          DO  ibord = 1 , Data%nbord , 1
             WRITE(6,*) Data%lgn(ibord), Data%Lgf1(ibord),&
                  & Data%lgf2(ibord), Data%lgc(ibord), Data%lgd(ibord)
          END DO
          STOP
       END IF
    END DO
 
    MaxLogfac =  Data%nbord   

    ALLOCATE( Data%Fr(MaxLogfac) )
    
    ! boucle sur le nombre de CL dans le cas test !
    DO  ibord = 1, Data%nbord , 1
       SELECT CASE(Data%lgc(ibord))
          CASE(4)  !Dirichlet 
             h    = UaFrl(1,ibord)
             u_x  = UaFrl(2,ibord)
             u_y  = UaFrl(3,ibord)
             
             UaFrl(1,ibord)  = h
             UaFrl(2,ibord)  = h*u_x 
             UaFrl(3,ibord)  = h*u_y
                          
             Data%Fr(ibord)%UA => UaFrl(:,ibord) 
             
             WRITE(6,*) " ibord = ", ibord, " Ua = ", UaFrl(:,ibord)
          END SELECT
       END DO


  END SUBROUTINE DonClimites
!***************************************************************************!

!***************************************************************************!
  SUBROUTINE Reprise_in(DATA, Mesh, Var, PbName )
    
    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var
    CHARACTER(LEN=*), INTENT(inout)     :: PbName
    
    INTEGER             :: is, k, ies, Ndim
    INTEGER             :: lPbName, lensd2
    CHARACTER(LEN=70)   :: filename

    WRITE(6,*)
    WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    WRITE(6,*) '::     debut de la preparation a la reprise         ::'
    WRITE(6,*) '::                                                  ::'
    WRITE(6,*) ':: ne pas oublier de copier .save dans .restart     ::'
    WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    WRITE(6,*)
    
    lPbName                          = INDEX(PbName,' ') - 1
    filename(1:lPbName)              = PbName(1:lPbName)
    filename(lPbName+1:lPbName+9)    = '.restart '
    lensd2                           = lPbName+9

    OPEN( UNIT=66, FILE=filename(1:lensd2), FORM='unformatted')
    REWIND(66)
  
    lensd2            = INDEX(filename,' ') - 1
    IF (Data%impre > 5 ) WRITE (6,*) "    ----------------------------------"
    IF (Data%impre > 5 ) WRITE(6, *) "   file name:  ", filename(1:lensd2)
    !     
    !     Lecture au format restart
    !     ------------------------    
    !
    READ(66) Data%T     , Data%KT
    !
    IF (Data%impre > 5 ) WRITE(6, *) " Temps =",Data%T , " Iteration = ", Data%KT
    IF (Data%impre > 5 ) WRITE(6, *) "    ----------------------------------"
    Ndim = Data%Ndim

    READ(66) ((Mesh%Coor(k,is), k = 1,Ndim), is = 1,Mesh%Npoint)
    READ(66) ((Var%Ua(ies,is), is = 1,Var%NCells), ies = 1,Data%Nvar)
    CLOSE(66)

    IF (Data%impre > 5 ) WRITE(6,*) ' Ua  min =  ' , MINVAL( Var%UA, DIM=2 )
    IF (Data%impre > 5 ) WRITE(6,*) ' Ua  max =  ' , MAXVAL( Var%UA, DIM=2 )
    IF (Data%impre > 5 ) WRITE(6,*) ' ..... Fin de la reprise .... '

    WRITE (6,*)
    WRITE (6,*) "----------------------------------"
    WRITE (6,*) "CONFIRMATION REPRISE"
    WRITE (6,*) "date de reprise =",Data%T
    WRITE (6,*) "date de fin     =",Data%Tmax
    WRITE (6,*)
    WRITE (6,*) "----------------------------------"
    WRITE (6,*) "CONFIRMATION REPRISE"
    WRITE (6,*) "iteration de reprise =",Data%KT
    WRITE (6,*) "nb d''iteration max  =",Data%KTmax
    WRITE (6,*) "----------------------------------"
    WRITE (6,*)

    ! verification
    IF(Data%T >= Data%Tmax) then
       WRITE(6,*) "Probleme entre la date de reprise et la date de fin"
       WRITE (6,*)
       STOP
    END IF

    IF(Data%KT >= Data%KTmax) then
       WRITE(6,*) "Probleme entre l''iteration de reprise et le nb d''iteration max"
       WRITE (6,*)
       STOP
    END IF

  END SUBROUTINE Reprise_in
!***************************************************************************!

!***************************************************************************!
  SUBROUTINE Reprise_out(DATA, Mesh, Var, PbName )
    
    TYPE(Donnees),     INTENT(inout)    :: DATA
    TYPE(MeshDef) ,    INTENT(inout)    :: Mesh
    TYPE(Variables),   INTENT(inout)    :: Var
    CHARACTER(LEN=*),  INTENT(inout)    :: PbName
    
    INTEGER             :: is, k, ies, Ndim
    INTEGER             :: lPbName, lensd2
    CHARACTER(LEN=70)   :: filename

    IF (Data%impre > 5 )  WRITE(6,*) ' ..... preparation de la reprise .... '

    lPbName                          = INDEX(PbName,' ') - 1
    filename(1:lPbName)              = PbName(1:lPbName)
    filename(lPbName+1:lPbName+6)    = '.save '
    lensd2                           = lPbName+6

    OPEN( UNIT=66, FILE=filename(1:lensd2), FORM='unformatted')
    REWIND(66)
    !!
    lensd2            = INDEX(filename,' ') - 1
    IF (Data%impre > 5 ) WRITE (6,*) "    ----------------------------------"
    IF (Data%impre > 5 ) WRITE(6, *) "   file name:  ", filename(1:lensd2)
    IF (Data%impre > 5 ) WRITE(6, *) "    ----------------------------------"
    !       
    !     Ecriture au format restart
    !     ------------------------    
    !
    WRITE(66) Data%T     , Data%KT
    !
    Ndim = Data%Ndim
    !
    WRITE(66) ((Mesh%Coor(k,is), k = 1,Ndim), is = 1,Mesh%Npoint)
    WRITE(66) ((Var%Ua(ies,is), is = 1,Var%NCells), ies = 1,Data%Nvar)
    CLOSE(66)
    IF (Data%impre > 5 ) WRITE(6,*) ' ..... Fin de la preparation a la reprise .... '


  END SUBROUTINE Reprise_out
!***************************************************************************!

!***************************************************************************!
SUBROUTINE CellVertexVtk(DATA, Mesh, Var, PbName, flag_topo)
    
    TYPE(Donnees)   , INTENT(in)     :: DATA
    TYPE(MeshDef)   , INTENT(in)     :: Mesh
    TYPE(Variables) , INTENT(in)     :: Var
    CHARACTER(LEN=70),INTENT(in)     :: PbName
    integer          ,intent(in)     :: flag_topo
   
    REAL(PR)          :: u_x, u_y
    CHARACTER(LEN=75) :: vtk="results/ "
    INTEGER           :: is, jt
    INTEGER           :: lensd3, lPbName


    IF (Data%Impre > 5 ) WRITE (6,FMT = *) " --------------------------------------------"
    IF (Data%Impre > 5 ) WRITE (6,FMT = *) " ---------      passage dans vtk   ----------"
    IF (Data%Impre > 5 ) WRITE (6,FMT = *) " --------------------------------------------"

    ! ECRITURE sous FICHIER vtk !

    lPbName                    = INDEX(PbName,' ') - 1
    vtk(9:lPbName+8)             = PbName(1:lPbName)
    vtk(lPbName+9:lPbName+13)   = '.vtk '
    lensd3                     = lPbName+13
    OPEN(UNIT=61,FILE=vtk(1:lensd3) )

    WRITE(61,'(A)')'# vtk DataFile Version 3.0'
    WRITE(61,'(A)')'# Solution du M1 TR'
    WRITE(61,'(A)')'ASCII'
    WRITE(61,'(A)')'DATASET UNSTRUCTURED_GRID'
    WRITE(61,'(A,I7,A)')'POINTS', Mesh%Npoint,'  float'

    DO is = 1,Mesh%Npoint
       if (flag_topo == 1) then
          WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') Mesh%coor(1,is), Mesh%coor(2,is),  DATA%Z(is)
       else
          WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') Mesh%coor(1,is), Mesh%coor(2,is),  max(1.E-08,Var%Ua(1,is)) + DATA%Z(is)
          end if
    END DO
    WRITE(61,'(A,1x,I7,1x,I6)')'CELLS',Mesh%Nelemt, 4*Mesh%Nelemt

    DO jt = 1,Mesh%Nelemt
       WRITE(61,'(I1,1x,I7,1x,I7,1x,I7)') 3, Mesh%Nu(1,jt)-1, Mesh%Nu(2,jt)-1, Mesh%Nu(3,jt)-1
    END DO

    WRITE(61,'(A,1x,I7)')'CELL_TYPES', Mesh%Nelemt
    DO jt = 1,Mesh%Nelemt
       WRITE(61,'(I1)') 5
    ENDDO

    WRITE(61,'(A,1x,I7)')'POINT_DATA',Mesh%Npoint
    WRITE(61,'(A)')'SCALARS hauteur float'
    WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

    DO is = 1,Mesh%Npoint
       WRITE(61,'(ES20.7)') max(1.E-08,Var%Ua(1,is))
    END DO

    WRITE(61,'(A)')'SCALARS topographie float'
    WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

    DO is = 1,Mesh%Npoint
       WRITE(61,'(ES20.7)') max(1.E-08,DATA%Z(is))
    END DO

    WRITE(61,'(A)')'SCALARS surface float'
    WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

    DO is = 1,Mesh%Npoint
       WRITE(61,'(ES20.7)') max(1.E-08,Var%Ua(1,is)) + DATA%Z(is)
    END DO

    WRITE(61,'(A)')'VECTORS vitesse float'

    DO is = 1,Mesh%Npoint
       u_x = vitesse(Var%Ua(1,is) , Var%Ua(2,is))
       u_y = vitesse(Var%Ua(1,is) , Var%Ua(3,is))
       WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') max(1.E-08,u_x), max(1.E-08,u_y), 0.
    END DO
   
    CLOSE(61)

  END SUBROUTINE CellVertexVtk
!***************************************************************!

END MODULE IO
