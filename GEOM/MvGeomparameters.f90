MODULE ParmGeom


  !!$ Module qui gere les caractéristiques des segments  !
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  !!$ listes des subroutines : 
  !!$ - GeomAires
  !!$ - CellVertexMvvno

CONTAINS
  
  !---------------------------------------------------------------------!
  SUBROUTINE GeomAires(Mesh)
    !---------------------------------------------
    !     role: calcule l'aire des cellules et des triangles 
    !
    !==============================================================
    ! Airs(i) ::  Aire de la cellule I 
    ! Airt(jt) :: Aire du triangle jt 
    !==============================================================
    !
    USE var_types
    ! 
    !.. Implicit Declarations .. 
    IMPLICIT NONE
 
    !.. Formal Arguments .. 
    TYPE(MeshDef) , INTENT(inout)      :: Mesh

    !.. Local variables .. 
    INTEGER                            :: Ns, Nt
    INTEGER, DIMENSION(:)  , POINTER   :: Ndegre
    INTEGER, DIMENSION(:,:), POINTER   :: Nu
    REAL(PR),DIMENSION(:,:), POINTER   :: coor
    REAL(PR),DIMENSION(:)  , POINTER   :: airt, airs
    
    !.. Local Scalars .. 
    INTEGER       :: is1,is2,is3,jt
    INTEGER       :: is4
    REAL(PR)      :: aire,x12,x13,x14,y12,y13,y14
     
    !.. Intrinsic Functions .. 
    INTRINSIC REAL
    ! 
    ! ... Executable Statements ...
   
    NULLIFY( Ndegre, Nu, Coor, Airt, Airs)

    Ns       = Mesh%Npoint
    Nt       = Mesh%Nelemt
    Ndegre   =>Mesh%Ndegre
    Nu       =>Mesh%Nu
    Coor     =>Mesh%Coor

    ALLOCATE(  Airt(1:Nt) )       ! aire des triangles !
    ALLOCATE(  Airs(1:Ns) )       ! aire des cellules !

    Airs = 0.0_PR

    ! boucle sur les elements !
    DO jt = 1,Nt
       is1 = Nu(1,jt)
       is2 = Nu(2,jt)
       is3 = Nu(3,jt)
       !
       x12 = Coor(1,is2) - Coor(1,is1)
       y12 = Coor(2,is2) - Coor(2,is1)
       !
       x13 = Coor(1,is3) - Coor(1,is1)
       y13 = Coor(2,is3) - Coor(2,is1)
       aire = ABS(0.5_PR * (x12*y13 - y12*x13))
       !
       is4 = -1
       IF (Ndegre(jt) == 4) THEN
          is4 = Nu(4,jt)
          !
          x14 = Coor(1,is4) - Coor(1,is1)
          y14 = Coor(2,is4) - Coor(2,is1)
          !
          aire = aire + ABS(0.5*(x14*y13-y14*x13))
       ENDIF
       !
       airt(jt) = aire
       !
       airs(is1) = airs(is1) + aire/REAL(Ndegre(jt),PR)
       airs(is2) = airs(is2) + aire/REAL(Ndegre(jt),PR)
       airs(is3) = airs(is3) + aire/REAL(Ndegre(jt),PR)
       !
       IF (Ndegre(jt) == 4) THEN
          IF (is4 == -1) THEN
             print*,'PROBLEME Dans MvGeomparameters/GeomAires/is4'
             stop
          END IF
          airs(is4) = airs(is4) + aire/REAL(Ndegre(jt),PR)
       ENDIF

    END DO
    
    NULLIFY( Mesh%Airt, Mesh%Airs)
    Mesh%Airt  => Airt
    Mesh%Airs  => Airs


  END SUBROUTINE GeomAires
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
  SUBROUTINE CellVertexMvvno( Mesh, Seg, Dt )
    !     -----------------------------------------------
    !     role: calcule les coordonnees de la normale sortante 
    !     et la longueur de chaque segment 
    !
    !==============================================================
    ! Vno(i,jseg) ::  Ieme coord de la normale sortante au segment jseg 
    ! Normext(i,jseg) ::  Ieme coord de la normale sortante unitaire au segment jseg 
    ! LenSeg(jseg) ::  longueur du segment jseg 
    !==============================================================
    !
    ! Modules utilises
    ! ----------------
    USE var_types 
  
    IMPLICIT NONE
    ! 
    ! Type et vocations des Variables d'appel
    !----------------------------------------
    TYPE(MeshDef) , INTENT(inout)            :: Mesh
    TYPE(MeshSeg) , INTENT(inout)            :: Seg
    REAL(PR)      , INTENT(in) , OPTIONAL    :: Dt 

    ! Variables locales scalaires
    !----------------------------
    INTEGER :: j, is1, is2, jdeg, compt_segtot
    INTEGER :: js1, js2, js3, js4, iseg, jt, k, k1
    INTEGER :: Nsegtot
    REAL(PR):: xija, yija, xGta, yGta, Dtl
    ! 
    ! Variables locales tableaux
    !----------------------------
    REAL(PR)   , DIMENSION(2)   :: Nij  
  
    INTEGER, DIMENSION(:)    , POINTER :: ndegre
    INTEGER, DIMENSION(:,:)  , POINTER :: Nu, Nutv
    INTEGER, DIMENSION(:,:,:), POINTER :: TabSeg
    REAL(PR),DIMENSION(:)    , POINTER :: LenSeg
    REAL(PR),DIMENSION(:,:)  , POINTER :: Coor, Vno, Normext
    ! 
    ! 
    ! Fonctions Internes
    INTRINSIC MOD
    ! 
    !--------------------------------------------------------------------
    !                Partie executable de la procedure
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !
    IF( .NOT. PRESENT(Dt) ) THEN
       WRITE(6,*) " passage dans MvCellVertexMvvno sans Dt "
    END IF

    NULLIFY( Ndegre, Nu, Coor, TabSeg, Nutv, Vno, Normext )

    Ndegre    =>Mesh%Ndegre
    Nu        =>Mesh%Nu
    Coor      =>Mesh%Coor
    Nsegtot   = Seg%Nsegmtot
    Nutv      =>Seg%Nutv
    TabSeg    =>Seg%TabSeg
   
    ALLOCATE(Normext(1:2,1:Nsegtot),Vno(1:2,1:Nsegtot),LenSeg(1:Nsegtot))

    IF( .NOT. PRESENT(Dt) ) THEN
       Dtl  = 1.0_PR
    ELSE
       Dtl  = Dt
    END IF

    Vno          = 0._PR
    LenSeg       = 0._PR
    Nij          = 0._PR
    compt_segtot = 1

    ! boucle sur les segments !
    DO iseg = 1,Nsegtot

          ! on recupere les cellules voisines du segment compt_segtot !
          is1 = Seg%TabSeg(1,compt_segtot,1)
          is2 = Seg%TabSeg(2,compt_segtot,1)

          ! calcul du milieu de l'arete is1-->is2 !
          xija = (coor(1,is1) + coor(1,is2) )*0.5_PR
          yija = (coor(2,is1) + coor(2,is2) )*0.5_PR
          
          DO j = 1, 2

             jt  = Nutv(j,iseg)
             IF ( jt == 0 ) CYCLE
             js1 = Nu(1,jt)   ; js2 = Nu(2,jt) ; js3 = Nu(3,jt)     
                       
             jdeg = Ndegre(jt)
             SELECT CASE (jdeg)
             CASE(4)
                js4 = Nu(4,jt) 
                xGta  = (   coor(1,js1) + coor(1,js2)&
                     &    + coor(1,js3) + coor(1,js4)  )*0.25_PR
                yGta  = (   coor(2,js1) + coor(2,js2)&
                     &    + coor(2,js3) + coor(2,js4) )*0.25_PR
                           
             CASE(3)  
                ! calcul des centres de gravite du triangle jt !
                IF(Seg%TabSeg(1,compt_segtot,2) == 0) THEN
                  
                   xGta  = (   coor(1,js1) + coor(1,js2)&
                        &    + coor(1,js3) )/3.0_PR
                   yGta  = (   coor(2,js1) + coor(2,js2)&
                        &    + coor(2,js3) )/3.0_PR
                   Nij  =    (/ (yGta-yija) ,  - (xGta-xija) /) 

                ELSE 
                   ! si c'est un segment frontiere, on calcule sa normale avec 
                   ! le premier de ses sommets !

                   xGta = coor(1,is1)
                   yGta = coor(2,is1)
                   Nij  =    (/ - (yGta-yija) ,  (xGta-xija) /)

                ENDIF

             CASE DEFAULT
                print*,'PROBLEME Dans MvGeomparameters/CellVertexMvvno'
                stop
                            
             END SELECT
             
             DO k = 1,jdeg
                js1   = NU(k,jt)
                k1    = MOD(k,jdeg) + 1
                js2   = NU(k1,jt)
                IF (js1==is2 .AND. js2==is1) THEN
                   Nij(1:2)   = - Nij(1:2)
                   EXIT
                END IF
             END DO

             vno(1:2,compt_segtot) = vno(1:2,compt_segtot) + Nij(1:2)
                              
          END DO
          !     Normale sortante de 1 -- vers ---> 2 (normalisee)
          !
          LenSeg(compt_segtot) = SQRT( Vno(1,compt_segtot)**2 + Vno(2,compt_segtot)**2 )
          Normext(1,compt_segtot) = Vno(1,compt_segtot)/LenSeg(compt_segtot)
          Normext(2,compt_segtot) = Vno(2,compt_segtot)/LenSeg(compt_segtot)
          
          compt_segtot = compt_segtot + 1
       
    END DO
    
    NULLIFY( Seg%Normext, Seg%Vno, Seg%LenSeg) 
    Seg%Normext   => Normext
    Seg%Vno       => Vno
    Seg%LenSeg    => LenSeg

    IF( .NOT. PRESENT(Dt) ) THEN
       WRITE(6,*) ' Fin de  CellVertexMvvno '
       CALL verif()
    END IF
    
    NULLIFY( Ndegre, Nu, Coor, TabSeg, Nutv )


  CONTAINS

    SUBROUTINE verif()
      REAL(PR), DIMENSION( 2, SIZE(Mesh%COOR, DIM=2) ) :: SNvol

      SNvol = 0 

      DO iseg = 1,Nsegtot
         is1 = TabSeg(1,iseg,1)
         is2 = TabSeg(2,iseg,1)
         SNvol(1:2, is1) = SNvol(1:2, is1) + vno(1:2,iseg)   
         SNvol(1:2, is2) = SNvol(1:2, is2) - vno(1:2,iseg)

      END DO

      WRITE(6,* ) " Min val sur les cellules ", MINVAL(SNvol) 
      WRITE(6,* ) " Max val sur les cellules ", MAXVAL(SNvol) 


    END SUBROUTINE verif

  END SUBROUTINE CellVertexMvvno
  !---------------------------------------------------------------------------!


END MODULE ParmGeom


