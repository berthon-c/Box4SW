MODULE SegmtsGeom

  !!$ Module qui gere la création du maillage vertex centered !
  !---------------------------------------------------------------!
  !!$ Auteurs : B. Nkonga , C. Berthon, R. Turpault, C. Sarazin
  !---------------------------------------------------------------!

  !!$ listes des subroutines : 
  !!$ - Segments2DP
  !!$ - WhichSeg
  !!$ - compute_voisin


CONTAINS
!------------------------------------------------------------------------------!
  SUBROUTINE Segments2DP(Mesh, Seg, DATA, Var)
    !---------------------------------------------
    !     role: creer tous les tableaux relatifs aux segments  
    !
    !==============================================================
    ! Nubo(k,jseg) ::  sommet numero K du segment jseg
    ! TabSeg(k,jseg,tag) :: numero du keme sommet du seg jseg et tag pour indiquer
    !                     si le point est a l'interieur ou sur la frontiere 
    !                     et sa contrainte (CL) s'il est sur la frontiere
    ! NumSegFr(jsegfr) :: numero global (dans TabSeg) du segment frontiere jsegfr
    ! Nusv(k,js) ::  numero du Keme sommet tel qu'avec le
    !                segment js, ils forment un element
    ! Nutv(jt,jseg) ::  numero du jteme (1er ou 2eme) element contenant le segment jseg
    ! NbVois(k)  ::  nombre de sommets voisins du sommet k 
    ! NbVois(k,i)::  numero du Ieme sommet voisin du sommet k
    !==============================================================
    !
    USE var_types

    IMPLICIT NONE

    ! ****************
    !   declarations
    ! ***************
    ! variables d'appel

    TYPE(MeshDef), INTENT(inout)        :: Mesh
    TYPE(MeshSeg), INTENT(out)          :: Seg
    TYPE(Donnees), INTENT(in)           :: DATA
    TYPE(Variables), INTENT(inout)      :: Var

    ! variables locales
    INTEGER                          :: Npoint, Nelemt
    INTEGER                          :: Nsegmt, NsegmtFr, Nsegmtot, Info
    INTEGER                          :: i, is, js, iv, jt, k, is1, is2
    INTEGER                          :: compt_seg, compt_segtot, compt_segfr
    INTEGER                          :: ismin, ismax, jseg
    INTEGER                          :: next, ntyp, noth

    ! variables locales
    INTEGER, DIMENSION(1:2)                     :: indj
    INTEGER, DIMENSION(:),   POINTER            :: NbSegFr, Numsegfr, NbVois, NbVoit
    INTEGER, DIMENSION(:,:), POINTER            :: Nubo, Nusv, Nutv
    INTEGER, DIMENSION(:,:,:), POINTER          :: TabSeg
    TYPE(cellule), POINTER                      :: NewCell, PtCell, PtCellPred
    TYPE(cellule), DIMENSION(:), POINTER        :: hashtable
    INTEGER, DIMENSION(SIZE(Mesh%Coor, DIM=2))  :: NSupVois 

    ! *************************** !
    !       Initialisations       !
    ! *************************** !

    Npoint   = Mesh%Npoint
    Nelemt   = Mesh%Nelemt
    Nsegmtfr = Mesh%Nsegmtfr

    NULLIFY(Nubo, NbVois, TabSeg, Nusv, Nutv)

    IF( Data%Impre > 2 )  WRITE(6,*) "Npoint, Nelemt, NsegmtFr",&
         &Mesh%Npoint, Mesh%Nelemt, Mesh%NsegmtFr

    ! *******************************************
    !   construction de hashtable, NSupVois, Nsegmt
    ! *******************************************

    ALLOCATE(hashtable(1:Npoint),STAT=Info)
    WRITE(6,*) " Apres Hashtable Info =", Info

    DO is = 1, Npoint
       hashtable(is)%val = is
       NULLIFY(hashtable(is)%suiv)
    END DO

    NSupVois(1:Npoint) = 0
    Nsegmt = 0

    ! construction de la table des segments 

    NULLIFY(PtCell, PtCellPred )

    DO jt = 1, Nelemt

       DO k = 1, Mesh%Ndegre(jt)

          ntyp = Mesh%Ndegre(jt)
          next = MOD(k, ntyp) + 1
          is1  = Mesh%Nu(   k, jt)
          is2  = Mesh%Nu(next, jt)

          ismin = MIN(is1,is2)
          ismax = MAX(is1,is2)

          ! recherche de la position d'insertion de ismax dans
          ! le tableau hashtable(ismin). on utilise une
          ! recherche sequentielle simple.
          ! initialisation des pointeurs
          PtCell     => hashtable(ismin)%suiv
          PtCellPred => hashtable(ismin)

          ! on recherche la place du nouveau voisin
          DO
             IF (.NOT. ASSOCIATED(PtCell)) EXIT
             IF (PtCell%val .GT. ismax)    EXIT
             ! on a trouve la place, on arrete

             PtCellPred => PtCell
             PtCell     => PtCell%suiv
          END DO

          ! on teste si le sommet n'est pas deja la
          IF (PtCellPred%val .LT. ismax) THEN

             ! creation de la nouvelle cellule
             ALLOCATE(NewCell)
             NewCell%val  = ismax

             ! insertion de la nouvelle cellule
             NewCell%suiv    => PtCell
             PtCellPred%suiv => NewCell

             ! on a rajoute un voisin, donc on incremente
             ! NSupVois(ismin) et Nsegmt
             NSupVois(ismin) = NSupVois(ismin) + 1
             Nsegmt = Nsegmt + 1

          END IF
       END DO

    END DO
   
    ! nombre de segment total !
    Nsegmtot = Nsegmt + NsegmtFr

    ! le tableau temporaire hashtable est constitue, 
    ! on le recopie dans le tableau Nubo

    ! **************************** !
    !   creation du tableau Nubo   !
    ! **************************** !

    IF( Data%Impre > 2 ) WRITE(6,*) " Stockage du tableau Nubo, Nsegmt = ",Nsegmt

    ! allocation memoire
    ALLOCATE(Nubo(1:2,1:Nsegmt), NbVois(1:Npoint))
    
    NbVois(1:Npoint)   = 0
    k           = 0
      
    ! boucle sur les sommets du maillage !
    DO is = 1, Npoint

       ! si le sommet is a un(des) voisin(s)
       IF (NSupVois(is) .NE. 0) THEN
          PtCell => hashtable(is)%suiv

          ! on prend tous les voisins du sommet is
          DO iv = 1, NSupVois(is)
             js            = PtCell%val
             NbVois(is)    = NbVois(is) + 1
             NbVois(js)    = NbVois(js) + 1
             k             = k +1
             Nubo(1,k) = is
             Nubo(2,k) = PtCell%val
        
             ! on pointe sur le voisin suivant
             PtCell => PtCell%suiv
          END DO

       END IF
    END DO 

    !maintenant que la numerotation dans Nubo est correcte,
    !on construit les tableaux TabSeg, Nusv, Nutv

    ! ***************************************
    !   construction du tableau TabSeg
    ! ***************************************

    IF( Data%Impre > 2 ) WRITE(6,*) " Stockage du tableau TabSeg, Nsegmtot = ",Nsegmtot

    ALLOCATE(TabSeg(1:2,1:Nsegmtot,1:2),Numsegfr(1:(Nsegmtfr)),NbSegFr(1:Npoint))
    TabSeg(1:2,1:Nsegmtot,1:2) = 0
    NbSegFr(1:Npoint)          = 0
    Numsegfr(1:(Nsegmtfr))     = 0
    
    compt_seg    = 1  ! compteur pour les segments internes !
    compt_segfr  = 1  ! compteur pour les segments frontiere !
    compt_segtot = 1  ! compteur pour tous les segments !

    ! boucle sur les cellules !
    DO is = 1 , Npoint-1
       
       ! on parcourt tous les segments qui ont le sommet is comme voisin !
       DO WHILE( is == Nubo(1,compt_seg))

          ! si le point est sur la frontiere !
          IF(Mesh%Typoint(is).NE.0) THEN 
             ! si les 2 sommets sont sur la frontiere et qu'on n'a pas depasse 4 segments fr!
             IF((Mesh%Typoint(Nubo(2,compt_seg)).GT.0).and.&
                  &(NbSegFr(Nubo(1,compt_seg)).LT.4).and.&
                  &(NbSegFr(Nubo(2,compt_seg)).LT.4)) THEN
                
                ! segment interieur !
                TabSeg(1,compt_segtot,1) = Nubo(1,compt_seg)
                TabSeg(2,compt_segtot,1) = Nubo(2,compt_seg)
                
                compt_seg    = compt_seg    + 1
                compt_segtot = compt_segtot + 1

                ! segments frontieres !
                TabSeg(1,compt_segtot,1) = Nubo(1,compt_seg-1)
                TabSeg(2,compt_segtot,1) = Nubo(2,compt_seg-1)
                TabSeg(1,compt_segtot,2) = -1
                TabSeg(2,compt_segtot,2) = MIN(Mesh%Typoint(Nubo(1,compt_seg-1)),Mesh%Typoint(Nubo(2,compt_seg-1))) 
                Numsegfr(compt_segfr)    = compt_segtot
                
                compt_segtot = compt_segtot + 1
                compt_segfr  = compt_segfr  + 1

                TabSeg(1,compt_segtot,1) = Nubo(2,compt_seg-1)
                TabSeg(2,compt_segtot,1) = Nubo(1,compt_seg-1)
                TabSeg(1,compt_segtot,2) = -1
                TabSeg(2,compt_segtot,2) = MIN(Mesh%Typoint(Nubo(1,compt_seg-1)),Mesh%Typoint(Nubo(2,compt_seg-1))) 
                Numsegfr(compt_segfr)    = compt_segtot

                compt_segtot = compt_segtot + 1
                compt_segfr  = compt_segfr  + 1

                NbSegFr(Nubo(1,compt_seg-1)) = NbSegFr(Nubo(1,compt_seg-1)) + 2
                NbSegFr(Nubo(2,compt_seg-1)) = NbSegFr(Nubo(2,compt_seg-1)) + 2
             ELSE
                ! segment interieur !
                TabSeg(1,compt_segtot,1) = Nubo(1,compt_seg)
                TabSeg(2,compt_segtot,1) = Nubo(2,compt_seg)
                
                compt_seg    = compt_seg    + 1
                compt_segtot = compt_segtot + 1
             ENDIF

          ELSE
             ! segment interieur !
             TabSeg(1,compt_segtot,1) = Nubo(1,compt_seg)
             TabSeg(2,compt_segtot,1) = Nubo(2,compt_seg)

             compt_seg    = compt_seg    + 1
             compt_segtot = compt_segtot + 1
          ENDIF

          if ( compt_seg>Nsegmt) exit ! modif a verifier

       ENDDO
          
    ENDDO

    ! modification a verifier
    NULLIFY(Seg%Nubo,Seg%Nusv,Seg%Nutv,Seg%NbVois,Seg%NbVoit)   ! modification a verifier
    Seg%Nsegmt   =  Nsegmt
    Seg%Nsegmtot =  Nsegmtot
    Seg%Nubo     => Nubo
    Seg%NbVois   => NbVois
    Seg%NbVoit   => NbVoit
    Seg%Numsegfr => Numsegfr
    Seg%TabSeg   => TabSeg
    ! modification a verifier

    CALL compute_voisin(Seg, Var)

    ! **************************************** !
    !   construction des tableaux Nutv, Nusv   !
    ! **************************************** !

     IF( Data%Impre > 2 ) WRITE(6,*) " Stockage des tableaux Nusv et Nutv  "

    ! initialisation des tableaux 
     ALLOCATE(Nusv(2,Nsegmtot), Nutv(2,Nsegmtot))
     Nusv(1:2,1:Nsegmtot) = 0
     Nutv(1:2,1:Nsegmtot) = 0

     ! boucle sur les elements !
     DO jt = 1, Nelemt

        ! boucle sur les sommets de l'element !
        DO k = 1, Mesh%Ndegre(jt)  ! 

          ntyp = Mesh%Ndegre(jt)
          next = MOD(   k,ntyp) + 1
          noth = MOD(next,ntyp) + 1
          is1 = MIN(Mesh%Nu(k,jt),Mesh%Nu(next,jt))
          is2 = MAX(Mesh%Nu(k,jt),Mesh%Nu(next,jt))

          ! si les deux sommets sont des sommets frontieres !
          IF ((Mesh%Typoint(is1).NE.0).AND.(Mesh%Typoint(is2).NE.0)) THEN 
             ! on cherche le/les segments qui coupent l'arete [is1;is2]
             indj = WhichSeg(1, Nsegmtot, TabSeg, is1, is2, -1)
             DO i = 1,2
                
                jseg= indj(i)

                IF( jseg .NE. 0 ) THEN
                   IF (Nutv(1,jseg) .EQ. 0) THEN
                      Nutv(1,jseg) = jt
                      Nusv(1,jseg) = Mesh%Nu(noth,jt)
                   ELSE
                      IF (Nutv(2,jseg) .EQ. 0) THEN
                         Nutv(2,jseg) = jt
                         Nusv(2,jseg) = Mesh%Nu(noth,jt)
                      ELSE
                         WRITE(6,*) " Probleme dans la construction des Nutv "
                         STOP
                      END IF
                   END IF
                ENDIF
             ENDDO

             ! on cherche le/les segments qui coupent l'arete [is2;is1]
             indj = WhichSeg(1, Nsegmtot, TabSeg, is2, is1, -1)
             DO i = 1,2
                jseg= indj(i)
                
                IF( jseg .NE. 0 ) THEN
                   IF (Nutv(1,jseg) .EQ. 0) THEN
                      Nutv(1,jseg) = jt
                      Nusv(1,jseg) = Mesh%Nu(noth,jt)
                   ELSE
                      IF (Nutv(2,jseg) .EQ. 0) THEN
                         Nutv(2,jseg) = jt
                         Nusv(2,jseg) = Mesh%Nu(noth,jt)
                      ELSE
                         WRITE(6,*) " Probleme dans la construction des Nutv "
                         STOP
                      END IF
                   END IF
                ENDIF
             ENDDO
             
          ! si un des sommets n'est pas sur la frontiere !
          ELSE 
             indj = WhichSeg(1, Nsegmtot, TabSeg, MIN(is1,is2),MAX(is1,is2),0)
             jseg= indj(1)
                       
             IF( jseg == 0 ) CYCLE
             IF (Nutv(1,jseg) .EQ. 0) THEN
                Nutv(1,jseg) = jt
                Nusv(1,jseg) = Mesh%Nu(noth,jt)
             ELSE
                IF (Nutv(2,jseg) .EQ. 0) THEN
                   Nutv(2,jseg) = jt
                   Nusv(2,jseg) = Mesh%Nu(noth,jt)
                ELSE
                   WRITE(6,*) " Probleme dans la construction des Nutv "
                   STOP
                END IF
             END IF
          ENDIF
             
       END DO
    END DO
 
    ! Nettoyage memoire !

    DO is = 1, Npoint
       IF (NSupVois(is).NE.0) THEN
          PtCellPred => hashtable(is)%suiv
          PtCell => PtCellPred%suiv
          DO
             DEALLOCATE(PtCellPred)
             IF (.NOT.ASSOCIATED(PtCell)) EXIT
             PtCellPred => PtCell
             PtCell => PtCell%suiv
          END DO
       END IF
    END DO

    NULLIFY(Seg%Nubo,Seg%Nusv,Seg%Nutv,Seg%NbVois,Seg%NbVoit)
    Seg%Nsegmt   =  Nsegmt
    Seg%Nsegmtot =  Nsegmtot
    Seg%Nubo     => Nubo
    Seg%Nusv     => Nusv
    Seg%Nutv     => Nutv
    Seg%NbVois   => NbVois
    Seg%NbVoit   => NbVoit
    Seg%Numsegfr => Numsegfr
    Seg%TabSeg   => TabSeg

    IF( Data%Impre > 2 )WRITE(6,*) " Sortie de la procedure Segments2DP"

  CONTAINS
   
    !-----------------------------------------------------------!
    ! fonction qui recherche le numero du segment qui coupe 
    ! l'arete [is1;is2] du maillage primal !
    FUNCTION WhichSeg(idebut, ifin, TabSeg, is1, is2, fr)

      ! variables d'appel
      INTEGER, INTENT(in)                :: is1, is2, idebut, ifin, fr
      INTEGER, DIMENSION(:,:,:), POINTER :: TabSeg
       
      ! variable de sortie 
      INTEGER, DIMENSION(1:2)            :: WhichSeg

      ! variables locales
      INTEGER                            :: TheSeg, debut, fin
      INTEGER                            :: incr, ismin, ismax

      ! initialisation des variables
      WhichSeg = 0
      ismin = is1
      ismax = is2
      debut = idebut
      fin   = ifin
      jseg  = 0
      incr  = 0

      IF( ifin == 0 ) RETURN
      
      ! recherche du premier point is1 dans TabSeg
      DO jseg = debut,fin

         ! fin de la recherche, on sort
         IF ((TabSeg(1,jseg,1).EQ.ismin)) then ! .OR.(jseg.EQ.jseg2)) THEN 
            TheSeg = jseg 
            EXIT
         END IF

      END DO

      ! si ismin est dans TabSeg, on cherche ismax
      IF (TabSeg(1,jseg,1).EQ.ismin) THEN
         
         ! ruse pour forcer jseg a diminuer ou augmenter
         IF (TabSeg(2,jseg,1) .GT. ismax) THEN
            incr = - 1
         ELSE
            incr = 1
         END IF
         
         ! recherche de ismax
         DO
            IF ((TabSeg(1,jseg,1).EQ.ismin).AND.(TabSeg(2,jseg,1).EQ.ismax)) THEN
               WhichSeg(1) = jseg
               EXIT
            END IF

            jseg = jseg + incr

            IF ( jseg .LT. idebut .OR. jseg .GT. Ifin ) THEN
               WhichSeg(1) =  0
               EXIT
            END IF
         END DO
      ELSE
         WhichSeg(1) = 0
      END IF
              
      ! si les deux points de l'arete [is1,is2] sont sur la frontiere,
      ! on recherche un eventuel segment frontiere sur cette arete !
      IF (fr == -1) THEN
         jseg = jseg + incr
         
         ! recherche de ismax dans la direction de incr (+- 1)
         DO
            if (jseg>Ifin) exit  ! modif a verifier
            
            IF ((TabSeg(1,jseg,1).EQ.ismin).AND.(TabSeg(2,jseg,1) .EQ. ismax)) THEN
               WhichSeg(2) = jseg
               EXIT
            END IF
            
            jseg = jseg + incr
            IF ( jseg .LT. idebut .OR. jseg .GT. Ifin ) THEN
               WhichSeg(2) =  0
                EXIT
            END IF
         END DO
         
         ! si on n'a pas trouve de 2eme segment dans la 1ere partie, 
         ! on parcourt la deuxieme !
         IF(WhichSeg(2) ==  0) THEN
            incr = -incr
            jseg = WhichSeg(1) + incr
            ! recherche de ismax
            DO
               if (jseg<Idebut) exit  ! modif a verifier

               IF ((TabSeg(1,jseg,1).EQ.ismin).AND.(TabSeg(2,jseg,1) .EQ. ismax)) THEN
                  WhichSeg(2) = jseg
                  EXIT
               END IF
               
               jseg = jseg + incr
               IF ( jseg .LT. idebut .OR. jseg .GT. Ifin ) THEN
                  WhichSeg(2) =  0
                  EXIT
               END IF
            END DO
         ENDIF

      ENDIF
         
    END FUNCTION WhichSeg
    !-----------------------------------------------------------!

  END SUBROUTINE Segments2DP
  !---------------------------------------------------------------------------!


  !---------------------------------------------------------------------------!
  ! fonction qui renvoie le nombre de sommets voisins de chaque sommet
  ! et leurs numéros !
  SUBROUTINE compute_voisin(Seg, Var)

    USE var_types

    IMPLICIT NONE

    ! variables d'appel
    TYPE(MeshSeg),   INTENT(in)         :: Seg
    TYPE(Variables), INTENT(inout)      :: Var

    ! variables locales
    INTEGER    :: iseg, is1, is2, Nseg, Nsegtot, i

    Nsegtot   = Seg%Nsegmtot
    Nseg      = Seg%Nsegmt
 
    ALLOCATE(Var%NbVois( 1:Var%NCells ))
    Var%NbVois = 0

    DO iseg = 1,Nseg
       is1 = Seg%Nubo(1,iseg)
       is2 = Seg%Nubo(2,iseg)       

       i   = MAX(is1,is2)
       is1 = MIN(is1,is2)
       is2 = i

       Var%NbVois(is1) = Var%NbVois(is1)+1
       Var%NbVois(is2) = Var%NbVois(is2)+1
    END DO

    ALLOCATE( Var%NuVois(1:Var%NCells, 1:MAXVAL(Var%NbVois)) )

    Var%NuVois = 0        
    Var%NbVois = 0

    DO iseg = 1,Nseg
       is1 = Seg%Nubo(1,iseg)
       is2 = Seg%Nubo(2,iseg)

       i   = MAX(is1,is2)
       is1 = MIN(is1,is2)
       is2 = i

       Var%NbVois(is1) = Var%NbVois(is1)+1
       Var%NbVois(is2) = Var%NbVois(is2)+1

       Var%NuVois(is1,Var%NbVois(is1)) = is2
       Var%NuVois(is2,Var%NbVois(is2)) = is1
    END DO


  END SUBROUTINE compute_voisin
  !---------------------------------------------------------------------------!


END MODULE SegmtsGeom
  
           
           





