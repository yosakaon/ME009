PROGRAM puissance_iteree
!	Calcul de valeurs propres
!       par la methode de puissance iteree et deflation
!       a compiler avec des reels 8 (double precision)
  USE v_propre
  IMPLICIT NONE

  INTEGER			:: n

  REAL, DIMENSION (:, :), ALLOCATABLE	:: A, L, U
  REAL, DIMENSION (:), ALLOCATABLE	    :: vp
  INTEGER			:: i,j
  REAL				:: lambda
  PRINT *, "Entrez n"
  READ *, n
  ALLOCATE (A(n,n))
  ALLOCATE (L(n,n))
  ALLOCATE (U(n,n))
  allocate(vp(n))
  A=0;
   do j=1, n
     do i=1, n
       if (i==j) then
       A(i,j) = 4;
     end if
     end do
   end do
   A(1,2) = -1;
   A(2,1)=-1;
   A(2,3)=-1;
   A(3,2)=-1;
   A(3,4)=-1;
   A(4,3)=-5;
  ! Determination des 1ers vecteur/valeur propre vp/lambda
  i = 1
   CALL puis_iter(A, vp, lambda)
   !MVP (:,i) = vp		! stockage des vecteurs propres
   !valp (i) = lambda		! et des valeurs propres associees
  PRINT *,'Valeur propre no ', i,' = ', lambda
  PRINT *,'Vecteur propre no ', i,' : ', vp



  ! Determination des n-1 vecteurs/valeurs propres restants
  ! en utilisant la méthode de la deflation
  DO i=2, n
  call deflation(A,vp,lambda)
  call puis_iter(A,vp,lambda)
 ! MVP (:,i) = vp		! stockage des vecteurs propres
  !valp (i) = lambda		! et des valeurs propres associees
  PRINT *,'Valeur propre no ', i,' = ', lambda
  PRINT *,'Vecteur propre no ', i,' : ', vp
  END DO



  ! PUISSANCE ITEREE INVERSE
  PRINT *
  PRINT *,'Puissance iteree inverse'
  ! Initialisation de A ici matrice de Vandermonde A FAIRE

  call decomp(A,L,U)
  ! Puissance iteree inverse
  CALL puis_iter_inv(L, U, vp, lambda) ! A FAIRE
  PRINT *,'Valeur propre no ', n,' = ', lambda
  PRINT *,'Vecteur propre no ', n,' : ', vp

  DEALLOCATE (A)
  DEALLOCATE (vp)

END PROGRAM puissance_iteree
