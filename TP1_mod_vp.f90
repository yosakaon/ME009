MODULE v_propre
  IMPLICIT NONE
  REAL, PARAMETER 			:: eps = 1.e-6

  CONTAINS

  ! Procedure de factorisation A=LU
  SUBROUTINE decomp(A,L,U)
    REAL, DIMENSION (:, :), INTENT(IN)    :: A
    REAL, DIMENSION (:, :), INTENT(OUT)   :: L, U
    INTEGER         :: i, j, k, n
    REAL            :: s

    n = SIZE(A,1)
    DO j=1,n
        L(j,j)=1.
        DO i=1,j
            s=0
            DO k=1,i-1
                s=s+L(i,k)*U(k,j)
            END DO
            U(i,j)=A(i,j)-s
        END DO
        DO i=j+1,n
            s=0
            DO k=1,j-1
                s=s+L(i,k)*U(k,j)
            END DO
            L(i,j)=(A(i,j)-s)/U(j,j)
        END DO
    END DO
  END SUBROUTINE decomp

  SUBROUTINE descente(L,y,b)
    REAL, DIMENSION (:, :), INTENT(IN)   :: L
    REAL, DIMENSION (:), INTENT(IN)      :: b
    REAL, DIMENSION (:), INTENT(OUT)     :: y
    INTEGER         :: i, j, n
    REAL            :: s

    n = SIZE(L,1)
    DO i = 1,n
         s = 0.
         DO j = 1, i-1
              s = s + L(i,j) * y(j)
         END DO
         y(i)= (b(i) - s) / L(i,i)
    END DO
  END SUBROUTINE descente

  SUBROUTINE remonte(U,x,y)
  REAL, DIMENSION (:, :), INTENT(IN)    :: U
  REAL, DIMENSION (:), INTENT(IN)       :: y
  REAL, DIMENSION (:), INTENT(OUT)      :: x
  INTEGER         :: i, j, n
  REAL            :: s

  n = SIZE(U,1)
  DO i = n, 1, -1
       s = 0.
       DO j= i+1, n
           s = s + U(i,j) * x(j)
       END DO
       x(i)= (y(i) - s) / U(i,i)
  END DO
  END SUBROUTINE remonte

  ! Affichage d'une matrice A de dimensions nl x nc
  SUBROUTINE affiche (M)
    REAL, DIMENSION (:, :), INTENT(IN)    :: M
    INTEGER         :: i, j, nl, nc

    nl = SIZE(M,1)
    nc = SIZE(M,2)
    DO i=1, nl
        WRITE (*,fmt='(30E12.4)') (M(i,j), j=1, nc)
    END DO
    PRINT *
  END SUBROUTINE affiche

  ! Affichage sur une ligne d'un vecteur x de dimension n
  SUBROUTINE affiche_vec (x)
    REAL, DIMENSION (:), INTENT(IN)    :: x
    INTEGER         :: i,n

    n = SIZE(x)
    WRITE (*,fmt='(30E12.4)') (x(i), i=1, n)
    PRINT *
  END SUBROUTINE affiche_vec

  ! Puissance iteree sur la matrice M
  SUBROUTINE puis_iter(M, vp, lambda)

    REAL, DIMENSION (:,:), INTENT(IN) 	:: M
    REAL, DIMENSION (:), INTENT(OUT)   	:: vp
    REAL, INTENT(INOUT)                 	:: lambda
    REAL, DIMENSION (:), ALLOCATABLE	:: u
    REAL 				:: anclambda
    INTEGER 			:: n, i
    real :: norme
    integer :: k
    integer :: kmax = 100
    lambda = 0.
    n = SIZE(M,1)
    ALLOCATE (u(n))
    CALL RANDOM_NUMBER(u)
    vp = 1.
    vp = vp/SQRT(SUM(vp*vp)) ! normalisation
    do k=2, kmax
    u = MATMUL(M,vp)
    anclambda = lambda
    lambda = dot_product(vp,u)
    norme = sqrt(sum(u**2))
    vp = u / norme
    if ((abs(lambda-anclambda)) < eps) exit
  end do
    DEALLOCATE (u)
    print*, "nombre d'itérations", k
  END SUBROUTINE puis_iter

  ! Calcul de la nouvelle matrice M pour la deflation a l'aide
  ! des vecteur/valeur propre vp/lambda
  SUBROUTINE deflation(M, vp, lambda)

    REAL, DIMENSION (:,:), INTENT(INOUT) :: M
    REAL, DIMENSION (:), INTENT(IN)   	 :: vp
    REAL, INTENT(IN)                 	 :: lambda
    real, allocatable, dimension(:) :: y
    real, allocatable, dimension(:,:) :: At
    REAL :: lambdat = 0.
    integer :: i,j,n
    n = SIZE(M,1)
    allocate(At(n,n))
    allocate(y(n))
    At = transpose(M)
    call puis_iter(At,y,lambdat)
    do i = 1, n
      do j=1, n
    M(i,j) = M(i,j) - lambda*((vp(i)*y(j))/dot_product(vp,y))
  end do
end do
  END SUBROUTINE deflation



  ! Puissance iteree inverse sur une matrice
  ! passee sous forme de sa decomposition L U
  SUBROUTINE puis_iter_inv(L, U, vp, lambda)
    REAL, DIMENSION(:,:), INTENT(IN) 	:: L, U
    REAL, DIMENSION(:), INTENT(OUT)   	:: vp
    REAL, INTENT(OUT)                 	:: lambda
    real, allocatable, dimension(:) :: y
    real, allocatable, dimension(:) :: uk
    real :: mu, norme, ancmu
    integer :: n, k, kmax
    kmax = 1000
    n = size(U,1)
    allocate(y(n))
    allocate(uk(n))
    vp = 1.
    vp = vp/SQRT(SUM(vp*vp)) ! normalisation
    ! AJOUTER LES DECLARATIONS NECESSAIRES
    do k=1, kmax
    call descente(L,y,vp)
    call remonte(U,uk,y)
    mu = dot_product(vp,uk)
    ancmu = mu
    norme = sqrt(sum(uk**2))
    vp = uk / norme
    if ((abs(mu-ancmu)) < eps) exit
    end do
    DEALLOCATE (uk)

  END SUBROUTINE puis_iter_inv



END MODULE v_propre
