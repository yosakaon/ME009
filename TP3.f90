PROGRAM sor_points

IMPLICIT NONE

! DECLARATIONS DES VARIABLES (A COMPLETER)
! PARAMETRES DE LA RESOLUTION
INTEGER, PARAMETER  :: nx=23, ny=19, nn=nx*ny, kmax=5000, kech = 10
REAL, PARAMETER     :: h=1.5/(nx+1), eps=1.e-8, omega=1.97
INTEGER             :: ix, iy, i, k, mx, my, ij, j
real :: gs, gn, gw, norme_res, somme, QH1, QV1, QH2,QV2,QH3,QV3

! VECTEURS QUI COMPOSENT LA MATRICE A DU SYSTEME (on ne travaille qu'avec ces vecteurs, pas avec A, voir TD2 !) ET SECOND MEMBRE b
REAL, DIMENSION(nn) :: As, Aw, Ap, Ae, An, b, res
! VECTEUR UTILISES POUR LA RESOLUTION SOR
REAL, DIMENSION(-nx+1:nn+nx) :: Q, dQ
! CHAMP DE TEMPERATURE 2D (pour tracé isovaleurs, lors du post-traitement !)
REAL, DIMENSION(0:nx+1,0:ny+1) :: T
! VECTEUR x
REAL, DIMENSION(0:nx+1) :: x
! VECTEUR y
REAL, DIMENSION(0:ny+1) :: y, ge


! OUVERTURE DES FICHIERS POUR ENREGISTREMENT DES RESULTATS
OPEN(10,file='isoT_sor.dat')
OPEN(11,file='res_sor.out')
PRINT *,' nx ny',nx,ny

! REMPLISSAGE DES VECTEURS As, Aw, Ap, Ae, An, et b, voir TD2
!
Ap=4;
As=-1;
Aw=-1;
Ae=-1;
An=-1;
b=0;
gw=0.;
do i=1, ny+1
ge=h*ny;
end do
do i=1, nx
  As(i)=-1;
  b(i)=b(i);
end do
do i=1, ny
  Aw(1+(i-1)*nx) = 0;
  b(1+(i-1)*nx) = b(1+(i-1)*nx) + gw;
end do
do i = 1, ny
  Ae(i*nx)=0;
  b(i*nx)=b(i*nx)+ge(i*h);
end do
do i=1, nx
  An(i+(ny-1)*nx)=-1;
  b(i+(ny-1)*nx)=b(i+(ny-1)*nx);
end do
! AFFICHAGE (POUR VERIFICATION) SEULEMENT SI nx = 5 et ny = 4
IF(nx==5) THEN

PRINT *,' vecteur As', As

PRINT *,' vecteur Aw', Aw

PRINT *,' vecteur Ap', Ap

PRINT *,' vecteur Ae', Ae

PRINT *,' vecteur An', An

PRINT *,' vecteur b' , b
ENDIF

! RESOLUTION PAR S.O.R (voir TD3)
! INITIALISATIONS A FAIRE AVANT ENTREE DANS LA BOUCLE
Q=0.;
dQ=0.;
! DEBUT BOUCLE CONDITIONNELLE
! A éCRIRE
do i=1, nx
  Q(i+nx*ny)=-Q(i+nx*ny-nx); ! nord
  Q(i)=0;
end do
do k=1, kmax
  do i=1, nn
    res(i)=b(i)-As(i)*Q(i-nx)-Aw(i)*Q(i-1)-Ap(i)*Q(i)-Ae(i)*Q(i+1)-An(i)*Q(i+nx);
    !print*, res
  end do
  norme_res=sqrt(sum(res**2));
  print*, norme_res;
  if ((norme_res)<eps) exit
    do i=1, nn
      dQ(i)=(omega/Ap(i))*(res(i)-As(i)*dQ(i-nx)-Aw(i)*dQ(i-1));
  end do
  Q=Q+dQ;

! POST-TRAITEMENT
! ECRITURE DES RESULTATS DANS UN FICHIER
! ET AFFICHAGE EVENTUEL A l'ECRAN POUR CONTROLE DES ITERATIONS (évolution du residu)
!IF(MOD(k,kech)==0) THEN
WRITE(11,fmt='(I6,E14.5)') k, norme_res ! norm_res (réel) désigne ici la norme du résidu
									   ! k (entier) désigne ici le nombre d'itérations
!ENDIF
end do


! FIN BOUCLE CONDITIONNELLE
! FIN RESOLUTION S.O.R

PRINT *,'convergence en',k,'iterations'

WRITE(11,fmt='(I6,E14.5)') k, norme_res ! norm_res désigne ici la norme du résidu

! POST-TRAITEMENT : ISOTHERMES DE TEMPERATURE ET CALCUL DES FLUX QH ET QV
! APPEL DE LA PROCEDURE POUR PASSER DU VECTEUR SOLUTION Q AU CHAMP 2D T(x_i,y_j), voir TD3
CALL ecritureT(Q, T, x, y, nx, ny, nn, h)
!
!

! CALCUL DES FLUX QH et QV SELON LES SECTIONS MEDIANES (voir TD3)

mx=(nx+1)/2;
my=(ny+1)/2;
somme=0.;
do j=0, ny
somme=somme+(-T(mx+1,j)-T(mx+1,j+1)+T(mx-1,j)+T(mx-1,j+1));
end do
if (nx==5) then
QH1=1/4*somme;
else if (nx==95) then
QH2=1/4*somme;
else
  QH3=1/4*somme;
end if
somme=0.;
do i=0, nx
  somme=somme+(-T(i,my+1)-T(i+1,my+1)+T(i,my-1)+T(i+1,my-1))
end do
if (nx==5) then
QV1=1/4*somme;
else if (nx==95) then
QV2=1/4*somme;
else
QV3=1/4*somme;
end if

! etude de convergence en maillage à faire

END PROGRAM sor_points

SUBROUTINE ecritureT(Q, T, x, y, nx, ny, nn, h)
! DECLARATION DES VARIABLES A COMPLETER
INTEGER :: nx, ny, nn, ix, iy
REAL, DIMENSION(-nx+1:nn+nx) :: Q
REAL, DIMENSION(0:nx+1,0:ny+1) :: T
REAL, DIMENSION(0:nx+1) :: x
! VECTEUR y
REAL, DIMENSION(0:ny+1) :: y

ix=h*nx;
iy=h*ny;
! initialiser x et y
do i=0, nx
  x(i)=i*h
end do
do j=0, ny
  y(j)=j*h
end do
!  conditions limites sur T
do i=0, nx
  T(i,0)=gs;
  T(i,ny+1)=gn;
end do
do j=0, ny
  T(0,j)=gw;
  T(nx+1,j)=ge;
end do

do i=1, nx
  do j=1, ny
    T(i,j)=Q(i+(j-1)*nx);
  end do
end do
do i=0, nx
  do j=0, ny
WRITE(10,fmt='(E14.5,E14.5,E14.5)') x(i),y(j),T(i,j);

end do
    write(10,fmt=*) CHAR(13);
end do
END SUBROUTINE

SUBROUTINE convergenceMaillage(QH1, QV1, QH2, QV2, QH3, QV3)
end subroutine
