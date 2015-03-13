MODULE MODULE_THICKNESS_INTE_GFD_BERCOVICI

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINE THICKNESS_NEWTON_SOLVER

  SUBROUTINE  THICKNESS_INTE_GFD_BERCOVICI(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,z,F_err,theta)

    !*****************************************************************
    !Solve for the thickness in the thickenss evolution equation using the Newton
    ! method
    !*****************************************************************
    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION , DIMENSION(:,:), INTENT(INOUT) :: H,P
    DOUBLE PRECISION, DIMENSIOn(:,:), INTENT(IN) :: T,BL,Ts
    DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: dist,ray

    !Parametre du model
    DOUBLE PRECISION , INTENT(IN) :: Dt,Dr,theta

    !Nombre sans dimensions
    DOUBLE PRECISION , INTENT(IN) :: el,grav,sigma,nu,delta0
    INTEGER, INTENT(IN) :: M, z
    DOUBLE PRECISION , INTENT(INOUT) :: F_err

    !Variable du sous programmes
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,S
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Hm,qa

    DOUBLE PRECISION :: U
    INTEGER :: ndyke,i,N,Size,code
    INTEGER :: err1,col,algo1
    LOGICAL :: cho

    ! Taille de la grille sur laquelle on fait l'inversion
    ndyke=sigma/Dr
    CHO=COUNT(H(:,1)>0.D0)<ndyke
    SELECT CASE (CHO)
    CASE(.TRUE.)
       N = ndyke+3
    CASE(.FALSE.)
       N = COUNT(H(:,1)>0.d0)
    END SELECT

    ! Caracterisation du flux
    ALLOCATE(qa(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur alloc flux';STOP
    END IF

    DO i = 1,N,1
       U = 2.d0/(sigma)**4.
       IF (i<ndyke+1) THEN
          qa(i) = U*(sigma**2.-dist(i)**2.)
       ELSE 
          qa(i) = 0.d0
       END IF
    END DO


    ! Remplissage matrix
    ALLOCATE(a(1:N),b(1:N),c(1:N),S(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur allocation dans coeff Temperature'; STOP
    END IF

    CALL THICKNESS_MATRIX_FILL(a,b,c,S,N,M,H,P,T,BL,Ts,Dt,Dr,dist,ray,el,grav,nu,delta0,qa)

    a(1)=0
    c(N)=0

    !Inversion de la matrice
    ALLOCATE(Hm(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF
    CALL TRIDIAG(a,b,c,S,N,Hm)
    DO i=1,N,1
       H(i,3)=Hm(i)
    END DO 

    H(N+1:M,3) = delta0

    ! Calcul du soeuil F_err
    IF (DOT_PRODUCT(H(:,2),H(:,2)) == 0D0) THEN
       F_err = ABS(MAXVAL((Hm(:))))
    ELSE
       Size = COUNT(H(:,2)>1D-6)
       F_err = ABS(MAXVAL((H(:Size,3)-H(:Size,2))/H(:Size,2)))
    ENDIF

    DEALLOCATE(Hm,a,b,c,S,qa)

  END SUBROUTINE THICKNESS_INTE_GFD_BERCOVICI
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE THICKNESS
  !-------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE THICKNESS_MATRIX_FILL(a,b,c,S,N,M,H,P,T,BL,Ts,Dt,Dr,dist,ray,el,grav,nu,delta0,qa)

    !*****************************************************************
    ! Give the jacobian coeficient a1,b1,c1
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: a,b,c,S
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,BL,Ts
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray,qa

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr 
    INTEGER ,INTENT(IN) :: N,M

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,nu,delta0

    ! Parametre pour le sous programme
    DOUBLE PRECISION :: Agrav,h_a3,Bgrav,h_b3
    INTEGER :: i,col,algo1,err1

    col =2
    ! Remplissage de la matrice Jacobienne
    DO i=1,N,1
       IF1: IF (i .NE. N) THEN
          Agrav = grav*Dt*(ray(i)/(dist(i)*Dr**2))
          ! IF (H(i,col) == 0 .AND. H(i+1,col) ==0) THEN
          !    h_a3 =0.d0
          ! ENDIF

          h_a3 = (0.5d0*(H(i,col)+H(i+1,col)))**3
          !h_a3 =0.5*(H(i+1,col)**3+H(i,col)**3)
          
       ENDIF IF1

       IF2: IF (i .NE. 1) THEN
          Bgrav=Dt*grav*(ray(i-1)/(dist(i)*Dr**2))
          ! IF (H(i,col) == 0 .AND. H(i-1,col) ==0) THEN
          !    h_a3 =0.d0
          ! ENDIF
          ! h_b3 = (H(i-1,col)+H(i,col))**3
          h_b3 = (0.5d0*(H(i,col)+H(i-1,col)))**3
       ENDIF IF2


       IF3:IF (i==1) THEN
          a(i) = 0.d0
          b(i) = 1.d0+Agrav*h_a3
          c(i) = -Agrav*h_a3
          S(i) = H(i,1)+qa(i)*Dt

       ELSEIF (i == N) THEN
          a(i) = -Bgrav*h_b3
          b(i) = 1.d0+Bgrav*h_b3
          c(i) = 0.d0
          S(i) = H(i,1)+qa(i)*Dt
       ELSE
          a(i) = -Bgrav*h_b3
          b(i) = 1.d0+Bgrav*h_b3+Agrav*h_a3
          c(i) = -Agrav*h_a3
          S(i) = H(i,1)+qa(i)*Dt
       END IF IF3

    END DO

  END SUBROUTINE THICKNESS_MATRIX_FILL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE NONA DIAG
  !-------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  NONA_DIAGO(N,Hm,a,b,c,d,e,f,g,k,l,S)

    !*****************************************************************
    ! Solves for a vector Hm of length N the nano diagonal linear set
    ! M Hm = S, where A, B, C, D, E, F, G, K and  L  are the three main 
    ! diagonals of matrix M(N,N), the other terms are 0.
    ! S is the right side vector.
    !*****************************************************************
    IMPLICIT NONE

    INTEGER , INTENT(IN) :: N
    INTEGER :: i
    INTEGER :: err3,err4
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: a,b,c,d,e,f,g,k,l,S
    DOUBLE PRECISION, DIMENSION(:),INTENT(INOUT) :: Hm
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: zeta,alpha,beta,mu,xi,lambda,eta,omega,gamma

    AllOCATE(zeta(1:N),alpha(1:N),beta(1:N),mu(1:N),xi(1:N),stat=err3)
    ALLOCATE(lambda(1:N),eta(1:N),omega(1:N),gamma(1:N),stat=err4)

    IF (err3>1 .OR. err4>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur P,Q'; STOP
    END IF


    zeta(1)=b(1)
    alpha(1)=c(1)
    beta(1)=d(1)
    mu(1)=e(1)
    xi(1)=f(1)/mu(1)
    lambda(1)=g(1)/mu(1)
    eta(1)=k(1)/mu(1)
    omega(1)=l(1)/mu(1)
    gamma(1)=S(1)/mu(1)

    zeta(2)=b(2)
    alpha(2)=c(2)
    beta(2)=d(2)
    mu(2)=e(2)-xi(1)*beta(2)
    xi(2)=(f(2)-lambda(1)*beta(2))/mu(2)
    lambda(2)=(g(2)-eta(1)*beta(2))/mu(2)
    eta(2)=(k(2)-omega(1)*beta(2))/mu(2)
    omega(2)=l(2)/mu(2)
    gamma(2)=(S(2)-beta(2)*gamma(1))/mu(2)

    zeta(3)=b(3)
    alpha(3)=c(3)
    beta(3)=d(3)-xi(1)*alpha(3)
    mu(3)=e(3)-lambda(1)*alpha(3)-xi(2)*beta(3)
    xi(3)=(f(3)-eta(1)*alpha(3)-lambda(2)*beta(3))/mu(3)
    lambda(3)=(g(3)-omega(1)*alpha(3)-eta(2)*beta(3))/mu(3)
    eta(3)=(k(3)-omega(2)*beta(3))/mu(3)
    omega(3)=l(3)/mu(3)
    gamma(3)=(S(3)-alpha(3)*gamma(1)-beta(3)*gamma(2))/mu(3)

    zeta(4)=b(4)
    alpha(4)=c(4)-xi(1)*zeta(4)
    beta(4)=d(4)-lambda(1)*zeta(4)-xi(2)*alpha(4)
    mu(4)=e(4)-eta(1)*zeta(4)-lambda(2)*alpha(4)-xi(3)*beta(4)
    xi(4)=(f(4)-omega(1)*zeta(4)-eta(2)*alpha(4)-lambda(3)*beta(4))/mu(4)
    lambda(4)=(g(4)-omega(2)*alpha(4)-eta(3)*beta(4))/mu(4)
    eta(4)=(k(4)-omega(3)*beta(4))/mu(4)
    omega(4)=l(4)/mu(4)
    gamma(4)=(S(4)-zeta(4)*gamma(1)-alpha(4)*gamma(2)-beta(4)*gamma(3))/mu(4)

    DO i=5,N

       zeta(i)=b(i)-a(i)*xi(i-4)
       alpha(i)=c(i)-a(i)*lambda(i-4)-xi(i-3)*zeta(i)
       beta(i)=d(i)-a(i)*eta(i-4)-lambda(i-3)*zeta(i)-alpha(i)*xi(i-2)
       mu(i)=e(i)-a(i)*omega(i-4)-zeta(i)*eta(i-3)-lambda(i-2)*alpha(i)-beta(i)*xi(i-1)
       xi(i)=(f(i)-omega(i-3)*zeta(i)-eta(i-2)*alpha(i)-lambda(i-1)*beta(i))/mu(i)
       lambda(i)=(g(i)-alpha(i)*omega(i-2)-eta(i-1)*beta(i))/mu(i)
       eta(i)=(k(i)-omega(i-1)*beta(i))/mu(i)
       omega(i)=l(i)/mu(i)
       gamma(i)=(S(i)-a(i)*gamma(i-4)-zeta(i)*gamma(i-3)-alpha(i)*gamma(i-2)-beta(i)*gamma(i-1))/mu(i)

    END DO

    Hm(N)=gamma(N)
    Hm(N-1)=gamma(N-1)-xi(N-1)*Hm(N)
    Hm(N-2)=gamma(N-2)-lambda(N-2)*Hm(N)-xi(N-2)*Hm(N-1)
    Hm(N-3)=gamma(N-3)-eta(N-3)*Hm(N)-lambda(N-3)*Hm(N-1)-xi(N-3)*Hm(N-2)

    DO i=N-4,1,-1
       Hm(i)=gamma(i)-xi(i)*Hm(i+1)-lambda(i)*Hm(i+2)-eta(i)*Hm(i+3)-omega(i)*Hm(i+4)
    END DO


    DEALLOCATE(zeta,alpha,beta,mu,xi,lambda,eta,omega,gamma)

  END SUBROUTINE NONA_DIAGO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE TRIDIAG
  !-------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TRIDIAG(A,B,C,S,N,U)
    !*****************************************************************
    ! Solves for a vector U of length N the tridiagonal linear set
    ! M U = R, where A, B and C are the three main diagonals of matrix
    ! M(N,N), the other terms are 0. R is the right side vector.
    !*****************************************************************
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: A,B,C,S
    DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: U
    INTEGER, INTENT(IN) :: N
    
    INTEGER :: CODE
    DOUBLE PRECISION, DIMENSION(N) :: GAM
    DOUBLE PRECISION :: BET
    INTEGER :: j

    BET = B(1)

    IF (BET == 0.D0) THEN
       PRINT*,'ERROR TRIDIAG'
       STOP
    ENDIF
    U(1) = S(1)/BET
    DO J=2,N                    !Decomposition and forward substitution
       GAM(j)=C(j-1)/BET
       BET=B(j)-A(j)*GAM(j)

       IF(BET.EQ.0.D0) THEN            !Algorithm fails
          PRINT*,'ERRORTRIDIAG2',j,N
          STOP
       END IF
       U(j)=(S(j)-A(j)*U(j-1))/BET
    END DO

    DO j=N-1,1,-1                     !Back substitution
       U(j)=U(j)-GAM(j+1)*U(j+1)
    END DO

    CODE=0
    RETURN
  END SUBROUTINE TRIDIAG


END MODULE MODULE_THICKNESS_INTE_GFD_BERCOVICI

