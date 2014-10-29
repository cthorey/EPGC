MODULE MODULE_THERMAL_NEWTON_OLD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ce module utilise résout l'équation de la chaleur en intégrant sur toute 
!!!!!!!!!! l'épaissuer contrziarment à l'approche utilisé ensuite pour copier
!!!!!!!!!! Balmforth. Ne pas oublie de changer la def de Xi dans le main.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  SUBROUTINE  THERMAL_NEWTON_SOLVER_OLD(Xi,H,P,T,BL,Dt,Dr,dist,ray,Mtot,sigma,nu,Pe,delta0,el,grav,theta,F_err)

    IMPLICIT NONE

    DOUBLE PRECISION , DIMENSION(:,:), INTENT(IN) :: H,P
    DOUBLE PRECISION , DIMENSION(:,:), INTENT(INOUT) :: Xi,T,BL
    DOUBLE PRECISION , INTENT(IN) :: Dt,Dr
    DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION , INTENT(IN) :: sigma,nu,Pe,delta0,el,grav,theta
    DOUBLE PRECISION, INTENT(INOUT) :: F_Err
    INTEGER, INTENT(IN) :: Mtot

    DOUBLE PRECISION :: Ite,U
    INTEGER :: ndyke,i,N
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: Xi_guess,Xi_tmps
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,d,e,f,g,k,l,S
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a1,b1,c1,d1,e1,f1,g1,k1,l1
    INTEGER :: err1,err2,col,Size
    DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: Xi_m,qa
    LOGICAL :: cho

    ! Taille de la grille
    ndyke=sigma/Dr
    CHO=COUNT(H(:,1)>delta0)<ndyke
    SELECT CASE (CHO)
    CASE(.TRUE.)
       N=ndyke +10   ! Cas ou on donne pas de profile initiale...
    CASE(.FALSE.)
       N=COUNT(H(:,1)>delta0)+10
    END SELECT

    ! Remplissage du flux au centre
    ALLOCATE(qa(1:N+1),stat=err2)
    IF (err2>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF
    DO i=1,N,1
       U=2.d0/(sigma)**4
       IF (i<ndyke+1) THEN
          qa(i)=U*(sigma**2-dist(i)**2)
       ELSE
          qa(i)=0.d0
       END IF
    END DO

    ! Remplissage di Xi_tmps et XI_gess
    ALLOCATE(Xi_tmps(1:N),Xi_guess(1:N),stat=err2)
    IF (err2>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF
    col=1
    CALL TEMPERATURE(Xi_tmps,col,N,Xi,H,T,BL,P,dist,ray,Dr,nu,Pe,qa,delta0,el,grav)
    col=2
    CALL TEMPERATURE(Xi_guess,col,N,Xi,H,T,BL,P,dist,ray,Dr,nu,Pe,qa,delta0,el,grav)

    ! Remplissage de la matrice jacobienne
    ALLOCATE(a1(1:N),b1(1:N),c1(1:N),d1(1:N),e1(1:N),stat=err1)
    ALLOCATE(f1(1:N),g1(1:N),k1(1:N),l1(1:N),stat=err2)
    IF (err1>1 .OR. err2>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF
    CALL JACOBI_TEMPERATURE(a1,b1,c1,d1,e1,f1,g1,k1,l1,N,H,T,P,BL,Dr,dist,ray,nu,Pe,qa,delta0,el,grav)

    ! Matrice à inverser
    ALLOCATE(a(1:N),b(1:N),c(1:N),d(1:N),e(1:N),stat=err1)
    ALLOCATE(f(1:N),g(1:N),k(1:N),l(1:N),S(1:N),stat=err2)
    IF (err1>1 .OR. err2>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF

    DO i=1,N,1
       a(i)=-theta*Dt*a1(i)
       b(i)=-theta*Dt*b1(i)
       c(i)=-theta*Dt*c1(i)
       d(i)=-theta*Dt*d1(i)
       e(i)=1.d0-theta*Dt*e1(i)
       f(i)=-theta*Dt*f1(i)
       g(i)=-theta*Dt*g1(i)
       k(i)=-theta*Dt*k1(i)
       l(i)=-theta*Dt*l1(i)
       S(i)=Xi(i,1)-Xi(i,2)+theta*Dt*Xi_guess(i)+(1-theta)*Dt*Xi_tmps(i)
    END DO

    a(1)=0;a(2)=0;a(3)=0;a(4)=0
    b(1)=0;b(2)=0;b(3)=0
    c(1)=0;c(2)=0
    d(1)=0
    l(N-3)=0;l(N-2)=0;l(N-1)=0;l(N)=0
    k(N-2)=0;k(N-1)=0;k(N)=0
    g(N-1)=0;g(N)=0
    f(N)=0

    ! Inversion de la matrice à l'aide de l'algo nona

    ALLOCATE(Xi_m(1:N),stat=err2)
    IF (err2>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF

    CALL NONA_DIAGO(N,Xi_m,a,b,c,d,e,f,g,k,l,S)

    DO i=1,N,1
       Xi(i,3)=Xi_m(i)+Xi(i,2)
    END DO
    
    ! Reconstruction des variable T et BL
    CALL Xi_Split(Xi,H,BL,T,N)

    ! Calcule de l'erreur
    IF (DOT_PRODUCT(Xi(:,2),Xi(:,2)) == 0D0) THEN
       F_err = ABS(MAXVAL(Xi_m(:)))
    ELSE
       Size = COUNT(Xi(:,2)>1D-10)
       F_err = ABS(MAXVAL(((Xi(:Size,3)-Xi(:Size,2))/Xi(:Size,2))))
    ENDIF

    ! Deallocation
    DEALLOCATE(Xi_m,qa,a,b,c,d,e,f,g,k,l,S)
    DEALLOCATE(Xi_guess,Xi_tmps,a1,b1,c1,d1,e1,f1,g1,k1,l1)

  END SUBROUTINE THERMAL_NEWTON_SOLVER_OLD

  SUBROUTINE  TEMPERATURE(Te,col,N,Xi,H,T,BL,P,dist,ray,Dr,nu,Pe,qa,delta0,el,grav)


    !### INTIALISATION DES VARIABLES ###!

    IMPLICIT NONE

    DOUBLE PRECISION  ,DIMENSION(:)  , INTENT(INOUT)          :: Te
    INTEGER           ,INTENT(IN)                             :: col,N
    DOUBLE PRECISION  ,DIMENSION(:,:), INTENT(IN)             :: Xi,H,T,BL,P
    DOUBLE PRECISION  ,INTENT(IN)                             :: Dr 
    DOUBLE PRECISION  ,DIMENSION(:), INTENT(IN)               :: qa,dist,ray
    DOUBLE PRECISION  ,INTENT(IN)                             :: nu,Pe,delta0,el,grav

    DOUBLE PRECISION                                          :: Ai,Bi
    DOUBLE PRECISION                                          :: h_a,h_a2,h_a3,delta_a,delta_a2,eta_a,zeta_a,Omega_a,T_a
    DOUBLE PRECISION                                          :: Psi_a,Sigma_a
    DOUBLE PRECISION                                          :: h_b,h_b2,h_b3,delta_b,delta_b2,eta_b,zeta_b,Omega_b,T_b
    DOUBLE PRECISION                                          :: Psi_b,Sigma_b
    DOUBLE PRECISION                                          :: loss
    INTEGER                                                   :: i


    !### REMPLISSAGE DE T ###!

    DO i=1,N,1
       IF1:IF  (i .NE. 1) THEN
          Bi=(ray(i-1)/(dist(i)*Dr))
          h_b=0.5d0*(H(i,3)+H(i-1,3))
          h_b2=0.5d0*(H(i,3)**2+H(i-1,3)**2)
          h_b3=0.5d0*(H(i,3)**3+H(i-1,3)**3)
          delta_b=0.5d0*(BL(i,col)+BL(i-1,col))
          delta_b2=0.5d0*(BL(i,col)**2+BL(i-1,col)**2)
          eta_b=(grav*(H(i,3)-H(i-1,3))+el*(P(i,3)-P(i-1,3)))/Dr
          T_b=0.5d0*(T(i,col)+T(i-1,col))
          zeta_b=h_b*(sqrt(5.d0)*h_b+sqrt(2.d0)*delta_b)/(sqrt(2.d0))
          Psi_b=(delta_b2/3.d0)*(2*h_b+(8.d0/35.d0)*delta_b)-2*h_b3
          Sigma_b=eta_b*(((sqrt(10.d0)-3.d0)/5.d0)*nu*h_b*zeta_b)*T_b
          Omega_b=eta_b*((3.d0)/(5.d0)*nu*(zeta_b-delta_b2))
       ENDIF IF1
       
       IF2:IF (i .NE. N) THEN
          Ai=(ray(i)/(dist(i)*Dr))
          h_a=0.5d0*(H(i+1,3)+H(i,3))
          h_a2=0.5d0*(H(i+1,3)**2+H(i,3)**2)
          h_a3=0.5d0*(H(i+1,3)**3+H(i,3)**3)
          delta_a=0.5d0*(BL(i+1,col)+BL(i,col))
          delta_a2=0.5d0*(BL(i+1,col)**2+BL(i,col)**2)
          eta_a=(grav*(H(i+1,3)-H(i,3))+el*(P(i+1,3)-P(i,3)))/Dr
          T_a=0.5d0*(T(i+1,col)+T(i,col))
          zeta_a=h_a*(sqrt(5.d0)*h_a+sqrt(2.d0)*delta_a)/(sqrt(2.d0))
          Psi_a=(delta_a2/3.d0)*(2*h_a+(8.d0/35.d0)*delta_a)-2*h_a3
          Sigma_a=eta_a*(((sqrt(10.d0)-3.d0)/5.d0)*nu*h_a*zeta_a)*T_a
          Omega_a=eta_a*((3.d0)/(5.d0)*nu*(zeta_a-delta_a2))            
       END IF IF2

       loss=(4*Pe*T(i,col))/BL(i,col)
       
       IF3: IF (i==1) THEN
          Te(i)=qa(i)-loss+Ai*Omega_a*Xi(i,col)+Ai*Sigma_a
       ELSEIF (i==N) THEN
          Te(i)=qa(i)-loss-Bi*Omega_b*Xi(i-1,col)-Bi*Sigma_b
        ELSE
           Te(i)=Ai*Omega_a*Xi(i,col)-Bi*Omega_b*Xi(i-1,col)&
               &+Ai*Sigma_a-Bi*Sigma_b &
               &-loss+qa(i)
          ! print*,i,(4*Phi*T(i,col))/BL(i,col)
        END IF IF3

    ENDDO
    
  END SUBROUTINE TEMPERATURE

 SUBROUTINE  JACOBI_TEMPERATURE(a,b,c,d,e,f,g,k,l,N,H,T,P,BL,Dr,dist,ray,nu,Pe,qa,delta0,el,grav)

    !### INTIALISATION DES VARIABLES ###!

    IMPLICIT NONE

    DOUBLE PRECISION  ,DIMENSION(:)  , INTENT(INOUT)          :: a,b,c,d,e,f,g,k,l
    INTEGER           ,INTENT(IN)                             :: N
    DOUBLE PRECISION  ,DIMENSION(:,:), INTENT(IN)             :: H,T,P,BL
    DOUBLE PRECISION  ,INTENT(IN)                             :: Dr 
    DOUBLE PRECISION  ,DIMENSION(:), INTENT(IN)               :: qa,dist,ray
    DOUBLE PRECISION  ,INTENT(IN)                             :: nu,Pe,delta0,el,grav

    DOUBLE PRECISION                                          :: Ai,Bi
    DOUBLE PRECISION                                          :: h_a,h_a2,delta_a,delta_a2,eta_a,zeta_a,Omega_a,T_a
    DOUBLE PRECISION                                          :: h_b,h_b2,delta_b,delta_b2,eta_b,zeta_b,Omega_b,T_b
    INTEGER                                                   :: i,col


    !### REMPLISSAGE DE LA MATRICE JACOBIENNE ###!
    
    col=2

    DO i=1,N,1
       IF1: IF (i .NE. 1) THEN
          Bi=(ray(i-1)/(dist(i)*Dr))
          h_b=0.5d0*(H(i,3)+H(i-1,3))
          h_b2=0.5d0*(H(i,3)**2+H(i-1,3)**2)
          delta_b=0.5d0*(BL(i,col)+BL(i-1,col))
          delta_b2=0.5d0*(BL(i,col)**2+BL(i-1,col)**2)
          eta_b=(grav*(H(i,3)-H(i-1,3))+el*(P(i,3)-P(i-1,3)))/Dr
          T_b=0.5d0*(T(i,col)+T(i-1,col))
          zeta_b=h_b*(sqrt(5.d0)*h_b+sqrt(2.d0)*delta_b)/(sqrt(2.d0))
          Omega_b=eta_b*((3.d0)/(5.d0)*nu*(zeta_b-delta_b2))
       ENDIF IF1
       IF2: IF (i .NE. N) THEN
          Ai=(ray(i)/(dist(i)*Dr))
          h_a=0.5d0*(H(i+1,3)+H(i,3))
          h_a2=0.5d0*(H(i+1,3)**2+H(i,3)**2)
          delta_a=0.5d0*(BL(i+1,col)+BL(i,col))
          delta_a2=0.5d0*(BL(i+1,col)**2+BL(i,col)**2)
          eta_a=(grav*(H(i+1,3)-H(i,3))+el*(P(i+1,3)-P(i,3)))/Dr
          T_a=0.5d0*(T(i+1,col)+T(i,col))
          zeta_a=h_a*(sqrt(5.d0)*h_a+sqrt(2.d0)*delta_a)/(sqrt(2.d0))
          Omega_a=eta_a*((3.d0)/(5.d0)*nu*(zeta_a-delta_a2))
       END IF IF2

       IF3:IF (i==1) THEN
          a(i)=0; b(i)=0; c(i)=0;d(i)=0
          e(i)=Ai*Omega_a
          f(i)=0;g(i)=0;k(i)=0;l(i)=0
       ELSEIF (i==N) THEN
          a(i)=0; b(i)=0; c(i)=0;
          d(i)=-Bi*Omega_b;
          e(i)=0.d0;
          f(i)=0;g(i)=0;k(i)=0;l(i)=0  
       ELSE
          a(i)=0.d0;b(i)=0.d0;c(i)=0.d0
          d(i)=-Bi*Omega_b
          e(i)=Ai*Omega_a
          f(i)=0;g(i)=0;k(i)=0;l(i)=0
       END IF IF3
    END DO

  END SUBROUTINE JACOBI_TEMPERATURE

  SUBROUTINE Xi_Split(Xi,H,BL,T,N)

    IMPLICIT NONE

    DOUBLE PRECISION  ,DIMENSION(:,:)  , INTENT(IN)             :: H
    DOUBLE PRECISION  ,DIMENSION(:,:)  , INTENT(INOUT)          :: T,BL,Xi
    INTEGER           ,INTENT(IN)                               :: N

    INTEGER                                                     :: i


!!! RECONSTRUCTION DU VECTEURS COUCHE LIMITES

!!! Test sur XI 
    DO i=1,N
       IF (Xi(i,3)>H(i,3) ) THEN
          print*,'erreur ds Xi',i,Xi(i,3),H(i,3)
          Xi(i,3)=H(i,3)
       ENDIF
    END DO

    DO i=1,N
       IF ((2.d0*H(i,3))/(3.d0)< Xi(i,3) .AND. Xi(i,3)<H(i,3)) THEN
          BL(i,3)=((3.d0)/(2.d0))*(H(i,3)-Xi(i,3))
          T(i,3)=1
       ELSEIF (Xi(i,3)<(2.d0/(3.d0))*H(i,3)) THEN
          BL(i,3)=H(i,3)/2
          T(i,3)=(3.d0*Xi(i,3))/(2*H(i,3))
       END IF
    END DO
    
  END SUBROUTINE Xi_Split


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

END MODULE MODULE_THERMAL_NEWTON_OLD
