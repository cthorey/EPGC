
MODULE MODULE_THERMAL_SKIN_NEWTON

  USE MODULE_MOBILITY
  
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Equation : Heat equation
!!!!!!!!  Model : Skin theory (Balmforth2004)
!!!!!!!!  Schema : Newton_Rhapsod method
!!!!!!!!  Model couplage : Bercovici rheology (1995) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE THERMAL_SKIN_NEWTON(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,&
       &M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps,Rheology,ERROR_CODE)

    !*****************************************************************
    ! Solve for the parameter Xi, and split in Temperature and thermal layer
    ! from   evolution equation using the Newton
    ! method
    !*****************************************************************
    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: H,P
    DOUBLE PRECISION , DIMENSION(:,:), INTENT(INOUT) :: Xi,T,BL,Ts
    DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: dist,ray

    !Parametre du model
    DOUBLE PRECISION , INTENT(IN) :: Dt,Dr,theta,tmps

    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology
    
    !Nombre sans dimensions
    DOUBLE PRECISION , INTENT(IN) :: sigma,nu,Pe,psi,delta0,el,grav,N1
    INTEGER, INTENT(IN) :: M, z
    DOUBLE PRECISION , INTENT(INOUT) :: F_err

    !Variable du sous programmes
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: Xi_guess,Xi_tmps
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,S
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a1,b1,c1
    DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: Xi_m

    DOUBLE PRECISION :: U
    INTEGER :: i,ndyke,N,Size
    INTEGER :: err1,col
    LOGICAL :: CHO

    ndyke=sigma/Dr
    DO i =1,M,1
       IF (H(i,1)<delta0) THEN
          N = i-1;EXIT
       ENDIF
    ENDDO
    N=M-1
    ! Calcule de f tmps n et n+
    ALLOCATE(Xi_tmps(1:N),Xi_guess(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur alloc xi_tmps-xi_guess'; STOP
    END IF

    col=1
    CALL TEMPERATURE_BALMFORTH(Xi_tmps,col,N,Xi,H,T,Ts,BL,P,dist,ray,Dr,nu,Pe,&
         &delta0,el,grav,N1,tmps+Dt,psi,Dt,Rheology,ERROR_CODE)
    col=2
    CALL TEMPERATURE_BALMFORTH(Xi_guess,col,N,Xi,H,T,Ts,BL,P,dist,ray,Dr,nu,Pe,&
         &delta0,el,grav,N1,tmps+Dt,psi,Dt,Rheology,ERROR_CODE)

    ! Jacobienner
    ALLOCATE(a1(1:N),b1(1:N),c1(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur allocation dans coeff Temperature'; STOP
    END IF

    CALL JACOBI_TEMPERATURE_BALMFORTH(a1,b1,c1,N,H,BL,T,Ts,Xi,P,Dr,Dt,dist,ray,&
         &nu,Pe,delta0,el,grav,Rheology,ERROR_CODE)

    !Systeme a inverser
    ALLOCATE(a(1:N),b(1:N),c(1:N),S(1:N),stat= err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF

    DO i=1,N,1
       a(i)=-theta*Dt*a1(i)
       b(i)=1D0-theta*Dt*b1(i)
       c(i)=-theta*Dt*c1(i)
       S(i)=(Xi(i,1)-Xi(i,2))+theta*Dt*Xi_guess(i)+(1-theta)*Dt*Xi_tmps(i)
    END DO

    a(1)=0
    c(N)=0

    !Inversion de la matrice
    ALLOCATE(Xi_m(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF

    CALL TRIDIAG(a,b,c,S,N,Xi_m)

    DO i=1,N,1
       Xi(i,3)=Xi_m(i)+Xi(i,2)
       IF (Xi(i,3) >H(i,3)/2D0) THEN
          Xi(i:,3) = H(i:,3)/2D0
          PRINT*,'Xi trop grand'
          EXIT
       ENDIF
    END DO

    ! Separation variables
    ! CALL XI_SPLIT(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)
    CALL XI_SPLIT_BALMFORTH(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)

    ! Calcule de l'erreur
    IF (DOT_PRODUCT(Xi(:,3),Xi(:,3))==DOT_PRODUCT(H(:,3)/2D0,H(:,3)/2D0)) THEN
       F_err = 0D0 ! Cas ou le refroidissemnt est trop important et tout devient nulle
    ELSEIF (DOT_PRODUCT(Xi(:,2),Xi(:,2)) == 0D0) THEN
       F_err = ABS(MAXVAL(Xi_m(:)))
    ELSE
       Size = COUNT(Xi(:,2)>1D-10)
       F_err = ABS(MAXVAL(((Xi(:Size,3)-Xi(:Size,2))/Xi(:Size,2))))
    ENDIF

    DEALLOCATE(Xi_m,a,b,c,S)
    DEALLOCATE(Xi_guess,Xi_tmps,a1,b1,c1)

  END SUBROUTINE THERMAL_SKIN_NEWTON

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


  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE TRIDIAG
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE TRIDIAG(A,B,C,S,N,U)
    !*****************************************************************
    ! Solves for a vector U of length N the tridiagonal linear set
    ! M U = R, where A, B and C are the three main diagonals of matrix
    ! M(N,N), the other terms are 0. R is the right side vector.
    !*****************************************************************
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: A,B,C,S
    DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: U
    INTEGER, INTENT(IN) :: N

    INTEGER :: CODE
    DOUBLE PRECISION, DIMENSION(N) :: GAM
    DOUBLE PRECISION :: BET
    INTEGER :: j

    IF(B(1) .EQ. 0.D0) THEN
       CODE=1
       RETURN
    END IF

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
          PRINT*,'ERRORTRIDIAG'
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


  SUBROUTINE XI_SPLIT(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)

    !*****************************************************************
    ! Solve for T and BL from xi for balmofrth with Ts=0
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: BL,T,Xi,Ts

    !Dimensionless parameter
    DOUBLE PRECISION, INTENT(IN) :: delta0,N1,Pe,el

    !Parametre du model
    DOUBLE PRECISION :: Dt,tmps
    INTEGER, INTENT(IN) :: N

    ! Parametre pour le sous programme
    INTEGER :: i
    DOUBLE PRECISION, PARAMETER :: pi=3.14159265
    DOUBLE PRECISION :: Xit,Tss,beta,d1,d2,Dr

    ! Separation des variables
    DO i=1,N
       Xit = H(i,3)/6.d0
       IF (Xi(i,3) <= Xit) THEN
          T(i,3) = 1.d0
          BL(i,3) = 3.d0*Xi(i,3)
       ELSEIF (Xi(i,3)> Xit) THEN
          T(i,3)= 3.d0/2.D0-(3.d0*Xi(i,3)/H(i,3))
          BL(i,3)=H(i,3)/2.d0
       ENDIF
    END DO
  END SUBROUTINE XI_SPLIT

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE X_SPLIT_Balmforth
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE XI_SPLIT_BALMFORTH(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)

    !*****************************************************************
    ! Solve for T and BL from Xi deriving Ts directly here
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: BL,T,Xi,Ts

    !Dimensionless parameter
    DOUBLE PRECISION, INTENT(IN) :: delta0,N1,Pe,el

    !Parametre du model
    DOUBLE PRECISION :: Dt,tmps
    INTEGER, INTENT(IN) :: N

    ! Parametre pour le sous programme
    INTEGER :: i
    DOUBLE PRECISION, PARAMETER :: pi=3.14159265
    DOUBLE PRECISION :: Xit,Tss,beta

    ! Separation des variables
    DO i=1,N
       beta = N1*Pe**(-0.5d0)/(sqrt(pi*(tmps+Dt)))
       Xit = beta*H(i,3)**2/(6.d0*beta*H(i,3)+24.d0)

       IF (Xi(i,3) <= Xit) THEN
          Ts(i,3) = 3.d0*beta/4.d0*Xi(i,3)&
               &-sqrt(3.d0)/4.d0*sqrt(beta*Xi(i,3)*(3.d0*Xi(i,3)*beta+8.d0))+1.d0
          T(i,3) = 1.d0
          BL(i,3) = 1/(Ts(i,3)*beta)*(2.d0-2.d0*Ts(i,3)) 
       ELSEIF (Xi(i,3)> Xit) THEN
          Ts(i,3) =(-12.d0*Xi(i,3)+6.d0*H(i,3))/((beta*H(i,3)+6.d0)*H(i,3))
          BL(i,3) = H(i,3)/2.d0
          T(i,3) = Ts(i,3)/4.d0*(beta*H(i,3)+4.d0)
       ENDIF
    END DO

  END SUBROUTINE XI_SPLIT_BALMFORTH

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE TEMPERATURE
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE TEMPERATURE_BALMFORTH(f,col,N,Xi,H,T,Ts,BL,P,dist,ray,Dr,nu,Pe,&
       &delta0,el,grav,N1,tmps,psi,Dt,Rheology,ERROR_CODE)

    !*****************************************************************
    ! Give the vector f
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: f
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,Xi,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dr 
    INTEGER ,INTENT(IN) :: col,N
    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: nu,Pe,delta0,el,grav,N1,tmps,psi,Dt

    ! Parametre pour le sous programme
    DOUBLE PRECISION, PARAMETER :: pi=3.14159265

    DOUBLE PRECISIOn :: omega_a,sigma_a
    DOUBLE PRECISIOn :: omega_b,sigma_b
    DOUBLE PRECISIOn :: Ai,Bi

    DOUBLE PRECISION :: loss,beta
    INTEGER :: i,Na

    ! Remplissage de f

    DO i=1,N,1   

       IF1:IF (i .NE. 1) THEN
          CALL fBi_thermal(ray,dist,Dr,i,Bi)
          CALL fomega_b(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,omega_b)
          CALL fsigma_b(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,sigma_b)
          
       ENDIF IF1

       IF2: IF (i .NE. N) THEN
          CALL fAi_thermal(ray,dist,Dr,i,Ai)
          CALL fomega_a(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,omega_a)
          CALL fsigma_a(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,sigma_a)
       END IF IF2

       ! beta = N1*Pe**(-0.5d0)/(sqrt(pi*(tmps+Dt)))
       ! loss = Pe*beta*psi*Ts(i,col)
       ! print*,'loss Old',loss
       loss = 2*Pe*psi*(T(i,col)-Ts(i,col))/BL(i,col)
       ! print*,i,'loss New',loss

       IF4: IF (i==1) THEN
          f(i)=loss+Ai*Omega_a*Xi(i,col)+Ai*Sigma_a
       ELSEIF (i==N) THEN
          f(i)=loss-Bi*Omega_b*Xi(i-1,col)-Bi*Sigma_b
       ELSE
          f(i)=Ai*Omega_a*Xi(i,col)&
               &-Bi*Omega_b*Xi(i-1,col)&
               &+Ai*Sigma_a-Bi*Sigma_b &
               &+loss
       END IF IF4
       

    END DO

  END SUBROUTINE TEMPERATURE_BALMFORTH

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE JACOBIENNE TEMPERATURE
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE JACOBI_TEMPERATURE_BALMFORTH(a,b,c,N,H,BL,T,Ts,Xi,P,Dr,Dt,dist,ray,nu,Pe,&
       &delta0,el,grav,Rheology,ERROR_CODE)

    !*****************************************************************
    ! Give the jacobian coeficient a1,b1,c1
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: a,b,c
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,BL,T,Ts,Xi,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dr,Dt
    INTEGER ,INTENT(IN) :: N
    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology
    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: nu,Pe,delta0,el,grav

    ! Parametre pour le sous programme

    DOUBLE PRECISIOn :: omega_a,sigma_a
    DOUBLE PRECISIOn :: omega_b,sigma_b
    DOUBLE PRECISIOn :: Ai,Bi

    DOUBLE PRECISION :: loss

    INTEGER :: i,col

    ! Remplissage de la matrice Jacobienne

    col=2

    DO i=1,N,1
       IF1:IF (i .NE. 1) THEN
          CALL fBi_thermal(ray,dist,Dr,i,Bi)
          CALL fomega_b(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,omega_b)
          
       ENDIF IF1
       IF2: IF (i .NE. N) THEN
          CALL fAi_thermal(ray,dist,Dr,i,Ai)
          CALL fomega_a(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,omega_a)
       END IF IF2

       IF3:IF (i==1) THEN
          a(i)=0.d0
          b(i)=Ai*Omega_a
          c(i)=0.d0
       ELSEIF (i==N) THEN
          a(i)=-Bi*Omega_b
          b(i)=0.d0
          c(i)=0.d0
       ELSE
          a(i) = -Bi*Omega_b
          b(i) = Ai*Omega_a
          c(i)=0.d0
       END IF IF3

    ENDDO
  END SUBROUTINE JACOBI_TEMPERATURE_BALMFORTH


END MODULE MODULE_THERMAL_SKIN_NEWTON
