MODULE MODULE_THERMAL_GFD_GRAVI_ARRHENIUS

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Subroutine pour résoudre eqaution de la chaleur à l'aide d'une 
!!!!!!! schema difference fini

  SUBROUTINE THERMAL_GFD_SOLVER_ARRHENIUS(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)

    !*****************************************************************
    ! Solve for the parameter Xi, and split in Temperature and thermal layer
    ! from  the evolution euqation using a center difference in space and a general 
    ! theta scheme is timre
    !*****************************************************************
    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: H,P
    DOUBLE PRECISION , DIMENSION(:,:), INTENT(INOUT) :: Xi,T,BL,Ts
    DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: dist,ray

    !Parametre du model
    DOUBLE PRECISION , INTENT(IN) :: Dt,Dr,theta,tmps

    !Nombre sans dimensions
    DOUBLE PRECISION , INTENT(IN) :: sigma,nu,Pe,psi,delta0,el,grav,N1
    INTEGER, INTENT(IN) :: M, z
    DOUBLE PRECISION , INTENT(INOUT) :: F_err

    !Variable du sous programmes
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: Xi_tmps
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,d,e,f,g,k,l,S
    DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: Xi_m

    DOUBLE PRECISION :: U
    INTEGER :: i,ndyke,N,Size
    INTEGER :: err1,col
    LOGICAL :: CHO


    ! Taille de la grille
    N = COUNT(H(:,3)>0.D0)

    !Systeme a inverser
    ALLOCATE(a(1:N),b(1:N),c(1:N),S(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur allocation dans coeff Temperature'; STOP
    END IF
     CALL FILLING_MATRIX(a,b,c,S,N,H,BL,T,Ts,Xi,P,Dr,dist,ray,nu,&
          &N1,Pe,delta0,el,grav,tmps,Dt,theta,psi)

    a(1)=0
    c(N)=0

    !Inversion de la matrice
    ALLOCATE(Xi_m(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF
    CALL TRIDIAG(a,b,c,s,N,Xi_m)
    DO i=1,N,1
       Xi(i,3)=Xi_m(i)
       IF (Xi(i,3) >H(i,3)/2D0) THEN
          Xi(i:,3) = H(i:,3)/2D0
          EXIT
       ELSEIF (XI(i,3)<0.D0) THEN
          Xi(i,3) =0.D0
       ENDIF
    END DO

    ! Separation variables
    CALL XI_SPLIT(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)

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

  END SUBROUTINE THERMAL_GFD_SOLVER_ARRHENIUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Routine pour separer les variables
!!!!!!! 

SUBROUTINE XI_SPLIT(Xi,T,BL,Ts,H,N,delta0,Dt,tmps,N1,Pe,el)

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


  END SUBROUTINE XI_SPLIT

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
  
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE FILLING MATRIX
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE FILLING_MATRIX(a,b,c,S,N,H,BL,T,Ts,Xi,P,Dr,dist,ray,nu,N1,Pe,delta0,el,grav,tmps,Dt,theta,psi)

    !*****************************************************************
    ! Give the jacobian coeficient a1,b1,c1
    !*****************************************************************
  IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: a,b,c,S
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,Xi,T,Ts,P,BL
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dr 
    INTEGER ,INTENT(IN) :: N

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: nu,Pe,delta0,el,grav,N1,tmps,Dt,theta,psi

    ! Parametre pour le sous programme
    DOUBLE PRECISION, PARAMETER :: pi=3.14159265

    DOUBLE PRECISION :: h_a,delta_a,delta_a2,eta_a,Ai,T_a,delta_a3
    DOUBLE PRECISIOn :: omega_a,sigma_a

    DOUBLE PRECISION :: h_b,delta_b,delta_b2,eta_b,Bi,T_b,delta_b3
    DOUBLE PRECISIOn :: omega_b,sigma_b,Ts_a,Ts_b,Ds_b,Ds_a
    DOUBLE PRECISION :: loss,beta,Crys
    INTEGER :: i,col

    ! Remplissage de f

    col = 2

    DO i=1,N,1  
       IF1:IF (i .NE. 1) THEN
          eta_b=(grav*(H(i,3)-H(i-1,3)))/Dr
          Bi=(ray(i-1)/(dist(i)*Dr))*Dt
          h_b=0.5d0*(H(i,3)+H(i-1,3))
          delta_b=0.5d0*(BL(i,col)+BL(i-1,col))
          delta_b2=0.5d0*(BL(i,col)**2+BL(i-1,col)**2)
          delta_b3=0.5d0*(BL(i,col)**3+BL(i-1,col)**3)
          T_b = 0.5d0*(T(i,col)+T(i-1,col))
          Ts_b = 0.5d0*(Ts(i,col)+Ts(i-1,col))
          Ds_b = T_b-Ts_b

          IF (nu < 1D0) THEN
              omega_b = 3*sqrt(pi)*Ds_b**(-1.5d0)*delta_b2*eta_b*nu*nu**(-T_b)&
             & *(-log(nu))**(-1.5d0)*erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*1.0/Ds_b&
             & *delta_b2*eta_b*nu*nu**(-T_b)*1.0/(-log(nu)) - 3*1.0/Ds_b*&
             & delta_b*eta_b*h_b*nu*nu**Ds_b*nu**(-T_b)*1.0/(-log(nu)) + 3*1.0/&
             & Ds_b*delta_b*eta_b*h_b*nu*nu**(-T_b)*1.0/(-log(nu))

             sigma_b = (3.0d0/2.0d0)*sqrt(pi)*Ds_b**(-1.5d0)*delta_b3*eta_b*&
             & nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*erf(sqrt(Ds_b)*sqrt(-log(nu&
             & ))) - 1.0/Ds_b*delta_b3*eta_b*nu*nu**Ds_b*nu**(-T_b)*(-log(nu))&
             & **(-2.0d0) - 2*1.0/Ds_b*delta_b3*eta_b*nu*nu**(-T_b)*(-log(nu))&
             & **(-2.0d0) - 1.0/Ds_b*delta_b2*eta_b*h_b*nu*nu**Ds_b*nu**(-T_b)&
             & *(-log(nu))**(-2.0d0) + 1.0/Ds_b*delta_b2*eta_b*h_b*nu*nu**(&
             & -T_b)*(-log(nu))**(-2.0d0) + sqrt(pi)*Ds_b**(-0.5d0)*delta_b3*&
             & eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*log(nu)*erf(sqrt(Ds_b)*&
             & sqrt(-log(nu))) - 2*delta_b3*eta_b*nu*nu**(-T_b)*(-log(nu))**(&
             & -2.0d0)*log(nu) + delta_b2*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**&
             & (-2.0d0)*log(nu)
          ELSEIF (nu .EQ. 1D0) THEN
             omega_b = -2*delta_b**2*eta_b + 3*delta_b*eta_b*h_b
             sigma_b = (7.0d0/15.0d0)*Ds_b*delta_b**3*eta_b - 1.0d0/2.0d0*Ds_b*&
                  & delta_b**2*eta_b*h_b
          ENDIF
       ENDIF IF1

       IF2: IF (i .NE. N) THEN
          eta_a=(grav*(H(i+1,3)-H(i,3)))/Dr
          Ai=(ray(i)/(dist(i)*Dr))*Dt
          h_a=0.5d0*(H(i+1,3)+H(i,3))
          delta_a=0.5d0*(BL(i+1,col)+BL(i,col))
          delta_a2 = 0.5d0*(BL(i+1,col)**2+BL(i,col)**2)
          delta_a3 = 0.5d0*(BL(i+1,col)**3+BL(i,col)**3)
          T_a = 0.5d0*(T(i,col)+T(i+1,col))
          Ts_a = 0.5d0*(Ts(i,col)+Ts(i+1,col))
          Ds_a = T_a-Ts_a
          
          IF (nu < 1D0) THEN
             omega_a = 3*sqrt(pi)*Ds_a**(-1.5d0)*delta_a2*eta_a*nu*nu**(-T_a)&
             & *(-log(nu))**(-1.5d0)*erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*1.0/Ds_a&
             & *delta_a2*eta_a*nu*nu**(-T_a)*1.0/(-log(nu)) - 3*1.0/Ds_a*&
             & delta_a*eta_a*h_a*nu*nu**Ds_a*nu**(-T_a)*1.0/(-log(nu)) + 3*1.0/&
             & Ds_a*delta_a*eta_a*h_a*nu*nu**(-T_a)*1.0/(-log(nu))

             sigma_a = (3.0d0/2.0d0)*sqrt(pi)*Ds_a**(-1.5d0)*delta_a3*eta_a*&
             & nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*erf(sqrt(Ds_a)*sqrt(-log(nu&
             & ))) - 1.0/Ds_a*delta_a3*eta_a*nu*nu**Ds_a*nu**(-T_a)*(-log(nu))&
             & **(-2.0d0) - 2*1.0/Ds_a*delta_a3*eta_a*nu*nu**(-T_a)*(-log(nu))&
             & **(-2.0d0) - 1.0/Ds_a*delta_a2*eta_a*h_a*nu*nu**Ds_a*nu**(-T_a)&
             & *(-log(nu))**(-2.0d0) + 1.0/Ds_a*delta_a2*eta_a*h_a*nu*nu**(&
             & -T_a)*(-log(nu))**(-2.0d0) + sqrt(pi)*Ds_a**(-0.5d0)*delta_a3*&
             & eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*log(nu)*erf(sqrt(Ds_a)*&
             & sqrt(-log(nu))) - 2*delta_a3*eta_a*nu*nu**(-T_a)*(-log(nu))**(&
             & -2.0d0)*log(nu) + delta_a2*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**&
             & (-2.0d0)*log(nu)
          ELSEIF (nu .EQ. 1D0) THEN
             omega_a = -2*delta_a2*eta_a + 3*delta_a*eta_a*h_a
             sigma_a = (7.0d0/15.0d0)*Ds_a*delta_a3*eta_a - 1.0d0/2.0d0*Ds_a*&
                  & delta_a2*eta_a*h_a
          ENDIF

       END IF IF2
       
       beta = N1*Pe**(-0.5d0)/(sqrt(pi*(tmps+Dt)))
       loss = Pe*beta*Ts(i,col)*Dt*psi
       IF (ABS(N-i)<2) THEN
          loss = loss/2.0
       ENDIF

       IF3:IF (i==1) THEN
          a(i) = 0.d0
          b(i) = 1D0-Ai*Omega_a
          c(i) = 0.d0
          S(i) = Xi(i,1)+Ai*Sigma_a+loss
       ELSEIF (i == N) THEN
          a(i) = Bi*Omega_b
          b(i) = 1D0
          c(i) = 0.d0
          S(i) = Xi(i,1)-Bi*Sigma_b+loss
       ELSE
          a(i) = Bi*Omega_b
          b(i) = 1D0-Ai*Omega_a
          c(i) = 0.d0
          S(i) = Xi(i,1)+Ai*Sigma_a-Bi*Sigma_b+loss
       END IF IF3
    ENDDO
  END SUBROUTINE FILLING_MATRIX
 

END MODULE MODULE_THERMAL_GFD_GRAVI_ARRHENIUS
