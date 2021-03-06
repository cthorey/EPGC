MODULE MODULE_THICKNESS_SKIN_NEWTON_GRAVI
  USE MODULE_MOBILITY
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINE THICKNESS_NEWTON_SOLVER

  SUBROUTINE  THICKNESS_SKIN_NEWTON_GRAVI(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma&
       &,nu,delta0,z,F_err,theta,Rheology,ERROR_CODE)

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

    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology

    !Variable du sous programmes
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coeff
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fguess,ftmps
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,d,e,f,g,k,l,S
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a1,b1,c1,d1,e1,f1,g1,k1,l1
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Hm,qa

    DOUBLE PRECISION :: U,Test_Eta,Test_Eta2
    INTEGER :: ndyke,i,N,Size,code
    INTEGER :: err1,col,algo1
    LOGICAL :: cho

    ! Taille de la grille sur laquelle on fait l'inversion
    ndyke=sigma/Dr
    CHO=COUNT(H(:,2)>0.D0)<ndyke
    SELECT CASE (CHO)
    CASE(.TRUE.)
       N = ndyke+1
    CASE(.FALSE.)
       N = COUNT(H(:,2)>1D-10)
       Test_Eta = ABS((H(N,2)-H(N-1,2))/Dr)
       Test_Eta2 =H(N,2)/Dr
       ! print*,Test_Eta,Test_Eta2
       IF (Test_Eta2 >Test_Eta) THEN
          N = N+1
       ENDIF
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


    ! Calcule coefficient pression elastique

    ALLOCATE(Coeff(1:N,7),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur alloc dans coeff du systeme'; STOP
    END IF


    ! Calcule de f tmps n et n+1
    ALLOCATE(ftmps(1:N),fguess(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur de f';STOP
    END IF

    col=1
    CALL THICKNESS(ftmps,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,dist,ray,qa,&
         &el,grav,nu,delta0,Rheology,ERROR_CODE)
    col=2
    CALL THICKNESS(fguess,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,dist,ray,qa,&
         &el,grav,nu,delta0,Rheology,ERROR_CODE)

    ! Jacbienne
    ALLOCATE(a1(1:N),b1(1:N),c1(1:N),d1(1:N),e1(1:N), stat = err1)
    ALLOCATE(f1(1:N),g1(1:N),k1(1:N),l1(1:N), stat = err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF
    
    CALL JACOBI_THICKNESS(a1,b1,c1,d1,e1,f1,g1,k1,l1,N,M,H,P,T,BL,Ts,&
         Coeff,Dt,Dr,dist,ray,el,grav,nu,delta0,Rheology,ERROR_CODE)

    !Systeme a inverser
    ALLOCATE(a(1:N),b(1:N),c(1:N),d(1:N),e(1:N),stat = err1)
    ALLOCATE(f(1:N),g(1:N),k(1:N),l(1:N),S(1:N),stat = err1)
    IF (err1>1) THEN
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
       S(i)=H(i,1)-H(i,2)+theta*Dt*fguess(i)+(1-theta)*Dt*ftmps(i)
    END DO

    d(1)=0.d0
    f(N) =0.d0

    !Inversion de la matrice
    ALLOCATE(Hm(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF

    CALL TRIDIAG(d,e,f,S,N,Hm)

    DO i=1,N,1
       H(i,3)=Hm(i)+H(i,2)
    END DO


    ! Calcul du soeuil F_err
    IF (DOT_PRODUCT(H(:,2),H(:,2)) == 0D0) THEN
       F_err = ABS(MAXVAL((Hm(:))))
    ELSE
       Size = COUNT(H(:,2)>0D0)
       F_err = ABS(MAXVAL((H(:Size,3)-H(:Size,2))/H(:Size,2)))
    ENDIF

    DEALLOCATE(Hm,a,b,c,d,e,f,g,k,l,S,Coeff,qa)
    DEALLOCATE(fguess,ftmps,a1,b1,c1,d1,e1,f1,g1,k1,l1)



  END SUBROUTINE THICKNESS_SKIN_NEWTON_GRAVI

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE THICKNESS
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  SUBROUTINE   THICKNESS(f,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,dist,ray,&
       &qa,el,grav,nu,delta0,Rheology,ERROR_CODE)

    !*****************************************************************
    ! Give the vector f 
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: f
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,BL,Ts
    DOUBLE PRECISION  ,DIMENSION(:,:), INTENT(INOUT) :: P,Coeff
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: qa
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr 
    INTEGER ,INTENT(IN) :: col,N,M

    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,nu,delta0

    ! Parametre pour le sous programme
    DOUBLE PRECISION :: phi_a,phi_b
    DOUBLE PRECISION :: Ai,Bi
    DOUBLE PRECISION :: Ael,Bel,Agrav,Bgrav
    INTEGER :: i,err1,algo1

    ! Remplissage de f

    !### REMPLISSAGE DE f ###!

    DO i=1,N,1
       IF (i .NE. N) THEN
          CALL fAi_thickness(ray,dist,Dr,i,Ai)
          CALL fPhi_A(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,phi_a)
       ENDIF

       IF (i .NE. 1) THEN
          CALL fBi_thickness(ray,dist,Dr,i,Bi)
          CALL fPhi_B(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
              &,nu,Rheology,ERROR_CODE,phi_b)
       ENDIF

       Ael = el*Ai
       Agrav = grav*Ai
       Bel = el*Bi
       Bgrav = grav*Bi

       IF (i==1) THEN
          f(i)=Ael*phi_a*(P(2,col)-P(1,col))+Agrav*phi_a*(H(2,col)-H(1,col))&
               &+qa(i)
       ELSEIF (i==N) THEN
          f(i)=-Bel*phi_b*(P(i,col)-P(i-1,col))-Bgrav*phi_b*(H(i,col)-H(i-1,col))&
               &+qa(i)
       ELSE
          f(i)=Ael*phi_a*(P(i+1,col)-P(i,col))-Bel*phi_b*(P(i,col)-P(i-1,col))&
               &+Agrav*phi_a*(H(i+1,col)-H(i,col))-Bgrav*phi_b*(H(i,col)-H(i-1,col))&
               &+qa(i)
       END IF

    END DO

  END SUBROUTINE THICKNESS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE THICKNESS
  !-------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE JACOBI_THICKNESS(a,b,c,d,e,f,g,k,l,N,M,H,P,T,BL,Ts,Coeff,&
       &Dt,Dr,dist,ray,el,grav,nu,delta0,Rheology,ERROR_CODE)

    !*****************************************************************
    ! Give the jacobian coeficient a1,b1,c1
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: a,b,c,d,e,f,g,k,l
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,BL,Ts
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: P,Coeff
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    ! Prametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr 
    INTEGER ,INTENT(IN) :: N,M

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,nu,delta0
   !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology
    
    ! Parametre pour le sous programme
    DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: alpha,beta,gamma,lambda,kappa,delta,epsilonn
    DOUBLE PRECISION ::Ael,Bel,Agrav,Bgrav,Ai,Bi
    DOUBLE PRECISION :: phi_a,phi_b,dphib_dhi,dphib_dhi1,dphia_dhi,dphia_dhi1
    DOUBLE PRECISION :: H1,H2,P1,P2
    INTEGER :: i,col,algo1,err1

    ! Allocation + remplissage pression
    ALLOCATE(alpha(1:N),beta(1:N),gamma(1:N),lambda(1:N),stat=err1)
    ALLOCATE(kappa(1:N),delta(1:N),epsilonn(1:N),stat=err1)
    IF (err1>1) THEN !On teste si nos tableau sont bien alloués
       PRINT*,"Erreur ds alloc alpha,beta.." ;STOP
    END IF

    algo1=2;col=2

    alpha=Coeff(:,1)
    beta=Coeff(:,2)
    gamma=Coeff(:,3)
    lambda=Coeff(:,4)
    kappa=Coeff(:,5)
    delta=Coeff(:,6)
    epsilonn=Coeff(:,7)

    ! Remplissage de la matrice Jacobienne
    DO i=1,N,1

       IF1: IF (i .NE. N) THEN
          CALL fAi_thickness(ray,dist,Dr,i,Ai)
          CALL fPhi_A(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,phi_a)
          CALL fdPhi_A(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,dphia_dhi,dphia_dhi1)

          P1=P(i+1,2)-P(i,2); P2=P(i,2)-P(i-1,2)
          H1=H(i+1,2)-H(i,2); H2=H(i,2)-H(i-1,2)
       ENDIF IF1

       IF2: IF (i .NE. 1) THEN

          CALL fBi_thickness(ray,dist,Dr,i,Bi)
          CALL fPhi_B(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,phi_b)
          CALL fDPhi_B(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
               &,nu,Rheology,ERROR_CODE,dphib_dhi,dphib_dhi1)
                    
          P2=P(i,2)-P(i-1,2)
          H2=H(i,2)-H(i-1,2)
       ENDIF IF2

       Ael = el*Ai
       Agrav = grav*Ai
       Bel = el*Bi
       Bgrav = grav*Bi
    
       IF3: IF (i==1) THEN
          a(i)=0;b(i)=0;c(i)=0;d(i)=0
          e(i)=Agrav*(dphia_dhi*H1-phi_a)
          f(i)=Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=0.d0
          k(i)=0.d0
          l(i)=0.d0
          ELSEIF (i ==N) THEN
             a(i)=0.d0
             b(i)=0.d0
             c(i)=0.d0
             d(i)=Bgrav*(-dphib_dhi1*H2+phi_b)
             e(i)=-Bgrav*(dphib_dhi*H2+phi_b)
             f(i)=0.d0
             g(i)=0.D0
             k(i)=0.d0
             l(i)=0.d0
       ELSE
          a(i)=0.d0
          b(i)=0.d0
          c(i)=0.d0
          d(i)=Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=0.D0
          k(i)=0.d0
          l(i)=0.d0
       END IF IF3

    END DO

    DEALLOCATE(alpha,beta,gamma,lambda,kappa,delta,epsilonn)
  END SUBROUTINE JACOBI_THICKNESS

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
  
END MODULE MODULE_THICKNESS_SKIN_NEWTON_GRAVI

