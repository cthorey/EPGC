MODULE MODULE_THICKNESS_SKIN_NEWTON

  USE MOBILITY_THICKNESS_SKIN_RHEOLOGY
  
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINE THICKNESS_NEWTON_SOLVER

  SUBROUTINE  THICKNESS_SKIN_NEWTON(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,&
       &gam,Inter_Q,z,F_err,theta,tmps,Rheology,ERROR_CODE)

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
    DOUBLE PRECISION, INTENT(IN) :: tmps

    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology
    
    !Nombre sans dimensions
    DOUBLE PRECISION , INTENT(IN) :: el,grav,sigma,nu,delta0,gam,Inter_Q
    INTEGER, INTENT(IN) :: M, z
    DOUBLE PRECISION , INTENT(INOUT) :: F_err

    !Variable du sous programmes
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coeff
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fguess,ftmps
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a,b,c,d,e,f,g,k,l,S
    DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: a1,b1,c1,d1,e1,f1,g1,k1,l1
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

    InterInjectionRate:IF (mod(tmps,Inter_Q)<Inter_Q/2D0) THEN
       DO i = 1,N,1
          U = 2.d0/(sigma)**4.
          Flux:IF (i<ndyke+1) THEN
             qa(i) = U*(1-gam*H(i,2))*(sigma**2.-dist(i)**2.)
          ELSE 
             qa(i) = 0.d0
          END IF Flux
       END DO
    ELSE
       qa(:)=0D0
    ENDIF InterInjectionRate
    
    ! Calcule coefficient pression elastique

    ALLOCATE(Coeff(1:N,7),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur alloc dans coeff du systeme'; STOP
    END IF

    algo1=1;col=1
    CALL PRESSURE_ORDER_4(M,N,Dr,P,H,dist,Coeff,algo1,col,el,delta0)

    ! Calcule de f tmps n et n+1
    ALLOCATE(ftmps(1:N),fguess(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur de f';STOP
    END IF

    col=1
    CALL THICKNESS(ftmps,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,&
         &dist,ray,qa,el,grav,nu,delta0,Rheology,ERROR_CODE)
    col=2
    CALL THICKNESS(fguess,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,&
         &dist,ray,qa,el,grav,nu,delta0,Rheology,ERROR_CODE)

    ! Jacbienne
    ALLOCATE(a1(1:N),b1(1:N),c1(1:N),d1(1:N),e1(1:N), stat = err1)
    ALLOCATE(f1(1:N),g1(1:N),k1(1:N),l1(1:N), stat = err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans coeff du systeme'; STOP
    END IF
    
    CALL JACOBI_THICKNESS(a1,b1,c1,d1,e1,f1,g1,k1,l1,N,M,H,P,T,BL,&
         &Ts,Coeff,Dt,Dr,dist,ray,el,grav,nu,delta0,Rheology,ERROR_CODE)

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


    a(1)=0.d0;a(2)=0.d0;a(3)=0.d0;a(4)=0.d0
    b(1)=0.d0;b(2)=0.d0;b(3)=0.d0
    c(1)=0.d0;c(2)=0.d0
    d(1)=0.d0
    l(N-3)=0.d0;l(N-2)=0.d0;l(N-1)=0.d0;l(N)=0.d0
    k(N-2)=0.d0;k(N-1)=0.d0;k(N)=0.d0
    g(N-1)=0.d0;g(N)=0.d0
    f(N)=0.d0

    !Inversion de la matrice
    ALLOCATE(Hm(1:N),stat=err1)
    IF (err1>1) THEN
       PRINT*, 'Erreur d''allocation dans vecteur Hm'; STOP
    END IF

    CALL NONA_DIAGO(N,Hm,a,b,c,d,e,f,g,k,l,S)

    DO i=1,N,1
       H(i,3)=Hm(i)+H(i,2)
    END DO

    H(N+1:M,3) = delta0

    ! Calcul de la pression
    algo1=2;col=3
    CALL PRESSURE_ORDER_4(M,N,Dr,P,H,dist,Coeff,algo1,col,el,delta0)

    ! Calcul du soeuil F_err
    IF (DOT_PRODUCT(H(:,2),H(:,2)) == 0D0) THEN
       F_err = ABS(MAXVAL((Hm(:))))
    ELSE
       Size = COUNT(H(:,2)>1D-6)
       F_err = ABS(MAXVAL((H(:Size,3)-H(:Size,2))/H(:Size,2)))
    ENDIF

    DEALLOCATE(Hm,a,b,c,d,e,f,g,k,l,S,Coeff,qa)
    DEALLOCATE(fguess,ftmps,a1,b1,c1,d1,e1,f1,g1,k1,l1)

  END SUBROUTINE THICKNESS_SKIN_NEWTON

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE PRESSURE
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------

  SUBROUTINE  PRESSURE_ORDER_4(M,N,Dr,P,H,dist,Coeff,algo1,col,el,delta0)

    !*****************************************************************
    ! Calcul the coefficient of a 7 stencil elastic
    !*****************************************************************

    !-------------------------------------------------------------------------------------
    !  Variables

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: H
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: Coeff,P
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist

    ! Parametres model
    INTEGER , INTENT(IN) :: M,N,col,algo1
    DOUBLE PRECISION, INTENT(IN) :: Dr,el,delta0

    ! Variables subroutines
    INTEGER :: i,err1
    DOUBLE PRECISION :: p1,p2,p3,p4

    !-------------------------------------------------------------------------------------
    !  Contains

    SELECT CASE(algo1)
    CASE(1)

       DO i=1,N,1
          p1=1.d0/dist(i)**3
          p2=-1.d0/(dist(i)**2)
          p3=2.d0/dist(i)
          p4=1.d0

          Coeff(i,1)=(1.d0/(24.d0*Dr**4))*(-4.d0*p4+3.d0*p3*Dr)
          Coeff(i,2)=(1.d0/(24.d0*Dr**4))*(48.d0*p4-24.d0*p3*Dr-2.d0*p2*Dr**2+2.d0*p1*Dr**3)
          Coeff(i,3)=(1.d0/(24.d0*Dr**4))*(-156.d0*p4+39.d0*p3*Dr+32.d0*p2*Dr**2-16.d0*p1*Dr**3)
          Coeff(i,4)=(1.d0/(24.d0*Dr**4))*(224.d0*p4-60.d0*p2*Dr**2)
          Coeff(i,5)=(1.d0/(24.d0*Dr**4))*(-156.d0*p4-39*p3*Dr+32.d0*p2*Dr**2+16.d0*p1*Dr**3)
          Coeff(i,6)=(1.d0/(24.d0*Dr**4))*(48.d0*p4+24.d0*p3*Dr-2.d0*p2*Dr**2-2.d0*p1*Dr**3)
          Coeff(i,7)=(1.d0/(24.d0*Dr**4))*(-4.d0*p4-3.d0*p3*Dr)

       END DO

    CASE(2)
       DO i=1,M,1
          IF (i==1) THEN
             P(i,col)=el*(Coeff(i,1)*H(3,col)+Coeff(i,2)*H(2,col)+Coeff(i,3)*H(1,col)+Coeff(i,4)*H(1,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*H(i+2,col)+Coeff(i,7)*H(i+3,col))
          ELSEIF (i==2) THEN
             P(i,col)=el*(Coeff(i,1)*H(2,col)+Coeff(i,2)*H(1,col)+Coeff(i,3)*H(1,col)+Coeff(i,4)*H(2,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*H(i+2,col)+Coeff(i,7)*H(i+3,col))
          ELSEIF (i==3) THEN
             P(i,col)=el*(Coeff(i,1)*H(1,col)+Coeff(i,2)*H(1,col)+Coeff(i,3)*H(2,col)+Coeff(i,4)*H(3,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*H(i+2,col)+Coeff(i,7)*H(i+3,col))
          ELSEIF (i==N-2) THEN
             P(i,col)=Coeff(i,1)*H(i-3,col)+Coeff(i,2)*H(i-2,col)+Coeff(i,3)*H(i-1,col)+Coeff(i,4)*H(i,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*H(i+2,col)+Coeff(i,7)*delta0
          ELSEIF (i==N-1) THEN
             P(i,col)=Coeff(i,1)*H(i-3,col)+Coeff(i,2)*H(i-2,col)+Coeff(i,3)*H(i-1,col)+Coeff(i,4)*H(i,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*delta0+Coeff(i,7)*delta0
          ELSEIF (i==N) THEN
             P(i,col)=Coeff(i,1)*H(i-3,col)+Coeff(i,2)*H(i-2,col)+Coeff(i,3)*H(i-1,col)+Coeff(i,4)*H(i,col)+Coeff(i,5)*delta0&
                  &+Coeff(i,6)*delta0+Coeff(i,7)*delta0
          ELSE
             P(i,col)=Coeff(i,1)*H(i-3,col)+Coeff(i,2)*H(i-2,col)+Coeff(i,3)*H(i-1,col)+Coeff(i,4)*H(i,col)+Coeff(i,5)*H(i+1,col)&
                  &+Coeff(i,6)*H(i+2,col)+Coeff(i,7)*H(i+3,col)
          END IF
          IF (el==0) THEN
             P(i,col)=0
          END IF
       END DO

    END SELECT

  END SUBROUTINE PRESSURE_ORDER_4

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  !  SUBROUTINE THICKNESS
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  SUBROUTINE   THICKNESS(f,col,N,M,H,P,T,BL,Ts,Coeff,Dt,Dr,&
       &dist,ray,qa,el,grav,nu,delta0,Rheology,ERROR_CODE)

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

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,nu,delta0

    !Parametre a transletre
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Rheology
    
    ! Parametre pour le sous programme
    DOUBLE PRECISION :: h_a,h_b,h_a2,h_b2,h_a3,h_b3,T_a,T_b
    DOUBLE PRECISION :: delta_a,delta_b,delta_a2,delta_b2,delta_a3,delta_b3
    DOUBLE PRECISION :: phi_a,phi_b
    DOUBLE PRECISION :: Ael,Bel,Agrav,Bgrav
    DOUBLE PRECISION :: Ts_a,Ts_b,Delta_T_b,Delta_T_a
    INTEGER :: i,err1,algo1

    ! Remplissage de f

    algo1=2
    CALL  PRESSURE_ORDER_4(M,N,Dr,P,H,dist,Coeff,algo1,col,el,delta0)

    !### REMPLISSAGE DE f ###!

    DO i=1,N,1
       IF (i .NE. N) THEN
          Ael=el*(ray(i)/(dist(i)*Dr**2))
          Agrav=grav*(ray(i)/(dist(i)*Dr**2))
          h_a=0.5d0*(H(i+1,col)+H(i,col))
          h_a2=0.5d0*(H(i+1,col)**2+H(i,col)**2)
          h_a3=0.5d0*(H(i+1,col)**3+H(i,col)**3)
          delta_a=0.5d0*(BL(i+1,3)+BL(i,3))
          delta_a2=0.5d0*(BL(i+1,3)**2+BL(i,3)**2)
          delta_a3=0.5d0*(BL(i+1,3)**3+BL(i,3)**3)
          T_a=0.5d0*(T(i,3)+T(i+1,3))
          Ts_a = 0.5d0*(Ts(i,3)+Ts(i+1,3))
          Delta_T_a = T_a -Ts_a

          CALL fPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,phi_a,nu,Rheology,ERROR_CODE)
       ENDIF

       IF (i .NE. 1) THEN
          Bel=el*(ray(i-1)/(dist(i)*Dr**2))
          Bgrav=grav*(ray(i-1)/(dist(i)*Dr**2))
          h_b=0.5d0*(H(i,col)+H(i-1,col))
          h_b2=0.5d0*(H(i,col)**2+H(i-1,col)**2)
          h_b3=0.5d0*(H(i,col)**3+H(i-1,col)**3)
          delta_b=0.5d0*(BL(i,3)+BL(i-1,3))
          delta_b2=0.5d0*(BL(i,3)**2+BL(i-1,3)**2)
          delta_b3=0.5d0*(BL(i,3)**3+BL(i-1,3)**3)
          T_b=0.5d0*(T(i,3)+T(i-1,3))
          Ts_b = 0.5d0*(Ts(i,3)+Ts(i-1,3))
          Delta_T_b = T_b - Ts_b

          CALL fPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,phi_b,nu,Rheology,ERROR_CODE)

       ENDIF

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

  SUBROUTINE JACOBI_THICKNESS(a,b,c,d,e,f,g,k,l,N,M,H,P,T,BL,Ts,&
       &Coeff,Dt,Dr,dist,ray,el,grav,nu,delta0,Rheology,ERROR_CODE)

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
    DOUBLE PRECISION ::Ael,Bel,Agrav,Bgrav,h_a,h_b,h_a2,h_b2,h_a3,h_b3,T_a,T_b
    DOUBLE PRECISION ::delta_a,delta_b,delta_a2,delta_b2,delta_a3,delta_b3
    DOUBLE PRECISION :: Ts_a,Ts_b,Delta_T_b,Delta_T_a
    DOUBLE PRECISION :: phi_a,phi_b,dphib_dhi,dphib_dhi1,dphia_dhi,dphia_dhi1
    DOUBLE PRECISION :: H1,H2,P1,P2,hi,hi2,hia,hib,hia2,hib2
    INTEGER :: i,col,algo1,err1

    ! Allocation + remplissage pression
    ALLOCATE(alpha(1:N),beta(1:N),gamma(1:N),lambda(1:N),stat=err1)
    ALLOCATE(kappa(1:N),delta(1:N),epsilonn(1:N),stat=err1)
    IF (err1>1) THEN !On teste si nos tableau sont bien alloués
       PRINT*,"Erreur ds alloc alpha,beta.." ;STOP
    END IF

    algo1=2;col=2
    CALL  PRESSURE_ORDER_4(M,N,Dr,P,H,dist,Coeff,algo1,col,el,delta0)

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
          Ael=el*(ray(i)/(dist(i)*Dr**2))
          Agrav=grav*(ray(i)/(dist(i)*Dr**2))
          h_a=0.5d0*(H(i+1,col)+H(i,col))
          h_a2=0.5d0*(H(i+1,col)**2+H(i,col)**2)
          h_a3=0.5d0*(H(i+1,col)**3+H(i,col)**3)
          delta_a=0.5d0*(BL(i+1,3)+BL(i,3))
          delta_a2=0.5d0*(BL(i+1,3)**2+BL(i,3)**2)
          delta_a3=0.5d0*(BL(i+1,3)**3+BL(i,3)**3)
          T_a=0.5d0*(T(i,3)+T(i+1,3))
          Ts_a = 0.5d0*(Ts(i,3)+Ts(i+1,3))
          Delta_T_a = T_a -Ts_a
          hia = H(i+1,col)
          hi = H(i,col)

          CALL fPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,phi_a,nu,Rheology,ERROR_CODE)
          CALL fdPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,hi,hia,&
               &dphia_dhi1,dphia_dhi,nu,Rheology,ERROR_CODE)
          
          P1=P(i+1,2)-P(i,2); P2=P(i,2)-P(i-1,2)
          H1=H(i+1,2)-H(i,2); H2=H(i,2)-H(i-1,2)
       ENDIF IF1

       IF2: IF (i .NE. 1) THEN
          Bel=el*(ray(i-1)/(dist(i)*Dr**2))
          Bgrav=grav*(ray(i-1)/(dist(i)*Dr**2))
          h_b=0.5d0*(H(i,col)+H(i-1,col))
          h_b2=0.5d0*(H(i,col)**2+H(i-1,col)**2)
          h_b3=0.5d0*(H(i,col)**3+H(i-1,col)**3)
          delta_b=0.5d0*(BL(i,3)+BL(i-1,3))
          delta_b2=0.5d0*(BL(i,3)**2+BL(i-1,3)**2)
          delta_b3=0.5d0*(BL(i,3)**3+BL(i-1,3)**3)
          T_b=0.5d0*(T(i,3)+T(i-1,3))
          Ts_b = 0.5d0*(Ts(i,3)+Ts(i-1,3))
          Delta_T_b = T_b - Ts_b
          hib = H(i-1,col); hib2 = H(i-1,col)**2
          hi = H(i,col);hi2 = H(i,col)**2

          CALL fPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,phi_b,nu,Rheology,ERROR_CODE)
          CALL fDPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,hib,hi,&
               &dphib_dhi,dphib_dhi1,nu,Rheology,ERROR_CODE)
                    
          P2=P(i,2)-P(i-1,2)
          H2=H(i,2)-H(i-1,2)
       ENDIF IF2

       IF3: IF (i==1) THEN
          a(i)=0;b(i)=0;c(i)=0;d(i)=0
          e(i)=Ael*dphia_dhi*P1+Ael*phi_a*((gamma(i+1)+beta(i+1))-(lambda(i)+gamma(i)))&
               &+Agrav*(dphia_dhi*H1-phi_a)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*((lambda(i+1)+alpha(i+1))-(kappa(i)+beta(i)))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-(delta(i)+alpha(i)))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))
          l(i)=Ael*phi_a*epsilonn(i+1)

       ELSEIF (i==2) THEN
          a(i)=0;b(i)=0;c(i)=0
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*((beta(i+1)+alpha(i+1))-(gamma(i)+beta(i)))&
               &-Bel*phi_b*((gamma(i)+beta(i))-(lambda(i-1)+gamma(i-1)))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-(lambda(i)+alpha(i)))&
               &-Bel*phi_b*((lambda(i)+alpha(i))-(kappa(i-1)+beta(i-1)))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))&
               &-Bel*phi_b*(kappa(i)-(delta(i-1)+alpha(i-1)))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))-Bel*phi_b*epsilonn(i)
          l(i)=Ael*phi_a*epsilonn(i+1)
       ELSEIF (i==3) THEN
          a(i)=0;b(i)=0
          c(i)=Ael*phi_a*(alpha(i+1)-(beta(i)+alpha(i)))-Bel*phi_b*((beta(i)+alpha(i))-(gamma(i-1)+beta(i-1)))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-(lambda(i-1)+alpha(i-1)))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))-Bel*phi_b*epsilonn(i)
          l(i)=Ael*phi_a*epsilonn(i+1)
       ELSEIF (i==4) THEN
          a(i)=0
          b(i)=-Ael*phi_a*alpha(i)-Bel*phi_b*(alpha(i)-(beta(i-1)+alpha(i-1)))
          c(i)=Ael*phi_a*(alpha(i+1)-beta(i))-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))-Bel*phi_b*epsilonn(i)
          l(i)=Ael*phi_a*epsilonn(i+1)

       ELSEIF (i==N-3) THEN
          a(i)=Bel*phi_b*alpha(i-1)
          b(i)=-Ael*phi_a*alpha(i)-Bel*phi_b*(alpha(i)-beta(i-1))
          c(i)=Ael*phi_a*(alpha(i+1)-beta(i))-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))-Bel*phi_b*epsilonn(i)
          l(i)=0
       ELSEIF (i==N-2) THEN
          a(i)=Bel*phi_b*alpha(i-1)
          b(i)=-Ael*phi_a*alpha(i)-Bel*phi_b*(alpha(i)-beta(i-1))
          c(i)=Ael*phi_a*(alpha(i+1)-beta(i))-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=0;l(i)=0
       ELSEIF (i==N-1) THEN
          a(i)=Bel*phi_b*alpha(i-1)
          b(i)=-Ael*phi_a*alpha(i)-Bel*phi_b*(alpha(i)-beta(i-1))
          c(i)=Ael*phi_a*(alpha(i+1)-beta(i))-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=0.d0;k(i)=0.d0;l(i)=0.d0

       ELSEIF (i==N) THEN
          a(i)=Bel*phi_b*alpha(i-1)
          b(i)=-Bel*phi_b*(alpha(i)-beta(i-1))
          c(i)=-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=-Bel*dphib_dhi*P2-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=0.d0;g(i)=0.d0;k(i)=0.d0;l(i)=0.d0
       ELSE
          a(i)=Bel*phi_b*alpha(i-1)
          b(i)=-Ael*phi_a*alpha(i)-Bel*phi_b*(alpha(i)-beta(i-1))
          c(i)=Ael*phi_a*(alpha(i+1)-beta(i))-Bel*phi_b*(beta(i)-gamma(i-1))
          d(i)=-Bel*dphib_dhi1*P2+Ael*phi_a*(beta(i+1)-gamma(i))-Bel*phi_b*(gamma(i)-lambda(i-1))&
               &+Bgrav*(-dphib_dhi1*H2+phi_b)
          e(i)=Ael*dphia_dhi*P1-Bel*dphib_dhi*P2+Ael*phi_a*(gamma(i+1)-lambda(i))-Bel*phi_b*(lambda(i)-kappa(i-1))&
               &+Agrav*(dphia_dhi*H1-phi_a)-Bgrav*(dphib_dhi*H2+phi_b)
          f(i)=Ael*dphia_dhi1*P1+Ael*phi_a*(lambda(i+1)-kappa(i))-Bel*phi_b*(kappa(i)-delta(i-1))&
               &+Agrav*(dphia_dhi1*H1+phi_a)
          g(i)=Ael*phi_a*(kappa(i+1)-delta(i))-Bel*phi_b*(delta(i)-epsilonn(i-1))
          k(i)=Ael*phi_a*(delta(i+1)-epsilonn(i))-Bel*phi_b*epsilonn(i)
          l(i)=Ael*phi_a*epsilonn(i+1)
       END IF IF3

    END DO

    DEALLOCATE(alpha,beta,gamma,lambda,kappa,delta,epsilonn)
  END SUBROUTINE JACOBI_THICKNESS

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


  SUBROUTINE TRIDIAG(A,B,C,S,N,U,CODE)
    !*****************************************************************
    ! Solves for a vector U of length N the tridiagonal linear set
    ! M U = R, where A, B and C are the three main diagonals of matrix
    ! M(N,N), the other terms are 0. R is the right side vector.
    !*****************************************************************
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: A,B,C,S
    DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: U
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(OUT) :: CODE

    DOUBLE PRECISION, DIMENSION(N) :: GAM
    DOUBLE PRECISION :: BET
    INTEGER :: j

    IF(B(1).EQ.0.D0) THEN
       CODE=1
       RETURN
    END IF

    BET = B(1)
    U(1) = S(1)/BET
    DO J=2,N                    !Decomposition and forward substitution
       GAM(j)=C(j-1)/BET
       BET=B(j)-A(j)*GAM(j)
       IF(BET.EQ.0.D0) THEN            !Algorithm fails
          CODE=2
          RETURN
       END IF
       U(j)=(S(j)-A(j)*U(j-1))/BET
    END DO

    DO J=N-1,1,-1                     !Back substitution
       U(J)=U(J)-GAM(J+1)*U(J+1)
    END DO

    CODE=0
    RETURN
  END SUBROUTINE TRIDIAG

END MODULE MODULE_THICKNESS_SKIN_NEWTON

