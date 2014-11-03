MODULE MODULE_COMPLEMENTAIRE

  USE MODULE_INTEGRATION

  
CONTAINS

  SUBROUTINE  STRESS_ELASTIC_FIELD(Srr,Stt,H,dist,Dr,Mtot)

    IMPLICIT NONE

    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist
    DOUBLE PRECISION ,INTENT(IN) :: Dr
    INTEGER ,INTENT(IN) :: Mtot

    DOUBLE PRECISION ,DIMENSION(:), INTENT(INOUT) :: Srr,Stt

    INTEGER :: i
    DOUBLE PRECISION :: nuu

    nuu=0.25
    DO i=1,Mtot-2
       IF (i==1) THEN
          Srr(i)=-(1.d0/(12*Dr**2*dist(i)))*((nuu*Dr-dist(i))*H(i+1,3)+(16*dist(i)-8*nuu*Dr)*H(i,3)&
               &-30*dist(i)*H(i,3)+(16*dist(i)+8*nuu*Dr)*H(i+1,3)-(dist(i)+nuu*Dr)*H(i+2,3))
          Stt(i)=-(1.d0/(12*Dr**2*dist(i)))*((Dr-nuu*dist(i))*H(i+1,3)+(16*nuu*dist(i)-8*Dr)*H(i,3)&
               &-30*nuu*dist(i)*H(i,3)+(16*nuu*dist(i)+8*Dr)*H(i+1,3)-(Dr+nuu*dist(i))*H(i+2,3))
       ELSEIF (i==2) THEN
          Srr(i)=-(1.d0/(12*Dr**2*dist(i)))*((nuu*Dr-dist(i))*H(i-1,3)+(16*dist(i)-8*nuu*Dr)*H(i-1,3)&
               &-30*dist(i)*H(i,3)+(16*dist(i)+8*nuu*Dr)*H(i+1,3)-(dist(i)+nuu*Dr)*H(i+2,3))
          Stt(i)=-(1.d0/(12*Dr**2*dist(i)))*((Dr-nuu*dist(i))*H(i-1,3)+(16*nuu*dist(i)-8*Dr)*H(i-1,3)&
               &-30*nuu*dist(i)*H(i,3)+(16*nuu*dist(i)+8*Dr)*H(i+1,3)-(Dr+nuu*dist(i))*H(i+2,3))
       ELSE
          Srr(i)=-(1.d0/(12*Dr**2*dist(i)))*((nuu*Dr-dist(i))*H(i-2,3)+(16*dist(i)-8*nuu*Dr)*H(i-1,3)&
               &-30*dist(i)*H(i,3)+(16*dist(i)+8*nuu*Dr)*H(i+1,3)-(dist(i)+nuu*Dr)*H(i+2,3))
          Stt(i)=-(1.d0/(12*Dr**2*dist(i)))*((Dr-nuu*dist(i))*H(i-2,3)+(16*nuu*dist(i)-8*Dr)*H(i-1,3)&
               &-30*nuu*dist(i)*H(i,3)+(16*nuu*dist(i)+8*Dr)*H(i+1,3)-(Dr+nuu*dist(i))*H(i+2,3))          
       END IF
    END DO

  END SUBROUTINE STRESS_ELASTIC_FIELD

  SUBROUTINE AVERAGE_QUANTITY(Xi,H,T,Ts,BL,dist,ray,Dt,Dr,el,grav,N1,Pe,Psi,nu,Tm,Vm,Mum,Phim,M,tmps,delta0,&
       &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,Tm01,Tm02,Tm05,Tm005)

    IMPLICIT NONE

    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: XI,H,T,Ts,BL
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,el,grav,N1,Pe,Psi,nu,tmps,delta0
    DOUBLE PRECISION ,INTENT(INOUT) :: Tm,Vm,Mum,Phim
    DOUBLE PRECISION ,INTENT(INOUT) :: Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005
    DOUBLE PRECISION ,INTENT(INOUT) :: Tm01,Tm02,Tm05,Tm005
    INTEGER, INTENT(IN) :: M

    DOUBLE PRECISION, EXTERNAL:: viscosity_1,viscosity_2,viscosity_3
    DOUBLE PRECISION :: hthetabar,a,beta,muPart1,muPart2,muPart3
    DOUBLE PRECISION :: hmubar,Phibar
    DOUBLE PRECISION :: abserr
    INTEGER :: i,N,ier,last
    INTEGER :: I01,I02,I05,I005

    DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v,ho
    common /arg/ delta,Thetas,Thetab,nu_v,ho

    N = 0
    DO i=1,M,1
       IF (H(i,3)-delta0>0.d0 .OR. N /= 0) THEN
          CYCLE
       ELSE
          N = i
       END IF
    ENDDO
    
    i01 = 0; I02 = 0; I05 = 0; I005 = 0
    DO i=1,N,1
       IF (dist(i)>(1D0-0.1D0)*dist(N) .AND. i01 == 0) THEN
          i01 = i
       ENDIF
       IF (dist(i)>(1D0-0.2D0)*dist(N) .AND. i02 == 0) THEN
          I02 = i
       ENDIF
       IF (dist(i)>(1D0-0.5D0)*dist(N) .AND. i05 == 0) THEN
          I05 = i
       ENDIF 
       IF (dist(i)>(1D0-0.05D0)*dist(N) .AND. i005 == 0) THEN
          I005 = i
       ENDIF
    ENDDO

    Vm01 = 0D0;Vm02 = 0D0;Vm05 = 0D0;Vm005 = 0D0
    Mum01 =0D0;Mum02 =0D0;Mum05 =0D0;Mum005 =0D0
    Tm01 = 0D0;Tm02 = 0D0;Tm05 = 0D0;Tm005 = 0D0
    
    DO i=1,N,1
       hthetabar = -2*(T(i,3)-Ts(i,3))/3.d0*BL(i,3)+T(i,3)*H(i,3)
       beta = (1.d0-nu)
       Thetas = Ts(i,3);Thetab = T(i,3);delta = BL(i,3);nu_v = nu;ho=H(i,3)
       CALL qxgs(Viscosity_1,0.d0,delta,1D-6,1D-3,muPart1,abserr,ier,10,last)
       CALL qxgs(Viscosity_2,delta,ho-delta,1D-6,1D-3,muPart2,abserr,ier,10,last)
       CALL qxgs(Viscosity_3,ho-delta,ho,1D-6,1D-3,muPart3,abserr,ier,10,last)
       hmubar = muPart1+muPart2+muPart3
       Phibar = -4.d0*Pe*((T(i,3)-Ts(i,3))/BL(i,3))

       IF (i ==1) THEN
          Tm = hthetabar*ray(i)**2
          Vm = H(i,3)*ray(i)**2
          Mum = hmubar*ray(i)**2
          Phim = Phibar*ray(i)**2
       ELSEIF (i == N) THEN
          Tm = Tm + hthetabar*(dist(N)**2-ray(N-1)**2)
          Vm = Vm +H(i,3)*(dist(N)**2-ray(N-1)**2)
          Mum = Mum + hmubar*(dist(N)**2-ray(N-1)**2)
          Phim = Phim+Phibar*(dist(N)**2-ray(N-1)**2)
          IF (i>=I01) THEN
             Vm01 = Vm01 +H(i,3)*(dist(N)**2-ray(N-1)**2)
             Mum01 = Mum01 + hmubar*(dist(N)**2-ray(N-1)**2)
             Tm01 = Tm01 + hthetabar*(dist(N)**2-ray(N-1)**2)
          ENDIF
          IF (i>=I02) THEN
             Vm02 = Vm02 +H(i,3)*(dist(N)**2-ray(N-1)**2)
             Mum02 = Mum02 + hmubar*(dist(N)**2-ray(N-1)**2)
             Tm02 = Tm02 + hthetabar*(dist(N)**2-ray(N-1)**2)
          ENDIF
          IF (i>=I05) THEN
             Vm05 = Vm05 +H(i,3)*(dist(N)**2-ray(N-1)**2)
             Mum05 = Mum05 + hmubar*(dist(N)**2-ray(N-1)**2)
             Tm05 = Tm05 + hthetabar*(dist(N)**2-ray(N-1)**2)
          ENDIF
          IF (i>=I005) THEN
             Vm005 = Vm005 +H(i,3)*(dist(N)**2-ray(N-1)**2)
             Mum005 = Mum005 + hmubar*(dist(N)**2-ray(N-1)**2)
             Tm005 = Tm005 + hthetabar*(dist(N)**2-ray(N-1)**2)
          ENDIF
       ELSE 
          Tm = Tm + hthetabar*(ray(i)**2-ray(i-1)**2)
          Vm = Vm +H(i,3)*(ray(i)**2-ray(i-1)**2)
          Mum = Mum + hmubar*(ray(i)**2-ray(i-1)**2)
          Phim = Phim+Phibar*(ray(i)**2-ray(i-1)**2)
          IF (i>=I01) THEN
             Vm01 = Vm01 +H(i,3)*(dist(i)**2-ray(i-1)**2)
             Mum01 = Mum01 + hmubar*(dist(i)**2-ray(i-1)**2)
             Tm01 = Tm01 + hthetabar*(dist(i)**2-ray(i-1)**2)
          ENDIF
          IF (i>=I02) THEN
             Vm02 = Vm02 +H(i,3)*(dist(i)**2-ray(i-1)**2)
             Mum02 = Mum02 + hmubar*(dist(i)**2-ray(i-1)**2)
             Tm02 = Tm02 + hthetabar*(dist(i)**2-ray(i-1)**2)
          ENDIF
          IF (i>=I05) THEN
             Vm05 = Vm05 +H(i,3)*(dist(i)**2-ray(i-1)**2)
             Mum05 = Mum05 + hmubar*(dist(i)**2-ray(i-1)**2)
             Tm05 = Tm05 + hthetabar*(dist(i)**2-ray(i-1)**2)
          ENDIF
          IF (i>=I005) THEN
             Vm005 = Vm005 +H(i,3)*(dist(i)**2-ray(i-1)**2)
             Mum005 = Mum005 + hmubar*(dist(i)**2-ray(i-1)**2)
             Tm005 = Tm005 + hthetabar*(dist(i)**2-ray(i-1)**2)
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE AVERAGE_QUANTITY

  SUBROUTINE TRACKING_FRONT(Xi,H,T,Ts,BL,dist,ray,Dt,Dr,el,grav,N1,Pe,Psi,nu,tmps,delta0,&
       &Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e)

    IMPLICIT NONE

    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: XI,H,T,Ts,BL
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,el,grav,N1,Pe,Psi,nu,tmps,delta0
    DOUBLE PRECISION ,INTENT(INOUT) :: Fr_d_R,Fr_d_T,Fr_d_Mu,Mu_e
    DOUBLE PRECISION ,INTENT(INOUT) :: Fr_001_R,Fr_001_T,Fr_001_Mu

    DOUBLE PRECISION, EXTERNAL:: viscosity_1,viscosity_2,viscosity_3
    DOUBLE PRECISION :: a,beta,muPart1,muPart2,muPart3
    DOUBLE PRECISION :: Vm,Tm,Mum,tbar,mbar
    DOUBLE PRECISION , DIMENSION(:),ALLOCATABLE  :: hmubar,hthetabar
    DOUBLE PRECISION :: abserr
    INTEGER :: i,N,ier,last,err1,N001

    DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v,ho
    common /arg/ delta,Thetas,Thetab,nu_v,ho

    N = 0
    DO i=1,COUNT(H(:,3)>0),1
       IF (H(i,3)-delta0>0.d0 .OR. N /= 0) THEN
          CYCLE
       ELSE
          N = i
       END IF
    ENDDO
    
    N001 = 0
    DO i=1,COUNT(H(:,3)>0),1
       IF (H(i,3)-0.01>0.d0 .OR. N001 /= 0) THEN
          CYCLE
       ELSE
          N001 = i
       END IF
    ENDDO
    
    ALLOCATE(hmubar(N),hthetabar(N), stat=err1)
    IF (err1>1) THEN
       PRINT*,'ERREUR endallocation tabelaux dans le main'
       STOP
    ENDIF

    DO i=1,N,1
       hthetabar = -2*(T(i,3)-Ts(i,3))/3.d0*BL(i,3)+T(i,3)*H(i,3)
       beta = (1.d0-nu)
       Thetas = Ts(i,3);Thetab = T(i,3);delta = BL(i,3);nu_v = nu;ho=H(i,3)
       CALL qxgs(Viscosity_1,0.d0,delta,1D-6,1D-3,muPart1,abserr,ier,10,last)
       CALL qxgs(Viscosity_2,delta,ho-delta,1D-6,1D-3,muPart2,abserr,ier,10,last)
       CALL qxgs(Viscosity_3,ho-delta,ho,1D-6,1D-3,muPart3,abserr,ier,10,last)
       hmubar(i) = muPart1+muPart2+muPart3
       hthetabar(i) = -2*(T(i,3)-Ts(i,3))/3.d0*BL(i,3)+T(i,3)*H(i,3)
    ENDDO

    Mu_e = 0.38*delta0**(-1D0/11D0)*(H(1,3)/(tmps**(8D0/22D0)))**(11D0/2D0)
    
    ! Premier cas definit le front avec delta0
    Mum = hmubar(N)*(dist(N)**2-ray(N-1)**2)
    Vm = H(N,3)*(dist(N)**2-ray(N-1)**2)
    Tm = hthetabar(N)*(dist(N)**2-ray(N-1)**2)
    mbar = Mum/Vm
    tbar = Tm/Vm
    Fr_d_R = dist(N)
    Fr_d_T = 0
    Fr_d_Mu = 1D0/nu

    DO i=N,1,-1
       IF (mbar<Mu_e) THEN
          Fr_d_R = dist(i)
          Fr_d_T = tbar
          Fr_d_Mu = mbar
          EXIT
       ELSE
          IF (i/=1) THEN
             Tm = Tm + hthetabar(i)*(ray(i)**2-ray(i-1)**2)
             Vm = Vm +H(i,3)*(ray(i)**2-ray(i-1)**2)
             Mum = Mum + hmubar(i)*(ray(i)**2-ray(i-1)**2)
             mbar = Mum/Vm
             tbar = Tm/Vm
          ELSE
             Tm = Tm+hthetabar(i)*ray(i)**2
             Vm = Vm+H(i,3)*ray(i)**2
             Mum = Mum+hmubar(i)*ray(i)**2
             mbar = Mum/Vm
             tbar = Tm/Vm
          ENDIF
       ENDIF
    ENDDO
    IF (mbar>Mu_e) THEN
        Fr_d_R = 0.d0
        Fr_d_T = tbar
        Fr_d_Mu = mbar
     ENDIF

    ! Deuxieme cas, on definit le front avec 0.001
! Premier cas definit le front avec delta0
    
    IF (N001>1) THEN
       Mum = hmubar(N001)*(ray(N001)**2-ray(N001-1)**2)
       Vm = H(N001,3)*(ray(N001)**2-ray(N001-1)**2)
       Tm = hthetabar(N001)*(ray(N001)**2-ray(N001-1)**2)
       mbar = Mum/Vm
       tbar = Tm/Vm
       Fr_001_R = dist(N001)
       Fr_001_T = 0
       Fr_001_Mu = 1D0/nu
       DO i=N001,1,-1
          IF (mbar<Mu_e) THEN
             Fr_001_R = dist(i)
             Fr_001_T = tbar
             Fr_001_Mu = mbar
             EXIT
          ELSE
             IF (i/=1) THEN
                Tm = Tm + hthetabar(i)*(ray(i)**2-ray(i-1)**2)
                Vm = Vm +H(i,3)*(ray(i)**2-ray(i-1)**2)
                Mum = Mum + hmubar(i)*(ray(i)**2-ray(i-1)**2)
                mbar = Mum/Vm
                tbar = Tm/Vm
             ELSE
                Tm = Tm+hthetabar(i)*ray(i)**2
                Vm = Vm+H(i,3)*ray(i)**2
                Mum = Mum+hmubar(i)*ray(i)**2
             ENDIF
          ENDIF
       ENDDO
       IF (mbar>Mu_e) THEN
          Fr_001_R = 0.d0
          Fr_001_T = tbar
          Fr_001_Mu = mbar
       ENDIF
    ELSE
        Fr_001_R = Fr_d_R
        Fr_001_T = Fr_d_T 
        Fr_001_Mu = Fr_d_Mu
     ENDIF
    
  
DEALLOCATE(hmubar,hthetabar)
       
  ENDSUBROUTINE TRACKING_FRONT


  END MODULE MODULE_COMPLEMENTAIRE

  DOUBLE PRECISION FUNCTION Viscosity_1(x)
    IMPLICIT NONE  
    DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v,ho
    common /arg/ delta,Thetas,Thetab,nu_v,ho
    DOUBLE PRECISION :: x
    
       Viscosity_1  = 1.d0/(nu_v+(1.d0-nu_v)*(Thetab-(Thetab-Thetas)*(1.d0-x/delta)**2))
  END FUNCTION Viscosity_1

  DOUBLE PRECISION FUNCTION Viscosity_2(x)
    IMPLICIT NONE  
    DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v,ho
    common /arg/ delta,Thetas,Thetab,nu_v,ho
    DOUBLE PRECISION :: x
    
       Viscosity_2  = 1.d0/(nu_v+(1.d0-nu_v)*Thetab)
  END FUNCTION Viscosity_2

  DOUBLE PRECISION FUNCTION Viscosity_3(x)
    IMPLICIT NONE  
    DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v,ho
    common /arg/ delta,Thetas,Thetab,nu_v,ho
    DOUBLE PRECISION :: x
    
       Viscosity_3  = 1.d0/(nu_v+(1.d0-nu_v)*(Thetab-(Thetab-Thetas)*(1.d0-(ho-x)/delta)**2))
  END FUNCTION Viscosity_3
