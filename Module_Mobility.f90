MODULE MODULE_MOBILITY
CONTAINS
  
  SUBROUTINE CoeffA_THickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a)

    IMPLICIT NONE

    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn, INTENT(INOUT) :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn, INTENT(INOUT) :: Ts_a,Ts_a2,Ts_a3,Ds_a

    h_a=0.5d0*(H(i+1,col)+H(i,col))
    h_a2=0.5d0*(H(i+1,col)**2+H(i,col)**2)
    h_a3=0.5d0*(H(i+1,col)**3+H(i,col)**3)

    delta_a=0.5d0*(BL(i+1,3)+BL(i,3))
    delta_a2=0.5d0*(BL(i+1,3)**2+BL(i,3)**2)
    delta_a3=0.5d0*(BL(i+1,3)**3+BL(i,3)**3)

    T_a = 0.5d0*(T(i,3)+T(i+1,3))
    T_a2 = 0.5d0*(T(i,3)**2+T(i+1,3)**2)
    T_a3 = 0.5d0*(T(i,3)**3+T(i+1,3)**3)

    Ts_a = 0.5d0*(Ts(i,3)+Ts(i+1,3))
    Ts_a2 = 0.5d0*(Ts(i,3)**2+Ts(i+1,3)**2)
    Ts_a3 = 0.5d0*(Ts(i,3)**3+Ts(i+1,3)**3)

    Ds_a = T_a-Ts_a

  END SUBROUTINE CoeffA_THickness

  
  SUBROUTINE CoeffA_Thermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a,eta_a)

    IMPLICIT NONE

    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn, INTENT(INOUT) :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn, INTENT(INOUT) :: Ts_a,Ts_a2,Ts_a3,Ds_a
    DOUBLE PRECISION, INTENT(INOUT) :: eta_a

    eta_a=(grav*(H(i+1,3)-H(i,3))+el*(P(i+1,3)-P(i,3)))/Dr

    h_a=0.5d0*(H(i+1,3)+H(i,3))
    h_a2=0.5d0*(H(i+1,3)**2+H(i,3)**2)
    h_a3=0.5d0*(H(i+1,3)**3+H(i,3)**3)

    delta_a=0.5d0*(BL(i+1,col)+BL(i,col))
    delta_a2=0.5d0*(BL(i+1,col)**2+BL(i,col)**2)
    delta_a3=0.5d0*(BL(i+1,col)**3+BL(i,col)**3)

    T_a = 0.5d0*(T(i,col)+T(i+1,col))
    T_a2 = 0.5d0*(T(i,col)**2+T(i+1,col)**2)
    T_a3 = 0.5d0*(T(i,col)**3+T(i+1,col)**3)

    Ts_a = 0.5d0*(Ts(i,col)+Ts(i+1,col))
    Ts_a2 = 0.5d0*(Ts(i,col)**2+Ts(i+1,col)**2)
    Ts_a3 = 0.5d0*(Ts(i,col)**3+Ts(i+1,col)**3)

    Ds_a = T_a-Ts_a

  END SUBROUTINE CoeffA_Thermal

  SUBROUTINE CoeffB_Thickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b)

    IMPLICIT NONE

    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn, INTENT(INOUT) :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn, INTENT(INOUT) :: Ts_b,Ts_b2,Ts_b3,Ds_b

    h_b=0.5d0*(H(i-1,col)+H(i,col))
    h_b2=0.5d0*(H(i-1,col)**2+H(i,col)**2)
    h_b3=0.5d0*(H(i-1,col)**3+H(i,col)**3)

    delta_b=0.5d0*(BL(i-1,3)+BL(i,3))
    delta_b2=0.5d0*(BL(i-1,3)**2+BL(i,3)**2)
    delta_b3=0.5d0*(BL(i-1,3)**3+BL(i,3)**3)

    T_b = 0.5d0*(T(i,3)+T(i-1,3))
    T_b2 = 0.5d0*(T(i,3)**2+T(i-1,3)**2)
    T_b3 = 0.5d0*(T(i,3)**3+T(i-1,3)**3)

    Ts_b = 0.5d0*(Ts(i,3)+Ts(i-1,3))
    Ts_b2 = 0.5d0*(Ts(i,3)**2+Ts(i-1,3)**2)
    Ts_b3 = 0.5d0*(Ts(i,3)**3+Ts(i-1,3)**3)

    Ds_b = T_b-Ts_b
    
  END SUBROUTINE CoeffB_Thickness
  
  SUBROUTINE CoeffB_Thermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b,eta_b)

    IMPLICIT NONE

    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn, INTENT(INOUT) :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn, INTENT(INOUT) :: Ts_b,Ts_b2,Ts_b3,Ds_b
    DOUBLE PRECISION, INTENT(INOUT) :: eta_b

    eta_b=(grav*(H(i,3)-H(i-1,3))+el*(P(i,3)-P(i-1,3)))/Dr

    h_b=0.5d0*(H(i-1,3)+H(i,3))
    h_b2=0.5d0*(H(i-1,3)**2+H(i,3)**2)
    h_b3=0.5d0*(H(i-1,3)**3+H(i,3)**3)

    delta_b=0.5d0*(BL(i-1,col)+BL(i,col))
    delta_b2=0.5d0*(BL(i-1,col)**2+BL(i,col)**2)
    delta_b3=0.5d0*(BL(i-1,col)**3+BL(i,col)**3)

    T_b = 0.5d0*(T(i,col)+T(i-1,col))
    T_b2 = 0.5d0*(T(i,col)**2+T(i-1,col)**2)
    T_b3 = 0.5d0*(T(i,col)**3+T(i-1,col)**3)

    Ts_b = 0.5d0*(Ts(i,col)+Ts(i-1,col))
    Ts_b2 = 0.5d0*(Ts(i,col)**2+Ts(i-1,col)**2)
    Ts_b3 = 0.5d0*(Ts(i,col)**3+Ts(i-1,col)**3)

    Ds_b = T_b-Ts_b
    
  END SUBROUTINE CoeffB_Thermal

  SUBROUTINE fAi_thermal(ray,dist,Dr,i,Ai)
    IMPLICIT NONE
    
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION , INTENT(IN) :: Dr
    INTEGER , INTENT(IN) :: i
    DOUBLE PRECISION , INTENT(INOUT) :: Ai

    Ai=(ray(i)/(dist(i)*Dr))
    
  END SUBROUTINE fAi_thermal

  SUBROUTINE fBi_thermal(ray,dist,Dr,i,Bi)
    IMPLICIT NONE
    
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION , INTENT(IN) :: Dr
    INTEGER , INTENT(IN) :: i
    DOUBLE PRECISION , INTENT(INOUT) :: Bi

    Bi = (ray(i-1)/(dist(i)*Dr))
    
  END SUBROUTINE fBi_thermal

  SUBROUTINE fAi_thickness(ray,dist,Dr,i,Ai)
    IMPLICIT NONE
    
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION , INTENT(IN) :: Dr
    INTEGER , INTENT(IN) :: i
    DOUBLE PRECISION , INTENT(INOUT) :: Ai

    Ai=(ray(i)/(dist(i)*Dr**2))
    
  END SUBROUTINE fAi_thickness

  SUBROUTINE fBi_thickness(ray,dist,Dr,i,Bi)
    IMPLICIT NONE
    
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray
    DOUBLE PRECISION , INTENT(IN) :: Dr
    INTEGER , INTENT(IN) :: i
    DOUBLE PRECISION , INTENT(INOUT) :: Bi

    Bi = (ray(i-1)/(dist(i)*Dr**2))
    
  END SUBROUTINE fBi_thickness

  SUBROUTINE fomega_a(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,omega_a)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn :: Ts_a,Ts_a2,Ts_a3,Ds_a
    DOUBLE PRECISION :: eta_a

    DOUBLE PRECISION, INTENT(INOUT) :: omega_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffA_Thermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a,eta_a)

    IF (Rheology == 0) THEN
       omega_a = (7.0d0/5.0d0)*T_a*delta_a**2*eta_a*nu - 7.0d0/5.0d0*T_a*&
            & delta_a**2*eta_a - 3.0d0/2.0d0*T_a*delta_a*eta_a*h_a*nu + (3.0d0/&
            & 2.0d0)*T_a*delta_a*eta_a*h_a + (3.0d0/5.0d0)*Ts_a*delta_a**2*&
            & eta_a*nu - 3.0d0/5.0d0*Ts_a*delta_a**2*eta_a - 3.0d0/2.0d0*Ts_a*&
            & delta_a*eta_a*h_a*nu + (3.0d0/2.0d0)*Ts_a*delta_a*eta_a*h_a - 2*&
            & delta_a**2*eta_a*nu + 3*delta_a*eta_a*h_a*nu

       ! omega_a = (eta_a*delta_a)/10.d0*(nu*(-20.d0*delta_a+30.d0*h_a)+&
       !      &(1.d0-nu)*(6.d0*Ds_a*delta_a-15.d0*Ds_a*h_a-20.d0*T_a*delta_a+30.d0*T_a*h_a))

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          omega_a = 3*sqrt(pi)*Ds_a**(-1.5d0)*delta_a2*eta_a*nu*nu**(-T_a)&
               & *(-log(nu))**(-1.5d0)*erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*1.0/Ds_a&
               & *delta_a2*eta_a*nu*nu**(-T_a)*1.0/(-log(nu)) - 3*1.0/Ds_a*&
               & delta_a*eta_a*h_a*nu*nu**Ds_a*nu**(-T_a)*1.0/(-log(nu)) + 3*1.0/&
               & Ds_a*delta_a*eta_a*h_a*nu*nu**(-T_a)*1.0/(-log(nu))

       ELSEIF (nu==1.0) THEN
          omega_a = -2.0d0*delta_a2*eta_a + 3.0d0*delta_a*eta_a*h_a
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fomega_A

  SUBROUTINE fsigma_a(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,sigma_a)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn :: Ts_a,Ts_a2,Ts_a3,Ds_a
    DOUBLE PRECISION :: eta_a

    DOUBLE PRECISION, INTENT(INOUT) :: sigma_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffA_Thermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a,eta_a)

    IF (Rheology == 0) THEN
       sigma_a = -38.0d0/105.0d0*T_a**2*delta_a**3*eta_a*nu + (38.0d0/&
            & 105.0d0)*T_a**2*delta_a**3*eta_a + (1.0d0/3.0d0)*T_a**2*delta_a**&
            & 2*eta_a*h_a*nu - 1.0d0/3.0d0*T_a**2*delta_a**2*eta_a*h_a + (9.0d0&
            & /35.0d0)*T_a*Ts_a*delta_a**3*eta_a*nu - 9.0d0/35.0d0*T_a*Ts_a*&
            & delta_a**3*eta_a - 1.0d0/6.0d0*T_a*Ts_a*delta_a**2*eta_a*h_a*nu +&
            & (1.0d0/6.0d0)*T_a*Ts_a*delta_a**2*eta_a*h_a + (7.0d0/15.0d0)*T_a*&
            & delta_a**3*eta_a*nu - 1.0d0/2.0d0*T_a*delta_a**2*eta_a*h_a*nu + (&
            & 11.0d0/105.0d0)*Ts_a**2*delta_a**3*eta_a*nu - 11.0d0/105.0d0*Ts_a&
            & **2*delta_a**3*eta_a - 1.0d0/6.0d0*Ts_a**2*delta_a**2*eta_a*h_a*&
            & nu + (1.0d0/6.0d0)*Ts_a**2*delta_a**2*eta_a*h_a - 7.0d0/15.0d0*&
            & Ts_a*delta_a**3*eta_a*nu + (1.0d0/2.0d0)*Ts_a*delta_a**2*eta_a*&
            & h_a*nu


       ! sigma_a = (-1.d0/210.d0)*Ds_a*delta_a2*eta_a*(nu*(-98.d0*delta_a+105.d0*h_a)+&
       !      &(1-nu)*(22.d0*Ds_a*delta_a-35.d0*Ds_a*h_a-98.d0*T_a*delta_a+105.d0*T_a*h_a))
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
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

       ELSEIF (nu==1.0) THEN
          sigma_a = (7.0d0/15.0d0)*Ds_a*delta_a3*eta_a - 1.0d0/2.0d0*Ds_a*&
               & delta_a2*eta_a*h_a
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_A

    SUBROUTINE fomega_b(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,omega_b)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn :: Ts_b,Ts_b2,Ts_b3,Ds_b
    DOUBLE PRECISION :: eta_b

    DOUBLE PRECISION, INTENT(INOUT) :: omega_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffB_THermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b,eta_b)


    IF (Rheology == 0) THEN
       omega_b = (7.0d0/5.0d0)*T_b*delta_b**2*eta_b*nu - 7.0d0/5.0d0*T_b*&
            & delta_b**2*eta_b - 3.0d0/2.0d0*T_b*delta_b*eta_b*h_b*nu + (3.0d0/&
            & 2.0d0)*T_b*delta_b*eta_b*h_b + (3.0d0/5.0d0)*Ts_b*delta_b**2*&
            & eta_b*nu - 3.0d0/5.0d0*Ts_b*delta_b**2*eta_b - 3.0d0/2.0d0*Ts_b*&
            & delta_b*eta_b*h_b*nu + (3.0d0/2.0d0)*Ts_b*delta_b*eta_b*h_b - 2*&
            & delta_b**2*eta_b*nu + 3*delta_b*eta_b*h_b*nu

       ! omega_b = (eta_b*delta_b)/10.d0*(nu*(-20.d0*delta_b+30.d0*h_b)+&
       !      &(1.d0-nu)*(6.d0*Ds_b*delta_b-15.d0*Ds_b*h_b-20.d0*T_b*delta_b+30.d0*T_b*h_b))

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
              omega_b = 3*sqrt(pi)*Ds_b**(-1.5d0)*delta_b2*eta_b*nu*nu**(-T_b)&
             & *(-log(nu))**(-1.5d0)*erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*1.0/Ds_b&
             & *delta_b2*eta_b*nu*nu**(-T_b)*1.0/(-log(nu)) - 3*1.0/Ds_b*&
             & delta_b*eta_b*h_b*nu*nu**Ds_b*nu**(-T_b)*1.0/(-log(nu)) + 3*1.0/&
             & Ds_b*delta_b*eta_b*h_b*nu*nu**(-T_b)*1.0/(-log(nu))        

       ELSEIF (nu==1.0) THEN
          omega_b = -2.0d0*delta_b2*eta_b + 3.0d0*delta_b*eta_b*h_b
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FOmega_B


  SUBROUTINE fsigma_b(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,sigma_b)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn :: Ts_b,Ts_b2,Ts_b3,Ds_b
    DOUBLE PRECISION :: eta_b

    DOUBLE PRECISION, INTENT(INOUT) :: sigma_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffB_Thermal(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b,eta_b)

    
    IF (Rheology == 0) THEN
       sigma_b = -38.0d0/105.0d0*T_b**2*delta_b**3*eta_b*nu + (38.0d0/&
            & 105.0d0)*T_b**2*delta_b**3*eta_b + (1.0d0/3.0d0)*T_b**2*delta_b**&
            & 2*eta_b*h_b*nu - 1.0d0/3.0d0*T_b**2*delta_b**2*eta_b*h_b + (9.0d0&
            & /35.0d0)*T_b*Ts_b*delta_b**3*eta_b*nu - 9.0d0/35.0d0*T_b*Ts_b*&
            & delta_b**3*eta_b - 1.0d0/6.0d0*T_b*Ts_b*delta_b**2*eta_b*h_b*nu +&
            & (1.0d0/6.0d0)*T_b*Ts_b*delta_b**2*eta_b*h_b + (7.0d0/15.0d0)*T_b*&
            & delta_b**3*eta_b*nu - 1.0d0/2.0d0*T_b*delta_b**2*eta_b*h_b*nu + (&
            & 11.0d0/105.0d0)*Ts_b**2*delta_b**3*eta_b*nu - 11.0d0/105.0d0*Ts_b&
            & **2*delta_b**3*eta_b - 1.0d0/6.0d0*Ts_b**2*delta_b**2*eta_b*h_b*&
            & nu + (1.0d0/6.0d0)*Ts_b**2*delta_b**2*eta_b*h_b - 7.0d0/15.0d0*&
            & Ts_b*delta_b**3*eta_b*nu + (1.0d0/2.0d0)*Ts_b*delta_b**2*eta_b*&
            & h_b*nu


       ! sigma_b = (-1.d0/210.d0)*Ds_b*delta_b2*eta_b*(nu*(-98.d0*delta_b+105.d0*h_b)+&
       !      &(1-nu)*(22.d0*Ds_b*delta_b-35.d0*Ds_b*h_b-98.d0*T_b*delta_b+105.d0*T_b*h_b))
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
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
       ELSEIF (nu==1.0) THEN
             sigma_b = (7.0d0/15.0d0)*Ds_b*delta_b**3*eta_b - 1.0d0/2.0d0*Ds_b*&
                  & delta_b**2*eta_b*h_b
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_B

    SUBROUTINE fPhi_A(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,phi_a)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn :: Ts_a,Ts_a2,Ts_a3,Ds_a

    DOUBLE PRECISION, INTENT(INOUT) :: phi_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffA_Thickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a)

    
    IF (Rheology == 0) THEN
       phi_a = (4.0d0/5.0d0)*T_a*delta_a**3*nu - 4.0d0/5.0d0*T_a*delta_a&
            & **3 - 2*T_a*delta_a**2*h_a*nu + 2*T_a*delta_a**2*h_a + 2*T_a*&
            & delta_a*h_a**2*nu - 2*T_a*delta_a*h_a**2 - T_a*h_a**3*nu + T_a*&
            & h_a**3 - 4.0d0/5.0d0*Ts_a*delta_a**3*nu + (4.0d0/5.0d0)*Ts_a*&
            & delta_a**3 + 2*Ts_a*delta_a**2*h_a*nu - 2*Ts_a*delta_a**2*h_a - 2&
            & *Ts_a*delta_a*h_a**2*nu + 2*Ts_a*delta_a*h_a**2 + h_a**3*nu
       
       ! phi_a = (nu+(1.d0-nu)*T_a)*h_a3-(2.d0/5.d0)*(1.d0-nu)*&
       !      &Delta_T_a*(2.d0*delta_a3-5.d0*delta_a2*h_a+5.d0*delta_a*h_a2)
       
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          phi_a = 6*sqrt(pi)*Ds_a**(-1.5d0)*delta_a3*nu*nu**(-T_a)*(-log(&
               & nu))**(-1.5d0)*erf(sqrt(Ds_a)*sqrt(-log(nu))) + 12*1.0/Ds_a*&
               & delta_a3*nu*nu**Ds_a*nu**(-T_a)*1.0/(-log(nu)) - 24*1.0/Ds_a*&
               & delta_a3*nu*nu**(-T_a)*1.0/(-log(nu)) - 12*1.0/Ds_a*delta_a2*&
               & h_a*nu*nu**Ds_a*nu**(-T_a)*1.0/(-log(nu)) + 12*1.0/Ds_a*delta_a2&
               & *h_a*nu*nu**(-T_a)*1.0/(-log(nu)) - 12*sqrt(pi)*Ds_a**(-0.5d0)*&
               & delta_a3*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
               & Ds_a)*sqrt(-log(nu))) + 12*sqrt(pi)*Ds_a**(-0.5d0)*delta_a2*h_a&
               & *nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(Ds_a)*sqrt(&
               & -log(nu))) - 3*sqrt(pi)*Ds_a**(-0.5d0)*delta_a*h_a**2*nu*nu**(&
               & -T_a)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(Ds_a)*sqrt(-log(nu&
               & ))) + 8*delta_a3*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) - 12*&
               & delta_a2*h_a*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) + 6*delta_a*&
               & h_a**2*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) - h_a**3*nu*nu**(-T_a&
               & )*1.0/(-log(nu))*log(nu)
       ELSEIF (nu==1.0) THEN
          phi_a = 1.0d0*h_a**3
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_A

    SUBROUTINE  fdPhi_A(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,dphia_dhi,dphia_dhi1)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3
    DOUBLE PRECISIOn :: Ts_a,Ts_a2,Ts_a3,Ds_a
    DOUBLE PRECISION :: hi,hia
    
    DOUBLE PRECISION, INTENT(INOUT) :: dphia_dhi,dphia_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffA_Thickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3,&
       &Ts_a,Ts_a2,Ts_a3,Ds_a)

     hia = H(i+1,col)
     hi = H(i,col)
    
    IF (Rheology == 0) THEN
       
          dphia_dhi1 = -T_a*delta_a**2*nu + T_a*delta_a**2 + T_a*delta_a*hi*&
               & nu - T_a*delta_a*hi + T_a*delta_a*hia*nu - T_a*delta_a*hia -&
               & 3.0d0/8.0d0*T_a*hi**2*nu + (3.0d0/8.0d0)*T_a*hi**2 - 3.0d0/4.0d0*&
               & T_a*hi*hia*nu + (3.0d0/4.0d0)*T_a*hi*hia - 3.0d0/8.0d0*T_a*hia**2&
               & *nu + (3.0d0/8.0d0)*T_a*hia**2 + Ts_a*delta_a**2*nu - Ts_a*&
               & delta_a**2 - Ts_a*delta_a*hi*nu + Ts_a*delta_a*hi - Ts_a*delta_a*&
               & hia*nu + Ts_a*delta_a*hia + (3.0d0/8.0d0)*hi**2*nu + (3.0d0/4.0d0&
               & )*hi*hia*nu + (3.0d0/8.0d0)*hia**2*nu
          dphia_dhi = -T_a*delta_a**2*nu + T_a*delta_a**2 + T_a*delta_a*hi*&
               & nu - T_a*delta_a*hi + T_a*delta_a*hia*nu - T_a*delta_a*hia -&
               & 3.0d0/8.0d0*T_a*hi**2*nu + (3.0d0/8.0d0)*T_a*hi**2 - 3.0d0/4.0d0*&
               & T_a*hi*hia*nu + (3.0d0/4.0d0)*T_a*hi*hia - 3.0d0/8.0d0*T_a*hia**2&
               & *nu + (3.0d0/8.0d0)*T_a*hia**2 + Ts_a*delta_a**2*nu - Ts_a*&
               & delta_a**2 - Ts_a*delta_a*hi*nu + Ts_a*delta_a*hi - Ts_a*delta_a*&
               & hia*nu + Ts_a*delta_a*hia + (3.0d0/8.0d0)*hi**2*nu + (3.0d0/4.0d0&
               & )*hi*hia*nu + (3.0d0/8.0d0)*hia**2*nu

         ! dphia_dhi1=(3.d0/2.d0)*(nu+(1-nu)*T_a)*hia**2-(1-nu)*Delta_T_a*(-delta_a2+delta_a*hia) ! d(phi_i+1/2)/d(h_i+1)
         ! dphia_dhi=(3.d0/2.d0)*(nu+(1-nu)*T_a)*hi**2-(1-nu)*Delta_T_a*(-delta_a2+delta_a*hi)     ! d(phi_i+1/2)/d(h_i)
         
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
         dphia_dhi1 = -12*1.0/Ds_a*delta_a2*nu*nu**Ds_a*nu**(-T_a)*1.0/(&
                  & -log(nu)) + 12*1.0/Ds_a*delta_a2*nu*nu**(-T_a)*1.0/(-log(nu)) +&
                  & 12*sqrt(pi)*Ds_a**(-0.5d0)*delta_a2*nu*nu**(-T_a)*(-log(nu))**(&
                  & -1.5d0)*log(nu)*erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_a&
                  & **(-0.5d0)*delta_a*hi*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*&
                  & erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_a**(-0.5d0)*&
                  & delta_a*hia*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
                  & Ds_a)*sqrt(-log(nu))) - 12*delta_a2*nu*nu**(-T_a)*1.0/(-log(nu&
                  & ))*log(nu) + 12*delta_a*hi*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) +&
                  & 12*delta_a*hia*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) - 3*hi**2*nu*&
                  & nu**(-T_a)*1.0/(-log(nu))*log(nu) - 6*hi*hia*nu*nu**(-T_a)*1.0/(&
                  & -log(nu))*log(nu) - 3*hia**2*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu)

             dphia_dhi = -12*1.0/Ds_a*delta_a2*nu*nu**Ds_a*nu**(-T_a)*1.0/(&
                  & -log(nu)) + 12*1.0/Ds_a*delta_a2*nu*nu**(-T_a)*1.0/(-log(nu)) +&
                  & 12*sqrt(pi)*Ds_a**(-0.5d0)*delta_a2*nu*nu**(-T_a)*(-log(nu))**(&
                  & -1.5d0)*log(nu)*erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_a&
                  & **(-0.5d0)*delta_a*hi*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*&
                  & erf(sqrt(Ds_a)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_a**(-0.5d0)*&
                  & delta_a*hia*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
                  & Ds_a)*sqrt(-log(nu))) - 12*delta_a2*nu*nu**(-T_a)*1.0/(-log(nu&
                  & ))*log(nu) + 12*delta_a*hi*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) +&
                  & 12*delta_a*hia*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu) - 3*hi**2*nu*&
                  & nu**(-T_a)*1.0/(-log(nu))*log(nu) - 6*hi*hia*nu*nu**(-T_a)*1.0/(&
                  & -log(nu))*log(nu) - 3*hia**2*nu*nu**(-T_a)*1.0/(-log(nu))*log(nu)
             
       ELSEIF (nu==1.0) THEN
          dphia_dhi1 = 3*hi**2 + 6*hi*hia + 3*hia**2
          dphia_dhi = 3*hi**2 + 6*hi*hia + 3*hia**2
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FdPhi_A

   SUBROUTINE  fPhi_B(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,phi_b)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn :: Ts_b,Ts_b2,Ts_b3,Ds_b

    DOUBLE PRECISION, INTENT(INOUT) :: phi_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffB_Thickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b)

    
    IF (Rheology == 0) THEN
       
       phi_b = (4.0d0/5.0d0)*T_b*delta_b**3*nu - 4.0d0/5.0d0*T_b*delta_b&
            & **3 - 2*T_b*delta_b**2*h_b*nu + 2*T_b*delta_b**2*h_b + 2*T_b*&
            & delta_b*h_b**2*nu - 2*T_b*delta_b*h_b**2 - T_b*h_b**3*nu + T_b*&
            & h_b**3 - 4.0d0/5.0d0*Ts_b*delta_b**3*nu + (4.0d0/5.0d0)*Ts_b*&
            & delta_b**3 + 2*Ts_b*delta_b**2*h_b*nu - 2*Ts_b*delta_b**2*h_b - 2&
            & *Ts_b*delta_b*h_b**2*nu + 2*Ts_b*delta_b*h_b**2 + h_b**3*nu
       
       ! phi_b=(nu+(1.d0-nu)*T_b)*h_b3-(2.d0/5.d0)*(1.d0-nu)*&
       !      &Delta_T_b*(2.d0*delta_b3-5.d0*delta_b2*h_b+5.d0*delta_b*h_b2)
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          phi_b = 6*sqrt(pi)*Ds_b**(-1.5d0)*delta_b3*nu*nu**(-T_b)*(-log(&
               & nu))**(-1.5d0)*erf(sqrt(Ds_b)*sqrt(-log(nu))) + 12*1.0/Ds_b*&
               & delta_b3*nu*nu**Ds_b*nu**(-T_b)*1.0/(-log(nu)) - 24*1.0/Ds_b*&
               & delta_b3*nu*nu**(-T_b)*1.0/(-log(nu)) - 12*1.0/Ds_b*delta_b2*&
               & h_b*nu*nu**Ds_b*nu**(-T_b)*1.0/(-log(nu)) + 12*1.0/Ds_b*delta_b2&
               & *h_b*nu*nu**(-T_b)*1.0/(-log(nu)) - 12*sqrt(pi)*Ds_b**(-0.5d0)*&
               & delta_b3*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
               & Ds_b)*sqrt(-log(nu))) + 12*sqrt(pi)*Ds_b**(-0.5d0)*delta_b2*h_b&
               & *nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(Ds_b)*sqrt(&
               & -log(nu))) - 3*sqrt(pi)*Ds_b**(-0.5d0)*delta_b*h_b**2*nu*nu**(&
               & -T_b)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(Ds_b)*sqrt(-log(nu&
               & ))) + 8*delta_b3*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) - 12*&
               & delta_b2*h_b*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) + 6*delta_b*&
               & h_b**2*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) - h_b**3*nu*nu**(-T_b&
               & )*1.0/(-log(nu))*log(nu)
       ELSEIF (nu==1.0) THEN
          phi_b = 1.0d0*h_b**3
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_B
  
  SUBROUTINE  fDPhi_B(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i&
       &,nu,Rheology,ERROR_CODE,dphib_dhi,dphib_dhi1)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION , INTENT(IN) :: Dr,Dt,el,grav
    INTEGER, INTENT(IN) :: col,i
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(IN) :: H,T,Ts,BL,P
    DOUBLE PRECISION ,DIMENSION(:), INTENT(IN) :: dist,ray

    DOUBLE PRECISIOn :: h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3
    DOUBLE PRECISIOn :: Ts_b,Ts_b2,Ts_b3,Ds_b

    DOUBLE PRECISION :: hi,hib

    DOUBLE PRECISION, INTENT(INOUT) :: dphib_dhi,dphib_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    CALL CoeffB_Thickness(H,T,Ts,BL,P,col,dist,ray,Dr,Dt,el,grav,i,&
       &h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3,&
       &Ts_b,Ts_b2,Ts_b3,Ds_b)

    hib = H(i-1,col)
    hi = H(i,col)
    
    IF (Rheology == 0) THEN
       
       dphib_dhi = -T_b*delta_b**2*nu + T_b*delta_b**2 + T_b*delta_b*hi*&
            & nu - T_b*delta_b*hi + T_b*delta_b*hib*nu - T_b*delta_b*hib -&
            & 3.0d0/8.0d0*T_b*hi**2*nu + (3.0d0/8.0d0)*T_b*hi**2 - 3.0d0/4.0d0*&
            & T_b*hi*hib*nu + (3.0d0/4.0d0)*T_b*hi*hib - 3.0d0/8.0d0*T_b*hib**2&
            & *nu + (3.0d0/8.0d0)*T_b*hib**2 + Ts_b*delta_b**2*nu - Ts_b*&
            & delta_b**2 - Ts_b*delta_b*hi*nu + Ts_b*delta_b*hi - Ts_b*delta_b*&
            & hib*nu + Ts_b*delta_b*hib + (3.0d0/8.0d0)*hi**2*nu + (3.0d0/4.0d0&
            & )*hi*hib*nu + (3.0d0/8.0d0)*hib**2*nu
       dphib_dhi1 = -T_b*delta_b**2*nu + T_b*delta_b**2 + T_b*delta_b*hi*&
            & nu - T_b*delta_b*hi + T_b*delta_b*hib*nu - T_b*delta_b*hib -&
            & 3.0d0/8.0d0*T_b*hi**2*nu + (3.0d0/8.0d0)*T_b*hi**2 - 3.0d0/4.0d0*&
            & T_b*hi*hib*nu + (3.0d0/4.0d0)*T_b*hi*hib - 3.0d0/8.0d0*T_b*hib**2&
            & *nu + (3.0d0/8.0d0)*T_b*hib**2 + Ts_b*delta_b**2*nu - Ts_b*&
            & delta_b**2 - Ts_b*delta_b*hi*nu + Ts_b*delta_b*hi - Ts_b*delta_b*&
            & hib*nu + Ts_b*delta_b*hib + (3.0d0/8.0d0)*hi**2*nu + (3.0d0/4.0d0&
            & )*hi*hib*nu + (3.0d0/8.0d0)*hib**2*nu
       
       ! dphib_dhi=(3.d0/2.d0)*(nu+(1-nu)*T_b)*hi**2-(1-nu)*Delta_T_b*(-delta_b2+delta_b*hi)     ! d(phi_i-1/2)/d(h_i)
       ! dphib_dhi1=(3.d0/2.d0)*(nu+(1-nu)*T_b)*hib**2-(1-nu)*Delta_T_b*(-delta_b2+delta_b*hib) ! d(phi_i-1/2)/d(h_i-1)
       
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
           dphib_dhi = -12*1.0/Ds_b*delta_b2*nu*nu**Ds_b*nu**(-T_b)*1.0/(&
             & -log(nu)) + 12*1.0/Ds_b*delta_b2*nu*nu**(-T_b)*1.0/(-log(nu)) +&
             & 12*sqrt(pi)*Ds_b**(-0.5d0)*delta_b2*nu*nu**(-T_b)*(-log(nu))**(&
             & -1.5d0)*log(nu)*erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_b&
             & **(-0.5d0)*delta_b*hi*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*&
             & erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_b**(-0.5d0)*&
             & delta_b*hib*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
             & Ds_b)*sqrt(-log(nu))) - 12*delta_b2*nu*nu**(-T_b)*1.0/(-log(nu&
             & ))*log(nu) + 12*delta_b*hi*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) +&
             & 12*delta_b*hib*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) - 3*hi**2*nu*&
             & nu**(-T_b)*1.0/(-log(nu))*log(nu) - 6*hi*hib*nu*nu**(-T_b)*1.0/(&
             & -log(nu))*log(nu) - 3*hib**2*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu)

             dphib_dhi1 = -12*1.0/Ds_b*delta_b2*nu*nu**Ds_b*nu**(-T_b)*1.0/(&
             & -log(nu)) + 12*1.0/Ds_b*delta_b2*nu*nu**(-T_b)*1.0/(-log(nu)) +&
             & 12*sqrt(pi)*Ds_b**(-0.5d0)*delta_b2*nu*nu**(-T_b)*(-log(nu))**(&
             & -1.5d0)*log(nu)*erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_b&
             & **(-0.5d0)*delta_b*hi*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*&
             & erf(sqrt(Ds_b)*sqrt(-log(nu))) - 6*sqrt(pi)*Ds_b**(-0.5d0)*&
             & delta_b*hib*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*log(nu)*erf(sqrt(&
             & Ds_b)*sqrt(-log(nu))) - 12*delta_b2*nu*nu**(-T_b)*1.0/(-log(nu&
             & ))*log(nu) + 12*delta_b*hi*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) +&
             & 12*delta_b*hib*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu) - 3*hi**2*nu*&
             & nu**(-T_b)*1.0/(-log(nu))*log(nu) - 6*hi*hib*nu*nu**(-T_b)*1.0/(&
             & -log(nu))*log(nu) - 3*hib**2*nu*nu**(-T_b)*1.0/(-log(nu))*log(nu)
             
       ELSEIF (nu == 1.0) THEN
          dphib_dhi1 = 3*hi**2 + 6*hi*hib + 3*hib**2
          dphib_dhi = 3*hi**2 + 6*hi*hib + 3*hib**2
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FDPhi_B


    
END MODULE MODULE_MOBILITY

