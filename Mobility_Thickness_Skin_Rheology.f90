MODULE MOBILITY_THICKNESS_SKIN_RHEOLOGY

CONTAINS
    
  SUBROUTINE fPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,delta_a2,delta_a3,h_a2,h_a3,Delta_T_a,&
       &phi_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ael,Agrav,h_a,delta_a,T_a,Ts_a
    DOUBLE PRECISION, INTENT(IN) :: delta_a2,delta_a3,h_a2,h_a3,Delta_T_a
    DOUBLE PRECISION, INTENT(INOUT) :: phi_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    
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
          phi_a = -12.0d0*sqrt(pi)*T_a*delta_a**3*nu*nu**(-T_a)*(-log(nu))**&
               & (-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(&
               & T_a - Ts_a)) + 8.0d0*T_a*delta_a**3*nu*nu**(-T_a)*1.0/(-log(nu))*&
               & 1.0/(T_a - Ts_a)*log(nu) + 12.0d0*sqrt(pi)*T_a*delta_a**2*h_a*nu*&
               & nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*&
               & erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 12.0d0*T_a*delta_a**2*h_a*&
               & nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) - 3.0d0*&
               & sqrt(pi)*T_a*delta_a*h_a**2*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(&
               & T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a&
               & )) + 6.0d0*T_a*delta_a*h_a**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(&
               & T_a - Ts_a)*log(nu) - 1.0d0*T_a*h_a**3*nu*nu**(-T_a)*1.0/(-log(nu&
               & ))*1.0/(T_a - Ts_a)*log(nu) + 12.0d0*sqrt(pi)*Ts_a*delta_a**3*nu*&
               & nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*&
               & erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 8.0d0*Ts_a*delta_a**3*nu*&
               & nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) - 12.0d0*sqrt(&
               & pi)*Ts_a*delta_a**2*h_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a -&
               & Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) +&
               & 12.0d0*Ts_a*delta_a**2*h_a*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a&
               & - Ts_a)*log(nu) + 3.0d0*sqrt(pi)*Ts_a*delta_a*h_a**2*nu*nu**(-T_a&
               & )*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(&
               & -log(nu))*sqrt(T_a - Ts_a)) - 6.0d0*Ts_a*delta_a*h_a**2*nu*nu**(&
               & -T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 1.0d0*Ts_a*h_a**3&
               & *nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 12.0d0*&
               & delta_a**3*nu*nu**(-Ts_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a) + 6.0d0&
               & *sqrt(pi)*delta_a**3*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a -&
               & Ts_a)**(-1.5d0)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 24.0d0*&
               & delta_a**3*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a) - 12.0d0&
               & *delta_a**2*h_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a) +&
               & 12.0d0*delta_a**2*h_a*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a -&
               & Ts_a)
       ELSEIF (nu==1.0) THEN
          phi_a = 1.0d0*h_a**3
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_A

  SUBROUTINE  fdPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,hi,hia,Delta_T_a,delta_a2,&
       &dphia_dhi1,dphia_dhi,nu,Rheology,ERROR_CODE)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ael,Agrav,h_a,T_a,delta_a,Ts_a,Delta_T_a
    DOUBLE PRECISION, INTENT(IN) :: hi,hia,delta_a2

    DOUBLE PRECISION, INTENT(INOUT) :: dphia_dhi,dphia_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    
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
          dphia_dhi = 6.0d0*sqrt(pi)*T_a*delta_a**2*nu*nu**(-T_a)*(-log(nu))&
          & **(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt&
          & (T_a - Ts_a)) - 6.0d0*T_a*delta_a**2*nu*nu**(-T_a)*1.0/(-log(nu))&
          & *1.0/(T_a - Ts_a)*log(nu) - 1.5d0*sqrt(pi)*T_a*delta_a*hi*nu*nu**&
          & (-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(&
          & sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*T_a*delta_a*hi*nu*nu**(&
          & -T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) - 1.5d0*sqrt(pi)*&
          & T_a*delta_a*hia*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**&
          & (-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*T_a&
          & *delta_a*hia*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu&
          & ) - 0.375d0*T_a*hi**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a -&
          & Ts_a)*log(nu) - 0.75d0*T_a*hi*hia*nu*nu**(-T_a)*1.0/(-log(nu))*&
          & 1.0/(T_a - Ts_a)*log(nu) - 0.375d0*T_a*hia**2*nu*nu**(-T_a)*1.0/(&
          & -log(nu))*1.0/(T_a - Ts_a)*log(nu) - 6.0d0*sqrt(pi)*Ts_a*delta_a&
          & **2*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log&
          & (nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 6.0d0*Ts_a*delta_a**2&
          & *nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 1.5d0*&
          & sqrt(pi)*Ts_a*delta_a*hi*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a&
          & - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) -&
          & 3.0d0*Ts_a*delta_a*hi*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a -&
          & Ts_a)*log(nu) + 1.5d0*sqrt(pi)*Ts_a*delta_a*hia*nu*nu**(-T_a)*(&
          & -log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(&
          & nu))*sqrt(T_a - Ts_a)) - 3.0d0*Ts_a*delta_a*hia*nu*nu**(-T_a)*1.0&
          & /(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 0.375d0*Ts_a*hi**2*nu*nu**&
          & (-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 0.75d0*Ts_a*hi*&
          & hia*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) +&
          & 0.375d0*Ts_a*hia**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)&
          & *log(nu) - 6.0d0*delta_a**2*nu*nu**(-Ts_a)*1.0/(-log(nu))*1.0/(&
          & T_a - Ts_a) + 6.0d0*delta_a**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(&
          & T_a - Ts_a)


          dphia_dhi1 = 6.0d0*sqrt(pi)*T_a*delta_a**2*nu*nu**(-T_a)*(-log(nu&
          & ))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*&
          & sqrt(T_a - Ts_a)) - 6.0d0*T_a*delta_a**2*nu*nu**(-T_a)*1.0/(-log(&
          & nu))*1.0/(T_a - Ts_a)*log(nu) - 1.5d0*sqrt(pi)*T_a*delta_a*hi*nu*&
          & nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*&
          & erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*T_a*delta_a*hi*nu*nu&
          & **(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) - 1.5d0*sqrt(pi)&
          & *T_a*delta_a*hia*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)&
          & **(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*&
          & T_a*delta_a*hia*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log&
          & (nu) - 0.375d0*T_a*hi**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a -&
          & Ts_a)*log(nu) - 0.75d0*T_a*hi*hia*nu*nu**(-T_a)*1.0/(-log(nu))*&
          & 1.0/(T_a - Ts_a)*log(nu) - 0.375d0*T_a*hia**2*nu*nu**(-T_a)*1.0/(&
          & -log(nu))*1.0/(T_a - Ts_a)*log(nu) - 6.0d0*sqrt(pi)*Ts_a*delta_a&
          & **2*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log&
          & (nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 6.0d0*Ts_a*delta_a**2&
          & *nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 1.5d0*&
          & sqrt(pi)*Ts_a*delta_a*hi*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a&
          & - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) -&
          & 3.0d0*Ts_a*delta_a*hi*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a -&
          & Ts_a)*log(nu) + 1.5d0*sqrt(pi)*Ts_a*delta_a*hia*nu*nu**(-T_a)*(&
          & -log(nu))**(-1.5d0)*(T_a - Ts_a)**(-1.5d0)*log(nu)*erf(sqrt(-log(&
          & nu))*sqrt(T_a - Ts_a)) - 3.0d0*Ts_a*delta_a*hia*nu*nu**(-T_a)*1.0&
          & /(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 0.375d0*Ts_a*hi**2*nu*nu**&
          & (-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) + 0.75d0*Ts_a*hi*&
          & hia*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)*log(nu) +&
          & 0.375d0*Ts_a*hia**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(T_a - Ts_a)&
          & *log(nu) - 6.0d0*delta_a**2*nu*nu**(-Ts_a)*1.0/(-log(nu))*1.0/(&
          & T_a - Ts_a) + 6.0d0*delta_a**2*nu*nu**(-T_a)*1.0/(-log(nu))*1.0/(&
          & T_a - Ts_a)
       ELSEIF (nu==1.0) THEN
          dphia_dhi1 = 0.375d0*hi**2 + 0.75d0*hi*hia + 0.375d0*hia**2
          dphia_dhi = 0.375d0*hi**2 + 0.75d0*hi*hia + 0.375d0*hia**2
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FdPhi_A

  SUBROUTINE  fPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,delta_b2,delta_b3,h_b2,h_b3,Delta_T_b,&
       &phi_b,nu,Rheology,ERROR_CODE)


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bel,Bgrav,h_b,delta_b,T_b,Ts_b
    DOUBLE PRECISION, INTENT(IN) :: delta_b2,delta_b3,h_b2,h_b3,Delta_T_b
    
    DOUBLE PRECISION, INTENT(INOUT) :: phi_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0
    
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
          phi_b = -12.0d0*sqrt(pi)*T_b*delta_b**3*nu*nu**(-T_b)*(-log(nu))**&
          & (-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(&
          & T_b - Ts_b)) + 8.0d0*T_b*delta_b**3*nu*nu**(-T_b)*1.0/(-log(nu))*&
          & 1.0/(T_b - Ts_b)*log(nu) + 12.0d0*sqrt(pi)*T_b*delta_b**2*h_b*nu*&
          & nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*&
          & erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 12.0d0*T_b*delta_b**2*h_b*&
          & nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) - 3.0d0*&
          & sqrt(pi)*T_b*delta_b*h_b**2*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(&
          & T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b&
          & )) + 6.0d0*T_b*delta_b*h_b**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(&
          & T_b - Ts_b)*log(nu) - 1.0d0*T_b*h_b**3*nu*nu**(-T_b)*1.0/(-log(nu&
          & ))*1.0/(T_b - Ts_b)*log(nu) + 12.0d0*sqrt(pi)*Ts_b*delta_b**3*nu*&
          & nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*&
          & erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 8.0d0*Ts_b*delta_b**3*nu*&
          & nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) - 12.0d0*sqrt(&
          & pi)*Ts_b*delta_b**2*h_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b -&
          & Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) +&
          & 12.0d0*Ts_b*delta_b**2*h_b*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b&
          & - Ts_b)*log(nu) + 3.0d0*sqrt(pi)*Ts_b*delta_b*h_b**2*nu*nu**(-T_b&
          & )*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(&
          & -log(nu))*sqrt(T_b - Ts_b)) - 6.0d0*Ts_b*delta_b*h_b**2*nu*nu**(&
          & -T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 1.0d0*Ts_b*h_b**3&
          & *nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 12.0d0*&
          & delta_b**3*nu*nu**(-Ts_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b) + 6.0d0&
          & *sqrt(pi)*delta_b**3*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b -&
          & Ts_b)**(-1.5d0)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 24.0d0*&
          & delta_b**3*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b) - 12.0d0&
          & *delta_b**2*h_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b) +&
          & 12.0d0*delta_b**2*h_b*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b -&
          & Ts_b)
       ELSEIF (nu==1.0) THEN
          phi_b = 1.0d0*h_b**3
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_B
  
  SUBROUTINE  fDPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,hib,hi,Delta_T_b,delta_b2,&
       &dphib_dhi,dphib_dhi1,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bel,Bgrav,h_b,T_b,delta_b,Delta_T_b
    DOUBLE PRECISION, INTENT(IN) :: Ts_b,hi,hib,delta_b2

    DOUBLE PRECISION, INTENT(INOUT) :: dphib_dhi,dphib_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0

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
          dphib_dhi = 6.0d0*sqrt(pi)*T_b*delta_b**2*nu*nu**(-T_b)*(-log(nu))&
          & **(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt&
          & (T_b - Ts_b)) - 6.0d0*T_b*delta_b**2*nu*nu**(-T_b)*1.0/(-log(nu))&
          & *1.0/(T_b - Ts_b)*log(nu) - 1.5d0*sqrt(pi)*T_b*delta_b*hi*nu*nu**&
          & (-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(&
          & sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*T_b*delta_b*hi*nu*nu**(&
          & -T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) - 1.5d0*sqrt(pi)*&
          & T_b*delta_b*hib*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**&
          & (-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*T_b&
          & *delta_b*hib*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu&
          & ) - 0.375d0*T_b*hi**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b -&
          & Ts_b)*log(nu) - 0.75d0*T_b*hi*hib*nu*nu**(-T_b)*1.0/(-log(nu))*&
          & 1.0/(T_b - Ts_b)*log(nu) - 0.375d0*T_b*hib**2*nu*nu**(-T_b)*1.0/(&
          & -log(nu))*1.0/(T_b - Ts_b)*log(nu) - 6.0d0*sqrt(pi)*Ts_b*delta_b&
          & **2*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log&
          & (nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 6.0d0*Ts_b*delta_b**2&
          & *nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 1.5d0*&
          & sqrt(pi)*Ts_b*delta_b*hi*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b&
          & - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) -&
          & 3.0d0*Ts_b*delta_b*hi*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b -&
          & Ts_b)*log(nu) + 1.5d0*sqrt(pi)*Ts_b*delta_b*hib*nu*nu**(-T_b)*(&
          & -log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(&
          & nu))*sqrt(T_b - Ts_b)) - 3.0d0*Ts_b*delta_b*hib*nu*nu**(-T_b)*1.0&
          & /(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 0.375d0*Ts_b*hi**2*nu*nu**&
          & (-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 0.75d0*Ts_b*hi*&
          & hib*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) +&
          & 0.375d0*Ts_b*hib**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)&
          & *log(nu) - 6.0d0*delta_b**2*nu*nu**(-Ts_b)*1.0/(-log(nu))*1.0/(&
          & T_b - Ts_b) + 6.0d0*delta_b**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(&
          & T_b - Ts_b)
     
          dphib_dhi1 = 6.0d0*sqrt(pi)*T_b*delta_b**2*nu*nu**(-T_b)*(-log(nu&
          & ))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*&
          & sqrt(T_b - Ts_b)) - 6.0d0*T_b*delta_b**2*nu*nu**(-T_b)*1.0/(-log(&
          & nu))*1.0/(T_b - Ts_b)*log(nu) - 1.5d0*sqrt(pi)*T_b*delta_b*hi*nu*&
          & nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*&
          & erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*T_b*delta_b*hi*nu*nu&
          & **(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) - 1.5d0*sqrt(pi)&
          & *T_b*delta_b*hib*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)&
          & **(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*&
          & T_b*delta_b*hib*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log&
          & (nu) - 0.375d0*T_b*hi**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b -&
          & Ts_b)*log(nu) - 0.75d0*T_b*hi*hib*nu*nu**(-T_b)*1.0/(-log(nu))*&
          & 1.0/(T_b - Ts_b)*log(nu) - 0.375d0*T_b*hib**2*nu*nu**(-T_b)*1.0/(&
          & -log(nu))*1.0/(T_b - Ts_b)*log(nu) - 6.0d0*sqrt(pi)*Ts_b*delta_b&
          & **2*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log&
          & (nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 6.0d0*Ts_b*delta_b**2&
          & *nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 1.5d0*&
          & sqrt(pi)*Ts_b*delta_b*hi*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b&
          & - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) -&
          & 3.0d0*Ts_b*delta_b*hi*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b -&
          & Ts_b)*log(nu) + 1.5d0*sqrt(pi)*Ts_b*delta_b*hib*nu*nu**(-T_b)*(&
          & -log(nu))**(-1.5d0)*(T_b - Ts_b)**(-1.5d0)*log(nu)*erf(sqrt(-log(&
          & nu))*sqrt(T_b - Ts_b)) - 3.0d0*Ts_b*delta_b*hib*nu*nu**(-T_b)*1.0&
          & /(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 0.375d0*Ts_b*hi**2*nu*nu**&
          & (-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) + 0.75d0*Ts_b*hi*&
          & hib*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)*log(nu) +&
          & 0.375d0*Ts_b*hib**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(T_b - Ts_b)&
          & *log(nu) - 6.0d0*delta_b**2*nu*nu**(-Ts_b)*1.0/(-log(nu))*1.0/(&
          & T_b - Ts_b) + 6.0d0*delta_b**2*nu*nu**(-T_b)*1.0/(-log(nu))*1.0/(&
          & T_b - Ts_b)
       ELSEIF (nu == 1.0) THEN
          dphib_dhi = 0.375d0*hi**2 + 0.75d0*hi*hib + 0.375d0*hib**2
          dphib_dhi1 = 0.375d0*hi**2 + 0.75d0*hi*hib + 0.375d0*hib**2
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FDPhi_B


END MODULE MOBILITY_THICKNESS_SKIN_RHEOLOGY
