MODULE MOBILITY_THERMAL_SKIN_RHEOLOGY

CONTAINS

  SUBROUTINE fomega_a(Ai,h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3&
       &,Ts_a,Ts_a2,Ts_a3,Ds_a,eta_a,omega_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0

    DOUBLE PRECISION, INTENT(IN) :: Ai,eta_a
    DOUBLE PRECISION, INTENT(IN) :: h_a,h_a2,h_a3
    DOUBLE PRECISION, INTENT(IN) :: delta_a,delta_a2,delta_a3
    DOUBLE PRECISION, INTENT(IN) :: T_a,T_a2,T_a3
    DOUBLE PRECISION, INTENT(IN) :: Ts_a,Ts_a2,Ts_a3
    DOUBLE PRECISION, INTENT(IN) :: Ds_a

    DOUBLE PRECISION, INTENT(INOUT) :: omega_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    IF (Rheology == 0) THEN
       omega_a = (7.0d0/5.0d0)*T_a*delta_a**2*eta_a*nu - 7.0d0/5.0d0*T_a*&
            & delta_a**2*eta_a - 3.0d0/2.0d0*T_a*delta_a*eta_a*h_a*nu + (3.0d0/&
            & 2.0d0)*T_a*delta_a*eta_a*h_a + (3.0d0/5.0d0)*Ts_a*delta_a**2*&
            & eta_a*nu - 3.0d0/5.0d0*Ts_a*delta_a**2*eta_a - 3.0d0/2.0d0*Ts_a*&
            & delta_a*eta_a*h_a*nu + (3.0d0/2.0d0)*Ts_a*delta_a*eta_a*h_a - 2*&
            & delta_a**2*eta_a*nu + 3*delta_a*eta_a*h_a*nu

       omega_a = (7.0d0/5.0d0)*T_a*delta_a2*eta_a*nu - 7.0d0/5.0d0*T_a*&
            & delta_a2*eta_a - 3.0d0/2.0d0*T_a*delta_a*eta_a*h_a*nu + (3.0d0/&
            & 2.0d0)*T_a*delta_a*eta_a*h_a + (3.0d0/5.0d0)*Ts_a*delta_a2*&
            & eta_a*nu - 3.0d0/5.0d0*Ts_a*delta_a2*eta_a - 3.0d0/2.0d0*Ts_a*&
            & delta_a*eta_a*h_a*nu + (3.0d0/2.0d0)*Ts_a*delta_a*eta_a*h_a - 2*&
            & delta_a2*eta_a*nu + 3*delta_a*eta_a*h_a*nu

       ! omega_a = (eta_a*delta_a)/10.d0*(nu*(-20.d0*delta_a+30.d0*h_a)+&
       !      &(1.d0-nu)*(6.d0*Ds_a*delta_a-15.d0*Ds_a*h_a-20.d0*T_a*delta_a+30.d0*T_a*h_a))

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          omega_a = -6.0d0*sqrt(pi)*T_a**2*delta_a**2*eta_a*nu*nu**(-T_a)*(&
               & -log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(&
               & nu))*sqrt(T_a - Ts_a)) + 3.0d0*sqrt(pi)*T_a**2*delta_a*eta_a*h_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)&
               & *erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 12.0d0*sqrt(pi)*T_a*Ts_a*&
               & delta_a**2*eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)&
               & **(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 6.0d0*&
               & sqrt(pi)*T_a*Ts_a*delta_a*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(&
               & -1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(&
               & T_a - Ts_a)) + 6.0d0*T_a*delta_a**2*eta_a*nu*nu**(-Ts_a)*1.0/(&
               & -log(nu))*(T_a - Ts_a)**(-2.0d0) + 3.0d0*sqrt(pi)*T_a*delta_a**2*&
               & eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*&
               & erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 12.0d0*T_a*delta_a**2*&
               & eta_a*nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) - 3.0d0&
               & *T_a*delta_a*eta_a*h_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)&
               & **(-2.0d0) + 3.0d0*T_a*delta_a*eta_a*h_a*nu*nu**(-T_a)*1.0/(-log(&
               & nu))*(T_a - Ts_a)**(-2.0d0) - 6.0d0*sqrt(pi)*Ts_a**2*delta_a**2*&
               & eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*&
               & log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*sqrt(pi)*&
               & Ts_a**2*delta_a*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a&
               & - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) -&
               & 6.0d0*Ts_a*delta_a**2*eta_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a -&
               & Ts_a)**(-2.0d0) - 3.0d0*sqrt(pi)*Ts_a*delta_a**2*eta_a*nu*nu**(&
               & -T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*erf(sqrt(-log(&
               & nu))*sqrt(T_a - Ts_a)) + 12.0d0*Ts_a*delta_a**2*eta_a*nu*nu**(&
               & -T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) + 3.0d0*Ts_a*delta_a*&
               & eta_a*h_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) -&
               & 3.0d0*Ts_a*delta_a*eta_a*h_a*nu*nu**(-T_a)*1.0/(-log(nu))*(T_a -&
               & Ts_a)**(-2.0d0) - 6.0d0*delta_a**2*eta_a*nu*nu**(-Ts_a)*1.0/(T_a&
               & - Ts_a)/log(nu) - 6.0d0*sqrt(pi)*delta_a**2*eta_a*nu*nu**(-T_a)*&
               & sqrt(-log(nu))*(T_a - Ts_a)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a&
               & - Ts_a))/log(nu) + 6.0d0*delta_a**2*eta_a*nu*nu**(-T_a)*1.0/(T_a&
               & - Ts_a)/log(nu) + 3.0d0*sqrt(pi)*delta_a*eta_a*h_a*nu*nu**(-T_a)*&
               & sqrt(-log(nu))*(T_a - Ts_a)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a&
               & - Ts_a))/log(nu)

       ELSEIF (nu==1.0) THEN
          omega_a = -2.0d0*delta_a**2*eta_a + 3.0d0*delta_a*eta_a*h_a
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fomega_A

  SUBROUTINE fsigma_a(Ai,h_a,h_a2,h_a3,delta_a,delta_a2,delta_a3,T_a,T_a2,T_a3&
       &,Ts_a,Ts_a2,Ts_a3,Ds_a,eta_a,sigma_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0

    DOUBLE PRECISION, INTENT(IN) :: Ai,eta_a
    DOUBLE PRECISION, INTENT(IN) :: h_a,h_a2,h_a3
    DOUBLE PRECISION, INTENT(IN) :: delta_a,delta_a2,delta_a3
    DOUBLE PRECISION, INTENT(IN) :: T_a,T_a2,T_a3
    DOUBLE PRECISION, INTENT(IN) :: Ts_a,Ts_a2,Ts_a3
    DOUBLE PRECISION, INTENT(IN) :: Ds_a

    DOUBLE PRECISION, INTENT(INOUT) :: sigma_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE

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

          sigma_a = -38.0d0/105.0d0*T_a2*delta_a3*eta_a*nu + (38.0d0/&
               & 105.0d0)*T_a2*delta_a3*eta_a + (1.0d0/3.0d0)*T_a2*delta_a2&
               & *eta_a*h_a*nu - 1.0d0/3.0d0*T_a2*delta_a2*eta_a*h_a + (9.0d0&
               & /35.0d0)*T_a*Ts_a*delta_a3*eta_a*nu - 9.0d0/35.0d0*T_a*Ts_a*&
               & delta_a3*eta_a - 1.0d0/6.0d0*T_a*Ts_a*delta_a2*eta_a*h_a*nu +&
               & (1.0d0/6.0d0)*T_a*Ts_a*delta_a2*eta_a*h_a + (7.0d0/15.0d0)*T_a*&
               & delta_a3*eta_a*nu - 1.0d0/2.0d0*T_a*delta_a2*eta_a*h_a*nu + (&
               & 11.0d0/105.0d0)*Ts_a2*delta_a3*eta_a*nu - 11.0d0/105.0d0*Ts_a2&
               & *delta_a3*eta_a - 1.0d0/6.0d0*Ts_a2*delta_a2*eta_a*h_a*&
               & nu + (1.0d0/6.0d0)*Ts_a2*delta_a2*eta_a*h_a - 7.0d0/15.0d0*&
               & Ts_a*delta_a3*eta_a*nu + (1.0d0/2.0d0)*Ts_a*delta_a2*eta_a*&
               & h_a*nu

       ! sigma_a = (-1.d0/210.d0)*Ds_a*delta_a2*eta_a*(nu*(-98.d0*delta_a+105.d0*h_a)+&
       !      &(1-nu)*(22.d0*Ds_a*delta_a-35.d0*Ds_a*h_a-98.d0*T_a*delta_a+105.d0*T_a*h_a))
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          sigma_a = -6.0d0*sqrt(pi)*T_a**3*delta_a**3*eta_a*nu*nu**(-T_a)*(&
               & -log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)**2*erf(sqrt(&
               & -log(nu))*sqrt(T_a - Ts_a)) - 4.0d0*sqrt(pi)*T_a**3*delta_a**3*&
               & eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*&
               & log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 3.0d0*sqrt(pi)*T_a&
               & **3*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a&
               & - Ts_a)**(-2.5d0)*log(nu)**2*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a&
               & )) + 2.0d0*sqrt(pi)*T_a**3*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*(&
               & -log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(&
               & nu))*sqrt(T_a - Ts_a)) + 12.0d0*sqrt(pi)*T_a**2*Ts_a*delta_a**3*&
               & eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*&
               & log(nu)**2*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 6.0d0*sqrt(pi)*&
               & T_a**2*Ts_a*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(&
               & T_a - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a&
               & )) - 6.0d0*sqrt(pi)*T_a**2*Ts_a*delta_a**2*eta_a*h_a*nu*nu**(-T_a&
               & )*(-log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)**2*erf(sqrt&
               & (-log(nu))*sqrt(T_a - Ts_a)) - 3.0d0*sqrt(pi)*T_a**2*Ts_a*delta_a&
               & **2*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(&
               & -2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 6.0d0*T_a&
               & **2*delta_a**3*eta_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*(T_a -&
               & Ts_a)**(-2.0d0)*log(nu) + 4.0d0*T_a**2*delta_a**3*eta_a*nu*nu**(&
               & -Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) + 3.0d0*sqrt(pi)*T_a&
               & **2*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a -&
               & Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) -&
               & 12.0d0*T_a**2*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.0d0)&
               & *(T_a - Ts_a)**(-2.0d0)*log(nu) + 2.0d0*sqrt(pi)*T_a**2*delta_a**&
               & 3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)&
               & *erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 8.0d0*T_a**2*delta_a**3*&
               & eta_a*nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) - 2.0d0&
               & *T_a**2*delta_a**2*eta_a*h_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*&
               & (T_a - Ts_a)**(-2.0d0)*log(nu) - 2.0d0*T_a**2*delta_a**2*eta_a*&
               & h_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) + 3.0d0*&
               & T_a**2*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(-2.0d0)*(&
               & T_a - Ts_a)**(-2.0d0)*log(nu) + 2.0d0*T_a**2*delta_a**2*eta_a*h_a&
               & *nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) - 6.0d0*sqrt&
               & (pi)*T_a*Ts_a**2*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(&
               & -2.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)**2*erf(sqrt(-log(nu))*sqrt&
               & (T_a - Ts_a)) + 3.0d0*sqrt(pi)*T_a*Ts_a**2*delta_a**2*eta_a*h_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)&
               & **2*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 6.0d0*T_a*Ts_a*delta_a&
               & **3*eta_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(&
               & -2.0d0)*log(nu) - 2.0d0*T_a*Ts_a*delta_a**3*eta_a*nu*nu**(-Ts_a)*&
               & 1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) - 3.0d0*sqrt(pi)*T_a*Ts_a*&
               & delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a - Ts_a)&
               & **(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 12.0d0*&
               & T_a*Ts_a*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))**(-2.0d0)*(T_a&
               & - Ts_a)**(-2.0d0)*log(nu) - 1.0d0*sqrt(pi)*T_a*Ts_a*delta_a**3*&
               & eta_a*nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*&
               & erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) + 4.0d0*T_a*Ts_a*delta_a**3*&
               & eta_a*nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) + 1.0d0&
               & *T_a*Ts_a*delta_a**2*eta_a*h_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0&
               & )*(T_a - Ts_a)**(-2.0d0)*log(nu) + 1.0d0*T_a*Ts_a*delta_a**2*&
               & eta_a*h_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) -&
               & 3.0d0*T_a*Ts_a*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*(-log(nu))**(&
               & -2.0d0)*(T_a - Ts_a)**(-2.0d0)*log(nu) - 1.0d0*T_a*Ts_a*delta_a**&
               & 2*eta_a*h_a*nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) -&
               & 1.0d0*T_a*delta_a**3*eta_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*(&
               & T_a - Ts_a)**(-2.0d0) + 2.0d0*T_a*delta_a**3*eta_a*nu*nu**(-Ts_a)&
               & *1.0/(T_a - Ts_a)/log(nu) + 1.5d0*sqrt(pi)*T_a*delta_a**3*eta_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*erf(&
               & sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 2.0d0*T_a*delta_a**3*eta_a*nu*&
               & nu**(-T_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(-2.0d0) + 2.0d0*&
               & sqrt(pi)*T_a*delta_a**3*eta_a*nu*nu**(-T_a)*sqrt(-log(nu))*(T_a -&
               & Ts_a)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a))/log(nu) -&
               & 2.0d0*T_a*delta_a**3*eta_a*nu*nu**(-T_a)*1.0/(T_a - Ts_a)/log(nu&
               & ) - 1.0d0*T_a*delta_a**2*eta_a*h_a*nu*nu**(-Ts_a)*(-log(nu))**(&
               & -2.0d0)*(T_a - Ts_a)**(-2.0d0) + 1.0d0*T_a*delta_a**2*eta_a*h_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(-2.0d0) - 1.0d0&
               & *sqrt(pi)*T_a*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*sqrt(-log(nu))*(&
               & T_a - Ts_a)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a))/log(nu&
               & ) - 2.0d0*sqrt(pi)*Ts_a**3*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(&
               & nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*&
               & sqrt(T_a - Ts_a)) + 1.0d0*sqrt(pi)*Ts_a**3*delta_a**2*eta_a*h_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*log(nu)&
               & *erf(sqrt(-log(nu))*sqrt(T_a - Ts_a)) - 2.0d0*Ts_a**2*delta_a**3*&
               & eta_a*nu*nu**(-Ts_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) -&
               & 1.0d0*sqrt(pi)*Ts_a**2*delta_a**3*eta_a*nu*nu**(-T_a)*(-log(nu))&
               & **(-1.5d0)*(T_a - Ts_a)**(-2.5d0)*erf(sqrt(-log(nu))*sqrt(T_a -&
               & Ts_a)) + 4.0d0*Ts_a**2*delta_a**3*eta_a*nu*nu**(-T_a)*1.0/(-log(&
               & nu))*(T_a - Ts_a)**(-2.0d0) + 1.0d0*Ts_a**2*delta_a**2*eta_a*h_a*&
               & nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(-2.0d0)*log(nu&
               & ) + 1.0d0*Ts_a**2*delta_a**2*eta_a*h_a*nu*nu**(-Ts_a)*1.0/(-log(&
               & nu))*(T_a - Ts_a)**(-2.0d0) - 1.0d0*Ts_a**2*delta_a**2*eta_a*h_a*&
               & nu*nu**(-T_a)*1.0/(-log(nu))*(T_a - Ts_a)**(-2.0d0) + 1.0d0*Ts_a*&
               & delta_a**3*eta_a*nu*nu**(-Ts_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)&
               & **(-2.0d0) - 2.0d0*Ts_a*delta_a**3*eta_a*nu*nu**(-Ts_a)*1.0/(T_a&
               & - Ts_a)/log(nu) - 1.5d0*sqrt(pi)*Ts_a*delta_a**3*eta_a*nu*nu**(&
               & -T_a)*(-log(nu))**(-2.5d0)*(T_a - Ts_a)**(-2.5d0)*erf(sqrt(-log(&
               & nu))*sqrt(T_a - Ts_a)) + 2.0d0*Ts_a*delta_a**3*eta_a*nu*nu**(-T_a&
               & )*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(-2.0d0) - 2.0d0*sqrt(pi)*&
               & Ts_a*delta_a**3*eta_a*nu*nu**(-T_a)*sqrt(-log(nu))*(T_a - Ts_a)**&
               & (-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a))/log(nu) + 2.0d0*&
               & Ts_a*delta_a**3*eta_a*nu*nu**(-T_a)*1.0/(T_a - Ts_a)/log(nu) +&
               & 1.0d0*Ts_a*delta_a**2*eta_a*h_a*nu*nu**(-Ts_a)*(-log(nu))**(&
               & -2.0d0)*(T_a - Ts_a)**(-2.0d0) - 1.0d0*Ts_a*delta_a**2*eta_a*h_a*&
               & nu*nu**(-T_a)*(-log(nu))**(-2.0d0)*(T_a - Ts_a)**(-2.0d0) + 1.0d0&
               & *sqrt(pi)*Ts_a*delta_a**2*eta_a*h_a*nu*nu**(-T_a)*sqrt(-log(nu))*&
               & (T_a - Ts_a)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_a - Ts_a))/log(&
               & nu)
       ELSEIF (nu==1.0) THEN
          sigma_a = 0.466666666666665d0*T_a*delta_a**3*eta_a - 0.5d0*T_a*&
               &delta_a**2*eta_a*h_a - 0.466666666666667d0*Ts_a*delta_a**3*eta_a&
               & + 0.5d0*Ts_a*delta_a**2*eta_a*h_a
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_A


  SUBROUTINE fomega_b(Bi,h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3&
       &,Ts_b,Ts_b2,Ts_b3,Ds_b,eta_b,omega_b,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0

    DOUBLE PRECISION, INTENT(IN) :: Bi,eta_b
    DOUBLE PRECISION, INTENT(IN) :: h_b,h_b2,h_b3
    DOUBLE PRECISION, INTENT(IN) :: delta_b,delta_b2,delta_b3
    DOUBLE PRECISION, INTENT(IN) :: T_b,T_b2,T_b3
    DOUBLE PRECISION, INTENT(IN) :: Ts_b,Ts_b2,Ts_b3
    DOUBLE PRECISION, INTENT(IN) :: Ds_b
    
    DOUBLE PRECISION, INTENT(INOUT) :: omega_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    IF (Rheology == 0) THEN
          omega_b = (7.0d0/5.0d0)*T_b*delta_b**2*eta_b*nu - 7.0d0/5.0d0*T_b*&
               & delta_b**2*eta_b - 3.0d0/2.0d0*T_b*delta_b*eta_b*h_b*nu + (3.0d0/&
               & 2.0d0)*T_b*delta_b*eta_b*h_b + (3.0d0/5.0d0)*Ts_b*delta_b**2*&
               & eta_b*nu - 3.0d0/5.0d0*Ts_b*delta_b**2*eta_b - 3.0d0/2.0d0*Ts_b*&
               & delta_b*eta_b*h_b*nu + (3.0d0/2.0d0)*Ts_b*delta_b*eta_b*h_b - 2*&
               & delta_b**2*eta_b*nu + 3*delta_b*eta_b*h_b*nu

          omega_b = (7.0d0/5.0d0)*T_b*delta_b2*eta_b*nu - 7.0d0/5.0d0*T_b*&
               & delta_b2*eta_b - 3.0d0/2.0d0*T_b*delta_b*eta_b*h_b*nu + (3.0d0/&
               & 2.0d0)*T_b*delta_b*eta_b*h_b + (3.0d0/5.0d0)*Ts_b*delta_b2*&
               & eta_b*nu - 3.0d0/5.0d0*Ts_b*delta_b2*eta_b - 3.0d0/2.0d0*Ts_b*&
               & delta_b*eta_b*h_b*nu + (3.0d0/2.0d0)*Ts_b*delta_b*eta_b*h_b - 2*&
               & delta_b2*eta_b*nu + 3*delta_b*eta_b*h_b*nu
       
       ! omega_b = (eta_b*delta_b)/10.d0*(nu*(-20.d0*delta_b+30.d0*h_b)+&
       !      &(1.d0-nu)*(6.d0*Ds_b*delta_b-15.d0*Ds_b*h_b-20.d0*T_b*delta_b+30.d0*T_b*h_b))

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          omega_b = -6.0d0*sqrt(pi)*T_b**2*delta_b**2*eta_b*nu*nu**(-T_b)*(&
               & -log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(&
               & nu))*sqrt(T_b - Ts_b)) + 3.0d0*sqrt(pi)*T_b**2*delta_b*eta_b*h_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)&
               & *erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 12.0d0*sqrt(pi)*T_b*Ts_b*&
               & delta_b**2*eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)&
               & **(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 6.0d0*&
               & sqrt(pi)*T_b*Ts_b*delta_b*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(&
               & -1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(&
               & T_b - Ts_b)) + 6.0d0*T_b*delta_b**2*eta_b*nu*nu**(-Ts_b)*1.0/(&
               & -log(nu))*(T_b - Ts_b)**(-2.0d0) + 3.0d0*sqrt(pi)*T_b*delta_b**2*&
               & eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*&
               & erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 12.0d0*T_b*delta_b**2*&
               & eta_b*nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) - 3.0d0&
               & *T_b*delta_b*eta_b*h_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)&
               & **(-2.0d0) + 3.0d0*T_b*delta_b*eta_b*h_b*nu*nu**(-T_b)*1.0/(-log(&
               & nu))*(T_b - Ts_b)**(-2.0d0) - 6.0d0*sqrt(pi)*Ts_b**2*delta_b**2*&
               & eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*&
               & log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*sqrt(pi)*&
               & Ts_b**2*delta_b*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b&
               & - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) -&
               & 6.0d0*Ts_b*delta_b**2*eta_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b -&
               & Ts_b)**(-2.0d0) - 3.0d0*sqrt(pi)*Ts_b*delta_b**2*eta_b*nu*nu**(&
               & -T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*erf(sqrt(-log(&
               & nu))*sqrt(T_b - Ts_b)) + 12.0d0*Ts_b*delta_b**2*eta_b*nu*nu**(&
               & -T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) + 3.0d0*Ts_b*delta_b*&
               & eta_b*h_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) -&
               & 3.0d0*Ts_b*delta_b*eta_b*h_b*nu*nu**(-T_b)*1.0/(-log(nu))*(T_b -&
               & Ts_b)**(-2.0d0) - 6.0d0*delta_b**2*eta_b*nu*nu**(-Ts_b)*1.0/(T_b&
               & - Ts_b)/log(nu) - 6.0d0*sqrt(pi)*delta_b**2*eta_b*nu*nu**(-T_b)*&
               & sqrt(-log(nu))*(T_b - Ts_b)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b&
               & - Ts_b))/log(nu) + 6.0d0*delta_b**2*eta_b*nu*nu**(-T_b)*1.0/(T_b&
               & - Ts_b)/log(nu) + 3.0d0*sqrt(pi)*delta_b*eta_b*h_b*nu*nu**(-T_b)*&
               & sqrt(-log(nu))*(T_b - Ts_b)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b&
               & - Ts_b))/log(nu)
       ELSEIF (nu==1.0) THEN
          omega_b = -2.0d0*delta_b**2*eta_b + 3.0d0*delta_b*eta_b*h_b
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FOmega_B


  SUBROUTINE fsigma_b(Bi,h_b,h_b2,h_b3,delta_b,delta_b2,delta_b3,T_b,T_b2,T_b3&
       &,Ts_b,Ts_b2,Ts_b3,Ds_b,eta_b,sigma_b,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979d0

    DOUBLE PRECISION, INTENT(IN) :: Bi,eta_b
    DOUBLE PRECISION, INTENT(IN) :: h_b,h_b2,h_b3
    DOUBLE PRECISION, INTENT(IN) :: delta_b,delta_b2,delta_b3
    DOUBLE PRECISION, INTENT(IN) :: T_b,T_b2,T_b3
    DOUBLE PRECISION, INTENT(IN) :: Ts_b,Ts_b2,Ts_b3
    DOUBLE PRECISION, INTENT(IN) :: Ds_b

    DOUBLE PRECISION, INTENT(INOUT) :: sigma_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE

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

          sigma_b = -38.0d0/105.0d0*T_b2*delta_b3*eta_b*nu + (38.0d0/&
               & 105.0d0)*T_b2*delta_b3*eta_b + (1.0d0/3.0d0)*T_b2*delta_b2&
               & *eta_b*h_b*nu - 1.0d0/3.0d0*T_b2*delta_b2*eta_b*h_b + (9.0d0&
               & /35.0d0)*T_b*Ts_b*delta_b3*eta_b*nu - 9.0d0/35.0d0*T_b*Ts_b*&
               & delta_b3*eta_b - 1.0d0/6.0d0*T_b*Ts_b*delta_b2*eta_b*h_b*nu +&
               & (1.0d0/6.0d0)*T_b*Ts_b*delta_b2*eta_b*h_b + (7.0d0/15.0d0)*T_b*&
               & delta_b3*eta_b*nu - 1.0d0/2.0d0*T_b*delta_b2*eta_b*h_b*nu + (&
               & 11.0d0/105.0d0)*Ts_b2*delta_b3*eta_b*nu - 11.0d0/105.0d0*Ts_b2&
               & *delta_b3*eta_b - 1.0d0/6.0d0*Ts_b2*delta_b2*eta_b*h_b*&
               & nu + (1.0d0/6.0d0)*Ts_b2*delta_b2*eta_b*h_b - 7.0d0/15.0d0*&
               & Ts_b*delta_b3*eta_b*nu + (1.0d0/2.0d0)*Ts_b*delta_b2*eta_b*&
               & h_b*nu
          
       ! sigma_b = (-1.d0/210.d0)*Ds_b*delta_b2*eta_b*(nu*(-98.d0*delta_b+105.d0*h_b)+&
       !      &(1-nu)*(22.d0*Ds_b*delta_b-35.d0*Ds_b*h_b-98.d0*T_b*delta_b+105.d0*T_b*h_b))
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       nutest:IF (nu /= 1.0) THEN
          sigma_b = -6.0d0*sqrt(pi)*T_b**3*delta_b**3*eta_b*nu*nu**(-T_b)*(&
               & -log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)**2*erf(sqrt(&
               & -log(nu))*sqrt(T_b - Ts_b)) - 4.0d0*sqrt(pi)*T_b**3*delta_b**3*&
               & eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*&
               & log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 3.0d0*sqrt(pi)*T_b&
               & **3*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b&
               & - Ts_b)**(-2.5d0)*log(nu)**2*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b&
               & )) + 2.0d0*sqrt(pi)*T_b**3*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*(&
               & -log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(&
               & nu))*sqrt(T_b - Ts_b)) + 12.0d0*sqrt(pi)*T_b**2*Ts_b*delta_b**3*&
               & eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*&
               & log(nu)**2*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 6.0d0*sqrt(pi)*&
               & T_b**2*Ts_b*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(&
               & T_b - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b&
               & )) - 6.0d0*sqrt(pi)*T_b**2*Ts_b*delta_b**2*eta_b*h_b*nu*nu**(-T_b&
               & )*(-log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)**2*erf(sqrt&
               & (-log(nu))*sqrt(T_b - Ts_b)) - 3.0d0*sqrt(pi)*T_b**2*Ts_b*delta_b&
               & **2*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(&
               & -2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 6.0d0*T_b&
               & **2*delta_b**3*eta_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*(T_b -&
               & Ts_b)**(-2.0d0)*log(nu) + 4.0d0*T_b**2*delta_b**3*eta_b*nu*nu**(&
               & -Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) + 3.0d0*sqrt(pi)*T_b&
               & **2*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b -&
               & Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) -&
               & 12.0d0*T_b**2*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.0d0)&
               & *(T_b - Ts_b)**(-2.0d0)*log(nu) + 2.0d0*sqrt(pi)*T_b**2*delta_b**&
               & 3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)&
               & *erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 8.0d0*T_b**2*delta_b**3*&
               & eta_b*nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) - 2.0d0&
               & *T_b**2*delta_b**2*eta_b*h_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*&
               & (T_b - Ts_b)**(-2.0d0)*log(nu) - 2.0d0*T_b**2*delta_b**2*eta_b*&
               & h_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) + 3.0d0*&
               & T_b**2*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(-2.0d0)*(&
               & T_b - Ts_b)**(-2.0d0)*log(nu) + 2.0d0*T_b**2*delta_b**2*eta_b*h_b&
               & *nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) - 6.0d0*sqrt&
               & (pi)*T_b*Ts_b**2*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(&
               & -2.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)**2*erf(sqrt(-log(nu))*sqrt&
               & (T_b - Ts_b)) + 3.0d0*sqrt(pi)*T_b*Ts_b**2*delta_b**2*eta_b*h_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)&
               & **2*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 6.0d0*T_b*Ts_b*delta_b&
               & **3*eta_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(&
               & -2.0d0)*log(nu) - 2.0d0*T_b*Ts_b*delta_b**3*eta_b*nu*nu**(-Ts_b)*&
               & 1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) - 3.0d0*sqrt(pi)*T_b*Ts_b*&
               & delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b - Ts_b)&
               & **(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 12.0d0*&
               & T_b*Ts_b*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))**(-2.0d0)*(T_b&
               & - Ts_b)**(-2.0d0)*log(nu) - 1.0d0*sqrt(pi)*T_b*Ts_b*delta_b**3*&
               & eta_b*nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*&
               & erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) + 4.0d0*T_b*Ts_b*delta_b**3*&
               & eta_b*nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) + 1.0d0&
               & *T_b*Ts_b*delta_b**2*eta_b*h_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0&
               & )*(T_b - Ts_b)**(-2.0d0)*log(nu) + 1.0d0*T_b*Ts_b*delta_b**2*&
               & eta_b*h_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) -&
               & 3.0d0*T_b*Ts_b*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*(-log(nu))**(&
               & -2.0d0)*(T_b - Ts_b)**(-2.0d0)*log(nu) - 1.0d0*T_b*Ts_b*delta_b**&
               & 2*eta_b*h_b*nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) -&
               & 1.0d0*T_b*delta_b**3*eta_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*(&
               & T_b - Ts_b)**(-2.0d0) + 2.0d0*T_b*delta_b**3*eta_b*nu*nu**(-Ts_b)&
               & *1.0/(T_b - Ts_b)/log(nu) + 1.5d0*sqrt(pi)*T_b*delta_b**3*eta_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*erf(&
               & sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 2.0d0*T_b*delta_b**3*eta_b*nu*&
               & nu**(-T_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(-2.0d0) + 2.0d0*&
               & sqrt(pi)*T_b*delta_b**3*eta_b*nu*nu**(-T_b)*sqrt(-log(nu))*(T_b -&
               & Ts_b)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b))/log(nu) -&
               & 2.0d0*T_b*delta_b**3*eta_b*nu*nu**(-T_b)*1.0/(T_b - Ts_b)/log(nu&
               & ) - 1.0d0*T_b*delta_b**2*eta_b*h_b*nu*nu**(-Ts_b)*(-log(nu))**(&
               & -2.0d0)*(T_b - Ts_b)**(-2.0d0) + 1.0d0*T_b*delta_b**2*eta_b*h_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(-2.0d0) - 1.0d0&
               & *sqrt(pi)*T_b*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*sqrt(-log(nu))*(&
               & T_b - Ts_b)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b))/log(nu&
               & ) - 2.0d0*sqrt(pi)*Ts_b**3*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(&
               & nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)*erf(sqrt(-log(nu))*&
               & sqrt(T_b - Ts_b)) + 1.0d0*sqrt(pi)*Ts_b**3*delta_b**2*eta_b*h_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*log(nu)&
               & *erf(sqrt(-log(nu))*sqrt(T_b - Ts_b)) - 2.0d0*Ts_b**2*delta_b**3*&
               & eta_b*nu*nu**(-Ts_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) -&
               & 1.0d0*sqrt(pi)*Ts_b**2*delta_b**3*eta_b*nu*nu**(-T_b)*(-log(nu))&
               & **(-1.5d0)*(T_b - Ts_b)**(-2.5d0)*erf(sqrt(-log(nu))*sqrt(T_b -&
               & Ts_b)) + 4.0d0*Ts_b**2*delta_b**3*eta_b*nu*nu**(-T_b)*1.0/(-log(&
               & nu))*(T_b - Ts_b)**(-2.0d0) + 1.0d0*Ts_b**2*delta_b**2*eta_b*h_b*&
               & nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(-2.0d0)*log(nu&
               & ) + 1.0d0*Ts_b**2*delta_b**2*eta_b*h_b*nu*nu**(-Ts_b)*1.0/(-log(&
               & nu))*(T_b - Ts_b)**(-2.0d0) - 1.0d0*Ts_b**2*delta_b**2*eta_b*h_b*&
               & nu*nu**(-T_b)*1.0/(-log(nu))*(T_b - Ts_b)**(-2.0d0) + 1.0d0*Ts_b*&
               & delta_b**3*eta_b*nu*nu**(-Ts_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)&
               & **(-2.0d0) - 2.0d0*Ts_b*delta_b**3*eta_b*nu*nu**(-Ts_b)*1.0/(T_b&
               & - Ts_b)/log(nu) - 1.5d0*sqrt(pi)*Ts_b*delta_b**3*eta_b*nu*nu**(&
               & -T_b)*(-log(nu))**(-2.5d0)*(T_b - Ts_b)**(-2.5d0)*erf(sqrt(-log(&
               & nu))*sqrt(T_b - Ts_b)) + 2.0d0*Ts_b*delta_b**3*eta_b*nu*nu**(-T_b&
               & )*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(-2.0d0) - 2.0d0*sqrt(pi)*&
               & Ts_b*delta_b**3*eta_b*nu*nu**(-T_b)*sqrt(-log(nu))*(T_b - Ts_b)**&
               & (-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b))/log(nu) + 2.0d0*&
               & Ts_b*delta_b**3*eta_b*nu*nu**(-T_b)*1.0/(T_b - Ts_b)/log(nu) +&
               & 1.0d0*Ts_b*delta_b**2*eta_b*h_b*nu*nu**(-Ts_b)*(-log(nu))**(&
               & -2.0d0)*(T_b - Ts_b)**(-2.0d0) - 1.0d0*Ts_b*delta_b**2*eta_b*h_b*&
               & nu*nu**(-T_b)*(-log(nu))**(-2.0d0)*(T_b - Ts_b)**(-2.0d0) + 1.0d0&
               & *sqrt(pi)*Ts_b*delta_b**2*eta_b*h_b*nu*nu**(-T_b)*sqrt(-log(nu))*&
               & (T_b - Ts_b)**(-0.5d0)*erf(sqrt(-log(nu))*sqrt(T_b - Ts_b))/log(&
               & nu)
       ELSEIF (nu==1.0) THEN
          sigma_b = 0.466666666666665d0*T_b*delta_b**3*eta_b - 0.5d0*T_b*&
               &delta_b**2*eta_b*h_b - 0.466666666666667d0*Ts_b*delta_b**3*eta_b&
               & + 0.5d0*Ts_b*delta_b**2*eta_b*h_b
       ELSE
          ERROR_CODE = 1
       ENDIF nutest
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_B

END MODULE MOBILITY_THERMAL_SKIN_RHEOLOGY
