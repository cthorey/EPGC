MODULE MOBILITY_THERMAL_SKIN_RHEOLOGY

CONTAINS
  
  SUBROUTINE fomega_a(Ai,h_a,delta_a,T_a,Ts_a,eta_a,omega_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ai,h_a,delta_a,T_a,Ts_a,eta_a

    DOUBLE PRECISION, INTENT(INOUT) :: omega_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
    IF (Rheology == 0) THEN
          omega_a = (7.0d0/5.0d0)*T_a*delta_a**2*eta_a*nu - 7.0d0/5.0d0*T_a*&
               & delta_a**2*eta_a - 3.0d0/2.0d0*T_a*delta_a*eta_a*h_a*nu + (3.0d0/&
               & 2.0d0)*T_a*delta_a*eta_a*h_a + (3.0d0/5.0d0)*Ts_a*delta_a**2*&
               & eta_a*nu - 3.0d0/5.0d0*Ts_a*delta_a**2*eta_a - 3.0d0/2.0d0*Ts_a*&
               & delta_a*eta_a*h_a*nu + (3.0d0/2.0d0)*Ts_a*delta_a*eta_a*h_a - 2*&
               & delta_a**2*eta_a*nu + 3*delta_a*eta_a*h_a*nu

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fomega_A

  SUBROUTINE fsigma_a(Ai,h_a,delta_a,T_a,Ts_a,eta_a,sigma_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ai,h_a,delta_a,T_a,Ts_a,eta_a

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

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_A


    SUBROUTINE fOmega_b(Bi,h_b,delta_b,T_b,Ts_b,eta_b,omega_b,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bi,h_b,delta_b,T_b,Ts_b,eta_b

    DOUBLE PRECISION, INTENT(INOUT) :: omega_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
    IF (Rheology == 0) THEN
          omega_b = (7.0d0/5.0d0)*T_b*delta_b**2*eta_b*nu - 7.0d0/5.0d0*T_b*&
               & delta_b**2*eta_b - 3.0d0/2.0d0*T_b*delta_b*eta_b*h_b*nu + (3.0d0/&
               & 2.0d0)*T_b*delta_b*eta_b*h_b + (3.0d0/5.0d0)*Ts_b*delta_b**2*&
               & eta_b*nu - 3.0d0/5.0d0*Ts_b*delta_b**2*eta_b - 3.0d0/2.0d0*Ts_b*&
               & delta_b*eta_b*h_b*nu + (3.0d0/2.0d0)*Ts_b*delta_b*eta_b*h_b - 2*&
               & delta_b**2*eta_b*nu + 3*delta_b*eta_b*h_b*nu

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FOmega_B


  SUBROUTINE fsigma_b(Bi,h_b,delta_b,T_b,Ts_b,sigma_b,eta_b,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bi,h_b,delta_b,T_b,Ts_b,eta_b

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

    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE Fsigma_B

END MODULE MOBILITY_THERMAL_SKIN_RHEOLOGY
