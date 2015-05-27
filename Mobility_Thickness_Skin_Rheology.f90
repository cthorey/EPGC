MODULE MOBILITY_THICKNESS_SKIN_RHEOLOGY

CONTAINS
  
  SUBROUTINE fPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,phi_a,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ael,Agrav,h_a,delta_a,T_a,Ts_a

    DOUBLE PRECISION, INTENT(INOUT) :: phi_a
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
    IF (Rheology == 0) THEN
       phi_a = (4.0d0/5.0d0)*T_a*delta_a**3*nu - 4.0d0/5.0d0*T_a*delta_a&
            & **3 - 2*T_a*delta_a**2*h_a*nu + 2*T_a*delta_a**2*h_a + 2*T_a*&
            & delta_a*h_a**2*nu - 2*T_a*delta_a*h_a**2 - T_a*h_a**3*nu + T_a*&
            & h_a**3 - 4.0d0/5.0d0*Ts_a*delta_a**3*nu + (4.0d0/5.0d0)*Ts_a*&
            & delta_a**3 + 2*Ts_a*delta_a**2*h_a*nu - 2*Ts_a*delta_a**2*h_a - 2&
            & *Ts_a*delta_a*h_a**2*nu + 2*Ts_a*delta_a*h_a**2 + h_a**3*nu       
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_A

  SUBROUTINE  fdPhi_A(Ael,Agrav,h_a,delta_a,T_a,Ts_a,hi,hia,&
       &dphia_dhi1,dphia_dhi,nu,Rheology,ERROR_CODE)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Ael,Agrav,h_a,T_a,delta_a,Ts_a
    DOUBLE PRECISION, INTENT(IN) :: hi,hia

    DOUBLE PRECISION, INTENT(INOUT) :: dphia_dhi,dphia_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
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
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FdPhi_A

  SUBROUTINE  fPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,phi_b,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bel,Bgrav,h_b,delta_b,T_b,Ts_b

    DOUBLE PRECISION, INTENT(INOUT) :: phi_b
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
    IF (Rheology == 0) THEN
       
       phi_b = (4.0d0/5.0d0)*T_b*delta_b**3*nu - 4.0d0/5.0d0*T_b*delta_b&
          & **3 - 2*T_b*delta_b**2*h_b*nu + 2*T_b*delta_b**2*h_b + 2*T_b*&
          & delta_b*h_b**2*nu - 2*T_b*delta_b*h_b**2 - T_b*h_b**3*nu + T_b*&
          & h_b**3 - 4.0d0/5.0d0*Ts_b*delta_b**3*nu + (4.0d0/5.0d0)*Ts_b*&
          & delta_b**3 + 2*Ts_b*delta_b**2*h_b*nu - 2*Ts_b*delta_b**2*h_b - 2&
          & *Ts_b*delta_b*h_b**2*nu + 2*Ts_b*delta_b*h_b**2 + h_b**3*nu
      
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FPhi_B
  
  SUBROUTINE  fDPhi_B(Bel,Bgrav,h_b,delta_b,T_b,Ts_b,hib,hi,&
       &dphib_dhi,dphib_dhi1,nu,Rheology,ERROR_CODE)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Rheology
    DOUBLE PRECISION , INTENT(IN) :: nu
    
    DOUBLE PRECISION, INTENT(IN) :: Bel,Bgrav,h_b,T_b,delta_b
    DOUBLE PRECISION, INTENT(IN) :: Ts_b,hi,hib

    DOUBLE PRECISION, INTENT(INOUT) :: dphib_dhi,dphib_dhi1
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    
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
      
    ELSEIF (Rheology == 1) THEN
       ERROR_CODE = 1
    ELSEIF (Rheology == 2) THEN
       ERROR_CODE = 1
    ELSE
       ERROR_CODE = 1
    ENDIF
  END SUBROUTINE FDPhi_B


END MODULE MOBILITY_THICKNESS_SKIN_RHEOLOGY
