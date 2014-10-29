MODULE MODULE_CONSERVATION

CONTAINS

  SUBROUTINE ENERGY_CONSERVATION(H,BL,T,Ts,Xi,E_c,Phi_s,Pe,M,tmps,Dt,dist,ray,k)

    !*****************************************************************
    ! Control the energy conservation at each time step
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:) ,INTENT(IN) :: Xi,Ts,H,T,BL
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(INOUT) :: E_c,Phi_s
    INTEGER ,INTENT(IN) :: M,k

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: Pe,tmps,Dt

    !Parametre du sous programme
    DOUBLE PRECISION :: E_n,D_p,D_c
    INTEGER :: i
    
    E_n = Xi(1,3)
    Phi_s = ((T(1,3)-Ts(1,3))/BL(1,3))
    DO i=2,M,-1
       E_n = E_n + Xi(i,3)*(ray(i)**2-ray(i-1)**2)
       Phi_s = Phi_s+((T(i,3)-Ts(i,3))/BL(i,3))*(ray(i)**2-ray(i-1)**2)
    ENDDO
    E_n = E_n + Xi(M,3)*(dist(M)**2-ray(M-1)**2)
    Phi_s = Phi_s +((T(M,3)-Ts(M,3))/BL(M,3))*(dist(M)**2-ray(M-1)**2)
    
    D_p = (E_n-E_c)
    D_c = (1.D0-4.D0*Pe*Phi_s)*Dt
    E_c = E_n
  END SUBROUTINE ENERGY_CONSERVATION

SUBROUTINE ENERGY_CONSERVATION_2(H,BL,T,Ts,Xi,E_c,Phi_s,Pe,M,tmps,Dt,dist,ray,k)

    !*****************************************************************
    ! Control the energy conservation at each time step pour le code balmforth
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:) ,INTENT(IN) :: Xi,Ts,H,T,BL
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(INOUT) :: E_c,Phi_s
    INTEGER ,INTENT(IN) :: M,k

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: Pe,tmps,Dt

    !Parametre du sous programme
    DOUBLE PRECISION :: E_n,D_p,D_c,hthetabar
    INTEGER :: i
    
    hthetabar = -2*(T(1,3)-Ts(1,3))/3.d0*BL(1,3)+T(1,3)*H(1,3)
    E_n = hthetabar
    Phi_s = ((T(1,3)-Ts(1,3))/BL(1,3))
    DO i=2,M,-1
       hthetabar = -2*(T(i,3)-Ts(i,3))/3.d0*BL(i,3)+T(i,3)*H(i,3)
       E_n = E_n + hthetabar*(ray(i)**2-ray(i-1)**2)
       Phi_s = Phi_s+((T(i,3)-Ts(i,3))/BL(i,3))*(ray(i)**2-ray(i-1)**2)
    ENDDO
    hthetabar = -2*(T(M,3)-Ts(M,3))/3.d0*BL(M,3)+T(M,3)*H(M,3)
    E_n = E_n + hthetabar*(dist(M)**2-ray(M-1)**2)
    Phi_s = Phi_s +((T(M,3)-Ts(M,3))/BL(M,3))*(dist(M)**2-ray(M-1)**2)
    
    D_p = E_n
    D_c = (1.D0-4.D0*Pe*Phi_s)*Dt+E_c
    E_c = E_n
    ! PRINT*,'Conservation chaleur',tmps,D_p-D_c,D_p
  END SUBROUTINE ENERGY_CONSERVATION_2

END MODULE MODULE_CONSERVATION

