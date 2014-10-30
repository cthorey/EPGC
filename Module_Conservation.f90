MODULE MODULE_CONSERVATION

CONTAINS

  SUBROUTINE MASS_CONSERVATION(H,Dt,dist,ray,k,BV_a,BV_b,V_t1,V_t2,delta0)

    !*****************************************************************
    ! Control the mass conservation at each time step pour le code balmforth
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:) ,INTENT(IN) :: H
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(INOUT) :: BV_a,BV_b
    DOUBLE PRECISION ,INTENT(INOUT) :: V_t1,V_t2
    INTEGER ,INTENT(IN) :: k

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: Dt,delta0

    !Parametre du sous programme
    DOUBLE PRECISION :: A_t1,Int_t1,A_t2,Int_t2,Dr
    INTEGER :: i,N

    N = COUNT(H(:,1)>delta0)
    Dr = ray(1)

    ! temps t
    A_t1 = H(1,1)
    Int_t1 = A_t1*Dr**2
    DO i=2,N-1
       A_t1 = H(i,1)
       Int_t1 = Int_t1+A_t1*(ray(i)**2-ray(i-1)**2)
    ENDDO
    A_t1 = H(N,1)
    Int_t1 = Int_t1+A_t1*(dist(N)**2-ray(N-1)**2)
    
    ! temps t+dt
    A_t2 = H(1,3)
    Int_t2 = A_t2*Dr**2
    DO i=2,N-1
       A_t2 = H(i,3)
       Int_t2 = Int_t2+A_t2*(ray(i)**2-ray(i-1)**2)
    ENDDO
    A_t2 = H(N,3)
    Int_t2 = Int_t2+A_t2*(dist(N)**2-ray(N-1)**2)

    ! Volume de l'energie au temps t+Dt
    BV_a = Int_t2
    ! Volume de l'intrusion au temp t + ce qu'on a rajoute
    BV_b = Int_t1+ Dt
    ! Volume temps t
    V_t1 = Int_t1
    ! Volume temps t+dt
    V_t2 = Int_t2

  END SUBROUTINE MASS_CONSERVATION

  SUBROUTINE ENERGY_CONSERVATION(H,BL,T,Ts,Pe,Dt,dist,ray,k,psi,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,delta0)

    !*****************************************************************
    ! Control the energy conservation at each time step pour le code balmforth
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:) ,INTENT(IN) :: Ts,H,T,BL
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(INOUT) :: BE_a,BE_b
    DOUBLE PRECISION ,INTENT(INOUT) :: En_t1,En_t2,Phi_s,Phi_l
    INTEGER ,INTENT(IN) :: k

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: Pe,Dt,psi,delta0

    !Parametre du sous programme
    DOUBLE PRECISION :: tbar_t1,tbar_t2,A_1_t1,A_1_t2,A_2,A_3
    DOUBLE PRECISION :: Int_1_t1,Int_1_t2,Int_2,Int_3
    DOUBLE PRECISION :: E_ta,E_tb,Dr
    INTEGER :: i,N
    
    N=0
    DO i=1,COUNT(H(:,3)>0D0),1
       IF (H(i,3)-delta0>0.d0 .OR. N /= 0) THEN
          CYCLE
       ELSE
          N = i
       END IF
    ENDDO

    Dr = ray(1)

    ! Calcule au temps t
    tbar_t1 = -(2*(T(1,1)-Ts(1,1))*BL(1,1))/(3.d0*H(1,1))+T(1,1)
    ! Calcule au temps t
    tbar_t2 = -(2*(T(1,3)-Ts(1,3))*BL(1,3))/(3.d0*H(1,3))+T(1,3)

    A_1_t1 = tbar_t1*H(1,1)
    A_1_t2 = tbar_t2*H(1,3)
    A_2 = (tbar_t2*H(1,3)-tbar_t1*H(1,1))/Dt
    A_3 = ((T(1,1)-Ts(1,1))/BL(1,1))

    Int_1_t1 = A_1_t1*Dr**2
    Int_1_t2 = A_1_t2*Dr**2
    Int_2 = A_2*Dr**2
    Int_3 = A_3*Dr**2
    
    DO i=2,N-1
       tbar_t1 = -(2*(T(i,1)-Ts(i,1))*BL(i,1))/(3.d0*H(i,1))+T(i,1)
       tbar_t2 = -(2*(T(i,3)-Ts(i,3))*BL(i,3))/(3.d0*H(i,3))+T(i,3)

       A_1_t1 = tbar_t1*H(1,1)
       A_1_t2 = tbar_t2*H(1,3)
       A_2 = (tbar_t2*H(1,3)-tbar_t1*H(1,1))/Dt
       A_3 = ((T(1,1)-Ts(1,1))/BL(1,1))

       Int_1_t1 = Int_1_t1 + A_1_t1*(ray(i)**2-ray(i-1)**2)
       Int_1_t2 = Int_1_t2 + A_1_t2*(ray(i)**2-ray(i-1)**2)
       Int_2 = Int_2 + A_2*(ray(i)**2-ray(i-1)**2)
       Int_3 = Int_3 + A_3*(ray(i)**2-ray(i-1)**2)
    ENDDO
    
    tbar_t1 = -(2*(T(N,1)-Ts(N,1))*BL(N,1))/(3.d0*H(N,1))+T(N,1)
    tbar_t2 = -(2*(T(N,3)-Ts(N,3))*BL(N,3))/(3.d0*H(N,3))+T(N,3)
    A_1_t1 = tbar_t1*H(N,1)
    A_1_t2 = tbar_t2*H(N,3)
    A_2 = (tbar_t2*H(N,3)-tbar_t1*H(N,1))/Dt
    A_3 = ((T(N,1)-Ts(N,1))/BL(N,1))

    Int_1_t1 = Int_1_t1 + A_1_t1*(dist(N)**2-ray(N-1)**2)
    Int_1_t2 = Int_1_t2 + A_1_t2*(dist(N)**2-ray(N-1)**2)
    Int_2 = Int_2 + A_2*(dist(N)**2-ray(N-1)**2)
    Int_3 = Int_3 + A_3*(dist(N)**2-ray(N-1)**2)

    ! Bilan d'energie
    
    !E_t2a: Energie Billan calculer a partir de l'intrusin au temps t+dt
    BE_a = Int_1_t2
    !E_t2b: Energie Bilan calculer a partir de l'intrusion aut temps t + gane perdu
    BE_b = Int_1_t1+(1D0-psi*Int_2-4D0*Pe*Int_3)*Dt
    ! Energie dans l'intrusion au temps t
    En_t1 = Int_1_t1
    ! Energie dans l'intrusion au temps t+dt
    En_t2 = Int_1_t2
    ! Energie sources en J s
    Phi_s = (1D0-psi*Int_2)
    ! Energie lost en J s
    Phi_l = 4D0*Pe*Int_3
    
    
    ! PRINT*,'Conservation chaleur',tmps,D_p-D_c,D_p
  END SUBROUTINE ENERGY_CONSERVATION

END MODULE MODULE_CONSERVATION

