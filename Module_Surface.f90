MODULE MODULE_SURFACE_TEMPERATURE

CONTAINS

  SUBROUTINE SURFACE_CONDUCTION_SOLVER(H,BL,T,Ts,N1,Pe,delta0,M,tmps,Dt)

    !*****************************************************************
    ! Solve for Ts in the case where only the convection is needed 
    !*****************************************************************

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:) ,INTENT(INOUT) :: Ts
    DOUBLE PRECISION , DIMENSION(:,:) ,INTENT(IN) :: H,T,BL

    ! Parametre du model
    INTEGER ,INTENT(IN) :: M

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: N1,Pe,delta0,tmps,Dt

    !Parametre du sous programme
    INTEGER :: size,i,N
    DOUBLE PRECISION :: beta
    DOUBLE PRECISION, PARAMETER :: pi=3.14159265

    N = COUNT(H(:,3)>delta0) 
    DO i = 1,N,1
       beta = N1*Pe**(-0.5d0)/(sqrt(pi*tmps))
       Ts(i,3) = 2.d0*T(i,3)/(beta*BL(i,3)+2.d0)
    ENDDO

  END SUBROUTINE SURFACE_CONDUCTION_SOLVER

 
END MODULE MODULE_SURFACE_TEMPERATURE

