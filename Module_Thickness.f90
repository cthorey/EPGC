MODULE MODULE_THICKNESS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  IMPORTATIONS DES MODULES

USE MODULE_THICKNESS_NEWTON
USE MODULE_THICKNESS_GFD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINES

CONTAINS

  SUBROUTINE  THICKNESS_SOLVER(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,eps_1,ERROR_CODE)

    IMPLICIT NONE
    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: H,P
    DOUBLE PRECISION, DIMENSIOn(:,:), INTENT(IN) :: T,Ts,BL
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,eps_1
    INTEGER ,INTENT(IN) :: M
    INTEGER, INTENT(INOUT) :: ERROR_CODE

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,sigma,nu,delta0
 
    ! Parametre du sous programme
    INTEGER                                                                            :: z
    DOUBLE PRECISION                                                             :: F_errt,F_err,theta

    ! Schema utilise
    theta = 1.d0

    F_err=20; F_errt=20; z=0
    THICKNESS_ITERATION: DO 
       CALL  THICKNESS_NEWTON_SOLVER(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,z,F_err,theta)
       z=z+1
       IF ( F_err>F_errt) THEN
          PRINT*,'Erreur_Ite_Epais',F_err,F_errt
       END IF
       IF (z>50000 .OR. F_Err>1D30) THEN
          ERROR_CODE = 1
          EXIT
       ENDIF
       IF (F_err<eps_1) EXIT
       H(:,2) = H(:,3); H(:,3) = delta0
       F_errt = F_err
    END DO THICKNESS_ITERATION

  END SUBROUTINE THICKNESS_SOLVER

END MODULE MODULE_THICKNESS
