MODULE MODULE_THERMAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  IMPORTATIONS DES MODULES
USE MODULE_THERMAL_NEWTON
USE MODULE_THERMAL_GFD
USE MODULE_SURFACE_TEMPERATURE
USE MODULE_THERMAL_NEWTON_OLD
USE MODULE_THERMAL_NEWTON_CRYS_2
USE MODULE_THERMAL_NEWTON_INT_EPAISSEUR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINES
CONTAINS

  SUBROUTINE THERMAL_SOLVER(Xi,H,T,Ts,BL,P,M,dist,ray,sigma,nu,Pe,psi,delta0,el,grav,Dr,Dt,eps_1,k,N1,tmps,ERROR_CODE)

    IMPLICIT NONE
    ! Tableaux
    DOUBLE PRECISION , DIMENSIOn(:,:) , INTENT(INOUT) :: H,P
    DOUBLE PRECISION , DIMENSION(:,:) ,INTENT(INOUT) :: Xi,T,Ts,BL
    DOUBLE PRECISION , DIMENSION(:) ,INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,eps_1,N1,tmps
    INTEGER ,INTENT(IN) :: M,k
    INTEGER ,INTENT(INOUT) :: ERROR_CODE

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: sigma,nu,Pe,psi,delta0,el,grav
 
    ! Parametre du sous programme
    INTEGER :: z,N,ndyke,i,Size
    DOUBLE PRECISION :: F_Errt,F_err,theta,F2,F2t,s
    LOGICAL :: cho

    ! Schema utilise
    theta = 1D0
    F_err = 20; F_errt = 20
    z =0

    THERMAL_ITERATION: DO 
       CALL THERMAL_NEWTON_SOLVER(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)
       z=z+1
       IF (F_err>F_errt) THEN
          ! PRINT*,tmps,z,'Erreur_Ite_Temp',F_err,F_errt
       ENDIF
       IF (z>20000) THEN
          ERROR_CODE = 1
          EXIT
       ENDIF
       IF (F_err<eps_1) EXIT

       Xi(:,2) = Xi(:,3)
       T(:,2) = T(:,3)
       Ts(:,2) = Ts(:,3)
       BL(:,2) = BL(:,3)
       
       F_errt = F_err
      
    END DO THERMAL_ITERATION

  END SUBROUTINE THERMAL_SOLVER
 
END MODULE MODULE_THERMAL
