MODULE MODULE_THERMAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  IMPORTATIONS DES MODULES
  
USE MODULE_THERMAL_SKIN_NEWTON_BERCOVICI
USE MODULE_THERMAL_SKIN_GFD_BERCOVICI
USE MODULE_THERMAL_INTE_NEWTON_BERCOVICI
USE MODULE_THERMAL_SKIN_NEWTON_ROSCOE
USE MODULE_THERMAL_SKIN_NEWTON_ARRHENIUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINES
CONTAINS

  SUBROUTINE THERMAL_SOLVER(Xi,H,T,Ts,BL,P,M,dist,ray,sigma,nu,Pe,psi,delta0,el,grav,gam,Inter_Q,Dr,Dt,eps_1,k,N1,tmps,&
       &ERROR_CODE,Model,Schema,Rheology)

    IMPLICIT NONE
    ! Tableaux
    DOUBLE PRECISION , DIMENSIOn(:,:) , INTENT(INOUT) :: H,P
    DOUBLE PRECISION , DIMENSION(:,:) ,INTENT(INOUT) :: Xi,T,Ts,BL
    DOUBLE PRECISION , DIMENSION(:) ,INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,eps_1,N1,tmps
    INTEGER ,INTENT(IN) :: M,k
    INTEGER ,INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Model,Schema,Rheology
    
    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: sigma,nu,Pe,psi,delta0,el,grav,gam,Inter_Q
 
    ! Parametre du sous programme
    INTEGER :: z,N,ndyke,i,Size
    DOUBLE PRECISION :: F_Errt,F_err,theta,F2,F2t,s
    LOGICAL :: cho

    ! Schema utilise
    theta = 1D0
    F_err = 100; F_errt = 100
    z =0

    THERMAL_ITERATION: DO
       ! CALL SUBROUTINE
       IF (Model == 1 .AND. Schema == 0 .AND. Rheology == 0) THEN
          CALL  THERMAL_SKIN_NEWTON_BERCOVICI(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)
       ELSEIF (Model == 1 .AND. Schema == 0 .AND. Rheology == 1) THEN
          CALL  THERMAL_SKIN_NEWTON_ROSCOE(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)
       ELSEIF (Model == 1 .AND. Schema == 0 .AND. Rheology == 2) THEN
          CALL  THERMAL_SKIN_NEWTON_ARRHENIUS(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)
       ELSEIF (Model == 1 .AND. Schema == 1 .AND. Rheology == 0) THEN
          CALL  THERMAL_SKIN_GFD_BERCOVICI(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,delta0,el,grav,N1,F_err,z,tmps)
       ELSEIF (Model == 0 .AND. Schema == 1 .AND. Rheology == 0) THEN
          CALL  THERMAL_INTE_NEWTON_BERCOVICI(Xi,H,P,T,Ts,BL,Dt,Dr,theta,dist,ray,M,sigma,nu,Pe,psi,&
               &delta0,el,grav,N1,gam,Inter_Q,F_err,z,tmps)
       ELSE
          PRINT*,'PAS DE MODULE THICKNESS  CORRESPONDANT IMPLEMENTE ENCORE'
          ERROR_CODE = 1
       ENDIF

       ! ITERATIVE PROCEDURE
       
       z=z+1
       ! PRINT*,F_err,z,'Temperature'
       IF (F_err>F_errt) THEN
          PRINT*,tmps,z,'Erreur_Ite_Temp',F_err,F_errt
       ENDIF
       IF (z>2000) THEN
          PRINT*,z
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
