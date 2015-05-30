MODULE MODULE_THICKNESS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  IMPORTATIONS DES MODULES

  USE MODULE_THICKNESS_INTE_GFD_BERCOVICI
  USE MODULE_THICKNESS_SKIN_NEWTON
  USE MODULE_THICKNESS_SKIN_NEWTON_BERCOVICI
  USE MODULE_THICKNESS_SKIN_NEWTON_GRAVI
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SUBROUTINES

CONTAINS

  SUBROUTINE  THICKNESS_SOLVER(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,gam&
       &,Inter_Q,eps_1,ERROR_CODE,Model,Schema,Rheology,tmps)

    IMPLICIT NONE
    ! Tableaux
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: H,P
    DOUBLE PRECISION, DIMENSIOn(:,:), INTENT(IN) :: T,Ts,BL
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: dist,ray

    ! Parametre du model
    DOUBLE PRECISION ,INTENT(IN) :: Dt,Dr,eps_1
    DOUBLE PRECISION, INTENT(IN) :: tmps
    INTEGER ,INTENT(IN) :: M
    INTEGER, INTENT(INOUT) :: ERROR_CODE
    INTEGER, INTENT(IN) :: Model,Schema,Rheology

    ! Nombre sans dimension
    DOUBLE PRECISION ,INTENT(IN) :: el,grav,sigma,nu,delta0,gam,Inter_Q
 
    ! Parametre du sous programme
    INTEGER                                                                            :: z
    DOUBLE PRECISION                                                             :: F_errt,F_err,theta

    ! Schema utilise
    theta = 1.d0

    F_err=20; F_errt=20; z=0
    THICKNESS_ITERATION: DO
       
       IF (Model == 1 .AND. Schema == 0 .AND. el == 1D0) THEN
          CALL  THICKNESS_SKIN_NEWTON(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,&
               &gam,Inter_Q,z,F_err,theta,tmps,Rheology,ERROR_CODE)
       ELSEIF (Model == 1 .AND. Schema ==0 .AND. el == 0D0) THEN
          CALL THICKNESS_SKIN_NEWTON_GRAVI(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma&
               &,nu,delta0,z,F_err,theta,Rheology,ERROR_CODE)
          
       ELSE
          PRINT*,'PAS DE MODULE CORRESPONDANT IMPLEMENTE ENCORE THICKNESS'
          PRINT*,'MODEL =',Model,'SCHEMA =',Schema,'Rheology =',Rheology
          ERROR_CODE = 1
       ENDIF

       ! ITERATIVE PROCEDURE
       z=z+1
       IF ( F_err>F_errt) THEN
          ! PRINT*,'Erreur_Ite_Epais',F_err,F_errt
       END IF
       IF (z>50000 .OR. ERROR_CODE == 1) THEN
          ERROR_CODE = 1
          PRINT*,z
          EXIT
       ENDIF

       IF (F_err<eps_1) EXIT
       H(:,2) = H(:,3); H(:,3) = delta0
       F_errt = F_err
    END DO THICKNESS_ITERATION

  END SUBROUTINE THICKNESS_SOLVER

  SUBROUTINE ADVANCE_RADIUS_TEST(H,P,T,BL,Ts,Xi,R_Intrusion)
    
    ! Thcheuqe is lintrusion avant dans le regime gravi
    IMPLICIT NONE
    ! Tableaux
    
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: H,P,T,Ts,BL,Xi
    INTEGER, INTENT(INOUT) :: R_Intrusion

    IF (R_Intrusion /= COUNT(H(:,3)>0.d0)) THEN
       IF (COUNT(H(:,3)-R_Intrusion>0.d0)>1) THEN
          BL(R_Intrusion:COUNT(H(:,3)>0.d0),1) = 1D-4
          BL(R_Intrusion:COUNT(H(:,3)>0.d0),2) = 1D-4
          T(R_Intrusion:COUNT(H(:,3)>0.d0),1) = 1D0
          T(R_Intrusion:COUNT(H(:,3)>0.d0),2) = 1D0
          Xi(R_Intrusion:COUNT(H(:,3)>0.d0),1) = 1D-4/3D0
          Xi(R_Intrusion:COUNT(H(:,3)>0.d0),2) = 1D-4/3D0
       ELSE
          BL(COUNT(H(:,3)>0.d0),1) = BL(COUNT(H(:,3)>0.d0)-1,1)
          BL(COUNT(H(:,3)>0.d0),2) = BL(COUNT(H(:,3)>0.d0)-1,2)
          T(COUNT(H(:,3)>0.d0),1) = T(COUNT(H(:,3)>0.d0)-1,1)
          T(COUNT(H(:,3)>0.d0),2) = T(COUNT(H(:,3)>0.d0)-1,2)
          Xi(COUNT(H(:,3)>0.d0),1) = Xi(COUNT(H(:,3)>0.d0)-1,1)
          Xi(COUNT(H(:,3)>0.d0),2) = Xi(COUNT(H(:,3)>0.d0)-1,2)
       ENDIF
       R_Intrusion = COUNT(H(:,3)>0.d0)
    ENDIF
  END SUBROUTINE ADVANCE_RADIUS_TEST

  
END MODULE MODULE_THICKNESS
