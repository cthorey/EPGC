PROGRAM MAIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  IMPORTATION DES MODULES

  USE MODULE_INITIALISATION
  USE MODULE_OUTPUT
  USE MODULE_THICKNESS
  USE MODULE_THERMAL
  USE MODULE_CONSERVATION
  USE MODULE_COMPLEMENTAIRE
  USE MODULE_INTEGRATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  DECLARATION DES VARIABLES

  IMPLICIT NONE
  ! Tableau pour stocker les donne
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H,T,Xi,BL,Ts,P
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Srr,Stt
  DOUBLE PRECISION , DIMENSION(:),ALLOCATABLE  :: hmubar,hthetabar,ubar

  ! Tableux pour definir la grille
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dist,ray

  ! Parametre du model
  INTEGER :: M
  DOUBLE PRECISION :: tmps_m,Dt,Dr,tmps,tmps_n
  INTEGER :: sample,Init
  INTEGER :: k,k1,k2,z,compteur,i

  ! SCHEMA A UTILISER
  INTEGER Model,T_Schema,H_Schema,Rheology
  
  ! Nombre sans dimension
  DOUBLE PRECISION :: el,grav,delta0,sigma,nu,Pe,psi,N1,gam

  ! Conservation quantityes
  ! Mass
  DOUBLE PRECISION  :: BV_a,BV_b,V_t1,V_t2
  ! Energie
  DOUBLE PRECISION  :: BE_a,BE_b ! Voir module conservation pour definition
  DOUBLE PRECISION  :: En_t1,En_t2,Phi_s,Phi_l
  
  !Average quantity
  DOUBLE PRECISION :: Tm,Vm,Mum,Phim
  DOUBLE PRECISION :: Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005
  DOUBLE PRECISION :: Tm01,Tm02,Tm05,Tm005
  ! Seuil de confiance
  DOUBLE PRECISION :: eps_1

  !Condition sur boulce + compteur
  DOUBLE PRECISION :: F1,F1t,F2,F2t
  INTEGER :: Ite_Glob

  ! Trackjing front routine
  DOUBLE PRECISION :: Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu
  DOUBLE PRECISION :: Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H
  DOUBLE PRECISION :: Fr_005_R,Fr_005_T,Fr_005_Mu
  DOUBLE PRECISION :: Mu_e

  !Necessaire au bon fonctionnement
  INTEGER :: cas,size
  INTEGER :: ERROR_CODE
  INTEGER :: Naa

  ! String variable pour les nom fichiers
  CHARACTER(LEN=300) :: Output_racine,Input_Racine,Input_Data_Name
  CHARACTER(LEN=300) :: Format_NSD,Format_O
  CHARACTER(LEN=300) :: Format_NSD_Init_0,Format_NSD_Init_1
  CHARACTER(LEN=300) :: Format_Input_Data,Format_RV,Format_Backup
  CHARACTER(LEN=300) :: Format_M
  CHARACTER(LEN=300) :: Format_NF,NF,Root_Code !File name output
  CHARACTER(LEN=300) :: NF_Name
  ! Associaion de varialbe
  DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v
  common /arg/ delta,Thetas,Thetab,nu_v

  ! Parametre pour erreur d'allocaiton

  INTEGER :: err1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  DEBUT DU PROGRAMME

  ! Call subroutine const dans le module CONSTANTE
  CALL  CONSTANTE(M,tmps_m,Dt,Dr,sample,el,grav,delta0,sigma,nu,Pe,psi,N1,gam,&
       &eps_1,Format_O,Format_NSD,Init,Input_Racine,Output_Racine,Input_Data_Name&
       &,Format_NSD_Init_0,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup,&
       &NF,Format_NF,Root_Code,Model,T_Schema,H_Schema,Rheology)
  
  ! Allocation des tableaux
  ALLOCATE(H(1:M,4),Xi(1:M,4),Ts(1:M,4),BL(1:M,4),T(1:M,4),&
       &P(1:M,4),ray(1:M),dist(1:M), stat=err1)
  ALLOCATE(Srr(1:M),Stt(1:M), stat = err1)
  ALLOCATE(hmubar(1:M),hthetabar(1:M),ubar(1:M), stat = err1)
  IF (err1>1) THEN
     PRINT*,'ERREUR endallocation tabelaux dans le main'
     STOP
  ENDIF

  ! Calll initialisation dans le cmodule Module_Initialisation
  CALL INITIALISATION(Format_O,Format_NSD,M,H,T,Ts,Xi,BL,P,dist,ray,k,k1,k2,z,tmps,&
       &Dt,Dr,eps_1,el,grav,delta0,sigma,nu,Pe,Psi,N1,gam,sample,Init,compteur,tmps_m,&
       &Input_Data_Name,Input_racine,Output_Racine,NF,Format_NF,Root_Code&
       &,Format_NSD_Init_0,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup&
       &,Model,T_Schema,H_Schema,Rheology)

  ! Debut de la boucle sur le temps
  tmps_n = 0D0
  ERROR_CODE = 0
  TEMPS: DO WHILE (tmps<tmps_m)
     
     IF (tmps>1D-1) THEN !Attention a l'ordre, il faut pas qu'il sort de la boucle avant d'avoir pu changer !
        Dt = 1D-6
     ELSEIF (tmps >10D0) THEN
        Dt = 1D-5
     ENDIF

     OPEN(unit=3,file='Output.txt')

     Cas=0
     CALL OUTPUT(Format_O,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup,&
          &NF,Format_NF,Root_Code,&
          &Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,&
          &Tm01,Tm02,Tm05,Tm005,&
          &Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e,&
          &Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H,hmubar,hthetabar,ubar,&
          &Fr_005_R,Fr_005_T,Fr_005_Mu,tmps_n)


     H(:,1) = H(:,3); H(:,4) = H(:,1)
     Xi(:,1) = Xi(:,3); Xi(:,4) = Xi(:,1)
     T(:,1) = T(:,3); BL(:,1)=BL(:,3)
     Ts(:,1) = Ts(:,3); Ts(:,4) = Ts(:,1)

     Ite_Glob = 0
     F1t = 20; F1 = 20
     F2t = 20; F2 = 20

     ITERATION_GLOBALE: DO

        H(:,2) = H(:,3); Xi(:,2) = Xi(:,3); T(:,2) = T(:,3); BL(:,2) = BL(:,3); Ts(:,2)= Ts(:,3)
        
        ! Module Epaisseur
        CALL THICKNESS_SOLVER(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,gam,eps_1,&
             &ERROR_CODE,Model,H_Schema,Rheology)
        IF (ERROR_CODE == 1) THEN
           WRITE(NF_Name,Format_NF),Root_Code,NF,'_BUG'
           OPEN(unit =1,file=NF_Name, action ="write",status ="replace")
           WRITE(1,*)'Le CODE A PLANTE DANS LE SOLVER THICKNESS',tmps
           PRINT*,'Le CODE A PLANTE DANS LE SOLVER THICKNESS',tmps
           CLOSE(1)
           STOP
        END IF
        F1 = ABS(MAXVAL((H(:,4)-H(:,3))/(H(:,3))))

        ! Module heat transport
        CALL  THERMAL_SOLVER(Xi,H,T,Ts,BL,P,M,dist,ray,sigma,nu,Pe,psi,delta0,el,grav,gam,Dr,Dt,eps_1,k,N1,tmps,&
             &ERROR_CODE,Model,T_Schema,Rheology)
        IF (ERROR_CODE == 1) THEN
           WRITE(NF_Name,Format_NF),Root_Code,NF,'_BUG'
           OPEN(unit =1,file=NF_Name, action ="write",status ="replace")
           WRITE(1,*)'Le CODE A PLANTE DANS LE THERMAL SOLVER',tmps
           PRINT*,'Le CODE A PLANTE DANS LE THERMAL THICKNESS',tmps
           CLOSE(1)
           STOP
        END IF
        Size = COUNT(Xi(:,3)>1D-10)
        F2 = ABS(MAXVAL((Xi(:Size,4)-Xi(:Size,3))/(Xi(:Size,3))))


        !Condition d'arret  ou de sorti de boucle
        IF (Ite_Glob>20000) THEN
           WRITE(NF_Name,Format_NF),Root_Code,NF,'_BUG'
           OPEN(unit =1,file=NF_Name, action ="write",status ="replace")
           WRITE(1,*)'Le CODE A PLANTE DANS LE MAIN BOUCLE',tmps
           PRINT*,'Le CODE A PLANTE DANS LE MAIN BOUCLE',tmps
           CLOSE(1)
           STOP
        END IF
        IF (F1<eps_1 .AND. F2<eps_1 ) EXIT

        F1t = F1; H(:,4) = H(:,3)
        F2t = F2; Xi(:,4) = Xi(:,3)
        Ite_Glob = Ite_Glob+1

     END DO ITERATION_GLOBALE

     CALL MASS_CONSERVATION(H,Dt,dist,ray,k,BV_a,BV_b,V_t1,V_t2,delta0)
     CALL ENERGY_CONSERVATION(H,BL,T,Ts,Pe,Dt,dist,ray,k,psi,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,delta0,sigma)
     CALL STRESS_ELASTIC_FIELD(Srr,Stt,H,dist,Dr,M)
     CALL AVERAGE_QUANTITY(Xi,H,T,Ts,BL,dist,ray,Dt,Dr,el,grav,N1,Pe,Psi,nu,Tm,Vm,Mum,Phim,M,tmps,delta0,&
       &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,Tm01,Tm02,Tm05,Tm005)
     CALL TRACKING_FRONT(Xi,H,T,Ts,BL,dist,ray,P,Dt,Dr,el,grav,N1,Pe,Psi,nu,tmps,delta0,&
       &Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e,&
       &Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H,hmubar,hthetabar,ubar,&
       &Fr_005_R,Fr_005_T,Fr_005_Mu)

     Cas = 1
     CALL OUTPUT(Format_O,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup,&
          &NF,Format_NF,Root_Code,&
          &Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,&
          &Tm01,Tm02,Tm05,Tm005,&
          &Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e,&
          &Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H,hmubar,hthetabar,ubar,&
          &Fr_005_R,Fr_005_T,Fr_005_Mu,tmps_n)

     ! On incremente les compteurs et le temps
     k = k+1

     tmps_n = tmps
     tmps = tmps+Dt
     print*,tmps

  END DO TEMPS

  DEALLOCATE(H,T,Xi,BL,Ts)
  DEALLOCATE(dist,ray)

END PROGRAM MAIN

