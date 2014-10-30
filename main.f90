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

  ! Tableux pour definir la grille
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dist,ray

  ! Parametre du model
  INTEGER :: M
  DOUBLE PRECISION :: tmps_m,Dt,Dr,tmps
  INTEGER :: sample,Init
  INTEGER :: k,k1,k2,z,compteur,i

  ! Nombre sans dimension
  DOUBLE PRECISION :: el,grav,delta0,sigma,nu,Pe,psi,N1

  ! Conservation quantityes
  ! Mass
  DOUBLE PRECISION  :: BV_a,BV_b,V_t1,V_t2
  ! Energie
  DOUBLE PRECISION  :: BE_a,BE_b ! Voir module conservation pour definition
  DOUBLE PRECISION  :: En_t1,En_t2,Phi_s,Phi_l
  
  !Average quantity
  DOUBLE PRECISION :: Tm,Vm,Mum,Phim
  DOUBLE PRECISION :: Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005
  ! Seuil de confiance
  DOUBLE PRECISION :: eps_1

  !Condition sur boulce + compteur
  DOUBLE PRECISION :: F1,F1t,F2,F2t
  INTEGER :: Ite_Glob

  !Necessaire au bon fonctionnement
  INTEGER :: cas,size

  ! String variable pour les nom fichiers
  CHARACTER(LEN=300) :: Output_racine,Input_Racine,Input_Data_Name
  CHARACTER(LEN=300) :: Format_NSD,Format_O
  CHARACTER(LEN=300) :: Format_NSD_Init_0,Format_NSD_Init_1
  CHARACTER(LEN=300) :: Format_Input_Data,Format_RV,Format_Backup
  CHARACTER(LEN=300) :: Format_M

  ! Associaion de varialbe
  DOUBLE PRECISION :: delta,Thetas,Thetab,nu_v
  common /arg/ delta,Thetas,Thetab,nu_v

  ! Parametre pour erreur d'allocaiton

  INTEGER :: err1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  DEBUT DU PROGRAMME

  ! Call subroutine const dans le module CONSTANTE
  CALL   CONSTANTE(M,tmps_m,Dt,Dr,sample,el,grav,delta0,sigma,nu,Pe,Psi,N1,&
       &eps_1,Format_O,Format_NSD,Init,Input_Racine,Output_Racine,Input_Data_Name&
       &,Format_NSD_Init_0,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup)

  ! Allocation des tableaux
  ALLOCATE(H(1:M,4),Xi(1:M,4),Ts(1:M,4),BL(1:M,4),T(1:M,4),&
       &P(1:M,4),ray(1:M),dist(1:M), stat=err1)
  ALLOCATE(Srr(1:M),Stt(1:M), stat = err1)
  IF (err1>1) THEN
     PRINT*,'ERREUR endallocation tabelaux dans le main'
     STOP
  ENDIF

  ! Calll initialisation dans le cmodule Module_Initialisation
  CALL INITIALISATION(Format_O,Format_NSD,M,H,T,Ts,Xi,BL,P,dist,ray,k,k1,k2,z,tmps,&
       &Dt,Dr,eps_1,el,grav,delta0,sigma,nu,Pe,psi,N1,sample,Init,compteur,tmps_m,&
       &Input_Data_Name,Input_racine,Output_Racine&
       &,Format_NSD_Init_0,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup)

  ! Debut de la boucle sur le temps

  TEMPS: DO WHILE (tmps<tmps_m)

     OPEN(unit=3,file='Output.txt')

     Cas=0
     CALL OUTPUT(Format_O,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup&
          &,Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l)


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
        CALL THICKNESS_SOLVER(H,P,T,BL,Ts,Dt,Dr,M,dist,ray,el,grav,sigma,nu,delta0,eps_1)
        F1 = ABS(MAXVAL((H(:,4)-H(:,3))/(H(:,3))))

        ! Module heat transport
        CALL  THERMAL_SOLVER(Xi,H,T,Ts,BL,P,M,dist,ray,sigma,nu,Pe,psi,delta0,el,grav,Dr,Dt,eps_1,k,N1,tmps)
        Size = COUNT(Xi(:,3)>1D-10)
        F2 = ABS(MAXVAL((Xi(:Size,4)-Xi(:Size,3))/(Xi(:Size,3))))


        !Condition d'arret  ou de sorti de boucle
        IF (Ite_Glob>20000) THEN
           PRINT*,'Erreur Ite Global',F1,F2
           STOP
        END IF

        ! print*,'Ite_Glob',F1,F2,Ite_Glob
        IF (F1<eps_1 .AND. F2<eps_1 ) EXIT

        F1t = F1; H(:,4) = H(:,3)
        F2t = F2; Xi(:,4) = Xi(:,3)
        Ite_Glob = Ite_Glob+1

     END DO ITERATION_GLOBALE

     CALL MASS_CONSERVATION(H,Dt,dist,ray,k,BV_a,BV_b,V_t1,V_t2,delta0)
     CALL ENERGY_CONSERVATION(H,BL,T,Ts,Pe,Dt,dist,ray,k,psi,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,delta0)
     CALL STRESS_ELASTIC_FIELD(Srr,Stt,H,dist,Dr,M)
     CALL AVERAGE_QUANTITY(Xi,H,T,Ts,BL,dist,ray,Dt,Dr,el,grav,N1,Pe,Psi,nu,Tm,Vm,Mum,Phim,M,tmps,delta0,&
       &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005)
     
     ! print*,k,Mum/Vm,Tm/Vm,Vm,tmps
     ! Ecriture des donne dans un fichier
     Cas = 1
     CALL OUTPUT(Format_O,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup&
          &,Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l)

     ! On incremente les compteurs et le temps
     k = k+1
     ! PRINT*,'ite',k,tmps,Mum/Vm,Mum01/Vm01,Mum005/Vm005,Mum05/Vm05
     print*,'ite',k,Phi_l,Phi_s
     tmps = tmps+Dt

  END DO TEMPS

  DEALLOCATE(H,T,Xi,BL,Ts)
  DEALLOCATE(dist,ray)

END PROGRAM MAIN

