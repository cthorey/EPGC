MODULE MODULE_OUTPUT

CONTAINS

  SUBROUTINE OUTPUT(Format,time_frame,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup,&
          &NF,Format_NF,Root_Code,&
          &Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,&
          &Tm01,Tm02,Tm05,Tm005,&
          &Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e,&
          &Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H,hmubar,hthetabar,ubar,&
          &Fr_005_R,Fr_005_T,Fr_005_Mu,tmps_n,R_Intrusion)

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: H,T,Ts,Xi,BL,P
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: dist,ray,Srr,Stt,hmubar,hthetabar,ubar
    DOUBLE PRECISION ,DIMENSION(:) ,INTENT(IN) :: time_frame

    ! Parametre du model
    DOUBLE PRECISION , INTENT(INOUT) :: tmps,delta0,Dt,tmps_n
    DOUBLE PRECISION , INTENT(INOUT) :: Phim,Vm,Tm,Mum
    DOUBLE PRECISION , INTENT(INOUT) :: Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005
    DOUBLE PRECISION , INTENT(INOUT) :: Tm01,Tm02,Tm05,Tm005
    DOUBLE PRECISION , INTENT(INOUT) :: Fr_d_R,Fr_d_T,Fr_d_Mu,Fr_001_R,Fr_001_T,Fr_001_Mu,Mu_e
    DOUBLE PRECISION , INTENT(INOUT) :: Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H
    DOUBLE PRECISION , INTENT(INOUT) :: Fr_005_R,Fr_005_T,Fr_005_Mu
    DOUBLE PRECISION ,INTENT(INOUT) :: BE_a,BE_b
    DOUBLE PRECISION ,INTENT(INOUT) :: En_t1,En_t2,Phi_s,Phi_l
    DOUBLE PRECISION ,INTENT(INOUT) :: BV_a,BV_b
    DOUBLE PRECISION ,INTENT(INOUT) :: V_t1,V_t2

    INTEGER , INTENT(INOUT) :: k,k1,k2,z,M,compteur,R_Intrusion

    ! Format files
    CHARACTER(LEN=300) , INTENT(IN) :: Format,Output_Racine
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_RV,Format_Backup
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_NF,NF,Root_Code
    LOGICAL :: FILE_EXISTS
    INTEGER , INTENT(IN) :: Cas,sample

    !Variable specitfique au sous programmes
    INTEGER :: i,comp_n,comp_o,samp
    CHARACTER(LEN=300) :: Output_Name_F,NF_Name
    CHARACTER(LEN=300) :: Data_file,Format_Data
    DOUBLE PRECISION :: R,hmax,power,tmps_Avant


    SELECT CASE (Cas)

!!! BACKUP FILE
    CASE(0)

!!!Compteur
       IF (tmps >= time_frame(compteur)) THEN
          compteur = compteur +1
       ENDIF

       WRITE(NF_Name,Format_NF),Root_Code,NF,'.txt'
       OPEN(unit =2,file=NF_Name, action ="write",status ="replace")
       WRITE(2,*)'tmps',tmps,'compteur',k,'k1',k1,'k2',k2,'comp_o',comp_o,'\n'
       WRITE(2,*)'comp_n',comp_n,'samp',samp,'power',power,'Dt',Dt
       CLOSE(2)

       WRITE(NF_Name,Format_NF),Root_Code,NF,'_GRANDEUR.txt'
       OPEN(unit =3,file=NF_Name, action ="write",status ="replace")
       WRITE(3,*)'tmps',tmps,'H',H(1,3),'T',T(1,3),'Dt',tmps-tmps_n
       CLOSE(3)

!!! Backup
       IF (compteur == k1) THEN
          WRITE(Output_Name_F,Format_Backup),Output_Racine,'Backup_',compteur,'.dat'
          OPEN(1,file=Output_Name_F,status='replace')
          DO i=1,M,1
             WRITE(1,Format)k,k1,k2,z,compteur,M,tmps,dist(i),ray(i)&
                  &,H(i,1),H(i,2),H(i,3),Xi(i,1),Xi(i,2),Xi(i,3),T(i,1),T(i,2),T(i,3)&
                  &,BL(i,1),BL(i,2),BL(i,3),Ts(i,1),Ts(i,2),Ts(i,3),P(i,1),P(i,2),P(i,3),R_Intrusion
          END DO
          CLOSE(1)
          k1=k1+1
       END IF

!!! DATA FILE
    CASE(1)
       IF (compteur==k2) THEN
          print*,'Write data',tmps,k
          ! Data pour chaque point de la grille
          WRITE(Data_File,Format_RV)Output_Racine,'RV_',compteur,'.dat'
          OPEN(unit=2,file=Data_File,status='replace')
          Format_Data='(54(D30.24,2X))'
          
          R = 0.d0
          DO i =1,M,1
             IF (H(i,3)<=delta0) THEN
                R = dist(i-1);EXIT
             ENDIF
          ENDDO
          
          ! Header
          WRITE(2,'(54(A,2X))')'tm', 'dist', 'H',&
               &'Te','BL','Xi','Ts','P','Srr','Stt',&
               &'hmubar','hthetabar','ubar',&
               &'R','Phi','Vm','Tm','Mum',&
               &'Vm01','Mum01','Vm02','Mum02','Vm05','Mum05',&
               &'Vm005','Mum005',&
               &'BV_a','BV_b','V_t1','V_t2','BE_a','BE_b','En_t1',&
               &'En_t2','Phi_s','Phi_l','Tm01','Tm02','Tm05','Tm005',&
               &'Fr_d_R','Fr_d_T','Fr_d_Mu','Fr_001_R','Fr_001_T','Fr_001_Mu','Mu_e',&
               &'Fr_Mu_R','Fr_Mu_T','Fr_Mu_Mu','Fr_Mu_H',&
               &'Fr_005_R','Fr_005_T','Fr_005_Mu'

          DO i=1,M,1
             IF (H(i,3) == delta0) EXIT
             WRITE(2,Format_Data)tmps, dist(i), H(i,3),&
                  &T(i,3), BL(i,3), Xi(i,3), Ts(i,3), P(i,3), Srr(i), Stt(i),&
                  &hmubar(i),hthetabar(i),ubar(i),&
                  &R, Phim, Vm, Tm, Mum,&
                  &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
                  &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l,&
                  &Tm01, Tm02, Tm05, Tm005,&
                  &Fr_d_R, Fr_d_T, Fr_d_Mu, Fr_001_R, Fr_001_T, Fr_001_Mu, Mu_e,&
                  &Fr_Mu_R,Fr_Mu_T,Fr_Mu_Mu,Fr_Mu_H,&
                  &Fr_005_R,Fr_005_T,Fr_005_Mu
          END DO
          CLOSE(2)
         
          !On augmente le compteut
          k2=k2+1
       END IF

!!! DATA 
    END SELECT

  END SUBROUTINE OUTPUT

END MODULE MODULE_OUTPUT
