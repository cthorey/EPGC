MODULE MODULE_OUTPUT

CONTAINS

  SUBROUTINE  OUTPUT(Format,Dt,M,H,T,Xi,BL,Ts,P,dist,ray,k,k1,k2,z,compteur,tmps,&
          &Output_Racine,delta0,Cas,sample,Format_RV,Format_Backup&
          &,Phim,Vm,Tm,Mum,Srr,Stt,&
          &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
          &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l)

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: H,T,Ts,Xi,BL,P
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: dist,ray,Srr,Stt

    ! Parametre du model
    DOUBLE PRECISION , INTENT(INOUT) :: tmps,delta0,Dt
    DOUBLE PRECISION , INTENT(INOUT) :: Phim,Vm,Tm,Mum
    DOUBLE PRECISION , INTENT(INOUT) :: Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005
    DOUBLE PRECISION ,INTENT(INOUT) :: BE_a,BE_b
    DOUBLE PRECISION ,INTENT(INOUT) :: En_t1,En_t2,Phi_s,Phi_l
    DOUBLE PRECISION ,INTENT(INOUT) :: BV_a,BV_b
    DOUBLE PRECISION ,INTENT(INOUT) :: V_t1,V_t2

    INTEGER , INTENT(INOUT) :: k,k1,k2,z,M,compteur

    ! Format files
    CHARACTER(LEN=300) , INTENT(IN) :: Format,Output_Racine
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_RV,Format_Backup
    LOGICAL :: FILE_EXISTS
    INTEGER , INTENT(IN) :: Cas,sample

    !Variable specitfique au sous programmes
    INTEGER :: i,comp_n,comp_o,samp
    CHARACTER(LEN=300) :: Output_Name_F
    CHARACTER(LEN=300) :: Data_file,Format_Data
    DOUBLE PRECISION :: R,hmax,power


    SELECT CASE (Cas)

!!! BACKUP FILE
    CASE(0)

!!!Compteur
       IF (tmps == 0D0) THEN
          samp = 1
       ELSE
          power = FLOOR(LOG10(tmps))
          samp = 10**power/Dt
       ENDIF

       comp_o = k/samp
       comp_n = (k+1)/samp

       IF (comp_o /= comp_n ) THEN
          compteur = compteur +1
       ENDIF

!!! Backup
       IF (compteur == k1) THEN
          WRITE(Output_Name_F,Format_Backup),Output_Racine,'Backup_',compteur,'.dat'
          OPEN(1,file=Output_Name_F,status='replace')
          DO i=1,M,1
             WRITE(1,Format)k,k1,k2,z,compteur,M,tmps,dist(i),ray(i)&
                  &,H(i,1),H(i,2),H(i,3),Xi(i,1),Xi(i,2),Xi(i,3),T(i,1),T(i,2),T(i,3)&
                  &,BL(i,1),BL(i,2),BL(i,3),Ts(i,1),Ts(i,2),Ts(i,3),P(i,1),P(i,2),P(i,3)
          END DO
          CLOSE(1)
          k1=k1+1
       END IF

!!! DATA FILE
    CASE(1)
       IF (compteur==k2) THEN
          ! Data pour chaque point de la grille
          WRITE(Data_File,Format_RV)Output_Racine,'RV_',compteur,'.dat'
          OPEN(unit=2,file=Data_File,status='replace')
          Format_Data='(33(D30.24,2X))'
          
          R = 0.d0
          DO i=1,M,1
             IF (H(i,3)-delta0>0.d0 .OR. R /= 0.d0) THEN
                CYCLE
             ELSE
                R = dist(i)
             END IF
          ENDDO

          ! Header
          WRITE(2,'(33(A,2X))')'tm', 'dist', 'H',&
               &'Te','BL','Xi','Ts','P','Srr','Stt',&
               &'R','Phi','Vm','Tm','Mum',&
               &'Vm01','Mum01','Vm02','Mum02','Vm05','Mum05',&
               &'Vm005','Mum005',&
               &'BV_a','BV_b','V_t1','V_t2','BE_a','BE_b','En_t1',&
               &'En_t2','Phi_s','Phi_l'

          DO i=1,M,1
             IF (H(i,3) == delta0) EXIT
             WRITE(2,Format_Data)tmps,dist(i),H(i,3)&
                  &,T(i,3),BL(i,3),Xi(i,3),Ts(i,3),P(i,3),Srr(i),Stt(i)&
                  &,R,Phim,Vm,Tm,Mum,&
                  &Vm01,Mum01,Vm02,Mum02,Vm05,Mum05,Vm005,Mum005,&
                  &BV_a,BV_b,V_t1,V_t2,BE_a,BE_b,En_t1,En_t2,Phi_s,Phi_l
          END DO
          CLOSE(2)
         
          !On augmente le compteut
          k2=k2+1
       END IF

!!! DATA 
    END SELECT

  END SUBROUTINE OUTPUT

END MODULE MODULE_OUTPUT
