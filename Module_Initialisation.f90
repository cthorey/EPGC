MODULE MODULE_INITIALISATION

CONTAINS

  SUBROUTINE  CONSTANTE(M,tmps_m,Dt,Dr,sample,el,grav,delta0,sigma,nu,Pe,psi,N1,&
       &eps_1,Format_O,Format_NSD,Init,Input_Racine,Output_Racine,Input_Data_Name&
       &,Format_NSD_Init_0,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  DECLARATION DES VARIABLES 

    IMPLICIT NONE

    ! Parametre du model
    INTEGER, INTENT(INOUT) :: M,sample,Init
    DOUBLE PRECISION, INTENT(INOUT) :: Dt,Dr,tmps_m

    ! Nombre sans dimension
    DOUBLE PRECISION , INTENT(OUT) :: el,grav,delta0,sigma,nu,Pe,Psi,N1

    ! Seuil de confiance
    DOUBLE PRECISION, INTENT(INOUT) :: eps_1

    ! String variable pour les nom fichiers
    CHARACTER(LEN=300), INTENT(INOUT)  :: Output_Racine,Format_O,Format_NSD
    CHARACTER(LEN=300) , INTENT(INOUT) :: Input_Racine,Input_Data_Name
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_NSD_Init_0,Format_NSD_Init_1
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_Input_Data,Format_RV,Format_Backup
    CHARACTER(LEN=300) :: Output_Name_NSD,Format_Racine
    CHARACTER(LEN=3) :: I_Racine,O_Racine
    CHARACTER(LEN=2) :: I_D_R,NSD

    CHARACTER(LEN=64) :: Root

    ! Parametre du model
    tmps_m = 1D32  
    M = 2000    
    Dt = 1D-6
    Dr = 1D-2
    eps_1 = 1D-4

    ! Nombre sans dimensions
    el = 1D0
    grav = 0D0
    delta0 = 1D-2
    sigma = 2D-2
    nu = 1D-3
    Pe = 1D-2
    psi = 0D0
    N1 = 1D5

    ! Variable pour l'outxsput
    sample = (Dt)/(Dt)
    Init = 1
    Input_Data_Name =  'Backup_0000064.dat'
    Root = '/Users/thorey/Documents/These/Projet/Refroidissement&
         &/Skin_Model/'
    Input_Racine = Root//'Code/Code_ELAS/Test/Run/'
    Output_Racine = Root//'Code/Code_ELAS/Test/Run/'

    ! Ecriture du nom des fichiers
    write(I_Racine,'(I3)'),len(trim(Input_Racine))
    write(O_Racine,'(I3)'),len(trim(Output_Racine))
    write(I_D_R,'(I2)'),len(trim(Input_Data_Name))
    write(NSD,'(I2)'),len(trim('NbSsDim.txt'))

    Format_NSD_Init_1='(a'//I_Racine//',a'//NSD//')'
    Format_Input_Data='(a'//I_Racine//',a'//I_D_R//')'
    Format_NSD_Init_0='(a'//I_Racine//',a'//NSD//')'
    Format_RV='(a'//O_Racine//',a3,i7.7,a4)'
    Format_Backup='(a'//O_Racine//',a7,i7.7,a4)'


    ! Format d'ecritur du fichier de sorti

    Format_O = "(I10,2X,I10,2X,I10,2X,I10,2X,I5,2X,I10,2X,20(D20.14&
         &,2X))"

    ! Format d'ecriture pour le fichier ou seront ecrit les nombres
    ! sans dimensions
    Format_NSD = "(D20.14,2X,D20.14,2X,D20.14,2X,D20.14,2X,D20.14,2X&
         &,D20.14,2X,D20.14,2X,D20.14,2X,I5,2X,D20.14,2X&
         &,D20.14,2X,D20.14,2X)"

    IF (Init==0) THEN
       WRITE(Output_Name_NSD,Format_NSD_Init_0),Output_Racine&
            &,'NbSsDim.txt'
       OPEN(unit=1,file=Output_Name_NSD,status='replace')
       WRITE(1,'(24(A,2X))')'el', 'grav', 'delta0', 'sigma', 'nu', 'Pe','Psi',&
            &'N1','M','Dt','Dr','eps'
       WRITE(1,Format_NSD),el,grav,delta0,sigma,nu,Pe,Psi,N1,M,Dt,Dr&
            &,eps_1
       CLOSE(1)
    END IF

  END SUBROUTINE CONSTANTE

  SUBROUTINE INITIALISATION(Format_O,Format_NSD,M,H,T,Ts,Xi,BL,P,dist&
       &,ray,k,k1,k2,z,tmps, Dt,Dr,eps_1,el,grav,delta0,sigma,nu,Pe&
       &,Psi,N1,sample,Init,compteur,tmps_m, Input_Data_Name&
       &,Input_racine,Output_Racine ,Format_NSD_Init_0&
       &,Format_NSD_Init_1,Format_Input_Data,Format_RV,Format_Backup)

    IMPLICIT NONE

    ! Tableaux
    DOUBLE PRECISION ,DIMENSION(:,:), INTENT(INOUT) :: H,T,Xi,BL,Ts,P
    DOUBLE PRECISION ,DIMENSION(:) , INTENT(INOUT) :: dist,ray

    ! Parametre necessaire au bon deroulement des operation
    INTEGER , INTENT(INOUT) :: k,k1,k2,z,compteur,Init
    INTEGER , INTENT(INOUT) :: sample
    DOUBLE PRECISION, INTENT(INOUT) :: tmps

    ! Nombre sans dimension
    DOUBLE PRECISION , INTENT(INOUT) :: el,grav,delta0,sigma,nu,Pe&
         &,Psi,N1

    ! Parameter du model
    DOUBLE PRECISION , INTENT(INOUT) :: Dt,Dr,eps_1,tmps_m
    INTEGER, INTENT(INOUT) :: M

    ! String variable pour les nom fichiers
    CHARACTER(LEN=300), INTENT(INOUT)  :: Output_Racine,Format_O&
         &,Format_NSD
    CHARACTER(LEN=300) , INTENT(INOUT) :: Input_Racine,Input_Data_Name
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_NSD_Init_0&
         &,Format_NSD_Init_1
    CHARACTER(LEN=300) , INTENT(INOUT) :: Format_Input_Data,Format_RV&
         &,Format_Backup

    !Variable specitif au sous programme
    CHARACTER(LEN=300) :: Input_Name_NSD,Input_Data
    LOGICAL :: FILE_EXISTS
    INTEGER :: i

    SELECT CASE (Init)

    CASE(0)
       CALL CONSTANTE(M,tmps_m,Dt,Dr,sample,el,grav,delta0,sigma,nu&
            &,Pe,Psi,N1, eps_1,Format_O,Format_NSD,Init,Input_Racine&
            &,Output_Racine,Input_Data_Name ,Format_NSD_Init_0&
            &,Format_NSD_Init_1,Format_Input_Data,Format_RV&
            &,Format_Backup)

       compteur = 0 
       dist = 0;ray = 0
       DO i = 1,M,1
          dist(i) = (i-0.5d0)*Dr
          ray(i) = dist(i)+0.5d0*Dr
       END DO
       H = delta0
       T = 0d0
       ! DO i = 1,10
       !    H(i,:) = (dist(10)-dist(i)**2)**2
       !    BL(i,:) = H(i,:)/2.d0
       !    T(i,:)  = H(i,:)/H(1,1)
       ! ENDDO
       Ts = 0.d0
       P = 0.d0
       BL = H/(2.d0) 
       Xi = BL/3.D0*(3.D0-2*T-Ts)
       !Xi = T*H-2.d0/3.d0*T*BL
       k = 0;k1 = 1;k2 = 1;z = 0;tmps = 0

    CASE(1)
       WRITE(Input_Name_NSD,Format_NSD_Init_1)Input_Racine,'NbSsDim.txt'
       print*,Input_Name_NSD
       INQUIRE(FILE = Input_Name_NSD, EXIST = FILE_EXISTS)
       IF (FILE_EXISTS) THEN
          OPEN(1,file = Input_Name_NSD)
          READ(1,*) ! Read the header but don't do nothing with it
          READ(1,Format_NSD),el,grav,delta0,sigma,nu,Pe,Psi,N1,M,Dt,Dr,eps_1
          CLOSE(1)
       ELSEIF( .NOT. FILE_EXISTS) THEN
          PRINT*,'ERREUR: PAS DE FICHIER AVEC LES NOMBRES SANS DIMENSIONS. RECOMMENCER LA SIMU DEPUIS LE DEPART'
          STOP
       ENDIF

       WRITE(Input_Data,Format_Input_Data)Input_Racine,Input_Data_Name
       print*,Input_Data
       INQUIRE(FILE = Input_Data, EXIST = FILE_EXISTS)  
       IF (FILE_EXISTS) THEN
          OPEN(1,file = Input_Data)
          DO i = 1,M,1
             READ(1,Format_O)k,k1,k2,z,compteur,M,tmps,dist(i),ray(i)&
                  &,H(i,1),H(i,2),H(i,3),Xi(i,1),Xi(i,2),Xi(i,3),T(i,1),T(i,2),T(i,3)&
                  &,BL(i,1),BL(i,2),BL(i,3),Ts(i,1),Ts(i,2),Ts(i,3),P(i,1),P(i,2),P(i,3)
          END DO
          CLOSE(1)
          
          k1  =  k1+1
          k2  =  k2+1

       ELSEIF( .NOT. FILE_EXISTS) THEN
          PRINT*,'ERREUR: PAS DE FICHIER INPUT. CHOISIR LE MODE INIT = 0 AND TRY AGAIN'
          STOP
       ENDIF

    END SELECT

  END SUBROUTINE INITIALISATION
     
   END MODULE MODULE_INITIALISATION
      
