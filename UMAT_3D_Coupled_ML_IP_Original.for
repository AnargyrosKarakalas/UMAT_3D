!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
!                                                                              !
!     Subroutine for the Implementation of Thermomechanical Constitutive Model !
!     for Shape Memory Alloys                                                  !
!                                                                              !
!     AUTHORS:                                                                 !
!     Dr. Anargyros A. Karakalas                                               !
!                                                                              !
!     Dimension: 3-D                                                           !
!     Version IP Scaling - Original v3.1                                       !
!                                                                              !
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW!
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,
     & NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      !O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O!
      
      !O-O-O-O-O- | UMAT Command Lines  required by Default | -O-O-O-O-O!
      
      CHARACTER*80 CNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     & DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     & SPREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),
     & DFGRD1(3,3),JSTEP(4)
      
      !O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O-O!
      
      !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
      
      ! USER Defined Variables
      
      ! Integer Variables
      INTEGER Elastic_Response, Debug, i
      
      ! Double Precision Variables
      REAL*8 E_A, E_M, PR, alpha_A, alpha_M, MVF_initial,
     & dRPL_dTemperature
      
      ! Double Precision Arrays
      REAL*8 Stress_previous(NTENS,1),
     & S_Iso(NTENS,NTENS), S_sma(NTENS,NTENS), C_sma(NTENS,NTENS),
     & alpha_Coef(NTENS,1), alpha_SMA(NTENS,1), Strain(NTENS,1),
     & dSigma_dEpsilon(NTENS,NTENS), dSigma_dTemperature(NTENS,1),
     & Stress_Current(NTENS,1),
     & dRPL_dEpsilon(1,NTENS),
     & Model_Parameters(7), Initial_State(8), SMA_State(31)
      
C:-----------------------------------------------------------------------!
      
      Elastic_Response = INT(PROPS(33))
      Debug = INT(PROPS(38))
      
C:-----------------------------------------------------------------------!
      
C: Print UMAT Initialization Messages
      IF (Debug == 1) THEN
          ! Write Initialization Messages
          WRITE(*,*)
          WRITE(*,*)"############### BEGINNING of UMAT ###############"
          WRITE(*,*)
          WRITE(*,*)"Step Number - JSTEP(1)",JSTEP(1)
          WRITE(*,*)"Increment Number - KINC",KINC
          WRITE(*,*)"Element Number - NOEL",NOEL
          WRITE(*,*)"Integration Point Number - NPT",NPT
          WRITE(*,*)
      ENDIF
      
C:-----------------------------------------------------------------------!
      
      IF (Elastic_Response == 0) THEN
C: Call the Subroutine for the Calculation of the Model Parameters 
          CALL Model_Parameters_Calculation(STATEV,TEMP,NSTATV,
     & PROPS,NPROPS,JSTEP,KINC,
     & Model_Parameters)
C: Call the Subroutine that Stores the required variables at the Initial State
          CALL Store_Initial_State(STATEV,TEMP,NSTATV,NTENS,PROPS,NPROPS,
     & JSTEP,KINC,Initial_State)
C: Call the Subroutine that describes the State of the Material
          CALL State_Assignment(STRESS,STATEV,NTENS,NSTATV,PROPS,
     & NPROPS,JSTEP,KINC,Initial_State,SMA_State)
C: Call the Subroutine for the Determination of Transformation Direction
          CALL Transformation_Direction(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,KINC,
     & Model_Parameters,Initial_State,SMA_State)
C: Call the Subroutine for Updating the Internal State of the Material
          CALL Update_Internal_State(STATEV,STRAN,DSTRAN,TEMP,
     & DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
C: Call the Subroutine to Update the STRESS Vector required by the Abaqus Global Solver
          CALL FE_STRESS_Calculation(STRESS,
     & STRAN,DSTRAN,TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Initial_State,SMA_State)
C: Call the Subroutine for providing the Partial Derivatives required by the Abaqus Global Solver
          CALL FE_Derivatives(STATEV,DDSDDE,DDSDDT,DRPLDE,
     & DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,NTENS,NSTATV,PROPS,
     & NPROPS,Model_Parameters,Initial_State,SMA_State)
C: Call the Subroutine to Update the RPL required by the Abaqus Global Solver
          CALL FE_RPL_Calculation(STRESS,STATEV,RPL,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
C: Call the Subroutine to Correct the Internal State variables if required
          CALL Internal_State_Correction(NTENS,SMA_State)
C: Call the Subroutine to Update the STATE variables
          CALL Update_STATEVs(STRESS,STATEV,TEMP,DTEMP,
     & NTENS,NSTATV,Initial_State,SMA_State)
          
      ELSE
          
C: Define the Elastic Properties for the Pure Phases
          E_A = PROPS(2)
          E_M = PROPS(3)
          PR = PROPS(4) ! Poisson's Ratio
C: Define the Thermal Expansion Coefficients of the Pure Phases
          alpha_A = PROPS(6)
          alpha_M = PROPS(7)
C: Read the Martensitic Volume Fraction at the Initial State
          MVF_initial = PROPS(31)
C: Define the previous Stress Vector
          DO i=1,NTENS
              Stress_previous(i,1) = STRESS(i)
          ENDDO
C: Define the Compliance Coefficients of a Generally Isotropic Elastic Material [6 x 6]
          S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
          S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
          S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
          S_Iso(4,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) ,0.0_8, 0.0_8 /)
          S_Iso(5,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
          S_Iso(6,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
C: Calculate the Compliance for the given State [6 x 6]
          S_sma = 
     & ((1.0_8/E_A)+(MVF_initial*((1.0_8/E_M)-(1.0_8/E_A))))*S_Iso
C: Calculate the Partial Derivative of Stress with respect to Total Strain [6 x 6]
          dSigma_dEpsilon = C_sma
C: Calculate the Stifness for the given State [6 x 6]
          CALL Matrix_Inversion(S_sma,C_sma,NTENS)
C: Define the Thermal Expansion Matrix Coefficients [6 x 1]
          alpha_Coef(1,1) = 1.0_8
          alpha_Coef(2,1) = 1.0_8
          alpha_Coef(3,1) = 1.0_8
          alpha_Coef(4,1) = 0.0_8
          alpha_Coef(5,1) = 0.0_8
          alpha_Coef(6,1) = 0.0_8
C: Calculate the Thermal Expansion Matrix Coefficients for the given State [6 x 1]
          alpha_SMA = 
     & (alpha_A+(MVF_initial*(alpha_M-alpha_A)))*alpha_Coef
C: Calculate the current Strain Vector [6 x 1]
          DO i=1,NTENS
              Strain(i,1) = STRAN(i)+DSTRAN(i)
          ENDDO
C: Calculate the current Stress Vector [6 x 1]
          Stress_Current = 
     & MATMUL(C_sma,(Strain-(alpha_SMA*(TEMP+DTEMP))))
C: Calculate the Partial Derivative of Stress with respect to Temperature [6 x 1]
          dSigma_dTemperature = -MATMUL(C_sma,alpha_SMA)
C: Calculate the Partial Derivative of the Volumetric Heat Generation per Unit Time with respect to Total Strain [1 x 6]
          dRPL_dEpsilon =
     & -(MATMUL(TRANSPOSE(alpha_SMA)*(TEMP+DTEMP),C_sma)/DTIME)
C: Calculate the Partial Derivative of Heat with respect to the Temperature [SCALAR]
          dRPL_dTemperature = (MINVAL(MATMUL(MATMUL(
     & TRANSPOSE(alpha_SMA),C_sma),alpha_SMA*(TEMP+DTEMP)))/DTIME)-
     & (MINVAL(MATMUL(TRANSPOSE(alpha_SMA),
     & (Stress_current-Stress_previous)))/DTIME)
C: Define the STRESS quantity required by the Abaqus Global Solver
          DO i=1,NTENS
              STRESS(i) = Stress_current(i,1)
          ENDDO
C: Define the DDSDDE quantity required by the Abaqus Global Solver          
          DDSDDE = dSigma_dEpsilon
C: Define the DDSDDT quantity required by the Abaqus Global Solver
          DO i=1,NTENS
              DDSDDT(i) = dSigma_dTemperature(i,1)
          ENDDO
C: Define the DRPLDE quantity required by the Abaqus Global Solver
          DO i=1,NTENS
              DRPLDE(i) = dRPL_dEpsilon(1,i)
          ENDDO
C: Define the DRPLDT quantity required by the Abaqus Global Solver
          DRPLDT = dRPL_dTemperature
      ENDIF
      
C:-----------------------------------------------------------------------!
      
      RETURN
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
C: Call this Subroutine to Calculate the required Model Parameters
      SUBROUTINE Model_Parameters_Calculation(STATEV,TEMP,NSTATV,
     & PROPS,NPROPS,JSTEP,KINC,
     & Model_Parameters)
      
C------------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NSTATV, NPROPS, KINC
      
C: Integer Arrays
      INTEGER*4 JSTEP(4)
      
C: Double Precision Variables
      REAL*8 TEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 CA_eq_CM, Temperature_Smoothing, Debug, i
      
C: Double Precision Variables
      REAL*8 Density, E_A, E_M, T_Calibration, alpha_A, alpha_M,
     & Sigma_Calibration, C_A, C_M, ESH_A, ESH_M,
     & Mf_input, Ms_input, As_input, Af_input, H_sat, k_t, Backstress_11,
     & n1, n2, n3, n4, Mf_0, Ms_0, As_0, Af_0, T_initial, S_A, S_M,
     & Delta_S, Delta_a, Delta_c,
     & Mf_smooth, Ms_smooth, As_smooth, Af_smooth,
     & Sigma_Eff_Calibration,
     & H_cur_calibration, dHcurrent_dSigma_calibration,
     & Lambda_calibration, dLambda_dSigma_calibration,
     & Sigma_0, Sigma_effective_0, H_current_0,
     & rDs0, D, a1, a2, a3, rDu0, Y0_t
      
C: Double Precision Arrays
C: Output
      REAL*8 Model_Parameters(7)
      
C---------------------------------------------------------------------------|
! Assign the values to the required Variables based on the Input Properties !
C---------------------------------------------------------------------------|
      
C: Define the Density of the Material
      Density = PROPS(1)
      
C: Define the Elastic Properties for the Pure Phases
      E_A = PROPS(2)
      E_M = PROPS(3)
      
C: Define the Calibration Temperature Level
      T_Calibration = PROPS(5)
      
C: Define the Thermal Expansion Coefficients of the Pure Phases
      alpha_A = PROPS(6)
      alpha_M = PROPS(7)
      
C: Define the Calibration Stress Level for C_A and C_M
      Sigma_Calibration = PROPS(8)
      
C: Define the Inclination of Transformation Boundary Surfaces
      C_A = PROPS(9)
      C_M = PROPS(10)
      
C: Define the Effective Specific Heat - ESH of the Pure Phases
      ESH_A = PROPS(11)
      ESH_M = PROPS(12)
      
C: Define the Transformation Temperature Thresholds at Zero Stress Level
      Mf_input = PROPS(13)
      Ms_input = PROPS(14)
      As_input = PROPS(15)
      Af_input = PROPS(16)
      
C: Define the Maximum Transformation Strain Function Parameters
      H_sat = PROPS(17) ! Maximum Saturated Transformation Strain
      k_t = PROPS(18) ! Exponent of Maximum Transformation Strain Function
      
C: Define the Backstress calibrated value at XX Direction
      Backstress_11 = PROPS(19)
      IF (Backstress_11 == 0.0_8) THEN
C: Use a very small value in order to avoid any singularities in the Lamda Tensor when STRAN=DSTRAN=0.0
          Backstress_11 = 1.0e-24
      ENDIF
      
C: Define the Transformation Hardening Function Parameters
      n1 = PROPS(25)
      n2 = PROPS(26)
      n3 = PROPS(27)
      n4 = PROPS(28)
      
C: Read the Index for the Definition/Calculation of the Stress Influence Coefficients
C: If CA_eq_CM = 0 then C_A and C_M obtain the User Input values
C: If CA_eq_CM = 1 then C_A equals C_M which equals the User Input value
C: If CA_eq_CM = 1 then C_A equals C_M which equals the average value of the Input values of C_A and C_M
      CA_eq_CM = INT(PROPS(34))
      
C: Read the value of the Index regarding the Use of Smoothing Temperature Thresholds
C: If Temperature_Smoothing = 0 then the User Input values for the Temperature Thresholds at Zero Stress Level are Considered
C: If Temperature_Smoothing = 1 then the Smoothed values for the Temperature Thresholds at Zero Stress Level are Considered
      Temperature_Smoothing = INT(PROPS(35))
      
C: Read the value of Debug variable
C: If Debug = 0 WRITE statements are deactivated
C: If Debug = 1 WRITE statements are activated
      Debug = INT(PROPS(38))
      
C------------------------------------------------------------------------------------|
! Assign the values to the required Variables based on the Step and Increment Number !
C------------------------------------------------------------------------------------|
      
      IF (JSTEP(1) == 1 .AND. KINC <= 1) THEN
C: Set the value for the Initial Temperature
          T_initial = TEMP
      ELSE
C: Set the value for the Initial Temperature
          T_initial = STATEV(31)
      ENDIF
      
C-----------------------------------------------------------------|
! Perform the required Calculations based on the Input Properties !
C-----------------------------------------------------------------|
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
      
C: Calculate the Difference between Martensite and Austenite Properties
      
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
      
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
      
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C: Based on the Stress Influence Coefficients Index Define their values
      IF (CA_eq_CM == 0) THEN
          C_A = PROPS(9)
          C_M = PROPS(10)
      ELSEIF (CA_eq_CM == 1) THEN
          C_A = C_M
      ELSEIF (CA_eq_CM == 2) THEN
          C_A = (C_A+C_M)/2.0_8
          C_M = (C_A+C_M)/2.0_8
      ENDIF
      
C: Based on the Temperature Smoothing Index Define the Transformation Temperature Thresholds
      IF (Temperature_Smoothing == 1) THEN
C: Calculate the Modified Transformation Temperature Thresholds considering Smooth Hardening
          Ms_smooth = 
     &((Ms_input/2.0_8)*(1.0_8+((2.0_8**(-n1))*(n1+1.0_8))+((2.0_8**
     &(-n2))*(n2-1.0_8)))/((n1*(2.0_8**(-n1)))+(n2*(2.0_8**(-n2)))))+
     &((Mf_input/2.0_8)*(-1.0_8+((2.0_8**(-n1))*(n1-1.0_8))+((2.0_8**
     &(-n2))*(n2+1.0_8)))/((n1*(2.0_8**(-n1)))+(n2*(2.0_8**(-n2)))))
          Mf_smooth = 
     &((Ms_input/2.0_8)*(-1.0_8+((2.0_8**(-n1))*(n1+1.0_8))+((2.0_8**
     &(-n2))*(n2-1.0_8)))/((n1*(2.0_8**(-n1)))+(n2*(2.0_8**(-n2)))))+
     &((Mf_input/2.0_8)*(1.0_8+((2.0_8**(-n1))*(n1-1.0_8))+((2.0_8**
     &(-n2))*(n2+1.0_8)))/((n1*(2.0_8**(-n1)))+(n2*(2.0_8**(-n2)))))
          As_smooth = 
     &((As_input/2.0_8)*(1.0_8+((2.0_8**(-n3))*(n3-1.0_8))+((2.0_8**
     &(-n4))*(n4+1.0_8)))/((n3*(2.0_8**(-n3)))+(n4*(2.0_8**(-n4)))))+
     &((Af_input/2.0_8)*(-1.0_8+((2.0_8**(-n3))*(n3+1.0_8))+((2.0_8**
     &(-n4))*(n4-1.0_8)))/((n3*(2.0_8**(-n3)))+(n4*(2.0_8**(-n4)))))
          Af_smooth = 
     &((As_input/2.0_8)*(-1.0_8+((2.0_8**(-n3))*(n3-1.0_8))+((2.0_8**
     &(-n4))*(n4+1.0_8)))/((n3*(2.0_8**(-n3)))+(n4*(2.0_8**(-n4)))))+
     &((Af_input/2.0_8)*(1.0_8+((2.0_8**(-n3))*(n3+1.0_8))+((2.0_8**
     &(-n4))*(n4-1.0_8)))/((n3*(2.0_8**(-n3)))+(n4*(2.0_8**(-n4)))))
C: Source:
C: Use of a Ni60Ti Shape Memory Alloy for Active Jet Engine Chevron Application: II. Experimentally Validated Numerical Analysis
C: (Hartl,Mooney,Lagoudas,Calkins,Mabe)
          Mf_0 = Mf_smooth
          Ms_0 = Ms_smooth
          As_0 = As_smooth
          Af_0 = Af_smooth
      ELSE
          Mf_0 = Mf_input
          Ms_0 = Ms_input
          As_0 = As_input
          Af_0 = Af_input
      ENDIF
      
C-----------------------------------------|
C Calculate the required Model Parameters !
C-----------------------------------------|
      
C: For the Calibration Process the 1-D version of the Model is Used - X-Direction
      
C: Calculate the Effective Stress at the Calibration Stress Level
      Sigma_Eff_Calibration = Sigma_Calibration+Backstress_11
C: Calculate the Maximum Transformation Strain at the Calibration Stress Level
      H_cur_calibration = H_sat*(1.0_8-EXP(-k_t*Sigma_Eff_Calibration))
C: Calculate the Derivative of the Maximum Transformation Strain Function at the Calibration Stress Level
      dHcurrent_dSigma_calibration = 
     & H_sat*k_t*EXP(-k_t*Sigma_Eff_Calibration)
C: Calculate the Transformation Tensor
      Lambda_calibration = H_cur_calibration
C: Calculate the Derivative of the Transformation Tensor with respect to the Stress
      dLambda_dSigma_calibration = dHcurrent_dSigma_calibration
      
C: Calculate Model parameter: rDs0
      rDs0 = ((-2.0_8)*(C_A*C_M)*((Sigma_Calibration*Delta_S)+
     & (Delta_a*(T_Calibration-T_initial))+Lambda_calibration+
     & (Sigma_Eff_Calibration*dLambda_dSigma_calibration))/
     & (C_M+C_A))-((Sigma_Calibration*Delta_a)+
     & (Density*Delta_c*LOG(T_Calibration/T_initial)))
C: Calculate Model parameter: D
      D = (C_M-C_A)*(
     & (Sigma_Calibration*Delta_S)+(Delta_a*(T_Calibration-T_initial))+
     & Lambda_calibration+(Sigma_Eff_Calibration*
     & dLambda_dSigma_calibration))/((C_M+C_A)*(Lambda_calibration+(
     & Sigma_Eff_Calibration*dLambda_dSigma_calibration)))
C: Calculate Model parameter: alpha 1
      a1 = (rDs0*(Mf_0-Ms_0))-
     & (Density*Delta_c*((Ms_0-T_initial)-(Ms_0*LOG(Ms_0/T_initial))))+
     & (Density*Delta_c*((Mf_0-T_initial)-(Mf_0*LOG(Mf_0/T_initial))))
C: Calculate Model parameter: alpha 2
      a2 = rDs0*(As_0-Af_0)+(Density*Delta_c*((As_0-T_initial)-
     & (As_0*LOG(As_0/T_initial))))-(Density*Delta_c*((Af_0-T_initial)-
     & (Af_0*LOG(Af_0/T_initial))))
C: Calculate Model parameter: alpha 3
      a3 = ((-a1/4.0_8)*(1.0_8+(1.0_8/(n1+1.0_8))-(1.0_8/(n2+1.0_8))))+
     & ((a2/4.0_8)*(1.0_8+(1.0_8/(n3+1.0_8))-(1.0_8/(n4+1.0_8))))
      
C: Define Zero Stress Level value
      Sigma_0 = 0.0_8
C: Calculate the Effective Stress at Zero Stress Level
      Sigma_effective_0 = (2.0_8/3.0_8)*(Sigma_0+Backstress_11)
C: Calculate the Maximum Transformation Strain at the Zero Stress Level
      H_current_0 = H_sat*(1.0_8-EXP(-k_t*Sigma_effective_0))
      
C: Calculate Model parameter: rDu0
      rDu0 = ((rDs0/2.0_8)*(Ms_0+Af_0))+(Backstress_11*H_current_0)+
     &(((Density*Delta_c*((Ms_0-T_initial)-(Ms_0*LOG(Ms_0/T_initial))))-
     &(Density*Delta_c*((Af_0-T_initial)-(Af_0*LOG(Af_0/T_initial)))))/
     &2.0_8)
C: Calculate Model parameter: Y0
      Y0_t = ((rDs0/2.0_8)*(Ms_0-Af_0))-a3-(D*Backstress_11*H_current_0)+
     &(((Density*Delta_c*((Ms_0-T_initial)-(Ms_0*LOG(Ms_0/T_initial))))+
     &(Density*Delta_c*((Af_0-T_initial)-(Af_0*LOG(Af_0/T_initial)))))/
     &2.0_8)
      
C: Source:
C: Constitutive Model for the Numerical Analysis of Phase Transformation in Polycrystalline Shape Memory Alloys
C: (Lagoudas,Hartl,Chemisky,Machado,Popov)
      
C--------------------------------------------------------------------------------------|
C Assign the values of the Model Parameters to the Components of the Respective Matrix !
C--------------------------------------------------------------------------------------|
      
      Model_Parameters(1) = rDs0
      Model_Parameters(2) = D
      Model_Parameters(3) = a1
      Model_Parameters(4) = a2
      Model_Parameters(5) = a3
      Model_Parameters(6) = rDu0
      Model_Parameters(7) = Y0_t
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (JSTEP(1) == 1 .AND. KINC <= 1) THEN
      !    IF (Debug == 1) THEN
      !        ! Write the values of Model Parameters during First Increment of the First Step
      !        WRITE(*,*)
      !        WRITE(*,*)"Mf_0",Mf_0
      !        WRITE(*,*)"Ms_0",Ms_0
      !        WRITE(*,*)"As_0",As_0
      !        WRITE(*,*)"Af_0",Af_0
      !        WRITE(*,*)"############# Model  Parameters #############"
      !        WRITE(*,*)"rDs0",rDs0
      !        WRITE(*,*)"D",D
      !        WRITE(*,*)"a1",a1
      !        WRITE(*,*)"a2",a2
      !        WRITE(*,*)"a3",a3
      !        WRITE(*,*)"rDu0",rDu0
      !        WRITE(*,*)"Y0_t",Y0_t
      !        WRITE(*,*)"#############################################"
      !        WRITE(*,*)
      !    ENDIF
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Store_Initial_State(STATEV,TEMP,NSTATV,NTENS,PROPS,NPROPS,
     & JSTEP,KINC,Initial_State)
      
C------------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NSTATV, NTENS, NPROPS, KINC
      
C: Integer Arrays
      INTEGER*4 JSTEP(4)
      
C: Double Precision Variables
      REAL*8 TEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 et_0_direction, Debug, i
      
C: Double Precision Variables
      REAL*8 H_sat, MVF_initial, T_initial
      
C: Double Precision Arrays
      REAL*8 et_initial(NTENS,1),
C: Output
     & Initial_State(8)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
C: Read the Manimum Transformation Strain Saturation value
      H_sat = PROPS(17)
C: Read the Martensitic Volume Fraction at the Initial State
      MVF_initial = PROPS(31)
C: Set the Initial value for the Martensitic Volume Fraction
      IF (MVF_initial <= 0.0001_8) THEN
          MVF_initial = 0.0_8
      ELSEIF (MVF_initial >= 0.9999_8) THEN
          MVF_initial = 1.0_8
      ENDIF
C: Read the Transformation Strain Direction at the Initialization of the Analysis
      et_0_direction = INT(PROPS(32))
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C: Calculate the Components of the Initial Transformation Strain Vector
      IF (MVF_initial > 0.0001_8) THEN
          IF (ABS(et_0_direction) == 1) THEN
              et_initial(1,1) = MVF_initial*H_sat
              et_initial(2,1) = -MVF_initial*H_sat/2.0_8
              et_initial(3,1) = -MVF_initial*H_sat/2.0_8
              et_initial(4,1) = 0.0_8
              et_initial(5,1) = 0.0_8
              et_initial(6,1) = 0.0_8
              et_initial = et_initial*SIGN(1,et_0_direction)
          ELSEIF (ABS(et_0_direction) == 2) THEN
              et_initial(1,1) = -MVF_initial*H_sat/2.0_8
              et_initial(2,1) = MVF_initial*H_sat
              et_initial(3,1) = -MVF_initial*H_sat/2.0_8
              et_initial(4,1) = 0.0_8
              et_initial(5,1) = 0.0_8
              et_initial(6,1) = 0.0_8
              et_initial = et_initial*SIGN(1,et_0_direction)
          ELSEIF (ABS(et_0_direction) == 3) THEN
              et_initial(1,1) = -MVF_initial*H_sat/2.0_8
              et_initial(2,1) = -MVF_initial*H_sat/2.0_8
              et_initial(3,1) = MVF_initial*H_sat
              et_initial(4,1) = 0.0_8
              et_initial(5,1) = 0.0_8
              et_initial(6,1) = 0.0_8
              et_initial = et_initial*SIGN(1,et_0_direction)
          ELSEIF (ABS(et_0_direction) == 4) THEN
              et_initial(1,1) = 0.0_8
              et_initial(2,1) = 0.0_8
              et_initial(3,1) = 0.0_8
              et_initial(4,1) = MVF_initial*H_sat*SQRT(3.0_8)
              et_initial(5,1) = 0.0_8
              et_initial(6,1) = 0.0_8
              et_initial = et_initial*SIGN(1,et_0_direction)
          ELSEIF (ABS(et_0_direction) == 5) THEN
              et_initial(1,1) = 0.0_8
              et_initial(2,1) = 0.0_8
              et_initial(3,1) = 0.0_8
              et_initial(4,1) = 0.0_8
              et_initial(5,1) = MVF_initial*H_sat*SQRT(3.0_8)
              et_initial(6,1) = 0.0_8
              et_initial = et_initial*SIGN(1,et_0_direction)
          ELSEIF (ABS(et_0_direction) == 6) THEN
              et_initial(1,1) = 0.0_8
              et_initial(2,1) = 0.0_8
              et_initial(3,1) = 0.0_8
              et_initial(4,1) = 0.0_8
              et_initial(5,1) = 0.0_8
              et_initial(6,1) = MVF_initial*H_sat*SQRT(3.0_8)
              et_initial = et_initial*SIGN(1,et_0_direction)
          ENDIF
      ELSE
          DO i=1,NTENS
              et_initial(i,1) = 0.0_8
          ENDDO
      ENDIF
      
      IF (JSTEP(1) == 1 .AND. KINC <= 1) THEN
C: Set the value for the Initial Temperature
          T_initial = TEMP
      ELSE
C: Set the value for the Initial Temperature
          T_initial = STATEV(31)
      ENDIF
      
C: Define the Initial value for the Martensitic Volume Fraction that will be used in the required Calculations
      IF (MVF_initial < 0.0001_8) THEN
          MVF_initial = 0.0001_8
      ELSEIF (MVF_initial > 0.9999_8) THEN
          MVF_initial = 0.9999_8
      ENDIF
      
      
C---------------------------------------------------------------------------------------------|
C Assign the values of the Initial State Variables to the Components of the Respective Matrix |
C---------------------------------------------------------------------------------------------|
      
      Initial_State(1) = T_initial
      Initial_State(2) = MVF_initial
      Initial_State(3) = et_initial(1,1)
      Initial_State(4) = et_initial(2,1)
      Initial_State(5) = et_initial(3,1)
      Initial_State(6) = et_initial(4,1)
      Initial_State(7) = et_initial(5,1)
      Initial_State(8) = et_initial(6,1)
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Initial_State Components"
      !    DO i=1,8
      !        WRITE(*,*)Initial_State(i)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE State_Assignment(STRESS,STATEV,NTENS,NSTATV,PROPS,
     & NPROPS,JSTEP,KINC,Initial_State,SMA_State)
      
C------------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NSTATV, NTENS, NPROPS, KINC
      
C: Integer Arrays
      INTEGER*4 JSTEP(4)
      
C: Double Precision Variables
      
C: Double Precision Arrays
      REAL*8 STRESS(NTENS), STATEV(NSTATV), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 Debug, i, TO_index, Direction_FWD, Direction_REV,
     & Completion_FWD, Completion_REV,
     & Reversal_FWD_REV, Reversal_REV_FWD,
     & MVF_correction_index
      
C: Double Precision Variables
      REAL*8 T_initial, MVF_initial, MVF_previous, MVF_current,
     & MVF_reversal_FWD_REV, MVF_reversal_REV_FWD, Phi_FWD, Phi_REV
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), et_initial(NTENS,1),
     & et_previous(NTENS,1), et_current(NTENS,1),
     & Stress_previous(NTENS,1), Lambda_REV(NTENS,1),
     & SMA_State(31)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C-----------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State Data |
C-----------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      MVF_initial = Initial_State(2)
      DO i=1,NTENS
          et_initial(i,1) = Initial_State(2+i)
      ENDDO
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Increment Number |
C---------------------------------------------------------------------------|
      
      IF (JSTEP(1) == 1 .AND. KINC <= 1) THEN
          
C: Set the value of the Martensitic Volume Fraction during the previous Increment
          MVF_previous = MVF_initial
C: Set the value of the Martensitic Volume Fraction for the current Increment
          MVF_current = MVF_initial
C: Set the values of the Transformation Strain Vector during the previous Increment
          et_previous = et_initial
C: Set the values of the Transformation Strain Vector for the current Increment
          et_current = et_initial
C: Set the values of the Stress Vector during the previous Increment
          DO i=1,NTENS
              Stress_previous(i,1) = STRESS(i)
          ENDDO
C: Set the Transformation Occurence Index equal to Zero
          TO_index = 0
C: Set the initial value for the Forward Transformation Direction Index
          Direction_FWD = 0
C: Set the initial value for the Reverse Transformation Direction Index
          Direction_REV = 0
          
          IF(MVF_current <= 0.0001_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 0
C: Set the Index for the Reverse Transformation Completion
              Completion_REV = 1
          ELSEIF(MVF_current >= 0.9999_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 1
              ! Reverse Transformation has not been Completed
              Completion_REV = 0
          ELSEIF(MVF_current > 0.0001_8 .OR. MVF_initial < 0.9999_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 0
C: Set the Index for the Reverse Transformation Completion
              Completion_REV = 0
          ENDIF
          
C: Set the initial value for the Forward to Reverse Transformation Reversal Index
          Reversal_FWD_REV = 0
C: Set the initial value for the Reverse to Forward Transformation Reversal Index
          Reversal_REV_FWD = 0
C: The value of the MVF at the Last Transformation Reversal Point from Forward to Reverse Transformation
          MVF_reversal_FWD_REV = 0.9999_8
C: The value of the MVF at the Last Transformation Reversal Point from Reverse to Forward Transformation
          MVF_reversal_REV_FWD = 0.0001_8
C: Set the value of the Transformation Function of the Forward Transformation at the end of the previous Increment
          Phi_FWD = 0.0_8
C: Set the value of the Transformation Function of the Reverse Transformation at the end of the previous Increment
          Phi_REV = 0.0_8
C: Set the value for the MVF Correction Index
          MVF_correction_index = 0
C: Set the value of the Transformation Tensor for the Reverse Transformation
          IF (MVF_current <= 0.0001_8) THEN
              Lambda_REV = 0.0_8
          ELSE
              Lambda_REV = et_current/MVF_current
          ENDIF
          
      ELSE
          
C: Set the value of the Martensitic Volume Fraction during the previous Increment
          MVF_previous = STATEV(1)
C: Set the value of the Martensitic Volume Fraction for the current Increment
          MVF_current = STATEV(1)
C: Set the values of the Transformation Strain Vector during the previous Increment
          DO i=1,NTENS
              et_previous(i,1) = STATEV(1+i)
          ENDDO
C: Set the values of the Transformation Strain Vector for the current Increment
          DO i=1,NTENS
              et_current(i,1) = STATEV(1+i)
          ENDDO
C: Set the values of the Stress Vector during the previous Increment
          DO i=1,NTENS
              Stress_previous(i,1) = STRESS(i)
          ENDDO
C: Set the Transformation Occurence Index equal to Zero
          TO_index = 0
C: Set the initial value for the Forward Transformation Direction Index
          Direction_FWD = INT(STATEV(15))
C: Set the initial value for the Reverse Transformation Direction Index
          Direction_REV = INT(STATEV(16))
          
          IF(MVF_current <= 0.0001_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 0
C: Set the Index for the Reverse Transformation Completion
              Completion_REV = 1
          ELSEIF(MVF_current >= 0.9999_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 1
              ! Reverse Transformation has not been Completed
              Completion_REV = 0
          ELSEIF(MVF_current > 0.0001_8 .OR. MVF_initial < 0.9999_8)THEN
C: Set the Index for the Forward Transformation Completion
              Completion_FWD = 0
C: Set the Index for the Reverse Transformation Completion
              Completion_REV = 0
          ENDIF
          
C: Set the initial value for the Forward to Reverse Transformation Reversal Index
          Reversal_FWD_REV = INT(STATEV(19))
C: Set the initial value for the Reverse to Forward Transformation Reversal Index
          Reversal_REV_FWD = INT(STATEV(20))
C: The value of the MVF at the Last Transformation Reversal Point from Forward to Reverse Transformation
          MVF_reversal_FWD_REV = STATEV(21)
C: The value of the MVF at the Last Transformation Reversal Point from Reverse to Forward Transformation
          MVF_reversal_REV_FWD = STATEV(22)
C: Set the value of the Transformation Function of the Forward Transformation at the end of the previous Increment
          Phi_FWD = STATEV(29)
C: Set the value of the Transformation Function of the Reverse Transformation at the end of the previous Increment
          Phi_REV = STATEV(30)
C: Set the value for the MVF Correction Index
          MVF_correction_index = 0
C: Set the value of the Transformation Tensor for the Reverse Transformation
          IF (MVF_current <= 0.0001_8) THEN
              Lambda_REV = 0.0_8
          ELSE
              DO i=1,NTENS
                  Lambda_REV(i,1) = STATEV(22+i)
              ENDDO
          ENDIF
          
      ENDIF
      
C-----------------------------------------------------------------------------------------------------------------------|
C Assign the values of the Variables that describe the State of the Material to the Components of the Respective Matrix |
C-----------------------------------------------------------------------------------------------------------------------|

      SMA_State(1) = MVF_current
      DO i=1,NTENS
              SMA_State(1+i) = et_current(i,1)
      ENDDO
      DO i=1,NTENS
              SMA_State(7+i) = Stress_previous(i,1)
      ENDDO
      SMA_State(14) = TO_index
      SMA_State(15) = Direction_FWD
      SMA_State(16) = Direction_REV
      SMA_State(17) = Completion_FWD
      SMA_State(18) = Completion_REV
      SMA_State(19) = Reversal_FWD_REV
      SMA_State(20) = Reversal_REV_FWD
      SMA_State(21) = MVF_reversal_FWD_REV
      SMA_State(22) = MVF_reversal_REV_FWD
      DO i=1,NTENS
          SMA_State(22+i) = Lambda_REV(i,1)
      ENDDO
      SMA_State(29) = Phi_FWD
      SMA_State(30) = Phi_REV
      SMA_State(31) = MVF_correction_index
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"SMA_State Components"
      !    DO i=1,30
      !        WRITE(*,*)i,SMA_State(i)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
      
C------------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS, NPROPS
      
C: Double Precision Arrays
      REAL*8 PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 E_A, E_M, PR, alpha_A, alpha_M, ESH_A, ESH_M, S_A, S_M,
     & Delta_S, Delta_a, Delta_c,
     & MVF_current, 
     & S_value, alpha_value, c_current
      
C: Double Precision Arrays
      REAL*8 SMA_State(31),
     & S_Iso(NTENS,NTENS), alpha_Iso(NTENS,1),
     & S_current(NTENS,NTENS), a_current(NTENS,1)
      
11    FORMAT(6(E24.8E3))
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4) ! Poisson's Ratio
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
      
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      
C:   |------------------------------------------------------------------|
      
C: Define the Compliance Coefficients of a Generally Isotropic Elastic Material [6 x 6]
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Calculate the Compliance for the Current Increment [6 x 6]
      S_value = (S_A+(MVF_current*Delta_S))
      S_current = S_value*S_Iso
      
C: Define the Thermal Expansion Matrix Coefficients for a Generally Isotropic Material [6 x 1]
      alpha_Iso(1,1) = 1.0_8
      alpha_Iso(2,1) = 1.0_8
      alpha_Iso(3,1) = 1.0_8
      alpha_Iso(4,1) = 0.0_8
      alpha_Iso(5,1) = 0.0_8
      alpha_Iso(6,1) = 0.0_8
      
C: Calculate the Thermal Expansion Coefficient for the Current Increment [6 x 1]
      alpha_value = (alpha_A+(MVF_current*Delta_a))
      a_current = alpha_value*alpha_Iso
      
C: Calculate the Value of the Effective Specific Heat for the Current Increment
      c_current = (ESH_A+(MVF_current*Delta_c))
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"S_current"
      !    DO i=1,NTENS
      !        WRITE(*,FMT=11)S_current(i,:)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
      
C------------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS, NPROPS
      
C: Double Precision
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Double Precision Variables
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 E_A, E_M, PR, alpha_A, alpha_M,ESH_A, ESH_M, T_initial,
     & MVF_current, c_current
      
C: Double Precision Arrays
      REAL*8 SMA_State(31), Initial_State(8), et_initial(NTENS,1),
     & et_current(NTENS,1),Strain(NTENS,1),
     & S_current(NTENS,NTENS), a_current(NTENS,1), S_Inverted(NTENS,NTENS),
C: Output
     & Stress_Current(NTENS,1)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4) ! Poisson's Ratio
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State of the Material |
C----------------------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      DO i=1,NTENS
          et_initial(i,1) = Initial_State(2+i)
      ENDDO
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      
C-------------------------------------------------------------|
C Calculate the Total Strain provided by Abaqus Global Solver |
C-------------------------------------------------------------|
      
      Do i=1,NTENS
          Strain(i,1) = STRAN(i)+DSTRAN(i)
      ENDDO
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Material Properties
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
C: Invert the Compliance to Calculate the Stiffness Matrix
      CALL Matrix_Inversion(S_current,S_Inverted,NTENS)
      
C-----------------------------------------------|
C Calculate the Stress based on Input variables |
C-----------------------------------------------|
      
      Stress_Current = MATMUL(S_Inverted,
     & (Strain-(a_current*(TEMP+DTEMP-T_initial))-
     & (et_current-et_initial)))
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    DO i=1,NTENS
      !        WRITE(*,*)"Current Stress Components"
      !        WRITE(*,*)Stress_Current(i,1)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 MVF_current,
C: Output
     & Sigma_Eff_VM
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31), BackStress(NTENS,1),
     & et_current(NTENS,1), I_6x6(NTENS,NTENS), Stress_Current(NTENS,1),
     & Sigma_Eff_Mean(NTENS,1), dSigma_Eff_VM_dSigma_Eff(1,NTENS),
     & dSigma_Eff_Dev_dSigma_Eff(NTENS,NTENS),
     & dSigma_Eff_Dev_dSigma(NTENS,NTENS),
     & dSigma_Eff_dSigma(NTENS,NTENS),
C: Output
     & Sigma_Eff(NTENS,1), Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS)
      
11    FORMAT(6(E24.8E3))
      
C----------------------------------------------------------------------|
C Define any Auxiliary Matrices required for the particular Subroutine |
C----------------------------------------------------------------------|
      
C: Define an Identity Matrix - [6 x 6]
      I_6x6(1,:) = (/ 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(2,:) = (/ 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(3,:) = (/ 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8 /)
      I_6x6(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8 /)
      I_6x6(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8 /)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
C: Define the Components of the BackStress Matrix [6 x 1]
      DO i=1,NTENS
          IF (PROPS(18+i) == 0.0_8) THEN
              BackStress(i,1) = 1.0e-24
          ELSE
              BackStress(i,1) = PROPS(18+i)
          ENDIF
      ENDDO
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
      
C-----------------------------------------------|
C Calculate the Effective Stress Vector [6 x 1] |
C-----------------------------------------------|
      Sigma_Eff = Stress_Current + BackStress
      
C-------------------------------------------------------------------------|
C Calculate the Effective Mean Stress / Hydrostatic Stress Vector [6 x 1] |
C-------------------------------------------------------------------------|
      Sigma_Eff_Mean(1,1) =
     & ((Sigma_Eff(1,1)+Sigma_Eff(2,1)+Sigma_Eff(3,1))/3.0_8)
      Sigma_Eff_Mean(2,1) =
     & ((Sigma_Eff(1,1)+Sigma_Eff(2,1)+Sigma_Eff(3,1))/3.0_8)
      Sigma_Eff_Mean(3,1) =
     & ((Sigma_Eff(1,1)+Sigma_Eff(2,1)+Sigma_Eff(3,1))/3.0_8)
      Sigma_Eff_Mean(4,1) = 0.0_8
      Sigma_Eff_Mean(5,1) = 0.0_8
      Sigma_Eff_Mean(6,1) = 0.0_8
      
C----------------------------------------------------------|
C Calculate the Effective Deviatoric Stress Vector [6 x 1] |
C----------------------------------------------------------|
      Sigma_Eff_Dev = Sigma_Eff - Sigma_Eff_Mean
      
C--------------------------------------------------------------|
C Calculate the Effective Von Mises Equivalent Stress [SCALAR] |
C--------------------------------------------------------------|
      Sigma_Eff_VM = SQRT(((1.0_8/2.0_8)*(
     & ((Sigma_Eff(1,1)-Sigma_Eff(2,1))**2.0_8)+
     & ((Sigma_Eff(2,1)-Sigma_Eff(3,1))**2.0_8)+
     & ((Sigma_Eff(3,1)-Sigma_Eff(1,1))**2.0_8)))+
     & (3.0_8*
     & ((Sigma_Eff(4,1)**2.0_8)+
     & (Sigma_Eff(5,1)**2.0_8)+
     & (Sigma_Eff(6,1)**2.0_8))))
      
C:-----------------------------------------------------------------------!
C:                  Calculation of Partial Derivatives                  :| 
C:-----------------------------------------------------------------------!
      
C: Calculate the Partial Derivative of the Effective Von Misses Stress with respect to the Effective Stress [1 x 6]
      dSigma_Eff_VM_dSigma_Eff(1,1) = (1.0_8/(2.0_8*Sigma_Eff_VM))*(
     & (+2.0_8*Sigma_Eff(1,1))+
     & (-1.0_8*Sigma_Eff(2,1))+
     & (-1.0_8*Sigma_Eff(3,1)))
      dSigma_Eff_VM_dSigma_Eff(1,2) = (1.0_8/(2.0_8*Sigma_Eff_VM))*(
     & (-1.0_8*Sigma_Eff(1,1))+
     & (+2.0_8*Sigma_Eff(2,1))+
     & (-1.0_8*Sigma_Eff(3,1)))
      dSigma_Eff_VM_dSigma_Eff(1,3) = (1.0_8/(2.0_8*Sigma_Eff_VM))*(
     & (-1.0_8*Sigma_Eff(1,1))+
     & (-1.0_8*Sigma_Eff(2,1))+
     & (+2.0_8*Sigma_Eff(3,1)))
      dSigma_Eff_VM_dSigma_Eff(1,4) = (3.0_8*Sigma_Eff(4,1))/Sigma_Eff_VM
      dSigma_Eff_VM_dSigma_Eff(1,5) = (3.0_8*Sigma_Eff(5,1))/Sigma_Eff_VM
      dSigma_Eff_VM_dSigma_Eff(1,6) = (3.0_8*Sigma_Eff(6,1))/Sigma_Eff_VM
      
C: Calculate the Partial Derivative of the Effective Deviatoric Stress with respect to the Effective Stress [6 x 6]
      dSigma_Eff_Dev_dSigma_Eff(1,1) = +(2.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(1,2) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(1,3) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(1,4) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(1,5) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(1,6) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(2,1) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(2,2) = +(2.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(2,3) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(2,4) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(2,5) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(2,6) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(3,1) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(3,2) = -(1.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(3,3) = +(2.0_8/3.0_8)
      dSigma_Eff_Dev_dSigma_Eff(3,4) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(3,5) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(3,6) = 0.0_8
      dSigma_Eff_Dev_dSigma_Eff(4,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8 /)
      dSigma_Eff_Dev_dSigma_Eff(5,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8 /)
      dSigma_Eff_Dev_dSigma_Eff(6,:) = 
     & (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8 /)
      
C: Calculate the Partial Derivative of the Effective Stress with respect to the Stress [6 x 6]
      dSigma_Eff_dSigma = MATMUL(dSigma_Eff_Dev_dSigma_Eff,I_6x6)
      
C: Calculate the Partial Derivative of the Effective Deviatoric Stress with respect to the Stress [6 x 6]
      dSigma_Eff_Dev_dSigma =
     & MATMUL(dSigma_Eff_Dev_dSigma_Eff,dSigma_Eff_dSigma)
      
C: Calculate the Partial Derivative of the Effective Von Mises Stress with respect to the Stress [1 x 6]
      dSigma_Eff_VM_dSigma = 
     & MATMUL(dSigma_Eff_VM_dSigma_Eff,dSigma_Eff_dSigma)
      
C: Calculate the Partial Derivative of the Effective Deviatoric Stress divided by the Effective Von Mises Stress with respect to the Stress [6 x 6]
      dSigma_Eff_Dev_over_Sigma_VM_dSigma = 
     & ((dSigma_Eff_Dev_dSigma*Sigma_Eff_VM)-
     & MATMUL(Sigma_Eff_Dev,dSigma_Eff_VM_dSigma))/
     & (Sigma_Eff_VM**2.0_8)
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Sigma Effective values"
      !    DO i=1,NTENS
      !        WRITE(*,*)Sigma_Eff(i,1)
      !    ENDDO
      !    WRITE(*,*)
      !    WRITE(*,*)"Sigma Effective Mean values"
      !    DO i=1,NTENS
      !        WRITE(*,*)Sigma_Eff_Mean(i,1)
      !    ENDDO
      !    WRITE(*,*)
      !    WRITE(*,*)"Sigma Effective Von Mises value"
      !    WRITE(*,*)Sigma_Eff_VM
      !    WRITE(*,*)
      !    WRITE(*,*)"SHAPE of Sigma Effective Deviatoric"
      !    WRITE(*,*)SHAPE(Sigma_Eff_Dev)
      !    WRITE(*,*)"Sigma Effective Deviatoric values"
      !    DO i=1,NTENS
      !        WRITE(*,11)Sigma_Eff_Dev(i,:)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE H_current_Calculation(STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & H_current,dH_current_dSigma)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 H_sat, k_t, Sigma_Eff_VM,
     & H_current, dH_current_dSigma_Eff_VM
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31), Sigma_Eff(NTENS,1),
     & Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_Dev_dSigma(NTENS,NTENS), dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS),
     & dH_current_dSigma(1,NTENS)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      H_sat = PROPS(17) ! Maximum Saturated Transformation Strain
      k_t = PROPS(18) ! Exponent of Maximum Transformation Strain Function
      Debug = INT(PROPS(38))
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Stress and its Derivatives
      CALL Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
      
C------------------------------------------------------|
C Calculation of Maximum Current Transformation Strain |
C------------------------------------------------------|
      
      H_current = H_sat*(1.0_8-EXP(-k_t*Sigma_Eff_VM))
      
C-----------------------------------------------------------------------------------------------------------|
C Calculation of the Derivative of the Maximum Current Transformation Strain with respect to Stress [1 x 6] |
C-----------------------------------------------------------------------------------------------------------|
      
      dH_current_dSigma_Eff_VM = H_sat*k_t*EXP(-k_t*Sigma_Eff_VM)
      dH_current_dSigma = dH_current_dSigma_Eff_VM*dSigma_Eff_VM_dSigma
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"H_current",H_current
      !    WRITE(*,*)"dH_current_dSigma"
      !    DO i=1,NTENS
      !        WRITE(*,*)dH_current_dSigma(1,i)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C:   |------------------------------------------------------------------|
      
C: Source: Phenomenological Modeling of Induced Transformation Anisotropy in Shape Memory Alloy Actuators
C: (Hartl,Solomou,Lagoudas,Saravanos)
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 Sigma_Eff_VM, H_current
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31),
     & dH_current_dSigma(1,NTENS),
     & Sigma_Eff(NTENS,1), Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_Dev_dSigma(NTENS,NTENS),
     & dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS),
C: Output
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1)
      
C: A is used in Format in Order to Add Text String before the Variable Values
10    FORMAT(A,6(E16.8))
11    FORMAT(6(E24.8E3))
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Debug = INT(PROPS(38))
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Stress and its Derivatives
      CALL Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
C: Call the Subroutine for the Calculation of H_current and its Derivatives 
      CALL H_current_Calculation(STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & H_current,dH_current_dSigma)
      
C--------------------------------------------------------------------------------|
C ---------------- F O R W A R D    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      
C: Calculate the values of the Transformation Tensor for the Forward Transformation [6 x 1]
      Lambda_FWD = (3.0_8/2.0_8)*(H_current*Sigma_Eff_Dev/Sigma_Eff_VM)
C: Transform the Tensorial values to Vector
      Lambda_FWD(1,1) = Lambda_FWD(1,1)
      Lambda_FWD(2,1) = Lambda_FWD(2,1)
      Lambda_FWD(3,1) = Lambda_FWD(3,1)
      Lambda_FWD(4,1) = Lambda_FWD(4,1)*2.0_8
      Lambda_FWD(5,1) = Lambda_FWD(5,1)*2.0_8
      Lambda_FWD(6,1) = Lambda_FWD(6,1)*2.0_8
      
C--------------------------------------------------------------------------------|
C ---------------- R E V E R S E    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      
C: Calculate the values of the Transformation Tensor for the Reverse Transformation
      DO i=1,NTENS
          Lambda_REV(i,1) = STATEV(22+i)
      ENDDO
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Lambda_FWD"
      !    WRITE(*,*)"SHAPE",SHAPE(Lambda_FWD)
      !    DO i=1,NTENS
      !        WRITE(*,FMT=11)Lambda_FWD(i,:)
      !    ENDDO
      !    WRITE(*,*)"Lambda_REV"
      !    WRITE(*,*)"SHAPE",SHAPE(Lambda_REV)
      !    DO i=1,NTENS
      !        WRITE(*,FMT=11)Lambda_REV(i,:)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Hardening_Function(NTENS,PROPS,NPROPS,Model_Parameters,
     & SMA_State,ft_FWD,ft_REV,dft_dMVF_FWD,dft_dMVF_REV)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NPROPS
      
C: Double Precision Arrays
      REAL*8 PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 ML_Activation, Debug, i
      
C: Double Precision Variables
      REAL*8 n1, n2, n3, n4, a1, a2, a3,
     & MVF_current, MVF_reversal_FWD_REV, MVF_reversal_REV_FWD,
     & MVF_fwd, MVF_rev,
     & ft_FWD, ft_REV, dft_dMVF_FWD, dft_dMVF_REV
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), SMA_State(31)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      n1 = PROPS(25)
      n2 = PROPS(26)
      n3 = PROPS(27)
      n4 = PROPS(28)
C: Read the Index for the consideration of Partial Transformation Loops in Transformation Hardening Function
C: If Minor_Loop_Consideration = 0 Partial Transformation Loops are Considered
C: If Minor_Loop_Consideration = 1 Partial Transformation Loops are not Considered
      ML_Activation = INT(PROPS(37))
      Debug = INT(PROPS(38))
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Model Parameters |
C---------------------------------------------------------------------------|
      
      a1 = Model_Parameters(3)
      a2 = Model_Parameters(4)
      a3 = Model_Parameters(5)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      MVF_reversal_FWD_REV = SMA_State(21)
      MVF_reversal_REV_FWD = SMA_State(22)
      
C-------------------------------------|
C Calculate the Scaled Values for MVF |
C-------------------------------------|
      
      MVF_fwd = (MVF_current-MVF_reversal_REV_FWD)/
     & (1.0_8-MVF_reversal_REV_FWD)
      
      IF (MVF_fwd <= 0.0_8) THEN
          MVF_fwd = 0.0001_8
      ELSEIF (MVF_fwd >= 1.0_8) THEN
          MVF_fwd = 0.9999_8
      ENDIF
      
      MVF_rev = (1.0_8/MVF_reversal_FWD_REV)*MVF_current
      
      IF (MVF_rev >= 1.0_8) THEN
          MVF_rev = 0.9999_8
      ELSEIF (MVF_rev <= 0.0_8) THEN
          MVF_rev = 0.0001_8
      ENDIF
      
C------------------------------------------------------------------------------------------|
C Calculate the values of the Hardening Function Expressions and their Partial Derivatives |
C------------------------------------------------------------------------------------------|
      
      IF (ML_Activation == 0) THEN
C: Calculate the Transformation Hardening Function for the Forward Transformation
          ft_FWD = ((1.0_8/2.0_8)*a1*
     & (1.0_8+(MVF_current**n1)-((1.0_8-MVF_current)**n2)))+a3
C: Calculate the Transformation Hardening Function for the Reverse Transformation
          ft_REV = ((1.0_8/2.0_8)*a2*
     & (1.0_8+(MVF_current**n3)-((1.0_8-MVF_current)**n4)))-a3
C: Calculate the Derivative of the Hardening Fuction with respect to the MVF for the Forward Transformation
          dft_dMVF_FWD = (1.0_8/2.0_8)*a1*
     & ((n1*(MVF_current**(n1-1.0_8)))+(n2*((1.0_8-MVF_current)**(n2-1.0_8))))
C: Calculate the Derivative of the Hardening Fuction with respect to the MVF for the Reverse Transformation
          dft_dMVF_REV = (1.0_8/2.0_8)*a2*
     & ((n3*(MVF_current**(n3-1.0_8)))+(n4*((1.0_8-MVF_current)**(n4-1.0_8))))
      ELSE
C: Calculate the Transformation Hardening Function for the Forward Transformation
          ft_FWD = ((1.0_8/2.0_8)*a1*
     & (1.0_8+(MVF_fwd**n1)-((1.0_8-MVF_fwd)**n2)))+a3
C: Calculate the Transformation Hardening Function for the Reverse Transformation
          ft_REV = ((1.0_8/2.0_8)*a2*
     & (1.0_8+(MVF_rev**n3)-((1.0_8-MVF_rev)**n4)))-a3
C: Calculate the Derivative of the Hardening Fuction with respect to the MVF for the Forward Transformation
          dft_dMVF_FWD = (1.0_8/2.0_8)*a1*(1.0_8/(1.0_8-MVF_reversal_REV_FWD))*
     & ((n1*(MVF_fwd**(n1-1.0_8)))+(n2*((1.0_8-MVF_fwd)**(n2-1.0_8))))
C: Calculate the Derivative of the Hardening Fuction with respect to the MVF for the Reverse Transformation
          dft_dMVF_REV = (1.0_8/2.0_8)*a2*(1.0_8/MVF_reversal_FWD_REV)*
     & ((n3*(MVF_rev**(n3-1.0_8)))+(n4*((1.0_8-MVF_rev)**(n4-1.0_8))))
      ENDIF
      
C: Source: Modeling of Partial Transformation Cycles of SMAs with a Modified Hardening Function
C: (Karakalas,Machairas,Solomou,Saravanos)
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"MVF_current",MVF_current
      !    WRITE(*,*)"MVF_reversal_FWD_REV",MVF_reversal_FWD_REV
      !    WRITE(*,*)"MVF_reversal_REV_FWD",MVF_reversal_REV_FWD
      !    WRITE(*,*)"MVF_fwd",MVF_fwd
      !    WRITE(*,*)"MVF_rev",MVF_rev
      !    WRITE(*,*)"ft_FWD",ft_FWD
      !    WRITE(*,*)"ft_REV",ft_REV
      !    WRITE(*,*)"dft_dMVF_FWD",dft_dMVF_FWD
      !    WRITE(*,*)"dft_dMVF_REV",dft_dMVF_REV
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Transformation_Function(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State,
     & Phi_FWD,Phi_REV)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      REAL*8 Density, E_A, E_M, PR, alpha_A, alpha_M, ESH_A, ESH_M,
     & S_A, S_M, Delta_S, Delta_a, Delta_c,
     & rDs0, D, rDu0, Y0_t, T_initial,
     & ft_FWD, ft_REV, dft_dMVF_FWD, dft_dMVF_REV,
     & Y_FWD, Y_REV, Pi_FWD, Pi_FWD_B, Pi_REV, Phi_FWD, Phi_REV
            
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31),
     & BackStress(NTENS,1),
     & Sigma_Eff(NTENS,1), S_Iso(NTENS,NTENS), alpha_Coef(NTENS,1),
     & Stress_current(NTENS,1),
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Density = PROPS(1)
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4)
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
C: Define the Components of the BackStress Matrix [6 x 1]
      DO i=1,NTENS
          IF (PROPS(18+i) == 0.0_8) THEN
              BackStress(i,1) = 1.0e-24
          ELSE
              BackStress(i,1) = PROPS(18+i)
          ENDIF
      ENDDO
      Debug = INT(PROPS(38))
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Model Parameters |
C---------------------------------------------------------------------------|
      
      rDs0 = Model_Parameters(1)
      D = Model_Parameters(2)
      rDu0 = Model_Parameters(6)
      Y0_t = Model_Parameters(7)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State of the Material |
C----------------------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      
C:   |------------------------------------------------------------------|
      
C: Define the Compliance Coefficients of a Generally Isotropic Elastic Material [6 x 6]
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Define the Thermal Expansion Matrix Coefficients
      alpha_Coef(1,1) = 1.0_8
      alpha_Coef(2,1) = 1.0_8
      alpha_Coef(3,1) = 1.0_8
      alpha_Coef(4,1) = 0.0_8
      alpha_Coef(5,1) = 0.0_8
      alpha_Coef(6,1) = 0.0_8
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
C: Call the Subroutine for the Calculation of Transformation Tensor
      CALL Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
C: Call the Subroutine for the Calculation of Hardening Function and its Derivatives
      CALL Hardening_Function(NTENS,PROPS,NPROPS,Model_Parameters,
     & SMA_State,ft_FWD,ft_REV,dft_dMVF_FWD,dft_dMVF_REV)
      
C-----------------------------------------------|
C Calculate the Effective Stress Vector [6 x 1] |
C-----------------------------------------------|
      
      Sigma_Eff = Stress_Current + BackStress
      
C----------------------------------------------------------------------------------------------|
C Calculate the Critical Thermodynamic Driving Force for each Transformation Direction [6 x 1] |
C----------------------------------------------------------------------------------------------|
      
C: Critical Thermodynamic Driving Force for the Forward Transformation
      Y_FWD = Y0_t+(D*MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_FWD)))
C: Critical Thermodynamic Driving Force for the Reverse Transformation
      Y_REV = Y0_t+(D*MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_REV)))
      
C------------------------------------------------------------------------------------|
C Calculate the Total Thermodynamic Force for each Transformation Direction [SCALAR] |
C------------------------------------------------------------------------------------|
      
C: Total Thermodynamic Force of the Forward Transformation
      Pi_FWD = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_FWD))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_FWD
      
C: Total Thermodynamic Force of the Reverse Transformation
      Pi_REV = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_REV))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_REV
      
C------------------------------------------------------------------------------|
C Calculate the Transformation Function value for each Transformation [SCALAR] |
C------------------------------------------------------------------------------|
      
C: Transformation Function value for the Forward Transformation
      Phi_FWD = +Pi_FWD-Y_FWD
C: Transformation Function value for the Reverse Transformation
      Phi_REV = -Pi_REV-Y_REV
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Y_FWD",Y_FWD
      !    WRITE(*,*)"Y_REV",Y_REV
      !    WRITE(*,*)"Pi_FWD",Pi_FWD
      !    WRITE(*,*)"Pi_REV",Pi_REV
      !    WRITE(*,*)"Phi_FWD",Phi_FWD
      !    WRITE(*,*)"Phi_REV",Phi_REV
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Transformation_Direction(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,KINC,
     & Model_Parameters,Initial_State,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS, KINC
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, i, TO_index_pre,
     & Completion_FWD, Completion_REV, Direction_FWD, Direction_REV,
     & Direction_Index, TO_index,
     & Reversal_FWD_REV, Reversal_REV_FWD
      
C: Double Precision Variables
      REAL*8 Tolerance_Phi,
     & Phi_FWD, Phi_REV,
     & Phi_FWD_pre, Phi_REV_pre,
     & Delta_Phi_FWD, Delta_Phi_REV,
     & MVF_reversal_FWD_REV, MVF_reversal_REV_FWD
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Tolerance_Phi = PROPS(29)
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      TO_index_pre = INT(SMA_State(14))
      Completion_FWD = INT(SMA_State(17))
      Completion_REV = INT(SMA_State(18))
      Reversal_FWD_REV = INT(SMA_State(19))
      Reversal_REV_FWD = INT(SMA_State(20))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      Phi_FWD_pre = STATEV(29)
      Phi_REV_pre = STATEV(30)
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Calculate the current Transformation Function values
      CALL Transformation_Function(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State,
     & Phi_FWD,Phi_REV)
      
C:   |------------------------------------------------------------------|
      
C: If the values of Phi for both Transformations are less than or equal the Tolerance value then they are set to 0.0
      IF (ABS(Phi_FWD_pre) <= Tolerance_Phi) THEN
          Phi_FWD_pre = 0.0_8
      ENDIF
      IF (ABS(Phi_REV_pre) <= Tolerance_Phi) THEN
          Phi_REV_pre = 0.0_8
      ENDIF
      IF (ABS(Phi_FWD) <= Tolerance_Phi) THEN
          Phi_FWD = 0.0_8
      ENDIF
      IF (ABS(Phi_REV) <= Tolerance_Phi) THEN
          Phi_REV = 0.0_8
      ENDIF
      
C------------------------------------------------------------------------------------------------|
C Calculate the Difference between the Current and the Previous value of Transformation Function |
C------------------------------------------------------------------------------------------------|
      
      Delta_Phi_FWD = Phi_FWD-Phi_FWD_pre
      Delta_Phi_REV = Phi_REV-Phi_REV_pre
      
C-----------------------------------------------|
C Determine the Direction of the Transformation |
C-----------------------------------------------|
      
      !IF(Delta_Phi_FWD > 0.0_8)THEN
      IF((Delta_Phi_FWD > 0.0_8).AND.(Delta_Phi_REV <= 0.0_8))THEN
          ! Direction towards Forward Transformation
          Direction_FWD = 1
          Direction_REV = 0
          Direction_Index = +1
      !ELSEIF(Delta_Phi_REV > 0.0_8)THEN
      ELSEIF((Delta_Phi_FWD <= 0.0_8).AND.(Delta_Phi_REV > 0.0_8))THEN
          ! Direction towards Reverse Transformation
          Direction_FWD = 0
          Direction_REV = 1
          Direction_Index = -1
      ELSE
          Direction_FWD = INT(STATEV(15))
          Direction_REV = INT(STATEV(16))
          IF((Delta_Phi_FWD==0.0_8).AND.(Delta_Phi_REV==0.0_8))THEN
              Direction_Index = 0
          ELSE
              Direction_Index = 2
          ENDIF
      ENDIF
      
C------------------------------------|
C Determine Transformation Occurence |
C------------------------------------|
      
      IF ((Phi_FWD > 0.0_8) .AND. (Delta_Phi_FWD > 0.0_8) .AND. 
     & (Completion_FWD == 0)) THEN
          ! Forward Transformation
          ! Assign a value to Transformation Occurence Index
          TO_index = +1
      ELSEIF ((Phi_REV > 0.0_8) .AND. (Delta_Phi_REV > 0.0_8) .AND. 
     & (Completion_REV == 0)) THEN
          ! Reverse Transformation
          ! Assign a value to Transformation Direction Index
          TO_index = -1
      ELSE
          ! Elastic Response
          ! Assign a value to Transformation Direction Index
          TO_index = 0
      ENDIF
      
C-----------------------------------------------------------------|
C Update the proper contents in the current State of the Material |
C-----------------------------------------------------------------|
      
      SMA_State(14) = TO_index
      SMA_State(15) = Direction_FWD
      SMA_State(16) = Direction_REV
      SMA_State(29) = Phi_FWD
      SMA_State(30) = Phi_REV
      
C------------------------------------------------------------------|
C Determine if there is a Reversal in the Transformation Direction |
C------------------------------------------------------------------|
      
      CALL Transformation_Reversal_Detection(STATEV,NSTATV,
     & PROPS,NPROPS,SMA_State)
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Completion_FWD",Completion_FWD
      !    WRITE(*,*)"Completion_REV",Completion_REV
      !    WRITE(*,*)"Direction_FWD",Direction_FWD
      !    WRITE(*,*)"Direction_REV",Direction_REV
      !    WRITE(*,*)"Delta_Phi_FWD",Delta_Phi_FWD
      !    WRITE(*,*)"Delta_Phi_REV",Delta_Phi_REV
      !    WRITE(*,*)"TO_index",TO_index
      !    WRITE(*,*)"Phi_FWD",Phi_FWD
      !    WRITE(*,*)"Phi_REV",Phi_REV
      !    WRITE(*,*)"Direction_FWD",Direction_FWD
      !    WRITE(*,*)"Direction_REV",Direction_REV
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Transformation_Reversal_Detection(STATEV,NSTATV,
     & PROPS,NPROPS,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NSTATV, NPROPS
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, Direction_FWD, Direction_REV,
     & Direction_FWD_pre, Direction_REV_pre,
     & Reversal_FWD_REV, Reversal_REV_FWD
      
C: Double Precision Variables
      REAL*8 MVF_current, MVF_reversal_FWD_REV, MVF_reversal_REV_FWD
      
C: Double Precision Arrays
      REAL*8 SMA_State(31)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      Direction_FWD = INT(SMA_State(15))
      Direction_REV = INT(SMA_State(16))
      Reversal_FWD_REV = INT(SMA_State(19))
      Reversal_REV_FWD = INT(SMA_State(20))
      MVF_reversal_FWD_REV = SMA_State(21)
      MVF_reversal_REV_FWD = SMA_State(22)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_current = STATEV(1)
      Direction_FWD_pre = INT(STATEV(15))
      Direction_REV_pre = INT(STATEV(16))
      
C--------------------------------------------------|
C Detect Transformation Reversal Points if present |
C--------------------------------------------------|
      
      IF (Direction_REV_pre == 1 .AND. Direction_FWD == 1) THEN
C: Set the Index for Transformation Reversal from the Forward to Reverse
          Reversal_FWD_REV = 0
C: Set the Index for Transformation Reversal from the Reverse to Forward
          Reversal_REV_FWD = 1
C: Update the MVF at the Reversal Point from the Reverse to Forward Transformation
          MVF_reversal_REV_FWD = MVF_current
      ENDIF
      
      IF (Direction_FWD_pre == 1 .AND. Direction_REV == 1) THEN
C: Set the Index for Transformation Reversal from the Forward to Reverse
          Reversal_FWD_REV = 1
C: Set the Index for Transformation Reversal from the Reverse to Forward
          Reversal_REV_FWD = 0
C: Update the MVF at the Reversal Point from the Forward to Reverse Transformation
          MVF_reversal_FWD_REV = MVF_current
      ENDIF
      
C-----------------------------------------------------------------|
C Update the proper contents in the current State of the Material |
C-----------------------------------------------------------------|
      
      SMA_State(19) = Reversal_FWD_REV
      SMA_State(20) = Reversal_REV_FWD
      SMA_State(21) = MVF_reversal_FWD_REV
      SMA_State(22) = MVF_reversal_REV_FWD
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Update_Internal_State(STATEV,STRAN,DSTRAN,TEMP,
     & DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER Debug, TO_index
      
C: Double Precision Variables
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
C: Read the Index to Determine Transformation Occurence
      TO_index = INT(SMA_State(14))
      
C----------------------------------------------------------------------------------------|
C If Transformation Occurs then Call the Subroutine to initiate Return Mapping Algorithm |
C----------------------------------------------------------------------------------------|
      
      IF ((TO_index == +1).OR.(TO_index == -1)) THEN
          !IF (Debug == 1) THEN
          !    WRITE(*,*)
          !    WRITE(*,*)"Transformation Occurence"
          !    WRITE(*,*)
          !ENDIF
          CALL Return_Mapping_Algorithm(STATEV,STRAN,DSTRAN,TEMP,
     & DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
      ELSE
          !IF (Debug == 1) THEN
          !    WRITE(*,*)
          !    WRITE(*,*)"No Transformation Occurs - Thermoelastic Response"
          !    WRITE(*,*)
          !ENDIF
      ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Transformation_Strain_Residual(STATEV,STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & R_FWD,R_REV)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER Debug, i
      
C: Double Precision Variables
      REAL*8 MVF_current, MVF_previous
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31), I_6x6(NTENS,NTENS),
     & et_current(NTENS,1), et_previous(NTENS,1),
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1),
     & R_FWD(NTENS,1), R_REV(NTENS,1)
      
11    FORMAT(6(E24.8E3))
      
C----------------------------------------------------------------------|
C Define any Auxiliary Matrices required for the particular Subroutine |
C----------------------------------------------------------------------|
      
C: Define an Identity Matrix - [6 x 6]
      I_6x6(1,:) = (/ 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(2,:) = (/ 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(3,:) = (/ 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8 /)
      I_6x6(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8 /)
      I_6x6(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8 /)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_previous = STATEV(1)
      DO i=1,NTENS
          et_previous(i,1) = STATEV(1+i)
      ENDDO
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Transformation Tensor
      CALL Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
      
C--------------------------------------------------------------------|
C Calculate the current values of the Transformation Strain Residual |
C--------------------------------------------------------------------|
      
C: Calculate the current values of the Transformation Strain Residual for the Forward Transformation [6 x 1]
      R_FWD = (-et_current+et_previous)+
     & (Lambda_FWD*(MVF_current-MVF_previous))
      
C: Calculate the current values of the Transformation Strain Residual for the Reverse Transformation [6 x 1]
      R_REV = (-et_current+et_previous)+
     & (Lambda_REV*(MVF_current-MVF_previous))
      
C:--------------------------------------------------------------------------------|
!     P R I N T   C O M M A N D S   for   D E B U G G I N G   M E S S A G E S     !
C:--------------------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"R_FWD"
      !    WRITE(*,*)"SHAPE",SHAPE(R_FWD)
      !    DO i=1,NTENS
      !        WRITE(*,FMT=11)R_FWD(i,:)
      !    ENDDO
      !    WRITE(*,*)"R_REV"
      !    WRITE(*,*)"SHAPE",SHAPE(R_REV)
      !    DO i=1,NTENS
      !        WRITE(*,FMT=11)R_REV(i,:)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Return_Mapping_Algorithm(STATEV,STRAN,DSTRAN,TEMP,
     & DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, TO_index, i,
     & RMA_Increment_Index, RMA_Exit, Convergence, MVF_correction_index,
     & Completion_FWD, Completion_REV
      
C: Double Precision Variables
      REAL*8 Tolerance_Phi, Tolerance_R, MVF_current,
     & Phi_FWD, Phi_REV, dPhi_dMVF, Delta_MVF, MVF_pre_correction
            
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31),
     & et_current(NTENS,1),
     & R_FWD(NTENS,1), R_REV(NTENS,1),
     & dR_det(NTENS,NTENS), dR_dMVF(NTENS,1), dPhi_det(1,NTENS),
     & Newton_Raphson_Coefficients(NTENS+1,NTENS+1),
     & Newton_Raphson_Equation_Values(NTENS+1,1),
     & Inverse(NTENS+1,NTENS+1), Newton_Raphson_Deltas(NTENS+1,1),
     & Delta_et(NTENS,1), et_pre_correction(NTENS,1), RMA_Residual(2,1)
      
11    FORMAT(7(E24.8E3))
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Tolerance_Phi = PROPS(29)
      Tolerance_R = PROPS(29)
C: Read the Index for Enabling Debbugin Messages
      Debug = INT(PROPS(38))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_current = STATEV(1)
      DO i=1,NTENS
          et_current(i,1) = STATEV(1+i)
      ENDDO
      Completion_FWD = INT(STATEV(17))
      Completion_REV = INT(STATEV(18))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      TO_index = INT(SMA_State(14))
      
C:-----------------------------------------------------------------------!
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Return Mapping Subroutine Input"
      !    WRITE(*,*)"MVF_current",MVF_current
      !    WRITE(*,*)"TO_index",TO_index
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine to Calculate the Updated values of the Transformation Strain Residual
      CALL Transformation_Strain_Residual(STATEV,STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & R_FWD,R_REV)
C: Calculate the current Transformation Function values
      CALL Transformation_Function(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State,
     & Phi_FWD,Phi_REV)
      
C:   |------------------------------------------------------------------|
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Values of R and Phi prior to RMA Initialization"
      !    IF (TO_index == +1) THEN
      !        DO i=1,NTENS
      !            WRITE(*,*)"R_FWD",R_FWD(i,1)
      !        ENDDO
      !        WRITE(*,*)"Phi_FWD",Phi_FWD
      !    ENDIF
      !    IF (TO_index == -1) THEN
      !        DO i=1,NTENS
      !            WRITE(*,*)"R_REV",R_REV(i,1)
      !        ENDDO
      !        WRITE(*,*)"Phi_REV",Phi_REV
      !    ENDIF
      !    WRITE(*,*)
      !ENDIF
      
C:   |------------------------------------------------------------------|
      
C: Set Return Mapping Algorithm Index equal to Zero
      RMA_Increment_Index = 0
C: Set the Index for stopping RMA iterations equal to Zero
      RMA_Exit = 0
C: Set the Convergence Index equal to Zero
      Convergence = 0
      
C--------------------------------------------------------------------------------|
C ---------------- F O R W A R D    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      
      IF (TO_index == +1) THEN
          DO WHILE ((Convergence == 0).AND.(RMA_Exit == 0))
C: Print the value of the current RMA Increment
              IF (Debug == 1) THEN
                  WRITE(*,*)
                  WRITE(*,*)"RMA_Increment_Index",RMA_Increment_Index
                  WRITE(*,*)
              ENDIF
C: Call the Subroutine for the Calculation of the Partial Derivatives required by RMA
              CALL RMA_Derivatives(STATEV,STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & NSTATV,PROPS,NPROPS,Model_Parameters,Initial_State,SMA_State,
     & dR_det,dR_dMVF,dPhi_det,dPhi_dMVF)
C: Initialize the Matrix that contains the Coefficients of the required Derivatives for the Newton-Raphson Method
              Newton_Raphson_Coefficients = 0.0_8
C: Assemble the Coefficients Matrix for Newton-Raphson [7 x 7]
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(i,1:NTENS) = 
     & dR_det(i,1:NTENS)
              ENDDO
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(NTENS+1,i) = dPhi_det(1,i)
              ENDDO
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(i,NTENS+1) = dR_dMVF(i,1)
              ENDDO
              Newton_Raphson_Coefficients(NTENS+1,NTENS+1) = dPhi_dMVF
C: Assemble the matrix containing the Calculated Current values of the equations used in Newton-Raphson
              DO i=1,NTENS
                  Newton_Raphson_Equation_Values(i,1) = R_FWD(i,1)
              ENDDO
              Newton_Raphson_Equation_Values(NTENS+1,1) = Phi_FWD
C: Initialize the Inverse Matrix
              Inverse = 0.0_8
C: Perform Matrix Inversion
              CALL Matrix_Inversion
     & (Newton_Raphson_Coefficients,Inverse,NTENS+1)
C: Calculate the Updated Correnction
              Newton_Raphson_Deltas = 
     & -MATMUL(Inverse,Newton_Raphson_Equation_Values)
C: Separate the Correction of Transformation Strain values
              Delta_et(:,1) = Newton_Raphson_Deltas(1:NTENS,1)
C: Separate the Correction of MVF
              Delta_MVF = Newton_Raphson_Deltas(NTENS+1,1)
C: Calculate the current-updated values of the Transfromation Strain
              et_current = et_current+Delta_et
C: Calculate the current-updated value of the Martensitic Volume Fraction
              MVF_current = MVF_current+Delta_MVF
C: Store the calculated values of Transformation Strain Vector prior to any Correnctions
              et_pre_correction = et_current
C: Store the calculated value of MVF prior to any Correnctions
              MVF_pre_correction = MVF_current
C: Check for Violation of Martensitic Volume Fraction Boundary values
              IF ((MVF_current < 0.0001_8) .OR. 
     & (MVF_current > 0.9999_8)) THEN
C: Update Martensitic Volume Fraction correction Index
                  MVF_correction_index = 1
                  IF (MVF_current < 0.0001_8) THEN
                      !IF (Debug == 1) THEN
                      !    WRITE(*,*)
                      !    WRITE(*,*)"MVF_current < 0"
                      !    WRITE(*,*)"Reverse Transformation Completed"
                      !    WRITE(*,*)
                      !ENDIF
C: Set the MVF to the Proper Boundary Value
                      MVF_current = 0.0001_8
C: Update the Indeces regarding Transformation Completion
                      !Completion_FWD = 0
                      !Completion_REV = 1
C: Update the values of the Transformation Strain Vector
                      ! TO BE DEFINED !
                  ELSEIF (MVF_current > 0.9999_8) THEN
                      !IF (Debug == 1) THEN
                      !    WRITE(*,*)
                      !    WRITE(*,*)"MVF_current > 1"
                      !    WRITE(*,*)"Forward Transformation Completed"
                      !    WRITE(*,*)
                      !ENDIF
C: Set the MVF to the Proper Boundary Value
                      MVF_current = 0.9999_8
C: Update the Indeces regarding Transformation Completion
                      !Completion_FWD = 1
                      !Completion_REV = 0
C: Update the values of the Transformation Strain Vector
                      ! TO BE DEFINED !
                  ENDIF
              ENDIF
C: Update the State of the Material
              SMA_State(1) = MVF_current
              DO i=1,NTENS
                  SMA_State(1+i) = et_current(i,1)
              ENDDO
              !SMA_State(17) = Completion_FWD
              !SMA_State(18) = Completion_REV
C: Call the Subroutine to Calculate the Updated values of the Transformation Strain Residual
              CALL Transformation_Strain_Residual(STATEV,STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & R_FWD,R_REV)
C: Call the Subroutine to Calculate the Updated values of the Transformation Function
              CALL Transformation_Function(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State,
     & Phi_FWD,Phi_REV)
C: Determine the Maximum Transformation Strain Residual
              RMA_Residual(1,1) = MAXVAL(ABS(R_FWD(1:NTENS,1)))
C: Determine the Residual of the Transformation Function
              RMA_Residual(2,1) = ABS(Phi_FWD)
C: Check if Convergence has been reached
              IF ((RMA_Residual(1,1) > Tolerance_R) .OR. 
     & (RMA_Residual(2,1) > Tolerance_Phi)) THEN
                  Convergence = 0
                  IF (Debug == 1) THEN
                      WRITE(*,*)
                      WRITE(*,*)"Convergence has not been Reached"
                      WRITE(*,*)
                  ENDIF
              ELSE
                  Convergence = 1
                  IF (Debug == 1) THEN
                      WRITE(*,*)
                      WRITE(*,*)"Convergence has been Reached"
                      WRITE(*,*)
                  ENDIF
              ENDIF
C: Update RMA Increment Index
              RMA_Increment_Index = RMA_Increment_Index+1        
              IF (RMA_Increment_Index >= 10) THEN
                  RMA_Exit = 1
                  !IF (Debug == 1) THEN
                  !    WRITE(*,*)
                  !    WRITE(*,*)"RMA_Increment_Index > 10"
                  !    WRITE(*,*)"Delta_MVF",Delta_MVF
                  !    WRITE(*,*)"RMA_Residual(1,1)",RMA_Residual(1,1)
                  !    WRITE(*,*)"RMA_Residual(2,1)",RMA_Residual(2,1)
                  !    WRITE(*,*)
                  !ENDIF
              ENDIF
          ENDDO
C: Update the State of the Material
          SMA_State(1) = MVF_current
          DO i=1,NTENS
              SMA_State(1+i) = et_current(i,1)
          ENDDO
          !SMA_State(17) = Completion_FWD
          !SMA_State(18) = Completion_REV
          !SMA_State(21) = MVF_current
          DO i=1,NTENS
              SMA_State(22+i) = et_current(i,1)/MVF_current
          ENDDO
          SMA_State(29) = Phi_FWD
          SMA_State(30) = Phi_REV
          SMA_State(31) = MVF_correction_index
          !IF (Debug == 1) THEN
          !    WRITE(*,*)"RMA_Residual(1,1)",RMA_Residual(1,1)
          !    WRITE(*,*)"RMA_Residual(2,1)",RMA_Residual(2,1)
          !    WRITE(*,*)"MVF_current",MVF_current
          !    WRITE(*,*)"et_current"
          !    DO i=1,NTENS
          !        WRITE(*,*)et_current(i,1)
          !    ENDDO
          !ENDIF
C--------------------------------------------------------------------------------|
C ---------------- R E V E R S E    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
          
      ELSEIF (TO_index == -1) THEN
          DO WHILE ((Convergence == 0).AND.(RMA_Exit == 0))
C: Print the value of the current RMA Increment
              IF (Debug == 1) THEN
                  WRITE(*,*)
                  WRITE(*,*)"RMA_Increment_Index",RMA_Increment_Index
                  WRITE(*,*)
              ENDIF
C: Call the Subroutine for the Calculation of the Partial Derivatives required by RMA
              CALL RMA_Derivatives(STATEV,STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & NSTATV,PROPS,NPROPS,Model_Parameters,Initial_State,SMA_State,
     & dR_det,dR_dMVF,dPhi_det,dPhi_dMVF)
C: Initialize the Matrix that contains the Coefficients of the required Derivatives for the Newton-Raphson Method
              Newton_Raphson_Coefficients = 0.0_8
C: Assemble the Coefficients Matrix for Newton-Raphson [7 x 7]
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(i,1:NTENS) = 
     & dR_det(i,1:NTENS)
              ENDDO
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(NTENS+1,i) = dPhi_det(1,i)
              ENDDO
              DO i=1,NTENS
                  Newton_Raphson_Coefficients(i,NTENS+1) = dR_dMVF(i,1)
              ENDDO
              Newton_Raphson_Coefficients(NTENS+1,NTENS+1) = dPhi_dMVF
C: Assemble the matrix containing the Calculated Current values of the equations used in Newton-Raphson
              DO i=1,NTENS
                  Newton_Raphson_Equation_Values(i,1) = R_REV(i,1)
              ENDDO
              Newton_Raphson_Equation_Values(NTENS+1,1) = Phi_REV
C: Initialize the Inverse Matrix
              Inverse = 0.0_8
C: Perform Matrix Inversion
              CALL Matrix_Inversion
     & (Newton_Raphson_Coefficients,Inverse,NTENS+1)
C: Calculate the Updated Correnction
              Newton_Raphson_Deltas = 
     & -MATMUL(Inverse,Newton_Raphson_Equation_Values)
C: Separate the Correction of Transformation Strain values
              Delta_et(:,1) = Newton_Raphson_Deltas(1:NTENS,1)
C: Separate the Correction of MVF
              Delta_MVF = Newton_Raphson_Deltas(NTENS+1,1)
C: Calculate the current-updated values of the Transfromation Strain
              et_current = et_current+Delta_et
C: Calculate the current-updated value of the Martensitic Volume Fraction
              MVF_current = MVF_current+Delta_MVF
C: Store the calculated values of Transformation Strain Vector prior to any Correnctions
              et_pre_correction = et_current
C: Store the calculated value of MVF prior to any Correnctions
              MVF_pre_correction = MVF_current
C: Check for Violation of Martensitic Volume Fraction Boundary values
              IF ((MVF_current < 0.0001_8) .OR. 
     & (MVF_current > 0.9999_8)) THEN
C: Update Martensitic Volume Fraction correction Index
                  MVF_correction_index = 1
                  IF (MVF_current < 0.0001_8) THEN
                      !IF (Debug == 1) THEN
                      !    WRITE(*,*)
                      !    WRITE(*,*)"MVF_current < 0"
                      !    WRITE(*,*)"Reverse Transformation Completed"
                      !    WRITE(*,*)
                      !ENDIF
C: Set the MVF to the Proper Boundary Value
                      MVF_current = 0.0001_8
C: Update the Indeces regarding Transformation Completion
                      !Completion_FWD = 0
                      !Completion_REV = 1
C: Update the values of the Transformation Strain Vector
                      ! TO BE DEFINED !
                  ELSEIF (MVF_current > 0.9999_8) THEN
                      !IF (Debug == 1) THEN
                      !    WRITE(*,*)
                      !    WRITE(*,*)"MVF_current > 1"
                      !    WRITE(*,*)"Forward Transformation Completed"
                      !    WRITE(*,*)
                      !ENDIF
C: Set the MVF to the Proper Boundary Value
                      MVF_current = 0.9999_8
C: Update the Indeces regarding Transformation Completion
                      !Completion_FWD = 1
                      !Completion_REV = 0
C: Update the values of the Transformation Strain Vector
                      ! TO BE DEFINED !
                  ENDIF
              ENDIF
C: Update the State of the Material
              SMA_State(1) = MVF_current
              DO i=1,NTENS
                  SMA_State(1+i) = et_current(i,1)
              ENDDO
              !SMA_State(17) = Completion_FWD
              !SMA_State(18) = Completion_REV
C: Call the Subroutine to Calculate the Updated values of the Transformation Function
              CALL Transformation_Function(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State,
     & Phi_FWD,Phi_REV)
C: Call the Subroutine to Calculate the Updated values of the Transformation Strain Residual
              CALL Transformation_Strain_Residual(STATEV,STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & R_FWD,R_REV)
C: Determine the Maximum Transformation Strain Residual
              RMA_Residual(1,1) = MAXVAL(ABS(R_REV(:,1)))
C: Determine the Residual of the Transformation Function
              RMA_Residual(2,1) = ABS(Phi_REV)
C: Check if Convergence has been reached
              IF ((RMA_Residual(1,1) > Tolerance_R) .OR. 
     & (RMA_Residual(2,1) > Tolerance_Phi)) THEN
                  Convergence = 0
                  IF (Debug == 1) THEN
                      WRITE(*,*)
                      WRITE(*,*)"Convergence has not been Reached"
                      WRITE(*,*)
                  ENDIF
              ELSE
                  Convergence = 1
                  IF (Debug == 1) THEN
                      WRITE(*,*)
                      WRITE(*,*)"Convergence has been Reached"
                      WRITE(*,*)
                  ENDIF
              ENDIF
C: Update RMA Increment Index
              RMA_Increment_Index = RMA_Increment_Index+1        
              IF (RMA_Increment_Index >= 10) THEN
                  RMA_Exit = 1
                  !IF (Debug == 1) THEN
                  !    WRITE(*,*)
                  !    WRITE(*,*)"RMA_Increment_Index > 10"
                  !    WRITE(*,*)"Delta_MVF",Delta_MVF
                  !    WRITE(*,*)"RMA_Residual(1,1)",RMA_Residual(1,1)
                  !    WRITE(*,*)"RMA_Residual(2,1)",RMA_Residual(2,1)
                  !    WRITE(*,*)
                  !ENDIF
              ENDIF
          ENDDO
C: Update the State of the Material
          SMA_State(1) = MVF_current
          DO i=1,NTENS
              SMA_State(1+i) = et_current(i,1)
          ENDDO
          !SMA_State(17) = Completion_FWD
          !SMA_State(18) = Completion_REV
          SMA_State(29) = Phi_FWD
          SMA_State(30) = Phi_REV
          SMA_State(31) = MVF_correction_index
          !IF (Debug == 1) THEN
          !    WRITE(*,*)"RMA_Residual(1,1)",RMA_Residual(1,1)
          !    WRITE(*,*)"RMA_Residual(2,1)",RMA_Residual(2,1)
          !    WRITE(*,*)"MVF_current",MVF_current
          !    WRITE(*,*)"et_current"
          !    DO i=1,NTENS
          !        WRITE(*,*)et_current(i,1)
          !    ENDDO
          !ENDIF
      ELSE
          SMA_State(1) = MVF_current
          DO i=1,NTENS
              SMA_State(1+i) = et_current(i,1)
          ENDDO
          SMA_State(31) = MVF_correction_index
      ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE RMA_Derivatives(STATEV,STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & NSTATV,PROPS,NPROPS,Model_Parameters,Initial_State,SMA_State,
     & dR_det,dR_dMVF,dPhi_det,dPhi_dMVF)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER Debug, TO_index, i
      
C: Double Precision Variables
      REAL*8 E_A, E_M, PR, alpha_A, alpha_M, ESH_A, ESH_M, S_A, S_M,
     & Delta_S, Delta_a, Delta_c, D, T_initial,
     & c_current, H_sat, k_t, MVF_current, MVF_previous,
     & Sigma_Eff_VM, H_current, dH_current_dSigma_Eff_VM,
     & dY_dMVF_FWD, dY_dMVF_REV,
     & ft_FWD, ft_REV, dft_dMVF_FWD, dft_dMVF_REV,
     & dPhi_dMVF_FWD, dPhi_dMVF_FWD_B, dPhi_dMVF_REV,
     & dPhi_dMVF
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31),
     & I_6x6(NTENS,NTENS), Backstress(NTENS,1),
     & et_initial(NTENS,1), et_current(NTENS,1), Strain(NTENS,1),
     & S_Iso(NTENS,NTENS), C_Iso(NTENS,NTENS), alpha_Coef(NTENS,1),
     & S_current(NTENS,NTENS), a_current(NTENS,1), S_Inverted(NTENS,NTENS),
     & Stress_current(NTENS,1),
     & Sigma_Eff(NTENS,1), Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS),
     & dSigma_det(NTENS,NTENS), dS_dMVF(NTENS,NTENS), dC_dMVF(NTENS,NTENS),
     & dSigma_dMVF(NTENS,1), dH_current_dSigma(1,NTENS),
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1),
     & dLambda_dSigma_FWD(NTENS,NTENS),
     & dLambda_det_FWD(NTENS,NTENS), dLambda_dMVF_FWD(NTENS,1),
     & dLambda_dSigma_REV(NTENS,NTENS),
     & dLambda_det_REV(NTENS,NTENS), dLambda_dMVF_REV(NTENS,1),
     & dY_dSigma_FWD(1,NTENS), dY_det_FWD(1,NTENS),
     & dY_dSigma_REV(1,NTENS), dY_det_REV(1,NTENS),
     & dR_det_FWD(NTENS,NTENS), dR_dMVF_FWD(NTENS,1),
     & dR_det_REV(NTENS,NTENS), dR_dMVF_REV(NTENS,1),
     & dPhi_det_FWD(1,NTENS), dPhi_det_REV(1,NTENS),
     & dR_det(NTENS,NTENS), dR_dMVF(NTENS,1), dPhi_det(1,NTENS)
            
C: A is used in Format in Order to Add Text String before the Variable Values
10    FORMAT(A,6(E16.8))
11    FORMAT(6(E24.8E3))
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4) ! Poisson's Ratio
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
      H_sat = PROPS(17) ! Maximum Saturated Transformation Strain
      k_t = PROPS(18) ! Exponent of Maximum Transformation Strain Function
C: Define the Components of the BackStress Matrix [6 x 1]
      DO i=1,NTENS
          IF (PROPS(18+i) == 0.0_8) THEN
              BackStress(i,1) = 1.0e-24
          ELSE
              BackStress(i,1) = PROPS(18+i)
          ENDIF
      ENDDO
      Debug = INT(PROPS(38))
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Model Parameters |
C---------------------------------------------------------------------------|
      
      D = Model_Parameters(2)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State of the Material |
C----------------------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      DO i=1,NTENS
          et_initial(i,1) = Initial_State(2+i)
      ENDDO
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      TO_index = INT(SMA_State(14))
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"RMA Derivatives Calculation"
      !    WRITE(*,*)"MVF_current",MVF_current
      !    DO i=1,NTENS
      !        WRITE(*,*)"et_current",et_current(i,1)
      !    ENDDO
      !    WRITE(*,*)"TO_index",TO_index
      !    WRITE(*,*)
      !ENDIF
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_previous = STATEV(1)
      
C----------------------------------------------------------------------|
C Define any Auxiliary Matrices required for the particular Subroutine |
C----------------------------------------------------------------------|
      
C: Define an Identity Matrix - [6 x 6]
      I_6x6(1,:) = (/ 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(2,:) = (/ 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(3,:) = (/ 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      I_6x6(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8, 0.0_8 /)
      I_6x6(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8, 0.0_8 /)
      I_6x6(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 1.0_8 /)
      
C:   |------------------------------------------------------------------|
      
C: Define the Compliance Matrix Coefficients
      
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Perform Matrix Inversion
      CALL Matrix_Inversion(S_Iso,C_Iso,NTENS)
      
C: Redefine the Original Matrix prior to the Inversion
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Define the Thermal Expansion Matrix Coefficients
      alpha_Coef(1,1) = 1.0_8
      alpha_Coef(2,1) = 1.0_8
      alpha_Coef(3,1) = 1.0_8
      alpha_Coef(4,1) = 0.0_8
      alpha_Coef(5,1) = 0.0_8
      alpha_Coef(6,1) = 0.0_8
      
C-------------------------------------------------------------|
C Calculate the Total Strain provided by Abaqus Global Solver |
C-------------------------------------------------------------|
      
      Do i=1,NTENS
          Strain(i,1) = STRAN(i)+DSTRAN(i)
      ENDDO
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Material Properties
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
C: Invert the Compliance to Calculate the Stiffness Matrix
      CALL Matrix_Inversion(S_current,S_Inverted,NTENS)
C: Call the Subroutine for the Calculation of Effective Material Properties to redefine S_current
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
      
C----------------------------------------------------------------------------------|
C Calculate the Partial Derivatives of Stress (Stress is a function of et and MVF) |
C----------------------------------------------------------------------------------|
      
C: Calculate the Partial Derivative of Stress with respect to the Transformation Strain [6 x 6]
      dSigma_det = -S_Inverted
      
C: Calculate the Partial Derivative of the Compliance Matrix with respect to the MVF [6 x 6]
      dS_dMVF = ((E_A-E_M)/(E_A*E_M))*S_Iso
      
C: Calculate the Partial Derivative of the Stiffness Matrix with respect to the MVF [6 x 6]
      dC_dMVF = MATMUL(MATMUL(-S_Inverted,dS_dMVF),S_Inverted)
      
C: Calculate the Partial Derivative of Stress with respect to the Martensitic Volume Fraction [6 x 1]
      dSigma_dMVF = MATMUL(dC_dMVF,(Strain-
     & (a_current*(TEMP+DTEMP-T_initial))-(et_current-et_initial)))+
     & (MATMUL(S_Inverted,Delta_a*alpha_coef)*(TEMP+DTEMP-T_initial))
      
C:   |------------------------------------------------------------------|
      
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
C: Call the Subroutine for the Calculation of Effective Stress and its Derivatives
      CALL Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
C: Call the Subroutine for the Calculation of H_current and its Derivatives
      CALL H_current_Calculation(STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & H_current,dH_current_dSigma)
C: Call the Subroutine for the Calculation of Transformation Tensor
      CALL Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
C: Call the Subroutine for the Calculation of Hardening Function and its Derivatives
      CALL Hardening_Function(NTENS,PROPS,NPROPS,Model_Parameters,
     & SMA_State,ft_FWD,ft_REV,dft_dMVF_FWD,dft_dMVF_REV)
      
C:   |------------------------------------------------------------------|
      
C: Calculate the RMA Derivatives based on the Current Direction of the Transformation
      
C--------------------------------------------------------------------------------|
C ---------------- F O R W A R D    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      IF (TO_index == +1) THEN
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
          dLambda_dSigma_FWD = (3.0_8/2.0_8)*
     & ((dSigma_Eff_Dev_over_Sigma_VM_dSigma*H_current)+
     & MATMUL((Sigma_Eff_Dev/Sigma_Eff_VM),dH_current_dSigma))
C: Calculate the Vectorized Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
              dLambda_dSigma_FWD(1,:) = dLambda_dSigma_FWD(1,:)
              dLambda_dSigma_FWD(2,:) = dLambda_dSigma_FWD(2,:)
              dLambda_dSigma_FWD(3,:) = dLambda_dSigma_FWD(3,:)
              dLambda_dSigma_FWD(4,:) = dLambda_dSigma_FWD(4,:)*2.0_8
              dLambda_dSigma_FWD(5,:) = dLambda_dSigma_FWD(5,:)*2.0_8
              dLambda_dSigma_FWD(6,:) = dLambda_dSigma_FWD(6,:)*2.0_8
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Transformation Strain [6 x 6]
          dLambda_det_FWD = MATMUL(dLambda_dSigma_FWD,dSigma_det)
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Martensitic Volume Fraction [6 x 1]
          dLambda_dMVF_FWD = MATMUL(dLambda_dSigma_FWD,dSigma_dMVF)
C: Partial Derivative Critical Thermodynamic Driving Force with respect to Stress [1 x 6]
          dY_dSigma_FWD = (D*MATMUL(TRANSPOSE(Sigma_Eff),
     & dLambda_dSigma_FWD))+(D*TRANSPOSE(Lambda_FWD))
C: Partial Derivative Critical Thermodynamic Driving Force with respect to Transformation Strain [1 x 6]
          dY_det_FWD = MATMUL(dY_dSigma_FWD,dSigma_det)
C: Partial Derivative Critical Thermodynamic Driving Force with respect to MVF [SCALAR]
          dY_dMVF_FWD = MINVAL(MATMUL(dY_dSigma_FWD,dSigma_dMVF))
C: Calculate the Derivative of R function with respect to the Transformation Strain [6 x 6]
          dR_det_FWD = 
     & -I_6x6+(dLambda_det_FWD*(MVF_current-MVF_previous))
C: Calculate the Derivative of R function with respect to the Martensitic Volume Fraction [6 x 1]
          dR_dMVF_FWD = 
     & (dLambda_dMVF_FWD*(MVF_current-MVF_previous))+Lambda_FWD
C: Calculate the Derivative of Transformation function with respect to the Transformation Strain [1 x 6]
          dPhi_det_FWD = MATMUL(TRANSPOSE(Sigma_Eff),dLambda_det_FWD)+
     & MATMUL(TRANSPOSE(Lambda_FWD),dSigma_det)+
     & MATMUL(MATMUL(TRANSPOSE(Stress_current),
     & TRANSPOSE(Delta_S*S_Iso)),dSigma_det)+MATMUL(
     & TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial),dSigma_det)-
     & dY_det_FWD
C: Calculate the Derivative of Transformation function with respect to the MVF [SCALAR]
          dPhi_dMVF_FWD =
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),dLambda_dMVF_FWD))+
     & MINVAL(MATMUL(TRANSPOSE(Lambda_FWD),dSigma_dMVF))+
     & MINVAL(MATMUL(MATMUL(TRANSPOSE(Stress_current),
     & TRANSPOSE(Delta_S*S_Iso)),dSigma_dMVF))+
     & MINVAL(MATMUL(TRANSPOSE(Delta_a*alpha_Coef)*
     & (TEMP+DTEMP-T_initial),dSigma_dMVF))-dft_dMVF_FWD-
     & dY_dMVF_FWD
          
C:-----------------------------------------------------------------------!
          dR_det = dR_det_FWD
          dR_dMVF = dR_dMVF_FWD
          dPhi_det = dPhi_det_FWD
          dPhi_dMVF = dPhi_dMVF_FWD
          
          !IF (Debug == 1) THEN
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_dSigma_FWD"
          !    DO i=1,NTENS
          !        WRITE(*,11)dLambda_dSigma_FWD(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_det_FWD"
          !    DO i=1,NTENS
          !        WRITE(*,11)dLambda_det_FWD(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_dMVF_FWD"
          !    DO i=1,NTENS
          !        WRITE(*,*)dLambda_dMVF_FWD(i,1)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dY_dSigma_FWD"
          !    WRITE(*,11)dY_dSigma_FWD(1,:)
          !    WRITE(*,*)
          !    WRITE(*,*)"dR_det"
          !    DO i=1,NTENS
          !        WRITE(*,11)dR_det(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dR_dMVF"
          !    DO i=1,NTENS
          !        WRITE(*,*)dR_dMVF(i,1)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dPhi_det_FWD"
          !    DO i=1,NTENS
          !        WRITE(*,*)dPhi_det_FWD(1,i)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dPhi_dMVF_FWD"
          !    WRITE(*,*)dPhi_dMVF_FWD
          !    WRITE(*,*)
          !ENDIF
          
C:-----------------------------------------------------------------------!
      ENDIF
      
C--------------------------------------------------------------------------------|
C ---------------- R E V E R S E    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      IF (TO_index == -1) THEN
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
          dLambda_dSigma_REV = 0.0_8
C: Calculate the Derivative of the Transformation Tensor with respect to the Transformation Strain [6 x 6]
          dLambda_det_REV = 0.0_8
C: Calculate the Derivative of the Transformation Tensor with respect to the Martensitic Volume Fraction [6 x 1]
          dLambda_dMVF_REV = 0.0_8
C: Partial Derivative Critical Thermodynamic Driving Force with respect to Stress [1 x 6]
          dY_dSigma_REV = (D*TRANSPOSE(Lambda_REV))
C: Partial Derivative Critical Thermodynamic Driving Force with respect to Transformation Strain [1 x 6]
          dY_det_REV = MATMUL(dY_dSigma_REV,dSigma_det)
C: Partial Derivative Critical Thermodynamic Driving Force with respect to MVF [SCALAR]
          dY_dMVF_REV = MINVAL(MATMUL(dY_dSigma_REV,dSigma_dMVF))
C: Calculate the Derivative of R function with respect to the Transformation Strain [6 x 6]
          dR_det_REV = -I_6x6
C: Calculate the Derivative of R function with respect to the Martensitic Volume Fraction [6 x 1]
          dR_dMVF_REV = Lambda_REV
C: Calculate the Derivative of Transformation function with respect to the Transformation Strain [1 x 6]
          dPhi_det_REV = -MATMUL(TRANSPOSE(Lambda_REV),dSigma_det)-
     & MATMUL(MATMUL(TRANSPOSE(Stress_current),
     & TRANSPOSE(Delta_S*S_Iso)),dSigma_det)-MATMUL(
     & TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial),dSigma_det)-
     & dY_det_REV
C: Calculate the Derivative of Transformation function with respect to the MVF [SCALAR]  
          dPhi_dMVF_REV = 
     & -MINVAL(MATMUL(TRANSPOSE(Lambda_REV),dSigma_dMVF))-
     & MINVAL(MATMUL(MATMUL(TRANSPOSE(Stress_current),
     & TRANSPOSE(Delta_S*S_Iso)),dSigma_dMVF))-MINVAL(MATMUL(TRANSPOSE(
     & Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial),dSigma_dMVF))+
     & dft_dMVF_REV-dY_dMVF_REV
C:-----------------------------------------------------------------------!
          dR_det = dR_det_REV
          dR_dMVF = dR_dMVF_REV
          dPhi_det = dPhi_det_REV
          dPhi_dMVF = dPhi_dMVF_REV
          
          !IF (Debug == 1) THEN
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_dSigma_REV"
          !    DO i=1,NTENS
          !        WRITE(*,11)dLambda_dSigma_REV(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_det_REV"
          !    DO i=1,NTENS
          !        WRITE(*,11)dLambda_det_REV(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dLambda_dMVF_REV"
          !    DO i=1,NTENS
          !        WRITE(*,*)dLambda_dMVF_REV(i,1)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dY_dSigma_REV"
          !    WRITE(*,11)dY_dSigma_REV(1,:)
          !    WRITE(*,*)
          !    WRITE(*,*)"dR_det"
          !    DO i=1,NTENS
          !        WRITE(*,11)dR_det(i,:)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dR_dMVF"
          !    DO i=1,NTENS
          !        WRITE(*,*)dR_dMVF(i,1)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dPhi_det_REV"
          !    DO i=1,NTENS
          !        WRITE(*,*)dPhi_det_REV(1,i)
          !    ENDDO
          !    WRITE(*,*)
          !    WRITE(*,*)"dPhi_dMVF_REV"
          !    WRITE(*,*)dPhi_dMVF_REV
          !    WRITE(*,*)
          !ENDIF
          
C:-----------------------------------------------------------------------!
      ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE FE_STRESS_Calculation(STRESS,
     & STRAN,DSTRAN,TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Initial_State,SMA_State)
          
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRESS(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     & PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 Debug, i
      
C: Double Precision Variables
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31), Stress_current(NTENS,1)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Debug = INT(PROPS(38))
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
      
C----------------------------------------------------------------------------------------------|
C Provide the required Volumetric Heat Generation to   A B A Q U S   G L O B A L   S O L V E R |
C----------------------------------------------------------------------------------------------|
      
      DO i=1,NTENS
          STRESS(i) = Stress_current(i,1)
      ENDDO
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE FE_Derivatives(STATEV,DDSDDE,DDSDDT,DRPLDE,
     & DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,NTENS,NSTATV,PROPS,
     & NPROPS,Model_Parameters,Initial_State,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 RPL, DRPLDT, TIME, DTIME, TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     & DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     & PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER TO_index, TM_Coupling, Debug, i
      
C: Double Precision Variables
      REAL*8 Density, E_A, E_M, PR, alpha_A, alpha_M, ESH_A, ESH_M,
     & S_A, S_M, Delta_S, Delta_a, Delta_c, rDs0, D, rDu0, T_initial,
     & MVF_current, MVF_previous, c_current, H_current, Sigma_Eff_VM,
     & ft_FWD, ft_REV, dft_dMVF_FWD, dft_dMVF_REV, Pi_FWD, Pi_REV,
     & dpi_dTemperature_FWD, dpi_dTemperature_REV,
     & dpi_dMVF_FWD, dpi_dMVF_REV,
     & dPhi_dTemperature_FWD, dPhi_dTemperature_REV,
     & dPhi_dMVF_FWD, dPhi_dMVF_REV,
     & dHeat_dTemperature_FWD, dHeat_dTemperature_REV,
     & dHeat_dMVF_FWD, dHeat_dMVF_REV,
     & dRPL_dTemperature_FWD, dRPL_dTemperature_REV,
     & dRPL_dTemperature_TE
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31),
     & et_initial(NTENS,1), et_current(NTENS,1),
     & Stress_previous(NTENS,1), Stress_Current(NTENS,1),
     & S_Iso(NTENS,NTENS), alpha_Coef(NTENS,1),
     & S_current(NTENS,NTENS), a_current(NTENS,1),
     & S_Inverted(NTENS,NTENS),
     & Sigma_Eff(NTENS,1), Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS),
     & dH_current_dSigma(1,NTENS),
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1),
     & dLambda_dSigma_FWD(NTENS,NTENS), dLambda_dSigma_REV(NTENS,NTENS),
     & dpi_dSigma_FWD(1,NTENS), dpi_dSigma_REV(1,NTENS),
     & dY_dSigma_FWD(1,NTENS), dY_dSigma_REV(1,NTENS),
     & dPhi_dSigma_FWD(1,NTENS), dPhi_dSigma_REV(1,NTENS),
     & OC1(NTENS,NTENS), OC2(NTENS,1),
     & IC_01_FWD(NTENS,NTENS), IC_01_REV(NTENS,NTENS),
     & IC_02_FWD(NTENS,NTENS), IC_02_REV(NTENS,NTENS),
     & dSigma_dEpsilon_FWD(NTENS,NTENS),
     & dSigma_dEpsilon_REV(NTENS,NTENS),
     & dHeat_dSigma_FWD(1,NTENS), dHeat_dSigma_REV(1,NTENS),
     & IC_03_FWD(NTENS,1), IC_03_REV(NTENS,1),
     & dSigma_dTemperature_FWD(NTENS,1),
     & dSigma_dTemperature_REV(NTENS,1),
     & dRPL_dEpsilon_FWD(1,NTENS), dRPL_dEpsilon_REV(1,NTENS),
     & dSigma_dEpsilon_TE(NTENS,NTENS), dSigma_dTemperature_TE(NTENS,1),
     & dRPL_dEpsilon_TE(1,NTENS)
      
11    FORMAT(6(E24.8E3))
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Density = PROPS(1) ! Density of the material
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4) ! Poisson's Ratio
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
      TM_Coupling = INT(PROPS(36))
      Debug = INT(PROPS(38))
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Model Parameters |
C---------------------------------------------------------------------------|
      
      rDs0 = Model_Parameters(1)
      D = Model_Parameters(2)
      rDu0 = Model_Parameters(6)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State of the Material |
C----------------------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      DO i=1,NTENS
          et_initial(i,1) = Initial_State(2+i)
      ENDDO
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      TO_index = INT(SMA_State(14))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_previous = STATEV(1)
      DO i=1,NTENS
          Stress_previous(i,1) = STATEV(7+i)
      ENDDO
      
C----------------------------------------------------------------------|
C Define any Auxiliary Matrices required for the particular Subroutine |
C----------------------------------------------------------------------|
      
C: Define the Compliance Matrix Coefficients
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Define the Thermal Expansion Matrix Coefficients
      alpha_Coef(1,1) = 1.0_8
      alpha_Coef(2,1) = 1.0_8
      alpha_Coef(3,1) = 1.0_8
      alpha_Coef(4,1) = 0.0_8
      alpha_Coef(5,1) = 0.0_8
      alpha_Coef(6,1) = 0.0_8
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Material Properties
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
C: Invert the Compliance to Calculate the Stiffness Matrix
      CALL Matrix_Inversion(S_current,S_Inverted,NTENS)
C: Call the Subroutine for the Calculation of Effective Material Properties to redefine S_current
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
      
C--------------------------------------------------------------------------------|
C Perform the required Calculations based on the Direction of the Transformation |
C--------------------------------------------------------------------------------|
      
      IF ((TO_index == +1).OR.(TO_index == -1)) THEN
          
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
          
C: Call the Subroutine for the Calculation of Effective Stress and its Derivatives
          CALL Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
C: Call the Subroutine for the Calculation of H_current and its Derivatives
          CALL H_current_Calculation(STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & H_current,dH_current_dSigma)
C: Call the Subroutine for the Calculation of Transformation Tensor
          CALL Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
C: Call the Subroutine for the Calculation of Hardening Function and its Derivatives
          CALL Hardening_Function(NTENS,PROPS,NPROPS,Model_Parameters,
     & SMA_State,ft_FWD,ft_REV,dft_dMVF_FWD,dft_dMVF_REV)
          
C--------------------------------------------------------------------------------|
C ---------------- F O R W A R D    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
          
          IF (TO_index == +1) THEN
              
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
              dLambda_dSigma_FWD = (3.0_8/2.0_8)*
     & ((dSigma_Eff_Dev_over_Sigma_VM_dSigma*H_current)+
     & MATMUL((Sigma_Eff_Dev/Sigma_Eff_VM),dH_current_dSigma))
C: Calculate the Vectorized Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
              dLambda_dSigma_FWD(1,:) = dLambda_dSigma_FWD(1,:)
              dLambda_dSigma_FWD(2,:) = dLambda_dSigma_FWD(2,:)
              dLambda_dSigma_FWD(3,:) = dLambda_dSigma_FWD(3,:)
              dLambda_dSigma_FWD(4,:) = dLambda_dSigma_FWD(4,:)*2.0_8
              dLambda_dSigma_FWD(5,:) = dLambda_dSigma_FWD(5,:)*2.0_8
              dLambda_dSigma_FWD(6,:) = dLambda_dSigma_FWD(6,:)*2.0_8
C: Calculate the Total Thermodynamic Force [SCALAR]
              Pi_FWD = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_FWD))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_FWD
              
C: Calculate the Partial Derivative of Total Transformation Energy with respect to Sigma [1 x 6]
              dpi_dSigma_FWD = 
     & MATMUL(TRANSPOSE(Sigma_Eff),dLambda_dSigma_FWD)+
     & TRANSPOSE(Lambda_FWD)+
     & MATMUL(TRANSPOSE(Stress_current),TRANSPOSE(Delta_S*S_Iso))+
     & (TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial))
C: Calculate the Partial Derivative of Total Transformation Energy with respect to Temperature [SCALAR]
              dpi_dTemperature_FWD = 
     & MINVAL(MATMUL(TRANSPOSE(Stress_current),(Delta_a*alpha_Coef)))+
     & (Density*Delta_c*LOG((TEMP+DTEMP)/T_initial))+rDs0
C: Calculate the Partial Derivative of Total Transformation Energy with respect to Martensitic Volume Fraction [SCALAR]
              dpi_dMVF_FWD = -dft_dMVF_FWD
C: Calculate the Partial Derivative Critical Thermodynamic Driving Force with respect to Stress [1 x 6]
              dY_dSigma_FWD = (D*MATMUL(TRANSPOSE(Sigma_Eff),
     & dLambda_dSigma_FWD))+(D*TRANSPOSE(Lambda_FWD))
C: Calculate the Partial Derivative of Transformation Function with respect to Sigma [1 x 6]
              dPhi_dSigma_FWD = dpi_dSigma_FWD-dY_dSigma_FWD
C: Calculate the Partial Derivative of Transformation Function with respect to Temperature [SCALAR]
              dPhi_dTemperature_FWD = dpi_dTemperature_FWD
C: Calculate the Partial Derivative of Transformation Function with respect to Martensitic Volume Fraction [SCALAR]
              dPhi_dMVF_FWD = dpi_dMVF_FWD
              
C:   |------------------------------------------------------------------|
              
C: Calculate the Quantity OC1 [6 x 6]
              OC1 = S_current+((MVF_current-MVF_previous)*
     & dLambda_dSigma_FWD)
C: Calculate the Quantity OC2 [6 x 1]
              OC2 = MATMUL((Delta_S*S_Iso),Stress_current)+
     & ((Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial))+Lambda_FWD
              
C----------------------------------------|
C Calculate the Jacobian of the Material |
C----------------------------------------|
              
C: Perform the Intermediate Calculation 01 [6 x 6]
              IC_01_FWD = 
     & (OC1-MATMUL(OC2,(dPhi_dSigma_FWD/dPhi_dMVF_FWD)))
C: Perform the Intermediate Calculation 02 [6 x 6]
              CALL Matrix_Inversion(IC_01_FWD,IC_02_FWD,NTENS)
C: Calculate the Partial Derivative of Stress with respect to Total Strain [6 x 6]
              dSigma_dEpsilon_FWD = IC_02_FWD
C: Perform the Intermediate Calculation 03 [6 x 1]
              IC_03_FWD = 
     & (a_current-(OC2*(dPhi_dTemperature_FWD/dPhi_dMVF_FWD)))
C: Calculate the Variation of Stress Increments with respect to the Temperature [6 x 1]
              dSigma_dTemperature_FWD = -MATMUL(IC_02_FWD,IC_03_FWD)
              
C-------------------------------------------|
C Calculate the Partial Derivatives of Heat |
C-------------------------------------------|
              
C: Calulate the Partial Derivative of Heat with respect to Sigma [1 x 6]
              dHeat_dSigma_FWD = 
     & ((dpi_dSigma_FWD-(TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP)))*
     & ((MVF_current-MVF_previous)/DTIME))-
     & ((TRANSPOSE(a_current)*(TEMP+DTEMP))/DTIME)
C: Calulate the Partial Derivative of Heat with respect to Sigma [SCALAR]
              dHeat_dTemperature_FWD =((dpi_dTemperature_FWD-
     & MINVAL(MATMUL(TRANSPOSE(Stress_current),(Delta_a*alpha_Coef)))-
     & (Density*Delta_c*(LOG((TEMP+DTEMP)/T_initial)+1))-rDs0)*
     & ((MVF_current-MVF_previous)/DTIME))-
     & MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & a_current))
C: Calulate the Partial Derivative of Heat with respect to MVF [SCALAR]
              dHeat_dMVF_FWD = 
     & ((Pi_FWD-(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & (Delta_a*alpha_Coef*(TEMP+DTEMP)))))-(Density*Delta_c*
     & (TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))-(rDs0*(TEMP+DTEMP)))/
     & DTIME)+(dpi_dMVF_FWD*((MVF_current-MVF_previous)/DTIME))-
     & (MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & (Delta_a*alpha_Coef)*(TEMP+DTEMP))))
              
C:   |------------------------------------------------------------------|
              
C: Calculate the Partial Derivative of the Volumetric Heat Generation per Unit Time with respect to Total Strain [1 x 6]
              dRPL_dEpsilon_FWD = 
     & MATMUL(dHeat_dSigma_FWD,dSigma_dEpsilon_FWD)-
     & (MATMUL(dHeat_dMVF_FWD*(dPhi_dSigma_FWD/dPhi_dMVF_FWD),
     & dSigma_dEpsilon_FWD))
C: Calculate the Partial Derivative of Volumetric Heat Generation per Unit Time with respect to Temperature [SCALAR]
              dRPL_dTemperature_FWD = 
     & MINVAL(MATMUL(dHeat_dSigma_FWD,dSigma_dTemperature_FWD))-
     & (dHeat_dMVF_FWD*(MINVAL(MATMUL((dPhi_dSigma_FWD/dPhi_dMVF_FWD),
     & dSigma_dTemperature_FWD))+(dPhi_dTemperature_FWD/dPhi_dMVF_FWD)))
     & +dHeat_dTemperature_FWD
              
C:   |------------------------------------------------------------------|
              
              !IF (Debug == 1) THEN
              !    WRITE(*,*)
              !    WRITE(*,*)"dHeat_dSigma_FWD"
              !    DO i = 1,NTENS
              !        WRITE(*,*)dHeat_dSigma_FWD(1,i)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dSigma_dEpsilon_FWD"
              !    DO i = 1,NTENS
              !        WRITE(*,11)dHeat_dSigma_FWD(i,:)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dHeat_dMVF_FWD",dHeat_dMVF_FWD
              !    WRITE(*,*)
              !    WRITE(*,*)"dPhi_dSigma_FWD"
              !    DO i = 1,NTENS
              !        WRITE(*,11)dPhi_dSigma_FWD(i,:)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dPhi_dMVF_FWD",dPhi_dMVF_FWD
              !    WRITE(*,*)
              !ENDIF
              
C:   |------------------------------------------------------------------|
              
          ENDIF
          
C--------------------------------------------------------------------------------|
C ---------------- R E V E R S E    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
          
          IF (TO_index == -1) THEN
C: Calculate the Partial Derivative of the Transformation Tensor with respect to the Stress [6 x 6]
              dLambda_dSigma_REV = 0.0_8
C: Calculate the Total Thermodynamic Force [SCALAR]
              Pi_REV = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_REV))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_REV
C: Calculate the Partial Derivative of Total Transformation Energy with respect to Sigma [1 x 6]
              dpi_dSigma_REV = 
     & MATMUL(TRANSPOSE(Sigma_Eff),dLambda_dSigma_REV)+
     & TRANSPOSE(Lambda_REV)+
     & MATMUL(TRANSPOSE(Stress_current),TRANSPOSE(Delta_S*S_Iso))+
     & (TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial))
C: Calculate the Partial Derivative of Total Transformation Energy with respect to Temperature [SCALAR]
              dpi_dTemperature_REV = 
     & MINVAL(MATMUL(TRANSPOSE(Stress_current),(Delta_a*alpha_Coef)))+
     & (Density*Delta_c*LOG((TEMP+DTEMP)/T_initial))+rDs0
C: Calculate the Partial Derivative Critical Thermodynamic Driving Force with respect to Stress [1 x 6]
              dY_dSigma_REV = (D*MATMUL(TRANSPOSE(Sigma_Eff),
     & dLambda_dSigma_REV))+(D*TRANSPOSE(Lambda_REV))
C: Calculate the Partial Derivative of Total Transformation Energy with Respect to Martensitic Volume Fraction [SCALAR]
              dpi_dMVF_REV = -dft_dMVF_REV
C: Calculate the Partial Derivative of Transformation Function with respect to Sigma [1 x 6]
              dPhi_dSigma_REV = -dpi_dSigma_REV-dY_dSigma_REV
C: Calculate the Partial Derivative of Transformation Function with respect to Temperature [SCALAR]
              dPhi_dTemperature_REV = -dpi_dTemperature_REV
C: Calculate the Partial Derivative of Transformation Function with respect to Martensitic Volume Fraction [SCALAR]
              dPhi_dMVF_REV = -dpi_dMVF_REV
              
C:   |------------------------------------------------------------------|
              
C: Calculate the Quantity OC1 [6 x 6]
              OC1 = S_current+((MVF_current-MVF_previous)*
     & dLambda_dSigma_REV)
              
C: Calculate the Quantity OC2 [6 x 1]
              OC2 = MATMUL((Delta_S*S_Iso),Stress_current)+
     & ((Delta_a*alpha_Coef)*(TEMP+DTEMP-T_initial))+Lambda_REV
              
C----------------------------------------|
C Calculate the Jacobian of the Material |
C----------------------------------------|
              
C: Perform the Intermediate Calculation 01 [6 x 6]
              IC_01_REV = 
     & (OC1-MATMUL(OC2,(dPhi_dSigma_REV/dPhi_dMVF_REV)))
C: Perform the Intermediate Calculation 02 [6 x 6]
              CALL Matrix_Inversion(IC_01_REV,IC_02_REV,NTENS)
C: Calculate the Partial Derivative of Stress with respect to Total Strain [6 x 6]
              dSigma_dEpsilon_REV = IC_02_REV
C: Perform the Intermediate Calculation 03 [6 x 1]
              IC_03_REV = 
     & (a_current-(OC2*(dPhi_dTemperature_REV/dPhi_dMVF_REV)))
C: Calculate the Variation of Stress Increments with respect to the Temperature [6 x 1]
              dSigma_dTemperature_REV = -MATMUL(IC_02_REV,IC_03_REV)
              
C-------------------------------------------|
C Calculate the Partial Derivatives of Heat |
C-------------------------------------------|
              
C: Calulate the Partial Derivative of Heat with respect to Sigma [1 x 6]
              dHeat_dSigma_REV = 
     & ((dpi_dSigma_REV-(TRANSPOSE(Delta_a*alpha_Coef)*(TEMP+DTEMP)))*
     & ((MVF_current-MVF_previous)/DTIME))-
     & ((TRANSPOSE(a_current)*(TEMP+DTEMP))/DTIME)
C: Calulate the Partial Derivative of Heat with respect to Sigma [SCALAR]
              dHeat_dTemperature_REV = ((dpi_dTemperature_REV-
     & MINVAL(MATMUL(TRANSPOSE(Stress_current),(Delta_a*alpha_Coef)))-
     & (Density*Delta_c*(LOG((TEMP+DTEMP)/T_initial)+1))-rDs0)*
     & ((MVF_current-MVF_previous)/DTIME))-
     & MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & a_current))
C: Calulate the Partial Derivative of Heat with respect to MVF [SCALAR]
              dHeat_dMVF_REV = 
     & ((Pi_REV-(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & (Delta_a*alpha_Coef*(TEMP+DTEMP)))))-(Density*Delta_c*
     & (TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))-(rDs0*(TEMP+DTEMP)))/
     & DTIME)+(dpi_dMVF_REV*((MVF_current-MVF_previous)/DTIME))-
     & (MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & (Delta_a*alpha_Coef)*(TEMP+DTEMP))))
              
C:   |------------------------------------------------------------------|
              
C: Calculate the Partial Derivative of the Volumetric Heat Generation per Unit Time with respect to Total Strain [1 x 6]
              dRPL_dEpsilon_REV = 
     & MATMUL(dHeat_dSigma_REV,dSigma_dEpsilon_REV)-
     & (MATMUL(dHeat_dMVF_REV*(dPhi_dSigma_REV/dPhi_dMVF_REV),
     & dSigma_dEpsilon_REV))
C: Calculate the Partial Derivative of Volumetric Heat Generation per Unit Time with respect to Temperature [SCALAR]
              dRPL_dTemperature_REV = 
     & MINVAL(MATMUL(dHeat_dSigma_REV,dSigma_dTemperature_REV))-
     & (dHeat_dMVF_REV*(MINVAL(MATMUL((dPhi_dSigma_REV/dPhi_dMVF_REV),
     & dSigma_dTemperature_REV))+(dPhi_dTemperature_REV/dPhi_dMVF_REV)))
     & +dHeat_dTemperature_REV
              
C:   |------------------------------------------------------------------|
              
              !IF (Debug == 1) THEN
              !    WRITE(*,*)
              !    WRITE(*,*)"dHeat_dSigma_REV"
              !    DO i = 1,NTENS
              !        WRITE(*,*)dHeat_dSigma_REV(1,i)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dSigma_dEpsilon_REV"
              !    DO i = 1,NTENS
              !        WRITE(*,11)dHeat_dSigma_REV(i,:)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dHeat_dMVF_REV",dHeat_dMVF_REV
              !    WRITE(*,*)
              !    WRITE(*,*)"dPhi_dSigma_REV"
              !    DO i = 1,NTENS
              !        WRITE(*,11)dPhi_dSigma_REV(i,:)
              !    ENDDO
              !    WRITE(*,*)
              !    WRITE(*,*)"dPhi_dMVF_REV",dPhi_dMVF_REV
              !    WRITE(*,*)
              !ENDIF
              
C:   |------------------------------------------------------------------|
              
          ENDIF
          
      ELSEIF (TO_index == 0) THEN
C: No Transformation Takes Place - Thermo-elastic Response
          
C: Calculate the Partial Derivative of Stress with respect to the Total Strain [6 x 6]
          dSigma_dEpsilon_TE = S_Inverted
C: Calculate the Partial Derivative of Stress with respect to the Temperature [6 x 1]
          dSigma_dTemperature_TE = -MATMUL(S_Inverted,a_current)
C: Calculate the Partial Derivative of Heat with respect to the Total Strain [1 x 6]
          dRPL_dEpsilon_TE =
     & -(MATMUL(TRANSPOSE(a_current)*(TEMP+DTEMP),S_Inverted)/DTIME)
C: Calculate the Partial Derivative of Heat with respect to the Temperature [SCALAR]
          dRPL_dTemperature_TE = (MINVAL(MATMUL(MATMUL(
     & TRANSPOSE(a_current),S_Inverted),a_current*(TEMP+DTEMP)))/DTIME)-
     & (MINVAL(MATMUL(TRANSPOSE(a_current),
     & (Stress_current-Stress_previous)))/DTIME)
          
      ELSE
          
C: Print the proper ERROR Message
          WRITE(*,*)"E R R O R", "Transformation Occurence N/A"
          
      ENDIF
      
C: ------------------------------------------------------------------------------------- :|
C: Provide the required Partial Derivatives to   A B A Q U S   G L O B A L   S O L V E R
C: ------------------------------------------------------------------------------------- :|
      
C--------------------------------------------------------------------------------|
C ---------------- F O R W A R D    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      
      IF (TO_index == +1) THEN
C: Provide the Material Jacobian Matrix
          DO i=1,NTENS
              DDSDDE(i,:) = dSigma_dEpsilon_FWD(i,:)
          ENDDO
C: Provide the Variation of the Stress Increments with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DDSDDT(i) = dSigma_dTemperature_FWD(i,1)
              ENDDO
          ELSE
              DDSDDT = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Total Strain Increments
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DRPLDE(i) = dRPL_dEpsilon_FWD(1,i)
              ENDDO
          ELSE
              DRPLDE = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DRPLDT = dRPL_dTemperature_FWD
          ELSE
              DRPLDT = 0.0_8
          ENDIF
      ENDIF
      
C--------------------------------------------------------------------------------|
C ---------------- R E V E R S E    T R A N S F O R M A T I O N ---------------- |
C--------------------------------------------------------------------------------|
      
      IF (TO_index == -1) THEN
C: Provide the Material Jscobian Matrix
          DO i=1,NTENS
              DDSDDE(i,:) = dSigma_dEpsilon_REV(i,:)
          ENDDO
C: Provide the Variation of the Stress Increments with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DDSDDT(i) = dSigma_dTemperature_REV(i,1)
              ENDDO
          ELSE
              DDSDDT = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Total Strain Increments
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DRPLDE(i) = dRPL_dEpsilon_REV(1,i)
              ENDDO
          ELSE
              DRPLDE = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DRPLDT = dRPL_dTemperature_REV
          ELSE
              DRPLDT = 0.0_8
          ENDIF
      ENDIF
      
C--------------------------------------------------------------------------------|
C--------------------------------------------------------------------------------|
      
      IF (TO_index == 0) THEN
C: Provide the Material Jscobian Matrix
          DO i=1,NTENS
              DDSDDE(i,:) = dSigma_dEpsilon_TE(i,:)
          ENDDO
C: Provide the Variation of the Stress Increments with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DDSDDT(i) = dSigma_dTemperature_TE(i,1)
              ENDDO
          ELSE
              DDSDDT = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Total Strain Increments
          IF (TM_Coupling == 1) THEN
              DO i=1,NTENS
                  DRPLDE(i) = dRPL_dEpsilon_TE(1,i)
              ENDDO
          ELSE
              DRPLDE = 0.0_8
          ENDIF
C: Provide the Variation of the Volumetric Heat Generation per Unit Time with respect to the Temperature
          IF (TM_Coupling == 1) THEN
              DRPLDT = dRPL_dTemperature_TE
          ELSE
              DRPLDT = 0.0_8
          ENDIF
      ENDIF
      
      !IF (Debug == 1) THEN
      !    WRITE(*,*)
      !    WRITE(*,*)"Transformation Index ->",TO_index
      !    WRITE(*,*)
      !    WRITE(*,*)"DDSDDE"
      !    DO i = 1,NTENS
      !        WRITE(*,11)DDSDDE(i,:)
      !    ENDDO
      !    WRITE(*,*)
      !    WRITE(*,*)"DDSDDT"
      !    DO i = 1,NTENS
      !        WRITE(*,*)DDSDDT(i)
      !    ENDDO
      !    WRITE(*,*)
      !    WRITE(*,*)"DRPLDE"
      !    DO i = 1,NTENS
      !        WRITE(*,*)DRPLDE(i)
      !    ENDDO
      !    WRITE(*,*)
      !ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE FE_RPL_Calculation(STRESS,STATEV,RPL,
     & STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,NTENS,NSTATV,PROPS,NPROPS,
     & Model_Parameters,Initial_State,SMA_State)
       
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 NTENS, NSTATV, NPROPS
      
C: Double Precision Variables
      REAL*8 RPL, TIME, DTIME, TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRESS(NTENS), STATEV(NSTATV), STRAN(NTENS), DSTRAN(NTENS),
     & PROPS(NPROPS)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integer Variables
      INTEGER*4 TM_Coupling, Debug, TO_index, i
      
C: Double Precision Variables
      REAL*8 Density, E_A, E_M, PR, alpha_A, alpha_M, ESH_A, ESH_M,
     & S_A, S_M, Delta_S, Delta_a, Delta_c,
     & rDs0, rDu0, T_initial, MVF_previous, MVF_current, c_current,
     & Sigma_Eff_VM, ft_FWD, ft_REV, dft_dMVF_FWD, dft_dMVF_REV,
     & Pi_FWD, Pi_REV,
     & VHG
      
C: Double Precision Arrays
      REAL*8 Model_Parameters(7), Initial_State(8), SMA_State(31),
     & S_Iso(NTENS,NTENS), S_current(NTENS,NTENS),
     & alpha_Coef(NTENS,1), a_current(NTENS,1),
     & Stress_previous(NTENS,1), Stress_Current(NTENS,1),
     & Sigma_Eff(NTENS,1), Sigma_Eff_Dev(NTENS,1),
     & dSigma_Eff_VM_dSigma(1,NTENS),
     & dSigma_Eff_Dev_over_Sigma_VM_dSigma(NTENS,NTENS),
     & Lambda_FWD(NTENS,1), Lambda_REV(NTENS,1)
      
C------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Material Properties |
C------------------------------------------------------------------------------|
      
      Density = PROPS(1) ! Density of the material
      E_A = PROPS(2) ! Elastic Modulus of Austenite
      E_M = PROPS(3) ! Elastic Modulus of Martensite
      PR = PROPS(4) ! Poisson's Ratio
      alpha_A = PROPS(6) ! Thermal Expansion Coefficient of Austenite
      alpha_M = PROPS(7) ! Thermal Expansion Coefficient of Martensite
      ESH_A = PROPS(11) ! Effective Specific Heat of Austenite
      ESH_M = PROPS(12) ! Effective Specific Heat of Martensite
      TM_Coupling = INT(PROPS(36))
      Debug = INT(PROPS(38))
      
C: Calculate the Compliance values for Austenite and Martensite
      S_A = 1.0_8/E_A
      S_M = 1.0_8/E_M
C: Calculate the Difference between the Compliance of Martensite and Austenite - Delta S
      Delta_S = S_M-S_A
C: Calculate the Difference between the Thermal Expansion Coefficient of Martensite and Austenite - Delta alpha
      Delta_a = alpha_M-alpha_A
C: Calculate the Difference between the Effective Specific Heat of Martensite and Austenite - Delta c
      Delta_c = ESH_M-ESH_A
      
C---------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Model Parameters |
C---------------------------------------------------------------------------|
      
      rDs0 = Model_Parameters(1)
      rDu0 = Model_Parameters(6)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the Initial State of the Material |
C----------------------------------------------------------------------------------------|
      
      T_initial = Initial_State(1)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      TO_index = INT(SMA_State(14))
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the values of the State Variables |
C----------------------------------------------------------------------------------------|
      
      MVF_previous = STATEV(1)
      DO i=1,NTENS
          Stress_previous(i,1) = STATEV(7+i)
      ENDDO
      
C----------------------------------------------------------------------------|
C Perform the required Calculations based on the Phase State of the Material |
C----------------------------------------------------------------------------|
      
C: Define the Compliance Coefficients of a Generally Isotropic Elastic Material [6 x 6]
      S_Iso(1,:) = (/ 1.0_8, -PR, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(2,:) = (/ -PR, 1.0_8, -PR, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(3,:) = (/ -PR, -PR, 1.0_8, 0.0_8, 0.0_8, 0.0_8 /)
      S_Iso(4,:) = (/ 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8, 0.0_8 /)
      S_Iso(5,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR), 0.0_8 /)
      S_Iso(6,:) = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 2.0_8*(1.0_8+PR) /)
      
C: Define the Thermal Expansion Matrix Coefficients
      alpha_Coef(1,1) = 1.0_8
      alpha_Coef(2,1) = 1.0_8
      alpha_Coef(3,1) = 1.0_8
      alpha_Coef(4,1) = 0.0_8
      alpha_Coef(5,1) = 0.0_8
      alpha_Coef(6,1) = 0.0_8
      
C--------------------------------------------|
C Call the required user Defined Subroutines |
C--------------------------------------------|
      
C: Call the Subroutine for the Calculation of Effective Material Properties
      CALL Effective_Material_Properties(NTENS,PROPS,NPROPS,
     & SMA_State,
     & S_current,a_current,c_current)
C: Call the Subroutine for the Calculation of Stress
      CALL Stress_Calculation(STRAN,DSTRAN,TEMP,DTEMP,NTENS,
     & PROPS,NPROPS,Initial_State,SMA_State,
     & Stress_Current)
C: Call the Subroutine for the Calculation of Effective Stress and its Derivatives
      CALL Sigma_Effective_and_Derivatives(STRAN,DSTRAN,
     & TEMP,DTEMP,NTENS,PROPS,NPROPS,Initial_State,SMA_State,
     & Sigma_Eff,Sigma_Eff_Dev,Sigma_Eff_VM,
     & dSigma_Eff_VM_dSigma,dSigma_Eff_Dev_over_Sigma_VM_dSigma)
C: Call the Subroutine for the Calculation of Transformation Tensor
      CALL Transformation_Tensor(STATEV,STRAN,DSTRAN,TEMP,DTEMP,
     & NTENS,NSTATV,PROPS,NPROPS,Initial_State,SMA_State,
     & Lambda_FWD,Lambda_REV)
C: Call the Subroutine for the Calculation of Hardening Function and its Derivatives
      CALL Hardening_Function(NTENS,PROPS,NPROPS,Model_Parameters,
     & SMA_State,ft_FWD,ft_REV,dft_dMVF_FWD,dft_dMVF_REV)
      
C:   |------------------------------------------------------------------|
      
C: Calculate the Total Thermodynamic Force based on Transformation Direction
      
C: Forward Transformation
      IF (TO_index == +1) THEN
          Pi_FWD = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_FWD))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_FWD
      ENDIF
C: Reverse Transformation
      IF (TO_index == -1) THEN
          Pi_REV = 
     & MINVAL(MATMUL(TRANSPOSE(Sigma_Eff),Lambda_REV))+
     & ((1.0_8/2.0_8)*(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & MATMUL((Delta_S*S_Iso),Stress_current)))))+
     & (MINVAL(MATMUL(TRANSPOSE(Stress_current),Delta_a*alpha_Coef))*
     & (TEMP+DTEMP-T_initial))-
     & (Density*Delta_c*((TEMP+DTEMP-T_initial)-
     & ((TEMP+DTEMP)*LOG((TEMP+DTEMP)/T_initial))))+(rDs0*(TEMP+DTEMP))-
     & rDu0-ft_REV
      ENDIF
      
C:   |------------------------------------------------------------------|
      
C: Calculate the Volumetric Heat Generation based on Transformation Direction
      
      IF (TO_index == +1) THEN
          VHG = ((Pi_FWD-(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & (Delta_a*alpha_Coef)*(TEMP+DTEMP))))-
     & ((TEMP+DTEMP)*Density*Delta_c*LOG((TEMP+DTEMP)/T_initial))-
     & (rDs0*(TEMP+DTEMP)))*((MVF_current-MVF_previous)/DTIME))-
     & MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & (a_current)*(TEMP+DTEMP)))
      ENDIF
      
      IF (TO_index == -1) THEN
          VHG = ((Pi_REV-(MINVAL(MATMUL(TRANSPOSE(Stress_current),
     & (Delta_a*alpha_Coef)*(TEMP+DTEMP))))-
     & ((TEMP+DTEMP)*Density*Delta_c*LOG((TEMP+DTEMP)/T_initial))-
     & (rDs0*(TEMP+DTEMP)))*((MVF_current-MVF_previous)/DTIME))-
     & MINVAL(MATMUL(TRANSPOSE((Stress_current-Stress_previous)/DTIME),
     & (a_current)*(TEMP+DTEMP)))
      ENDIF
      
      IF (TO_index == 0) THEN
          VHG = 0.0_8
      ENDIF
      
C----------------------------------------------------------------------------------------------|
C Provide the required Volumetric Heat Generation to   A B A Q U S   G L O B A L   S O L V E R |
C----------------------------------------------------------------------------------------------|
      
      IF (TM_Coupling == 1) THEN
          RPL = VHG
      ELSE
          RPL = 0.0_8
      ENDIF
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Internal_State_Correction(NTENS,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER TO_index, MVF_correction_index,
     & Completion_FWD, Completion_REV, i
      
C: Double Precision Variables
      REAL*8 MVF_current
      
C: Double Precision Arrays
      REAL*8 SMA_State(31), et_current(NTENS,1)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      TO_index = INT(SMA_State(14))
      Completion_FWD = INT(SMA_State(17))
      Completion_REV = INT(SMA_State(18))
      MVF_correction_index = INT(SMA_State(31))
      
C--------------------------------------------------------------------------------------|
C Perform the required Corrections to the values of the Internal State of the Material |
C--------------------------------------------------------------------------------------|
      
      IF (MVF_correction_index == 1) THEN
          IF (MVF_current > 0.9999_8) THEN
C: Set the MVF to the Proper Boundary Value
              MVF_current = 0.9999_8
C: Update the Indeces regarding Transformation Completion
              Completion_FWD = 1
              Completion_REV = 0
          ELSEIF (MVF_current < 0.0001_8) THEN
C: Set the MVF to the Proper Boundary Value
              MVF_current = 0.0001_8
C: Update the Indeces regarding Transformation Completion
              Completion_FWD = 0
              Completion_REV = 1
C: Set the et_current to the Proper Boundary Values
              et_current = 0.0_8
          ENDIF
      ENDIF
      
      IF (TO_index == +1) THEN
          DO i=1,NTENS
              SMA_State(22+i) = et_current(i,1)/MVF_current
          ENDDO
      ENDIF
      
C: Update the State of the Material
      SMA_State(1) = MVF_current
      DO i=1,NTENS
          SMA_State(1+i) = et_current(i,1)
      ENDDO
      SMA_State(17) = Completion_FWD
      SMA_State(18) = Completion_REV
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Update_STATEVs(STRESS,STATEV,TEMP,DTEMP,
     & NTENS,NSTATV,Initial_State,SMA_State)
      
C:-----------------------------------------------------------------------!
      
      IMPLICIT NONE
      
C:   | ----------------------------------------------------------------- |
C:   |                       ABAQUS UMAT Variables                       |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 NTENS, NSTATV
      
C: Double Precision Variables
      REAL*8 TEMP, DTEMP
      
C: Double Precision Arrays
      REAL*8 STRESS(NTENS), STATEV(NSTATV)
      
C:   | ----------------------------------------------------------------- |
C:   |                 Subroutine USER Defined Variables                 |
C:   | ----------------------------------------------------------------- |
      
C: Integers
      INTEGER*4 i, Direction_FWD, Direction_REV,
     & Completion_FWD, Completion_REV,
     & Reversal_FWD_REV, Reversal_REV_FWD, TO_index
      
C: Double Precision Variables
      REAL*8 MVF_current, MVF_reversal_FWD_REV, MVF_reversal_REV_FWD,
     & Phi_FWD, Phi_REV
      
C: Double Precision Arrays
      REAL*8 Initial_State(8), SMA_State(31), et_current(NTENS,1),
     & Lambda_REV(NTENS,1)
      
C----------------------------------------------------------------------------------------|
C Assign the values to the required Variables based on the current State of the Material |
C----------------------------------------------------------------------------------------|
      
      MVF_current = SMA_State(1)
      DO i=1,NTENS
          et_current(i,1) = SMA_State(1+i)
      ENDDO
      TO_index = INT(SMA_State(14))
      Direction_FWD = INT(SMA_State(15))
      Direction_REV = INT(SMA_State(16))
      Completion_FWD = INT(SMA_State(17))
      Completion_REV = INT(SMA_State(18))
      Reversal_FWD_REV = INT(SMA_State(19))
      Reversal_REV_FWD = INT(SMA_State(20))
      MVF_reversal_FWD_REV = SMA_State(21)
      MVF_reversal_REV_FWD = SMA_State(22)
      DO i=1,NTENS
          Lambda_REV(i,1) = SMA_State(22+i)
      ENDDO
      Phi_FWD = SMA_State(29)
      Phi_REV = SMA_State(30)
      
C:   |------------------------------------------------------------------|
      
      STATEV(1) = MVF_current
      STATEV(2) = et_current(1,1)
      STATEV(3) = et_current(2,1)
      STATEV(4) = et_current(3,1)
      STATEV(5) = et_current(4,1)
      STATEV(6) = et_current(5,1)
      STATEV(7) = et_current(6,1)
      STATEV(8) = STRESS(1)
      STATEV(9) = STRESS(2)
      STATEV(10) = STRESS(3)
      STATEV(11) = STRESS(4)
      STATEV(12) = STRESS(5)
      STATEV(13) = STRESS(6)
      STATEV(14) = TO_index
      STATEV(15) = Direction_FWD
      STATEV(16) = Direction_REV
      STATEV(17) = Completion_FWD
      STATEV(18) = Completion_REV
      STATEV(19) = Reversal_FWD_REV
      STATEV(20) = Reversal_REV_FWD
      STATEV(21) = MVF_reversal_FWD_REV
      STATEV(22) = MVF_reversal_REV_FWD
      STATEV(23) = Lambda_REV(1,1)
      STATEV(24) = Lambda_REV(2,1)
      STATEV(25) = Lambda_REV(3,1)
      STATEV(26) = Lambda_REV(4,1)
      STATEV(27) = Lambda_REV(5,1)
      STATEV(28) = Lambda_REV(6,1)
      STATEV(29) = Phi_FWD
      STATEV(30) = Phi_REV
      STATEV(31) = Initial_State(1)
      
C--------------------------------------------------------------------------------|
C ------------------------ E N D    S U B R O U T I N E ------------------------ |
C--------------------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      SUBROUTINE Matrix_Inversion(a,c,n)
      
C ================================================================= !
C: Inverse matrix
C: Method: Based on Doolittle LU factorization for Ax=b
C: Alex G. December 2009
C:------------------------------------------------------------------!
C: input ...
C: a(n,n) - array of coefficients for matrix A
C: n      - dimension
C: output ...
C: c(n,n) - inverse matrix of A
C: comments ...
C: the original matrix a(n,n) will be destroyed during the calculation
C ================================================================= !
      
      IMPLICIT NONE
      
      INTEGER n
      DOUBLE PRECISION a(n,n), c(n,n)
      DOUBLE PRECISION L(n,n), U(n,n), b(n), d(n), x(n)
      DOUBLE PRECISION coeff
      INTEGER i, j, k
      
C:   |------------------------------------------------------------------|
      
C: Step 0: initialization for matrices L and U and b
C: Fortran 90/95 aloows such operations on matrices
      L=0.0_8
      U=0.0_8
      b=0.0_8
      
C:   |------------------------------------------------------------------|
      
C: Step 1: forward elimination
      DO k=1, n-1
          DO i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              DO j=k+1,n
                  a(i,j) = a(i,j)-coeff*a(k,j)
              ENDDO
          ENDDO
      ENDDO
      
C:   |------------------------------------------------------------------|
      
C: Step 2: prepare L and U matrices 
C: L matrix is a matrix of the elimination coefficient + the diagonal elements are 1.0
      DO i=1,n
          L(i,i) = 1.0
      ENDDO
      ! U matrix is the upper triangular part of A
      DO j=1,n
          DO i=1,j
              U(i,j) = a(i,j)
          ENDDO
      ENDDO
      
C:   |------------------------------------------------------------------|
      
C: Step 3: compute columns of the inverse matrix C
      DO k=1,n
          b(k)=1.0
          d(1) = b(1)
C: Step 3a: Solve Ld=b using the forward substitution
          DO i=2,n
              d(i)=b(i)
              DO j=1,i-1
                  d(i) = d(i) - L(i,j)*d(j)
              ENDDO
          ENDDO
C: Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          DO i = n-1,1,-1
              x(i) = d(i)
              DO j=n,i+1,-1
                  x(i)=x(i)-U(i,j)*x(j)
              ENDDO
              x(i) = x(i)/u(i,i)
          ENDDO
C: Step 3c: fill the solutions x(n) into column k of C
          DO i=1,n
              c(i,k) = x(i)
          ENDDO
          b(k)=0.0
      ENDDO
      
C:   |------------------------------------------------------------------|
      
      END SUBROUTINE
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
      
      
C------------------------------------------------------------------------|
C            D E F I N I T I O N S    of    F U N C T I O N S
C------------------------------------------------------------------------|
      
      
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
      
C: Call this Function to Calculate the Dot Product of Two Vectors
      FUNCTION VP(VA,VB,SIZE)
      
      IMPLICIT NONE
      
C: Integer Variables
      INTEGER*4 i, SIZE
      
C: Double Precision Variables
      REAL*8 VP
      
C: Double Precision Arrays
      REAL*8 VA(SIZE,1), VB(SIZE,1)
      
      VP = 0.0_8
      
      DO i=1,3
          VP = VP+(VA(i,1)*VB(i,1))
      ENDDO
      
      RETURN
      
      END FUNCTION VP
      
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!
C------------------------------------------------------------------------>
      !S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S-S!