!=======================================================================
! Below is the supporting code required by the seismic version of 
! UserPtfmLd
!=======================================================================

! Note: This code includes components derived from routines
! in RspMatch.  Additionally, functionality to create spectrum compatible
! motions depends on use of a system call to RspMatch for processing.
! Details on RspMatch can be found in :
! Atik and Ambrahamson (2010) "An Improved Method for Nonstationary 
! Spectral Matching." Earthquake Spectra, 26(3), 601-617


MODULE Seismic

! This MODULE stores seismic parameters.

USE                                 NWTC_Library

Integer(4)                          :: SeismicMode                                     ! Switch to identify the type of ground motion ( User-defined or synthethic motion)
Integer(4)                          :: nPass                                           ! Integer identifying the number of passes for spectral matching
Integer(4)                          :: MaxIter                                         ! Integer identifying the number of iterations for each pass for spectral matching
INTEGER(4)                          :: PtfmMotionType                                  ! Type of motion specified for the platform 1: acceleration, 2: velocity, and 3: displacement
INTEGER(4)                          :: numTimeSteps                                    ! The number of time steps in the whole simulation considering EQDelay
INTEGER(4)                          :: firstMotionIndex                                ! The index of the first acceleration
INTEGER(4)                          :: SynRandomSeedX(2)                               ! Random seeds for creating synthetic motion in X direction [-2147483648 to 2147483647]
INTEGER(4)                          :: SynRandomSeedY(2)                               ! Random seeds for creating synthetic motion in X direction [-2147483648 to 2147483647]
INTEGER(4)                          :: SynRandomSeedZ(2)                               ! Random seeds for creating synthetic motion in X direction [-2147483648 to 2147483647]
INTEGER(4)                          :: lastMotionIndex                                 ! The index of the last acceleration for use when detrending the data
LOGICAL                             :: InitSeismic = .TRUE.                            ! Switch to identify if the seismic parameters should be initialized.
LOGICAL                             :: BLineCorrection                                 ! Flag Specifying the need for Baseline Correction
LOGICAL                             :: TgtResponse                                     ! Flag Specifying the need for response spectral matching
LOGICAL                             :: PtfmArraysAllocated = .FALSE.
REAL(ReKi)                           :: SMTol                                           ! variable specifying the amount of tolerence needed for spectral matching
REAL(ReKi)                           :: MaxFreq                                         ! Indicates the maximum frequency up to which spectral matching is performed
REAL(ReKi)                           :: FreqMatch12(2)                                  ! Two frequencies which clarify the frequency range for spectral matching which are separated by commas
REAL(ReKi)                           :: SynDuration                                     ! Number specifying the duration of synthetic motion
REAL(ReKi)                           :: PtfmMotionFactor                                ! Factor applied to the platform motion files
REAL(ReKi)                           :: SynRMSAmp(3)                                    ! A value defining RMS amplitude for synthetic acceleration in the X, Y, and Z direction
REAL(ReKi)                          :: ActFreq
! Yang: 2018-5-21
REAL(ReKi)                           :: NatFreq                                         ! The first order natural Frequency of the wind turbine system (used to generate synthethic motions)
! End of the modification purpose
REAL(ReKi)                          :: EqDelay                                         ! Time to initiate seismic loading
! It has been defined in FAST_Mods.f90
!CHARACTER(128)                      :: SeismicFile                                 ! A string with the filename of the seismic file.
CHARACTER(128)                      :: PtfmXMotionFn                                   ! A string with the filename of the acceleration in the X direction
CHARACTER(128)                      :: PtfmYMotionFn                                   ! A string with the filename of the acceleration in the Y direction
CHARACTER(128)                      :: PtfmZMotionFn                                   ! A string with the filename of the acceleration in the Z direction
CHARACTER(128)                      :: SynInitRampFnX                                  ! A string with the filename of the initial ramping for synthetic acceleration in the X direction
CHARACTER(128)                      :: SynInitRampFnY                                  ! A string with the filename of the initial ramping for synthetic acceleration in the Y direction
CHARACTER(128)                      :: SynInitRampFnZ                                  ! A string with the filename of the initial ramping for synthetic acceleration in the Z direction
CHARACTER(128)                      :: SynFinalRampFnX                                 ! A string with the filename of the final ramping for synthetic acceleration in the X direction
CHARACTER(128)                      :: SynFinalRampFnY                                 ! A string with the filename of the final ramping for synthetic acceleration in the Y direction
CHARACTER(128)                      :: SynFinalRampFnZ                                 ! A string with the filename of the final ramping for synthetic acceleration in the Z direction
CHARACTER(128)                      :: TgtResFnX                                       ! A string with the filename of the Target Response Spectrum File in the X direction
CHARACTER(128)                      :: TgtResFnY                                       ! A string with the filename of the Target Response Spectrum File in the Y direction
CHARACTER(128)                      :: TgtResFnZ                                       ! A string with the filename of the Target Response Spectrum File in the Z direction

REAL(ReKi), ALLOCATABLE             :: PtfmAccel(:,:)                                  ! An array with the platform acceleration used to calculate base force
REAL(ReKi), ALLOCATABLE             :: PtfmVel(:,:)                                    ! An array with the platform velocity used to calculate base force
REAL(ReKi), ALLOCATABLE             :: PtfmDisp(:,:)                                   ! An array with the platform velocity used to calculate base force

CONTAINS

!=======================================================================
!++++++++++++++++++++++++++++++++++++++
!=======================================================================
      SUBROUTINE ReadSeismicFile()
      USE                                 Precision
      USE                                 General,ONLY:UnIn,SeismicFile
      USE                                 SimCont                                            ! Use simcont to get the length of the simulation
      USE                                 SysSubs

      IMPLICIT                            NONE

      INTEGER(4)                          :: IOS =0                                          ! I/O status returned from the read statement.
      INTEGER(4)                          :: numTimeSteps                                    ! The number of time steps in the simulation
      INTEGER(4)                          :: numTimeStepsData                                ! The number of time steps in the simulation

      INTEGER(4)                          :: Sttus                                           ! Status of allocation attempts.
      INTEGER(4)                          :: strLen1
      INTEGER(4)                          :: strLen2
      INTEGER(4)                          :: strLen3
      INTEGER(4)                          :: strLenInit
      INTEGER(4)                          :: strLenFinal
      INTEGER(4)                          :: strLenTgt 
      INTEGER(4)                          :: strLenSyn 
      INTEGER(4)                          :: result1 =0
      LOGICAL                             :: ex = .FALSE.                                    ! Used to see if files exist
      INTEGER(4)                          :: ex1 = 0                                         ! Flag for existance of x motion file
      INTEGER(4)                          :: ex2 = 0                                         ! Flag for existance of y motion file
      INTEGER(4)                          :: ex3 = 0                                         ! Flag for existance of z motion file

      INTEGER(4)                          :: motionSteps(3) = (/ 0, 0, 0 /)                  ! The maximum number of motion steps
      INTEGER(4)                          :: maxMotionSteps = 0                              ! The maximum number of motion steps
      REAL(ReKi), ALLOCATABLE             :: xMotionData(:)                                  ! The motion data in the x direction
      REAL(ReKi), ALLOCATABLE             :: yMotionData(:)                                  ! The motion data in the y direction
      REAL(ReKi), ALLOCATABLE             :: zMotionData(:)                                  ! The motion data in the z direction

      REAL(ReKi), ALLOCATABLE             :: rawMotionData(:,:)                              ! The motion data in the z direction
      REAL(ReKi), ALLOCATABLE             :: dispArray(:,:)
      REAL(ReKi), ALLOCATABLE             :: velArray(:,:)
      REAL(ReKi), ALLOCATABLE             :: accelArray(:,:)
      REAL(ReKi), ALLOCATABLE             :: correctedAccel(:,:)
      
      REAL(ReKi)                          :: motionStartTime = 0                             ! The time that the motion data starts at.  Currently ignored.
      
      INTEGER(4)                          :: tgtOutSteps(3) = (/ 0, 0, 0 /)                  ! The maximum number of motion steps in the motion 
      REAL(ReKi), ALLOCATABLE             :: xTgtAccel(:)                                    ! The acceleration data in the x direction corrected using RSPMatch
      REAL(ReKi), ALLOCATABLE             :: yTgtAccel(:)                                    ! The acceleration data in the y direction corrected using RSPMatch
      REAL(ReKi), ALLOCATABLE             :: zTgtAccel(:)                                    ! The acceleration data in the z direction corrected using RSPMatch
      REAL(ReKi), ALLOCATABLE             :: targetAccel(:,:)                                ! A matrix of all three acceleration components corrected using RSPMatch
      
      ! Input related variables
      CHARACTER(64)                       :: statusString                                    ! A string to store a status message in
      INTEGER(4)                          :: iStep                                           ! The current time step
      INTEGER(4)                          :: iDir                                            ! The current direction

!      REAL(ReKi)                          :: tempDamp
      INTEGER(4)                          :: I
      REAL(ReKi)                          :: velZeroStep(3)                                  ! The change in velocity to get to zero velocity in a second
      INTEGER(4)                          :: oneSecondSteps                                  ! The number of time steps in one second

      
      ! Open the FAST seismic input file:
      CALL OpenFInpFile(UnIn, SeismicFile )

      READ (UnIn,'(//)',IOSTAT=IOS)

      IF ( IOS < 0 )  THEN
            CALL PremEOF(SeismicFile, 'unused FAST seismic-file header' )
      ENDIF
      CALL ReadIVar(UnIn, SeismicFile, SeismicMode, 'SeismicMode', 'variable specifying the type of motion (user-defined or synthetic) for the platform' )
      IF ( ( SeismicMode <= 0 ) .OR. ( SeismicMode >= 3 ) ) THEN
            CALL ProgAbort(' SeismicMode takes values of 1 or 2')
      ENDIF
      CALL ReadCVar(UnIn, SeismicFile, PtfmXMotionFn, 'PtfmXMotionFn', 'file containing motion in the X direction' )
      CALL ReadCVar(UnIn, SeismicFile, PtfmYMotionFn, 'PtfmYMotionFn', 'file containing motion in the Y direction' )
      CALL ReadCVar(UnIn, SeismicFile, PtfmZMotionFn, 'PtfmZMotionFn', 'file containing motion in the Z direction' )
      CALL ReadIVar(UnIn, SeismicFile, PtfmMotionType, 'PtfmMotionType', 'variable specifing the type of motion for the platform' )
! Yang: Start for adding a judgement for the value of PtfmMotionType
      IF ( ( PtfmMotionType <= 0 ) .OR. ( PtfmMotionType >= 4 ) ) THEN
            CALL ProgAbort(' PtfmMotionType takes values of 1, 2 or 3')
      ENDIF
!Yang: end for adding a judgement for the value of PtfmMotionType
      CALL ReadRVar(UnIn, SeismicFile, PtfmMotionFactor, 'PtfmMotionFactor', 'The factor to multiply the input motion by to get to required units' )
! Yang: Start for reading the added parameters     
      CALL ReadRVar(UnIn, SeismicFile, NatFreq, 'NatFreq', 'The first order natural Frequency of the wind turbine system' )
! Yang: End for reading the added parameters    
      CALL ReadRVar(UnIn, SeismicFile, EqDelay, 'EqDelay', 'The initiation time for seismic loading' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedX(1), 'SynRandomSeedX1', 'Integer specifing the seed number for random number generation for X direction' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedX(2), 'SynRandomSeedX2', 'Integer specifing the seed number for random number generation for X direction' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedY(1), 'SynRandomSeedY1', 'Integer specifing the seed number for random number generation for Y direction' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedY(2), 'SynRandomSeedY2', 'Integer specifing the seed number for random number generation for Y direction' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedZ(1), 'SynRandomSeedZ1', 'Integer specifing the seed number for random number generation for Z direction' )
      CALL ReadIVar(UnIn, SeismicFile, SynRandomSeedZ(2), 'SynRandomSeedZ2', 'Integer specifing the seed number for random number generation for Z direction' )
      CALL ReadRVar(UnIn, SeismicFile, SynRMSAmp(1), 'SynRMSAmp1', 'The RMS amplitude used for synthetic motion in X direction in m/s^2' )
      CALL ReadRVar(UnIn, SeismicFile, SynRMSAmp(2), 'SynRMSAmp2', 'The RMS amplitude used for synthetic motion in Y direction in m/s^2' )
      CALL ReadRVar(UnIn, SeismicFile, SynRMSAmp(3), 'SynRMSAmp3', 'The RMS amplitude used for synthetic motion in Z direction in m/s^2' )
      CALL ReadRVar(UnIn, SeismicFile, SynDuration, 'SynDuration', 'Duration of the synthetic motion' )
      CALL ReadCVar(UnIn, SeismicFile, SynInitRampFnX, 'SynInitRampFnX', 'Initial ramping file name for synthetic motion in X direction' )
      CALL ReadCVar(UnIn, SeismicFile, SynInitRampFnY, 'SynInitRampFnY', 'Initial ramping file name for synthetic motion in Y direction' )
      CALL ReadCVar(UnIn, SeismicFile, SynInitRampFnZ, 'SynInitRampFnZ', 'Initial ramping file name for synthetic motion in Z direction' )
      CALL ReadCVar(UnIn, SeismicFile, SynFinalRampFnX, 'SynFinalRampFnX', 'Final ramping file name for synthetic motion in X direction' )
      CALL ReadCVar(UnIn, SeismicFile, SynFinalRampFnY, 'SynFinalRampFnY', 'Final ramping file name for synthetic motion in Y direction' )
      CALL ReadCVar(UnIn, SeismicFile, SynFinalRampFnZ, 'SynFinalRampFnZ', 'Final ramping file name for synthetic motion in Z direction' )
      CALL ReadLVar(UnIn, SeismicFile, BLineCorrection, 'BLineCorrection', 'variable specifying if motion should be baseline corrected' )
      CALL ReadLVar(UnIn, SeismicFile, TgtResponse, 'TgtResponse', 'variable specifying if spectral matching should be performed on motion' )
      CALL ReadCVar(UnIn, SeismicFile, TgtResFnX, 'TgtResFnX', 'Target response file name in X direction' )
      CALL ReadCVar(UnIn, SeismicFile, TgtResFnY, 'TgtResFnY', 'Target response file name in Y direction' )
      CALL ReadCVar(UnIn, SeismicFile, TgtResFnZ, 'TgtResFnZ', 'Target response file name in Z direction' )
      CALL ReadIVar(UnIn, SeismicFile, nPass, 'nPass', 'Integer specifying the number of passes in spectral matching' )
      CALL ReadIVar(UnIn, SeismicFile, MaxIter, 'MaxIter', 'Integer specifying maximum number of interations in each pass for spectral matching' )
      CALL ReadRVar(UnIn, SeismicFile, SMTol, 'SMTol', 'Tolerence used for spectral matching' )
      CALL ReadRVar(UnIn, SeismicFile, MaxFreq, 'MaxFreq', 'maximum frequency up to which spectral matching is performed' )
      READ (UnIn,*,IOSTAT=IOS)  ( FreqMatch12(I), I=1,2 )
      IF ( IOS < 0 ) THEN
            CALL PremEOF ( SeismicFile , 'FreqMatch12' )
      ELSEIF ( IOS > 0 ) THEN
            CALL WrScr1 (' Invalid numerical input for file "'//TRIM( SeismicFile )//'.')
            CALL ProgAbort(' FreqMatch12 should have two variables seperated by comma.')
      ENDIF

      ! That should be everything from the seismic file.  Go ahead and close it
      CLOSE ( UnIn )

      ! Convert the actuator frequency to radians per second.
      ActFreq= NatFreq*2.0*Pi*10.0
      ! divide actuator damping by 100
      !ActDamp= tempDamp/100.0
      ! divide Spectral Matching Tolerence by 100
      SMTol = SMTol/100.0
      ! Figure out how many time steps are in the simulation
      numTimeSteps = CEILING(TMax/DT)+1

      ! Allocate memory to store the data from the acceleration file
      ALLOCATE(PtfmAccel(3,numTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the PtfmAccel array.')
      ELSE
            PtfmArraysAllocated = .TRUE.
      ENDIF

      ALLOCATE(PtfmVel(3,numTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the PtfmVel array.')
      ELSE
            PtfmArraysAllocated = .TRUE.
      ENDIF
      ALLOCATE(PtfmDisp(3,numTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the PtfmDisp array.')
      ELSE
            PtfmArraysAllocated = .TRUE.
      ENDIF


      CALL ZeroFillArray(numTimeSteps)

      IF (SeismicMode == 1) THEN
            ! if the user has provided the input motion go ahead and read files 
            ! Check if files exist and abort program if they are printed incorect in the input file and 
            ! figure out the last direction for response match
            strLen1 = LEN_TRIM(PtfmXMotionFn)
            strLen2 = LEN_TRIM(PtfmYMotionFn)
            strLen3 = LEN_TRIM(PtfmZMotionFn)

            ! Check if a motion file was specified that the motion files exist otherwise abort
            CALL CheckFileExists( PtfmXMotionFn, 'X', ex1)
            CALL CheckFileExists( PtfmYMotionFn, 'Y', ex2)
            CALL CheckFileExists( PtfmZMotionFn, 'Z', ex3)

            IF (ex1 .EQ. 0 .AND. ex2 .EQ. 0 .AND. ex3 .EQ. 0) THEN
                  CALL ProgAbort ( ' No motion specified in any direction' )
            ENDIF
            ! if the x motion file exists then read it

            IF( ex1 ) THEN
                  CALL WrScr1(' Using '//TRIM(PtfmXMotionFn(1:strLen1))//' for X motion file')
                  CALL ReadMotionFile(PtfmXMotionFn, PtfmMotionFactor, xMotionData, motionStartTime)
                  motionSteps(1) = SIZE(xMotionData)
            ENDIF
            ! if the y motion file exists then read it

            IF( ex2 ) THEN
                  CALL WrScr1(' Using '//TRIM(PtfmYMotionFn(1:strLen2))//' for Y motion file')
                  CALL ReadMotionFile(PtfmYMotionFn, PtfmMotionFactor, yMotionData, motionStartTime)
                  motionSteps(2) = SIZE(yMotionData)
            ENDIF
            ! if the z motion file exists then read it

            IF( ex3 ) THEN
                  CALL WrScr1(' Using '//TRIM(PtfmZMotionFn(1:strLen3))//' for Z motion file')
                  CALL ReadMotionFile(PtfmZMotionFn, PtfmMotionFactor, zMotionData, motionStartTime)
                  motionSteps(3) = SIZE(zMotionData)
            ENDIF
      ELSEIF (SeismicMode == 2 .AND. PtfmMotionType .NE. 1) THEN
            CALL ProgAbort( 'PtfmMotionType must be 1 for synthetic motion.' )
      ELSEIF (SeismicMode == 2) THEN
            IF ( SynRMSAmp(1) .GT. 0.0 ) THEN
                  PtfmXMotionFn = 'SynMotionX.dat'
                  ex1 = 1
                  CALL SynMotionGen(SynRMSAmp(1), SynDuration, SynRandomSeedX, SynInitRampFnX, SynFinalRampFnX, ActFreq, DT, xMotionData)
                  motionSteps(1) = SIZE(xMotionData)
                  CALL WriteMotionFile(xMotionData, PtfmXMotionFn, DT)
            ELSE
                  ex1 = 0
            ENDIF
            IF ( SynRMSAmp(2) .GT. 0.0 ) THEN
                  PtfmYMotionFn = 'SynMotionY.dat'
                  ex2 = 1
                  CALL SynMotionGen(SynRMSAmp(2), SynDuration, SynRandomSeedY, SynInitRampFnY, SynFinalRampFnY, ActFreq, DT, yMotionData)
                  motionSteps(2) = SIZE(yMotionData)
                  CALL WriteMotionFile(yMotionData, PtfmYMotionFn, DT)
            ELSE
                  ex2 = 0
            ENDIF
            IF ( SynRMSAmp(3) .GT. 0.0 ) THEN
                  PtfmZMotionFn = 'SynMotionZ.dat'
                  ex3 = 1
                  CALL SynMotionGen(SynRMSAmp(3), SynDuration, SynRandomSeedZ, SynInitRampFnZ, SynFinalRampFnZ, ActFreq, DT, zMotionData)
                  motionSteps(3) = SIZE(zMotionData)
                  CALL WriteMotionFile(zMotionData, PtfmZMotionFn, DT)
            ELSE
                  ex3 = 0
            ENDIF
      ENDIF

      ! Figure out what direction has the largest number of points
      maxMotionSteps = MAXVAL(motionSteps)
      ALLOCATE(rawMotionData(3,maxMotionSteps), STAT=Sttus)
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the rawMotionData array.' )
      ENDIF
      
      ! Copy in the X data
      CALL getXYZMotion( rawMotionData, maxMotionSteps, xMotionData, motionSteps(1), 1)
      ! Copy in the Y data
      CALL getXYZMotion( rawMotionData, maxMotionSteps, yMotionData, motionSteps(2), 2)
      ! Copy in the Z data
      CALL getXYZMotion( rawMotionData, maxMotionSteps, zMotionData, motionSteps(3), 3)
      ! Calculate displacement, velocity, and acceleration from the given motion
      CALL getDVA( rawMotionData, maxMotionSteps, PtfmMotionType, dispArray, velArray, accelArray)
      
      
      ! If requested, baseline correct the motion.
      IF(BLineCorrection) THEN
            ! Allocate memory for the corrected acceleration
            ALLOCATE(correctedAccel(3,maxMotionSteps), STAT=Sttus)
            IF ( Sttus /= 0 ) THEN
                  CALL ProgAbort( ' Error allocating memory for the correctedAccel array.' )
            ENDIF
            ! Do baseline correction on each of the three directions
            DO iDir = 1,3
                  CALL BaselineCorrect(maxMotionSteps, iDir, accelArray, velArray, dispArray, correctedAccel)
            ENDDO
            ! Cleanup the disp, vel, and accel arrays without baseline correction
            DEALLOCATE(dispArray)
            DEALLOCATE(velArray)
            DEALLOCATE(accelArray)
            ! Calculate the disp, vel, and accel arrays with baseline correction
            CALL getDVA( correctedAccel, maxMotionSteps, 1, dispArray, velArray, accelArray)
            ! Cleanup the corrected acceleration because we are done using it
            DEALLOCATE(correctedAccel)
      ENDIF
      
      IF(TgtResponse) THEN
            DO iDir = 1,3
                  IF (iDir .EQ. 1 .AND. ex1 .EQ. 1) THEN
                        CALL TargetResonse(maxMotionSteps, iDir, accelArray, TgtResFnX, PtfmXMotionFn, xTgtAccel, tgtOutSteps(iDir))
                  ELSEIF (iDir .EQ. 2 .AND. ex2 .EQ. 1) THEN
                        CALL TargetResonse(maxMotionSteps, iDir, accelArray, TgtResFnY, PtfmYMotionFn, yTgtAccel, tgtOutSteps(iDir))
                  ELSEIF (iDir .EQ. 3 .AND. ex3 .EQ. 1) THEN
                        CALL TargetResonse(maxMotionSteps, iDir, accelArray, TgtResFnZ, PtfmZMotionFn, zTgtAccel, tgtOutSteps(iDir))
                  ENDIF
            ENDDO
            
            maxMotionSteps = MAXVAL(tgtOutSteps)
            
            ALLOCATE(targetAccel(3,maxMotionSteps), STAT=Sttus)
            IF ( Sttus /= 0 ) THEN
                  CALL ProgAbort( ' Error allocating memory for the targetAccel array.' )
            ENDIF
            
            DO iStep = 1,maxMotionSteps
                  IF (iStep .LT. tgtOutSteps(1)) THEN
                        targetAccel(1,iStep) = xTgtAccel(iStep)
                  ELSE
                        targetAccel(1,iStep) = 0
                  ENDIF
            ENDDO
            DO iStep = 1,maxMotionSteps
                  IF (iStep .LT. tgtOutSteps(2)) THEN
                        targetAccel(2,iStep) = yTgtAccel(iStep)
                  ELSE
                        targetAccel(2,iStep) = 0
                  ENDIF
            ENDDO
            DO iStep = 1,maxMotionSteps
                  IF (iStep .LT. tgtOutSteps(3)) THEN
                        targetAccel(3,iStep) = zTgtAccel(iStep)
                  ELSE
                        targetAccel(3,iStep) = 0
                  ENDIF
            ENDDO
            ! Cleanup the disp, vel, and accel arrays without target spectrum correction
            DEALLOCATE(dispArray)
            DEALLOCATE(velArray)
            DEALLOCATE(accelArray)
            IF (ex1 .EQ. 1) THEN
                  DEALLOCATE(xTgtAccel)
            ENDIF
            IF (ex2 .EQ. 1) THEN
                  DEALLOCATE(yTgtAccel)
            ENDIF
            IF (ex3 .EQ. 1) THEN
                  DEALLOCATE(zTgtAccel)
            ENDIF
            ! Calculate the disp, vel, and accel arrays with target spectrum correction
            CALL getDVA( targetAccel, maxMotionSteps, 1, dispArray, velArray, accelArray)
            DEALLOCATE(targetAccel)
      ENDIF
      
      
      firstMotionIndex = EqDelay/DT
      
      ! Figure out how much to remove from the velocity and displacement
      ! for each time step at the end of ground motion
      IF (firstMotionIndex+maxMotionSteps+oneSecondSteps .GT. numTimeSteps) THEN
            CALL ProgAbort ( 'The simulation is not long enough to include all seismic motion.' )
      ENDIF
      ! Copy the resulting displacement, velocity, and acceleration array into the appropriate variables
      DO iStep = 1,maxMotionSteps
            IF (iStep<= maxMotionSteps) THEN
                  DO iDir = 1,3
                        PtfmAccel(iDir,iStep+firstMotionIndex) = accelArray(iDir,iStep)
                        PtfmVel(iDir,iStep+firstMotionIndex)   = velArray(iDir,iStep)
                        PtfmDisp(iDir,iStep+firstMotionIndex)  = dispArray(iDir,iStep)
                  ENDDO
            ENDIF
      ENDDO
      ! Bring the velocity to zero over one second
      
      oneSecondSteps = 1/DT 
      DO iDir = 1,3
            velZeroStep(iDir) = velArray(iDir,maxMotionSteps)/oneSecondSteps
      ENDDO
      
      ! Set sensible values after the specified motion stops
      DO iStep = 1,oneSecondSteps
            DO iDir = 1,3
                  PtfmAccel(iDir,firstMotionIndex+maxMotionSteps+iStep) = 0
                  ! Ramp the velocity down
                  PtfmVel(iDir,firstMotionIndex+maxMotionSteps+iStep)   = velArray(iDir,maxMotionSteps)-iStep*velZeroStep(iDir)
                  ! Calcuate the resulting displacement based on the ramped velocity
                  PtfmDisp(iDir,firstMotionIndex+maxMotionSteps+iStep)  = PtfmDisp(iDir,firstMotionIndex+maxMotionSteps+iStep-1)+PtfmVel(iDir,firstMotionIndex+maxMotionSteps+iStep-1)*DT
            ENDDO
      ENDDO
      
      ! Set sensible values after the specified motion stops
      DO iStep = firstMotionIndex+maxMotionSteps+oneSecondSteps+1,numTimeSteps
            DO iDir = 1,3
                  PtfmAccel(iDir,iStep) = 0
                  ! Keep the final diplacement
                  PtfmDisp(iDir,iStep)  = PtfmDisp(iDir,iStep-1)
                  ! Set the velocity to zero
                  PtfmVel(iDir,iStep)   = 0
            ENDDO
      ENDDO
      
      ! Clean up memory and return.
      DEALLOCATE(rawMotionData)
      DEALLOCATE(dispArray)
      DEALLOCATE(velArray)
      DEALLOCATE(accelArray)
      END SUBROUTINE ReadSeismicFile
      
!=======================================================================
      SUBROUTINE getXYZMotion( motionArray, maxSteps, curData, maxDirSteps, dirIndex)

      IMPLICIT                            NONE
      REAL(ReKi), INTENT(INOUT)           :: motionArray(:,:)                                ! The number of time steps in the simulation
      INTEGER(4), INTENT(IN)              :: maxSteps
      REAL(ReKi), INTENT(IN)              :: curData(:)                                      ! The number of time steps in the simulation
      INTEGER(4), INTENT(IN)              :: maxDirSteps
      INTEGER(4), INTENT(IN)              :: dirIndex
      
      
      INTEGER(4)                          :: iTimeStep                                       ! The current time step
            

      DO iTimeStep=1,maxSteps
            IF (iTimeStep <= maxDirSteps) THEN
                  motionArray(dirIndex,iTimeStep) = curData(iTimeStep)
            ELSE
                  motionArray(dirIndex,iTimeStep) = 0
            ENDIF
      ENDDO

      END SUBROUTINE getXYZMotion
      
      !=======================================================================
      SUBROUTINE getDVA( motionArray, steps, mType, dArray, vArray, aArray)

      USE                                 SimCont,ONLY:DT
      IMPLICIT                            NONE
      
      REAL(ReKi), INTENT(IN)              :: motionArray(:,:)
      INTEGER(4), INTENT(IN)              :: steps
      INTEGER(4), INTENT(IN)              :: mType
      REAL(ReKi), ALLOCATABLE, INTENT(OUT) :: dArray(:,:)                                     ! The number of time steps in the simulation
      REAL(ReKi), ALLOCATABLE, INTENT(OUT) :: vArray(:,:)
      REAL(ReKi), ALLOCATABLE, INTENT(OUT) :: aArray(:,:)
      
      INTEGER(4)                          :: iStep
      INTEGER(4)                          :: iDir
      INTEGER(4)                          :: Sttus
      
      ALLOCATE(dArray(3,steps), STAT=Sttus)
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the dArray array.' )
      ENDIF
      ALLOCATE(vArray(3,steps), STAT=Sttus)
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the vArray array.' )
      ENDIF
      ALLOCATE(aArray(3,steps), STAT=Sttus)
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the aArray array.' )
      ENDIF
      DO iStep=1,steps
            DO iDir=1,3
                  IF (mType == 1) THEN
                        aArray(iDir,iStep) = motionArray(iDir,iStep)
                        IF ( iStep == 1) THEN
                              dArray(iDir,iStep) = 0
                              vArray(iDir,iStep) = 0
                        ELSE
                              vArray(iDir,iStep) = vArray(iDir,iStep-1) + DT/2*(aArray(iDir,iStep)+aArray(iDir,iStep-1))
                              dArray(iDir,iStep) = dArray(iDir,iStep-1) + DT/2*(vArray(iDir,iStep)+vArray(iDir,iStep-1))
                        ENDIF
                  ELSEIF (mType == 2) THEN
                        vArray(iDir,iStep) = motionArray(iDir,iStep)
                        IF ( iStep == 1) THEN
                              dArray(iDir,iStep) = 0
                              aArray(iDir,iStep) = 0
                        ELSE
                              aArray(iDir,iStep) = (vArray(iDir,iStep-1)-vArray(iDir,iStep))/DT
                              dArray(iDir,iStep) = dArray(iDir,iStep-1) + DT/2*(vArray(iDir,iStep)+vArray(iDir,iStep-1))
                        ENDIF
                  ELSEIF (mType == 3) THEN
                        dArray(iDir,iStep) = motionArray(iDir,iStep)
                        IF ( iStep == 1) THEN
                              vArray(iDir,iStep) = 0
                              aArray(iDir,iStep) = 0
                        ELSE
                              vArray(iDir,iStep) = (dArray(iDir,iStep-1)-dArray(iDir,iStep))/DT
                              aArray(iDir,iStep) = (vArray(iDir,iStep-1)-vArray(iDir,iStep))/DT
                        ENDIF
                  ELSE
                        CALL ProgAbort('PtfmMotionType must be 1, 2, or 3')
                  ENDIF
            ENDDO
      ENDDO
      
      END SUBROUTINE getDVA
!=======================================================================
      SUBROUTINE ReadMotionFile(fileName, fileFactor, motionData, motionStartTime)
      
      USE                                 SimCont                                            ! Use simcont to get the length of the simulation

      IMPLICIT                            NONE
      CHARACTER(128), INTENT(IN )         :: fileName
      REAL(ReKi), INTENT(IN )             :: fileFactor                                      ! Factor for scaling values in the data file
      REAL(ReKi), INTENT(OUT), ALLOCATABLE :: motionData(:)                                  ! Interpolated motion data
      REAL(ReKi), INTENT(OUT)             :: motionStartTime                                 ! The time that the motion starts
      
      REAL(ReKi), ALLOCATABLE             :: fileTimeData(:)
      REAL(ReKi), ALLOCATABLE             :: fileMotionData(:)
      INTEGER(4)                          :: numTimeStepsFile = 0                            ! The number of time steps in the data file
      INTEGER(4)                          :: numTimeStepsInt                                 ! The number of time steps after interpolation
      INTEGER(4)                          :: lastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call
      REAL(ReKi)                          :: motionTime                                      ! Duration of motion
      INTEGER(4)                          :: iTimeStep                                       ! The current time step
      INTEGER                             :: lineNumberData = 0
      REAL(ReKi)                          :: curFileValueData                                ! The motion value on the current line for raw motion file.
      REAL(ReKi)                          :: curTimeValueData                                ! The time value on the current line for raw motion file.
      INTEGER(4)                          :: finalMotionTime                                 ! The motion value on the current line.
      INTEGER(4)                          :: IOS = 0                                         ! I/O status returned from the read statement.
      LOGICAL                             :: EOF = .FALSE.                                   ! Status for reading raw motion file
      REAL(ReKi)                          :: curMotionValue                                  ! The value on the current line.
      REAL(ReKi)                          :: curTimeValue                                    ! The time value on the current line.
      INTEGER(4)                          :: Sttus
      INTEGER(4)                          :: UnInMtn = 34                                    ! File pointer for use in files to be written in TargetResonse subroutine
      CHARACTER(128)                      :: StatusString

      numTimeStepsFile = 0
      lineNumberData = 0
      IOS = 0   
      EOF = .FALSE.

      !-----------------------------------------------------------------------
      ! Figure out number of time steps in the data file  
      !-----------------------------------------------------------------------
      CALL OpenFInpFile(UnInMtn, fileName, IOS)
      
      READ(UnInMtn, *, IOSTAT=IOS) curTimeValueData, curFileValueData
      
      DO WHILE (IOS == 0)
            numTimeStepsFile = numTimeStepsFile + 1
            READ(UnInMtn, *, IOSTAT=IOS) curTimeValueData, curFileValueData
      ENDDO
      
      !-------------------------------------------------------------------------------------------------
      ! Rewind the file (to the beginning)
      !-------------------------------------------------------------------------------------------------
      REWIND( UnInMtn )
      !-----------------------------------------------------------------------
      ! Allocate Arrays for time and motion
      !-----------------------------------------------------------------------
      ALLOCATE(fileTimeData(numTimeStepsFile), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort (' Error allocating memory for the fileTimeData array.')
      ENDIF
      ALLOCATE(fileMotionData(numTimeStepsFile), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the fileMotionData array.' )
      ENDIF 

      !-----------------------------------------------------------------------
      !Reading input motion file 
      !-----------------------------------------------------------------------
      DO iTimeStep=1,numTimeStepsFile
            ! lineNumber will be set to -1 if we encounter the end of file
            IF (lineNumberData >= 0) THEN
                  READ(UnInMtn, *, IOSTAT=IOS) curTimeValueData,curFileValueData
                  IF (IOS == 0 .AND. iTimeStep <= numTimeStepsFile) THEN
                        lineNumberData = lineNumberData + 1
                        fileTimeData(iTimeStep) = curTimeValueData
                        fileMotionData(iTimeStep) = fileFactor*curFileValueData
                        ! Save the current motion for the current time step
                  ELSE
                        ! We ran out of input.  Hold the last motion value we had for all
                        ! the following time steps.
                        lineNumberData = -1
                  ENDIF
            ENDIF
            IF (iTimeStep == numTimeStepsFile) THEN
                  lineNumberData = -1
            ENDIF
      ENDDO

      ! Give a warning if we did not get to the end of the acceleration file     
      IF (lineNumberData >=0) THEN
            CALL WrScr1(' Warning : It appears like not all data was read from '//fileName)
            CALL WrScr1('')
      ENDIF
                       
      !-----------------------------------------------------------------------
      !Interpolat motion according to timestep in main file 
      !-----------------------------------------------------------------------
      motionStartTime = CEILING(fileTimeData(1)/DT)*DT
      motionTime      = fileTimeData(numTimeStepsFile)
      numTimeStepsInt = CEILING(motionTime/DT)+1

      ALLOCATE(motionData(numTimeStepsInt), STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
            CALL ProgAbort( 'Error allocating memory for the motionData array.' )
      ENDIF

      DO iTimeStep=1,numTimeStepsInt
            IF (iTimeStep == 1) THEN
                  motionData(iTimeStep) = fileMotionData(1)
            ELSE
                  motionData(iTimeStep) = InterpStpReal(motionStartTime+(iTimeStep-1)*DT, fileTimeData, fileMotionData, lastInd, numTimeStepsFile )
            ENDIF
      ENDDO

      CALL WrScr1(' Done reading motion file: '//fileName)
      
      ! Close the input file
      CLOSE ( UnInMtn )

      ! Deallocate arrays for next motion   

      DEALLOCATE(fileTimeData)
      DEALLOCATE(fileMotionData)
      
      END SUBROUTINE ReadMotionFile

      !=======================================================
      ! subroutine ZeroFillArray
      !=======================================================

      SUBROUTINE ZeroFillArray( numTimeSteps )

      IMPLICIT                            NONE
      INTEGER(4), INTENT(IN )             :: numTimeSteps                                    ! The number of time steps in the simulation
      INTEGER(4)                          :: iDir
      INTEGER(4)                          :: iTimeStep                                       ! The current time step
      !CHARACTER(64)                      :: StatusString                                    ! A string to store a status message in

      DO iDir=1,3 
            DO iTimeStep=1,numTimeSteps
                  PtfmAccel(iDir,iTimeStep) = 0
                  PtfmVel(iDir,iTimeStep)   = 0
                  PtfmDisp(iDir,iTimeStep)  = 0
            ENDDO
      ENDDO

      END SUBROUTINE ZeroFillArray

      !=======================================================
      ! subroutine CheckFileExists
      !=======================================================

      SUBROUTINE CheckFileExists( fileName, direction, result1)

      IMPLICIT                            NONE
      CHARACTER(128), INTENT(IN )         :: fileName
      CHARACTER(1), INTENT(IN )           :: direction
      INTEGER(4), INTENT(OUT )            :: result1

      LOGICAL                             :: ex = .FALSE.
      INTEGER(4)                          :: strLen

      strLen = LEN_TRIM( fileName )
      INQUIRE (file=fileName, exist=ex)
      IF (ex .EQ. .FALSE.) THEN
            IF  (strLen /= 0) THEN 
                  CALL ProgAbort( 'File defined in '//TRIM(direction)//' direction was not read properly.' )
            ELSE
                  CALL WrScr1( ' No motion is defined in '//TRIM(direction)//' direction.')   
                  result1 = 0
            ENDIF
      ELSE
            result1 = 1
      ENDIF      

      END SUBROUTINE CheckFileExists

      !=======================================================================
      !=======================================================================
      !Subroutine for baseline correction 
      !=======================================================================
      SUBROUTINE BaselineCorrect(npts, dirIndex, accel, vel, disp, correctedAccel)
      
      !USE                                 NWTC_Library
      USE                                 SimCont,ONLY:DT
      
      IMPLICIT NONE
      
      INTEGER(4), INTENT(IN)              :: npts                                            ! number of timesteps
      INTEGER(4), INTENT(IN)              :: dirIndex                                        ! The current direction (1, 2, or 3)
      REAL(ReKi), INTENT(IN)              :: accel(:,:)                                      ! Acceleration of the motion
      REAL(ReKi), INTENT(IN)              :: vel(:,:)                                        ! Velocity of the motion
      REAL(ReKi), INTENT(IN)              :: disp(:,:)                                       ! Displacement of motion
      REAL(ReKi), INTENT(INOUT)           :: correctedAccel(:,:)                             ! The corrected acceleration signal
      
      REAL(ReKi)                          :: te = 40                                         ! Precentage of timesteps at the end for applying tapered cosine function
      
      INTEGER                             :: nParam = 2                                      ! order of polynomial used
      
      REAL(ReKi), ALLOCATABLE             :: coef (:)                                        ! Coefficients used for baseline correction according to order of polynomial
      REAL(ReKi)                          :: t , tt , b1 , bPrime , bPrime2                  ! Local variables 
      REAL(ReKi)                          :: w , wPrime ,wPrime2                             ! Local variables
      REAL(ReKi)                          :: sum0 , sum1 , sum2                              ! Local variables
      INTEGER                             :: i , ii, j , n1                                  ! Local variables
      REAL(ReKi), ALLOCATABLE             :: b(:,:)
      REAL(ReKi), ALLOCATABLE             :: A(:,:)
      REAL(ReKi), ALLOCATABLE             :: AT(:,:)
      REAL(ReKi), ALLOCATABLE             :: ATA(:,:)
      REAL(ReKi), ALLOCATABLE             :: ATAI(:,:)
      REAL(ReKi), ALLOCATABLE             :: ATAIAT(:,:)
      REAL(ReKi), ALLOCATABLE             :: work (:)
      REAL(ReKi)                          :: EPS = 1.0e-11
      REAL(ReKi)                          :: DETER
      REAL(ReKi), ALLOCATABLE             :: xhat(:,:)
      INTEGER                             :: iStep
      REAL(ReKi)                          :: Y (200)
      INTEGER                             :: irow (200)
      INTEGER                             :: jcol (200)
      INTEGER                             :: jord (200)
      REAL(ReKi)                          :: MAXV, K , KM1, PIVOT , ISCAN,JSCAN , IROWK , JCOLK , AIJCK , IROWI, JCOLI 
      REAL(ReKi)                          :: INTCH , NM1 , IP1 , JTEMP, IROWJ, JCOLJ
      REAL(ReKi)                          :: INDIC = -1
      !=======================================================
      ! Allocate Arrays
      !=======================================================
      ALLOCATE (work(npts))
      ALLOCATE (coef (nParam))
      ALLOCATE (b(npts,1))
      ALLOCATE (A(npts,nParam))
      ALLOCATE (xhat(nParam,1))
      ALLOCATE (ATAIAT(nParam,npts))
      ALLOCATE (AT(nParam,npts))
      ALLOCATE (ATA(nParam,nParam))
      ALLOCATE (ATAI(nParam,nParam))
      
      
      !=======================================================
      ! Compute Baseline for uncorrected displacement
      ! Note: This code is derived from the baseline correction routine
      ! in RspMatch.  Details on RspMatch can be found in :
      ! Atik and Ambrahamson (2010) "An Improved Method for Nonstationary 
      ! Spectral Matching." Earthquake Spectra, 26(3), 601-617
      !=======================================================
      ! Load ATA and AT matrices (Ax=b) for baseline correction
      ! Load A first
      t = 0
      iStep = 10
      DO i=1,npts
            DO j=1,nparam
                  A(i,j) = 0
            ENDDO
      ENDDO
      DO i=1,npts,iStep
            t = (i-1)*DT
            tt = t*t
            DO j=1,nParam
                  A(i,j) = tt
                  tt = tt * t
            ENDDO
      ENDDO
      !!!!!!!!!!!!!!!!! Compute ATAIAT
      !! Transpose of A matrix
      DO i=1,npts
            DO j=1,nParam
                  AT(j,i) = A(i,j)
            ENDDO
      ENDDO
      !! Multiply AT by A to find ATA
      DO i = 1 , nParam
            DO j = 1, nParam
                  ATA(i,j) = 0
                  DO ii = 1, npts
                        ATA(i,j) = ATA(i,j) + AT(i,ii)*A(ii,j)
                  ENDDO
            ENDDO
      ENDDO
      !! FIND ATAI
      DO i=1,nParam
            DO j=1,nParam
                  ATAI(i,j) = ATA(i,j)
            ENDDO
      ENDDO
      !=======================================================
      ! Compute Inverse of ATA (Inverse of Singular Matrix)
      !=======================================================
      !        WHEN INDIC IS NEGATIVE, SIMUL COMPUTES THE INVERSE OF THE N BY
      !        N MATRIX A IN PLACE.  WHEN INDIC IS ZERO, SIMUL COMPUTES THE
      !        N SOLUTIONS X(1)...X(N) CORRESPONDING TO THE SET OF LINEAR
      !        EQUATIONS WITH AUGMENTED MATRIX OF COEFFICIENTS IN THE N BY
      !        N+1 ARRAY A AND IN ADDITION COMPUTES THE INVERSE OF THE
      !        COEFFICIENT MATRIX IN PLACE AS ABOVE.  IF INDIC IS POSITIVE,
      !        THE SET OF LINEAR EQUATIONS IS SOLVED BUT THE INVERSE IS NOT
      !        COMPUTED IN PLACE. THE GAUSS-JORDAN COMPLETE ELIMINATION METHOD
      !        IS EMPLOYED WITH THE MAXIMUM PIVOT STRATEGY.  ROW AND COLUMN
      !        SUBSCRIPTS OF SUCCESSIVE PIVOT ELEMENTS ARE SAVED IN ORDER IN
      !        THE IROW AND JCOL ARRAYS RESPECTIVELY.  K IS THE PIVOT COUNTER,
      !        PIVOT THE ALGEBRAIC VALUE OF THE PIVOT ELEMENT, MAX
      !        THE NUMBER OF COLUMNS IN A AND DETER THE DETERMINANT OF THE
      !        COEFFICIENTS MATRIX.  THE SOLUTIONS ARE COMPUTED IN THE (N+1)TH
      !        COLUMN OF A AND THEN UNSCRAMBLED AND PUT IN PROPER ORDER IN
      !        X(1)...X(N) USING THE PIVOT SUBSCRIPT INFORMATION AVAILABLE
      !        IN THE IROW AND JCOL ARRAYS.  THE SIGN OF THE DETERMINANT IS
      !        ADJUSTED, IF NECESSARY, BY DETERMINING IF AN ODD OR EVEN NUMBER
      !        OF PAIRWISE INTERCHANGES IS REQUIRED TO PUT THE ELEMENTS OF THE
      !        JORD ARRAY IN ASCENDING SEQUENCE WHERE JORD(IROW(I)) = JCOL(I).
      !        IF THE INVERSE IS REQUIRED, IT IS UNSCRAMBLED IN PLACE USING
      !        Y(1)...Y(N) AS TEMPORARY STORAGE.  THE VALUE OF THE DETERMINANT
      !        IS RETURNED AS THE VALUE OF THE FUNCTION.  SHOULD THE POTENTIAL
      !        PIVOT OF LARGEST MAGNITUDE BE SMALLER IN MAGNITUDE THAN EPS,
      !        THE MATRIX IS CONSIDERED TO BE SINGULAR AND A TRUE ZERO IS
      !        RETURNED AS THE VALUE OF THE FUNCTION.

            MAXV = nParam
            IF ( INDIC.GE.0 )  MAXV = nParam + 1

      !     ..... BEGIN ELIMINATION PROCEDURE .....
       5    DETER = 1.
            DO 18 K = 1, nParam
            KM1 = K - 1
      !     ..... SEARCH FOR THE PIVOT ELEMENT .....
            PIVOT = 0.
            DO 11 I = 1, nParam
            DO 11 J = 1, nParam
      !     ..... SCAN IROW AND JCOL ARRAYS FOR INVALID PIVOT SUBSCRIPTS .....
            IF ( K.EQ.1 ) GO TO 9
            DO 8 ISCAN = 1, KM1
            DO 8 JSCAN = 1, KM1
            IF ( I.EQ.IROW(ISCAN) ) GO TO 11
            IF ( J.EQ.JCOL(JSCAN) ) GO TO 11
       8    CONTINUE
       9    CONTINUE
            PIVOT = ATA(I,J)
            IROW(K) = I
            JCOL(K) = J
       11   CONTINUE

      !     ..... INSURE THAT SELECTED PIVOT IS LARGER THAN EPS .....
            IF ( ABS(PIVOT).GT.EPS ) GO TO 13
            DETER = 0.

      !     ..... UPDATE THE DETERMINANT VALUE .....
       13   IROWK = IROW(K)
            JCOLK = JCOL(K)
            DETER = DETER*PIVOT

      !     ..... NORMALIZE PIVOT ROW ELEMENTS .....
            DO 14 J = 1, MAXV
       14   ATA(IROWK,J) = ATA(IROWK,J)/PIVOT

      !     ..... CARRY OUT ELIMINATION AND  DEVELOP INVERSE .....
            ATA(IROWK,JCOLK) = 1./PIVOT
            DO 18 I = 1, nParam
            AIJCK = ATA(I,JCOLK)
            IF ( I.EQ.IROWK ) GO TO 18
            ATA(I,JCOLK) = - AIJCK/PIVOT
            DO 17 J = 1, MAXV
       17   IF ( J.NE.JCOLK ) ATA(I,J) = ATA(I,J) - AIJCK*ATA(IROWK,J)
       18   CONTINUE

      !     ..... ORDER SOLUTION VALUES (IF ANY) AND CREATE JORD ARRAY .....
            DO 20 I = 1, nParam
            IROWI = IROW(I)
            JCOLI = JCOL(I)
            JORD(IROWI) = JCOLI
       20   IF ( INDIC.GE.0 ) work(JCOLI) = ATA(IROWI,MAXV)

      !     ..... ADJUST SIGN OF DETERMINANT .....
            INTCH = 0
            NM1 = nParam - 1
            DO 22 I = 1, NM1
            IP1 = I + 1
            DO 22 J = IP1,nParam
            IF ( JORD(J).GE.JORD(I) ) GO TO 22
            JTEMP = JORD(J)
            JORD(J) = JORD(I)
            JORD(I) = JTEMP
            INTCH = INTCH + 1
       22   CONTINUE
            IF( INTCH/2*2.NE.INTCH ) DETER = - DETER

      !     ..... IF INDIC IS POSITIVE RETURN WITH RESULTS .....
            IF ( INDIC.LE.0 ) GO TO 26


      !     ..... IF INDIC IS NEGATIVE OR ZERO, UNSCRAMBLE THE INVERSE
      !           FIRST BY ROWS .....
       26   DO 28 J = 1, nParam
            DO 27 I = 1, nParam
            IROWI = IROW(I)
            JCOLI = JCOL(I)
       27   Y(JCOLI) = ATA(IROWI,J)
            DO 28 I = 1, nParam
       28   ATA(I,J) = Y(I)
      !     ..... THEN BY COLUMNS .....
            DO 30 I = 1, nParam
            DO 29 J = 1, nParam
            IROWJ = IROW(J)
            JCOLJ = JCOL(J)
       29   Y(IROWJ) = ATA(I,JCOLJ)
            DO 30 J = 1, nParam
       30   ATA(I,J) = Y(J)
      !=======================================================
      ! Form ATAIAT (multiply ATA by AT)
      !======================================================= 
      DO i = 1, nParam
            DO j = 1 , npts
                  ATAIAT(i,j) = 0
                  DO ii = 1, nParam
                        ATAIAT(i,j) = ATAIAT(i,j) + ATA(i,ii)*AT(ii,j)
                  ENDDO
            ENDDO
      ENDDO
      !=======================================================
      ! Load vector b
      !=======================================================
      DO i=1,nPts
            b(i,1) = disp(dirIndex,i)
      ENDDO
      !=======================================================
      ! Compute Xhat (multiply ATAIAT by b vector)
      !=======================================================                
      DO i = 1, nParam
            j = 1
            xhat(i,j) = 0
            DO ii = 1, npts
                  xhat(i,j) = xhat(i,j) + ATAIAT(i,ii)*b(ii,j)
            ENDDO
      ENDDO
      !=======================================================
      ! Coef parameters
      !=======================================================  
      DO i=1,nParam
            coef(i) = xHat(i,1)
      ENDDO
      !=======================================================
      ! Remove baseline and Apply Taper to end of displacement 
      ! trace (Forces D=0, V=0 at end of trace) 
      !======================================================= 
      t = 0
      n1 = npts - (npts * te) / 100
      DO i=1,npts
      
            !       Compute baseline
            sum1 = 0
            sum2 = 0
            sum0 = 0
            tt = 1
            DO j=1,nParam
                  sum0 = sum0 + coef(j)*tt*(t*t)
                  sum1 = sum1 + coef(j)*dble( j+1 ) * tt * t 
                  sum2 = sum2 + coef(j)*dble( (j+1)*j ) * tt
                  tt = tt * t
            ENDDO
            b1 = sum0
            bPrime = sum1
            bPrime2 = sum2
            !       Set taper (cosine bell)        
            IF ( i .le. n1 ) THEN
                  w = 1.
                  wPrime = 0.
                  wPrime2 = 0.
            ELSE
                  w = 0.5 * ( cos( pi*(i-n1)/(npts-n1) ) + 1 )
                  wPrime=-0.5 * pi/((npts-n1)*dt) * sin( pi*(i-n1)/(npts-n1) )
                  wPrime2 = -0.5 * (pi/((npts-n1)*dt))**2 *cos( pi*(i-n1)/(npts-n1) )
            ENDIF
            !       Compute corrected accelerogram
            correctedAccel(dirIndex,i) = (accel(dirIndex,i) - bPrime2)*w+ 2.*(vel(dirIndex,i) - bPrime)*wPrime + (disp(dirIndex,i) - b1)*wPrime2
            t = t + DT
      ENDDO

      DEALLOCATE (work)
      DEALLOCATE (coef)
      DEALLOCATE (b)
      DEALLOCATE (A)
      DEALLOCATE (xhat)
      DEALLOCATE (ATAIAT)
      DEALLOCATE (AT)
      DEALLOCATE (ATA)
      DEALLOCATE (ATAI)
      
      END SUBROUTINE BaselineCorrect
      
      !=======================================================
      ! subroutine RESPONSE SPECTRAL Matching Using RSPM09
      !=======================================================
      SUBROUTINE TargetResonse( npts, dirIndex, accel, tgtResFileName, motionFileName, matchedAccel, outStep)
      USE                                 SimCont,ONLY:DT
      USE                                 EnvCond,ONLY:Gravity
      
      IMPLICIT NONE
      
      INTEGER(4), INTENT(IN)              :: npts                                            ! number of timesteps
      INTEGER(4), INTENT(IN)              :: dirIndex                                        ! The current direction (1, 2, or 3)
      REAL(ReKi), INTENT(IN)              :: accel(:,:)                                      ! Acceleration of the motion
      CHARACTER(128), INTENT(IN)          :: tgtResFileName                                  ! A string with the filename of the Target Response Spectrum File in the X direction
      CHARACTER(128), INTENT(IN)          :: motionFileName                                  ! A string with the filename of the Target Response Spectrum File in the X direction
      REAL(ReKi), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: matchedAccel                     ! The corrected acceleration signal
      INTEGER(4), INTENT(OUT)             :: outStep                                         ! Time steps in the RSPMatch output file

      INTEGER(4)                          :: strLenTgt                                       ! String length of target response spectrum file path
      INTEGER(4)                          :: strLenMotion                                    ! String length of motion file path
      INTEGER(4)                          :: tgtLen                                          ! String length of response spectrum target file path
      INTEGER(4)                          :: motionLen                                       ! String length of motion file path
      INTEGER(4)                          :: accLen                                          ! String length of motion file path
      INTEGER(4)                          :: unmStrLen                                       ! String length of motion file path
      INTEGER(4)                          :: curLine                                         ! The current line in the output file
      INTEGER(4)                          :: curPos                                          ! The current position on the current line of the output file

      LOGICAL                             :: fileFound                                       ! Logical to test is files exist before trying to delete them
      CHARACTER(128)                      :: cmdString                                       ! String for storing system commands to be executed
      INTEGER(4)                          :: iTimeStep                                       ! Pointer for looping on time step
      INTEGER(4)                          :: UnInTGT = 31                                    ! File pointer for use in files to be written in TargetResonse subroutine
      INTEGER(4)                          :: I , J                                           ! Generic Variables 
      
      INTEGER(4)                          :: nFreq                                           ! Number of lines in the final acceleration file output (data is given in multiple columns which have to be turned into a vector)
      INTEGER(4)                          :: nDamp                                           ! Number of lines in the final acceleration file output (data is given in multiple columns which have to be turned into a vector)
      
      INTEGER(4)                          :: numLines                                        ! The number of lines in the RSPMatch output file
      INTEGER(4)                          :: nAdd                                            ! The number of points added by RSPMatch
      INTEGER(4)                          :: IOS =0                                          ! I/O status returned from the read statement.
      REAL(ReKi), ALLOCATABLE, DIMENSION(:) :: dampValues                                    ! Damping Values in target response spectra file
      REAL(ReKi), ALLOCATABLE, DIMENSION(:) :: freqValues                                    ! frequency Values in target response spectra file
      REAL(ReKi), ALLOCATABLE, DIMENSION(:) :: time1                                         ! time1 Values in target response spectra file
      REAL(ReKi), ALLOCATABLE, DIMENSION(:) :: time2                                         ! time2 Values in target response spectra file
      REAL(ReKi), ALLOCATABLE, DIMENSION(:,:) :: PSA                                         ! PSA values in target response spectra file
      REAL(ReKi), ALLOCATABLE, DIMENSION(:,:) :: PSAg                                        ! PSA Values in target response spectra file in g
      REAL(ReKi), ALLOCATABLE, DIMENSION(:,:) :: tgtMotionMat                                ! The data matrix of values in the RSPMatch output file
      REAL(ReKi)                          :: run1FreqMatch                                   ! Lower frequency for spectral matching in run 1 of RSPMatch
      REAL(ReKi)                          :: run2FreqMatch                                   ! Lower frequency for spectral matching in run 2 of RSPMatch
      REAL(ReKi)                          :: run3FreqMatch                                   ! Lower frequency for spectral matching in run 3 of RSPMatch
      REAL(ReKi)                          :: run4FreqMatch                                   ! Lower frequency for spectral matching in run 4 of RSPMatch
      REAL(ReKi)                          :: DT1                                             ! Timestep in the RSPMatch output file
      
      !=============================================================
      ! Make ACC Input file for RSPM09 from TGTMotion
      !============================================================= 
      strLenTgt = LEN_TRIM( tgtResFileName ) 
      strLenMotion = LEN_TRIM( motionFileName )-4

      ! Remove any files in folder that are from previous runs
      INQUIRE (file='ALL_FILES.INP', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('ALL_FILES.INP')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file='Run1.INP', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('Run1.INP')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file='Run2.INP', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('Run2.INP')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file='Run3.INP', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('Run3.INP')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file='Run4.INP', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('Run4.INP')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT.ACC', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT.ACC')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=tgtResFileName(1:strLenTgt)//'_TGT(g).tgt', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (tgtResFileName(1:strLenTgt)//'_TGT(g).tgt')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN1.acc', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN1.acc')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN2.acc', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN2.acc')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN3.acc', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN3.acc')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN4.acc', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN4.acc')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN1.rsp', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN1.rsp')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN2.rsp', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN2.rsp')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN3.rsp', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN3.rsp')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN4.rsp', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN4.rsp')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN1.unm', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN1.unm')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN2.unm', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN2.unm')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN3.unm', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN3.unm')
            CALL system(cmdString)
      ENDIF
      INQUIRE (file=motionFileName(1:strLenMotion)//'_TGT_RUN4.unm', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM (motionFileName(1:strLenMotion)//'_TGT_RUN4.unm')
            CALL system(cmdString)
      ENDIF
      ! end deleting previous files.

      OPEN(UnInTGT, file=motionFileName(1:strLenMotion)//'_TGT.ACC', action="WRITE", status="REPLACE")
      WRITE(UnInTGT, '(A)') 'ACC FILE USED FOR TARGET SPECTRUM MATCH PROGRAM'
      WRITE(UnInTGT, '(I11,X,F9.5,X,I1)') npts,  DT,  1     ! one being the added points at the start for padding to take place in every pass
      DO iTimeStep = 1, npts                                ! Loop through all time steps and write out the 
            WRITE(UnInTGT, '( 8e13.5)') accel(dirIndex,iTimeStep)/Gravity
      ENDDO
      CLOSE(UnInTGT)
      CALL WrScr1('Done Making ACC Input file for RSPM09')
      
      !=============================================================
      ! Open Target Response File and Convert PSA from m/s^2 to g for RSPM09
      !============================================================= 
      CALL OpenFInpFile( UnInTGT, tgtResFileName(1:strLenTgt) )
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,*) nFreq , nDamp
      ALLOCATE(dampValues(nDamp))
      ALLOCATE(freqValues(nFreq))
      ALLOCATE(time1(nFreq))
      ALLOCATE(time2(nFreq))
      ALLOCATE(PSA(nDamp,nFreq))
      ALLOCATE(PSAg(nDamp,nFreq))
      READ(UnInTGT,*) (dampValues(I),I=1,nDamp)
      DO I = 1, nFreq
            READ(UnInTGT,*) freqValues(I),time1(I), time2(I), (PSA(J,I),J = 1,nDamp)
      ENDDO
      CLOSE(UnInTGT)

      ! Convert PSA from m/s^2 into g units
      DO I = 1, nFreq
            DO J = 1, nDamp
                  PSAg(J,I) = PSA(J,I)/Gravity
            ENDDO
      ENDDO
      CALL WrScr1 ( 'Done opening target response file')
      ! Make input target response file in g units
      OPEN(UnInTGT,file=tgtResFileName(1:strLenTgt)//'_TGT(g).tgt', action="WRITE", status="REPLACE")
      WRITE(UnInTGT,'(A)') cmdString
      WRITE(UnInTGT,'(I10,X,I10)') nFreq , nDamp
      WRITE(UnInTGT,*) (dampValues(I), I=1, nDamp)
      DO I = 1, nFreq
            WRITE(UnInTGT,*) freqValues(I),time1(I), time2(I), (PSAg(J,I), J = 1, nDamp)
      ENDDO
      CLOSE(UnInTGT)
      
      CALL WrScr1( 'Done converting units for PSA')
      !=============================================================
      ! Create a text file with the input parameters to rspm09.exe so the user
      ! is not prompted for input.
      !============================================================= 
      OPEN(UnInTGT, file='runrspm.txt', action="WRITE", status="REPLACE")
      WRITE(UnInTGT,'(A)') 'ALL_FILES.INP'
      WRITE(UnInTGT,'(A)') ''
      CLOSE(UnInTGT)

      !=============================================================
      ! Make The Overall Input File used to Run each pass in RSPM09
      !============================================================= 
      OPEN(UnInTGT, file='ALL_FILES.INP', action="WRITE", status="REPLACE")
      WRITE(UnInTGT, '(I1)') nPass
      WRITE(UnInTGT, '(A)') 'Run1.INP'
      WRITE(UnInTGT, '(A)') 'Run2.INP'
      WRITE(UnInTGT, '(A)') 'Run3.INP'
      WRITE(UnInTGT, '(A)') 'Run4.INP'
      CLOSE(UnInTGT)
      CALL WrScr1('Done making All-Files')
      
      !=============================================================
      ! Check if file paths are readable by rspm
      !=============================================================
      tgtLen = LEN_TRIM(tgtResFileName(1:strLenTgt)//'_TGT(g).tgt')
      motionLen = LEN_TRIM(motionFileName(1:strLenMotion)//'_TGT.acc')
      accLen = LEN_TRIM(motionFileName(1:strLenMotion)//'_TGT_RUN1.acc')
      IF ((tgtLen .GE. 60) .OR. (motionLen .GE. 60) .OR. (accLen .GE. 60)) THEN
            CALL ProgAbort( 'Path length of the motion files and target spectrum must be shorter than 60 characters for the rspm program to work.' )
      ENDIF
      
      !=============================================================
      ! Make RUN1.INP FILE for the first Pass  (IF more passes >4 are needed, more files need to be created here but default is 4 for nPass)
      !=============================================================
      run1FreqMatch = FreqMatch12(1)*8                                      ! Make the first frequency match smaller for next pass 
      
      OPEN(UnInTGT, file='RUN1.INP', action="WRITE", status="REPLACE")
      WRITE(UnInTGT,'(I3)') MaxIter                                         ! Maximum number of iterations for 1st pass
      WRITE(UnInTGT,'(F6.3)') SMTol                                         ! Tolerance for maximum mismatch (in fraction of target)
      WRITE(UnInTGT,'(I1)')  1                                              ! convergence damping
      WRITE(UnInTGT,'(I1)')  7                                              ! model 7 which applies tapered cosine functions
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,F6.3,X,F6.3)') 1.25 , 0.25 , 1.0 , 4.0  ! alpha model, a1, a2, f1, f2 (not used for model 7)
      WRITE(UnInTGT,'(I2,X,F6.3)') 2, 0.0                                   ! scale falg, scale period to PGA  (=0 no, =1 yes, =2 yes but once) just for first pass
      WRITE(UnInTGT,'(I1)') 1                                               ! no interpolation required (Already interpolated motion according to DT)
      WRITE(UnInTGT,'(F7.4)') 1.0e-04                                       ! Minimum Eigenvalue
      WRITE(UnInTGT,'(I2)') 30                                              ! Group Size
      WRITE(UnInTGT,'(F6.3)') MaxFreq                                       ! Max Frequency
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,I1)') 0.0, 0.0, 4                       ! fBand, nPole (deactivated for model 7)
      WRITE(UnInTGT,'(I1)') 0                                               ! Mod PGA (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(I1,X,F6.3)') 0, 0.0                                   ! randomize target? (iSeed, ranFactor)
      WRITE(UnInTGT,'(F6.3,X,F6.3)') run1FreqMatch, FreqMatch12(2)          ! freqMatch
      WRITE(UnInTGT,'(I1)') 0                                               ! baseline cor flag (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(F6.3)') 1.0                                           ! scale factor
      WRITE(UnInTGT,'(A)') tgtResFileName(1:strLenTgt)//'_TGT(g).tgt'       ! Target Response File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT.acc'       ! Motion File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN1.acc'  ! Motion File created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN1.rsp'  ! Response Spectra File name created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN1.unm'  ! unmatched response spectrum file
      CLOSE(UnInTGT)
      CALL WrScr1( 'Done making RUN1.INP')

      !=============================================================
      ! Make RUN2.INP FILE for the second Pass
      !=============================================================
      run2FreqMatch = FreqMatch12(1)*4                                      ! Make the first frequency match smaller for next pass 

      OPEN(UnInTGT,file='RUN2.INP',action="WRITE",status="REPLACE")
      WRITE(UnInTGT,'(I3)') MaxIter                                         ! Maximum number of iterations for 1st pass
      WRITE(UnInTGT,'(F6.3)') SMTol                                         ! Tolerance for maximum mismatch (in fraction of target)
      WRITE(UnInTGT,'(I1)')  1                                              ! convergence damping
      WRITE(UnInTGT,'(I1)')  7                                              ! model 7 which applies tapered cosine functions
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,F6.3,X,F6.3)') 1.25 , 0.25 , 1.0 , 4.0  ! alpha model, a1, a2, f1, f2 (not used for model 7)
      WRITE(UnInTGT,'(I2,X,F6.3)') 0, 0.0                                   ! scale falg, scale period to PGA  (=0 no, =1 yes, =2 yes but once) just for first pass
      WRITE(UnInTGT,'(I1)') 1                                               ! no interpolation required (Already interpolated motion according to DT)
      WRITE(UnInTGT,'(F7.4)') 1.0e-04                                       ! Minimum Eigenvalue
      WRITE(UnInTGT,'(I2)') 30                                              ! Group Size
      WRITE(UnInTGT,'(F6.3)') MaxFreq                                       ! Max Frequency
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,I1)') 0.0, 0.0, 4                       ! fBand, nPole (deactivated for model 7)
      WRITE(UnInTGT,'(I1)') 0                                               ! Mod PGA (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(I1,X,F6.3)') 0, 0.0                                   ! randomize target? (iSeed, ranFactor)
      WRITE(UnInTGT,'(F6.3,X,F6.3)') run2FreqMatch, FreqMatch12(2)          ! freqMatch
      WRITE(UnInTGT,'(I1)') 0                                               ! baseline cor flag (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(F6.3)') 1.0                                           ! scale factor
      WRITE(UnInTGT,'(A)') tgtResFileName(1:strLenTgt)//'_TGT(g).tgt'       ! Target Response File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN1.acc'  ! Motion File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN2.acc'  ! Motion File created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN2.rsp'  ! Response Spectra File name created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN2.unm'  ! unmatched response spectrum file
      CLOSE(UnInTGT)
      CALL WrScr1('Done making RUN2.INP')

      !=============================================================
      ! Make RUN3.INP FILE for the third Pass
      !=============================================================
      run3FreqMatch = FreqMatch12(1)*2

      OPEN(UnInTGT,file='RUN3.INP',action="WRITE",status="REPLACE")
      WRITE(UnInTGT,'(I3)') MaxIter                                         ! Maximum number of iterations for 1st pass
      WRITE(UnInTGT,'(F6.3)') SMTol                                         ! Tolerance for maximum mismatch (in fraction of target)
      WRITE(UnInTGT,'(I1)')  1                                              ! convergence damping
      WRITE(UnInTGT,'(I1)')  7                                              ! model 7 which applies tapered cosine functions
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,F6.3,X,F6.3)') 1.25 , 0.25 , 1.0 , 4.0  ! alpha model, a1, a2, f1, f2 (not used for model 7)
      WRITE(UnInTGT,'(I2,X,F6.3)') 0, 0.0                                   ! scale falg, scale period to PGA  (=0 no, =1 yes, =2 yes but once) just for first pass
      WRITE(UnInTGT,'(I1)') 1                                               ! no interpolation required (Already interpolated motion according to DT)
      WRITE(UnInTGT,'(F7.4)') 1.0e-04                                       ! Minimum Eigenvalue
      WRITE(UnInTGT,'(I2)') 30                                              ! Group Size
      WRITE(UnInTGT,'(F6.3)') MaxFreq                                       ! Max Frequency
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,I1)') 0.0, 0.0, 4                       ! fBand, nPole (deactivated for model 7)
      WRITE(UnInTGT,'(I1)') 0                                               ! Mod PGA (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(I1,X,F6.3)') 0, 0.0                                   ! randomize target? (iSeed, ranFactor)
      WRITE(UnInTGT,'(F6.3,X,F6.3)') run3FreqMatch, FreqMatch12(2)          ! freqMatch
      WRITE(UnInTGT,'(I1)') 0                                               ! baseline cor flag (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(F6.3)') 1.0                                           ! scale factor
      WRITE(UnInTGT,'(A)') tgtResFileName(1:strLenTgt)//'_TGT(g).tgt'       ! Target Response File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN2.acc'  ! Motion File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN3.acc'  ! Motion File created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN3.rsp'  ! Response Spectra File name created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN3.unm'  ! unmatched response spectrum file
      CLOSE(UnInTGT)
      CALL WrScr1( 'Done making RUN3.INP')

      !=============================================================
      ! Make RUN4.INP FILE for the last Pass
      !=============================================================
      run4FreqMatch = FreqMatch12(1)

      OPEN (UnInTGT,file='RUN4.INP',action="WRITE",status="REPLACE")
      WRITE(UnInTGT,'(I3)') MaxIter                                         ! Maximum number of iterations for 1st pass
      WRITE(UnInTGT,'(F6.3)') SMTol                                         ! Tolerance for maximum mismatch (in fraction of target)
      WRITE(UnInTGT,'(I1)')  1                                              ! convergence damping
      WRITE(UnInTGT,'(I1)')  7                                              ! model 7 which applies tapered cosine functions
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,F6.3,X,F6.3)') 1.25 , 0.25 , 1.0 , 4.0  ! alpha model, a1, a2, f1, f2 (not used for model 7)
      WRITE(UnInTGT,'(I2,X,F6.3)') 0, 0.0                                   ! scale falg, scale period to PGA  (=0 no, =1 yes, =2 yes but once) just for first pass
      WRITE(UnInTGT,'(I1)') 1                                               ! no interpolation required (Already interpolated motion according to DT)
      WRITE(UnInTGT,'(F7.4)') 1.0e-04                                       ! Minimum Eigenvalue
      WRITE(UnInTGT,'(I2)') 30                                              ! Group Size
      WRITE(UnInTGT,'(F6.3)') MaxFreq                                       ! Max Frequency
      WRITE(UnInTGT,'(F6.3,X,F6.3,X,I1)') 0.0, 0.0, 4                       ! fBand, nPole (deactivated for model 7)
      WRITE(UnInTGT,'(I1)') 0                                               ! Mod PGA (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(I1,X,F6.3)') 0, 0.0                                   ! randomize target? (iSeed, ranFactor)
      WRITE(UnInTGT,'(F6.3,X,F6.3)') run4FreqMatch, FreqMatch12(2)          ! freqMatch
      WRITE(UnInTGT,'(I1)') 0                                               ! baseline cor flag (1=yes) (deactivated for model 7)
      WRITE(UnInTGT,'(F6.3)') 1.0                                           ! scale factor
      WRITE(UnInTGT,'(A)') tgtResFileName(1:strLenTgt)//'_TGT(g).tgt'       ! Target Response File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN3.acc'  ! Motion File name used for matching
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN4.acc'  ! Motion File created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN4.rsp'  ! Response Spectra File name created during first pass
      WRITE(UnInTGT,'(A)') motionFileName(1:strLenMotion)//'_TGT_RUN4.unm'  ! unmatched response spectrum file
      CLOSE(UnInTGT)
      CALL WrScr1( 'Done making RUN4.INP')
      !=============================================================
      ! CALL RSPM09 Program
      !=============================================================
      CALL system('rspm09.exe < runrspm.txt')
      !=============================================================
      ! See if Iterations were enough
      !=============================================================
      CALL OpenFInpFile( UnInTGT, MotionFileName(1:strLenMotion)//'_TGT_RUN4.unm' )
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(A)') cmdString
      unmStrLen = LEN_TRIM( cmdString )
      IF (unmStrLen .EQ. 59) THEN
            CALL ProgAbort ( 'Insufficient iterations to match the response.  Try increasing MaxIter in the seismic configuration file.' )
      ENDIF  
      CLOSE (UnInTGT)
      
      
      !=============================================================
      ! delete text file after you're done with rspm09
      !=============================================================
      INQUIRE(file='runrspm.txt', exist=fileFound)
      IF (fileFound) THEN
            WRITE(cmdString, '("del ", A)' ) TRIM ('runrspm.txt')
            CALL system(cmdString)
      ENDIF

      !=============================================================
      ! Get new motion Vector and convert it into m/s^2 again
      !=============================================================
      CALL OpenFInpFile( UnInTGT, MotionFileName(1:strLenmotion)//'_TGT_RUN4.acc' )
      READ(UnInTGT,'(A)') cmdString
      READ(UnInTGT,'(I10,X,F7.4,X,I10)') outStep , DT1 , nAdd
      
      numLines = CEILING(REAL(outStep)/5)
      ALLOCATE(tgtMotionMat(NumLines,5))
      
      ! Initialize the matrix
      DO I = 1, numLines
            DO J = 1, 5
                  tgtMotionMat(I,J) = 0
            ENDDO
      ENDDO
      ! Read the RSPMatch output file
      DO I=1,numLines
            READ(UnInTGT,*, IOSTAT=IOS)  (tgtMotionMat(I,J), J = 1,5)
      ENDDO
      ! Close the RSPMatch output file
      CLOSE (UnInTGT)
      ! Allocate the matched acceleration output
      ALLOCATE(matchedAccel(outStep))
      
      ! Making a vector from the acceleration data matrix
      DO J = 1,outStep
            curLine = CEILING(REAL(J)/5)
            curPos = MOD(J,5)
            IF (curPos .EQ. 0) THEN
                  curPos = 5
            ENDIF
            matchedAccel(J) = tgtMotionMat(curLine,curPos)*Gravity
      ENDDO
      
      ! Deallocate arrays for next direction
      DEALLOCATE(dampValues)
      DEALLOCATE(freqValues)
      DEALLOCATE(time1)
      DEALLOCATE(time2)
      DEALLOCATE(PSA)
      DEALLOCATE(PSAg)
      DEALLOCATE(tgtMotionMat)

      END SUBROUTINE TargetResonse

      !=======================================================================
      ! Synthetic Motion Subroutine
      !=======================================================================
      SUBROUTINE SynMotionGen( synRMSAmpValue, synDuration, randomSeeds, initRampFile, finalRampFile, ActFreq, timeStep, synAccel)
      USE                                 FFT_Module
      USE                                 Precision
      IMPLICIT                            NONE

      REAL(ReKi), INTENT(IN)              :: synRMSAmpValue                                  ! The synthetic motion duration in seconds
      REAL(ReKi), INTENT(IN)              :: synDuration                                     ! The synthetic motion duration in seconds
      INTEGER(4), INTENT(IN)              :: randomSeeds(2)                                  ! Seed to allow reproducible synthetic motion
      CHARACTER(128), INTENT(IN)          :: initRampFile                                    ! The filename of the initial ramp data
      CHARACTER(128), INTENT(IN)          :: finalRampFile                                   ! The filename of the final ramp data
      REAL(ReKi), INTENT(IN)              :: ActFreq                                         ! The actuator frequency
      REAL(ReKi), INTENT(IN)              :: timeStep
      REAL(ReKi), INTENT(OUT), ALLOCATABLE :: synAccel(:)                                    ! The synthetic acceleration time history
      
      INTEGER(4)                          :: Sttus                                           ! Status of allocation attempts.
      INTEGER(4)                          :: I                                               ! Generic index
      INTEGER(4)                          :: synTimeSteps                                    ! The number of time steps in the synthetic motion
      INTEGER(4)                          :: compTimeSteps                                   ! The number of time steps for the complex matrix
      REAL(ReKi)                          :: synDOmega                                       ! Frequency step (rad/s)
      REAL(ReKi)                          :: cornerFreq                                      ! Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 2/3 the Actuator Frequency
      REAL(ReKi)                          :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
      COMPLEX                             :: WGNC                                            ! Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component
      REAL(ReKi)                          :: WGNC_FactSyn                                    ! Factor used to scale the magnitude of the WGNC      as required by the discrete time inverse Fourier transform (-)
      REAL(ReKi)                          :: S2Sd_FactSyn                                    ! Factor used to scale the magnitude of the WaveS2Sdd as required by the discrete time inverse Fourier transform (-)
      REAL(ReKi)                          :: C1                                              ! Intermediate variable
      REAL(ReKi)                          :: C2                                              ! Intermediate variable
      REAL(ReKi)                          :: U1                                              ! First uniformly distributed random
      REAL(ReKi)                          :: U2                                              ! Second uniformly distributed random
      
      REAL(ReKi), ALLOCATABLE             :: SynCompMotion(:)                                ! Acceleration at the tower base (m-s^2)
      REAL(ReKi), ALLOCATABLE             :: AbsSynCompMotion(:)                             ! The absolute value of acceleration at the tower base (m-s^2)
      REAL(ReKi), ALLOCATABLE             :: FilteredSynCompMotion(:)                        ! The filtered value of acceleration at the tower base (m-s^2)
      COMPLEX, ALLOCATABLE                :: SynCompMotionC(:)                               ! Fourier transform of the acceleration at the tower base
      REAL(ReKi)                          :: synCompMotionRange(2)                           ! The min and max value to transform the random values to a uniform normal series
      REAL(ReKi)                          :: centerValue                                     ! The center value to transform the random values to a uniform normal series
      INTEGER(4)                          :: strLen
      REAL(ReKi)                          :: initRampStartTime
      INTEGER(4)                          :: initRampStartIndex = 1
      INTEGER(4)                          :: numInitRampPoints = 2
      REAL(ReKi), ALLOCATABLE             :: initRampData(:)                                 ! The initial ramp
      REAL(ReKi)                          :: finalRampStartTime
      INTEGER(4)                          :: finalRampStartIndex = 1
      INTEGER(4)                          :: numFinalRampPoints = 2
      REAL(ReKi), ALLOCATABLE             :: finalRampData(:)                                ! The initial ramp
      REAL(ReKi)                          :: CurSynRMSAmpValue                               ! The RMS value of the generated motion

      
      synTimeSteps = CEILING(synDuration/timeStep)+1
      
      !-----------------------------------------------------------------------
      ! Generate random numbers according to duration and random seeds 
      !-----------------------------------------------------------------------
      CALL WrScr(' Generating synthetic motion.')
      CALL RANDOM_SEED( PUT=randomSeeds(1:2) )
      IF ( MOD(synTimeSteps,2) == 1 ) THEN
            synTimeSteps = synTimeSteps + 1      ! Make timestep number even
      ENDIF
      
      compTimeSteps = MAX(synTimeSteps/2, 1 )                                    ! number of data for complex matrix
      synTimeSteps  = 2*PSF(compTimeSteps, 9 )                                   ! Convert to its primes 
      compTimeSteps = synTimeSteps/2+1                                           ! Update the value of numTimeStepsData2 based on the value needed for numTimeStepsData.
      
      synDOmega = TwoPi/(synTimeSteps*timeStep)                                        ! Compute the frequency step for synthetic motion calculations.
      cornerFreq = ActFreq*(2.0/3.0)
      Alpha = EXP((-timeStep)*(cornerFreq))
  
      ! Calculate the factors needed by the discrete time inverse Fourier
      ! transform in the calculations of the White Gaussian Noise (WGN) and
      ! the two-sided power spectral density of the time history per unit time:
      WGNC_FactSyn = SQRT( Pi/(synDOmega*timeStep) )                        ! This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN process has unit variance
      S2Sd_FactSyn = 1.0/timeStep                                     ! This factor is also needed by the discrete time inverse Fourier transform

      ! Allocate storage
      ALLOCATE(SynCompMotion(synTimeSteps) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
            CALL ProgAbort('Error allocating memory for the SynCompMotion array.')
      ENDIF
      ALLOCATE(AbsSynCompMotion(synTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
            CALL ProgAbort(' Error allocating memory for the AbsSynCompMotion array.')
      ENDIF
   
      ALLOCATE(SynCompMotionC(compTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
            CALL ProgAbort(' Error allocating memory for the SynCompMotionC array.')
      ENDIF
      ALLOCATE(FilteredSynCompMotion(synTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
            CALL ProgAbort('Error allocating memory for the FilteredSynCompMotion array.')
      ENDIF
      
      ALLOCATE(synAccel(synTimeSteps), STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
            CALL ProgAbort('Error allocating memory for the synAccel array.')
      ENDIF
      ! Compute the positive-frequency components (including zero) of the Fourier
      ! transforms of the synthetic motion
      DO I = 1,compTimeSteps                                    ! Loop through the positive frequency components (including zero) of the Fourier transforms
            ! Compute the Fourier transform of the realization of a White Gaussian Noise
            ! (WGN) time series process with unit variance:
            ! NOTE: For the time series process to be real with zero mean, the values at
            ! Omega == 0.0 and Omega == compTimeSteps must be zero.
            IF ( ( I == 1 ) .OR. ( I == compTimeSteps ) )  THEN   ! .TRUE. if ( Omega == 0.0 ) or ( Omega == compTimeSteps)
                  WGNC = (0.0,0.0)
            ELSE                                                        ! All other Omega
                  U1 = 0.0
                  DO WHILE ( U1 == 0.0 )
                        CALL RANDOM_NUMBER(U1)
                  ENDDO
                  CALL RANDOM_NUMBER(U2)
                  ! Compute intermediate variables:
                  C1 = SQRT( -2.0*LOG(U1) )
                  C2 = TwoPi*U2
                  ! Compute the unit normal randoms using boxmuller:
                  WGNC = CMPLX( C1*COS(C2), C1*SIN(C2) )
            ENDIF
            ! Compute the Fourier transform of the synthetic Motion
            SynCompMotionC(I) = ( WGNC_FactSyn*WGNC )*SQRT( TwoPi*( S2Sd_FactSyn ) )
      ENDDO
      
      ! Compute the inverse Fourier transforms to find the time-domain
      ! representations of the Synthetic Motion: 
      CALL InitFFT( synTimeSteps, .FALSE. )
      CALL ApplyFFT_cx( SynCompMotion(:), SynCompMotionC(:) )
      CALL ExitFFT
      
      !-----------------------------------------------------------
      ! Normalize random motion between -1 to 1
      !-----------------------------------------------------------
      ! Find maximum and minimum amplitudes
      synCompMotionRange(1) = SynCompMotion(1)
      synCompMotionRange(2) = SynCompMotion(1)
      DO I = 1,synTimeSteps  ! Loop through all time steps
            IF (SynCompMotion(I) .LT. synCompMotionRange(1) ) THEN
                  synCompMotionRange(1) = SynCompMotion(I)
            ELSEIF ( SynCompMotion(I) .GT. synCompMotionRange(2)  ) THEN
                  synCompMotionRange(2) = SynCompMotion(I)
            ENDIF
      ENDDO        
      centerValue = (synCompMotionRange(1) + synCompMotionRange(2))/2
      ! Bring Motion to zero
      DO I = 1,synTimeSteps  ! Loop through all time steps
            SynCompMotion(I) = SynCompMotion(I) - centerValue
      ENDDO
      DO I = 1,synTimeSteps  ! Loop through all time steps
            AbsSynCompMotion(I) = ABS(SynCompMotion(I))
      ENDDO
      synCompMotionRange(2)  = AbsSynCompMotion(1)
      DO I = 1,synTimeSteps  ! Loop through all time steps
            IF( AbsSynCompMotion(I) .GE. synCompMotionRange(2)) THEN 
                  synCompMotionRange(2) = AbsSynCompMotion(I)
            ENDIF  
      ENDDO
      DO I = 1,synTimeSteps  ! Loop through all time steps
            SynCompMotion(I) = SynCompMotion(I)/synCompMotionRange(2)
      ENDDO
      
      !-----------------------------------------------------------
      ! Filter the actuator frequency content
      !-----------------------------------------------------------
      FilteredSynCompMotion(1) = 0.0
      DO I = 2,synTimeSteps                                  ! Loop through all time steps
            FilteredSynCompMotion(I) = (1-Alpha)*SynCompMotion(I)+(Alpha)*(FilteredSynCompMotion(I-1))
      ENDDO
      
      ! Read the ramp files
      ! Check if the file is specified
      strLen = LEN_TRIM(initRampFile)
      IF (strLen .GT. 1) THEN
            CALL ReadMotionFile(initRampFile(1:strLen), 1.0, initRampData, initRampStartTime)
      ELSE ! Provide a dummy array of 1s if a ramp is not specified
            ALLOCATE(initRampData(2), STAT=Sttus )
            IF ( Sttus /= 0 ) THEN
                  CALL ProgAbort('Error allocating memory for the initRampData array.')
            ENDIF
            initRampData(1) = 1
            initRampData(2) = 1
            initRampStartTime = 0
      ENDIF
      
      strLen = LEN_TRIM(finalRampFile)
      IF (strLen .GT. 1) THEN
            CALL ReadMotionFile(finalRampFile(1:strLen), 1.0, finalRampData, finalRampStartTime)
      ELSE! Provide a dummy array of 1s if a ramp is not specified
            ALLOCATE(finalRampData(2), STAT=Sttus )
            IF ( Sttus /= 0 ) THEN
                  CALL ProgAbort('Error allocating memory for the finalRampData array.')
            ENDIF
            finalRampData(1) = 1
            finalRampData(2) = 1
            finalRampStartTime = (synTimeSteps-1)*timeStep
      ENDIF

      initRampStartIndex = MAX(FLOOR(initRampStartTime/timeStep),1)
      numInitRampPoints = SIZE(initRampData)
      finalRampStartIndex = MAX(CEILING(finalRampStartTime/timeStep),1)
      ! Truncate the final ramp if it goes beyond the number of data points
      numFinalRampPoints = MIN(SIZE(finalRampData),synTimeSteps - finalRampStartIndex+1)
      
      IF (initRampStartIndex + numInitRampPoints - 1 .GT. finalRampStartIndex) THEN
            CALL ProgAbort(' Error : Start time in final ramp should be greater than end time in initial ramp file')
      ENDIF
 
      IF (synTimeSteps .GT. finalRampStartIndex+numFinalRampPoints - 1) THEN
            CALL WrScr('WARNING: Synthetic data extends beyond the final ramp.')
      ENDIF
      
      DO I = 1,synTimeSteps                                  ! Loop through all time steps
            !synAccel(I) = 0
            IF (I .LT. initRampStartIndex) THEN
                  ! We are before the initial ramp
                  synAccel(I) = FilteredSynCompMotion(I)
            ELSE IF (I .LT. initRampStartIndex + numInitRampPoints) THEN
                  ! We are in the initial ramp
                  synAccel(I) = FilteredSynCompMotion(I)*initRampData(I-initRampStartIndex+1)
            ELSE IF (I .LT. finalRampStartIndex) THEN
                  ! We are between the two ramps
                  synAccel(I) = FilteredSynCompMotion(I)
            ELSE IF (I .LT. finalRampStartIndex + numFinalRampPoints) THEN
                  ! We are in the final ramp
                  synAccel(I) = FilteredSynCompMotion(I)*finalRampData(I-finalRampStartIndex+1)
            ELSE
                  ! We are beyond the final ramp
                  synAccel(I) = FilteredSynCompMotion(I)
            ENDIF
            ! Calculate the RMS value of the generated motion
            IF (I .EQ. 1 ) THEN
                  CurSynRMSAmpValue =  synAccel(I)**2
            ELSE
                  CurSynRMSAmpValue = CurSynRMSAmpValue + synAccel(I)**2
            ENDIF
      ENDDO
      
      CurSynRMSAmpValue = SQRT(CurSynRMSAmpValue/synTimeSteps)
      
      !----------------------------------------------------------------------------
      ! Scale motion to the desired RMS value according to SynRMSAmpValue
      !----------------------------------------------------------------------------
      DO I = 1,synTimeSteps
            synAccel(I) = synAccel(I) * (SynRMSAmpValue/CurSynRMSAmpValue)
      ENDDO
      ! Clean up memory that was allocated locally
      DEALLOCATE(SynCompMotion)
      DEALLOCATE(AbsSynCompMotion)
      DEALLOCATE(SynCompMotionC)
      DEALLOCATE(FilteredSynCompMotion)
      DEALLOCATE(initRampData)
      DEALLOCATE(finalRampData)
      
      END SUBROUTINE SynMotionGen
      !=======================================================================
      SUBROUTINE WriteMotionFile(motionData, fileName, timeStep)
      
      IMPLICIT                            NONE
      INTEGER(4)                          :: Sttus
      
      REAL(ReKi), INTENT(IN)              :: motionData(:)                               ! The time history to be written to file
      CHARACTER(128), INTENT(IN)          :: fileName                                    ! The filename of the initial ramp data
      REAL(ReKi), INTENT(IN)              :: timeStep                                    ! The time step for the motion data.
      INTEGER(4)                          :: numDataPoints                               ! The number of data points in 
      INTEGER(4)                          :: fineNameLength
      INTEGER(4)                          :: UnOut = 35
      INTEGER(4)                          :: I
      
      numDataPoints = SIZE(motionData)
      fineNameLength = LEN_TRIM(fileName)
      !----------------------------------------------------------------------------
      ! Create output file
      !----------------------------------------------------------------------------
      OPEN (UnOut,file= fileName(1:fineNameLength),action="WRITE",status="REPLACE")
      DO I = 1,numDataPoints                              ! Loop through all time steps
          WRITE (UnOut,*) (I-1)*timeStep, motionData(I)
      ENDDO                                               ! I - All time steps
      CLOSE (UnOut)
      END SUBROUTINE WriteMotionFile
      !=======================================================================
      SUBROUTINE CleanupSeismic()

      IMPLICIT                            NONE
      INTEGER(4)                          :: Sttus

      IF( PtfmArraysAllocated ) THEN
            DEALLOCATE(PtfmAccel, STAT=Sttus )
            IF ( Sttus /= 0 )  THEN
                  CALL ProgAbort ( ' Error deallocating memory for the PtfmAccel array.' )
            ELSE
                  PtfmArraysAllocated = .FALSE.
            ENDIF
            IF (ActFreq /= 0 ) THEN
                  DEALLOCATE(PtfmVel, STAT=Sttus )
                  IF ( Sttus /= 0 )  THEN
                        CALL ProgAbort ( ' Error deallocating memory for the PtfmVel array.' )
                  ELSE
                        PtfmArraysAllocated = .FALSE.
                  ENDIF
                  DEALLOCATE(PtfmDisp, STAT=Sttus )
                  IF ( Sttus /= 0 )  THEN
                        CALL ProgAbort ( ' Error deallocating memory for the PtfmDisp array.' )
                  ELSE
                        PtfmArraysAllocated = .FALSE.
                  ENDIF
            ENDIF
      ENDIF
      END SUBROUTINE CleanupSeismic
      
END MODULE Seismic
!=======================================================================
!=======================================================================
! End of the supporting code required by the seismic version of 
! UserPtfmLd
!=======================================================================