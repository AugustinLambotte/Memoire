CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:25Z creation; 2020-11-17T12:19:03Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8(   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8h   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9,   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         90   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    98   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9<   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9D   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9L   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9T   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9X   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9`   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9d   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9h   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9l   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :l   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   Ml   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  W�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   `h   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  b�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   j�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  m   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  u�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   }�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   �|   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �d   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �h   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �l   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �p   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �t   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �8   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �8   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �8   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �8Argo profile    3.1 1.2 19500101000000  20200828144325  20220127170408  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               MA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @�,���j11   @�,٘�k`@Q��:���0���R�K1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 7.1 dbar]                                                      @���A��A��AffA1��A@  AP  Ac33Aq��A�  A�  A�  A�  A�  A�  A�  A���A���A�  A���A���A�33A�  A�A�  B   B  B  B  B  BffB��B33B ��B#��B(ffB,ffB0  B3��B8  B=33B@  BD��BG��BK��BPffBS��BX  B\  B`ffBd  Bg��Bk��Bo��Bt  BxffB{��B�  B���B�33B�33B�ffB�33B�  B�33B�33B�33B�33B�33B�33B�33B�33B���B���B���B�33B���B�33B�  B���B���B�33B�33B�33B���B�  B�ffB���B�33B�  B�  B���B�33BǙ�B�ffB�33B�ffB���B�ffB���B�  B���B�  B�  B�B���B�33B�  C� C� C� C��C	ffCL�C��C��C� C� CffCffCffCL�CffCffC!ffC#� C%� C'��C)33C+L�C-L�C/L�C1ffC3ffC5� C733C9L�C;� C=33C?ffCA� CC33CEL�CG� CI33CKL�CM� CO33CQffCS� CU33CWL�CY� C[33C]L�C_� Ca� Cc�3CeffCgffCi��CkL�CmffCo��CqL�Cs� Cu�3CwffCyffC{��C}L�CffC���C��fC�� C�ٚC�� C��fC���C��3C�ٚC���C��3C���C�� C��3C�ٚC���C��3C��fC���C���C�� C�� C��3C��3C��fC��fC��fC��fC��fC��fC��fC��fC��3C��3C�� C�� C���C���C���C��fC��3C�� C���C�ٚC��3C�� C�� C���C��fC��3C�� C�ٚC��fC��3C�� C��3C���C���C��fC��fC��3C�� C���C��fC��3C�� C�ٚCó3CČ�CŦfCƳ3Cǌ�CȦfC�� Cʙ�C˦fC�� CͦfC�� C���Cг3Cь�Cҳ3C���CԦfCՙ�Cֳ3C���Cس3Cٙ�Cڌ�C۳3C���Cݳ3CަfCߌ�C�� C��fC���C���C�3C�fC晚C��C� C�3C�ٚC���C�� C��3C�fCC�C��C� C�3C��fC��fC�ٚC�ٚC���C���C�� C��3D 33DffD�fD� D9�D��DٚD	�D
Y�D� D� D� D  DffD��D��DFfD�3D� D��DL�D�3D�3D��D&fDs3D � D!�3D#&fD$� D%�fD&��D(33D)y�D*�fD+�3D-@ D.y�D/�3D0�3D2,�D3s3D4��D5��D7FfD8��D9� D:�3D<&fD=y�D>��D@fDAFfDB�fDC�fDE�DFL�DGy�DH�fDI��DK@ DL� DM��DO�DPFfDQs3DR� DTfDU9�DVs3DW�3DX��DZFfD[�fD\�fD^�D_9�D`l�Da� Dc�DdL�De�fDf� Dg��Di33Djs3Dk�3Dl��Dn@ Do��Dp��Dq�fDs33Dt�fDu�3Dw�DxFfDy�fDz��D{��D},�D~y�D� D�y�D�#3D��3D�` D�  D���D�@ D���D�|�D��D�� D�VfD���D��3D�@ D�ٚD�vfD�3D�� D�Y�D�fD��fD�FfD��3D�� D�  D���D�\�D�3D���D�6fD���D���D�&fD��fD�i�D���D��3D�9�D��3D�� D��D���D�\�D�  D��fD�9�D�� D�y�D� D���D�c3D�  D���D�<�D���D��3D�  D���D�VfD��fD��3D�33D��3D�vfD��D��3D�Y�D��fD��3D�C3D�� D�� D�  D���D�Y�D���D��fD�9�D�ٚD�y�D�fD���D�Y�D���D�� D�C3D�ٚD�|�D�  D��fD�P D���D��fD�C3D�� D�y�D�fD��3D�\�D�fD�� D�<�D�ٚD�� D��D¶fD�S3D��3Dē3D�6fD�ٚDƃ3D��Dǹ�D�c3D�3Dɠ D�@ D��fD�y�D�  D̶fD�` D�	�DΣ3D�C3D�� D�|�D��DѼ�D�\�D���DӜ�D�C3D��fD�|�D�  D��fD�` D�ٚ111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@���A��A��AffA1��A@  AP  Ac33Aq��A�  A�  A�  A�  A�  A�  A�  A���A���A�  A���A���A�33A�  A�A�  B   B  B  B  B  BffB��B33B ��B#��B(ffB,ffB0  B3��B8  B=33B@  BD��BG��BK��BPffBS��BX  B\  B`ffBd  Bg��Bk��Bo��Bt  BxffB{��B�  B���B�33B�33B�ffB�33B�  B�33B�33B�33B�33B�33B�33B�33B�33B���B���B���B�33B���B�33B�  B���B���B�33B�33B�33B���B�  B�ffB���B�33B�  B�  B���B�33BǙ�B�ffB�33B�ffB���B�ffB���B�  B���B�  B�  B�B���B�33B�  C� C� C� C��C	ffCL�C��C��C� C� CffCffCffCL�CffCffC!ffC#� C%� C'��C)33C+L�C-L�C/L�C1ffC3ffC5� C733C9L�C;� C=33C?ffCA� CC33CEL�CG� CI33CKL�CM� CO33CQffCS� CU33CWL�CY� C[33C]L�C_� Ca� Cc�3CeffCgffCi��CkL�CmffCo��CqL�Cs� Cu�3CwffCyffC{��C}L�CffC���C��fC�� C�ٚC�� C��fC���C��3C�ٚC���C��3C���C�� C��3C�ٚC���C��3C��fC���C���C�� C�� C��3C��3C��fC��fC��fC��fC��fC��fC��fC��fC��3C��3C�� C�� C���C���C���C��fC��3C�� C���C�ٚC��3C�� C�� C���C��fC��3C�� C�ٚC��fC��3C�� C��3C���C���C��fC��fC��3C�� C���C��fC��3C�� C�ٚCó3CČ�CŦfCƳ3Cǌ�CȦfC�� Cʙ�C˦fC�� CͦfC�� C���Cг3Cь�Cҳ3C���CԦfCՙ�Cֳ3C���Cس3Cٙ�Cڌ�C۳3C���Cݳ3CަfCߌ�C�� C��fC���C���C�3C�fC晚C��C� C�3C�ٚC���C�� C��3C�fCC�C��C� C�3C��fC��fC�ٚC�ٚC���C���C�� C��3D 33DffD�fD� D9�D��DٚD	�D
Y�D� D� D� D  DffD��D��DFfD�3D� D��DL�D�3D�3D��D&fDs3D � D!�3D#&fD$� D%�fD&��D(33D)y�D*�fD+�3D-@ D.y�D/�3D0�3D2,�D3s3D4��D5��D7FfD8��D9� D:�3D<&fD=y�D>��D@fDAFfDB�fDC�fDE�DFL�DGy�DH�fDI��DK@ DL� DM��DO�DPFfDQs3DR� DTfDU9�DVs3DW�3DX��DZFfD[�fD\�fD^�D_9�D`l�Da� Dc�DdL�De�fDf� Dg��Di33Djs3Dk�3Dl��Dn@ Do��Dp��Dq�fDs33Dt�fDu�3Dw�DxFfDy�fDz��D{��D},�D~y�D� D�y�D�#3D��3D�` D�  D���D�@ D���D�|�D��D�� D�VfD���D��3D�@ D�ٚD�vfD�3D�� D�Y�D�fD��fD�FfD��3D�� D�  D���D�\�D�3D���D�6fD���D���D�&fD��fD�i�D���D��3D�9�D��3D�� D��D���D�\�D�  D��fD�9�D�� D�y�D� D���D�c3D�  D���D�<�D���D��3D�  D���D�VfD��fD��3D�33D��3D�vfD��D��3D�Y�D��fD��3D�C3D�� D�� D�  D���D�Y�D���D��fD�9�D�ٚD�y�D�fD���D�Y�D���D�� D�C3D�ٚD�|�D�  D��fD�P D���D��fD�C3D�� D�y�D�fD��3D�\�D�fD�� D�<�D�ٚD�� D��D¶fD�S3D��3Dē3D�6fD�ٚDƃ3D��Dǹ�D�c3D�3Dɠ D�@ D��fD�y�D�  D̶fD�` D�	�DΣ3D�C3D�� D�|�D��DѼ�D�\�D���DӜ�D�C3D��fD�|�D�  D��fD�` D�ٚ111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?�Ĝ?�?}?���?�?�^5?��+?��^?��?��T?���?�33?�M�?�r�?�\)?��>�`B>H�9�O�;���Ͽ]�cS��e�T�l1�����Ĝ��ff���7��/��hs������w���j������
=���˿�X�������^��r���r����ۿ�33��������öF��KǿǍP���1'�Ǯ�Ƈ+�š˿ě��öF�����녿�-��-��J��n���n��°!�����Ͽ�Z��t������;d��V��{��/���������׿�`B������%��w�t9X�|푿�  ��r������{"ѿbJ�xQ쿈�u�~vɿa���Qhs�$�/�(r�����dZ�*���R񪿄Z������ƨ����y�#�h1'�\(��LI��H�9�0�׿'+�!G����񪾢MӾ�u:�o>~��?�7?��?'�?-��?<�?N��?Vȴ?bM�?n��?q�?m�h?co?\j?Tz�?T�j?S33?Q�?F�y?��>��>��?��?:��?B�\??|�?<j?49X?�?�?(�9?,I�?#o>�l�>��>���>��^>��^>��7>��\>w��>t�j>p��>j~�>aG�>M��>)��>$�/>+> Ĝ>O�>�-> Ĝ>��>C�=�G�=�hs=y�#=u=}�=��=�w=8Q�=49X=49X=8Q�='�=#�
=��=t�=+=o<��<�/<�j<�C�<t�;�`B    �49X�49X�T����9X��P�P�`�Y��Y��e`B�u��%���������1��vɽ\�ě��ȴ9�ȴ9�\�\�\�����E����{��^5�ě��������ě�����"ѽ����"ѽ�G���`B��h����   �o���+�
=q�V�V�V�
=q���J�o�V�1'���#��xս�S���S���S���G���G���;d��;d��G���/��"ѽ�/��/��"ѽ�/��"ѽ�
=������
=������������������������ȴ9�������ͽ��ͽ��`��������������
=�������ȴ9�\��Q콴9X��9X��-�� Ž� Ž�1���T���w������O߽���u�aG��]/�<j�0 Ž<j�0 Ž'�P��P��P���t��\)�L�ͽixսixսq���8Q�C���o;��
;ě��t���j��j��j��o�t��D��;o<o<e`B<�9X=o=#�
=e`B=Y�=8Q�<�/<��=0 �='�<���<�/<�9X;�`B;o    ��o�#�
��C���9X������h�o�t��49X�L�ͽ]/�ixսq���y�#��o��7L��hs���㽥�T�� Ž�Q콼j�ě��ȴ9���`��/��G���F���#�J�+�	7L�C��V�hs�t��z��u���"��)��.{�0 ž2-�6E��:^5�=p��B�\�G��I�^�S�ϾZ��\(��`A��gl��k��p�׾t�j�vȴ�z�H�}󶾁%�������˾����7L������O߾�V���`��n��������+���u��"Ѿ�/���R��;d��A����徥`B��l����þ�~�������D��{�� ž������F��?}��E���Q쾹�#���H��dZ��p���|��%�Õ�����Ƨ�ȴ9��=q��C���I���V��hs��t�����������b����ܬ��5?��;d�߾w��A���G���MӾ�Z��ff���~���D��{��� ž�&��&��-����Q���H��j��� A��Mӿ������y��9�	xտ	�^�
~��
���C����1�V���{����;� ſbN��`��`��`�&�&�hs�hs����녿녿-���j��+�b��#�"ѿ��/�5?��R��ۿ Ĝ�"�\�#���$��$���%��%`B�%`B�%`B�%��%��%��%��$�/111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111?�Ĝ?�?}?���?�?�^5?��+?��^?��?��T?���?�33?�M�?�r�?�\)?��>�`B>H�9�O�;���Ͽ]�cS��e�T�l1�����Ĝ��ff���7��/��hs������w���j������
=���˿�X�������^��r���r����ۿ�33��������öF��KǿǍP���1'�Ǯ�Ƈ+�š˿ě��öF�����녿�-��-��J��n���n��°!�����Ͽ�Z��t������;d��V��{��/���������׿�`B������%��w�t9X�|푿�  ��r�G�O�G�O��bJ�xQ쿈�u�~vɿa���Qhs�$�/�(r�����dZ�*���R񪿄Z������ƨ����y�#�h1'�\(��LI��H�9�0�׿'+�!G����񪾢MӾ�u:�o>~��?�7?��?'�?-��?<�?N��?Vȴ?bM�?n��?q�?m�h?co?\j?Tz�?T�j?S33?Q�?F�y?��>��>��?��?:��?B�\??|�?<j?49X?�?�?(�9?,I�?#o>�l�>��>���>��^>��^>��7>��\>w��>t�j>p��>j~�>aG�>M��>)��>$�/>+> Ĝ>O�>�-> Ĝ>��>C�=�G�=�hs=y�#=u=}�=��=�w=8Q�=49X=49X=8Q�='�=#�
=��=t�=+=o<��<�/<�j<�C�<t�;�`B    �49X�49X�T����9X��P�P�`�Y��Y��e`B�u��%���������1��vɽ\�ě��ȴ9�ȴ9�\�\�\�����E����{��^5�ě��������ě�����"ѽ����"ѽ�G���`B��h����   �o���+�
=q�V�V�V�
=q���J�o�V�1'���#��xս�S���S���S���G���G���;d��;d��G���/��"ѽ�/��/��"ѽ�/��"ѽ�
=������
=������������������������ȴ9�������ͽ��ͽ��`��������������
=�������ȴ9�\��Q콴9X��9X��-�� Ž� Ž�1���T���w������O߽���u�aG��]/�<j�0 Ž<j�0 Ž'�P��P��P���t��\)�L�ͽixսixսq���8Q�C���o;��
;ě��t���j��j��j��o�t��D��;o<o<e`B<�9X=o=#�
=e`B=Y�=8Q�<�/<��=0 �='�<���<�/<�9X;�`B;o    ��o�#�
��C���9X������h�o�t��49X�L�ͽ]/�ixսq���y�#��o��7L��hs���㽥�T�� Ž�Q콼j�ě��ȴ9���`��/��G���F���#�J�+�	7L�C��V�hs�t��z��u���"��)��.{�0 ž2-�6E��:^5�=p��B�\�G��I�^�S�ϾZ��\(��`A��gl��k��p�׾t�j�vȴ�z�H�}󶾁%�������˾����7L������O߾�V���`��n��������+���u��"Ѿ�/���R��;d��A����徥`B��l����þ�~�������D��{�� ž������F��?}��E���Q쾹�#���H��dZ��p���|��%�Õ�����Ƨ�ȴ9��=q��C���I���V��hs��t�����������b����ܬ��5?��;d�߾w��A���G���MӾ�Z��ff���~���D��{��� ž�&��&��-����Q���H��j��� A��Mӿ������y��9�	xտ	�^�
~��
���C����1�V���{����;� ſbN��`��`��`�&�&�hs�hs����녿녿-���j��+�b��#�"ѿ��/�5?��R��ۿ Ĝ�"�\�#���$��$���%��%`B�%`B�%`B�%��%��%��%��$�/111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oA��hA��\A�I�A�
=A�+A�?}A��uA��/B 
=B �B /B �FB �HB�B33Bv�B|�B�BVBVBr�B�%B�LB�Bz�B�BŢB!�BJ�B`BB�B��B�ZB��B�B-BK�B\)Bk�B�7B��B�9BB��B��B�fB�yB�B��B	B	hB	�B	%�B	-B	5?B	=qB	@�B	F�B	T�B	^5B	`BB	ffB	t�B	}�B	�DB	�\B	��B	��B	��B	��B	�-B	�jB	��B	��B	��B	�fB	�B	�B	�B	��B	�B	�B	��B	��B
PB
JB
�B
�B
'�B
$�B
49B
2-B
>wB
?}B
F�B
J�B
G�B
N�B
R�B
YB
n�B
}�B
�1B
��B
��B
��B
�-B
�LB
ƨB
��B
�yB
��BVB"�B@�BA�BF�BQ�BS�B^5BcTBdZBm�Bo�Bn�Bn�Bn�BjBhsBhsBcTBaHBS�BP�BdZBiyBaHBu�Bt�BjBaHBr�BiyBq�Br�Bp�BiyB_;B^5B`BB_;BaHBaHBaHBbNBbNBcTBdZBe`Be`BdZBdZBcTBq�Be`BdZBffBdZBe`BbNBcTBe`BdZBffBffBhsBiyBhsBiyBjBjBjBjBjBjBjBjBjBk�Bk�BjBk�Bl�BjBjBjBiyBk�BjBjBjBjBjBiyBjBiyBiyBiyBhsBiyBiyBiyBjBjBjBl�Bl�Bk�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bl�Bm�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bk�Bk�Bl�Bl�Bm�Bm�Bn�Bn�Bm�Bm�Bn�Bo�Bp�Bq�Bq�Bp�Bq�Bq�Bq�Bq�Bq�Bq�Br�Br�Bq�Br�Br�Br�Br�Bs�Br�Bs�Bs�Bs�Bs�Bt�Bt�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bv�Bw�Bv�Bw�Bw�Bw�Bv�Bw�Bx�Bx�By�Bz�B{�B{�B|�B}�B~�B~�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�+B�DB�JB�PB�DB�=B�DB�JB�JB�VB�\B�bB�hB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A���A��qA�Q�A�`A�3*A�G�A���A��YB TB �B 34B �bB �gB�B7UBz�B�!B�BZ>BZFBv�B�jB��B�SB+B�QB��B&BOBd�B�ZB�DB�BLB�B1iBP"B`�Bo�B��B�EB��B��B�B�JB��B��B��B�"B	xB	�B	 B	*CB	1qB	9�B	A�B	D�B	KB	Y_B	b�B	d�B	j�B	yB	�WB	��B	��B	��B	�B	�:B	�\B	��B	��B	��B	�'B	�WB	��B	��B	�B	�B	�$B	�B	�G�O�G�O�B
�B
�B
�B

B
,[B
)HB
8�B
6�B
B�B
C�B
KB
O+B
LB
SEB
W^B
]�B
sB
�aB
��B
��B
�B
�TB
��B
��B
�B
�nB
��B
�IB�B'EBD�BFBKBVbBXnBb�Bg�Bh�BrBtBsBsBsBn�Bl�Bl�Bg�Be�BXoBU\Bh�Bm�Be�Bz:By5Bn�Be�Bw'Bm�Bv#Bw&BuBm�Bc�Bb�Bd�Bc�Be�Be�Be�Bf�Bf�Bg�Bh�Bi�Bi�Bh�Bh�Bg�Bv"Bi�Bh�Bj�Bh�Bi�Bf�Bg�Bi�Bh�Bj�Bj�Bl�Bm�Bl�Bm�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bn�Bo�BqBn�Bn�Bn�Bm�Bo�Bn�Bn�Bn�Bn�Bn�Bm�Bn�Bm�Bm�Bm�Bl�Bm�Bm�Bm�Bn�Bn�Bn�BqBqBo�BqBqBqBqBrBrBp�BqBq BrBqBqBq BqBo�BqBqBo�Bo�BqBqBrBrBsBsBrBr	BsBtBuBv!Bv!BuBv!BvBv!BvBv!Bv!Bw&Bw&BvBw%Bw&Bw&Bw%Bx,Bw&Bx,Bx-Bx,Bx.By4By4Bx.By5By1By4By5By2By0By4By2Bz<Bz:Bz7Bz7B{AB|EB{@B|CB|EB|EB{?B|FB}IB}KB~SBWB�`B�]B�eB�lB�pB�rB�vB�yB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�B�B�B�
B�	B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�$B�B� B�#B�"B�"B�#B�#B�%B�!B�#B�#B�#B�"B�$B�*B�+B�)B�)B�)B�)B�*B�+B�1B�/B�.B�/B�/B�*B�/B�0B�/B�0B�/B�/B�0B�0B�-B�0B�1B�0B�6B�0B�1B�6B�6B�7B�5B�5B�.B�7B�8B�5B�.B�5B�/B�6B�4B�;B�6B�6B�>B�4B�4B�6B�<B�7B�9B�7B�;B�9B�;B�<B�<B�AB�=B�<B�<B�:B�8B�8B�:B�;B�BB�EB�;B�=B�9B�CB�;B�FB�FB�CB�>B�CB�AB�<B�CB�>B�DB�DB�:B�DB�DB�DB�BB�DB�BB�BB�BB�EB�EB�CB�CB�CB�CB�EB�CB�CB�CB�AB�DB�BB�JB�KB�BB�GB�HB�JB�HB�HB�KB�HB�KB�HB�JB�IB�IB�KB�IB�FB�FB�FB�HB�IB�FB�HB�IB�IB�FB�HB�IB�IB�IB�IB�GB�OB�OB�PB�QB�PB�MB�NB�NB�OB�TB�OB�PB�TB�LB�LB�NB�NB�OB�NB�NB�OB�O111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0043345 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040820220127170408  IF  ARFMCODA035h                                                                20200828144325                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144433  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200828144433  QCF$                G�O�G�O�G�O�0000000000004000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            @���D�ٚG�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170408  IP  PSAL            @���D�ٚG�O�                