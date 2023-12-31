CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  P   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:23Z creation; 2020-11-17T12:19:02Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        	@  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  C�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	@  F    PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  O@   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	@  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  Z�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  d   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  f`   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  o�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  q�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  {0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �p   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  �P   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �<   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �L   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �P   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �`   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �d   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �h   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �l   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828144323  20220127170401  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               )A   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @���So�1   @���A��@Sbr1ǲ�� 7��1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 5.8 dbar]                                                      @�33@�33A��A��A��A.ffAA��AQ��A`  AnffA�  A���A���A�ffA�  A���A�  A���A���A�33A���A�  A�33A���A���A�ffA�ffB��B��B��B  BffB��B  B ffB#��B(  B,��B0  B333B6ffB;��B@ffBC��BF��BK��BPffBS33BV��B[33B`  Bc33Bh  Bl��Bp  Bs33BvffB{��B�33B���B�ffB���B�33B���B���B�  B�ffB�33B���B���B�33B���B�ffB�33B���B���B�ffB���B���B�33B�  B���B���B�ffB�  B���B�ffB�33B�  B���B���B�ffB�33B���Bș�B�33B���Bҙ�B�33B�  Bߙ�B�ffB�  B���BB�ffB�  B���B�ffC��C� CffCL�C	33C�CffC��C�3C�3C�3C��C��C� C� C� C!� C#� C%ffC'ffC)L�C+L�C-L�C/L�C133C333C533C7�C933C;�C=�C?�CA�CC�CE�CG�CI33CK33CM33COL�CQL�CSL�CUL�CWffCYffC[ffC]ffC_ffCa� Cc� Ce� Cg��Ci��Ck�3Cm��Co��Cq��Cs33Cu33CwL�CyL�C{ffC}� C��C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��fC��3C�� C���C���C���C��fC��fC��fC��fC��3C��3C�� C�� C���C�ٚC��fC�� C���C��3C���C��fC��3C���C��fC��3C���C��fC��3C���C��fC�� C��fC�� C�ٚC�� C��fC���C��3C�ٚC�� C��3C���C�� C�ٚC���C��3C���C���C��3C�ٚC�� C��3C��fC���C�� C��fC�ٚC�� CƳ3CǦfCȌ�Cɳ3C��fC���C�� Cͳ3CΦfCϙ�CЌ�Cѳ3C��fC�ٚC�ٚC�ٚC���C�� C�� Cٳ3CڦfCۦfCܦfCݦfCަfCߦfC�fC�fC�fC�fC�3C�3C�3C�� C�� C���C���C�ٚC�fC� C��C��C�C�C�fC�fC��3C��3C�� C���C���C�ٚC�ffC�� D 33D�fD��D�3D33Dl�D��D��D
33Ds3D��D  D33DffD�3D�3D33Dy�D� D�D@ Dl�D� D��DL�D�fD �fD"  D#FfD$�fD%��D'3D(9�D)ffD*��D+��D-L�D.y�D/�fD0��D233D3y�D4��D6fD7FfD8�fD9�fD;fD<@ D=� D>�fD?�3DA,�DBy�DC�fDE  DF33DGy�DH�3DI�3DK9�DLy�DM� DO  DP,�DQs3DR��DT  DUL�DV� DW��DX��DZFfD[y�D\��D]� D_33D`�fDa� Db�3Dd,�DeffDf� Dg� Di&fDjl�Dk�3Dl��Dn@ Do�fDp�3Dr  Ds33Dtl�Du�fDv��DxS3Dyy�Dz�fD{��D}@ D~y�D�3D�|�D��D�� D�VfD�  D���D�C3D�� D�|�D��D���D�Y�D���D���D�<�D��3D�y�D� D���D�Y�D��fD��3D�33D��fD�y�D��D�� D�Y�D��3D���D�FfD��3D�|�D��D���D�Y�D��fD��fD�6fD�ٚD�|�D�  D��fD�\�D�  D���D�33D���D��fD�#3D���D�Y�D���D���D�9�D���D�|�D�#3D��fD�\�D���D��fD�33D���D���D�&fD��fD�\�D��3D��fD�@ D�ٚD�vfD�#3D��3D�` D�  D�� D�C3D��fD�|�D�#3D��fD�\�D�3D���D�33D�ٚD��3D�  D���D�Y�D��fD��fD�6fD��fD�vfD�3D��3D�S3D��fD���D�<�D��fD�|�D�#3D¼�D�\�D���Dę�D�9�D��fD�y�D��D��3D�` D���Dɜ�D�<�D���D�s3D��D�� D�VfD�� DΜ�D�FfD��fDІfD�&fD�ɚD�l�D���DӐ D�33D�� DՆfD�&fD��3D�` D�  Dؠ D�C3D��fDډ�D��D۶fD�` D���Dݙ�D�C3D��3D߀ D�#3D�fD�Y�D���D� D�C3D��fD�|�D�#3D幚D�S3D���D�3D�@ D�� D�|�D��D�� D�Y�D���D� D�9�D��fD�s3D� D��D�ffD�	�D� D�6fD���D� D��D��3D�\�D���D��fD�@ D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�33@�33A��A��A��A.ffAA��AQ��A`  AnffA�  A���A���A�ffA�  A���A�  A���A���A�33A���A�  A�33A���A���A�ffA�ffB��B��B��B  BffB��B  B ffB#��B(  B,��B0  B333B6ffB;��B@ffBC��BF��BK��BPffBS33BV��B[33B`  Bc33Bh  Bl��Bp  Bs33BvffB{��B�33B���B�ffB���B�33B���B���B�  B�ffB�33B���B���B�33B���B�ffB�33B���B���B�ffB���B���B�33B�  B���B���B�ffB�  B���B�ffB�33B�  B���B���B�ffB�33B���Bș�B�33B���Bҙ�B�33B�  Bߙ�B�ffB�  B���BB�ffB�  B���B�ffC��C� CffCL�C	33C�CffC��C�3C�3C�3C��C��C� C� C� C!� C#� C%ffC'ffC)L�C+L�C-L�C/L�C133C333C533C7�C933C;�C=�C?�CA�CC�CE�CG�CI33CK33CM33COL�CQL�CSL�CUL�CWffCYffC[ffC]ffC_ffCa� Cc� Ce� Cg��Ci��Ck�3Cm��Co��Cq��Cs33Cu33CwL�CyL�C{ffC}� C��C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��fC��3C�� C���C���C���C��fC��fC��fC��fC��3C��3C�� C�� C���C�ٚC��fC�� C���C��3C���C��fC��3C���C��fC��3C���C��fC��3C���C��fC�� C��fC�� C�ٚC�� C��fC���C��3C�ٚC�� C��3C���C�� C�ٚC���C��3C���C���C��3C�ٚC�� C��3C��fC���C�� C��fC�ٚC�� CƳ3CǦfCȌ�Cɳ3C��fC���C�� Cͳ3CΦfCϙ�CЌ�Cѳ3C��fC�ٚC�ٚC�ٚC���C�� C�� Cٳ3CڦfCۦfCܦfCݦfCަfCߦfC�fC�fC�fC�fC�3C�3C�3C�� C�� C���C���C�ٚC�fC� C��C��C�C�C�fC�fC��3C��3C�� C���C���C�ٚC�ffC�� D 33D�fD��D�3D33Dl�D��D��D
33Ds3D��D  D33DffD�3D�3D33Dy�D� D�D@ Dl�D� D��DL�D�fD �fD"  D#FfD$�fD%��D'3D(9�D)ffD*��D+��D-L�D.y�D/�fD0��D233D3y�D4��D6fD7FfD8�fD9�fD;fD<@ D=� D>�fD?�3DA,�DBy�DC�fDE  DF33DGy�DH�3DI�3DK9�DLy�DM� DO  DP,�DQs3DR��DT  DUL�DV� DW��DX��DZFfD[y�D\��D]� D_33D`�fDa� Db�3Dd,�DeffDf� Dg� Di&fDjl�Dk�3Dl��Dn@ Do�fDp�3Dr  Ds33Dtl�Du�fDv��DxS3Dyy�Dz�fD{��D}@ D~y�D�3D�|�D��D�� D�VfD�  D���D�C3D�� D�|�D��D���D�Y�D���D���D�<�D��3D�y�D� D���D�Y�D��fD��3D�33D��fD�y�D��D�� D�Y�D��3D���D�FfD��3D�|�D��D���D�Y�D��fD��fD�6fD�ٚD�|�D�  D��fD�\�D�  D���D�33D���D��fD�#3D���D�Y�D���D���D�9�D���D�|�D�#3D��fD�\�D���D��fD�33D���D���D�&fD��fD�\�D��3D��fD�@ D�ٚD�vfD�#3D��3D�` D�  D�� D�C3D��fD�|�D�#3D��fD�\�D�3D���D�33D�ٚD��3D�  D���D�Y�D��fD��fD�6fD��fD�vfD�3D��3D�S3D��fD���D�<�D��fD�|�D�#3D¼�D�\�D���Dę�D�9�D��fD�y�D��D��3D�` D���Dɜ�D�<�D���D�s3D��D�� D�VfD�� DΜ�D�FfD��fDІfD�&fD�ɚD�l�D���DӐ D�33D�� DՆfD�&fD��3D�` D�  Dؠ D�C3D��fDډ�D��D۶fD�` D���Dݙ�D�C3D��3D߀ D�#3D�fD�Y�D���D� D�C3D��fD�|�D�#3D幚D�S3D���D�3D�@ D�� D�|�D��D�� D�Y�D���D� D�9�D��fD�s3D� D��D�ffD�	�D� D�6fD���D� D��D��3D�\�D���D��fD�@ D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�n�@� �@��@�;d@�7L@��@�V@���@{S�@rM�@iG�@`��@^5?@Z�\@V@[C�@\��@UO�@IX@FV@B�H@>�@>ff@>�+@@�`@A�@?�;@@�`@@bN@:J@Bn�@C�m@E?}@IG�@H�9@KC�@M/@P  @O�@O�@Pb@P1'@P�`@Q��@P�u@PbN@Q&�@QX@P��@P �@N�R@M�T@PQ�@Nȴ@N��@N��@N��@N�@N��@N�R@Nȴ@M?}@L��@MO�@J~�@I��@Ix�@I7L@IG�@I�@H�@H �@G��@G|�@G\)@Fv�@D(�@D1@CdZ@B~�@A��@A7L@A��@B�@@�@@�u@@Q�@?��@?|�@?K�@?\)@?;d@>�y@>{@=�T@=?}@<j@<�@;t�@;t�@;dZ@;S�@;S�@;33@:��@:M�@:-@9�7@9�7@8��@8r�@81'@8 �@8b@8  @8r�@8��@8�@8Q�@8 �@7�;@7�@6��@7�P@7;d@6��@5��@5O�@4��@3�F@4�@2��@1��@1X@1X@0�@/��@/��@/l�@.�@.ȴ@.ȴ@.V@-�-@-?}@,9X@+t�@)�#@( �@(  @'��@&�@%��@&��@'l�@&@#�
@#33@"�!@"-@!hs@ �`@!hs@ Ĝ@ �@�@�@�y@/@�@I�@M�@�@hs@A�@�;@�@��@@�!@�@l�@�@K�@�w@t�@|�@V@I�@O�@�@�@�@5?@��@A�@�9@r�@ȴ@C�@��@ Ĝ@  �?��?��h?�?���?�M�?ꟾ?��`?�|�?��?�5??ݲ-?ܬ?�Q�?׍P?ؓu?�Q�?ش9?��?�z�?�1'?�1'?��T?��?��?��y?��?��?�;d?��?�b?ٙ�?ش9?�l�?ҏ\?�J?Ѓ?��?�p�?�~�?�X?�ff?ě�?���?���?��D?�ƨ?�X?���?��j?�o?��`?�5??���?�C�?���?�ȴ?���?���?���?��?���?�S�?���?��?��?���?��h?�/?���?�1?�1?��?��H?�=q?��9?��?�ȴ?��
?��!?�hs?��`?~��?� �?;d?�Ĝ?�S�?���?��?�z�?�`B?�9X?��?��?��\?|(�?wK�?s�F?o��?j=q?d�?b��?`A�?_�w?^v�?`A�?^v�?]�?[�m?Tz�?S33?V�+?Q�?Q&�?Rn�?P��?N��?L��?G�?D��?>��?<(�?<(�?:^5?4��?,I�?%��?   ?(�?�u?bN?
~�?��>�^5>��>�x�>�5?>�>��`>�=q>\>��>�?}>�{>�V>���>�G�>�(�>�>��;>���>�  >w��>n��>aG�>V>J��>@�>1&�>z�==���=Ƨ�=�1=�C�=ix�=L��=H�9=�w=#�
=��<�`B<T��:�o�49X��C���`B���C��<j�ixս�%��O߽�1��vɽȴ9������G���h���#�o���
=q�V�hs�n��������-�#�
�&�y�,1�/��333�6E��<j�?|�A�7�C���H�9�M��S�ϾW
=�Y��Z��]/�_;d�bMӾgl��ixվl�D�p�׾s�F�vȴ�y�#�z�H��  ��J����������$ݾ�1'��7L��=q��ƨ��O߾����bN��hs��녾��Ͼ������+��
=��b������"Ѿ�/���R��A���G�������MӾ�S���Z���T��l���r���xվ���1��V��{�� ž�-���F����ȴ���پ�X���H���m��p���|��  ���7��o�ě�����Ƨ�ȴ9��=q��I����;�����;���`����z�և+�ؓu�����(���5?�޸R��A���MӾ�Z��`B��ff������D�������33��F��9X������ȴ���پ�^5��j�����ۿ   � ��%�G��Mӿ��o��
��/�`B�ff�l������r��	7L�	�^�
=q�
�������
�����I���ͿV�V��ͿO߿V���\)� ſ�׿�׿&�-�����!��!�t��z��E��
=��P�b��u�����������X��#�^5��H�"ѿ��dZ�(��j�푿/�/�p���vɿ�ۿ|��w� A�� �1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@���@�n�@� �@��@�;d@�7L@��@�V@���@{S�@rM�@iG�@`��@^5?@Z�\@V@[C�@\��@UO�@IX@FV@B�H@>�@>ff@>�+@@�`@A�@?�;@@�`@@bN@:J@Bn�@C�m@E?}@IG�@H�9@KC�@M/@P  @O�@O�@Pb@P1'@P�`@Q��@P�u@PbN@Q&�@QX@P��@P �@N�R@M�T@PQ�@Nȴ@N��@N��@N��@N�@N��@N�R@Nȴ@M?}@L��@MO�@J~�@I��@Ix�@I7L@IG�@I�@H�@H �@G��@G|�@G\)@Fv�@D(�@D1@CdZ@B~�@A��@A7L@A��@B�@@�@@�u@@Q�@?��@?|�@?K�@?\)@?;d@>�y@>{@=�T@=?}@<j@<�@;t�@;t�@;dZ@;S�@;S�@;33@:��@:M�@:-@9�7@9�7@8��@8r�@81'@8 �@8b@8  @8r�@8��@8�@8Q�@8 �@7�;@7�@6��@7�P@7;d@6��@5��@5O�@4��@3�F@4�@2��@1��@1X@1X@0�@/��@/��@/l�@.�@.ȴ@.ȴ@.V@-�-@-?}@,9X@+t�@)�#@( �@(  @'��@&�@%��@&��@'l�@&@#�
@#33@"�!@"-@!hs@ �`@!hs@ Ĝ@ �@�@�@�y@/@�@I�@M�@�@hs@A�@�;@�@��@@�!@�@l�@�@K�@�w@t�@|�@V@I�@O�@�@�@�@5?@��@A�@�9@r�@ȴ@C�@��@ Ĝ@  �?��?��h?�?���?�M�?ꟾ?��`?�|�?��?�5??ݲ-?ܬ?�Q�?׍P?ؓu?�Q�?ش9?��?�z�?�1'?�1'?��T?��?��?��y?��?��?�;d?��?�b?ٙ�?ش9?�l�?ҏ\?�J?Ѓ?��?�p�?�~�?�X?�ff?ě�?���?���?��D?�ƨ?�X?���?��j?�o?��`?�5??���?�C�?���?�ȴ?���?���?���?��?���?�S�?���?��?��?���?��h?�/?���?�1?�1?��?��H?�=q?��9?��?�ȴ?��
?��!?�hs?��`?~��?� �?;d?�Ĝ?�S�?���?��?�z�?�`B?�9X?��?��?��\?|(�?wK�?s�F?o��?j=q?d�?b��?`A�?_�w?^v�?`A�?^v�?]�?[�m?Tz�?S33?V�+?Q�?Q&�?Rn�?P��?N��?L��?G�?D��?>��?<(�?<(�?:^5?4��?,I�?%��?   ?(�?�u?bN?
~�?��>�^5>��>�x�>�5?>�>��`>�=q>\>��>�?}>�{>�V>���>�G�>�(�>�>��;>���>�  >w��>n��>aG�>V>J��>@�>1&�>z�==���=Ƨ�=�1=�C�=ix�=L��=H�9=�w=#�
=��<�`B<T��:�o�49X��C���`B���C��<j�ixս�%��O߽�1��vɽȴ9������G���h���#�o���
=q�V�hs�n��������-�#�
�&�y�,1�/��333�6E��<j�?|�A�7�C���H�9�M��S�ϾW
=�Y��Z��]/�_;d�bMӾgl��ixվl�D�p�׾s�F�vȴ�y�#�z�H��  ��J����������$ݾ�1'��7L��=q��ƨ��O߾����bN��hs��녾��Ͼ������+��
=��b������"Ѿ�/���R��A���G�������MӾ�S���Z���T��l���r���xվ���1��V��{�� ž�-���F����ȴ���پ�X���H���m��p���|��  ���7��o�ě�����Ƨ�ȴ9��=q��I����;�����;���`����z�և+�ؓu�����(���5?�޸R��A���MӾ�Z��`B��ff������D�������33��F��9X������ȴ���پ�^5��j�����ۿ   � ��%�G��Mӿ��o��
��/�`B�ff�l������r��	7L�	�^�
=q�
�������
�����I���ͿV�V��ͿO߿V���\)� ſ�׿�׿&�-�����!��!�t��z��E��
=��P�b��u�����������X��#�^5��H�"ѿ��dZ�(��j�푿/�/�p���vɿ�ۿ|��w� A�� �1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��BB1BiyBĜB	;dB	iyB
�B
9XB
 �B
uB
�B
�B
9XB
iyB
�VB
��B
��B
��B
ŢB
�#B
�BB
�BDB�B%�B@�BH�BR�BW
Bp�B}�B�+B�uB��B��B�B�^BȴB��B��B��B��B�
B�#B�ZB�`B�TB�`B�mB�yB�mB�mB�B�B�B�B�B�B�B�B�B�yB�B�B�B�B�B�B�fB�B�B�B�B�B�B�sB�yB�mB�sB�fB�sB�`B�ZB�yB�yB�TB�B�mB�fB�mB�mB�sB�sB�yB�sB�yB�sB�yB�sB�sB�yB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�sB�mB�fB�mB�fB�fB�mB�yB�sB�`B�`B�ZB�TB�TB�HB�ZB�`B�`B�TB�NB�TB�NB�NB�BB�HB�;B�5B�;B�#B�#B�B�B�B�B�B��B�B�
B�5B�5B��B��B��B��B�
B�B�B�)B�5B�5B�5B�B�B�
B�B�B��B��B��B��B��BɺB��BB��BB��B�}B�wB��B��B��B��B��BƨBBȴBĜBǮBBÖB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBȴBɺBɺBǮBȴBǮBǮBƨBƨBƨBŢBÖBB��B�wB�qB�qB�}B�jB�dB�jB�jB�qB�qB�qB�jB�^B�jB�dB�^B�^B�^B�^B�XB�^B�dB�dB�jB�wB��B��BÖBÖBĜBƨBǮBȴBƨBŢBŢBĜBÖB��BBÖBĜBŢBƨBƨBǮBǮBŢBĜBǮBƨBƨBǮBǮBǮBǮBŢBĜBĜBBÖBB��B�}B�wB�jB�dB�XB�LB�LB�?B�9B�-B�-B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B��B�@B�Bo1B�ZB	A&B	o?B
aB
?%B
&�B
?B
QB
eB
?"B
oGB
�#B
��B
��B
��B
�pB
��B
�B
��BB!zB+�BFXBN�BX�B\�Bv{B��B�B�NB��B��B��B�8BΏBОBӯBչB��B��B��B�7B�=B�0B�=B�KB�UB�KB�LB�mB�hB�hB�tB�~B��B��B��B�|B�WB�bB�tB�mB�nB�oB�gB�CB�zB�nB�iB�bB�iB�bB�NB�XB�MB�NB�CB�NB�<B�6B�XB�VB�2B�\B�IB�DB�IB�KB�QB�RB�VB�RB�TB�OB�TB�PB�OB�VB�TB�[B�aB�cB�cB�gB�hB�jB�hB�hB�nB�nB�nB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�vB�vB�tB�vB�oB�oB�mB�nB�oB�tB�oB�oB�hB�ZB�dB�NB�KB�BB�KB�BB�BB�IB�TB�NB�>B�>B�6B�0B�0B�$B�7B�<B�>B�0B�'B�0B�'B�(B�B�#B�B�B�B��B��B��B��B��B��B��B��B��B��B�B�B��BМB��B��B��B��B��B�B�B�B�B��B��B��B��B��B��B��B��B��BҧBϖB�ZB�jB�]B�gB�cB�VB�OB�]B�]B�cB�bB�[B̀B�gBΎB�tB͆B�iB�oB�]BҧB��B��BչBֽB��B��B��BֽBնBӬBնBչBӭBշBѠBҧBЙBϑBϔBϓBΏBϒBϒBͅB΋B̓B͉B̀B�~B�B�{B�mB�hB�[B�PB�IB�KB�XB�AB�;B�@B�AB�IB�JB�HB�DB�7B�AB�<B�4B�2B�5B�2B�/B�5B�:B�:B�BB�OB�[B�_B�qB�oB�vB�B͆BΊB́B�yB�zB�tB�mB�aB�gB�oB�tB�|B�}B́B͇B̈́B�xB�tB̈́B�~B́B̈́B͆B͇BͅB�xB�uB�tB�hB�mB�gB�^B�RB�LB�BB�;B�/B�#B�$B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.005703 (+/- 0.01)                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219022022012717040120220127170401  IF  ARFMCODA035h                                                                20200828144323                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144407  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200828144407  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121902  IP  PSAL            @�33D��3G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170401  IP  PSAL            @�33D��3G�O�                