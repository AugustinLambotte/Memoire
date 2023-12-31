CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  I   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:29Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8    FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8(   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    88   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8H   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8P   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9    	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9,   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    90   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     94   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9T   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9t   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	$  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  C�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	$  FH   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  Ol   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	$  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  Z�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  d    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  fL   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  op   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  q�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  z�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  �P   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  �t   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	$  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �@   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �D   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �H   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �L   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �P   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �Argo profile    3.1 1.2 19500101000000  20200829092229  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               8A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��v�8�1   @��v�8�@RkA�k�*8��r�8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A�ffA�33A���A���A�ffA�ffA�33A�  A���A�ffA�  A���A�33BffB��BffB��B  B��B  B   B$ffB(  B,ffB0  B4  B8  B<  B@ffBDffBH  BJ��BO33BS33BX  B[33B`  BdffBhffBl  Bp��Bs33Bx  B|ffB�  B�  B�33B�33B�  B�33B�  B���B�  B�  B�  B�33B���B���B���B�33B�  B���B�ffB�33B���B���B�ffB�33B�  B�  B���B���B���B�ffB�ffB�33B�33B���B���Bƙ�Bș�B�33B�  Bҙ�B�33B���B�ffB�33B�  B���BB�ffB�ffB�33B�  CffC�fC�fC� C	  C  C�C�C33C33C33CL�C��C��C�3C��C!��C#�fC%�fC'ffC)  C+33C-33C/L�C133C333C5ffC7ffC9L�C;L�C=L�C?L�CA33CC33CE33CG33CI33CK33CM33CO33CQ33CS33CU33CW33CY33C[L�C]L�C_L�CaL�CcffCeffCg� Ci� Ck� Cm��CoL�CqL�CsffCu� Cw��CyL�C{ffC}� C33C��fC�� C���C��fC�� C���C��3C���C��fC���C��3C���C��3C���C�� C��fC���C�� C��fC���C�� C��fC�ٚC���C�� C��3C��fC��fC���C���C���C���C���C�� C�� C�� C�� C�� C���C���C���C��fC��3C�� C�� C��3C��3C��3C��3C��3C��fC��fC��fC��fC��3C��3C�� C�� C�� C���C���C��fC��3C��3C�� C���C�ٚCæfCĳ3C�� C�� C���CȦfCɳ3C�� C˙�C̦fCͳ3C�� C�ٚCг3Cь�CҦfCӳ3Cԙ�CզfC�� Cי�C�� C���Cڳ3Cۙ�Cܳ3C���C޳3Cߙ�C�3C���C�3C㙚C�� C�ٚC�� C�fC虚C� C�3C�ٚC���C��3C�fC��C�� C��fC�ٚC���C�� C��3C��fC���C���C��3C��fC��D L�D��D��D�DL�D��D��D	�D
S3D��D� D��D33Dy�D�3D��D@ Ds3D��D��D,�Ds3D�3D��D@ D��D ��D!��D#9�D$�fD%��D'fD(9�D)��D*�fD,fD-L�D.s3D/�3D1  D2,�D3` D4�3D6�D7FfD8� D9� D;  D<@ D=� D>��D?��DA@ DB�fDC�fDE�DF9�DG` DH��DI��DK33DLffDM�fDN�fDP,�DQl�DR�3DS��DUFfDVy�DW��DY  DZ@ D[y�D\�3D]�3D_9�D`� Da� Db��Dd9�De� Df�3Dg�fDi9�Dj��Dk� Dl��Dn33Dos3Dp�3Dq��Ds@ Dt�fDu��Dv��Dx9�Dy�fDz� D{��D}33D~s3D�3D�|�D�  D��fD�` D���D��3D�@ D�� D�|�D��D���D�` D�3D���D�@ D��3D�y�D�  D��fD�` D���D��3D�33D�� D�p D��D�� D�P D��3D��fD�9�D���D�vfD��D��3D�\�D���D���D�<�D�� D��3D��D��3D�Y�D�3D���D�9�D��3D��3D�  D��3D�VfD���D��3D�9�D�� D��fD�  D���D�VfD��3D�� D�<�D���D���D�,�D���D�` D��3D��3D�@ D�� D�y�D�  D�ɚD�ffD�fD��fD�FfD��fD���D�  D��fD�` D���D��3D�<�D�ٚD�vfD�#3D��3D�c3D�3D���D�9�D���D��3D�&fD���D�S3D���D��fD�@ D�ٚD�vfD�3D³3D�S3D��3Dē3D�6fD��fD�|�D�#3DǼ�D�\�D���Də�D�9�D���Dˀ D�#3D��fD�` D��fDΜ�D�FfD�� D�|�D��DѶfD�c3D�3DӖfD�6fD�ٚD�|�D�  DֶfD�\�D�3Dؠ D�9�D��fD�s3D� D۰ D�P D��3DݖfD�6fD���D߃3D�  D��D�\�D���D��D�@ D�� D�3D�&fD��D�S3D���D�3D�C3D�� D�3D�  D��3D�VfD���D��D�@ D��fD� D��D�3D�\�D�	�D�fD�C3D�� D�|�D��D�� D�ffD�  D���D�6fD��3D�p D�3D��fD�Y�D���D���D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A�ffA�33A���A���A�ffA�ffA�33A�  A���A�ffA�  A���A�33BffB��BffB��B  B��B  B   B$ffB(  B,ffB0  B4  B8  B<  B@ffBDffBH  BJ��BO33BS33BX  B[33B`  BdffBhffBl  Bp��Bs33Bx  B|ffB�  B�  B�33B�33B�  B�33B�  B���B�  B�  B�  B�33B���B���B���B�33B�  B���B�ffB�33B���B���B�ffB�33B�  B�  B���B���B���B�ffB�ffB�33B�33B���B���Bƙ�Bș�B�33B�  Bҙ�B�33B���B�ffB�33B�  B���BB�ffB�ffB�33B�  CffC�fC�fC� C	  C  C�C�C33C33C33CL�C��C��C�3C��C!��C#�fC%�fC'ffC)  C+33C-33C/L�C133C333C5ffC7ffC9L�C;L�C=L�C?L�CA33CC33CE33CG33CI33CK33CM33CO33CQ33CS33CU33CW33CY33C[L�C]L�C_L�CaL�CcffCeffCg� Ci� Ck� Cm��CoL�CqL�CsffCu� Cw��CyL�C{ffC}� C33C��fC�� C���C��fC�� C���C��3C���C��fC���C��3C���C��3C���C�� C��fC���C�� C��fC���C�� C��fC�ٚC���C�� C��3C��fC��fC���C���C���C���C���C�� C�� C�� C�� C�� C���C���C���C��fC��3C�� C�� C��3C��3C��3C��3C��3C��fC��fC��fC��fC��3C��3C�� C�� C�� C���C���C��fC��3C��3C�� C���C�ٚCæfCĳ3C�� C�� C���CȦfCɳ3C�� C˙�C̦fCͳ3C�� C�ٚCг3Cь�CҦfCӳ3Cԙ�CզfC�� Cי�C�� C���Cڳ3Cۙ�Cܳ3C���C޳3Cߙ�C�3C���C�3C㙚C�� C�ٚC�� C�fC虚C� C�3C�ٚC���C��3C�fC��C�� C��fC�ٚC���C�� C��3C��fC���C���C��3C��fC��D L�D��D��D�DL�D��D��D	�D
S3D��D� D��D33Dy�D�3D��D@ Ds3D��D��D,�Ds3D�3D��D@ D��D ��D!��D#9�D$�fD%��D'fD(9�D)��D*�fD,fD-L�D.s3D/�3D1  D2,�D3` D4�3D6�D7FfD8� D9� D;  D<@ D=� D>��D?��DA@ DB�fDC�fDE�DF9�DG` DH��DI��DK33DLffDM�fDN�fDP,�DQl�DR�3DS��DUFfDVy�DW��DY  DZ@ D[y�D\�3D]�3D_9�D`� Da� Db��Dd9�De� Df�3Dg�fDi9�Dj��Dk� Dl��Dn33Dos3Dp�3Dq��Ds@ Dt�fDu��Dv��Dx9�Dy�fDz� D{��D}33D~s3D�3D�|�D�  D��fD�` D���D��3D�@ D�� D�|�D��D���D�` D�3D���D�@ D��3D�y�D�  D��fD�` D���D��3D�33D�� D�p D��D�� D�P D��3D��fD�9�D���D�vfD��D��3D�\�D���D���D�<�D�� D��3D��D��3D�Y�D�3D���D�9�D��3D��3D�  D��3D�VfD���D��3D�9�D�� D��fD�  D���D�VfD��3D�� D�<�D���D���D�,�D���D�` D��3D��3D�@ D�� D�y�D�  D�ɚD�ffD�fD��fD�FfD��fD���D�  D��fD�` D���D��3D�<�D�ٚD�vfD�#3D��3D�c3D�3D���D�9�D���D��3D�&fD���D�S3D���D��fD�@ D�ٚD�vfD�3D³3D�S3D��3Dē3D�6fD��fD�|�D�#3DǼ�D�\�D���Də�D�9�D���Dˀ D�#3D��fD�` D��fDΜ�D�FfD�� D�|�D��DѶfD�c3D�3DӖfD�6fD�ٚD�|�D�  DֶfD�\�D�3Dؠ D�9�D��fD�s3D� D۰ D�P D��3DݖfD�6fD���D߃3D�  D��D�\�D���D��D�@ D�� D�3D�&fD��D�S3D���D�3D�C3D�� D�3D�  D��3D�VfD���D��D�@ D��fD� D��D�3D�\�D�	�D�fD�C3D�� D�|�D��D�� D�ffD�  D���D�6fD��3D�p D�3D��fD�Y�D���D���D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�������ٿ���ٿ��ٿ��ٿ����������ٿ����ٿ�����P��P��P��P��l���l���Kǿ�Kǿ�Kǿ�+��
=��
=��ȴ��+��+��$ݿ��T��`B���/���
��o��n����;�ۅ��?}��G���7L�������-��������I�^�$�� ���%�#�
=���>q��>��?>v�?Kƨ?i��?��?��m?�l�?�^5?z�H?h�9?���?�b?�ƨ?�%?�n�?���?�j?�G�?�+?У�?��T?׮?�?�"�?ߝ�?�V?�~�@��@	��@�?���?�G�?��?�hs?�S�?��?��?��m?�v�?��?�C�?�~�?�O�?��?��?��R?���?�;d?���?��H?�+?���?��?��?�  ?��H?�ȴ?��?���?���?�Q�?�?ӶF?У�?Η�?�j?��?�^5?�~�?ȴ9?Ƨ�?�bN?��?��H?�b?�E�?�33?���?���?�5??�V?��y?�o?�-?�bN?�r�?���?�`B?�S�?���?�|�?��?�I�?���?�+?�n�?�A�?�  ?|�?}�-?z^5?y�?wK�?w
=?t9X?rn�?p�`?o\)?nV?kC�?h1'?e�T?d�?b��?`�?\�?Y��?U?}?PbN?MO�?I��?Hr�?H��?H��?I7L?G�?D�/?;��?8b?49X?4��?8��?=/?>5??>5??=p�?;�m?9�#?8��?6?+?"�\?�?dZ?�m? A�?-�h?-��?,1?&ff?#��?#�
?#�
?$�?&$�?&ff?'�?'�?'+?&�y?'l�?'l�?-O�?1hs?0��?-��?,��?-V?)��?&��?"�\?!�7?�w?�?b??n�?��?
��?ff?Z? �>�^5>��#>�^5>���>�ȴ>� �>�`B>�~�>���>�ff>�;d>��>�^5>��>�o>��>��>�^5>�K�>�9X>�`B>�S�>�(�>���>p��>gl�>e`B>dZ>dZ>l�D>u>��>���>�o>�J>�  >s�F>Xb>,1>�w>
=q=�
==ȴ9=�{=� �=� �=���=���=���=�^5=���=H�9=ě�>V>	7L>t�>V==���=���=���=Ƨ�=�9X=�{=�Q�=��
=�Q�=�{=� �=���=��=���=�C�=y�#=P�`=@�=��<���<��
<�o;��
�t��D�����㼼j�����t��<j�Y���+��C���hs�������
���罩�署-��^5��^5��vɽȴ9����
=��S���xս�������	7L�I��hs��P��-� Ĝ�$�/�(�þ.{�1&�5?}�8Q�;dZ�=p��@��D���F��I�^�N��P�`�T���Z��]/�_;d�e`B�gl��j~��k��n���p�׾s�F�vȴ�x���y�#�{�m�~�۾�J��������1'���^��ƨ��O߾�������`��n����Ͼ�zᾔ�����+���P���u�����������㾜���/��/���R���w��A���Ĝ���徤Z���/���T��l����þ��羫���1��{����� ž������!��33��9X��?}��E����پ��H��푾�󶾾�۾�  ��J�Õ������$ݾǮ��1'������I�������;��hs��t��������׍P��b�����(���5?��;d��A����������Z���T���y���þ�~�����1��{�����׾����33��F���j��E���ȴ���پ�X���H���m��푾�p���vɾ�|� A�� Ĝ�%�J���S�����/�`B��˿ff��y���r���9�	�^�
~��
������ƨ�I���ͿO߿����V���\)�\)�bN��`����녿-��!�33��Ͽ�j��E��E��ȴ�Kǿ�P�b��u�X��#����"ѿdZ���(��j�/��-��5?��R��R��R��ۿvɿj��m�푿5?�5?��-��-���R�!%�"��#S��#�
�#o�#�
�%��%�˿%�˿$���$���&ff�&ff�&��%`B�%�˿&$ݿ&$ݿ&ff�&��&�y�'��'(1'�(�ÿ)7L�)xտ)�^�*=q�*~��*���+111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   �����ٿ���ٿ��ٿ��ٿ����������ٿ����ٿ�����P��P��P��P��l���l���Kǿ�Kǿ�Kǿ�+��
=��
=��ȴ��+��+��$ݿ��T��`B���/���
��o��n����;�ۅ��?}��G���7L�������-��������I�^�$�� ���%�#�
=���>q��>��?>v�?Kƨ?i��?��?��mG�O�G�O�?z�H?h�9?���?�b?�ƨ?�%?�n�?���?�j?�G�?�+?У�?��T?׮?�?�"�?ߝ�?�V?�~�@��@	��@�?���?�G�?��?�hs?�S�?��?��?��m?�v�?��?�C�?�~�?�O�?��?��?��R?���?�;d?���?��H?�+?���?��?��?�  ?��H?�ȴ?��?���?���?�Q�?�?ӶF?У�?Η�?�j?��?�^5?�~�?ȴ9?Ƨ�?�bN?��?��H?�b?�E�?�33?���?���?�5??�V?��y?�o?�-?�bN?�r�?���?�`B?�S�?���?�|�?��?�I�?���?�+?�n�?�A�?�  ?|�?}�-?z^5?y�?wK�?w
=?t9X?rn�?p�`?o\)?nV?kC�?h1'?e�T?d�?b��?`�?\�?Y��?U?}?PbN?MO�?I��?Hr�?H��?H��?I7L?G�?D�/?;��?8b?49X?4��?8��?=/?>5??>5??=p�?;�m?9�#?8��?6?+?"�\?�?dZ?�m? A�?-�h?-��?,1?&ff?#��?#�
?#�
?$�?&$�?&ff?'�?'�?'+?&�y?'l�?'l�?-O�?1hs?0��?-��?,��?-V?)��?&��?"�\?!�7?�w?�?b??n�?��?
��?ff?Z? �>�^5>��#>�^5>���>�ȴ>� �>�`B>�~�>���>�ff>�;d>��>�^5>��>�o>��>��>�^5>�K�>�9X>�`B>�S�>�(�>���>p��>gl�>e`B>dZ>dZ>l�D>u>��>���>�o>�J>�  >s�F>Xb>,1>�w>
=q=�
==ȴ9=�{=� �=� �=���=���=���=�^5=���=H�9=ě�>V>	7L>t�>V==���=���=���=Ƨ�=�9X=�{=�Q�=��
=�Q�=�{=� �=���=��=���=�C�=y�#=P�`=@�=��<���<��
<�o;��
�t��D�����㼼j�����t��<j�Y���+��C���hs�������
���罩�署-��^5��^5��vɽȴ9����
=��S���xս�������	7L�I��hs��P��-� Ĝ�$�/�(�þ.{�1&�5?}�8Q�;dZ�=p��@��D���F��I�^�N��P�`�T���Z��]/�_;d�e`B�gl��j~��k��n���p�׾s�F�vȴ�x���y�#�{�m�~�۾�J��������1'���^��ƨ��O߾�������`��n����Ͼ�zᾔ�����+���P���u�����������㾜���/��/���R���w��A���Ĝ���徤Z���/���T��l����þ��羫���1��{����� ž������!��33��9X��?}��E����پ��H��푾�󶾾�۾�  ��J�Õ������$ݾǮ��1'������I�������;��hs��t��������׍P��b�����(���5?��;d��A����������Z���T���y���þ�~�����1��{�����׾����33��F���j��E���ȴ���پ�X���H���m��푾�p���vɾ�|� A�� Ĝ�%�J���S�����/�`B��˿ff��y���r���9�	�^�
~��
������ƨ�I���ͿO߿����V���\)�\)�bN��`����녿-��!�33��Ͽ�j��E��E��ȴ�Kǿ�P�b��u�X��#����"ѿdZ���(��j�/��-��5?��R��R��R��ۿvɿj��m�푿5?�5?��-��-���R�!%�"��#S��#�
�#o�#�
�%��%�˿%�˿$���$���&ff�&ff�&��%`B�%�˿&$ݿ&$ݿ&ff�&��&�y�'��'(1'�(�ÿ)7L�)xտ)�^�*=q�*~��*���+111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�ZB�ZB�fB�B�B�B�B�B��B��B��B��B��B��B��B��B��B��BBBBBB%B%B1BDBVB\B{B�B�B�B"�B+B1'B:^BB�BN�By�B��B�qB�`B��B�BG�Bw�B��B��B		7B	5?B	q�B	�oB	�}B	�5B
%�B
.B
YB
l�B
�B
�`B
�{B
�}B
�-B
B
��B
��B
�B
�B
��B
��BBbB�B�B �B �B)�B&�B7LBH�B]/Bl�Bo�Bm�BhsBm�Bu�Bz�B�B�bB��B��B��B��B��B�B�3B�?B�XB�dB�wB��B�}B��B��BBÖBĜBŢBĜBBB��BB��B��B��B��B��B��BÖBĜBÖBBB�wB�wB�}B�wB�dB�jB�jB�qB�}B�dB�dB�XB�FB�LB�FB�?B�?B�?B�?B�?B�?B�FB�9B�3B�3B�-B�3B�3B�3B�3B�3B�FB�?B�?B�?B�?B�?B�9B�9B�9B�9B�9B�9B�9B�9B�9B�3B�3B�3B�3B�3B�?B�?B�?B�3B�-B�'B�-B�-B�FB�LB�FB�RB�LB�LB�RB�LB�LB�jB�3B�-B�'B�'B�FB�FB�RB�RB�^B�LB�LB�RB�RB�XB�^B�dB�dB�dB�jB�jB�qB�wB��B��B��B�}B��B�wB�wB�qB�qB�jB�dB�^B�XB�LB�RB�RB�FB�?B�?B�9B�9B�9B�9B�3B�3B�'B�3B�-B�'B�3B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�ZB�ZB�fB�B�B�B�B�B��B��B��B��B��B��B��B��B��B��BBBBBB%B%B1BDBVB\B{B�B�B�B"�B+B1'B:^BB�BN�By�B��B�qB�`B��B�BG�Bw�B��B��B		7B	5?B	q�B	�oB	�}B	�5B
%�B
.B
YB
l�B
�G�O�G�O�B
�}B
�-B
B
��B
��B
�B
�B
��B
��BBbB�B�B �B �B)�B&�B7LBH�B]/Bl�Bo�Bm�BhsBm�Bu�Bz�B�B�bB��B��B��B��B��B�B�3B�?B�XB�dB�wB��B�}B��B��BBÖBĜBŢBĜBBB��BB��B��B��B��B��B��BÖBĜBÖBBB�wB�wB�}B�wB�dB�jB�jB�qB�}B�dB�dB�XB�FB�LB�FB�?B�?B�?B�?B�?B�?B�FB�9B�3B�3B�-B�3B�3B�3B�3B�3B�FB�?B�?B�?B�?B�?B�9B�9B�9B�9B�9B�9B�9B�9B�9B�3B�3B�3B�3B�3B�?B�?B�?B�3B�-B�'B�-B�-B�FB�LB�FB�RB�LB�LB�RB�LB�LB�jB�3B�-B�'B�'B�FB�FB�RB�RB�^B�LB�LB�RB�RB�XB�^B�dB�dB�dB�jB�jB�qB�wB��B��B��B�}B��B�wB�wB�qB�qB�jB�dB�^B�XB�LB�RB�RB�FB�?B�?B�9B�9B�9B�9B�3B�3B�'B�3B�-B�'B�3B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092229                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092322  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092322  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A�ffD���G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172540  IP  PSAL            A�ffD���G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A�ffD���G�O�                