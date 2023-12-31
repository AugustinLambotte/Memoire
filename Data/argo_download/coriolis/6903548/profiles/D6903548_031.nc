CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  N   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:50Z creation; 2021-06-07T15:42:39Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        	8  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	8  F`   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	8  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  [    TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  dX   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  f�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  o�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  r0   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  {h   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	8  �x   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �\   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �l   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �p   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20190620115650  20210607154239  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @؇���>�1   @؇���>�@TW�S�x�@*�.HU��8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A0  AA��AQ��Aa��Ap  A���A���A���A�  A�33A�  A�33A�  A���A�  A�  A���A�33A�ffA�33A�ffB ffBffBffB33BffBffB��B��B��B$ffB(  B+��B/��B3��B7��B;33B?33BC33BG33BK33BO��BS��BX  B\  B`��Bc��Bh  BlffBo��Bt  BxffB{��B~��B���B�  B���B�ffB���B�33B���B���B�33B���B�ffB�  B���B���B�ffB�  B���B�ffB�33B�  B���B���B���B�ffB�33B�  B���B���B���B�ffB�33B�33B�  B�  B���BǙ�Bʙ�B�ffB�33B�33B�33B�33B�33B�33B�33B�33B�33B�33B�ffB�ffC33CL�CL�CL�C	ffCffCffC� C��C�3C��C�fC� C  C�C33C!L�C#ffC%ffC'��C)��C+�3C-��C/ffC1  C3�C533C7L�C9� C;��C=�3C?��CAffCC  CE�CG33CIL�CKffCM��CO�3CQ��CSffCU  CW�CY�C[L�C]ffC_� Ca��Cc�3Ce��CgffCi  Ck�Cm33CoL�CqffCs� Cu��Cw�3Cy�fC{ffC}  C�C��fC��fC��3C���C�ٚC��fC��3C�� C���C���C��fC��3C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��fC��3C�� C���C��3C��3C��3C��3C��3C�� C�� C���C���C��fC��3C�� C���C��fC��3C���C��fC�� C�ٚC��3C���C�� C�ٚC�� C��fC���C��3C�ٚC�� C��fC���C�� C��fC���C��3C��fC���C�� C��3C�ٚC���C�� Cų3CƦfCǙ�CȌ�Cɳ3C��fC��fC�ٚC�ٚC�ٚC���C���C���C���C�ٚC�ٚC�ٚC��fC��3Cس3Cـ Cڀ Cی�Cܙ�CݦfC޳3C�� C�ٚC��fC�3C� C䙚C�fC�� C�ٚC�3C��C�fC�� C�ٚC��fC�� C��C�fC�� C��fC�� C��C��fC�� C�ٚC�� C���C�s3C��3D 9�Ds3D�3D��D,�D�fD��D�3D
L�D� D�3D�D@ Ds3D��D  D33D��D��D� D&fDy�D� D�fD&fDffD ��D!��D#FfD$�fD%��D&��D(S3D)��D*�fD,�D-L�D.�3D/��D0ٚD29�D3��D4ٚD5��D7  D8ffD9�fD;  D<FfD=y�D>�3D?�3DA9�DB�fDC� DD�3DFFfDG� DH��DI��DK@ DL��DM��DN�fDP33DQs3DR�3DS�3DU,�DV�fDW� DX�3DZ33D[l�D\� D]��D_33D`s3Da��DcfDdFfDey�Df��Dg�3Di33Djs3Dk��Dl��Dn9�Do� Dp�fDq�3Ds&fDts3Du�fDv��Dx9�Dys3Dz��D{�fD}9�D~� D�fD�y�D��D���D�c3D�3D��3D�FfD��3D��3D�&fD��fD�ffD���D��3D�9�D��3D��3D��D��3D�\�D�	�D��fD�FfD��fD��fD�&fD��fD�ffD�	�D�� D�@ D��3D�y�D�fD���D�\�D�3D�� D�<�D�� D�vfD� D�� D�L�D�� D��3D�9�D�� D��fD�#3D���D�ffD�  D��fD�@ D��D�� D�fD���D�c3D�	�D�� D�9�D��3D�y�D�3D�� D�\�D�  D��fD�C3D�� D�� D�&fD�� D�Y�D�  D��fD�@ D�� D�|�D��D���D�Y�D��fD��3D�@ D��3D�|�D�#3D�� D�Y�D�  D���D�33D���D�y�D��D���D�` D�fD���D�6fD���D���D�#3D¼�D�\�D�  Dģ3D�FfD��D�|�D�3DǼ�D�VfD��fDɖfD�6fD��fD�vfD�&fD�� D�Y�D�  DΦfD�<�D��3D�s3D��DѼ�D�S3D���DӠ D�FfD���D�s3D�fDֹ�D�c3D�fDؖfD�9�D���D�s3D��D�ɚD�i�D�	�Dݣ3D�C3D�� D�|�D��D๚D�S3D���D�3D�@ D��3D�3D�&fD��fD�\�D��3D��D�<�D���D� D�fD��D�c3D���D�fD�<�D��3D�y�D�fD�3D�S3D�� D��D�FfD��D�|�D��D��3D�ffD���D��3D�9�D���D��fD�)�D��f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A0  AA��AQ��Aa��Ap  A���A���A���A�  A�33A�  A�33A�  A���A�  A�  A���A�33A�ffA�33A�ffB ffBffBffB33BffBffB��B��B��B$ffB(  B+��B/��B3��B7��B;33B?33BC33BG33BK33BO��BS��BX  B\  B`��Bc��Bh  BlffBo��Bt  BxffB{��B~��B���B�  B���B�ffB���B�33B���B���B�33B���B�ffB�  B���B���B�ffB�  B���B�ffB�33B�  B���B���B���B�ffB�33B�  B���B���B���B�ffB�33B�33B�  B�  B���BǙ�Bʙ�B�ffB�33B�33B�33B�33B�33B�33B�33B�33B�33B�33B�ffB�ffC33CL�CL�CL�C	ffCffCffC� C��C�3C��C�fC� C  C�C33C!L�C#ffC%ffC'��C)��C+�3C-��C/ffC1  C3�C533C7L�C9� C;��C=�3C?��CAffCC  CE�CG33CIL�CKffCM��CO�3CQ��CSffCU  CW�CY�C[L�C]ffC_� Ca��Cc�3Ce��CgffCi  Ck�Cm33CoL�CqffCs� Cu��Cw�3Cy�fC{ffC}  C�C��fC��fC��3C���C�ٚC��fC��3C�� C���C���C��fC��3C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��fC��3C�� C���C��3C��3C��3C��3C��3C�� C�� C���C���C��fC��3C�� C���C��fC��3C���C��fC�� C�ٚC��3C���C�� C�ٚC�� C��fC���C��3C�ٚC�� C��fC���C�� C��fC���C��3C��fC���C�� C��3C�ٚC���C�� Cų3CƦfCǙ�CȌ�Cɳ3C��fC��fC�ٚC�ٚC�ٚC���C���C���C���C�ٚC�ٚC�ٚC��fC��3Cس3Cـ Cڀ Cی�Cܙ�CݦfC޳3C�� C�ٚC��fC�3C� C䙚C�fC�� C�ٚC�3C��C�fC�� C�ٚC��fC�� C��C�fC�� C��fC�� C��C��fC�� C�ٚC�� C���C�s3C��3D 9�Ds3D�3D��D,�D�fD��D�3D
L�D� D�3D�D@ Ds3D��D  D33D��D��D� D&fDy�D� D�fD&fDffD ��D!��D#FfD$�fD%��D&��D(S3D)��D*�fD,�D-L�D.�3D/��D0ٚD29�D3��D4ٚD5��D7  D8ffD9�fD;  D<FfD=y�D>�3D?�3DA9�DB�fDC� DD�3DFFfDG� DH��DI��DK@ DL��DM��DN�fDP33DQs3DR�3DS�3DU,�DV�fDW� DX�3DZ33D[l�D\� D]��D_33D`s3Da��DcfDdFfDey�Df��Dg�3Di33Djs3Dk��Dl��Dn9�Do� Dp�fDq�3Ds&fDts3Du�fDv��Dx9�Dys3Dz��D{�fD}9�D~� D�fD�y�D��D���D�c3D�3D��3D�FfD��3D��3D�&fD��fD�ffD���D��3D�9�D��3D��3D��D��3D�\�D�	�D��fD�FfD��fD��fD�&fD��fD�ffD�	�D�� D�@ D��3D�y�D�fD���D�\�D�3D�� D�<�D�� D�vfD� D�� D�L�D�� D��3D�9�D�� D��fD�#3D���D�ffD�  D��fD�@ D��D�� D�fD���D�c3D�	�D�� D�9�D��3D�y�D�3D�� D�\�D�  D��fD�C3D�� D�� D�&fD�� D�Y�D�  D��fD�@ D�� D�|�D��D���D�Y�D��fD��3D�@ D��3D�|�D�#3D�� D�Y�D�  D���D�33D���D�y�D��D���D�` D�fD���D�6fD���D���D�#3D¼�D�\�D�  Dģ3D�FfD��D�|�D�3DǼ�D�VfD��fDɖfD�6fD��fD�vfD�&fD�� D�Y�D�  DΦfD�<�D��3D�s3D��DѼ�D�S3D���DӠ D�FfD���D�s3D�fDֹ�D�c3D�fDؖfD�9�D���D�s3D��D�ɚD�i�D�	�Dݣ3D�C3D�� D�|�D��D๚D�S3D���D�3D�@ D��3D�3D�&fD��fD�\�D��3D��D�<�D���D� D�fD��D�c3D���D�fD�<�D��3D�y�D�fD�3D�S3D�� D��D�FfD��D�|�D��D��3D�ffD���D��3D�9�D���D��fD�)�D��f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����푿��#��`B�Ұ!�����S��e�T�<(��ƨ�\���;��<�`B>��?�?7�P?N�?@�?6E�?�?��?�7>��?V?�+?j?=p�?MO�?Z�?��?�O�?��?�33?���?���?��u?��?�|�?°!?�z�?�?�V?�+?�R?�Z?���@ �@(�@�j@��@�+@@K�@!�@$��@%`B@&E�@&ȴ@(r�@-`B@,�@/��@0�`@0  @.v�@-��@-p�@/��@2-@4j@4��@6�y@7l�@<1@=��@<��@<�@<z�@=/@;33@7;d@5O�@4Z@3�m@1��@3"�@1�#@1�7@/�P@/��@4Z@7�@6v�@3t�@17L@1�#@17L@:��@;�F@;��@:�!@81'@7��@7�@4Z@2J@0�u@1�#@1��@/�@/;d@.ȴ@.$�@.��@.ȴ@-�@-/@/|�@/l�@/K�@-�@+"�@+�
@+��@*�!@)��@)�#@)�@)��@)�#@)��@(��@(r�@'�w@'K�@&@&@%@%�h@%V@%��@%@$�/@$��@#S�@"��@#"�@#S�@#��@#C�@!�#@!��@ ��@|�@�@��@ 1'@ Q�@ r�@!hs@!�@!hs@!�@!&�@ Ĝ@�@�w@�w@   @   @|�@@O�@Z@1@��@J@�7@1'@�w@l�@��@`B@@��@��@��@Q�@b@��@��@��@��@|�@|�@|�@;d@�@@��@I�@Z@Z@Z@Z@j@I�@(�@��@33@
=q@	�^@	G�@	�@	&�@	7L@	G�@	�7@	��@	�7@	�^@	�@
M�@
~�@
~�@
�\@
��@
�!@
M�@
M�@	��@	�@
�@
�@
=q@
n�@
=q@
�@	�#@	��@
n�@
��@
�\@
�!@
�!@
�\@
~�@
-@
�@	��@	&�@�@�u@1'@K�@+@+@�R@ff@5?@5?@5?@p�@9X@�
@t�@"�@�\@x�@%@ �@ A�?�|�?��?�{?�I�?���?�^5?��u?���?���?�?�!?�hs?��?�;d?�?�(�?�1?��H?�+?�?}?㕁?�bN?�^5?�-?�O�?�X?�`B?�  ?�C�?� �?���?���?��?�(�?���?�S�?� �?��?�o?�V?�9X?���?z�?q��?j~�?g+?b��?c�
?`  ?^v�?YX?R�!?J��?Fff?=p�?5?}?,�D?%��?#��?"��?�?1?��?   >�~�>�b>���>�$�>�%>��j>��h>��y>�z�>}�>["�>L��>C��>>v�>9X>�w>bN>   =�=���=���=m�h=<j=+<���<��
<u���
��t���1�����'@��m�h������1��Q�ě�������%�\)��u��w�1&�?|�M��Q녾["Ѿe`B�vȴ�|푾�J���˾��^��=q��bN���Ͼ��+������(���;d��G���S����
��`B���V���׾�ȴ��p����۾�%��o��+��O߾����bN��t���z���ؓu��"Ѿݲ-��;d��Ĝ��Ĝ���
��`B��l���xվ���1��h������!��33��F��9X��E���^5�����ۿ   �%�G��J�����
��/�ff�l��r���ÿ	��
=q���1��ͿV�{�����;�hs����n��t���Ͽ9X�9X��������+��+��+��+�ȴ��P��u��#�"ѿ��j���/�/���w�   � Ĝ�!�7�"Mӿ#o�#���#���#�
�#�
�$Z�$���$���%��&��&�y�'+�'l��'(1'�(r��)7L�)7L�)7L�)xտ)xտ)��+�+C��+C��+C��,1�,�Ϳ-V�-O߿-�h�-�h�-��.���.��/��/\)�/���0�׿1&�1���2-�2�!�2�2�3�Ͽ5?}�5�6�6�6E��6�+�6�+�6ȴ�6ȴ�7
=�7�ٿ8b�8Q�8Q�8�u�8���9��9��9��9��:^5�:^5�:^5�:�H�;dZ�;dZ�;dZ�;��;�m�<j�<��<푿=p��=�=�>vɿ>vɿ>vɿ>�ۿ>�ۿ?�w�@A��@��@Ĝ�@��@��@��@Ĝ�AG��AG��AG��A%�A%11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ��푿��#��`B�Ұ!�����S��e�T�<(��ƨ�\���;��<�`B>��?�?7�P?N�?@�?6E�?�?��?�7>��?V?�+?j?=p�?MO�?Z�?��?�O�?��?�33?���?���?��u?��?�|�?°!?�z�?�?�V?�+?�R?�Z?���@ �@(�@�j@��@�+@@K�@!�@$��@%`B@&E�@&ȴ@(r�@-`B@,�@/��@0�`@0  @.v�@-��@-p�@/��@2-@4j@4��@6�y@7l�@<1@=��@<��@<�@<z�@=/@;33@7;d@5O�@4Z@3�m@1��@3"�@1�#@1�7@/�P@/��@4Z@7�@6v�@3t�@17L@1�#@17L@:��@;�F@;��@:�!@81'@7��@7�@4Z@2J@0�u@1�#@1��@/�@/;d@.ȴ@.$�@.��@.ȴ@-�@-/@/|�@/l�@/K�@-�@+"�@+�
@+��@*�!@)��@)�#@)�@)��@)�#@)��@(��@(r�@'�w@'K�@&@&@%@%�h@%V@%��@%@$�/@$��@#S�@"��@#"�@#S�@#��@#C�@!�#@!��@ ��@|�@�@��@ 1'@ Q�@ r�@!hs@!�@!hs@!�@!&�@ Ĝ@�@�w@�w@   @   @|�@@O�@Z@1@��@J@�7@1'@�w@l�@��@`B@@��@��@��@Q�@b@��@��@��@��@|�@|�@|�@;d@�@@��@I�@Z@Z@Z@Z@j@I�@(�@��@33@
=q@	�^@	G�@	�@	&�@	7L@	G�@	�7@	��@	�7@	�^@	�@
M�@
~�@
~�@
�\@
��@
�!@
M�@
M�@	��@	�@
�@
�@
=q@
n�@
=q@
�@	�#@	��@
n�@
��@
�\@
�!@
�!@
�\@
~�@
-@
�@	��@	&�@�@�u@1'@K�@+@+@�R@ff@5?@5?@5?@p�@9X@�
@t�@"�@�\@x�@%@ �@ A�?�|�?��?�{?�I�?���?�^5?��u?���?���?�?�!?�hs?��?�;d?�?�(�?�1?��H?�+?�?}?㕁?�bN?�^5?�-?�O�?�X?�`B?�  ?�C�?� �?���?���?��?�(�?���?�S�?� �?��?�o?�V?�9X?���?z�?q��?j~�?g+?b��?c�
?`  ?^v�?YX?R�!?J��?Fff?=p�?5?}?,�D?%��?#��?"��?�?1?��?   >�~�>�b>���>�$�>�%>��j>��h>��y>�z�>}�>["�>L��>C��>>v�>9X>�w>bN>   =�=���=���=m�h=<j=+<���<��
<u���
��t���1�����'@��m�h������1��Q�ě�������%�\)��u��w�1&�?|�M��Q녾["Ѿe`B�vȴ�|푾�J���˾��^��=q��bN���Ͼ��+������(���;d��G���S����
��`B���V���׾�ȴ��p����۾�%��o��+��O߾����bN��t���z���ؓu��"Ѿݲ-��;d��Ĝ��Ĝ���
��`B��l���xվ���1��h������!��33��F��9X��E���^5�����ۿ   �%�G��J�����
��/�ff�l��r���ÿ	��
=q���1��ͿV�{�����;�hs����n��t���Ͽ9X�9X��������+��+��+��+�ȴ��P��u��#�"ѿ��j���/�/���w�   � Ĝ�!�7�"Mӿ#o�#���#���#�
�#�
�$Z�$���$���%��&��&�y�'+�'l��'(1'�(r��)7L�)7L�)7L�)xտ)xտ)��+�+C��+C��+C��,1�,�Ϳ-V�-O߿-�h�-�h�-��.���.��/��/\)�/���0�׿1&�1���2-�2�!�2�2�3�Ͽ5?}�5�6�6�6E��6�+�6�+�6ȴ�6ȴ�7
=�7�ٿ8b�8Q�8Q�8�u�8���9��9��9��9��:^5�:^5�:^5�:�H�;dZ�;dZ�;dZ�;��;�m�<j�<��<푿=p��=�=�>vɿ>vɿ>vɿ>�ۿ>�ۿ?�w�@A��@��@Ĝ�@��@��@��@Ĝ�AG��AG��AG��A%�A%11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B�?B��B��B��B��BE�Bx�B��B�)B'�B�hB��B<jB��B�B	B	YB	m�B	� B	��B	�^B	�HB	�B
+B
�B
'�B
S�B
jB
y�B
|�B
�B
�=B
ŢB
�B
�wB
�-B
��B
�B
��B
�BB
=B+B�B�B;dBB�B@�BF�BVBiyB]/BgmBk�Bp�Bs�Bu�Bz�B�+B�%B�PB�JB�JB�\B�PB��B��B��B��B��B��B�bB�?B�FB�?B�FB�?B�dB�9B�XB�!B�!B�?B�9B�'B�B�!B�!B�XB�wB�qB�dB�jB�^B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B��B�B�
B�B�
B��B��B�B�
B�B�B�B�B�B�B�B�B�B�#B�B�#B�;B�)B�#B�#B�#B�)B�#B�;B�/B�/B�)B�;B�#B�B�#B�)B�/B�5B�;B�;B�;B�BB�BB�5B�BB�;B�;B�BB�BB�BB�;B�/B�5B�#B�B�B�B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBƨBƨBŢBĜBĜBŢBĜBŢBŢBŢBŢBŢBƨBǮBǮBǮBǮBȴBȴBȴBǮBȴBǮBɺBȴBɺBɺBɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺB��B��B��BɺBɺBȴBɺBɺB��BɺBǮBȴBǮBɺBǮBŢBƨBƨBĜBƨBÖBÖBŢBÖBBBBB�}B�qB�jB�dB�dB�^B�^B�9B�-B�!B�B�B�B�B�B�B�9B�?B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�8B�WB��B��B�HB�,B�B�mBI.B|aB�JBߵB+|B��B�[B?�B�&B�B	�B	\�B	qB	��B	�DB	��B	��B	�B

�B
 8B
+|B
W�B
nB
}gB
�zB
��B
��B
�.B
��B
�B
��B
�qB
�<B
�_B
��B�B
�B&B!>B>�BFBDBJ4BY�BmB`�Bj�BoBt0BwBByOB~mB��B��B��B��B��B��B��B�B�B�DB�WB��B�iB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B�B�MB�~B�~B�~B�xB�xB�kB�_B�_B�YB�kB�qB�qB�xB�~B�~BׄB؊BِBۜBݩBܣBݩBׄBِBږBِBږB؊B؊BِBږBܣBݩBܣBݩBݩBܣBܣBݩBܣBޯBܣBޯB��BߵBޯBޯBޯBߵBޯB��B�B�BߵB��BޯBܣBޯBߵB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��BޯBݩBݩBݩBۜBږBِB؊B؊B�qB�qB�qB�_B�eB�_B�_B�YB�YB�_B�YB�YB�YB�YB�MB�MB�FB�@B�@B�@B�@B�@B�@B�@B�@B�@B�4B�4B�.B�(B�(B�.B�(B�.B�.B�.B�.B�.B�4B�:B�:B�:B�:B�@B�@B�@B�:B�@B�:B�FB�@B�FB�FB�FB�FB�FB�FB�MB�SB�YB�YB�YB�YB�_B�YB�_B�_B�YB�SB�YB�SB�SB�SB�YB�SB�SB�SB�SB�SB�SB�SB�MB�FB�MB�MB�MB�FB�FB�@B�FB�FB�MB�FB�:B�@B�:B�FB�:B�.B�4B�4B�(B�4B�"B�"B�.B�"B�B�B�B�B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�uB�iB�iB�iB�]B�]B�QB�QB�WB�WB�JB�QB�JB�DB�8B�2B�2B�2B�2B�2B�2B�2B�2B�8B�8B�8B�8B�,B�&B�,B�,B�2B�,B�,B�,B�2B�2B�2B�2B�>B�8B�8B�>B�>B�>B�>B�DB�>B�>B�>B�8B�>B�>B�>B�>B�>B�>B�DB�>B�>B�DB�>B�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�QB�QB�QB�QB�QB�QB�QB�QB�QB�QB�QB�QB�QB�QB�WB�QB�QB�QB�WB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�]B�]B�]B�]B�cB�]B�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�iB�iB�oB�oB�oB�iB�oB�oB�oB�oB�oB�oB�uB�oB�uB�uB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542392021060715423920210607154239  IF  ARFMCODA029d                                                                20190620115650                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115719  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.2                                                                 20190620115719  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154239  IP  PSAL            A0  D��fG�O�                