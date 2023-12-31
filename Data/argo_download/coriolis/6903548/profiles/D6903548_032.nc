CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  O   	N_HISTORY          N_CALIB          
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
resolution        =���   axis      Z        	<  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	<  Fd   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	<  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  [,   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  dh   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  f�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  o�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  rD   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  {�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �H   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	<  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �0   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �4   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �8   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �<   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �@   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �Argo profile    3.1 1.2 19500101000000  20190620115650  20210607154239  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL                A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @؈���>�1   @؈���>�@TS6H��6@*�f5��`8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A.ffA@  AS33Aa��Ap  A~ffA�  A���A���A���A���A�33A�ffA���A�  A���A�  A�33A�ffA�  A���A�33B33B33B33B��B��B��BffB��B"��B'��B,ffB/��B4  B8ffB;33B?33BC��BG��BK��BO��BT  BX  B\ffB_33Bc��BhffBk��Bp  Bu33BxffB{��B~��B���B�ffB�  B���B�ffB�33B���B�ffB�33B�  B���B�ffB�33B�33B���B���B���B���B�ffB�ffB�ffB�ffB�33B�33B�  B�  B�  B�  B�  B���B���B���B���BÙ�Bř�B�ffB�ffB�33B�33B�  B���Bߙ�B㙚B�ffB�ffB�ffB�ffB�33B�33B�33C� C��C��C��C	�3C�3C��C��C�fC� C  C�CL�CffC��C�3C!��C#ffC%  C'  C)�C+L�C-L�C/L�C1L�C3L�C5ffC7ffC9ffC;ffC=� C?33CAL�CCffCE33CGL�CI� CK33CML�CO� CQ33CSL�CU� CWL�CY� C[��C]ffC_33Ca� Cc�3Ce� CgL�Ci� Ck�3Cm� CoL�Cq�CsffCu�3Cw� CyL�C{33C}� C�3C���C��3C��fC���C�� C��3C�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C�� C��3C��fC��fC���C���C�� C��3C��fC���C�� C��3C��fC���C���C�� C��3C��3C��3C��fC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C�� C���C���C���C���C��fC��fC��3C�� C�� C���C�ٚC��fC��3C�� C�� C���C��fC³3C�� C�� Cų3CƦfCǦfCȦfC���C���C���C���C͙�CΦfCϦfCЦfCѳ3Cҳ3C�� C�� C�ٚC֦fC׳3C�� Cٙ�CڦfC�� Cܙ�CݦfC�� Cߙ�C�fC�� C♚C�fC�� C噚C�� C�ٚC���C�3CꙚC�� C�ٚC�� C�3C�fC�C��C�3C�ٚC���C�� C��3C��fC���C�� C�ffC��D 9�Dl�D� D�3DFfD�fD��D��D
33Dl�D�fD  DS3D�3D�3D3DY�Dy�D� D��D9�D��D��D�fD33D�fD �3D!� D#,�D$s3D%��D&�3D(&fD)ffD*� D+�fD-&fD.s3D/� D1  D29�D3� D4� D6fD7,�D8y�D9�fD:��D<9�D=��D>��D?��DA&fDBy�DC��DE  DF,�DG� DH��DI��DK,�DL� DM�3DO  DP,�DQy�DR�fDS��DU,�DVffDW��DY3DZL�D[� D\��D]�3D_,�D`l�Da�fDb�fDd&fDeffDf�fDg�3Di@ Dj��Dk��Dl�3Dn@ Do�3Dp�fDr  Ds@ Dt�fDu�fDw3Dx9�DyffDz�3D|fD}@ D~� D��D�� D�fD�� D�ffD�  D���D�9�D��fD�vfD�fD��fD�VfD��fD���D�<�D�ٚD�� D�&fD��3D�c3D�fD��fD�<�D��fD�|�D��D��3D�` D�3D��3D�FfD�ٚD�|�D�#3D���D�\�D�fD�� D�<�D�ٚD�y�D��D���D�Y�D���D���D�@ D��3D��fD�&fD�ɚD�\�D�� D��fD�<�D��fD�y�D��D���D�` D�  D���D�<�D�� D��fD��D��3D�Y�D�  D���D�6fD�� D�y�D�fD�� D�Y�D��fD�� D�9�D��fD��3D�  D���D�\�D���D���D�<�D�� D�s3D��D���D�S3D���D��3D�9�D��3D�� D�)�D�ɚD�i�D�	�D���D�0 D���D�|�D�fD³3D�P D�� Dē3D�6fD���DƆfD�#3D�� D�c3D��fDə�D�<�D��fDˀ D�)�D��fD�ffD�3DΣ3D�C3D��fDЉ�D��Dѳ3D�Y�D�fDӠ D�9�D�� DՃ3D�  Dּ�D�\�D���D؜�D�<�D���Dڀ D�#3D��fD�\�D��3Dݜ�D�FfD��3D�s3D�fD��D�S3D���D� D�<�D��fD�|�D�)�D��3D�\�D���D�3D�<�D��D�fD�,�D��fD�ffD�3D� D�@ D�� D� D�&fD��3D�\�D�fD�3D�<�D��3D�fD��D�� D�S3D���D���D�6fD���D���D�  D���D�0 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A.ffA@  AS33Aa��Ap  A~ffA�  A���A���A���A���A�33A�ffA���A�  A���A�  A�33A�ffA�  A���A�33B33B33B33B��B��B��BffB��B"��B'��B,ffB/��B4  B8ffB;33B?33BC��BG��BK��BO��BT  BX  B\ffB_33Bc��BhffBk��Bp  Bu33BxffB{��B~��B���B�ffB�  B���B�ffB�33B���B�ffB�33B�  B���B�ffB�33B�33B���B���B���B���B�ffB�ffB�ffB�ffB�33B�33B�  B�  B�  B�  B�  B���B���B���B���BÙ�Bř�B�ffB�ffB�33B�33B�  B���Bߙ�B㙚B�ffB�ffB�ffB�ffB�33B�33B�33C� C��C��C��C	�3C�3C��C��C�fC� C  C�CL�CffC��C�3C!��C#ffC%  C'  C)�C+L�C-L�C/L�C1L�C3L�C5ffC7ffC9ffC;ffC=� C?33CAL�CCffCE33CGL�CI� CK33CML�CO� CQ33CSL�CU� CWL�CY� C[��C]ffC_33Ca� Cc�3Ce� CgL�Ci� Ck�3Cm� CoL�Cq�CsffCu�3Cw� CyL�C{33C}� C�3C���C��3C��fC���C�� C��3C�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C�� C��3C��fC��fC���C���C�� C��3C��fC���C�� C��3C��fC���C���C�� C��3C��3C��3C��fC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C�� C���C���C���C���C��fC��fC��3C�� C�� C���C�ٚC��fC��3C�� C�� C���C��fC³3C�� C�� Cų3CƦfCǦfCȦfC���C���C���C���C͙�CΦfCϦfCЦfCѳ3Cҳ3C�� C�� C�ٚC֦fC׳3C�� Cٙ�CڦfC�� Cܙ�CݦfC�� Cߙ�C�fC�� C♚C�fC�� C噚C�� C�ٚC���C�3CꙚC�� C�ٚC�� C�3C�fC�C��C�3C�ٚC���C�� C��3C��fC���C�� C�ffC��D 9�Dl�D� D�3DFfD�fD��D��D
33Dl�D�fD  DS3D�3D�3D3DY�Dy�D� D��D9�D��D��D�fD33D�fD �3D!� D#,�D$s3D%��D&�3D(&fD)ffD*� D+�fD-&fD.s3D/� D1  D29�D3� D4� D6fD7,�D8y�D9�fD:��D<9�D=��D>��D?��DA&fDBy�DC��DE  DF,�DG� DH��DI��DK,�DL� DM�3DO  DP,�DQy�DR�fDS��DU,�DVffDW��DY3DZL�D[� D\��D]�3D_,�D`l�Da�fDb�fDd&fDeffDf�fDg�3Di@ Dj��Dk��Dl�3Dn@ Do�3Dp�fDr  Ds@ Dt�fDu�fDw3Dx9�DyffDz�3D|fD}@ D~� D��D�� D�fD�� D�ffD�  D���D�9�D��fD�vfD�fD��fD�VfD��fD���D�<�D�ٚD�� D�&fD��3D�c3D�fD��fD�<�D��fD�|�D��D��3D�` D�3D��3D�FfD�ٚD�|�D�#3D���D�\�D�fD�� D�<�D�ٚD�y�D��D���D�Y�D���D���D�@ D��3D��fD�&fD�ɚD�\�D�� D��fD�<�D��fD�y�D��D���D�` D�  D���D�<�D�� D��fD��D��3D�Y�D�  D���D�6fD�� D�y�D�fD�� D�Y�D��fD�� D�9�D��fD��3D�  D���D�\�D���D���D�<�D�� D�s3D��D���D�S3D���D��3D�9�D��3D�� D�)�D�ɚD�i�D�	�D���D�0 D���D�|�D�fD³3D�P D�� Dē3D�6fD���DƆfD�#3D�� D�c3D��fDə�D�<�D��fDˀ D�)�D��fD�ffD�3DΣ3D�C3D��fDЉ�D��Dѳ3D�Y�D�fDӠ D�9�D�� DՃ3D�  Dּ�D�\�D���D؜�D�<�D���Dڀ D�#3D��fD�\�D��3Dݜ�D�FfD��3D�s3D�fD��D�S3D���D� D�<�D��fD�|�D�)�D��3D�\�D���D�3D�<�D��D�fD�,�D��fD�ffD�3D� D�@ D�� D� D�&fD��3D�\�D�fD�3D�<�D��3D�fD��D�� D�S3D���D���D�6fD���D���D�  D���D�0 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����p���푿�푿�j�ؓu��|������D�����D���#�vȴ�>�ۿ�˾�;d=L��?��?�E�?�o?��/?�{@ Q�@��@v�@(�@ Ĝ@V@��@��@��@Z@
�H?�|�?ܬ?��9?�9X?��/?�bN?���?�1'?�?tz�?c�
?c��?b�\?|�?��?�&�?��H?��P?��-?��?̋D?�p�?�Z?�-?�n�?�j?��/?�o@|�@@�@  �@%/@!�^@!�#@'��@+33@0�9@4j@9�7@;"�@<��@=V@=��@>@;�m@9X@9hs@97L@9x�@9��@:M�@;��@;"�@;t�@;��@9x�@7�@2��@(��@&��@&E�@)�#@.��@/�@8�`@9��@9��@9��@8��@7�@4�/@3"�@-@*�H@'\)@$1@!hs@;d@ bN@!�#@'��@(��@*=q@*n�@+@+o@+ƨ@,�@,(�@,�/@,��@,�D@-`B@.��@/|�@1��@1��@17L@/;d@-��@-O�@-V@,�/@,��@,Z@,(�@+ƨ@*�H@*^5@*J@)��@)hs@)7L@(��@'�;@'
=@&ff@%�T@%�@$�j@$Z@#�m@#�
@#��@#S�@#dZ@#dZ@#t�@#dZ@#t�@#t�@#t�@#S�@"�@"�@!7L@ ��@��@;d@+@;d@\)@\)@�R@@?}@z�@��@@M�@��@x�@�@��@�`@��@Q�@  @  @��@$�@��@(�@M�@hs@Ĝ@��@G�@��@�@J@��@��@�7@G�@�@�9@�9@1'@  @�@��@ �@bN@�@�@bN@��@Z@I�@
�@
�@
�@
�@
�@��@(�@Z@z�@��@��@�j@j@dZ@
~�@	&�@	�7@	�#@
��@
~�@
n�@
�\@
��@o@o@
�@
��@
��@	�^@	X@	x�@	�#@��@�;@|�@v�@@�@��@�@?}@j@(�@1@�@�H@^5@�@ b?�j?��?�?��?�=q?�b?�+?��+?�$�?�
=?�ȴ?�-?�bN?��?�O�?�v�?�5??�h?��#?�7?��H?��y?ӕ�?�|�?�j?���?���?���?�`B?��!?���?��D?�K�?���?�n�?�v�?�ƨ?��?��
?�%?���?�ƨ?��?�t�?~v�?{"�?w��?q�?h�9?a%?W�P?Qhs?J=q?A%?7��?1��?,I�?#�
?;d?�#?z�?\)?1'?��>���>�ff>ۥ�>�bN>���>�ƨ>�o>���>��T>�S�>�hs>���>|�>u>hr�>G�>?|�>,1>%�T> Ĝ>$�=��=�F=�=�E�=�\)=49X=C�<�/;�`B�D����o�49X�T���� Ž�/��S����m�J�o�
=q�\)����-�)��7KǾA�7�E�˾I�^�Xb�ixվk��m�h�|푾������˾�$ݾ����;��;��n����P��"Ѿ�A���MӾ�`B���D�� ž����پ�X���H��p���|�\�š˾Ƨ�ȴ9��=q��O߾�녾�z����b����ۥ��/�޸R��A���G�����������/���/���T��xվ���h��&�����33��9X��E���ȴ���پ��H��푾�p��   � Ĝ�G��J������������˿�y�1'��9�	�^�	���C��1��D��h�V�������;��׿&�녿-�33��j��E��ȴ��ٿQ������X�^5��m�j�푿�-�;d��w� A��!%�"J�"�\�#�
�%��%�˿%�T�%�T�&ff�&�y�'(r��)7L�)�^�+C��+��+ƨ�+ƨ�,I��,I��,�D�-O߿.{�.V�.V�.��/��/��/\)�/���0 ſ0�`�1hs�1hs�1���1녿1녿2-�2n��2�!�2�2�333�333�3�F�3�F�49X�4�j�5?}�5��6�6E��6�+�6�+�6�+�6ȴ�6ȴ�6ȴ�6ȴ�7
=�7
=�7
=�7
=�7
=�7Kǿ7Kǿ7�P�7�P�7�P�7�ٿ8b�8���9��9��9��9X�9���:^5�:�H�;dZ�;��<(��<(��<j�=/�=/�=/�=p��=�=�=�>5?�>vɿ>�R�>�ۿ>�ۿ>�ۿ>�ۿ>�ۿ>�ۿ>��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��p���푿�푿�j�ؓu��|������D�����D���#�vȴ�>�ۿ�˾�;d=L��?��?�E�?�o?��/?�{@ Q�@��@v�@(�@ Ĝ@V@��@��@��@Z@
�H?�|�?ܬ?��9?�9X?��/?�bN?���?�1'?�?tz�?c�
?c��?b�\?|�?��?�&�?��H?��P?��-?��?̋D?�p�?�Z?�-?�n�?�j?��/?�o@|�@@�@  �@%/@!�^@!�#@'��@+33@0�9@4j@9�7@;"�@<��@=V@=��@>@;�m@9X@9hs@97L@9x�@9��@:M�@;��@;"�@;t�@;��@9x�@7�@2��@(��@&��@&E�@)�#@.��@/�@8�`@9��@9��@9��@8��@7�@4�/@3"�@-@*�H@'\)@$1@!hs@;d@ bN@!�#@'��@(��@*=q@*n�@+@+o@+ƨ@,�@,(�@,�/@,��@,�D@-`B@.��@/|�@1��@1��@17L@/;d@-��@-O�@-V@,�/@,��@,Z@,(�@+ƨ@*�H@*^5@*J@)��@)hs@)7L@(��@'�;@'
=@&ff@%�T@%�@$�j@$Z@#�m@#�
@#��@#S�@#dZ@#dZ@#t�@#dZ@#t�@#t�@#t�@#S�@"�@"�@!7L@ ��@��@;d@+@;d@\)@\)@�R@@?}@z�@��@@M�@��@x�@�@��@�`@��@Q�@  @  @��@$�@��@(�@M�@hs@Ĝ@��@G�@��@�@J@��@��@�7@G�@�@�9@�9@1'@  @�@��@ �@bN@�@�@bN@��@Z@I�@
�@
�@
�@
�@
�@��@(�@Z@z�@��@��@�j@j@dZ@
~�@	&�@	�7@	�#@
��@
~�@
n�@
�\@
��@o@o@
�@
��@
��@	�^@	X@	x�@	�#@��@�;@|�@v�@@�@��@�@?}@j@(�@1@�@�H@^5@�@ b?�j?��?�?��?�=q?�b?�+?��+?�$�?�
=?�ȴ?�-?�bN?��?�O�?�v�?�5??�h?��#?�7?��H?��y?ӕ�?�|�?�j?���?���?���?�`B?��!?���?��D?�K�?���?�n�?�v�?�ƨ?��?��
?�%?���?�ƨ?��?�t�?~v�?{"�?w��?q�?h�9?a%?W�P?Qhs?J=q?A%?7��?1��?,I�?#�
?;d?�#?z�?\)?1'?��>���>�ff>ۥ�>�bN>���>�ƨ>�o>���>��T>�S�>�hs>���>|�>u>hr�>G�>?|�>,1>%�T> Ĝ>$�=��=�F=�=�E�=�\)=49X=C�<�/;�`B�D����o�49X�T���� Ž�/��S����m�J�o�
=q�\)����-�)��7KǾA�7�E�˾I�^�Xb�ixվk��m�h�|푾������˾�$ݾ����;��;��n����P��"Ѿ�A���MӾ�`B���D�� ž����پ�X���H��p���|�\�š˾Ƨ�ȴ9��=q��O߾�녾�z����b����ۥ��/�޸R��A���G�����������/���/���T��xվ���h��&�����33��9X��E���ȴ���پ��H��푾�p��   � Ĝ�G��J������������˿�y�1'��9�	�^�	���C��1��D��h�V�������;��׿&�녿-�33��j��E��ȴ��ٿQ������X�^5��m�j�푿�-�;d��w� A��!%�"J�"�\�#�
�%��%�˿%�T�%�T�&ff�&�y�'(r��)7L�)�^�+C��+��+ƨ�+ƨ�,I��,I��,�D�-O߿.{�.V�.V�.��/��/��/\)�/���0 ſ0�`�1hs�1hs�1���1녿1녿2-�2n��2�!�2�2�333�333�3�F�3�F�49X�4�j�5?}�5��6�6E��6�+�6�+�6�+�6ȴ�6ȴ�6ȴ�6ȴ�7
=�7
=�7
=�7
=�7
=�7Kǿ7Kǿ7�P�7�P�7�P�7�ٿ8b�8���9��9��9��9X�9���:^5�:�H�;dZ�;��<(��<(��<j�=/�=/�=/�=p��=�=�=�>5?�>vɿ>�R�>�ۿ>�ۿ>�ۿ>�ۿ>�ۿ>�ۿ>��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�}B��BÖBŢB��B�;B��B(�BT�Br�B~�B��B�TB �BC�B�yB�B��BaHB�!B��B	/B	C�B	B�B	�+B	�B	��B	�B
  B
+B

=B
"�B
B
%B
{B
�B
�B
#�B
!�B
!�B
J�B
L�B
T�B
[#B
dZB
r�B
�B
�PB
�DB
ȴB
��B
�B
�sB
�sB
�sBBB1BVB�B{BH�BcTB`BBjBo�Bo�B^5B�DB�\B��B��B��B��B��B�B�B��B��B��B��B�B�B�B�'B�3B��B�-B�RB�-B�B��B��B��B��B�9B�?BB��B��BȴBɺBǮBŢB�qB�}B�?B�jB�!B�B�!B�9B�^BŢB�qB��B��B��B��B��B�B��B�B�B�#B�/B�HB�NB�B�B�B�B�yB�sB�sB�B�B�yB�yB�yB�mB�mB�mB�mB�mB�fB�mB�fB�`B�`B�ZB�TB�TB�TB�TB�TB�TB�TB�TB�ZB�ZB�ZB�`B�`B�ZB�`B�ZB�TB�NB�NB�HB�NB�HB�BB�HB�BB�;B�/B�5B�5B�)B�)B�#B�B�B�B�B�B�B�B�
B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȴBȴBȴBȴBɺBɺB��B��B��B��B��B��B��BɺBɺBǮBƨBȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBȴBɺBȴBɺBȴBȴBǮBȴBƨBĜBĜBĜBÖBÖB��BBĜBĜBÖBBB��B��B�}B��B�}B��B�qB�dB�RB�RB�RB�RB�RB�RB�LB�?B�3B�-B�3B�!B�B�!B�'B�!B�!B�B�!B�!B�!B�B�!B�B�!B�!B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�	B�B�"B�.B�YB��B�UB,�BX�Bv<B��B�]B��B$QBG"B�B��B�Bd�B��B�aB	2�B	G"B	FB	��B	��B	�uB	ۜB
�B

�B
�B
&]B
�B
	�B
B
!>B
#JB
'cB
%WB
%WB
NMB
PYB
X�B
^�B
g�B
v<B
��B
��B
��B
�@B
�eB
ِB
��B
��B
��B�B�B�B�BBBL@Bf�Bc�BnBs*Bs*Ba�B��B��B�B�iB�]B�8B��B��B��B��B��B�|B��B��B��B��B��B��B�cB��B��B��B��B�cB�8B�8B�B��B��B�B�MB�MB�@B�FB�:B�.B��B�	B��B��B��B��B��B��B��B�.B��B�SB�YB�kB�~BׄBِB؊BۜBܣBޯB�B��B��B�B�B�$B�B�B��B��B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��BߵBߵBޯBݩBݩBܣBܣBܣBܣBۜBږBږBۜBِBׄB�xB�eB�eB�eB�kB�kB�qB�xB�qB�kB�xB�kB�qB�xB�~B�eB�_B�eB�eB�eB�eB�kB�kB�kB�kB�YB�MB�SB�@B�@B�@B�@B�FB�FB�MB�SB�SB�SB�SB�MB�YB�FB�FB�:B�4B�@B�MB�MB�MB�MB�SB�YB�_B�_B�_B�YB�SB�YB�YB�YB�YB�MB�MB�SB�FB�FB�FB�FB�@B�FB�@B�FB�@B�@B�:B�@B�4B�(B�(B�(B�"B�"B�B�B�(B�(B�"B�B�B�B�B�	B�B�	B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�uB�uB�|B�|B�|B�uB�|B�oB�iB�cB�]B�]B�WB�WB�WB�QB�JB�DB�DB�JB�DB�>B�>B�>B�>B�>B�8B�>B�>B�>B�8B�8B�8B�>B�>B�DB�DB�>B�JB�DB�JB�DB�8B�2B�2B�2B�2B�2B�2B�8B�8B�8B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�DB�DB�DB�DB�DB�JB�JB�DB�JB�JB�JB�JB�JB�QB�QB�JB�QB�QB�QB�QB�QB�QB�QB�QB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�]B�]B�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�oB�oB�oB�oB�iB�iB�oB�oB�oB�oB�uB�oB�oB�uB�oB�oB�oB�oB�uB�uB�uB�uB�uB�|B�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B��B��B�|B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542392021060715423920210607154239  IF  ARFMCODA029d                                                                20190620115650                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115720  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.2                                                                 20190620115720  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154239  IP  PSAL            A.ffD�0 G�O�                