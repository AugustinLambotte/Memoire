CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  U   	N_HISTORY          N_CALIB          
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
resolution        =���   axis      Z        	T  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  D,   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	T  F�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	T  R0   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  [�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  g0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  p�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  |0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  �0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �8   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �<   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �@   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �D   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �H   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �Argo profile    3.1 1.2 19500101000000  20200829092229  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               5A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @���1   @���@R����gl�#�Զ�R8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      @���@�33A��A��A!��A1��AA��ANffA`  As33A~ffA�33A�  A���A���A���A�ffA���A�  A�  A�33A���A�  A���A���A�ffB ffB  B��B  B  B��B��B  B ffB#��B&��B+��B0��B4  B7��B:��B>ffBC��BH��BLffBP  BS��BW33BZ��B^ffBc��Bi33Bm33Bp��Bt��Bx��B|ffB�ffB�ffB�ffB�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB���B���B���B���B���B���B���B���B���B���B���B���B�  B�  B�  B�33B�33B�33B�33B�ffB�33B�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�B���B���B���C� C��C�3C��C	�fC� C  C33CL�CL�CffCffC� C��C�3C�3C!��C#ffC%  C'  C)�C+ffC-� C/� C1� C3� C5��C733C9L�C;� C=33C?L�CAffCC�CEL�CGffCI� CK��CML�COffCQ� CS33CUffCW��CYL�C[ffC]��C_L�CaffCc� Ce33CgffCi� CkL�Cm� Co�3CqffCs33CuffCw� CyL�C{�C}ffC��C��3C���C�� C�ٚC�� C��3C���C�� C��fC���C��3C��fC���C�� C��fC���C�� C��fC���C�� C��3C�ٚC���C�� C��3C��fC���C���C��3C��fC�ٚC���C�� C��3C���C���C���C��3C��fC�ٚC���C�� C��fC���C���C���C�� C��3C��fC��fC�ٚC�ٚC���C���C�� C��3C��3C��fC���C���C���C��3C��3C��fC�ٚC�ٚC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C�� Cͳ3C���C���C���C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C���C���C�ٚCަfCߦfC�3C�� C���C�3C�� C�ٚC�3C�� C�ٚC�3C�� C���C�fC��3C�� C���C�ٚC�3C�3C���C��fC�� C�ٚC��3C���C��3C�� C��3D 9�D� D��D��D,�D` D��D	�D
FfD� D� D  D9�D� D�fD��D33Ds3D��D��DFfD��D� D��D9�D��D � D!��D#9�D$s3D%�3D&��D(,�D)s3D*� D+�3D-  D.s3D/�fD1  D29�D3s3D4�3D5��D733D8s3D9��D;  D<L�D=y�D>� D?��DA9�DB�fDC��DD�3DF&fDG` DH� DI�fDK,�DLs3DM� DOfDP9�DQl�DR� DS�3DU,�DV�fDW�fDX��DZ33D[y�D\� D^�D_@ D`s3Da�fDb��Dd33De� Df�fDg�3Di&fDjy�Dk��DmfDnFfDo�fDp�fDrfDsL�Dts3Du� Dv��Dx@ Dys3Dz�fD{� D}  D~` D��D�y�D��D��3D�Y�D��fD���D�9�D��fD��3D�  D�� D�` D��3D���D�@ D�ٚD�s3D��D�ɚD�ffD�3D���D�9�D�ٚD�vfD�fD��fD�VfD���D�� D�FfD���D�s3D��D�� D�ffD���D���D�9�D�ٚD�y�D��D���D�` D�fD���D�6fD���D��fD�#3D�� D�` D���D���D�<�D�� D��3D�&fD���D�S3D���D��3D�<�D��fD�s3D�  D�ɚD�ffD�fD��fD�FfD��D���D�,�D�� D�S3D���D�� D�FfD���D�p D��D���D�S3D��fD�� D�9�D��3D�|�D��D���D�` D��3D���D�@ D��fD�� D�)�D��3D�` D���D��fD�33D�� D�p D� D° D�P D�� Dē3D�33D���D�|�D�fDǳ3D�P D�� Dɓ3D�6fD���DˆfD�  D̹�D�` D���DΓ3D�<�D�ٚD�vfD�  D�� D�` D�  DӖfD�9�D���DՀ D�&fD�� D�Y�D�  Dة�D�C3D��3Dڀ D�  D�� D�c3D�	�Dݜ�D�@ D��fD߃3D�&fD��D�VfD���D�fD�C3D�� D� D�  D�� D�\�D���D��D�<�D�� D�fD�  D�fD�\�D�fD� D�9�D��fD�s3D�3D� D�\�D�	�D�D�I�D���D� D�3D�� D�` D�  D���D�<�D���D�� D�#3D���D�S3D�  D�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @���@�33A��A��A!��A1��AA��ANffA`  As33A~ffA�33A�  A���A���A���A�ffA���A�  A�  A�33A���A�  A���A���A�ffB ffB  B��B  B  B��B��B  B ffB#��B&��B+��B0��B4  B7��B:��B>ffBC��BH��BLffBP  BS��BW33BZ��B^ffBc��Bi33Bm33Bp��Bt��Bx��B|ffB�ffB�ffB�ffB�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB���B���B���B���B���B���B���B���B���B���B���B���B�  B�  B�  B�33B�33B�33B�33B�ffB�33B�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�B���B���B���C� C��C�3C��C	�fC� C  C33CL�CL�CffCffC� C��C�3C�3C!��C#ffC%  C'  C)�C+ffC-� C/� C1� C3� C5��C733C9L�C;� C=33C?L�CAffCC�CEL�CGffCI� CK��CML�COffCQ� CS33CUffCW��CYL�C[ffC]��C_L�CaffCc� Ce33CgffCi� CkL�Cm� Co�3CqffCs33CuffCw� CyL�C{�C}ffC��C��3C���C�� C�ٚC�� C��3C���C�� C��fC���C��3C��fC���C�� C��fC���C�� C��fC���C�� C��3C�ٚC���C�� C��3C��fC���C���C��3C��fC�ٚC���C�� C��3C���C���C���C��3C��fC�ٚC���C�� C��fC���C���C���C�� C��3C��fC��fC�ٚC�ٚC���C���C�� C��3C��3C��fC���C���C���C��3C��3C��fC�ٚC�ٚC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C�� Cͳ3C���C���C���C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C���C���C�ٚCަfCߦfC�3C�� C���C�3C�� C�ٚC�3C�� C�ٚC�3C�� C���C�fC��3C�� C���C�ٚC�3C�3C���C��fC�� C�ٚC��3C���C��3C�� C��3D 9�D� D��D��D,�D` D��D	�D
FfD� D� D  D9�D� D�fD��D33Ds3D��D��DFfD��D� D��D9�D��D � D!��D#9�D$s3D%�3D&��D(,�D)s3D*� D+�3D-  D.s3D/�fD1  D29�D3s3D4�3D5��D733D8s3D9��D;  D<L�D=y�D>� D?��DA9�DB�fDC��DD�3DF&fDG` DH� DI�fDK,�DLs3DM� DOfDP9�DQl�DR� DS�3DU,�DV�fDW�fDX��DZ33D[y�D\� D^�D_@ D`s3Da�fDb��Dd33De� Df�fDg�3Di&fDjy�Dk��DmfDnFfDo�fDp�fDrfDsL�Dts3Du� Dv��Dx@ Dys3Dz�fD{� D}  D~` D��D�y�D��D��3D�Y�D��fD���D�9�D��fD��3D�  D�� D�` D��3D���D�@ D�ٚD�s3D��D�ɚD�ffD�3D���D�9�D�ٚD�vfD�fD��fD�VfD���D�� D�FfD���D�s3D��D�� D�ffD���D���D�9�D�ٚD�y�D��D���D�` D�fD���D�6fD���D��fD�#3D�� D�` D���D���D�<�D�� D��3D�&fD���D�S3D���D��3D�<�D��fD�s3D�  D�ɚD�ffD�fD��fD�FfD��D���D�,�D�� D�S3D���D�� D�FfD���D�p D��D���D�S3D��fD�� D�9�D��3D�|�D��D���D�` D��3D���D�@ D��fD�� D�)�D��3D�` D���D��fD�33D�� D�p D� D° D�P D�� Dē3D�33D���D�|�D�fDǳ3D�P D�� Dɓ3D�6fD���DˆfD�  D̹�D�` D���DΓ3D�<�D�ٚD�vfD�  D�� D�` D�  DӖfD�9�D���DՀ D�&fD�� D�Y�D�  Dة�D�C3D��3Dڀ D�  D�� D�c3D�	�Dݜ�D�@ D��fD߃3D�&fD��D�VfD���D�fD�C3D�� D� D�  D�� D�\�D���D��D�<�D�� D�fD�  D�fD�\�D�fD� D�9�D��fD�s3D�3D� D�\�D�	�D�D�I�D���D� D�3D�� D�` D�  D���D�<�D���D�� D�#3D���D�S3D�  D�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����𿆧��T������j��Mӿ|�z�H�x�u�j=q�gl��i7L�d�/�\��G+�)��)��#����1>]/>�-?&ff?^��?\(�?���?�j@@�-@&�R@Up�@dz�@l�@uO�@t�@z�\@t�D@s�F@oK�@z��@x�u@v�@t��@r~�@k@hbN@h�u@_
=@aX@i&�@m�@j�@l�@g�@d�j@_��@W+@Up�@P�u@L�D@J�@FE�@D�/@C�
@CS�@C@B�@Bn�@@r�@>E�@=�@=��@;S�@;o@;@9�^@8b@7l�@5p�@3��@2=q@1X@1�@0b@/+@.ȴ@.ff@-V@,I�@+�
@)�@(r�@(bN@(A�@(b@&v�@%�-@%V@$I�@!X@{@n�@�@?}@�@�j@��@n�@+@�F@	��@	�7@	x�@	7L@��@E�@��@=q?���?��?���?�A�?�w?��-?�1?��?��?�$�?�\)?�C�?֧�?�n�?θR?�(�?�C�?���?�E�?��/?���?�"�?��?���?��?��?���?�j?���?�X?�K�?�?��?��?�;d?�p�?��?�l�?���?�ff?��T?��?��?�-?�G�?�;d?�I�?�X?�`B?�33?�&�?~�R?y�#?w��?u?r-?p��?nV?k?e`B?c�
?co?co?^��?["�?U�?S33?S33?S�F?R-?O�;?N{?MV?KC�?I7L?G�?F�y?F$�?D��?DZ?BM�??;d?=�-?;�m?;dZ?9X?8�u?5�?49X?2�!?0��?/\)?.V?.{?-�h?,�D?,�D?-��?.{?-V?+�?(�9?$��?!��?|�?�?/?�?��?��?��?�#?�u?�P?
=?E�?��?
=??}?z�?�F?��?dZ?dZ?�#?��?
=?E�?�?
~�?$�?J>��m>�K�>�F>�D>��>��T>�`B>�S�>ܬ>ٙ�>׍P>�z�>��>�=q>��>��>��>���>�M�>�;d>��>���>���>�>��>�V>��>�=q>�1'>u>p��>m�h>gl�>aG�>\(�>V>Kƨ>I�^>H�9>D��>=p�>2->/�>(��>%�T>�>hs>$�=��#=�;d=�
==��=�v�=��w=��P=�+=q��=<j=\)<�`B<ě�<�t�<t�;�`B;��
    �ě��T����j���,1�0 Ž49X�<j�H�9��%��o�u�q���q���}󶽁%�����C���+������罺^5�Ƨ����l��o�+�	7L�C��V�n���+�����R�#�
�+�/��333�:^5�?|�B�\�D���E�˾J���L�;N��O�;�Q녾T���Xb�Z��]/�_;d�^5?�]/�\(��]/�\(��^5?�bMӾdZ�fff�gl��hr��j~��o���u�x���z�H�}󶾀  ���\���˾�+���9������I���O߾�����`��녾��Ͼ���
=���P���P�����������"Ѿ��㾜���5?���w��G���S���Z��`B��`B���T��ff��r����þ�~���1��V������ ž��׾�&龳33��9X����ȴ��Q쾹X��^5��푾��۾��7��J�Õ��š˾ɺ^��=q��C����;�V���;��hs���Ͼ����Ձ�׍P����ۥ��/�޸R��Ĝ�����`B���y����þ����1��h�����׾����F����KǾ�����^5��p���vɾ�vɾ�|�   �%��7�Mӿ��S����������`B�$ݿ�y��y�l��r��	7L�	�^�	��
=q��ƨ�I���ͿO߿��V��� ſ ſ ſbN��`�hs�녿n��n���t���Ͽz����?}���E���+�
=�
=�ȴ�ȴ��ٿ�u�X�������H�"ѿ�m���푿/�p���vɿvɿ�R��ۿ|�   � A�� A�� ��!%�!G��"J�"Mӿ"Mӿ#o�#S��#S��#S��#�
�$Z�%��%�T�&ff�'l��'��(1'�(r��(r��(�ÿ)7L�)7L�)�^�*=q�*���+�+C��+ƨ�,1�,I��,�D�-V�-V�-�h�-��-��.{�.{�.���.��.��.��/\)�/���/���/�;�0 ſ0 ſ0 ſ0bN�0bN111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ���𿆧��T������j��Mӿ|�z�H�x�u�j=q�gl��i7L�d�/�\��G+�)��)��#����1>]/>�-?&ff?^��?\(�?���?�j@@�-@&�R@Up�@dz�@l�@uO�@t�@z�\@t�D@s�F@oK�@z��@x�u@v�@t��@r~�@k@hbN@h�u@_
=@aX@i&�@m�@j�@l�@g�@d�j@_��@W+@Up�@P�u@L�D@J�@FE�@D�/@C�
@CS�@C@B�@Bn�@@r�@>E�@=�@=��@;S�@;o@;@9�^@8b@7l�@5p�@3��@2=q@1X@1�@0b@/+@.ȴ@.ff@-V@,I�@+�
@)�@(r�@(bN@(A�@(b@&v�@%�-@%V@$I�@!X@{@n�@�@?}@�@�j@��@n�@+@�F@	��@	�7@	x�@	7L@��@E�@��@=q?���?��?���?�A�?�w?��-?�1?��?��?�$�?�\)?�C�?֧�?�n�?θR?�(�?�C�?���?�E�?��/?���?�"�?��?���?��?��?���?�j?���?�X?�K�?�?��?��?�;d?�p�?��?�l�?���?�ff?��T?��?��?�-?�G�?�;d?�I�?�X?�`B?�33?�&�?~�R?y�#?w��?u?r-?p��?nV?k?e`B?c�
?co?co?^��?["�?U�?S33?S33?S�F?R-?O�;?N{?MV?KC�?I7L?G�?F�y?F$�?D��?DZ?BM�??;d?=�-?;�m?;dZ?9X?8�u?5�?49X?2�!?0��?/\)?.V?.{?-�h?,�D?,�D?-��?.{?-V?+�?(�9?$��?!��?|�?�?/?�?��?��?��?�#?�u?�P?
=?E�?��?
=??}?z�?�F?��?dZ?dZ?�#?��?
=?E�?�?
~�?$�?J>��m>�K�>�F>�D>��>��T>�`B>�S�>ܬ>ٙ�>׍P>�z�>��>�=q>��>��>��>���>�M�>�;d>��>���>���>�>��>�V>��>�=q>�1'>u>p��>m�h>gl�>aG�>\(�>V>Kƨ>I�^>H�9>D��>=p�>2->/�>(��>%�T>�>hs>$�=��#=�;d=�
==��=�v�=��w=��P=�+=q��=<j=\)<�`B<ě�<�t�<t�;�`B;��
    �ě��T����j���,1�0 Ž49X�<j�H�9��%��o�u�q���q���}󶽁%�����C���+������罺^5�Ƨ����l��o�+�	7L�C��V�n���+�����R�#�
�+�/��333�:^5�?|�B�\�D���E�˾J���L�;N��O�;�Q녾T���Xb�Z��]/�_;d�^5?�]/�\(��]/�\(��^5?�bMӾdZ�fff�gl��hr��j~��o���u�x���z�H�}󶾀  ���\���˾�+���9������I���O߾�����`��녾��Ͼ���
=���P���P�����������"Ѿ��㾜���5?���w��G���S���Z��`B��`B���T��ff��r����þ�~���1��V������ ž��׾�&龳33��9X����ȴ��Q쾹X��^5��푾��۾��7��J�Õ��š˾ɺ^��=q��C����;�V���;��hs���Ͼ����Ձ�׍P����ۥ��/�޸R��Ĝ�����`B���y����þ����1��h�����׾����F����KǾ�����^5��p���vɾ�vɾ�|�   �%��7�Mӿ��S����������`B�$ݿ�y��y�l��r��	7L�	�^�	��
=q��ƨ�I���ͿO߿��V��� ſ ſ ſbN��`�hs�녿n��n���t���Ͽz����?}���E���+�
=�
=�ȴ�ȴ��ٿ�u�X�������H�"ѿ�m���푿/�p���vɿvɿ�R��ۿ|�   � A�� A�� ��!%�!G��"J�"Mӿ"Mӿ#o�#S��#S��#S��#�
�$Z�%��%�T�&ff�'l��'��(1'�(r��(r��(�ÿ)7L�)7L�)�^�*=q�*���+�+C��+ƨ�,1�,I��,�D�-V�-V�-�h�-��-��.{�.{�.���.��.��.��/\)�/���/���/�;�0 ſ0 ſ0 ſ0bN�0bN111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�qB�?B��B�
B�TB��B\B$�B49BF�BL�BL�BJ�BYB`BB�JB�7B�bB��B�yB�yB  BB�B`BB�JB��B	>wB	�B	��B
H�B
� B
��B
��B
��B
�
B
�B�B�BR�BQ�BP�BO�BR�BN�BN�BO�Bo�Be`B�B�bB��B��B�dB�^BȴB��B��B�TB�sB�NB�B�B�B�B�B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�sB�`B�TB�;B�TB�BB�BB�NB�HB�BB�/B�#B�#B�)B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBƨBÖBÖBÖBÖB�}B�wB�}B�wB�qB�qB�qB�qB�qB�qB�qB�jB�dB�jB�^B�XB�^B�XB�XB�XB�XB�XB�XB�RB�RB�LB�FB�?B�?B�3B�9B�9B�3B�3B�3B�3B�3B�9B�3B�-B�3B�3B�-B�9B�'B�'B�-B�3B�9B�3B�3B�3B�3B�3B�9B�3B�9B�9B�3B�3B�3B�9B�3B�3B�3B�3B�3B�-B�-B�3B�-B�-B�9B�9B�3B�?B�FB�FB�LB�?B�9B�9B�3B�?B�?B�9B�9B�?B�?B�?B�?B�?B�?B�FB�?B�LB�LB�LB�LB�LB�RB�^B�XB�^B�XB�RB�RB�RB�3B�-B�3B�-B�-B�'B�'B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�qB�?B��B�
B�TB��B\B$�B49BF�BL�BL�BJ�BYB`BB�JB�7B�bB��B�yB�yB  BB�B`BB�JB��B	>wB	�B	��B
H�B
� B
��B
��B
��B
�
B
�B�B�BR�BQ�BP�BO�BR�BN�BN�BO�Bo�Be`B�B�bB��B��B�dB�^BȴB��B��B�TB�sB�NB�B�B�B�B�B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�sB�`B�TB�;B�TB�BB�BB�NB�HB�BB�/B�#B�#B�)B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBƨBÖBÖBÖBÖB�}B�wB�}B�wB�qB�qB�qB�qB�qB�qB�qB�jB�dB�jB�^B�XB�^B�XB�XB�XB�XB�XB�XB�RB�RB�LB�FB�?B�?B�3B�9B�9B�3B�3B�3B�3B�3B�9B�3B�-B�3B�3B�-B�9B�'B�'B�-B�3B�9B�3B�3B�3B�3B�3B�9B�3B�9B�9B�3B�3B�3B�9B�3B�3B�3B�3B�3B�-B�-B�3B�-B�-B�9B�9B�3B�?B�FB�FB�LB�?B�9B�9B�3B�?B�?B�9B�9B�?B�?B�?B�?B�?B�?B�FB�?B�LB�LB�LB�LB�LB�RB�^B�XB�^B�XB�RB�RB�RB�3B�-B�3B�-B�-B�'B�'B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092229                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092319  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092319  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            @���D�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172540  IP  PSAL            @���D�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            @���D�� G�O�                