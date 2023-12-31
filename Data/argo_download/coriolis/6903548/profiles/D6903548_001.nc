CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  O   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:48Z creation; 2021-06-07T15:42:38Z last update (BSH ARSQ software)    
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
_FillValue                  ,  �Argo profile    3.1 1.2 19500101000000  20190620115648  20210607154238  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @�qv�O��1   @�qw�u0�@TP}rc;�@-ٜ û�1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 9.3 dbar]                                                      AffA!��A0  A@  AQ��A^ffAl��A�  A�  A�33A���A���A�  A�  A�33A�ffA���A�ffA�33A�  A���A�  A�ffB   B��BffBffBffB  B33B33B��B$��B(ffB+��B133B4��B7��B;33B>��BDffBHffBM33BP  BTffBXffB\ffB`ffBd  Bg��BlffBpffBtffBx  B|  B�  B�  B�  B�  B�  B�  B�  B�33B���B���B�  B���B���B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B���B���B���B���B���B���B�  B�  B�  B�ffB���BÙ�B���B�  B���B�  B�33B���B�ffB���B�ffB�  B���BB�ffB�  B�  B�  C� CffCffCffC	ffCffCffCffC� C� C��C33CL�CffC� C��C!33C#L�C%� C'� C)�3C+ffC-ffC/��C133C3ffC5� C7L�C9� C;��C=L�C?� CA��CCL�CE� CG�3CIffCK� CM��CO��CQffCS33CUffCW� CYffC[33C]� C_�3CaffCc33CeffCg��CiffCk33CmffCo��CqffCs33CuffCw��CyffC{33C}� C�3C�� C��fC���C��3C���C�� C��fC���C��3C�ٚC���C�� C��3C���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��3C��fC�ٚC�ٚC���C���C�� C�� C�� C��3C��3C��3C��fC��fC��fC��fC���C���C���C���C���C���C���C���C���C���C��fC��fC��fC��fC��3C��3C�� C���C���C�ٚC��fC�� C�� C�� C¦fC���C���C���C���C���C���C���C�ٚC˦fC̦fCͳ3C�� C�ٚCЦfCѳ3C���CӦfC�� C�ٚCֳ3C׌�CئfC�� CڦfCی�CܦfC���C޳3CߦfC���C�3C��fC���C�� C�3C�fC癚C��C��C� C� C�3C�ٚC���C�� C�3C�fC�C�C��C���C�� C�� C��3C��3C��fC��D @ Ds3D� D�D9�Dl�D�3D	  D
33DffD�fD��DS3D��D� D��D9�Dy�D� DfDS3D� D�fD��D,�Ds3D � D"�D#9�D$l�D%��D&�3D(9�D)� D*��D+�3D-33D.s3D/��D1  D2L�D3� D4� D6  D79�D8y�D9� D;�D<@ D=s3D>�fD@  DAS3DB�fDC� DD�3DF33DGs3DH��DJ  DKFfDL�3DM� DN��DP&fDQffDR� DS� DU  DVffDW��DX��DZ@ D[��D\��D]��D_9�D`� Da�fDb��Dd&fDel�Df��DhfDi33DjffDk�3Dl��DnFfDo� Dp�3Dq��Ds9�Dt� Du�fDw3Dx@ Dyy�Dz�3D{��D},�D~s3D��D��3D��D��3D�Y�D��fD��3D�6fD��3D�vfD��D��3D�\�D��fD�� D�9�D��fD�vfD�fD���D�\�D�fD�� D�9�D��3D�p D��D���D�i�D���D�� D�9�D�� D��fD�  D���D�` D���D��3D�<�D��3D�|�D�fD��3D�\�D�	�D���D�33D��fD�y�D�  D���D�VfD�  D�� D�33D�ٚD�� D�fD�� D�i�D�  D���D�@ D��D��3D��D���D�VfD��3D�� D�<�D���D���D�,�D�ɚD�i�D�	�D���D�I�D���D�p D�3D���D�` D���D��3D�9�D�� D��fD�  D���D�VfD��3D�� D�0 D��3D�vfD��D���D�\�D�  D��fD�<�D��fD�� D�  D¼�D�\�D���DĜ�D�@ D��fD�|�D�fD�� D�c3D���Dɠ D�FfD�� D�|�D�fD̰ D�Y�D�fDΣ3D�@ D�� DЀ D�  D�� D�` D�  DӠ D�<�D�ٚD�|�D�  D��fD�c3D�  Dؠ D�C3D�ٚD�vfD��D��3D�\�D��fDݜ�D�9�D��fD߀ D�  D��3D�c3D���D��D�<�D�� D�3D�  D�� D�c3D�fD�3D�C3D��fD�3D�  D�� D�c3D���D왚D�6fD��3D�p D� D� D�P D��fD��D�C3D���D�y�D�#3D�� D�\�D�  D���D�@ D���D���D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 AffA!��A0  A@  AQ��A^ffAl��A�  A�  A�33A���A���A�  A�  A�33A�ffA���A�ffA�33A�  A���A�  A�ffB   B��BffBffBffB  B33B33B��B$��B(ffB+��B133B4��B7��B;33B>��BDffBHffBM33BP  BTffBXffB\ffB`ffBd  Bg��BlffBpffBtffBx  B|  B�  B�  B�  B�  B�  B�  B�  B�33B���B���B�  B���B���B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B���B���B���B���B���B���B�  B�  B�  B�ffB���BÙ�B���B�  B���B�  B�33B���B�ffB���B�ffB�  B���BB�ffB�  B�  B�  C� CffCffCffC	ffCffCffCffC� C� C��C33CL�CffC� C��C!33C#L�C%� C'� C)�3C+ffC-ffC/��C133C3ffC5� C7L�C9� C;��C=L�C?� CA��CCL�CE� CG�3CIffCK� CM��CO��CQffCS33CUffCW� CYffC[33C]� C_�3CaffCc33CeffCg��CiffCk33CmffCo��CqffCs33CuffCw��CyffC{33C}� C�3C�� C��fC���C��3C���C�� C��fC���C��3C�ٚC���C�� C��3C���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��3C��fC�ٚC�ٚC���C���C�� C�� C�� C��3C��3C��3C��fC��fC��fC��fC���C���C���C���C���C���C���C���C���C���C��fC��fC��fC��fC��3C��3C�� C���C���C�ٚC��fC�� C�� C�� C¦fC���C���C���C���C���C���C���C�ٚC˦fC̦fCͳ3C�� C�ٚCЦfCѳ3C���CӦfC�� C�ٚCֳ3C׌�CئfC�� CڦfCی�CܦfC���C޳3CߦfC���C�3C��fC���C�� C�3C�fC癚C��C��C� C� C�3C�ٚC���C�� C�3C�fC�C�C��C���C�� C�� C��3C��3C��fC��D @ Ds3D� D�D9�Dl�D�3D	  D
33DffD�fD��DS3D��D� D��D9�Dy�D� DfDS3D� D�fD��D,�Ds3D � D"�D#9�D$l�D%��D&�3D(9�D)� D*��D+�3D-33D.s3D/��D1  D2L�D3� D4� D6  D79�D8y�D9� D;�D<@ D=s3D>�fD@  DAS3DB�fDC� DD�3DF33DGs3DH��DJ  DKFfDL�3DM� DN��DP&fDQffDR� DS� DU  DVffDW��DX��DZ@ D[��D\��D]��D_9�D`� Da�fDb��Dd&fDel�Df��DhfDi33DjffDk�3Dl��DnFfDo� Dp�3Dq��Ds9�Dt� Du�fDw3Dx@ Dyy�Dz�3D{��D},�D~s3D��D��3D��D��3D�Y�D��fD��3D�6fD��3D�vfD��D��3D�\�D��fD�� D�9�D��fD�vfD�fD���D�\�D�fD�� D�9�D��3D�p D��D���D�i�D���D�� D�9�D�� D��fD�  D���D�` D���D��3D�<�D��3D�|�D�fD��3D�\�D�	�D���D�33D��fD�y�D�  D���D�VfD�  D�� D�33D�ٚD�� D�fD�� D�i�D�  D���D�@ D��D��3D��D���D�VfD��3D�� D�<�D���D���D�,�D�ɚD�i�D�	�D���D�I�D���D�p D�3D���D�` D���D��3D�9�D�� D��fD�  D���D�VfD��3D�� D�0 D��3D�vfD��D���D�\�D�  D��fD�<�D��fD�� D�  D¼�D�\�D���DĜ�D�@ D��fD�|�D�fD�� D�c3D���Dɠ D�FfD�� D�|�D�fD̰ D�Y�D�fDΣ3D�@ D�� DЀ D�  D�� D�` D�  DӠ D�<�D�ٚD�|�D�  D��fD�c3D�  Dؠ D�C3D�ٚD�vfD��D��3D�\�D��fDݜ�D�9�D��fD߀ D�  D��3D�c3D���D��D�<�D�� D�3D�  D�� D�c3D�fD�3D�C3D��fD�3D�  D�� D�c3D���D왚D�6fD��3D�p D� D� D�P D��fD��D�C3D���D�y�D�#3D�� D�\�D�  D���D�@ D���D���D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?���?���?ܬ?ؓu?�X?�=q?�;d?��?�\)?���?�ff?��?��?ꟾ?��?�%?�z�?�{?�Z?ļj?ۅ?�j?��?�?��
?�&�?��?�7L?� �?Ͳ-?͑h?�F?�1?�^?��m?�$�?���?�G�?֧�?�bN@ Ĝ?�1'?��+@`B@��@(Q�@/�;@3o@7\)@6$�@/��@ ��@(�@33@�@��@dZ@C�@1'@E�@�@S�?��h?�O�@-@�@hs@ �u@�@��@�`@(�@Q�@  @@�/@�@��@��@ ��@ Q�@C�@E�@+@�`@��@r�@Q�@b@ �@	X@��@�`@  @S�@J@=q@�!@=q@-@-@-@�\@$�@�y@��@�w@1'@�`@
=q@�m@
��@o@I�@O�@p�@�T@��@$�@v�@V@�@|�@�@�u@X@M�@S�@dZ@��@M�@^5@�!@��@��@t�@�F@�m@��@�j@��@z�@z�@j@Z@�
@"�@�!@J@%@�@l�@�y@@�-@(�@dZ@
�@	��@	�@	&�@	&�@	�7@
��@
M�@	�^@	��@	��@	��@
J@�`@Q�@Q�@	X@�`@��@
J@
M�@
J@
J@	�#@
M�@�@`B@�h@`B@O�@�@?}@�@�@?}@`B@O�@��@�@��@��@�@�/@�@��@�j@��@�@V@�@V@�@��@Z@Z@9X@1@1@1@�m@��@t�@��@��@dZ@t�@�F@9X@��@�j@�j@�/@�/@��@��@�j@��@�j@��@z�@��@��@�D@Z@�
@C�@"�@
�!@
^5@
-@
J@	��@�`@�@ �@  @��@K�@v�@$�@�T@?}@�@9X@��@�!@�H@%@ �@ A�?��;?�;d?��?��?�dZ?��?���?�ȴ?��+?�E�?�?�?��?�Z?�?�!?�J?���?��?���?�^5?陚?�?�J?��;?�v�?�/?�7L?���?�O�?�^5?ȓu?�1'?�+?���?��?�E�?�o?�{?�E�?���?�1?�b?�Z?�&�?�{?��?�-?}�-?q��?k�?c��?X��?I7L?F��?B��?5?(��?!�7?"�?��?�!?C�?�/?   >�E�>�`B>�n�>ȴ9>�ȴ>�&�>�{>�>�
=>�C�>���>vȴ>k�>]/>F��>1&�>�R>O�=��=�`B=���=��T=�7L=q��=0 �<���<���<D���D���T����t���h�8Q�T���u���w��vɽ�;d���پC����'-V�<j�I�^�P�`�Y��dZ�l�D�vȴ��  �����+��ƨ�����t���
=���-���R��G���r���~������h��������Q쾻�m��p����7��o��+�Ǯ�ɺ^��=q������O߾�bN��t�����z�׍P�ؓu��"Ѿݲ-�߾w��A������`B���T��ff��r���~���1��V��{�� ž�!��?}��ȴ��Q������X��dZ��j��푾�vɿ A��%��7���Z����T����9�	7L�	�^�	��	��
=q�
~��ƨ�����\)�����`�&�����!��F��Ͽ9X����E��ȴ�b�����������"ѿ"ѿdZ�����푿�vɿ�R�|�|��w��w� A��!%�!���"Mӿ#S��#���$�/�%�˿&$ݿ&�y�'��(�9�(�ÿ)�^�)7L�)��*���,�Ϳ-V�.���/��/���/�;�0 ſ0�`�1&�1hs�1���1녿1녿1녿2n��3�F�49X�4z�4z�49X�3�Ͽ49X�49X�4���5�6�+�6ȴ�7
=�7Kǿ7�P�7�ٿ8Q�9��9X�9X�9���:��:��:��:�H�;��;�m�;�m�<(��<(��<j�<j�<j�<j�<j�=/�=�-�=�=�>5?�>�ۿ?;d�?;d�?;d�?|�?�w�?�w�?�w�@��@Ĝ�A%�A%�AG��A�7�AG��A�7�A�7�A�7�A�7�A�7�A���A���BJ�BJ�BJ�BJ�BJ�BJ�BJ�BJ111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ?���?���?ܬ?ؓu?�X?�=q?�;d?��?�\)?���?�ff?��?��?ꟾ?��?�%?�z�?�{?�Z?ļj?ۅ?�j?��?�?��
?�&�?��?�7L?� �?Ͳ-?͑h?�F?�1?�^?��m?�$�?���?�G�G�O�G�O�@ Ĝ?�1'?��+@`B@��@(Q�@/�;@3o@7\)@6$�@/��@ ��@(�@33@�@��@dZ@C�@1'@E�@�@S�?��h?�O�@-@�@hs@ �u@�@��@�`@(�@Q�@  @@�/@�@��@��@ ��@ Q�@C�@E�@+@�`@��@r�@Q�@b@ �@	X@��@�`@  @S�@J@=q@�!@=q@-@-@-@�\@$�@�y@��@�w@1'@�`@
=q@�m@
��@o@I�@O�@p�@�T@��@$�@v�@V@�@|�@�@�u@X@M�@S�@dZ@��@M�@^5@�!@��@��@t�@�F@�m@��@�j@��@z�@z�@j@Z@�
@"�@�!@J@%@�@l�@�y@@�-@(�@dZ@
�@	��@	�@	&�@	&�@	�7@
��@
M�@	�^@	��@	��@	��@
J@�`@Q�@Q�@	X@�`@��@
J@
M�@
J@
J@	�#@
M�@�@`B@�h@`B@O�@�@?}@�@�@?}@`B@O�@��@�@��@��@�@�/@�@��@�j@��@�@V@�@V@�@��@Z@Z@9X@1@1@1@�m@��@t�@��@��@dZ@t�@�F@9X@��@�j@�j@�/@�/@��@��@�j@��@�j@��@z�@��@��@�D@Z@�
@C�@"�@
�!@
^5@
-@
J@	��@�`@�@ �@  @��@K�@v�@$�@�T@?}@�@9X@��@�!@�H@%@ �@ A�?��;?�;d?��?��?�dZ?��?���?�ȴ?��+?�E�?�?�?��?�Z?�?�!?�J?���?��?���?�^5?陚?�?�J?��;?�v�?�/?�7L?���?�O�?�^5?ȓu?�1'?�+?���?��?�E�?�o?�{?�E�?���?�1?�b?�Z?�&�?�{?��?�-?}�-?q��?k�?c��?X��?I7L?F��?B��?5?(��?!�7?"�?��?�!?C�?�/?   >�E�>�`B>�n�>ȴ9>�ȴ>�&�>�{>�>�
=>�C�>���>vȴ>k�>]/>F��>1&�>�R>O�=��=�`B=���=��T=�7L=q��=0 �<���<���<D���D���T����t���h�8Q�T���u���w��vɽ�;d���پC����'-V�<j�I�^�P�`�Y��dZ�l�D�vȴ��  �����+��ƨ�����t���
=���-���R��G���r���~������h��������Q쾻�m��p����7��o��+�Ǯ�ɺ^��=q������O߾�bN��t�����z�׍P�ؓu��"Ѿݲ-�߾w��A������`B���T��ff��r���~���1��V��{�� ž�!��?}��ȴ��Q������X��dZ��j��푾�vɿ A��%��7���Z����T����9�	7L�	�^�	��	��
=q�
~��ƨ�����\)�����`�&�����!��F��Ͽ9X����E��ȴ�b�����������"ѿ"ѿdZ�����푿�vɿ�R�|�|��w��w� A��!%�!���"Mӿ#S��#���$�/�%�˿&$ݿ&�y�'��(�9�(�ÿ)�^�)7L�)��*���,�Ϳ-V�.���/��/���/�;�0 ſ0�`�1&�1hs�1���1녿1녿1녿2n��3�F�49X�4z�4z�49X�3�Ͽ49X�49X�4���5�6�+�6ȴ�7
=�7Kǿ7�P�7�ٿ8Q�9��9X�9X�9���:��:��:��:�H�;��;�m�;�m�<(��<(��<j�<j�<j�<j�<j�=/�=�-�=�=�>5?�>�ۿ?;d�?;d�?;d�?|�?�w�?�w�?�w�@��@Ĝ�A%�A%�AG��A�7�AG��A�7�A�7�A�7�A�7�A�7�A���A���BJ�BJ�BJ�BJ�BJ�BJ�BJ�BJ111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B��B:^B�LB�B7LB��B�
B�;B�sB�B	�B	'�B	P�B	]/B	o�B	�PB	~�B	�hB	�FB	�qB	�B	�mB
uB
�B
#�B
�B
0!B
<jB
?}B
^5B
�+B
�uB
��B
�B
�B
�!B
�'B
�!B
�DB
�B
�;B
��BPB:^B]/Bn�B[#B�Bx�Bm�BbNBJ�BR�BQ�BVBW
BI�BffBe`BbNBu�B]/B^5BcTBdZBffBhsBjB\)Bz�Br�B{�B�B{�B|�B}�B�B|�B}�B� B� B�JB�VB�uB�{B�{B��B�{B��B��B��B��B��B�uB�uB�hB��B�{B�{B��B��B��B��B��B��B��B��B��B�B�B�B�-B�-B�?B�FB�LB�RB�RB�XB�^B�jB�qB�wB�wBBŢBŢBŢBŢBŢBŢBŢBǮBǮBǮBɺB��B��B��B��B��B��B��B��B��B��BɺBɺBǮBǮBǮBĜBƨBBÖB��B�wB�wB�qB�qB�wB�^B��B��B��B�}B��B��B��B��B�qB��B�jBBBÖBBBBBĜBȴBɺBɺBɺBɺBȴBɺBɺBɺB��BɺBɺB��BɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȴBǮBǮBŢBŢBĜBĜBÖBÖB��B��B�}B�}B�wB��B�}B�wB�qB�jB�dB�^B�XB�LB�LB�?B�?B�9B�3B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 BݩB�B=�B��B�0B:�B�BږB��B��B�B	B	+|B	TqB	`�B	s*B	��B	��B	��B	��B	��B	�B	��B
B
B
'cB
#JB
3�B
?�B
C	B
a�B
��B
�B
�uB
��B
��B
��B
��G�O�G�O�B
ِB
��B
�MB�B=�B`�Br$B^�B��B|aBqBe�BNMBV~BUxBY�BZ�BMFBi�Bh�Be�ByOB`�Ba�Bf�Bg�Bi�Bk�BnB_�B~mBv<BsB��BsB�zB��B��B�zB��B��B��B��B��B�B�B�B�B�B�B�B�B�B�,B�B�B��B�B�B�B�B�B�B�B�QB�QB�]B�iB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�.B�.B�.B�.B�.B�.B�.B�:B�:B�:B�FB�MB�YB�YB�YB�YB�SB�YB�SB�SB�SB�FB�FB�:B�:B�:B�(B�4B�B�"B�B�B�B��B��B�B��B�B�B�B�	B�B�B�B�B��B�B��B�B�B�"B�B�B�B�B�(B�@B�FB�FB�FB�FB�@B�FB�FB�FB�MB�FB�FB�MB�FB�FB�MB�MB�MB�MB�MB�MB�MB�MB�SB�SB�SB�SB�SB�MB�MB�SB�SB�MB�MB�MB�MB�MB�MB�SB�MB�SB�SB�_B�_B�eB�eB�eB�eB�eB�kB�eB�kB�kB�kB�kB�qB�qB�qB�qB�qB�kB�kB�qB�kB�eB�eB�eB�eB�_B�eB�eB�eB�eB�_B�_B�_B�_B�eB�eB�SB�_B�MB�_B�eB�_B�YB�_B�YB�YB�_B�YB�SB�SB�YB�SB�SB�SB�YB�YB�YB�_B�YB�SB�MB�@B�:B�:B�.B�.B�(B�(B�"B�"B�B�B�	B�	B�B�B�	B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�uB�oB�oB�iB�iB�iB�iB�iB�iB�iB�cB�iB�cB�]B�]B�WB�WB�WB�QB�JB�JB�JB�JB�JB�JB�JB�DB�DB�JB�>B�>B�>B�>B�>B�>B�>B�DB�>B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�>B�>B�>B�8B�8B�8B�8B�8B�8B�8B�8B�>B�>B�8B�8B�>B�>B�>B�>B�>B�>B�>B�DB�DB�DB�DB�DB�JB�JB�JB�QB�QB�JB�QB�JB�QB�QB�QB�QB�QB�QB�JB�QB�QB�QB�QB�QB�QB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�]B�]B�]B�]B�]B�cB�cB�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�cB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�oB�oB�iB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�|B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542382021060715423820210607154238  IF  ARFMCODA029d                                                                20190620115648                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115701  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115701  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154238  IP  PSAL            AffD��fG�O�                