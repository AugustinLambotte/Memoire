CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   =   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2013-04-08T10:43:04Z creation; 2019-06-03T17:52:39Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7(   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    78   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7<   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7@   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7P   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7`   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7p   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7x   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8(   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8,   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    80   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     84   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8T   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8X   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8\   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8|   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T      
resolution        >�EȠ�Q)        8�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�EȠ�Q)        8�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z         �  :    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  :�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z         �  ;4   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  <(   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���      �  <h   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  =\   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  >P   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  >�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  ?�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  ?�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  @�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  A�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  A�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  @  B�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  C    	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  D   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    DD   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    GD   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    JD   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  MD   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    Mp   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    Mt   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    Mx   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    M|   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  M�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    M�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    M�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    M�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         M�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         M�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        M�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    M�Argo profile    3.1 1.2 19500101000000  20130408104304  20190603175239  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL               �A   IF  29143153                        2C  D   APEX                            2573                            n/a                             846 @֍N�6P1   @֍P<���@Q��t�j�#yXbM�1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   A   A   @�ffAffAh  A�  A���A���BffB��B1��BE��BY��Bn  B�33B���B�ffB�  B�33B�  B�  Bƙ�B�  B�ffCL�C�CffC33C)L�C3��C=� CG33CQL�C[� Ce� Co��Cy33C��3C���C���C�� C���C��3C��fC��fC��3C��3C�s3C���C�� C CǙ�C�� CѦfCֳ3Cۀ C�s3C� C�fC�s3C�ffC���C�ٚ2111211222211222222222222222222111211221111222222222222222211   ����@���A  Ah  A���A���A���B	��B��B1��BE��BZ  BlffB���B�ffB�  B�33B�  B�  B���B�  B�ffB���C�CffC33C$L�C.��C8� CB33CLL�CV� C`� Cj��Ct33C~ffC��C��C�@ C�L�C�33C�&fC�&fC�33C�33C��3C��C�  C�  C��C�@ C�&fC�33C�  C��3C�  C�&fC��3C��fC�L�C�Y�1111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B
�B
�B
�B
�B
�B
�B
�B
�B
�B
��B
��B
��B
��B
��B
��B
��BBB+B
=B �B7LBK�BcTBn�Bo�Bq�Bq�Bp�Bq�Br�Bt�Bt�Bt�Bs�Bt�Bt�Bu�Bw�By�By�Bz�B� B~�B�B�B�=B�DB�B�+B�7B�PB�\B�hB�oB�uB�uB�uB�{B��B�{1111111111111111111111111111111111111111111111111111111111111   B
�B
�B
�B
�B
�B
�B
�B
�B
�B
��B
��B
��B
��B
��B
��B
��BBB+B
=B �B7LBK�BcTBn�Bo�Bq�Bq�Bp�Bq�Br�Bt�Bt�Bt�Bs�Bt�Bt�Bu�Bw�By�By�Bz�B� B~�B�B�B�=B�DB�B�+B�7B�PB�\B�hB�oB�uB�uB�uB�{B��B�{1111111111111111111111111111111111111111111111111111111111111   <��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
�� ſz�H�x�u�v�+�t���r-�kC��k��i��e��bMӿ`��_�w�_;d�["ѿVȴ�Rn��N{�J=q�G+�1���녾�
=��ƨ�5?}�8Q�'(�þ>vɾ:^5�333�%�T�����P�)��)���-���   ��l���;d��;d�'�O߼���D��=L��<�1��C��<j�#�
;ě�<��
<�9X<�/<�C�<T��;ě�:�o�#�
��`B1111111111111111111111111111111111111111111111111111111111111   �� ſz�H�x�u�v�+�t���r-�kC��k��i��e��bMӿ`��_�w�_;d�["ѿVȴ�Rn��N{�J=q�G+�1���녾�
=��ƨ�5?}�8Q�'(�þ>vɾ:^5�333�%�T�����P�)��)���-���   ��l���;d��;d�'�O߼���D��=L��<�1��C��<j�#�
;ě�<��
<�9X<�/<�C�<T��;ě�:�o�#�
��`B1111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071519192017090715192020170907151919  IF  CORTCOOA5.2 RTQCGL01                                                        20130411062113  QCP$PSAL            G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20130408204542  QCF$TEMP            G�O�G�O�G�O�4               IF  CORTCOOA5.2 RTQCGL01                                                        20130411061255  QCP$TEMP            G�O�G�O�G�O�                IF  ARGQCOAR1.0                                                                 20130408162805  QCP$                G�O�G�O�G�O�06B68           IF  ARGQCOAR1.0                                                                 20130408162805  QCF$                G�O�G�O�G�O�04000           IF      SCOO1.4                                                                 20130410150315  QC                  G�O�G�O�G�O�                IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            C�ٚC���@�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            C��fC���@�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            C���C��3@�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            Co��C[� @�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            B�33Bn  @�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            BffA���@�                  IF  ARGQSCOO1.4                                                                 20130408162757  CF  TEMP            A�  Aff@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            C�ٚC���@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            C�ffC��3?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            C��fC���@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            C�� C���?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            C���C��3@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            Cy33Cy33?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            Co��C[� @�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            CQL�B���?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            B�33Bn  @�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            BY��B��?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            BffA���@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            A���A���?�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            A�  Aff@�                  IF  ARGQSCOO1.4                                                                 20130408162759  CF  PSAL            @�ff@�ff?�                  IF  CORTCOOA5.2 RTQCGL01                                                        20130408205747  QCF$PSAL            G�O�G�O�G�O�4               IF  ARSQOW  1.0 CTD2016V1                                                       20170907151920  IP  PSAL            @�ffC�ٚG�O�                IF      COFC3.2                                                                 20190603175239                      G�O�G�O�G�O�                