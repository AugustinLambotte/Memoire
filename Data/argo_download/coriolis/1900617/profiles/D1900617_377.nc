CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   Q   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2016-09-04T08:37:04Z creation; 2019-06-03T17:53:39Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        D  :    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ;D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        D  ;�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  <�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     D  =0   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  >t   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ?�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  @   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  AP   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  A�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  B�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  D,   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  D�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  E�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  F   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  G\   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    G�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    J�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    M�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  P�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    P�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    P�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    P�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    P�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  P�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    Q   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    Q   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    Q   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Q,   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Q0   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Q4   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Q8Argo profile    3.1 1.2 19500101000000  20160904083704  20190603175339  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL              yA   IF  45123551                        2C  D   APEX                            2573                            n/a                             846 @��F�-��1   @��F�-��@P�/��w�)�I�^61   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   B   A   @���A��Ai��A�  A�  A���B��B��B/��BE��BY��Bm��B�  B���B�33B���B�  B�ffB�  B���B�ffB�33C��C33CL�C� C)33C3ffC=L�CGffCQ� C[ffCe  Co��Cy33C���C��fC�� C��3C��3C�s3C��fC���C�� C�� C�� C��3C��fC�� Cǳ3Č�Cљ�C֌�Cی�C���C��CꙚC�� C��3C��3C�� D�3D@ D� D	FfD�3DL�D��DS3D�3DFfDٚDL�D� D"9�D$��D'@ D)��D,l�D.�fD/�3222222222222222222222222222222222222222222222222222222222222222222222222222222222   �L��@���AffAd��A�ffA�33A�  B  B��B0��BD��BX��Bk33B�ffB���B�33B���B�  B���B�ffB�  B���B���C  C�CL�C$  C.33C8�CB33CLL�CV33C_��CjffCt  C}�fC��C��fC��C��C�ٚC��C��3C��fC��fC��fC��C��C�&fC��C��3C�  C��3C��3C�  C��3C�  C�&fC��C��C��fD �fD�3Ds3D��D
ffD  D� DfDffD��D��D  Ds3D ��D#l�D%�3D(� D+  D-��D.ff222222222222222222222222222222222222222222222222222222222222222222222222222222222   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�ZB!�B
ȴB
ǮB
��B
�bB�B  BoB&�B�BCbNAЃA�B�B�B�B�B�B�B�B"�B&�B,B49B<jB?}B?}BC�BH�BM�BP�B^5Bk�Bp�Bz�B~�B�B�B�+B�1B�1B�7B�7B�=    B�DB�JB�JB�JB�PB�PB�VB�VB�\B�\B�bB�bB�\B�VB�VB�VB�VB�VB�\B�\B�\B�\B�bB�bB�\B�bB�bB�bB�\B�bB�\B�bB�bB�bB�b441111111114441111111111111111111111111111111411111111111111111111111111111111111   G�O�G�O�B
ȴB
ǮB
��B
�bB�B  BoB&�B�G�O�G�O�G�O�B�B�B�B�B�B�B�B"�B&�B,B49B<jB?}B?}BC�BH�BM�BP�B^5Bk�Bp�Bz�B~�B�B�B�+B�1B�1B�7B�7B�=G�O�B�DB�JB�JB�JB�PB�PB�VB�VB�\B�\B�bB�bB�\B�VB�VB�VB�VB�VB�\B�\B�\B�\B�bB�bB�\B�bB�bB�bB�\B�bB�\B�bB�bB�bB�b441111111114441111111111111111111111111111111411111111111111111111111111111111111   G�O�G�O�<��
<��
<��
<��
<��
<��
<��
<��
<��
G�O�G�O�G�O�<��
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
G�O�<��
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
A\)A`BAdZA�\A�PA
A�@���@G�?�j?�n�?��?r-?Z�H?8��?#��?�R?�?�>�^5>ۥ�>�x�>�V>{�m>��>���>�n�>��9>]/>8Q�>"��=�;d=�\)>#�
>J��>)��>1&�>��>1'=���=�G�=�j=��
=�C�=aG�=t�<���<t�;�o�t���j�+�8Q�aG���O߽�����{��Q���`�   �t��'7KǾB�\�Kƨ�V�dZ�r�!�x����J���˾��^��\)��zᾙ���5?��Z��l����h������������111111111111111111111111111111111111111111111111111111111111111111111111111111111   A\)A`BAdZA�\A�PA
A�@���@G�?�j?�n�?��?r-?Z�H?8��?#��?�R?�?�>�^5>ۥ�>�x�>�V>{�m>��>���>�n�>��9>]/>8Q�>"��=�;d=�\)>#�
>J��>)��>1&�>��>1'=���=�G�=�j=��
=�C�=aG�=t�<���<t�;�o�t���j�+�8Q�aG���O߽�����{��Q���`�   �t��'7KǾB�\�Kƨ�V�dZ�r�!�x����J���˾��^��\)��zᾙ���5?��Z��l����h������������111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071551342017090715513420170907155134  IF  ARGQCOAR1.0                                                                 20160904081542  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20160904081542  QCF$                G�O�G�O�G�O�00CA40          IF  ARGQCOAR1.0                                                                 20160904081542  QCC$                G�O�G�O�G�O�00CA40          IF  CORTCOOA6.2 RTQCGL01                                                        20160905055901  QCF$TEMP            G�O�G�O�G�O�6               IF  CORTCOOA6.2 RTQCGL01                                                        20160905062339  QCF$PSAL            G�O�G�O�G�O�5               IF      SCOO0.19                                                                20160906173741  CF  PSAL            @���@���?�                  IF  CORTCOOA6.2 RTQCGL01                                                        20160907072833  QCP$PSAL            G�O�G�O�G�O�                IF  CORTCOOA6.2 RTQCGL01                                                        20160929072226  QCP$TEMP            G�O�G�O�G�O�                IF      SCOO0.19                                                                20161116153952  QC                  G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20170623001844  QCF$TEMP            G�O�G�O�G�O�6               IF  CODMCOOA6.2 DMQCGL01                                                        20170623004920  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO0.30                                                                20170907113059  CF  PSAL            BY��BY��@�                  IF      SCOO0.30                                                                20170907113059  CF  TEMP            C�� C�� @�                  IF      SCOO0.30                                                                20170907113059  CF  PSAL            B�33B�33@�                  IF      SCOO0.30                                                                20170907113059  CF  PSAL            C��3C��3@�                  IF      SCOO0.30                                                                20170907113059  CF  TEMP            Bm��B�  @�                  IF      SCOO0.30                                                                20170907113059  CF  PSAL            C�� C�� @�                  IF  ARSQOW  1.0 CTD2016V1                                                       20170907155134  IP  PSAL            @���D/�3G�O�                IF      COFC3.2                                                                 20190603175339                      G�O�G�O�G�O�                