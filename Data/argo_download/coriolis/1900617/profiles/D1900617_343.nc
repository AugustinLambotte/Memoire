CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   P   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2015-09-30T08:37:06Z creation; 2019-06-03T17:53:22Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        @  :    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  ;@   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        @  ;�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  <�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     @  =    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  >`   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  ?�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  ?�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  A0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  A�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  B�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  D    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  DP   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  E�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  E�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  G    SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    GP   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    JP   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    MP   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  PP   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    P|   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    P�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    P�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    P�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  P�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    P�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    P�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    P�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         P�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         P�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        P�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    P�Argo profile    3.1 1.2 19500101000000  20150930083706  20190603175322  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL              WA   IF  41857620                        2C  D   APEX                            2573                            n/a                             846 @�sG��1   @�sG��@Qc�����+St�j~�1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   A   A   @���AffA\��A�ffA���A陚B
  B��B2��BE��BZ  Bn  B�  B���B�33B�  B�  B�33B���Bƙ�B�  B���C� C  C  C� C)ffC3ffC=ffCG  CQ33C[  Ce  Co33Cy� C��3C�ffC���C�s3C�� C��fC��3C��fC��3C��3C�� C��fC��fC���Cǀ C̦fCљ�Cֳ3Cۀ C�3C�fCꙚCC��3C���C�s3D� D9�D�fD	` D��D@ D�3D9�D�3D9�D�3DFfD� D"9�D$� D'FfD)�3D,Y�D/  22222222222222222222222222222222222222222222222222222222222222222222222222222222�L��@�ffA	��Aa��A�33A�  A�ffB
  B  B0��BE33BY33Bm34B�33B���B���B���B���B�ffB�33BЙ�B�ffB���C��C��CL�C$33C.33C833CA��CL  CU��C_��Cj  CtL�C~34C���C��3C�ٚC�&fC��C��C��C��C��C�&fC��C��C�33C��fC��C�  C��C��fC��C��C�  C�  C��C��3C�ٚD �3D��Dy�D3D
� D�3DffD��D�fD��D�fD��Ds3D ��D#s3D%��D(�fD+�D-�322222222222222222222222222222222222222222222222222222222222222222222222222222222@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B
��B
�{B
�{B
�{B
�{B
�uB
�VB
jB�BH�BQ�BQ�BP�BM�BJ�BE�BA�B<jB:^B;dB7LB?}BM�BT�BO�BS�BN�B\)BffBe`B\)BbNBjBjBq�B{�B�B�1B�7B�=B�=B�7B�=B�7B�=B�=B�DB�DB�JB�JB�PB�VB�\B�VB�\B�\B�bB�hB�hB�oB�hB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�hB�hB�hB�oB�oB�oB�o11111111111111111111111111111111111111111111111111111111111111111111111111111111B
��B
�{B
�{B
�{B
�{B
�uB
�VB
jB�BH�BQ�BQ�BP�BM�BJ�BE�BA�B<jB:^B;dB7LB?}BM�BT�BO�BS�BN�B\)BffBe`B\)BbNBjBjBq�B{�B�B�1B�7B�=B�=B�7B�=B�7B�=B�=B�DB�DB�JB�JB�PB�VB�\B�VB�\B�\B�bB�hB�hB�oB�hB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�hB�hB�hB�oB�oB�oB�o11111111111111111111111111111111111111111111111111111111111111111111111111111111<��
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
@�t�@�t�@�t�@�|�@�t�@��@��@���?��D?��+?�V?���?s�F?ix�?V?:�H?%�?�
? �?M�>��w>���>��j?�>�\)>�"�>   >{�m>�"�>s�F=@�=+=�E�:�o=0 �=�->'�>��>��=��=��=u=��<#�
��9X����w�8Q�D����hs��E����
�� Žȴ9������S���S������l���l��hs����$�/�49X�:^5�P�`�e`B�u�}󶾆��C���bN���+�������㾠Ĝ���y��l�������11111111111111111111111111111111111111111111111111111111111111111111111111111111@�t�@�t�@�t�@�|�@�t�@��@��@���?��D?��+?�V?���?s�F?ix�?V?:�H?%�?�
? �?M�>��w>���>��j?�>�\)>�"�>   >{�m>�"�>s�F=@�=+=�E�:�o=0 �=�->'�>��>��=��=��=u=��<#�
��9X����w�8Q�D����hs��E����
�� Žȴ9������S���S������l���l��hs����$�/�49X�:^5�P�`�e`B�u�}󶾆��C���bN���+�������㾠Ĝ���y��l�������11111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071539462017090715394620170907153945  IF  ARGQCOAR1.0                                                                 20150930081537  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20150930081537  QCF$                G�O�G�O�G�O�008000          IF  ARGQCOAR1.0                                                                 20150930081537  QCC$                G�O�G�O�G�O�008000          IF  CORTCOOA6.2 RTQCGL01                                                        20151001150743  QCP$TEMP            G�O�G�O�G�O�                IF  CORTCOOA6.2 RTQCGL01                                                        20151001154536  QCP$PSAL            G�O�G�O�G�O�                IF  ARSQOW  1.0 CTD2016V1                                                       20170907153946  IP  PSAL            @���D/  G�O�                IF      COFC3.2                                                                 20190603175322                      G�O�G�O�G�O�                