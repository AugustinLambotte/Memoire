CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   P   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2016-08-25T10:37:04Z creation; 2019-06-03T17:53:38Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
_FillValue                    P�Argo profile    3.1 1.2 19500101000000  20160825103704  20190603175338  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL              xA   IF  45025845                        2C  D   APEX                            2573                            n/a                             846 @���*���1   @���*���@P���`A��)�����1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   B   A   @���A��Ai��A���Ař�A홚BffB  B/33BE��BX��Bk��B�33B�  B���B���B���B�  B�ffB�  B�  B���C33C� C�C  C)ffC3  C=ffCGffCQffCZ�fCe  Co33Cx�fC��3C��fC���C�� C�� C���C���C���C��fC�� C�� C���C���C�CǙ�C̳3Cљ�C֌�CۦfC���C� C��C�fC� C�� C��3D� D@ D�3D	@ D��D9�D�fDY�D� DY�D��DFfD� D"@ D$��D'l�D)�3D,y�D-  22222222222222222222222222222222222222222222222222222222222222222222222222222222����@y��AffAfffA�  A�  A�34B33BffB0��BD  BV��Bk��B���B�ffB�33B�ffB���B�  B���BЙ�B�ffB�  CL�C�fC��C$33C-��C833CB33CL33CU�3C_��Cj  Cs�3C~34C��C�  C��fC�&fC�  C�33C�  C��C�&fC�&fC�  C��3C�  C�  C��C�  C��3C��C��3C��fC��3C��C��fC�&fC��D �3D�3D�fD�3D
l�D��Dy�D�Ds3D�D� D��Ds3D �3D#� D&  D(�fD+,�D+�322222222222222222222222222222222222222222222222222222222222222222222222222222222@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B
�%B
�%B
�B
�+B
�B
n�B
ÖB
�BBBB+B
=A�1'A��/BhB{B�B�B�B �B$�B.B2-B8RB:^B<jBA�BG�BK�BS�BZBn�Bk�Bp�Bu�B}�B�    B�%B�+B�1B�=B�=B�DB�DB�DB�JB�JB�PB�VB�PB�VB�bB�bB�bB�\B�bB�hB�hB�hB�oB�hB�hB�oB�bB�bB�\B�\B�\B�\B�bB�bB�bB�\B�\B�bB�bB�bB�b11111111111114411111111111111111111111411111111111111111111111111111111111111111B
�%B
�%B
�B
�+B
�B
n�B
ÖB
�BBBB+B
=G�O�G�O�BhB{B�B�B�B �B$�B.B2-B8RB:^B<jBA�BG�BK�BS�BZBn�Bk�Bp�Bu�B}�B�G�O�B�%B�+B�1B�=B�=B�DB�DB�DB�JB�JB�PB�VB�PB�VB�bB�bB�bB�\B�bB�hB�hB�hB�oB�hB�hB�oB�bB�bB�\B�\B�\B�\B�bB�bB�bB�\B�\B�bB�bB�bB�b11111111111114411111111111111111111111411111111111111111111111111111111111111111<��
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
G�O�G�O�<��
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
<��
<��
<��
<��
<��
<��
A��A�A�AȴA��A�@�G�@"�H?��/?�$�?X��?:^5?.�?#o?��>��m>�>�/>��>�G�>�7L>���>��j>��>���>{�m>O�;>6E�>=p�>&�y>!��>8Q�>��>bN>�-=�/>�u>%�T>�=�x�=��=�%=<j='�<�/<��ͻD����o�C���h�'e`B���T��hs���w������o���
=q��u���.{�=p��J���^5?�p�׾�����1'��O߾�hs�������+���R��S���xվ��F���F���F���F11111111111111111111111111111111111111111111111111111111111111111111111111111111A��A�A�AȴA��A�@�G�@"�H?��/?�$�?X��?:^5?.�?#o?��>��m>�>�/>��>�G�>�7L>���>��j>��>���>{�m>O�;>6E�>=p�>&�y>!��>8Q�>��>bN>�-=�/>�u>%�T>�=�x�=��=�%=<j='�<�/<��ͻD����o�C���h�'e`B���T��hs���w������o���
=q��u���.{�=p��J���^5?�p�׾�����1'��O߾�hs�������+���R��S���xվ��F���F���F���F11111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071551202017090715512020170907155120  IF  ARGQCOAR1.0                                                                 20160825101538  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20160825101538  QCF$                G�O�G�O�G�O�00CA40          IF  ARGQCOAR1.0                                                                 20160825101538  QCC$                G�O�G�O�G�O�00CA40          IF  CORTCOOA6.2 RTQCGL01                                                        20160826070044  QCF$TEMP            G�O�G�O�G�O�5               IF  CORTCOOA6.2 RTQCGL01                                                        20160826072338  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO0.19                                                                20160826120424  QC                  G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20161226172837  QCP$TEMP            G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20170623001844  QCF$TEMP            G�O�G�O�G�O�5               IF  CODMCOOA6.2 DMQCGL01                                                        20170623004920  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO0.30                                                                20170907113108  CF  TEMP            B�33B���@�                  IF      SCOO0.30                                                                20170907113108  CF  PSAL            C�� C�� @�                  IF      SCOO0.30                                                                20170907113108  CF  PSAL            C���C���@�                  IF      SCOO0.30                                                                20170907113108  CF  PSAL            B���B���@�                  IF      SCOO0.30                                                                20170907113108  CF  PSAL            B�33B�33@�                  IF      SCOO0.30                                                                20170907113108  CF  TEMP            C���C�� @�                  IF  ARSQOW  1.0 CTD2016V1                                                       20170907155120  IP  PSAL            @���D-  G�O�                IF      COFC3.2                                                                 20190603175338                      G�O�G�O�G�O�                