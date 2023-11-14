CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   Q   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2016-07-06T07:37:07Z creation; 2019-06-03T17:53:36Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
_FillValue                    Q8Argo profile    3.1 1.2 19500101000000  20160706073707  20190603175336  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL              sA   IF  44505216                        2C  D   APEX                            2573                            n/a                             846 @׹F�֩&1   @׹F�֩&@Q	��l�D�*�+I�1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   B   B   @�33AffAi��A�  A�33A�ffBffB  B1��BF  BX  BnffB�33B�  B�33B�  B���B�  B���B�ffBڙ�BC�C�CL�C33C)�C3�C=33CGffCQ33C[  Ce� CoffCyffC�� C���C�s3C��3C��3C�� C�s3C�� C��3C��3C���C��3C���C�� Cǌ�C̀ Cь�Cֳ3C�� C���C�� C�� C� C� C���C�s3DٚDS3D�3D	Y�D� D33D�3D9�D��DY�D�fD9�D��D"9�D$�fD'L�D)�3D,FfD.��D/Ff222222222222222222222222222222222222222222222222222222222222222222222222222222222   �   @���A  AfffA�ffA���A�  B	��B33B1��BC��BZ  Bn  B���B�  B���B���B���B�ffB�33B�ffB�ffB�  C  C33C�C$  C.  C8�CBL�CL�CU�fC`ffCjL�CtL�C~ffC��C��fC�&fC�&fC�33C��fC��3C�&fC�&fC��C�&fC�@ C�33C�  C��3C�  C�&fC�33C�  C�33C�33C��3C��3C��C��fD �3D�D��D3D
y�D��Dl�D�3D�fD3D� D�3Ds3D �3D#� D&fD(l�D+  D-�fD.  222222222222222222222222222222222222222222222222222222222222222222222222222222222   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B
��B
ǮB
ƨB
�hB
�;B
�B
��B
��B
��B  BBDBJBPB	7B%BBBBBBBB+BuB �B&�B6FB9XBI�BT�BbNBo�Bw�B�B�B�+B�1B�7B�DB�DB�DB(�B�DB�DB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�VB�\B�\B�bB�\B�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�b        B�h114411111111111111111111111111111111111111411111111111111111111111111111111111441   B
��B
ǮG�O�G�O�B
�;B
�B
��B
��B
��B  BBDBJBPB	7B%BBBBBBBB+BuB �B&�B6FB9XBI�BT�BbNBo�Bw�B�B�B�+B�1B�7B�DB�DB�DG�O�B�DB�DB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�VB�\B�\B�bB�\B�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bG�O�G�O�B�h114411111111111111111111111111111111111111411111111111111111111111111111111111441   <��
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
G�O�G�O�<��
@�l�@�  @��/@���@ �?���?�Z?�$�?Y�#?6E�?)�^?(r�?&��?%�?�D>�M�>��7>�`B>�bN>��>["�>+=�=��=ix�>��=���>@�=���=@�=]/=��P=�t�=��>	7L=��=\=���=@�=o<t���1��P�]/��t�����ͽ��������پC����)��1&�;dZ�C���I�^�R�]/�hr��q���y�#��%�|푾����\)��녾�񪾗
=�������R��MӾ�ff�����h���!��KǾ��H��|�\��J114411111111111111111111111111111111111111111111111111111111111111111111111111111   @�l�@�  G�O�G�O�@ �?���?�Z?�$�?Y�#?6E�?)�^?(r�?&��?%�?�D>�M�>��7>�`B>�bN>��>["�>+=�=��=ix�>��=���>@�=���=@�=]/=��P=�t�=��>	7L=��=\=���=@�=o<t���1��P�]/��t�����ͽ��������پC����)��1&�;dZ�C���I�^�R�]/�hr��q���y�#��%�|푾����\)��녾�񪾗
=�������R��MӾ�ff�����h���!��KǾ��H��|�\��J114411111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071545132017090715451320170907154512  IF  ARGQCOAR1.0                                                                 20160706071602  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20160706071602  QCF$                G�O�G�O�G�O�00CA40          IF  ARGQCOAR1.0                                                                 20160706071602  QCC$                G�O�G�O�G�O�00C840          IF  CORTCOOA6.2 RTQCGL01                                                        20160707075147  QCF$TEMP            G�O�G�O�G�O�6               IF  CORTCOOA6.2 RTQCGL01                                                        20160707081510  QCP$PSAL            G�O�G�O�G�O�                IF  CORTCOOA6.2 RTQCGL01                                                        20160708071424  QCP$TEMP            G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20170622152430  QCF$TEMP            G�O�G�O�G�O�6               IF  CODMCOOA6.2 DMQCGL01                                                        20170622155240  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO0.30                                                                20170907113154  CF  PSAL            D)�3D)�3@�                  IF      SCOO0.30                                                                20170907113154  CF  PSAL            C��3C��3@�                  IF      SCOO0.30                                                                20170907113154  CF  PSAL            C�s3C�s3@�                  IF      SCOO0.30                                                                20170907113154  CF  TEMP            D)�3D,Ff@�                  IF      SCOO0.30                                                                20170907113154  CF  TEMP            C�� C��3@�                  IF  ARSQOW  1.0 CTD2016V1                                                       20170907154513  IP  PSAL            @�33D/FfG�O�                IF      COFC3.2                                                                 20190603175336                      G�O�G�O�G�O�                