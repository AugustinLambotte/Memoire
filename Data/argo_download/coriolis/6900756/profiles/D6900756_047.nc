CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2012-11-02T05:43:08Z creation; 2018-01-08T15:14:32Z last update (coriolis COFC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    6�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    6�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7    PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7@   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    7�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    7�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     7�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8$   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8D   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           8H   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8P   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >��	4E�   
_FillValue        A.�~            8T   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8\   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8d   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8l   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8p   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    8x   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    8|   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z         �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  :P   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z         �  :�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  ;L   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���      �  ;�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  <H   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  =   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  =D   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  >   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  >@   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  ?   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  ?�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  @   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  4  @�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o      �  A    HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    K$   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    K(   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    K,   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    K0   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  K4   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    Kt   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    K�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    K�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         K�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         K�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        K�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    K�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  A�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    A�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    D�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    G�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  J�Argo profile    3.1 1.2 19500101000000  20121102054308  20180108151432  6900756 DAP                                                             Andreas STERL                                                   PRES            TEMP            PSAL               /A   IF  27794626                        2C  D   APEX                            5790                            021009                          846 @�i�Z��1   @�i�io @O&$�/��A�=p��1   ARGOS   A   A   A   Primary sampling: discrete []                                                                                                                                                                                                                                      @�  A!��A���A�B ffBI��Bq��B�ffB�ffB�  B�33Bܙ�B�33C33C�3C�3C   C)��C4ffC>�CG�3CR� C\  Cf� Co�fCzffC��fC�&fC�@ C�@ C��3C�  C�&fC�L�C�L�C�  C��C��C�@ C�� C�&fD	� DfD"� D/�D;�3DH  DT�fD`�3Dm� 11111111111111111111111111111111111111111111111111  @ٙ�A&fgA�  A�  B!��BJ��Br��B�  B�  B���B���B�34B���C� C  C  C L�C)�gC4�3C>fgCH  CR��C\L�Cf��Cp33Cz�3C��C�L�C�ffC�ffC��C�&fC�L�C�s3C�s3C�&fC�33C�@ C�ffC��fC�L�D	�3D�D"�3D/,�D;�fDH3DT��DafDm�311111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�  @�b@�(�@� �@�I�@�9X@�9X@�1'@�-@��m@�7L@�x�@ǥ�@�{@��@�\)@�7L@�K�@���@��h@�%@��F@�J@�r�@��@��^@�&�@�1'@��@���@���@���@�%@�(�@�~�@���@�ff@���@���@��D@�j@�Q�@�I�@�9X@��#@�ƨ@���@�{11111111111111111111111111111111111111111111111111  @��@��@�  @�b@�(�@� �@�I�@�9X@�9X@�1'@�-@��m@�7L@�x�@ǥ�@�{@��@�\)@�7L@�K�@���@��h@�%@��F@�J@�r�@��@��^@�&�@�1'@��@���@���@���@�%@�(�@�~�@���@�ff@���@���@��D@�j@�Q�@�I�@�9X@��#@�ƨ@���@�{11111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBVBVB\B\BbB\B\BVBVBPB@�BN�BT�BYBZBYBXBVBQ�BM�BI�BL�BK�BH�BE�B@�B<jB9XB8RB5?B0!B,B(�B$�B"�B �B�B�BuBhBhBbB\B\BVBPBB�B�NB��11111111111111111111111111111111111111111111111111  B.B.B4B4B:B4B4B.B.B(B@[BN�BT�BX�BY�BX�BW�BU�BQ�BM�BI�BL�BK�BH�BEzB@[B<BB90B8*B5B/�B+�B(�B$�B"�B �B�B^BLB?B@B:B4B4B.B(B�B�B�&B��11111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Surface pressure = -0.3 dbar                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted by using pressure offset at the sea surface. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                     No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected (salinity adjusted for pressure offset). OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                      201504011553082015040115530820150401155308  IF  CORTCOOA5.2 RTQCGL01                                                        20121102160915  QCF$TEMP            G�O�G�O�G�O�5               IF  CORTCOOA5.2 RTQCGL01                                                        20121106170441  QCP$PSAL            G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20121106170003  QCP$TEMP            G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20121102161611  QCF$PSAL            G�O�G�O�G�O�5               IF  ARGQSCOO1.4                                                                 20121106115445  CF  TEMP            D`�3C�� ?�                  IF  ARGQSCOO1.4                                                                 20121106115445  CF  TEMP            ���ͽ���?�                  IF  ARGQSCOO1.4                                                                 20121106115445  CF  PSAL            ���ͽ���@@                  IF  ARGQCOAR1.0                                                                 20121102052758  QCP$                G�O�G�O�G�O�DEBFC           IF  ARGQCOAR1.0                                                                 20121102052758  QCF$                G�O�G�O�G�O�00900           IF  ARGQCOAR1.0                                                                 20141125165920  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20141125165920  QCF$                G�O�G�O�G�O�000000          IF  ARGQCOAR1.0                                                                 20141125165920  QCC$                G�O�G�O�G�O�000000          IF  ARGQCOAR1.0                                                                 20131014183752  QCP$                G�O�G�O�G�O�9EBFC           IF  ARGQCOAR1.0                                                                 20131014183752  QCF$                G�O�G�O�G�O�00000           IF  ARGQCOAR1.0                                                                 20131014183752  QCC$                G�O�G�O�G�O�00000           IF      SCOO1.4                                                                 20131115133333  QC                  G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20140818121343  QCF$PSAL            D/  Dz  G�O�6               IF  CODMCOOA6.2 DMQCGL01                                                        20140818114558  QCF$TEMP            D/  Dz  G�O�6               GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2013V01 + ARGO climatology 20150401155308  IP  PSAL            @�  Dm� G�O�                IF      COFC3.0                                                                 20180108151432                      G�O�G�O�G�O�                