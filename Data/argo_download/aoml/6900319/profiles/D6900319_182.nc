CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   /   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2015-05-04T22:31:17Z creation      
references        (http://www.argodatamgt.org/Documentation   comment       	free text      user_manual_version       3.03   Conventions       Argo-3.0 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      
_FillValue                    4�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    4�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    4�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    4�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    4�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    5    PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    5   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  5   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  5X   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  5�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        5�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    5�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    5�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     5�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    5�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    5�   PLATFORM_TYPE                     	long_name         Type of float      
_FillValue                     5�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                    6   FIRMWARE_VERSION                  	long_name         Instrument version     
_FillValue                    6,   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    6<   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T           6@   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    6H   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            6L   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           6T   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           6\   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    6d   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    6h   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    6p   CONFIG_MISSION_NUMBER                  	long_name         'Float's mission number for each profile    conventions       @0..N, 0 : launch mission (if exists), 1 : first complete mission   
_FillValue         ��        7p   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    7t   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    7x   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    7|   PRES         
      
   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z         �  7�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  8<   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���      �  8l   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  9(   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���      �  9X   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  :   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  :�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  ;    TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  ;�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  ;�   PSAL         
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  <�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  =d   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  =�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  0  >P   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �  >�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ?<   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ?l   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Bl   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    El   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    
_FillValue                  ,  Hl   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    H�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    H�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    H�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    H�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  H�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    H�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    H�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    H�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         I   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         I   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        I   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    IArgo profile    3.0 1.2 19500101000000  20150504223117  20170420153733  6900319 NAVY, Argo equivalent                                           CARL SZCZECHOWSKI                                               PRES            TEMP            PSAL               �A   AO  4941_112838_182                 2C  D   APEX                            6146            062608          846 @�N\Y� 1   @�N��� @Q�z�G����`A�71   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   A   A   @���A  As33A�33A�33A홚BffB ffB2  BF��BZffBp��B�ffB���B�  B�  B�  B�  B�33BǙ�Bљ�B���B�C�C  C�3C4  CH� Ca  Cy�3C�� C�L�C���C��C�� C�  C�@ C�  D	s3D�D(��D;s3DNS3Da  Dy��D���D���11111111111111111111111111111111111111111111111 @���A  As33A�33A�33A홚BffB ffB2  BF��BZffBp��B�ffB���B�  B�  B�  B�  B�33BǙ�Bљ�B���B�C�C  C�3C4  CH� Ca  Cy�3C�� C�L�C���C��C�� C�  C�@ C�  D	s3D�D(��D;s3DNS3Da  Dy��D���D���11111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��+@���@��\@��!@��R@���@��\@�n�@��\@��\@�v�@��\@�^5@���@�^5@�n�@��+@���@���@�X@�&�@���@���@�M�@�O�@a�#@K��@@�9@F�@&�+@^5?�P?׮?�dZ?���?d�?!�7>���>0 �<����,1�j~���xվ��D����Mӿ�`11111111111111111111111111111111111111111111111 @��+@���@��\@��!@��R@���@��\@�n�@��\@��\@�v�@��\@�^5@���@�^5@�n�@��+@���@���@�X@�&�@���@���@�M�@�O�@a�#@K��@@�9@F�@&�+@^5?�P?׮?�dZ?���?d�?!�7>���>0 �<����,1�j~���xվ��D����Mӿ�`11111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBYBYBYBYBYBYBYBYBYBXBXBXBZBXBW
BW
BW
BW
BS�BS�BS�BP�BM�BL�BE�B�B1BB�B��B�HB�B�B��B�wB�RB�B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111 BYBYBYBYBYBYBYBYBYBXBXBXBZBXBW
BW
BW
BW
BS�BS�BS�BP�BM�BL�BE�B�B1BB�B��B�HB�B�B��B�wB�RB�B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Surface pressure = 0 dbar                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted by using pressure offset at the sea surface. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                     No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected (salinity adjusted for pressure offset). OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                      201704201537342017042015373520170420153734  AO  ARCAADJP                                                                    20150504223117    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20150504223117  QCP$                G�O�G�O�G�O�0               AO  ARGQQCPL                                                                    20150504223117  QCF$                G�O�G�O�G�O�0               GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2016V01 + ARGO climatology 20170420153735  IP  PSAL            @���D���G�O�                