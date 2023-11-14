CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   Q   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       2012-12-06T11:49:16Z creation      
references        (http://www.argodatamgt.org/Documentation   comment           user_manual_version       2.4    Conventions       Argo-2.4 CF-1.6    featureType       trajectoryProfile         >   	DATA_TYPE                  	long_name         	Data type      
_FillValue                    2�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    2�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    2�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    2�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    3   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    3   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    3,   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  34   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  3t   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  3�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        3�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    3�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    3�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     3�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    4   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    4   INST_REFERENCE                    	long_name         Instrument type    conventions       Brand, type, serial number     
_FillValue                  @  4   FIRMWARE_VERSION                  	long_name         Instrument version     conventions           
_FillValue                    4X   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    4h   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T           4l   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    4t   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            4x   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           4�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           4�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    4�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    4�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    4�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    4�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    4�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    4�   PRES         
      	   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     D  5�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  6�   PRES_ADJUSTED            
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     D  7@   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  8�   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     D  8�   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  :   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ;`   TEMP_ADJUSTED            
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  ;�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  <�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  =L   PSAL         
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_practical_salinity   
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  >�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ?�   PSAL_ADJUSTED            
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  @(   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  Al   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  A�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  C   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    C4   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    F4   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    I4   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    
_FillValue                  ,  L4   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    L`   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    Ld   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    Lh   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    Ll   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  Lp   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    L�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    L�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    L�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         L�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         L�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        L�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    L�Argo profile    2.3 1.2 19500101000000  20121128054423  20130702144515  6901902 ARGO POLAND                                                     Waldemar WALCZOWSKI                                             PRES            PSAL            TEMP               1A   IF  28016469                        2C  D   NEMO Profiling Float                                                            860 @�p�N �1   @�p�N �@SF�����@$���D�1   GPS     B   B   B                                                                                                                                                                                                                                                                   =���>L��>L��>���?   ?   >���@���@ٙ�A`  A�  A噚B33BD  Bk33B�  B�33B�  B�ffB���B���C � C
��C�C�C(�3C2��C<�fCG�CP�fCZ� Cd� Cn��Cx�fC���C���C��fC��fC�Y�C���C���C�&fC��3C�Y�C���C�Y�CǙ�C�  C�Y�C�3C�L�D� D	L�D�3DFfD�fD"FfD(` D.�3D;33DG��DT&fD`� Dm33Dyl�D�#3D�c3D��fD��fD�#3D�\�D���D���D�,�D��D�l�D��3D�3D�ffDڣ3D�,�114114411111111111111111111111111111111111111111111111111111111111111111111111111   =���>L��G�O�>���?   G�O�G�O�@���@ٙ�A`  A�  A噚B33BD  Bk33B�  B�33B�  B�ffB���B���C � C
��C�C�C(�3C2��C<�fCG�CP�fCZ� Cd� Cn��Cx�fC���C���C��fC��fC�Y�C���C���C�&fC��3C�Y�C���C�Y�CǙ�C�  C�Y�C�3C�L�D� D	L�D�3DFfD�fD"FfD(` D.�3D;33DG��DT&fD`� Dm33Dyl�D�#3D�c3D��fD��fD�#3D�\�D���D���D�,�D��D�l�D��3D�3D�ffDڣ3D�,�114114411111111111111111111111111111111111111111111111111111111111111111111111111   @��@��G�O�@��@��G�O�G�O�@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@yG�@jM�@w|�@~��@|�@
=@{�@��@��@���@���@��`@���@���@��j@��/@���@�z�@�Z@�9X@�9X@��@|��@{@zM�@y�#@z�@z�\@y&�@w|�@w\)@w�@u�-@r�H@m��@iX@f�+@]V@\9X@Y�^@X��@Q�@IG�@D�j@<��@2�\@"�H@��?�V?ۅ?�b?�&�?�$�?D�/?,��>��>�X>t�j>+    ��l��~�۾�-��V��~���7�	�^���%`B�,I��49X�:���?;d�C�
�I�^�N���Rn��W
=�W
=�V�+144111441111111111111111111111111111111111111111111111111111111111111111111111111   @yG�G�O�G�O�@~��@|�G�O�G�O�G�O�@��@���@���@��`@���@���@��j@��/@���@�z�@�Z@�9X@�9X@��@|��@{@zM�@y�#@z�@z�\@y&�@w|�@w\)@w�@u�-@r�H@m��@iX@f�+@]V@\9X@Y�^@X��@Q�@IG�@D�j@<��@2�\@"�H@��?�V?ۅ?�b?�&�?�$�?D�/?,��>��>�X>t�j>+    ��l��~�۾�-��V��~���7�	�^���%`B�,I��49X�:���?;d�C�
�I�^�N���Rn��W
=�W
=�V�+144114441111111111111111111111111111111111111111111111111111111111111111111111111   ;oG�O�G�O�;o;oG�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B�Bp�B�7B��B��B�BB�B�B�B�B�B�B�B�+B�%B�%B�%B�%B�%B�%B�B�B� B~�B~�B�B�B�B� B� B~�B}�Bz�Bv�Bo�Bl�Be`BcTB`BB^5BVBM�BG�B@�B6FB&�B�B%B��B�`B�)B��BB�wB�-B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��144111441111111111111111111111111111111111111111111111111111111111111111111111111   B��G�O�G�O�B�7B��G�O�G�O�G�O�B�B�B�B�B�B�B�+B�%B�%B�%B�%B�%B�%B�B�B� B~�B~�B�B�B�B� B� B~�B}�Bz�Bv�Bo�Bl�Be`BcTB`BB^5BVBM�BG�B@�B6FB&�B�B%B��B�`B�)B��BB�wB�-B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��144114441111111111111111111111111111111111111111111111111111111111111111111111111   <#�
G�O�G�O�<#�
<#�
G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            PSAL            TEMP            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          201307021445152013070214451520130702144515  IF      SCOO1.4                                                                 20121206115254  QC                  G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20121128174921  QCF$PSAL            G�O�G�O�G�O�4               IF  CORTCOOA5.2 RTQCGL01                                                        20121128174219  QCF$TEMP            G�O�G�O�G�O�4               IF  ARGQCOAR1.0                                                                 20121128053502  QCP$                G�O�G�O�G�O�DEBFC           IF  ARGQCOAR1.0                                                                 20121128053502  QCF$                G�O�G�O�G�O�04100           GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2013V01 + ARGO climatology 20130702144515  IP  PSAL            =���D�,�G�O�                