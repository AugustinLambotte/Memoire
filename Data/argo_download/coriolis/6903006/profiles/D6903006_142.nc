CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:25Z creation; 2023-01-04T00:23:16Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_054b      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    :�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    :�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  ;   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  ;D   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  ;�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        ;�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     <   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     <8   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    <X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         <\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    <d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            <h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           <p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           <x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    <�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    <�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    <�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    <�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        =�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�       =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  D�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  E�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  I   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  I�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Mx   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N\   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Q�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Ud   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  VH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Y�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Z�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^4   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  a�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  f    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  g   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    w    HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    w   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    w   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    w   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  w   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    wP   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    w`   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    wd   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         wt   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         wx   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        w|   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    w�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  j�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    j�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    n�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    r�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  v�Argo profile    3.1 1.2 19500101000000  20230104002025  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�X����1   @�X����@RB<�����+���`�8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �b����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �tb�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��M���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��(�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����݀  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���8�x  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���t�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���>�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��"�8  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����X  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��+<M`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��M^o�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��q�,  999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990   A#33A.ffA@  AQ��A\��Al��A~ffA�ffA�33A�ffA�  A�33A�ffA�33A�  A�33A�ffA�  A�  A�33A�33A�33A�33B��B��B  BffB33B��BffB��B"��B'��B,ffB0  B333B8  B=33B@ffBD  BG33BS33Bg33B{��B�33B�  B�  B�ffB�33B�  B�33B�33B�ffB�ffB�ffB�  C 33C�fC	��C��C��C��C�3C"�fC(  C-  C233C733C;��C@��CE�fCK  CO��CT�3CY�fC_  Cc�fCi�Cn�Cs�Cx33C}�C�  C�s3C��3C�� C��C�� C�  C�ffC�ٚC�L�C���C�@ C��fC���C��fC�s3C�  C�ffC��fC�L�C��fC�� C��C���C�&fC��fC�  C�Y�C��3CɌ�C��3C�ffC��fCӀ C��C�s3C�ٚC݀ C��C� C�ٚC� C��3C�Y�C��3C�s3C��C���C�  C�� C��D L�D�fD��D��D9�D�fD��D��D
&fDl�D� D��D9�D�fD��D�3D9�Ds3D��D�fD33Dy�D� D�3D&fDy�D �3D"3D#S3D$��D%� D&��D(9�D)�fD*��D+�fD-33D.��D/�fD0�3D233D3� D4� D6fD79�D8y�D9� D:��D<33D=l�D>�3D?�3DA@ DBy�DC� DD��DF@ DGy�DH�fDI�3DKFfDL� DM�3DN��DP&fDQl�DR��DS��DU33DV� DW�3DY  DZ33D[l�D\�fD^  D_33D`l�Da�fDb��Dd,�De� Df�3Dg�3DiL�Dj�fDk��Dm�Dn9�Dol�Dp�3Dq�f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A#33A.ffA@  AQ��A\��Al��A~ffA�ffA�33A�ffA�  A�33A�ffA�33A�  A�33A�ffA�  A�  A�33A�33A�33A�33B��B��B  BffB33B��BffB��B"��B'��B,ffB0  B333B8  B=33B@ffBD  BG33BS33Bg33B{��B�33B�  B�  B�ffB�33B�  B�33B�33B�ffB�ffB�ffB�  C 33C�fC	��C��C��C��C�3C"�fC(  C-  C233C733C;��C@��CE�fCK  CO��CT�3CY�fC_  Cc�fCi�Cn�Cs�Cx33C}�C�  C�s3C��3C�� C��C�� C�  C�ffC�ٚC�L�C���C�@ C��fC���C��fC�s3C�  C�ffC��fC�L�C��fC�� C��C���C�&fC��fC�  C�Y�C��3CɌ�C��3C�ffC��fCӀ C��C�s3C�ٚC݀ C��C� C�ٚC� C��3C�Y�C��3C�s3C��C���C�  C�� C��D L�D�fD��D��D9�D�fD��D��D
&fDl�D� D��D9�D�fD��D�3D9�Ds3D��D�fD33Dy�D� D�3D&fDy�D �3D"3D#S3D$��D%� D&��D(9�D)�fD*��D+�fD-33D.��D/�fD0�3D233D3� D4� D6fD79�D8y�D9� D:��D<33D=l�D>�3D?�3DA@ DBy�DC� DD��DF@ DGy�DH�fDI�3DKFfDL� DM�3DN��DP&fDQl�DR��DS��DU33DV� DW�3DY  DZ33D[l�D\�fD^  D_33D`l�Da�fDb��Dd,�De� Df�3Dg�3DiL�Dj�fDk��Dm�Dn9�Dol�Dp�3Dq�f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����T�Ձ�Ձ�Ձ�ա˿ա˿ա˿Ձ�Ձ�Ձ�Ձ�Ձ�ա˿ա˿ա˿Ձ�ԛ���J�У׿�p��͑h��������|�У׿�Z�׍P�ش9�����C��ڟ���b���#��r���Kǿ��ÿ�n���b��ȴ�Õ���&鿿|����r����H��Ĝ��dZ�̋D�� ſ�����/��^5�ۅ���H��"ѿ�+��9X������/��ff�������H��p��Y�#�-��푾D�������ě�=�t�>�u>��9>ɺ^>���?I�?'�?9��?N{?NV?S�F?V�+?W�P?XQ�?W�P?V�+?W��?W�P?G�?,I�?E�?	7L??)x�?8b?BJ?D��?I�^?T9X?T9X?MV?H��??;d?6ȴ?(��?!G�??&�?��>��m>�>�x�>�;d>��j>�5?>�j>��H>�ȴ>�&�>���>��D>��y>�S�>��>��+>���>���>��>�O�>�O�>�C�>�1'>�+>���>��>}�>w��>l�D>hr�>cS�>]/>Xb>Q�>M��>L��>I�^>H�9>J��>J��>@�>>v�><j>7K�>)��>%�T>$�/>#�
>��>z�>	7L>   =�h=�G�=�
==��=�{=��w=���=�\)=q��=D��=t�<���<�1<u<#�
;D�����
�T����j�o�,1�T���u��t����{��-��^5�ȴ9����/��`B��F�   ���+�O߾���-��R���!���.{�0 ž1&�333�333�5?}�5?}�6E��7KǾ8Q�7KǾ7KǾ7KǾ7KǾ9X�<j�;dZ�<j�?|�@��A�7�D���E��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ���T�Ձ�Ձ�Ձ�ա˿ա˿ա˿Ձ�Ձ�Ձ�Ձ�Ձ�ա˿ա˿ա˿Ձ�ԛ���J�У׿�p��͑h��������|�У׿�Z�׍P�ش9�����C��ڟ���b���#��r���Kǿ��ÿ�n���b��ȴ�Õ���&鿿|����r����H��Ĝ��dZ�̋D�� ſ�����/��^5�ۅ���H��"ѿ�+��9X������/��ff�������H��p��Y�#�-��푾D�������ě�=�t�>�u>��9>ɺ^>���?I�?'�?9��?N{?NV?S�F?V�+?W�P?XQ�?W�P?V�+?W��?W�P?G�?,I�?E�?	7L??)x�?8b?BJ?D��?I�^?T9X?T9X?MV?H��??;d?6ȴ?(��?!G�??&�?��>��m>�>�x�>�;d>��j>�5?>�j>��H>�ȴ>�&�>���>��D>��y>�S�>��>��+>���>���>��>�O�>�O�>�C�>�1'>�+>���>��>}�>w��>l�D>hr�>cS�>]/>Xb>Q�>M��>L��>I�^>H�9>J��>J��>@�>>v�><j>7K�>)��>%�T>$�/>#�
>��>z�>	7L>   =�h=�G�=�
==��=�{=��w=���=�\)=q��=D��=t�<���<�1<u<#�
;D�����
�T����j�o�,1�T���u��t����{��-��^5�ȴ9����/��`B��F�   ���+�O߾���-��R���!���.{�0 ž1&�333�333�5?}�5?}�6E��7KǾ8Q�7KǾ7KǾ7KǾ7KǾ9X�<j�;dZ�<j�?|�@��A�7�D���E��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB33B6FB7LB7LB7LB8RB8RB:^B;dB>wB?}BA�BK�BT�BffBx�B�+B�uB��B��B�B�LB��BĜB��B��B�HB�ZB�sB�B�B�B��B��B��B	1B	JB	�B	uB	\B	�B	!�B	)�B	9XB	B�B	N�B	l�B	�%B	��B	�wB	ǮB	��B	�
B	�ZB	�B	��B
+B
!�B
-B
D�B
cTB
�B
��B
�XB
�B
�BBPB�B�B/B<jBK�BN�B_;Bk�Br�Bz�B}�B~�B�B�B�B�B�B�%B�B� By�Bt�Bq�Bw�B�B�+B�JB�VB�oB��B��B��B��B��B��B�{B�uB�oB�oB�uB�{B�uB��B��B�oB�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B33B6FB7LB7LB7LB8RB8RB:^B;dB>wB?}BA�BK�BT�BffBx�B�+B�uB��B��B�B�LB��BĜB��B��B�HB�ZB�sB�B�B�B��B��B��B	1B	JB	�B	uB	\B	�B	!�B	)�B	9XB	B�B	N�B	l�B	�%B	��B	�wB	ǮB	��B	�
B	�ZB	�B	��B
+B
!�B
-B
D�B
cTB
�B
��B
�XB
�B
�BBPB�B�B/B<jBK�BN�B_;Bk�Br�Bz�B}�B~�B�B�B�B�B�B�%B�B� By�Bt�Bq�Bw�B�B�+B�JB�VB�oB��B��B��B��B��B��B�{B�uB�oB�oB�uB�{B�uB��B��B�oB�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002025                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002316  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002316  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A#33Dq�fG�O�                