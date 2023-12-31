CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:24Z creation; 2023-01-04T00:23:03Z last update (coriolis COQC software)   
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
resolution        5�7�     �  =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Dx   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        l  ET   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        l  I�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     l  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  QP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  T�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  U�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Y   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  Y�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  ]L   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  `�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  a�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  e    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     l  e�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    u�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    u�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    u�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    u�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  u�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    v   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    v    HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    v$   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         v4   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         v8   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        v<   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    v@   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  iH   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    i�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    m�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    q�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  u�Argo profile    3.1 1.2 19500101000000  20230104002024  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               vA   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�
Y�q�1   @�
Y�q�@R�0�K^����t�8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �Q�?V   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �vx���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��m�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��$��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��7�H�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���8�x  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����'p  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����X  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���֩$  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��3�JT  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��lwؐ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��	{B`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��`�`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����	  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��M^o�  999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990 A!��A0  A@  ANffA^ffAl��A~ffA�  A���A���A�  A�  A�33A�  A�ffA�  A���A�  A�33A陚A�  A�ffB ffB  B��BffB��B  B  BffB   B$  B(  B,  B0��B3��B8  B<��B?33BC��BI33BTffBh  B|ffB���B�ffB���B�33B�ffB�ffB�ffB�ffB���B�33B�ffB�  C �CL�C	�fC�3C  C  C�fC#  C(�C-�C1��C6�fC<  C@��CE��CJ��CP�CT��CY��C_  CdL�Ci�Cm�fCs  Cx�C|�fC�ٚC�� C�&fC���C�&fC�� C�ٚC�� C��C�� C��fC�ffC��3C���C��fC�L�C�� C�s3C��C�� C��fC�L�C�  C�s3C�� C�@ C���C�L�C�ٚC�Y�C���C�Y�C�  CӦfC��Cؙ�C��C�s3C���C�Y�C��3C� C��C왚C��3C�L�C���C�Y�C��3C�� C��D 33DffD��D�3D,�D� D�fD�fD
33D� D��D  DS3D��D�fDfDFfD�fD�fDfDL�D��D�3D��D  Dy�D � D!��D#33D$l�D%��D'3D(S3D)�3D*��D+�3D-FfD.� D/��D0�3D2,�D3� D4�fD5��D733D8y�D9�fD:��D<S3D=��D>��D@fDA@ DBy�DC�3DD�3DF,�DGffDH�fDI��DKL�DL�3DM��DN�fDP33DQ� DR��DT  DU9�DVs3DW��DX��DZ33D[y�D\� D^�D_@ D`s3Da�3Db�3Dd9�Dey�Df�3Dg��DiFfDl9�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A!��A0  A@  ANffA^ffAl��A~ffA�  A���A���A�  A�  A�33A�  A�ffA�  A���A�  A�33A陚A�  A�ffB ffB  B��BffB��B  B  BffB   B$  B(  B,  B0��B3��B8  B<��B?33BC��BI33BTffBh  B|ffB���B�ffB���B�33B�ffB�ffB�ffB�ffB���B�33B�ffB�  C �CL�C	�fC�3C  C  C�fC#  C(�C-�C1��C6�fC<  C@��CE��CJ��CP�CT��CY��C_  CdL�Ci�Cm�fCs  Cx�C|�fC�ٚC�� C�&fC���C�&fC�� C�ٚC�� C��C�� C��fC�ffC��3C���C��fC�L�C�� C�s3C��C�� C��fC�L�C�  C�s3C�� C�@ C���C�L�C�ٚC�Y�C���C�Y�C�  CӦfC��Cؙ�C��C�s3C���C�Y�C��3C� C��C왚C��3C�L�C���C�Y�C��3C�� C��D 33DffD��D�3D,�D� D�fD�fD
33D� D��D  DS3D��D�fDfDFfD�fD�fDfDL�D��D�3D��D  Dy�D � D!��D#33D$l�D%��D'3D(S3D)�3D*��D+�3D-FfD.� D/��D0�3D2,�D3� D4�fD5��D733D8y�D9�fD:��D<S3D=��D>��D@fDA@ DBy�DC�3DD�3DF,�DGffDH�fDI��DKL�DL�3DM��DN�fDP33DQ� DR��DT  DU9�DVs3DW��DX��DZ33D[y�D\� D^�D_@ D`s3Da�3Db�3Dd9�Dey�Df�3Dg��DiFfDl9�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���� ſ�bN��녿������`�� ſ�vɿޗ������A���vɿ�5?�����/���m��ƨ�ۅ���H��"ѿ�=q�ؓu�׍P�׮��
=��l����y��`B�Ѓ���y���u�������j��Z��Mӿ�\)���!�[�m�CS��A���A�7�@  �3t���C�>8Q�?~�R?~�R?��?��9?�{?��?��7?���?�1?�r�?ļj?ǍP?�I�?��F?���?�1'?��?���?��u?�ȴ?��T?��?��?���?� �?�v�?�p�?��?�Q�?��u?�K�?�$�?�z�?�-?���?���?��?��?�Q�?�9X?���?�?}?�ff?��H?��?z^5?z��?qhs?R�!?Z^5?���?�X?��-?��h?��m?�+?��m?��/?st�?j��?`�?^v�?W�P?I7L?H��?f�y?d�?_;d?Z�H?X��?X�u?W�P?Rn�?I��?F��?BM�?:�?1�?+C�?(r�?�R?V>���>�h>��>��F>���>�->ɺ^>ɺ^>��>�^5>�Ĝ>�I�>��>�%>t�j>gl�>q��>`A�>F��>:^5>6E�>/�>)��>(��>,1>0 �>aG�>k�>O�;>B�\>.{>#�
>+=��==�x�=�l�=�=��#=��=��#=�=�Q�=�1=�C�=8Q�=t�=o<���<�t�<e`B<D��;�`B��o�ě��D����t���������t��#�
�8Q�P�`�]/�}�u�]/�u�����������{��E��ě��ȴ9�����;d��l����%�+�C��O߾bN�hs�hs�t��z��P�����-��-�2-111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 �� ſ�bN��녿������`�� ſ�vɿޗ������A���vɿ�5?�����/���m��ƨ�ۅ���H��"ѿ�=q�ؓu�׍P�׮��
=��l����y��`B�Ѓ���y���u�������j��Z��Mӿ�\)���!�[�m�CS��A���A�7�@  �3t���C�>8Q�?~�R?~�R?��?��9?�{?��?��7?���?�1?�r�?ļj?ǍP?�I�?��F?���?�1'?��?���?��u?�ȴ?��T?��?��?���?� �?�v�?�p�?��?�Q�?��u?�K�?�$�?�z�?�-?���?���?��?��?�Q�?�9X?���?�?}?�ff?��H?��?z^5?z��?qhs?R�!?Z^5?���?�X?��-?��h?��m?�+?��m?��/?st�?j��?`�?^v�?W�P?I7L?H��?f�y?d�?_;d?Z�H?X��?X�u?W�P?Rn�?I��?F��?BM�?:�?1�?+C�?(r�?�R?V>���>�h>��>��F>���>�->ɺ^>ɺ^>��>�^5>�Ĝ>�I�>��>�%>t�j>gl�>q��>`A�>F��>:^5>6E�>/�>)��>(��>,1>0 �>aG�>k�>O�;>B�\>.{>#�
>+=��==�x�=�l�=�=��#=��=��#=�=�Q�=�1=�C�=8Q�=t�=o<���<�t�<e`B<D��;�`B��o�ě��D����t���������t��#�
�8Q�P�`�]/�}�u�]/�u�����������{��E��ě��ȴ9�����;d��l����%�+�C��O߾bN�hs�hs�t��z��P�����-��-�2-111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�B	�B	�B	� B	�B	�B	{�B	�B	�B	�B	�B	�B	�B	�B	�B	�%B	�1B	�7B	�1B	�DB	��B	�hB	�oB	�uB	�uB	��B	��B	��B	��B	�?B	�qB	ɺB	��B	��B	�B	�HB	�B	��B
  B
B
%B
�B
ffB
ɺB+B49BA�B^5BdZBdZBn�Bt�B�%B��B�VB��B�hB�DB�=B�+B�DB�PB�JB�PB�\B�\B�\B�\B�VB�VB�VB�PB�DB�JB�JB�JB�JB�JB�VB�\B��B��B��B�{B��B��B��B��B�JB�DB�=B�%B~�B�1B��B��B�B�B�B�B��B��B��B�{B�uB�{B�uB�\B�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�hB�VB�JB�1B�1B�=B�bB�bB�bB�bB�VB�=B�DB�=B�7B�=B�=B�=B�7B�7B�=B�DB�=B�=B�PB�PB�bB�oB�uB�hB�bB�VB�VB�PB�PB�PB�VB�VB�\B�bB�oB�uB�bB�\B�\B�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�bB�bB�bB�bB�hB�hB�hB�oB�hB�oB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B	�B	�B	�B	� B	�B	�B	{�B	�B	�B	�B	�B	�B	�B	�B	�B	�%B	�1B	�7B	�1B	�DB	��B	�hB	�oB	�uB	�uB	��B	��B	��B	��B	�?B	�qB	ɺB	��B	��B	�B	�HB	�B	��B
  B
B
%B
�B
ffB
ɺB+B49BA�B^5BdZBdZBn�Bt�B�%B��B�VB��B�hB�DB�=B�+B�DB�PB�JB�PB�\B�\B�\B�\B�VB�VB�VB�PB�DB�JB�JB�JB�JB�JB�VB�\B��B��B��B�{B��B��B��B��B�JB�DB�=B�%B~�B�1B��B��B�B�B�B�B��B��B��B�{B�uB�{B�uB�\B�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�hB�VB�JB�1B�1B�=B�bB�bB�bB�bB�VB�=B�DB�=B�7B�=B�=B�=B�7B�7B�=B�DB�=B�=B�PB�PB�bB�oB�uB�hB�bB�VB�VB�PB�PB�PB�VB�VB�\B�bB�oB�uB�bB�\B�\B�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�bB�bB�bB�bB�hB�hB�hB�oB�hB�oB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002024                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002303  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002303  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A!��Dl9�G�O�                