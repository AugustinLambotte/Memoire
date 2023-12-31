CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:25Z creation; 2023-01-04T00:23:22Z last update (coriolis COQC software)   
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
_FillValue                  �  D�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        |  Ex   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        |  I�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  MP   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     |  N0   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  Q�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  U(   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  V   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Y�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  Zd   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  ]�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  a\   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  b<   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  e�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     |  f�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    v�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    v�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    v�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    v�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  v�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    v�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    v�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    v�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         w    HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         w   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        w   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    w   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  j   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    jT   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    nT   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    rT   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  vTArgo profile    3.1 1.2 19500101000000  20230104002025  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�%Zq�r1   @�%Zq�r@Q��ێ�t�+����8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��:��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����7   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���6�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���-��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���~K  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���Q)X  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��`T�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��~K�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���r�d  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��5��$  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���l��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���>2�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���Յ�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �� �Ř  9999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990 A&ffA1��AA��AQ��A`  AnffA���A�33A���A�33A�ffA���A���A���A�33A�ffA�33A���A�ffA�ffA�  A�33A�ffB��B  B��B  B��B  B33B ��B#��B(ffB,  B/��B4ffB8  B<  B?��BDffBHffBS33Bh  B~  B�  B���B���B���B���B�ffBÙ�B���B�ffB�ffB���B���B���C��C	��C�fC�fC�fC�fC#  C(  C-�C1��C6��C;�fCAffCE�fCJ��CO��CT�fCZ  C_33Cd�Ch��Cn  Cs33Cx�C}�C��C���C�&fC�� C�ٚC�ffC��C�� C��3C�s3C��3C���C��C�� C�  C�L�C�ٚC�s3C���C�ffC�  C�Y�C��3C���C��3C�Y�C�  CĦfC��C�ffC���C�s3C�&fCӦfC�33C�s3C�� C�L�C�ٚC�ffC��3C� C��C�s3C��fC�Y�C��3C���C��C���C��D 9�Dl�D��D�3D&fDy�D�3D��D
FfD�fD��D�3D9�D� D��D��D33Dl�D�fD�fD&fDl�D��D��D33Dy�D ��D!��D#@ D$�3D%��D'fD(S3D)� D*�3D+�3D-,�D.�fD/��D0��D2&fD3` D4��D63D7L�D8��D9��D;fD<L�D=��D>�3D@  DA,�DB� DC�3DD�fDF,�DGffDH��DJ�DK@ DLy�DM��DN��DP33DQs3DR�3DS�3DU33DVl�DW��DX�3DZ33D[y�D\�fD^fD_9�D`l�Da�3Db�3Dd9�De� Df�fDg�3Di&fDjy�Dk��Dm  DnY�Doff1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A&ffA1��AA��AQ��A`  AnffA���A�33A���A�33A�ffA���A���A���A�33A�ffA�33A���A�ffA�ffA�  A�33A�ffB��B  B��B  B��B  B33B ��B#��B(ffB,  B/��B4ffB8  B<  B?��BDffBHffBS33Bh  B~  B�  B���B���B���B���B�ffBÙ�B���B�ffB�ffB���B���B���C��C	��C�fC�fC�fC�fC#  C(  C-�C1��C6��C;�fCAffCE�fCJ��CO��CT�fCZ  C_33Cd�Ch��Cn  Cs33Cx�C}�C��C���C�&fC�� C�ٚC�ffC��C�� C��3C�s3C��3C���C��C�� C�  C�L�C�ٚC�s3C���C�ffC�  C�Y�C��3C���C��3C�Y�C�  CĦfC��C�ffC���C�s3C�&fCӦfC�33C�s3C�� C�L�C�ٚC�ffC��3C� C��C�s3C��fC�Y�C��3C���C��C���C��D 9�Dl�D��D�3D&fDy�D�3D��D
FfD�fD��D�3D9�D� D��D��D33Dl�D�fD�fD&fDl�D��D��D33Dy�D ��D!��D#@ D$�3D%��D'fD(S3D)� D*�3D+�3D-,�D.�fD/��D0��D2&fD3` D4��D63D7L�D8��D9��D;fD<L�D=��D>�3D@  DA,�DB� DC�3DD�fDF,�DGffDH��DJ�DK@ DLy�DM��DN��DP33DQs3DR�3DS�3DU33DVl�DW��DX�3DZ33D[y�D\�fD^fD_9�D`l�Da�3Db�3Dd9�De� Df�fDg�3Di&fDjy�Dk��Dm  DnY�Doff1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���������+��z῔9X���u���#���Ͽ�1��V���㿫���V���7��S������
=���D���ļj�Ƈ+��I���-��?}�֧���ա˿�z�Ұ!����-��%��  ���ۿ�O߿�I���1�ə����T�Õ���t����T����� ſ�p����T���Ϳ��^�ۿ/����&�q��=C�>���>��+<�1�49X>�$�>���>��H?'l�?<j?G+?E��?F$�?BM�?AG�?;��?;��?7K�?2n�?*��?&�y?%��?#�
?7
=?S��?Rn�?I7L?Hr�?AG�?F��?>��?2�!?*��?(r�?(�9?(1'?&��?   ?(�??�?V?�/>�>�{>��>�(�>�n�>��>ڟ�>��>�^5>u>�=�/=D��;�o���ͽH�9�D���@��49X�\)���0 ŽY��T����P�D����t����P���㽛�㽛�㽕���7L��o�L�ͽ�P��/��P��1��o�D���ě�<t�<o�o�ě���o�D����o:�o��o;�o<D��<T��<u<�o<�t�<�1<ě�<�`B<���<T��;�o;D��    ��o�D����o���
�#�
��o���
��1�ě��o�t��#�
�0 Ž8Q�]/�]/�}󶽑hs��t����P���
��-��j���������"ѽ�l�����h����ٽ��#���	7L�t���u��w�%�T�(�þ0 ž6E��7KǾ<j�A�7�D���G��Kƨ�P�`�R�Y��]/�`A��e`B�l�D�q���x���{�m��  ���7���\������������1'1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 �������+��z῔9X���u���#���Ͽ�1��V���㿫���V���7��S������
=���D���ļj�Ƈ+��I���-��?}�֧���ա˿�z�Ұ!����-��%��  ���ۿ�O߿�I���1�ə����T�Õ���t����T����� ſ�p����T���Ϳ��^�ۿ/����&�q��=C�>���>��+<�1�49X>�$�>���>��H?'l�?<j?G+?E��?F$�?BM�?AG�?;��?;��?7K�?2n�?*��?&�y?%��?#�
?7
=?S��?Rn�?I7L?Hr�?AG�?F��?>��?2�!?*��?(r�?(�9?(1'?&��?   ?(�??�?V?�/>�>�{>��>�(�>�n�>��>ڟ�>��>�^5>u>�=�/=D��;�o���ͽH�9�D���@��49X�\)���0 ŽY��T����P�D����t����P���㽛�㽛�㽕���7L��o�L�ͽ�P��/��P��1��o�D���ě�<t�<o�o�ě���o�D����o:�o��o;�o<D��<T��<u<�o<�t�<�1<ě�<�`B<���<T��;�o;D��    ��o�D����o���
�#�
��o���
��1�ě��o�t��#�
�0 Ž8Q�]/�]/�}󶽑hs��t����P���
��-��j���������"ѽ�l�����h����ٽ��#���	7L�t���u��w�%�T�(�þ0 ž6E��7KǾ<j�A�7�D���G��Kƨ�P�`�R�Y��]/�`A��e`B�l�D�q���x���{�m��  ���7���\������������1'1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBt�B�'B0!B�B  B�B[#B��B��B�BB�B7LBK�BgmB��BB5?BA�B]/B�+B��B�wB��B�B�BB�NB�B��B	\B	uB	�B	�B	%�B	'�B	(�B	,B	9XB	G�B	O�B	W
B	p�B	��B	�B	ɺB	�`B
B
,B
VB
{�B
�B
��B
�fB
�B
�yB
�B+B�B2-BG�BT�B]/B^5BcTBe`BffBiyBl�Bo�Bp�Bq�Br�Bt�Bw�B}�B�DB�DB�DB�PB�PB�bB�bB�PB�PB�\B�hB�uB�{B�{B��B��B��B�{B�uB�uB�hB�hB�bB�bB�uB��B��B�hB�VB�B� B|�Bz�Bz�Bx�By�By�Bz�B}�B� B~�B~�B~�B}�B� B|�B|�B|�B|�B|�B}�B~�B� B�B�B�%B�B�1B�=B�DB�DB�PB�PB�=B�JB�VB�PB�\B�\B�\B�hB�uB�{B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 Bt�B�'B0!B�B  B�B[#B��B��B�BB�B7LBK�BgmB��BB5?BA�B]/B�+B��B�wB��B�B�BB�NB�B��B	\B	uB	�B	�B	%�B	'�B	(�B	,B	9XB	G�B	O�B	W
B	p�B	��B	�B	ɺB	�`B
B
,B
VB
{�B
�B
��B
�fB
�B
�yB
�B+B�B2-BG�BT�B]/B^5BcTBe`BffBiyBl�Bo�Bp�Bq�Br�Bt�Bw�B}�B�DB�DB�DB�PB�PB�bB�bB�PB�PB�\B�hB�uB�{B�{B��B��B��B�{B�uB�uB�hB�hB�bB�bB�uB��B��B�hB�VB�B� B|�Bz�Bz�Bx�By�By�Bz�B}�B� B~�B~�B~�B}�B� B|�B|�B|�B|�B|�B}�B~�B� B�B�B�%B�B�1B�=B�DB�DB�PB�PB�=B�JB�VB�PB�\B�\B�\B�hB�uB�{B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002025                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002322  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002322  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A&ffDoffG�O�                