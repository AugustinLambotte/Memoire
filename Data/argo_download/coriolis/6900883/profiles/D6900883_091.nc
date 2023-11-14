CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   W   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2022-08-04T07:30:44Z creation; 2022-08-25T13:29:27Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_050n      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    9�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    9�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :    REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :$   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :4   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  :<   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  :|   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  :�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        :�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;    DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;(   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;,   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;0   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     ;P   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     ;p   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    ;�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           ;�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    ;�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            ;�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           ;�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           ;�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    ;�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    ;�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ;�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    ;�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     �  <�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  ?�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        \  ?�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  AD   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        \  A�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  B�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     \  CP   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  D�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  F   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  F`   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  G�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  H   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  Ip   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  J�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  K$   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  L�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     \  L�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    Z�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    Z�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    Z�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    Z�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  Z�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    Z�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    [   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    [   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         [    HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         [$   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        [(   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    [,   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  N4   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Nt   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Rt   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Vt   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  ZtArgo profile    3.1 1.2 19500101000000  20220804073044  20220825132927  6900883 AWI                                                             Gerd ROHARDT                                                    MTIME           PRES            TEMP            PSAL               [A   IF                                  2C  D   NEMO                            220                             23-May-2012                     860 @�%V""")1   @�%V""")@R�^�C.��%�n(���8   NONE        A   A   F   Primary sampling: discrete []                                                                                                                                                                                                                                      ��9���  ����5�  ��韫�  ���eC�  �����   ����Y�  ����   ����   ���5  ��ZC��  ������  �����`  ��"�9�  �����   ���%�p  ����c  ��Eȡ@  ���w�  ���a   ��"��  ��J��   ��(3��  ���Y �  �����  ��k���  ��OC�  ��ODP  �����   ����@  ��Eȡ  ���[g   ��Vٲ  �� �.x  ����8�  ��HpB�  ����p  ���Ϩ  ��$��H  ���%+  ���9D�  ��'O  ����ˀ  ��t���  �����  ���RL�  ��Er�  ��+��  ����ex  ��b�H  ��r��  �����H  ��T2(  ���[@  ����  ��)���  �����L  ���kT�  ��+�d�  ��}��  ���tn�  ��+�{�  ���F�   �����  ���}�  ������  ��~�/$  ��Q�n  ��+<Mx  �����  ����$�  ���I2�  ������  ��m:�  ��io<  ����k�  ������  ��R��r  ��F��  �Ʊ4  ��a"�  ��"�H  ���}'�  �͎��  ��I��X  �Ѐ0��  ��[�ާ  ��7_1�  000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 A  A1��AC33Ak33A�ffA�  A�33A�ffA�33A�33AᙚA�33B33B  B  B��B"��B)��B/33B:ffB@  BW��Bk��B��B�  B�33B���B���B�ffB���B�33B�ffB�ffC 33C
��CL�C�3C(ffC1�fC<�3CF� CP��CZ�3Cd�fCn�Cw��C���C��C�s3C�33C�  C��C�Y�C�&fC�33C���C�L�C�33C��C��fC�@ C��fC�@ C��3C�33D  D	33DY�D�3D�3D"�D(Y�D.�fD:� DG�3D`��Dy� D�Y�D��3D�C3D�� D�I�D�ٚD�C3D�ɚD�C3D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A��A+34A<��Ad��A�33A���A�  A�33A�  A�  A�fgA�  A�33B
ffBffB33B!33B(  B-��B8��B>ffBV  Bj  B~  B�33B�ffB�  B�  B���B���B�ffBי�B뙙B���C
34C�gCL�C(  C1� C<L�CF�CPfgCZL�Cd� Cm�4CwfgC�Y�C��gC�@ C�  C���C�ٚC�&gC��3C�  C�Y�C��C�  C��gC��3C��Cҳ3C��C�� C�  D�fD	�D@ D��D��D"  D(@ D.l�D:�fDG��D`� Dy�fD�L�D��fD�6fD��3D�<�D���D�6fD��D�6fD�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����t���J��%���-���ۿ����
=>.{?B�\?�n�?���?� �?���?���?�r�?ޗ�?�n�?�j?��y?���?�ƨ?�1?�`B?Ͳ-?�?�t�?�
=@7L@(�@#�m@1�@<�@AX@>�y@7�@41@1�@/��@,j@&�+@"��@E�@@�@�D?��?�J?ɺ^?��y?���?��\?��R?�X?�r�?�5??��?��j?lI�?Gl�?9��?��>�G�>ݲ->�ff>���>8Q�>;dZ>}�>��`>W
==���>)��=�h��w�J�cS�����ڟ����#��9�9X����!�7�+ƨ�6�+�;�m�A%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��t���J��%���-���ۿ����
=>.{?B�\?�n�?���?� �?���?���?�r�?ޗ�?�n�?�j?��y?���?�ƨ?�1?�`B?Ͳ-?�?�t�?�
=@7L@(�@#�m@1�@<�@AX@>�y@7�@41@1�@/��@,j@&�+@"��@E�@@�@�D?��?�J?ɺ^?��y?���?��\?��R?�X?�r�?�5??��?��j?lI�?Gl�?9��?��>�G�>ݲ->�ff>���>8Q�>;dZ>}�>��`>W
==���>)��=�h��w�J�cS�����ڟ����#��9�9X����!�7�+ƨ�6�+�;�m�A%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�BE�B_;B�?B�hB�FB�oB�B�BB�{BhsBB�B�B��Bv�BbNBT�BE�B�BBhBPB33Br�B�
B#�BQ�BaHB�BBaHB��B�#B�B1B��B�B�5B�TB��B��B�PB��B��B�7B�B��B�B,B  BuB%B�mBƨB��B�}B�9B�B��B�BuB"�B#�BuBL�BR�BO�B!�B{B>wB8RB%B;dBA�B@�B:^BB�BH�BI�BN�BP�BS�BS�BXBXB^5444444444444444444444444444444444444444444444444444444444444444444444444444444444444444 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�444444444444444444444444444444444444444444444444444444444444444444444444444444444444444 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�MTIME           PRES            TEMP            PSAL            Not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Not applicable                                                                                                                                                                                                                                                  Surface pressure = 0.4 dbar                                                                                                                                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Not applicable                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Severe sensor malfunction, measurements are bad and not adjustable                                                                                                                                                                                              20220825132927202208251329272022082513292720220825132927IF  ARFMCODA050n                                                                20220804073044                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220804073243  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC5.8                                                                 20220804073243  QCF$                G�O�G�O�G�O�0000000000008000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20220825132927  IP  PSAL            A  D���G�O�                