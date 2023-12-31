CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:24Z creation; 2023-01-04T00:23:08Z last update (coriolis COQC software)   
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
resolution        5�7�       =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  D�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  E�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  I(   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  J   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N|   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  R   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  U�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Vx   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Z   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Z�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^t   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  b    PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  fp   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  gT   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    wX   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    w\   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    w`   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    wd   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  wh   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    w�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    w�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    w�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         w�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         w�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        w�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    w�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  j�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    k    SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    o    SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    s    SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  w Argo profile    3.1 1.2 19500101000000  20230104002024  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               ~A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�Zq�r1   @�Zq�r@R������$E�>��8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    ?8EȠ   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �m3�a�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �~�Q�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���-��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��1M��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��=Ѻ�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��(�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���w`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��N��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��"�P  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��c]�T  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���$�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��hK�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��A��v  90999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990 A!��A,��AA��AS33A`  Ap  A���A�  A���A�  A���A�33A�33A�  A�  A�ffA�ffA�  A�  A�33A���A�  B33B��B  B��B��BffB33B  B ��B#��B(ffB,  B0  B4ffB8  B<ffB?��BDffBG33BT  Bg��B{��B���B���B���B�  B�33B�  B�  B�ffB�ffB���B�  B���B���C�fC
  C��C�fC�fC��C"�fC'�fC,�fC1�3C7  C<�C@�3CE��CJ�fCP  CU�CZ  C_�Cd  Ch�fCn  Cs�Cx�C}  C�  C���C��C�ffC��3C�s3C���C�ffC��3C�� C��C�Y�C�ٚC�ffC�ٚC�Y�C�ٚC�ffC�ٚC�ffC��3C�� C�ٚC�ffC��3C�s3C�  C�ffC�  CɌ�C��fC΀ C��CӀ C�  C�s3C�  C�ffC�ٚC� C�  C�Y�C�ٚC�s3C���C�s3C��C�s3C�ٚC�� C��D 9�Ds3D��D�fD  DffD��D��D
FfD�fD�fDfD,�Ds3D� D��D9�D�fD�3D  DL�Dy�D�fD�3D9�Dl�D ��D"fD#9�D$l�D%��D'fD(9�D)s3D*�fD+��D-L�D.�fD/�fD1  D2@ D3��D4��D5��D7&fD8ffD9�3D:�3D<33D=y�D>� D@fDA9�DBffDC��DD�3DF,�DGs3DH�3DI�3DK9�DL� DM��DN��DP@ DQy�DR� DS�3DU&fDVy�DW��DX��DZ&fD[s3D\�fD^fD_9�D`y�Da��Db��Dd@ De��Df��Dg�fDi33Dj� Dk�3Dl�fDn9�Do��Dp� Dq�3Ds9�Du�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A!��A,��AA��AS33A`  Ap  A���A�  A���A�  A���A�33A�33A�  A�  A�ffA�ffA�  A�  A�33A���A�  B33B��B  B��B��BffB33B  B ��B#��B(ffB,  B0  B4ffB8  B<ffB?��BDffBG33BT  Bg��B{��B���B���B���B�  B�33B�  B�  B�ffB�ffB���B�  B���B���C�fC
  C��C�fC�fC��C"�fC'�fC,�fC1�3C7  C<�C@�3CE��CJ�fCP  CU�CZ  C_�Cd  Ch�fCn  Cs�Cx�C}  C�  C���C��C�ffC��3C�s3C���C�ffC��3C�� C��C�Y�C�ٚC�ffC�ٚC�Y�C�ٚC�ffC�ٚC�ffC��3C�� C�ٚC�ffC��3C�s3C�  C�ffC�  CɌ�C��fC΀ C��CӀ C�  C�s3C�  C�ffC�ٚC� C�  C�Y�C�ٚC�s3C���C�s3C��C�s3C�ٚC�� C��D 9�Ds3D��D�fD  DffD��D��D
FfD�fD�fDfD,�Ds3D� D��D9�D�fD�3D  DL�Dy�D�fD�3D9�Dl�D ��D"fD#9�D$l�D%��D'fD(9�D)s3D*�fD+��D-L�D.�fD/�fD1  D2@ D3��D4��D5��D7&fD8ffD9�3D:�3D<33D=y�D>� D@fDA9�DBffDC��DD�3DF,�DGs3DH�3DI�3DK9�DL� DM��DN��DP@ DQy�DR� DS�3DU&fDVy�DW��DX��DZ&fD[s3D\�fD^fD_9�D`y�Da��Db��Dd@ De��Df��Dg�fDi33Dj� Dk�3Dl�fDn9�Do��Dp� Dq�3Ds9�Du�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����%��G���&��G���G���&��%��Ĝ��׿�Ĝ��Ĝ��׿����׿�׿� ſ� ſ� ſ� ſ�  �� ſ�  ��  �� ſ� ſ� ſ�  ���;�߾w������ۿ�V��5?����{�ݑh�܋D��j�܋D�܋D�܋D��1��ƨ�ۥ���m��I����Ϳ�I��ۅ���ÿ����I���p��q���&$ݾ�1'�q���D��>��9>��j?J=q?�X?���?�ff?��?��T?�?��u?�ff?�b?�?}?�t�?�S�?�1'?��w?�1?�M�?�M�?�hs?���?�bN?�|�?�~�?�dZ?�x�?���?���?�7L?���?��u?�9X?�  ?�p�?�K�?�t�?��;?��?�Q�?�7L?��+?��F?��9?���?�ȴ?�E�?�?}?��?���?���?��9?�1'?�r�?��9?�r�?��?�E�?�S�?���?��R?��?���?��9?��u?�ȴ?��?�z�?��?��
?��!?��\?��7?�%?��`?|j?z�?y�?xQ�?s�F?o�;?f$�?_�w?^�R?^v�?W��?L��?F$�?B�\??|�?9X?6ȴ?7��?:^5?7��?/��?*=q?(r�?*��?#�
? A�?��?�F?��?�>��>�ȴ>�Z>��>ۥ�>ݲ->�\)>��>��>�C�>�J>�j>���>� �>��#>��H>�9X>��D>���>���>w��>W
=>I�^>9X>-V>"��>V>o=�=�x�=��=�{=�+=q��=T��='�=+<��<�/<���<u<D��<t�;ě�:�o�D���49X�e`B��o��/�o�+�\)�'P�`�q����7L��O߽�t����w���罶E���-��^511111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��%��G���&��G���G���&��%��Ĝ��׿�Ĝ��Ĝ��׿����׿�׿� ſ� ſ� ſ� ſ�  �� ſ�  ��  �� ſ� ſ� ſ�  ���;�߾w������ۿ�V��5?����{�ݑh�܋D��j�܋D�܋D�܋D��1��ƨ�ۥ���m��I����Ϳ�I��ۅ���ÿ����I���p��q���&$ݾ�1'�q���D��>��9>��j?J=q?�X?���?�ff?��?��T?�?��u?�ff?�b?�?}?�t�?�S�?�1'?��w?�1?�M�?�M�?�hs?���?�bN?�|�?�~�?�dZ?�x�?���?���?�7L?���?��u?�9X?�  ?�p�?�K�?�t�?��;?��?�Q�?�7L?��+?��F?��9?���?�ȴ?�E�?�?}?��?���?���?��9?�1'?�r�?��9?�r�?��?�E�?�S�?���?��R?��?���?��9?��u?�ȴ?��?�z�?��?��
?��!?��\?��7?�%?��`?|j?z�?y�?xQ�?s�F?o�;?f$�?_�w?^�R?^v�?W��?L��?F$�?B�\??|�?9X?6ȴ?7��?:^5?7��?/��?*=q?(r�?*��?#�
? A�?��?�F?��?�>��>�ȴ>�Z>��>ۥ�>ݲ->�\)>��>��>�C�>�J>�j>���>� �>��#>��H>�9X>��D>���>���>w��>W
=>I�^>9X>-V>"��>V>o=�=�x�=��=�{=�+=q��=T��='�=+<��<�/<���<u<D��<t�;ě�:�o�D���49X�e`B��o��/�o�+�\)�'P�`�q����7L��O߽�t����w���罶E���-��^511111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�B	!�B	 �B	�B	�B	�B	!�B	 �B	 �B	 �B	!�B	 �B	!�B	!�B	!�B	"�B	"�B	#�B	#�B	"�B	"�B	"�B	"�B	"�B	#�B	"�B	"�B	#�B	#�B	&�B	&�B	(�B	)�B	)�B	+B	,B	/B	1'B	1'B	1'B	1'B	49B	9XB	;dB	<jB	=qB	D�B	E�B	F�B	J�B	M�B	YB	�?B	�B
�B
;dB
YB
y�B
��B
�dB
�B!�BB�BH�BI�B1'B49B=qBP�BVBYB[#BaHBjBp�B�B�1B�=B�=B�DB�DB�=B�=B�=B�JB�7B�7B�PB�VB�\B�PB�DB�=B�+B�+B�B�%B�VB�\B�VB�hB�oB��B��B��B��B��B��B�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�oB�oB�hB�uB��B��B�{B�uB�{B�{B�{B�uB�hB�hB�bB�hB�bB�\B�VB�VB�\B�\B�\B�PB�bB�uB�uB�oB�oB�uB��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B	�B	!�B	 �B	�B	�B	�B	!�B	 �B	 �B	 �B	!�B	 �B	!�B	!�B	!�B	"�B	"�B	#�B	#�B	"�B	"�B	"�B	"�B	"�B	#�B	"�B	"�B	#�B	#�B	&�B	&�B	(�B	)�B	)�B	+B	,B	/B	1'B	1'B	1'B	1'B	49B	9XB	;dB	<jB	=qB	D�B	E�B	F�B	J�B	M�B	YB	�?B	�B
�B
;dB
YB
y�B
��B
�dB
�B!�BB�BH�BI�B1'B49B=qBP�BVBYB[#BaHBjBp�B�B�1B�=B�=B�DB�DB�=B�=B�=B�JB�7B�7B�PB�VB�\B�PB�DB�=B�+B�+B�B�%B�VB�\B�VB�hB�oB��B��B��B��B��B��B�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�oB�oB�hB�uB��B��B�{B�uB�{B�{B�{B�uB�hB�hB�bB�hB�bB�\B�VB�VB�\B�\B�\B�PB�bB�uB�uB�oB�oB�uB��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002024                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002308  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002308  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A!��Du�G�O�                