CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:24Z creation; 2023-01-04T00:23:13Z last update (coriolis COQC software)   
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
resolution        5�7�       =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  D�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  E�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  I   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  J    PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  Nl   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Q�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  U|   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  V`   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Y�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Z�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^T   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  a�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  fH   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  g,   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    w,   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    w0   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    w4   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    w8   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  w<   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    w|   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    w�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    w�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         w�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         w�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        w�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    w�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  j�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    j�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    n�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    r�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  v�Argo profile    3.1 1.2 19500101000000  20230104002024  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�ڈ���1   @�ڈ���@Rn�i���(����^8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     ��E�    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �s|�@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���=   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��,���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���	0  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��ۗX  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���Y�p  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����/�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��[fǀ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����b�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��^З�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��<M^p  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���v�L  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��EȢ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��[f�~  0999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990  A!��A0  AA��AP  A\��AnffA|��A�ffA�33A�  A���A�ffA�33A�  A���A���A���A�33A�ffA�33A���A�  B ffBffBffB33B��BffB  B33B ffB$ffB'��B,ffB0  B4  B8��B<  B@��BDffBHffBTffBhffB|  B�33B�33B�  B���B�  B���B�33B�  B�  B�33B�  B���B���C  C
  C�C�C�fC�3C"��C(  C-�C2�C7�C<�CA33CF�CJ��CO��CU  CZ  C_33Cd  Ch��Cm�3Cr��Cw��C|��C���C�L�C�ٚC�Y�C�ٚC�s3C�  C�ffC��fC�s3C�  C���C��C���C�&fC�s3C�ٚC�s3C��fC�Y�C��C���C�  C�� C�  C�� C�  CĀ C��CɌ�C�ٚC�ffC�  Cӌ�C��C؀ C��3C�ffC�ٚC�Y�C�ٚC�s3C�ٚC�s3C��C�s3C��fC�Y�C��fC�s3C��D 9�Dl�D��D��D,�Dl�D��D�3D
@ D�fD��D�3D33Ds3D� D��D33Dy�D�fD�3D@ Dy�D��D�fD  D` D � D!�fD#,�D$� D%� D&��D(33D)s3D*��D+�3D-@ D.��D/�fD1  D2@ D3� D4� D5�fD733D8� D9�fD;�D<@ D=s3D>�fD?�3DA@ DBs3DC�fDD��DFFfDGy�DH��DJ  DKFfDLy�DM� DN�3DP&fDQs3DR�fDT  DU9�DVy�DW� DYfDZ9�D[s3D\�fD^fD_@ D`� Da��Dc  Dd&fDes3Df��Dg�fDi9�Dj�fDk��Dl�3Dn&fDoy�Dp�3Dr�Dt� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A!��A0  AA��AP  A\��AnffA|��A�ffA�33A�  A���A�ffA�33A�  A���A���A���A�33A�ffA�33A���A�  B ffBffBffB33B��BffB  B33B ffB$ffB'��B,ffB0  B4  B8��B<  B@��BDffBHffBTffBhffB|  B�33B�33B�  B���B�  B���B�33B�  B�  B�33B�  B���B���C  C
  C�C�C�fC�3C"��C(  C-�C2�C7�C<�CA33CF�CJ��CO��CU  CZ  C_33Cd  Ch��Cm�3Cr��Cw��C|��C���C�L�C�ٚC�Y�C�ٚC�s3C�  C�ffC��fC�s3C�  C���C��C���C�&fC�s3C�ٚC�s3C��fC�Y�C��C���C�  C�� C�  C�� C�  CĀ C��CɌ�C�ٚC�ffC�  Cӌ�C��C؀ C��3C�ffC�ٚC�Y�C�ٚC�s3C�ٚC�s3C��C�s3C��fC�Y�C��fC�s3C��D 9�Dl�D��D��D,�Dl�D��D�3D
@ D�fD��D�3D33Ds3D� D��D33Dy�D�fD�3D@ Dy�D��D�fD  D` D � D!�fD#,�D$� D%� D&��D(33D)s3D*��D+�3D-@ D.��D/�fD1  D2@ D3� D4� D5�fD733D8� D9�fD;�D<@ D=s3D>�fD?�3DA@ DBs3DC�fDD��DFFfDGy�DH��DJ  DKFfDLy�DM� DN�3DP&fDQs3DR�fDT  DU9�DVy�DW� DYfDZ9�D[s3D\�fD^fD_@ D`� Da��Dc  Dd&fDes3Df��Dg�fDi9�Dj�fDk��Dl�3Dn&fDoy�Dp�3Dr�Dt� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����\)��\)��;d��\)��;d��;d��\)��;d����ޗ��ݑh�ݲ-��푿��m��~���b��ȴ�և+��ȴ���ٿ�Q�ٙ��ۅ�ۥ�ۥ��1��1��(���1��ƨ�����o��Ĝ��7L�öF��\)��푿�(���C���b�����Ĝ��hs���ۿ�녿��������������5��Z��
=>>v�>�ƨ>� �?
~�?��=<j?2�?}�?���?�;d?�ƨ?��`?�(�?���?�hs?�$�?�X?ԛ�?�S�?�o?��?�-?У�?��?Η�?�A�?�dZ?��#?ļj?�E�?��/?�-?��?��7?�"�?���?�X?�K�?�&�?�33?�Z?�o?�%?�\)?�V?���?��7?��;?�A�?�M�?� �?��-?�dZ?�=q?���?�?}?�Ĝ?���?~�R?qhs?st�?m��?]/?Q��?/\)?5�?3�F?;��?PbN?L�D?E��?E��?>�R?49X?#��?E�?&�?
~�?�?��>�Q�>�l�>�>��>޸R>ٙ�>�"�>�O�>ě�>��`>޸R>׍P>և+>�;d>�=q>���>�j>��#>�ȴ>�9X>��#>��>��T>�A�>�/>��>��>�"�>��+>�n�>���>r�!>fff>]/>H�9>:^5>/�>-V>/�> Ĝ>�>\)>C�>+>1'>+>+=�l�=��=Ƨ�=��w=���=�hs=u=T��=8Q�=t�<�/<��
<�t�<o;o���
��t����+�8Q�T���8Q�,1�49X�8Q�D���L�ͽP�`�ixս�%��+��C����P�������
���
���置{��E���^5�ȴ9������`B�����ٽ��#�J1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ��\)��\)��;d��\)��;d��;d��\)��;d����ޗ��ݑh�ݲ-��푿��m��~���b��ȴ�և+��ȴ���ٿ�Q�ٙ��ۅ�ۥ�ۥ��1��1��(���1��ƨ�����o��Ĝ��7L�öF��\)��푿�(���C���b�����Ĝ��hs���ۿ�녿��������������5��Z��
=>>v�>�ƨ>� �?
~�?��=<j?2�?}�?���?�;d?�ƨ?��`?�(�?���?�hs?�$�?�X?ԛ�?�S�?�o?��?�-?У�?��?Η�?�A�?�dZ?��#?ļj?�E�?��/?�-?��?��7?�"�?���?�X?�K�?�&�?�33?�Z?�o?�%?�\)?�V?���?��7?��;?�A�?�M�?� �?��-?�dZ?�=q?���?�?}?�Ĝ?���?~�R?qhs?st�?m��?]/?Q��?/\)?5�?3�F?;��?PbN?L�D?E��?E��?>�R?49X?#��?E�?&�?
~�?�?��>�Q�>�l�>�>��>޸R>ٙ�>�"�>�O�>ě�>��`>޸R>׍P>և+>�;d>�=q>���>�j>��#>�ȴ>�9X>��#>��>��T>�A�>�/>��>��>�"�>��+>�n�>���>r�!>fff>]/>H�9>:^5>/�>-V>/�> Ĝ>�>\)>C�>+>1'>+>+=�l�=��=Ƨ�=��w=���=�hs=u=T��=8Q�=t�<�/<��
<�t�<o;o���
��t����+�8Q�T���8Q�,1�49X�8Q�D���L�ͽP�`�ixս�%��+��C����P�������
���
���置{��E���^5�ȴ9������`B�����ٽ��#�J1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBZBZB]/B]/B]/B\)B[#B[#BZB^5B`BBaHBffBn�B}�B�3B��B��B�ZB��B��BBhB�B$�B/B6FB7LB9XBA�BVBl�B{�B�hB��B�jBɺB��B�/B�5B�B��B	uB	1'B	e`B	� B	��B	�XB	��B	��B
%�B
B�B
m�B
�B
��B
�9B
�^B
�'B
�mBuB2-BS�BbNBk�BjBgmBffBjBv�By�Bz�B}�B~�B� B�B�B�B�+B�1B�B�B�7B�7B�7B�7B�B|�B|�B~�B}�B{�B~�B� B�B�B�B�B�1B�7B�DB�JB�VB�\B�\B�\B�VB�\B�VB�DB�+B�B�B�B�B� B|�Bv�Bz�B}�B�B�DB�DB�DB�PB�PB�=B�%B�B�B�B�B�B�B�B�%B�1B�+B�+B�7B�7B�=B�VB�hB�oB�oB��B�uB�uB�oB�oB�uB�{B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BZBZB]/B]/B]/B\)B[#B[#BZB^5B`BBaHBffBn�B}�B�3B��B��B�ZB��B��BBhB�B$�B/B6FB7LB9XBA�BVBl�B{�B�hB��B�jBɺB��B�/B�5B�B��B	uB	1'B	e`B	� B	��B	�XB	��B	��B
%�B
B�B
m�B
�B
��B
�9B
�^B
�'B
�mBuB2-BS�BbNBk�BjBgmBffBjBv�By�Bz�B}�B~�B� B�B�B�B�+B�1B�B�B�7B�7B�7B�7B�B|�B|�B~�B}�B{�B~�B� B�B�B�B�B�1B�7B�DB�JB�VB�\B�\B�\B�VB�\B�VB�DB�+B�B�B�B�B� B|�Bv�Bz�B}�B�B�DB�DB�DB�PB�PB�=B�%B�B�B�B�B�B�B�B�%B�1B�+B�+B�7B�7B�=B�VB�hB�oB�oB��B�uB�uB�oB�oB�uB�{B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002024                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002313  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002313  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A!��Dt� G�O�                