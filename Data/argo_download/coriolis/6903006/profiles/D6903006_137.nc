CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:25Z creation; 2023-01-04T00:23:14Z last update (coriolis COQC software)   
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
_FillValue                  8  v�Argo profile    3.1 1.2 19500101000000  20230104002025  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @������1   @������@Rg�Qtf��(�0�7nP8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     ?;N��   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �r��c�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �~io@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���ܺ�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��3�a�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��9D�`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���@H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���>2�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���N�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��Pt  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���lw�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���,_�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����l�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��
=p�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��x9��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���=�  0999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990  A$��A.ffA@  AQ��Aa��AnffA�  A�33A�33A���A�33A�  A�  A�  A�  A�  A���A�  A�ffA�33A�  A�  A�33B  B  BffB  B  B  B��B��B$ffB(  B+33B0  B4ffB8  B:��B>��BB��BHffBU33Bg33Bz��B�33B�ffB���B�  B�  B�ffBÙ�B͙�Bי�B���B뙚B���C �C  C
  C��C�3C�3C�3C"�fC'��C,� C1�fC7�C;��C@��CE��CJ� CO� CT� CY��C^�3Cc�3Ch��Cn  Cs33Cx�C|�fC��fC�L�C�  C�s3C��fC�s3C��fC�ffC��fC�ffC��fC�s3C�  C���C��3C�Y�C��fC�� C��C�s3C���C�ffC��C�� C��fC�ffC�ٚCČ�C��3C�ffC��CΌ�C��Cә�C��3C؀ C��C�s3C��fC�L�C�� C�s3C�33C�3C��3C�@ C���C���C��C�� C�  D FfD��D��D�fD9�D�fD� D��D
9�Dy�D��D�fD33D� D��D�3D33Ds3D�3D�3D33Dy�D� D  DFfD�fD ��D"�D#9�D$ffD%� D&��D(33D)s3D*�3D+��D-@ D.�fD/��D0��D2&fD3y�D4��D6  D733D8� D9��D;  D<33D=�fD>� D?��DA33DBl�DC�fDE  DFFfDG�fDH�fDJfDK@ DLs3DM�fDN�3DP33DQs3DR��DS�3DU9�DVs3DW�3DX��DZL�D[� D\�3D]�fD_9�D`��Da� Db�3Dd33Del�Df��DhfDi9�Djy�Dk��Dm  DnL�Do� Dp�3Dr  Dt�31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A$��A.ffA@  AQ��Aa��AnffA�  A�33A�33A���A�33A�  A�  A�  A�  A�  A���A�  A�ffA�33A�  A�  A�33B  B  BffB  B  B  B��B��B$ffB(  B+33B0  B4ffB8  B:��B>��BB��BHffBU33Bg33Bz��B�33B�ffB���B�  B�  B�ffBÙ�B͙�Bי�B���B뙚B���C �C  C
  C��C�3C�3C�3C"�fC'��C,� C1�fC7�C;��C@��CE��CJ� CO� CT� CY��C^�3Cc�3Ch��Cn  Cs33Cx�C|�fC��fC�L�C�  C�s3C��fC�s3C��fC�ffC��fC�ffC��fC�s3C�  C���C��3C�Y�C��fC�� C��C�s3C���C�ffC��C�� C��fC�ffC�ٚCČ�C��3C�ffC��CΌ�C��Cә�C��3C؀ C��C�s3C��fC�L�C�� C�s3C�33C�3C��3C�@ C���C���C��C�� C�  D FfD��D��D�fD9�D�fD� D��D
9�Dy�D��D�fD33D� D��D�3D33Ds3D�3D�3D33Dy�D� D  DFfD�fD ��D"�D#9�D$ffD%� D&��D(33D)s3D*�3D+��D-@ D.�fD/��D0��D2&fD3y�D4��D6  D733D8� D9��D;  D<33D=�fD>� D?��DA33DBl�DC�fDE  DFFfDG�fDH�fDJfDK@ DLs3DM�fDN�3DP33DQs3DR��DS�3DU9�DVs3DW�3DX��DZL�D[� D\�3D]�fD_9�D`��Da� Db�3Dd33Del�Df��DhfDi9�Djy�Dk��Dm  DnL�Do� Dp�3Dr  Dt�31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����˿��T��˿���˿�����˿�`B���������j�䛦��Z��9X��9X��9X������Ͽ��Ͽ�F��F��\�����%��;d��V��5?�܋D��j��j��I��܋D�ܬ��V�ݲ-��{�ޗ�����������1��Kǿ�l���+��ȴ��-�θR��$ݿ���9X�����r��ixվE��?�1'?�n�?� �?��y?��?�1'?��w?�?}?�z�?�|�?���?�+?���?��?�-?�
=?��?���?�K�?�?�
=?��P?�$�?�o?���?�?��+?���?��?��#?��y?��?��H?��?�X?�E�?��j?���?��R?�ƨ?�+?��?��?�&�?z�H?wK�?q�?m��?j��?j=q?ix�?d�/?c��?d�/?gl�?a%?[dZ?F$�?%��?�+?/?�?�9?K�?O�?�
>��y>��>և+>�5?>�5?>��;>�=q>�=q>Ǯ>�J>���>�?}>��F>��y>���>���>�M�>�/>�t�>�I�>�$�>��>{�m>r�!>m�h>ix�>^5?>Z�>L��>1&�>+>!��>��>n�>�=��=�
==��=ě�=�Q�=���=�hs=�O�=�t�=�1=���=�t�=aG�=L��=�w=+<���<�1<D��;��
;D��;o;o;o:�o    �ě��#�
�49X�D�������`B���L�ͽe`B�m�h��C���hs�������罧���^5�Ƨ���ͽ��`��G���xս����h��   �$ݾ	7L�V�n������-�$�/�+�.{�2-�49X�:^5�F��Kƨ�N��O�;�T��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ��˿��T��˿���˿�����˿�`B���������j�䛦��Z��9X��9X��9X������Ͽ��Ͽ�F��F��\�����%��;d��V��5?�܋D��j��j��I��܋D�ܬ��V�ݲ-��{�ޗ�����������1��Kǿ�l���+��ȴ��-�θR��$ݿ���9X�����r��ixվE��?�1'?�n�?� �?��y?��?�1'?��w?�?}?�z�?�|�?���?�+?���?��?�-?�
=?��?���?�K�?�?�
=?��P?�$�?�o?���?�?��+?���?��?��#?��y?��?��H?��?�X?�E�?��j?���?��R?�ƨ?�+?��?��?�&�?z�H?wK�?q�?m��?j��?j=q?ix�?d�/?c��?d�/?gl�?a%?[dZ?F$�?%��?�+?/?�?�9?K�?O�?�
>��y>��>և+>�5?>�5?>��;>�=q>�=q>Ǯ>�J>���>�?}>��F>��y>���>���>�M�>�/>�t�>�I�>�$�>��>{�m>r�!>m�h>ix�>^5?>Z�>L��>1&�>+>!��>��>n�>�=��=�
==��=ě�=�Q�=���=�hs=�O�=�t�=�1=���=�t�=aG�=L��=�w=+<���<�1<D��;��
;D��;o;o;o:�o    �ě��#�
�49X�D�������`B���L�ͽe`B�m�h��C���hs�������罧���^5�Ƨ���ͽ��`��G���xս����h��   �$ݾ	7L�V�n������-�$�/�+�.{�2-�49X�:^5�F��Kƨ�N��O�;�T��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBB��BBBBBB��BBǮBŢBǮB��B��B��B��B��B��B��B��B��B�BB�yB�BB1'BH�B\)Bo�B}�B�=B��B��B�B��B��B��B��B�)B�TB�B	  B	1'B	\)B	hsB	w�B	�=B	��B	��B	�dB	��B	�sB
�B
M�B
�B�B,B>wBL�BZBdZBe`Bp�Bs�Br�BffBdZBbNBbNBcTBiyBp�Bu�Bz�B|�B}�B� B�B�B~�B�%B�7B�VB�\B�VB�VB�JB�JB�=B�DB�DB�DB�PB�JB�JB�DB�JB�VB�PB�PB�JB�VB�VB�\B�bB�hB�uB�{B��B��B��B��B�\B�+B�+B�7B�=B�1B�VB�PB�VB�VB�=B�PB�bB�hB�oB�oB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BB��BBBBBB��BBǮBŢBǮB��B��B��B��B��B��B��B��B��B�BB�yB�BB1'BH�B\)Bo�B}�B�=B��B��B�B��B��B��B��B�)B�TB�B	  B	1'B	\)B	hsB	w�B	�=B	��B	��B	�dB	��B	�sB
�B
M�B
�B�B,B>wBL�BZBdZBe`Bp�Bs�Br�BffBdZBbNBbNBcTBiyBp�Bu�Bz�B|�B}�B� B�B�B~�B�%B�7B�VB�\B�VB�VB�JB�JB�=B�DB�DB�DB�PB�JB�JB�DB�JB�VB�PB�PB�JB�VB�VB�\B�bB�hB�uB�{B��B��B��B��B�\B�+B�+B�7B�=B�1B�VB�PB�VB�VB�=B�PB�bB�hB�oB�oB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002025                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002314  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002314  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A$��Dt�3G�O�                