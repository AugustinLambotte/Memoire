CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       06-May-2016 13:52:41Zcreation      
references        (http://www.argodatamgt.org/Documentation   comment       bThis netCDF file is generated using BODC's argoReader and netCDF writer software (argo@bodc.ac.uk)     user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    <H   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    <X   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    <\   REFERENCE_DATE_TIME                	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    <`   DATE_CREATION                  	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    <p   DATE_UPDATE                	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    <�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    <�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  <�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  <�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  =   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        =H   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    =L   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    =P   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     =T   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    =t   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    =x   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     =|   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     =�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     =�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    =�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   axis      T      
_FillValue        A.�~       
resolution        >�E�vQ�        =�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    =�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�E�vQ�        =�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   	valid_min         �V�        	valid_max         @V�        axis      Y      
_FillValue        @�i�            =�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    	valid_min         �f�        	valid_max         @f�        axis      X      
_FillValue        @�i�            =�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    >   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    >   VERTICAL_SAMPLING_SCHEME                   	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    >   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ?   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ?   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        axis      Z      
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     D  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  Fd   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  M�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     D  Zh   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  a�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  h�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  p4   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  r   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  s�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     D  u�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  |�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     D  �8   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �|   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �0   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �L   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �h   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �`   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �P   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �l   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225043602  20210225043602  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125469                          2C  D   APEX                            6229                            120210                          846 @�
�����1   @�
�����@P��z�H�3�~��"�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            )A   A   A   @,��@�  @�  A   A   A@  A`  A�  A�  A�  A�33A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D��Dy�D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D-��D.y�D/  D/� D0  D0� D1  D1�fD2  D2� D3  D3� D4  D4� D5  D5� D6fD6� D7  D7� D8  D8� D9  D9� D:  D:� D;fD;� D;��D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DKy�DL  DL�fDM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DUfDU� DU��DVy�DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[�fD\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Bz�Bz�Bz�By�By�By�Bx�Bv�Bo�B_;BA�B(�BB�BA�B<jB;dB9XB6FB5?B49B33B2-B33B1'B1'B33B2-B33B6FB6FB6FB6FB6FB49B7LB8RB8RB:^B:^B:^B:^B<jBA�BB�BB�B@�B@�BE�BE�BE�BE�BF�BG�BG�BH�BJ�BJ�BI�BH�BH�BI�BI�BJ�BK�BK�BK�BK�BK�BK�BK�BK�BK�BL�BL�BL�BL�BL�BL�BK�BK�BK�BL�BL�BL�BM�BO�BP�BQ�BQ�BQ�BQ�BR�BR�BR�BS�BS�BS�BT�BT�BT�BT�BT�BVBVBVBVBVBVBVBVBW
BVBVBW
BXBYBZB[#B\)B\)B\)B]/B^5B_;B_;B`BB`BB`BBaHBaHBbNBcTBcTBdZBdZBdZBffBgmBhsBiyBjBm�Bq�Bs�Bu�Bv�Bw�By�By�Bz�Bz�B{�B{�B|�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�+B�1B�7B�=B�=B�DB�JB�JB�JB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�bB�bB�bB�bB�hB�oB�oB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�I�@�A�@�9X@�(�@��m@��@�  @��F@�~�@�{@~E�@;dZ?У�?��m?��?��?�ƨ?�z�?o�?]�?S33?C�
?49X?/��?'+?/?�?�?�?-?��?\)?�?�?��? Ĝ>�j>��>��>>�~�>�ff>��>�b>�O�>�+>�x�>{�m>dZ>W
=>Kƨ>0 �>��=�;d=��T=���=��T=L��=\)<�`B<���<�/<��<��<��<��=o=o=C�=+=C�=\)=o=o<���<���<���<ě�<ě�<�1<���<���<�t�<�t�<u<e`B<�C�<�C�<�t�<�t�<���<���<��
<��
<��
<�1<�9X<�j<ě�<���<�/<�`B<�`B<�h<�`B<�9X<e`B<t�<t�<t�<e`B:�o�o�T��<#�
<�C�<���=C�=0 �=<j=@�=D��=L��=�7L=���=��-=��-=��-=��
=�{=�E�=�^5=\=Ƨ�=���=��=�/=�l�=�>$�>\)>�->5?}>F��>N�>_;d>bM�>ix�>n��>q��>s�F>w��>x��>y�#>��\>�o>�o>��>���>�$�>���>���>���>�+>�+>�+>��>�1'>�1'>��9>�7L>��9>���>���>�o>�%>�  >}�>w��>u>q��>m�h>k�>ix�>hr�>hr�>gl�>fff>e`B>e`B>dZ>cS�>cS�>aG�>`A�>^5?>^5?>]/>["�>W
=>T��>Q�>O�;>N�>M��>J��>I�^>H�9>E��>D��>A�7>>v�><j>49X>1&�>/�>.{>)��>$�/>"��>�R>�>��>��>�u>�>bN>V>	7L>+>J=���=�h=�;d=���=��=���=ȴ9=ě�=�j=�E�=� �=�{=�1=��
=���=��=�hs=�C�=�7L=�o=q��=aG�=]/=T��=@�=8Q�=8Q�=49X=0 �=,1='�=��=+<��<�<�/<�/<�/<���<�9X<��
<�1<�t�<e`B<49X<#�
<t�<o;�`B;�`B;�o    �D���D���o�49X�e`B��o��t��ě������o�C��C��C��C���P��w�'<j�@��Y��u�}󶽅���7L��t����㽡�����罴9X��j�\�ȴ9��
=��G���xս�������#�   �o�1'�V��+��-��R� Ĝ�$�/�,1�0 ž1&�333�333�6E��:^5�>vɾA�7�D���T���\(��["Ѿ["ѾZ��\(��]/�e`B�n���o���p�׾p�׾p�׾q���w�پvȴ�z�H�~�۾�o������$ݾ�+��������1'���9��=q���;�V��bN���`��hs��녾��Ͼ��u������"Ѿ�����-��5?���w��G���G����w��G���MӾ��徢�徣�
��Z���/���/���/���/���/��ff��l����þ�~���~��������1��{������-��?}��KǾ���������^5��j��p����۾�|��%���7��o��������$ݾƧ�Ƨ��+�Ǯ��1'��1'��1'��1'�ȴ9�ȴ9�ȴ9�ȴ9�ɺ^��V���`���`���`���Ͼ����և+�׍P��������"Ѿۥ�ܬ�ݲ-�ݲ-�ݲ-�ݲ-��/��/��/��/��/�ݲ-�ݲ-�ݲ-�ݲ-�ݲ-��5?�޸R��5?111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @(Q�@{�@�@�A�HA>�HA^�HA~�HA�p�A�p�A���A�p�A�p�A�p�A�p�A�p�B�RB�RB�RB�RB'�RB/�RB7�RB?�RBG�RBO�RBW�RB_�RBg�RBo�RBw�RB�RB��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)C�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C3�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��=C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
C��
D {�D ��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D�DuD��D{�D��D	{�D	��D
{�D
��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D {�D ��D!{�D!��D"{�D"��D#{�D#��D${�D$��D%{�D%��D&{�D&��D'{�D'��D({�D(��D){�D)��D*{�D*��D+{�D+��D,{�D,��D-{�D-�D.uD.��D/{�D/��D0{�D0��D1��D1��D2{�D2��D3{�D3��D4{�D4��D5{�D6�D6{�D6��D7{�D7��D8{�D8��D9{�D9��D:{�D;�D;{�D;�D<{�D<��D={�D=��D>{�D>��D?{�D?��D@{�D@��DA{�DA��DB{�DB��DC{�DC��DD{�DD��DE{�DE��DF{�DF��DG{�DG��DH{�DH��DI{�DI��DJ{�DJ��DKuDK��DL��DL��DM{�DM��DN{�DN��DO{�DO��DP{�DP��DQ{�DQ��DR{�DR��DS{�DS��DT{�DU�DU{�DU�DVuDV��DW{�DW��DX{�DX��DY{�DY��DZ{�DZ��D[��D[��D\{�D\��D]{�D]��D^{�D^��D_{�D_��D`{�D`��Da{�Da��Db{�Db��Dc{�Dc��Dd{�Dd��De{�De��Df{�Df��Dg{�Dg��Dh��Bz�Bz�Bz�Bz3By�By�ByHBx�BvtBv)BoPBdHBJHBGnBE�B>SB</B;B8�B6XB6B5B4#B2�B3B3TB2�B4�B6rB6�B6TB6_B6JB6ZB7�B8�B9SB:{B:�B:�B:�B=�BA�BC�BC;BC[BD�BF�BFFBF5BF�BG�BI�BIBH�BJ�BL7BJuBIBH�BI�BI�BJ�BK�BK�BK�BK�BK�BK�BK�BK�BK�BL�BM BL�BL�BL�BL�BK�BK�BK�BL�BL�BL�BM�BO�BP�BQ�BQ�BQ�BQ�BR�BR�BR�BS�BS�BS�BT�BT�BT�BT�BU BU�BVBVOBVgBVBBVBVBU�BW�BVtBV;BU�BW�BX�BY�BZ�B[�B\B\B]B]jB^�B_#B`AB`AB` Ba	BaBb4Bc!Bc=Bd5BdKBdBf#BgBg�BiBi�BlyBp�BsMBu Bv�Bw|By�By�Bz�Bz�B{�B{�B|rB�B� B��B��B� B� B�B�B�B�B�B� B�B�B�B�B�+B�iB�]B�@B�iB�UB�[B��B�dB�|B�~B�gB�iB�\B�UB�bB�cB�bB�WB�bB�bB�ZB�uB�jB�wB�aB�oB�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B� B��B��B��B�B��B��B��B��B��B��B��B��B��B�pB�B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�"B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�I�@�A�@�9X@�(�@��m@��@�  @��F@�~�@�{@~E�@;dZ?У�?��m?��?��?�ƨ?�z�?o�?]�?S33?C�
?49X?/��?'+?/?�?�?�?-?��?\)?�?�?��? Ĝ>�j>��>��>>�~�>�ff>��>�b>�O�>�+>�x�>{�m>dZ>W
=>Kƨ>0 �>��=�;d=��T=���=��T=L��=\)<�`B<���<�/<��<��<��<��=o=o=C�=+=C�=\)=o=o<���<���<���<ě�<ě�<�1<���<���<�t�<�t�<u<e`B<�C�<�C�<�t�<�t�<���<���<��
<��
<��
<�1<�9X<�j<ě�<���<�/<�`B<�`B<�h<�`B<�9X<e`B<t�<t�<t�<e`B:�o�o�T��<#�
<�C�<���=C�=0 �=<j=@�=D��=L��=�7L=���=��-=��-=��-=��
=�{=�E�=�^5=\=Ƨ�=���=��=�/=�l�=�>$�>\)>�->5?}>F��>N�>_;d>bM�>ix�>n��>q��>s�F>w��>x��>y�#>��\>�o>�o>��>���>�$�>���>���>���>�+>�+>�+>��>�1'>�1'>��9>�7L>��9>���>���>�o>�%>�  >}�>w��>u>q��>m�h>k�>ix�>hr�>hr�>gl�>fff>e`B>e`B>dZ>cS�>cS�>aG�>`A�>^5?>^5?>]/>["�>W
=>T��>Q�>O�;>N�>M��>J��>I�^>H�9>E��>D��>A�7>>v�><j>49X>1&�>/�>.{>)��>$�/>"��>�R>�>��>��>�u>�>bN>V>	7L>+>J=���=�h=�;d=���=��=���=ȴ9=ě�=�j=�E�=� �=�{=�1=��
=���=��=�hs=�C�=�7L=�o=q��=aG�=]/=T��=@�=8Q�=8Q�=49X=0 �=,1='�=��=+<��<�<�/<�/<�/<���<�9X<��
<�1<�t�<e`B<49X<#�
<t�<o;�`B;�`B;�o    �D���D���o�49X�e`B��o��t��ě������o�C��C��C��C���P��w�'<j�@��Y��u�}󶽅���7L��t����㽡�����罴9X��j�\�ȴ9��
=��G���xս�������#�   �o�1'�V��+��-��R� Ĝ�$�/�,1�0 ž1&�333�333�6E��:^5�>vɾA�7�D���T���\(��["Ѿ["ѾZ��\(��]/�e`B�n���o���p�׾p�׾p�׾q���w�پvȴ�z�H�~�۾�o������$ݾ�+��������1'���9��=q���;�V��bN���`��hs��녾��Ͼ��u������"Ѿ�����-��5?���w��G���G����w��G���MӾ��徢�徣�
��Z���/���/���/���/���/��ff��l����þ�~���~��������1��{������-��?}��KǾ���������^5��j��p����۾�|��%���7��o��������$ݾƧ�Ƨ��+�Ǯ��1'��1'��1'��1'�ȴ9�ȴ9�ȴ9�ȴ9�ɺ^��V���`���`���`���Ͼ����և+�׍P��������"Ѿۥ�ܬ�ݲ-�ݲ-�ݲ-�ݲ-��/��/��/��/��/�ݲ-�ݲ-�ݲ-�ݲ-�ݲ-��5?�޸R��5?111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�4<#�<#י<#�<#�<#�{<#�p<&z�<D��<Ț�=;v�=p�<L�#<<��<]�<*8�<)ѻ<4�N<+�}<',<*�<*�<$t�<%õ<&`�<#�X<$�<%I<#�<#��<#�<#׎<#�=<'0<$
�<$r<$�a<#��<#��<#�x<#��<$�<#��<$�<$%<)�.<0<$�-<$t<$�<%�<$��<&�<%&�<#�?<#�.<%e_<$4<#�y<#ٜ<#�u<#�{<#�z<#ס<#�u<#ح<#�}<#ڱ<#�<#ة<#؈<#�<#�S<#�l<#�Q<#ؗ<#׀<#�y<#�`<#ו<#�M<#�<#�P<#�><#�
<#�@<#׍<#�<#�w<#ذ<#�O<#�<#�r<#�t<#�	<#��<#�j<#��<#��<#��<#��<#�V<#ء<#�<#�C<#�|<#޿<#�^<#ב<#�u<$4<#��<#ܣ<$ɺ<#��<#��<$Z<$x<#��<#��<#��<#ܜ<$c�<#��<#��<#׈<#ׅ<#�f<#�E<#�H<#�h<#�<#��<#�g<#�9<#�A<#�<#��<$�<$�<$=�<$�)<$p�<$ �<$X�<#�<#�<#�S<#��<#�<#��<#�<#�;<$@<#�<#�i<#�<#�0<#ؼ<#ج<#�l<#�p<#؞<#�q<#�q<#ؼ<#��<#�I<#�<#ؘ<#�.<#��<#ۧ<#�<#ۋ<#׊<#�<#�r<#ױ<#۬<#�/<#�}<#ם<#�<#׀<#�<#�<#�<#�^<#�<#�<#�9<#מ<#�<#��<#ׁ<#�<#׸<#� <#��<#��<#מ<#�
<#�<#��<#�<#�<#ؗ<#�<#�"<#�-<#��<#�<#�<#ף<#�<#ۘ<#ߨ<#ף<#ې<#׬<#�<#�<#װ<#�Y<#�v<#ע<#ޘ<#׫<#�<#�}<#ߴ<#�-<#��<#�9<#�<#ף<#׽<#ۊ<#�O<#�<#�<#�<#ۃ<#ۮ<#�9<#׳<#�<#�<#ٷ<#ް<#ۅ<#�<#׾<#ޛ<#ח<#�V<#�<#�<#�<#�<#�B<#ޠ<#׵<#�<#��<#�R<#�^<#�<#ے<#�c<#�P<#� <#۰<#�'<#�<#�<#�<#�<#�K<#�<#۞<#�<#�9<#�c<#�a<#�.<#װ<#��<#� <#��<#�<#�<#�<#׈<#�f<#�d<#�l<#� <#�s<#��<#�<#�|<#��<#��<#�*<#��<#��<#��<#�V<#��<#�<#��<#�M<#ِ<#�/<#�~<#��<#۝<#��<#�y<#��<#ٟ<#�i<#�<#�<#�U<#�<#׻<#��<#�,<#ۤ<#�<#ב<#�D<#�2<#��<#ۻ<#�`<#�Q<$?<#�r<#�@<#�u<#�~<#צ<#�+<#�y<#�)<#�<#�<#�_<#�^<#�<#�<#�#<#�Y<#�+<#��<#�<#�"<#�q<#�<#�Z<#�<#�
<#��<#�<#ٛ<#��<#�<#�
<#�<#��<#��<#��<#�<#�&<#׭<#�<#�-<#�<#�l<#�$<#��<#ק<#�<#�Q<#ב<#�<#�<#�`<#�e<#�d<#�P<#�!<#��<#�5<#�<#�^<#�<#�<#�<#ۥ<#�~<#�K<#�7<#�^<#�^<#�B<#�%<#ۑ<#��<#�<#�<#�<#�<#�8<#׶<#ׯ<#ץ<#�<#�Z<#�<#�<#�<#�e<#�d<#�^<#�<#�b<#�d<#�T<#��<#��<#޼<#�K<#�A<#��<#��<#�2<#��<#�<#�=<#ې<#�<#ץ<#ל<#�^<#�d<#�l<#ؕ<#�i<#�d<#�d<#�]<#�<#�a<#�d<#�d<#�\<#�<#�<#�<#�};o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.07                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929155840201909291558402019092915584420190929155851202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@,��@,��@,��G�O�G�O�G�O�G�O�Dh� Dh� Dh� G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          