CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS     N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
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
resolution        ?�������     <  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  G\   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  O�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   W�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   Y�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   [�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     <  ^   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  f@   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  n|   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                   v�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                   x�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                   z�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     <  |�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  �$   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     <  �`   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �P   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �l   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �    HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �p   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225042951  20210225042951  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125452                          2C  D   APEX                            6229                            120210                          846 @���i 1   @���i @Qn��P�4�n��O�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            %A   A   A   @@  @�  @�  A   AffA@  A`  A�  A�  A�  A�  A�33A�  A�  A�33B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Ck�fCm�fCp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D�� D�  D�@ D�� D�ɚB
�'B
�!B
�!B
�!B
�!B
�!B
�'B
�!B
�!B
�'B
�-B
�'B
�3B
�-B
�3B
�?B
�LB
�RB
�RB
�XB
�XB
�dB
�dB
�qB
�}B
��B
ÖB
ĜB
ĜB
ĜB
ĜB
ĜB
ŢB
ŢB
ŢB
ŢB
ŢB
ŢB
ƨB
ǮB
ǮB
ǮB
ɺB
ɺB
��B
��B
��B
��B
��B
��B
�
B
�B
�B
�#B
�;B
�TB
��BhB!�B,B5?B[#BdZBe`Be`BdZBdZBgmBl�Br�Bs�Bv�Bx�Bz�B}�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�+B�+B�%B�+B�+B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�7B�1B�1B�1B�+B�+B�B�B�B�B�%B�+B�%B�B�B�B� B~�B� B~�B}�B}�B~�B~�B� B�B�B�%B�1B�1B�1B�1B�1B�1B�1B�1B�1B�7B�7B�1B�1B�7B�=B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�+B�+B�+B�1B�1B�1B�1B�1B�7B�7B�7B�7B�7B�7B�7B�7B�7B�=B�=B�=B�=B�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�JB�PB�PB�VB�\B�\B�\B�bB�bB�bB�bB�hB�hB�hB�oB�oB�uB�uB�uB�uB�uB�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@j@Z@Z@z�@�D@z�@z�@�j@�@��@z�@z�@�F@�m@�F@o@
��@
��@
�!@
��@
n�@
n�@
n�@
�\@
=q@	hs@	�@
M�@
^5@
�\@
�!@
��@
�H@
�H@
�@
�H@@C�@t�@33@@
��@	��@	hs@bN@�@ƨ?��w?�V?���?��u?���?�z�?�t�?�&�?��?��?�9X?���?�?щ7?�^5?��R?�(�?��H?���?�G�?�$�?k?a�7?c��?Xb?PbN?O�?S�F?V�+?X��?W�P?W
=?W��?U?T�j?S�F?O�;?L1?DZ??|�?;dZ?>��?E`B?E��?A��?A%??|�?=�-?<�?<�?<j?;dZ?9�#?8��?8b?7�P?2�?1hs?.V?,��?*=q?'+?$�/?#o?p�?(�?�H?�?b?E�?�!?bN?�?ƨ?	��?	�^?r�?+?S�?�
?��?J? A�>�v�>�dZ>���>�E�>�&�>>�1>�x�>��T>�Ĝ>�A�>ݲ->�t�>�n�>�bN>���>���>��!>�K�>�J>\>�v�>��T>�
=>�hs>���>���>�+>ix�>e`B>dZ>fff>hr�>p��>p��>���>�(�>���>�
=>�n�>��`>�\)>�V>���>���>��^>���>��>w��>t�j>��>��>}�>x��>u>s�F>cS�>Y�>\(�>\(�>Y�>Xb>P�`>L��>E��>2->%�T>�->�>��>�P>�u>�u>t�>hs>bN>bN>hs>bN>bN>\)>C�>	7L>�>J>%=��#=��#=��#=���=�=�F==�x�=�;d=�;d=�;d=�"�=��=��`=���=Ƨ�=\=��=�j=�9X=�9X=�9X=�9X=�9X=�9X=� �=�1=��=��T=��
=��=�Q�=�Q�=�^5=�^5=�^5=�^5=�^5=�^5=�^5=�Q�=�E�=� �=�{=� �=�1=���=��-=�t�=�C�=�C�=�C�=�C�=�C�=�C�=�O�=�O�=�O�=�C�=�7L=�+=��=�o=�o=��=��=�o=�o=�%=�%=}�=u=u=u=q��=ix�=e`B=e`B=e`B=e`B=aG�=Y�=T��=T��=L��=@�=8Q�=0 �=,1=#�
=�w=��=�P=t�=\)=+<�/<���<���<���<���<���<���<�j<�j<�j<�9X<�9X<�9X<�t�<u<49X<#�
<t�<t�<t�<t�;ě�;��
;D��;o:�o��o��o�ě��t��D����t����㼴9X��9X��j��j���ͼ�`B��`B���o��P�0 Ž8Q�@��]/�aG��q����o��������7L��t����-���罰 Ž�E��������������������/��S���S���`B��l���h���m�J�o�����$ݾ1'�\)�bN�hs�n��hs�t��z������u��-��w�!���"��&�y�',1�-V�,1�.{�/��0 ž0 ž333�7KǾ9X�>vɾA�7�B�\�G��J���Kƨ�Kƨ�M��O�;�P�`�R�T���Xb�^5?�aG��bMӾbMӾcS��cS��fff�ixվj~��l�D�n���q���w�پ~�۾��7��J���˾��𾇮���9���^���^��I�����\)��n����Ͼ�zᾕ�����
=���u����������/���w���徣S���S���Z���/��`B��ff��l���l���r���~����D������!���j��E���E����پ�X��dZ��푾��۾�o�Ƨ���;�\)��\)���;��bN��bN��bN���`��hs��녾�n���녾�녾�녾�녾�녾�녾�hs��녾�녾�녾�hs��hs��hs��hs��hs��hs��hs��hs��hs��hs��hs���`���`���`���`���`���`���`���`11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @@  @�  @�  A   AffA@  A`  A�  A�  A�  A�  A�33A�  A�  A�33B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Ck�fCm�fCp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D�� D�  D�@ D�� D�ɚB
�;B
�#B
�B
�B
�.B
�#B
��B
�B
�>B
�dB
�5B
��B
�B
�\B
��B
�rB
�NB
�iB
�dB
�~B
�_B
�`B
�WB
��B
�B
�3B
�NB
ĎB
�zB
ĆB
đB
ĊB
ŢB
řB
ŮB
ŌB
�uB
ŃB
��B
��B
��B
�MB
�B
ʀB
��B
�7B
ΞB
��B
ԇB
�"B
�B
وB
څB
�B
�B
��B
�+B�B*rB3�BC�B_sBeUBe�Bf�Bf�Bh�Bm�BnVBrfBu�Bx8ByBzB}hB�B�;B�B��B�jB�CB�IB��B��B��B� B��B�_B��B�B��B�@B�hB�wB�HB�'B�?B�WB�oB�ZB�QB�JB� B�{B��B�~B��B��B��B��B�9B�wB�qB��B�eB��B��B��B�{B��B��B�<B�iB�nB��B�B�?B�}B��B�fB�{B�pB�rB��B�rB�pB�qB��B��B�FB�vB�!B�OB�dB��B��B��B��B�)B�B��B�]B�B��B��BB�\B��B~,B~B~�B~�B�B��B�.B��B�]B�rB��B�YB�WB�KB�VB�5B�xB��B�XB��B�SB�^B��B�^B�uB�^B�WB��B��B�B�1B�UB�AB��B�fB��B�B��B��B�6B�JB�=B�B�'B�dB�CB�8B�+B� B�7B�,B�=B�bB�LB�cB�XB�CB�gB�9B�8B�DB�DB�DB�QB�_B�wB�@B�>B�VB�LB�mB�XB�bB�WB�OB�^B�tB�DB�FB�GB�JB�KB�cB�cB�cB�WB�UB�3B��B�SB�OB�\B�_B�_B�bB�bB�eB�rB�vB��B�zB�gB��B��B��B��B��B�zB�{B�{B�~B�~B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�
B�B��B��B�B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B�!B�B�QB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@j@Z@Z@z�@�D@z�@z�@�j@�@��@z�@z�@�F@�m@�F@o@
��@
��@
�!@
��@
n�@
n�@
n�@
�\@
=q@	hs@	�@
M�@
^5@
�\@
�!@
��@
�H@
�H@
�@
�H@@C�@t�@33@@
��@	��@	hs@bN@�@ƨ?��w?�V?���?��u?���?�z�?�t�?�&�?��?��?�9X?���?�?щ7?�^5?��R?�(�?��H?���?�G�?�$�?k?a�7?c��?Xb?PbN?O�?S�F?V�+?X��?W�P?W
=?W��?U?T�j?S�F?O�;?L1?DZ??|�?;dZ?>��?E`B?E��?A��?A%??|�?=�-?<�?<�?<j?;dZ?9�#?8��?8b?7�P?2�?1hs?.V?,��?*=q?'+?$�/?#o?p�?(�?�H?�?b?E�?�!?bN?�?ƨ?	��?	�^?r�?+?S�?�
?��?J? A�>�v�>�dZ>���>�E�>�&�>>�1>�x�>��T>�Ĝ>�A�>ݲ->�t�>�n�>�bN>���>���>��!>�K�>�J>\>�v�>��T>�
=>�hs>���>���>�+>ix�>e`B>dZ>fff>hr�>p��>p��>���>�(�>���>�
=>�n�>��`>�\)>�V>���>���>��^>���>��>w��>t�j>��>��>}�>x��>u>s�F>cS�>Y�>\(�>\(�>Y�>Xb>P�`>L��>E��>2->%�T>�->�>��>�P>�u>�u>t�>hs>bN>bN>hs>bN>bN>\)>C�>	7L>�>J>%=��#=��#=��#=���=�=�F==�x�=�;d=�;d=�;d=�"�=��=��`=���=Ƨ�=\=��=�j=�9X=�9X=�9X=�9X=�9X=�9X=� �=�1=��=��T=��
=��=�Q�=�Q�=�^5=�^5=�^5=�^5=�^5=�^5=�^5=�Q�=�E�=� �=�{=� �=�1=���=��-=�t�=�C�=�C�=�C�=�C�=�C�=�C�=�O�=�O�=�O�=�C�=�7L=�+=��=�o=�o=��=��=�o=�o=�%=�%=}�=u=u=u=q��=ix�=e`B=e`B=e`B=e`B=aG�=Y�=T��=T��=L��=@�=8Q�=0 �=,1=#�
=�w=��=�P=t�=\)=+<�/<���<���<���<���<���<���<�j<�j<�j<�9X<�9X<�9X<�t�<u<49X<#�
<t�<t�<t�<t�;ě�;��
;D��;o:�o��o��o�ě��t��D����t����㼴9X��9X��j��j���ͼ�`B��`B���o��P�0 Ž8Q�@��]/�aG��q����o��������7L��t����-���罰 Ž�E��������������������/��S���S���`B��l���h���m�J�o�����$ݾ1'�\)�bN�hs�n��hs�t��z������u��-��w�!���"��&�y�',1�-V�,1�.{�/��0 ž0 ž333�7KǾ9X�>vɾA�7�B�\�G��J���Kƨ�Kƨ�M��O�;�P�`�R�T���Xb�^5?�aG��bMӾbMӾcS��cS��fff�ixվj~��l�D�n���q���w�پ~�۾��7��J���˾��𾇮���9���^���^��I�����\)��n����Ͼ�zᾕ�����
=���u����������/���w���徣S���S���Z���/��`B��ff��l���l���r���~����D������!���j��E���E����پ�X��dZ��푾��۾�o�Ƨ���;�\)��\)���;��bN��bN��bN���`��hs��녾�n���녾�녾�녾�녾�녾�녾�hs��녾�녾�녾�hs��hs��hs��hs��hs��hs��hs��hs��hs��hs��hs���`���`���`���`���`���`���`���`11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#�<#ر<#�~<#�^<#�
<#ޔ<#�(<#�(<#�A<#�#<$�<#�<#�<#�J<#�<#�
<#�G<#��<#��<#�<#�,<#׽<#�	<$`<#��<#��<#��<#��<#��<#׎<#�Q<#�<#�k<#�Y<#��<#��<#�q<#݉<#��<#�<$!d<#�<$K�<$�/<(z4<)�)<$�^<$�<$�/<$�e<#�i<#� <$p<$f�<$ <(�<-d<Vu<O�<��<1�|<$�4<$�<$�D<(�H<1uN<>�I<&U<#��<'Co<%n~<#��<$l�<$<#��<#�_<#؍<#��<#��<#�V<#�7<$B<$D�<%n<$�]<$B�<$2�<$�<#�<$6�<#ۥ<#�	<#�<#��<#�<#�	<#ޱ<#�d<#�Y<#�I<#��<$h�<#�d<$o<#� <$8<$�<#��<#�<$��<#��<#�S<#�<#�P<#�<$7�<#�b<#�u<$�<#�B<#��<#��<#�<$8�<#�+<#ס<#�{<#��<#�e<#��<#�!<#��<$�<#��<#�<#�h<#� <#��<#�Y<#�Z<$�E<#�<#��<%�P<%a�<#�Y<#��<$�><#ם<$ }<'��<%w�<$�<$"�<#�s<#�<%��<#��<#׌<#�8<#��<#�<#�<&��<%s�<#��<#��<#��<#��<#�H<#�<#�<#�<#�V<#�<$�B<#��<#�f<$g�<#��<#ۖ<#�<#ۈ<#�<$J�<$�<#�<#�
<#� <#��<#�<#ߪ<#�<$%<$�<#��<#��<#�$<#��<#�k<#�<#��<#��<#׀<#�
<#�m<#�r<#�
<#�s<#� <#�=<#އ<#ہ<#�|<#��<#�<#�
<#׀<#׆<#׏<#�<#��<#�<#�<#�<#��<#ר<#�<#�+<#�@<#�<#�f<#�
<#�<#�
<#�<#�<#�
<#�<#��<#��<#��<#׋<#�k<#ٖ<#�<#�<#ב<#�
<#�<#�<#�
<#�
<#�<#�M<#ט<#�^<#�g<#�:<#؎<#�K<#�/<#�<#�B<#�<#�
<#�
<#�<#�<#�U<#�
<#�
<#׼<#�X<#׆<#מ<#�Q<#�
<#�w<#�
<#�x<#�<#�v<#�<#׋<#��<#�<#�<#׏<#��<#׆<#�
<#�
<#�
<#ׇ<#�L<#�_<#�<#��<#�0<#�<#�<#�{<#��<#׏<#׆<#׆<#׆<#א<#�/<#�<#�x<#�
<#�
<#�
<#�
<#׈<#��<#�<#�
<#�z<#�
<#�<#�F<#�|<#ݯ<#ס<#�}<#�
<#�
<#�<#��<#כ<#��<#׌<#׏<#��<#�@<#�B<#ڼ<#۱<#�<#ח<#��<#�<#�w<#�<#��<#�<#�<#��<#�1<#�+<#�#<#��<#�R<#��<#��<#ށ<#�<#�~<#�<#�<#�<#�K<#�K<#ڄ<#�R<#��<#��<#�w<#��<#�<#��<#�f<#�<#׽<#��<#��<#�b<#��<#��<#��<#�<#׉<#�A<#�<#׆<#׆<#�v<#�`<#��<#׉<#�~<#�<#�c<#�<#�<#��<#ת<#�=<#��<#�H<#׏<#׈<#��<#�f<#�~<#�<#�a<#�{<#�Q<#�P<#�Y<#׬<#�X<#�d<#׎<#�<#��<#��<#ט<#��<#�<#ۣ<#�F<#�v<#ׅ<#�<#�w<#�<#�0<#�4<#ל<#��<#�<#ۯ<#�q<#��<#��<#��<#�<#؛<#�u<#�g<#؆<#�<#�<#۟<#۷<#��<#�j<#ע<#��<#ד<#�<#�8<#�<#�2<#�.<#�R<#�J<#ק<#�<#��<#׏<#א<#��<#��<#�<#��<#ޔ<#��<#�r<#��<#�<#�/<#�<#�6<#ۄ<#ސ<#ۓ<#�S<#��<#�F<$�<#�[<#�<#�~<#�~<#�
<#�
<#�<#ׇ<#׆<#�v<#�s<#�
<#�
<#�
<#�
<#�
<#�n<#�s<#�
<#�
<#�v<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�{<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151525201909291515252019092915152920190929151535202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@@  @@  @@  G�O�G�O�G�O�G�O�D�ɚD�ɚD�ɚG�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          