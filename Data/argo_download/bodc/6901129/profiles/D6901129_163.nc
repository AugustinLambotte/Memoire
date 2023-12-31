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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  F�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  N�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Z�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  \�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  dx   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  lT   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  t0   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  v(   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  x    PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  z   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �`   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �|   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225041939  20210225041939  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125451                          2C  D   APEX                            6229                            120210                          846 @��~d���1   @��~d���@Qě��T�4�KƧ�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @333@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx�Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DKfDK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{fD{� B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
�B
�)B
�HB
�HB
�NB
�HB
�HB
�HB
�HB
�HB
�NB
�NB
�TB
�mB
�mB
�B
��B
�B
�B
��BBBJB�B%�B.B49B:^BM�BW
BXBZB[#BbNBl�Bn�Bn�Bn�Bm�Bl�BjBiyBl�Bl�Bk�Bk�Bk�Bl�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bo�Bo�Bp�Bq�Bq�Bq�Bq�Br�Br�Bs�Bt�Bt�Bt�Bu�Bv�Bv�Bw�Bw�Bw�Bw�Bw�Bx�By�Bz�Bz�Bz�B{�B{�B{�B{�B{�B|�B|�B}�B}�B~�B~�B~�B~�B~�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�1B�1B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�=B�=B�=B�=B�=B�=B�=B�DB�DB�=B�=B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�DB�JB�JB�JB�JB�JB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�VB�\B�\B�\B�bB�bB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@'�;@(  @(  @(  @'�@(  @(  @(  @( �@(b@(  @(  @(  @(b@(b@( �@(b@(  @'�@'�P@';d@'�@'�@&5?@#�
@!x�@ff@$�@�-@�j@��@j@�@(�@��@�^@Q�@�D@33@ff@��@b@�@��@33?��R?�K�?֧�?�A�?�o?�7L?���?�"�?��F?��?~�R?z^5?g+?Q�?NV?N��?N��?N��?NV?NV?;�m?$�/?#�
?"J?^5?9X?�?l�?��?�
?J? �?   >�p�>�X>�?}>��j>�&�>�{>�V>�~�>��>��T>��
>�;d>�"�>ڟ�>ٙ�>�
=>��>�\)>�O�>�=q>Ǯ>š�>��7>�Q�>�33>�33>��!>� �>�{>��h>��D>�>��>��>��T>�G�>��R>�5?>��->�/>��>�b>��>��>��`>�bN>��>�O�>�C�>��9>��>�$�>��>��\>~��>y�#>w��>s�F>n��>k�>l�D>e`B>bM�>`A�>\(�>Z�>Y�>V>T��>Q�>Kƨ>C��>B�\>@�><j>6E�>1&�>.{>+>'�>#�
> Ĝ>�R>��>�u>�+>z�>hs>I�>$�>o>%=��=��m=�=�=�F=�h=�x�=�=�l�=�`B=�S�=�/=���=��=��=ȴ9=Ƨ�=\=ě�=ě�=��=�v�=�j=�j=�j=�j=�^5=�Q�=�E�=�Q�=�Q�=�9X=�1=���=��T=��
=��T=��T=��
=���=��w=��w=��-=���=��=�t�=�t�=�\)=�O�=�O�=�C�=�C�=�7L=�7L=�+=�+=�+=�+=��=�o=�o=�o=�o=�%=}�=y�#=m�h=Y�=D��=@�=@�=8Q�=49X=#�
=\)=C�=C�=t�=��=�w='�=0 �=0 �=0 �=0 �=0 �=,1=��=\)<��<�/<�j<�9X<�9X<�1<�t�<�o<e`B<T��<T��<D��<#�
<#�
<o;��
;�o;�o;ě�;ě�;�`B;�o        �D����`B�o�t��D���u��t���t���1�ě��ě��ě���j��j��j��9X��j���ͼ��+�\)�#�
�'',1�0 Ž<j�@��@��@��D���D���<j�8Q�49X�'�w�,1�0 Ž,1�,1�',1�0 ŽH�9�T���]/�e`B�m�h�q���u��%�u�aG��ixսq���u�y�#�}󶽇+�}�u�u�ixսaG��e`B�ixսm�h�y�#��7L��\)��t���t����㽥�T����罴9X��9X��9X��^5�Ƨ�ě����`��"ѽ�"ѽ�"ѽ�/��G���S���xս�������m�%�o�o���+�O߾bN����+�������-��-��R�!���%�T�',1�2-�7KǾ8Q�9X�>vɾ?|�A�7�C���E�˾Kƨ�O�;�R�W
=�Y��]/�_;d�bMӾbMӾcS��e`B�j~��u�w�پy�#�{�m�}󶾀���o�������˾�����^��ƨ������`��t���zᾖ�+���P������R������Z���/��`B���y��r����þ�~������33��9X��?}������X���#���H��j��vɾ�|�\�������š˾�$ݾƧ��1'�ȴ9��7L�ɺ^��=q������I����;�V��\)��bN��bN��녾Ձ��b�ۥ�߾w��;d�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w��A���A���A���A���A���A���A���A���A���A���Ĝ��Ĝ��Ĝ��G�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @333@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx�Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DKfDK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{fD{� B
̠B
��B
��B
��B
˾B
��B
��B
˲B
��B
��B
��B
��B
˾B
��B
˿B
��B
��B
��B
�B
�B
�B
��B
ҞB
ԵB
��B
�\B
�~B
�B
��B
�fB
�kB
�"B
�B
�fB
��B
�jB
�
B
�B
��B
�B
�UB
��B
�B
�KBB�B�B�B*�B1�B7�BA�BP�BW�BY	B[
B^�Bf<Bm9Bn�Bn�Bn�Bm�Bl�Bm�Bm�Bl�Bl�Bl�Bl�Bm'BmUBm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bo�Bo�Bp�Bq�Bq�Bq�Bq�Br�Bs!BtBt�Bt�Bt�Bv'BwBwBxBxBxBx7Bx�ByPBy�Bz�B{B{B{�B|B|B|/B{�B}!B}^B~0B~BB	BCB.B:B�>B�4B�B�/B�3B�?B�JB�'B�2B�'B�MB�ZB�PB�,B�EB�SB�<B�B�nB�CB�7B�PB�8B�-B�CB�.B�FB�kB��B�/B�9B�SB�lB�aB�LB�KB�LB�WB�KB�AB�QB�QB�EB�EB�RB�jB�uB�RB�EB�DB�=B�XB�6B�EB�\B�PB�-B�NB�DB�EB�]B�hB�EB�QB�^B�HB�UB�1B�>B�UB�JB�KB�CB�BB�AB�LB�QB�PB�9B�EB�^B�uB�SB�]B�PB�8B�FB�QB�VB�UB�KB�WB�XB�qB�\B�QB�kB�`B�WB�bB�WB�bB�WB�bB�WB�YB�XB�eB�gB�]B�_B�aB�oB�oB�pB��B��B��B�vB�jB��B�vB��B��B�vB�hB�QB�TB�bB�XB�]B�uB�uB�uB�vB��B��B��B��B��B��B��B�|B��B��B��B��B��B�|B��B��B�}B��B��B��B��B�jB��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�vB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�	B��B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B�	B�B��B��B�B��B��B��B��B�	B��B��B��B��B��B��B��B��B��B��B��B�CB��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�.B�
B��B��B��B��B��B��B��B�.B�!B��B��B�B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@'�;@(  @(  @(  @'�@(  @(  @(  @( �@(b@(  @(  @(  @(b@(b@( �@(b@(  @'�@'�P@';d@'�@'�@&5?@#�
@!x�@ff@$�@�-@�j@��@j@�@(�@��@�^@Q�@�D@33@ff@��@b@�@��@33?��R?�K�?֧�?�A�?�o?�7L?���?�"�?��F?��?~�R?z^5?g+?Q�?NV?N��?N��?N��?NV?NV?;�m?$�/?#�
?"J?^5?9X?�?l�?��?�
?J? �?   >�p�>�X>�?}>��j>�&�>�{>�V>�~�>��>��T>��
>�;d>�"�>ڟ�>ٙ�>�
=>��>�\)>�O�>�=q>Ǯ>š�>��7>�Q�>�33>�33>��!>� �>�{>��h>��D>�>��>��>��T>�G�>��R>�5?>��->�/>��>�b>��>��>��`>�bN>��>�O�>�C�>��9>��>�$�>��>��\>~��>y�#>w��>s�F>n��>k�>l�D>e`B>bM�>`A�>\(�>Z�>Y�>V>T��>Q�>Kƨ>C��>B�\>@�><j>6E�>1&�>.{>+>'�>#�
> Ĝ>�R>��>�u>�+>z�>hs>I�>$�>o>%=��=��m=�=�=�F=�h=�x�=�=�l�=�`B=�S�=�/=���=��=��=ȴ9=Ƨ�=\=ě�=ě�=��=�v�=�j=�j=�j=�j=�^5=�Q�=�E�=�Q�=�Q�=�9X=�1=���=��T=��
=��T=��T=��
=���=��w=��w=��-=���=��=�t�=�t�=�\)=�O�=�O�=�C�=�C�=�7L=�7L=�+=�+=�+=�+=��=�o=�o=�o=�o=�%=}�=y�#=m�h=Y�=D��=@�=@�=8Q�=49X=#�
=\)=C�=C�=t�=��=�w='�=0 �=0 �=0 �=0 �=0 �=,1=��=\)<��<�/<�j<�9X<�9X<�1<�t�<�o<e`B<T��<T��<D��<#�
<#�
<o;��
;�o;�o;ě�;ě�;�`B;�o        �D����`B�o�t��D���u��t���t���1�ě��ě��ě���j��j��j��9X��j���ͼ��+�\)�#�
�'',1�0 Ž<j�@��@��@��D���D���<j�8Q�49X�'�w�,1�0 Ž,1�,1�',1�0 ŽH�9�T���]/�e`B�m�h�q���u��%�u�aG��ixսq���u�y�#�}󶽇+�}�u�u�ixսaG��e`B�ixսm�h�y�#��7L��\)��t���t����㽥�T����罴9X��9X��9X��^5�Ƨ�ě����`��"ѽ�"ѽ�"ѽ�/��G���S���xս�������m�%�o�o���+�O߾bN����+�������-��-��R�!���%�T�',1�2-�7KǾ8Q�9X�>vɾ?|�A�7�C���E�˾Kƨ�O�;�R�W
=�Y��]/�_;d�bMӾbMӾcS��e`B�j~��u�w�پy�#�{�m�}󶾀���o�������˾�����^��ƨ������`��t���zᾖ�+���P������R������Z���/��`B���y��r����þ�~������33��9X��?}������X���#���H��j��vɾ�|�\�������š˾�$ݾƧ��1'�ȴ9��7L�ɺ^��=q������I����;�V��\)��bN��bN��녾Ձ��b�ۥ�߾w��;d�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w�߾w��A���A���A���A���A���A���A���A���A���A���Ĝ��Ĝ��Ĝ��G�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#�<#�
<#�j<#�g<#�
<#�<#ء<#�l<#�t<#�
<#�
<#�q<#�<#�e<#י<#�}<#ר<#��<#��<#�y<#�.<$5�<&8�<&9K<'��<#�m<#�<$,�<#�h<#�l<#� <#��<$˿<$2�<$�I<)u<$�<-kX</�G<$!<$�<(4.<#��<*��<T�O<A!a<R~�<5�j</^�<-�$<K��<)�|<$2�<$�C<$��<-�"</z2<$3<#�p<#�
<#�
<#׌<#��<,��<1d�<#�4<#�<%m�<$��<%�l<$S_<#�<#�<#�I<#�#<#�<#�<#�(<#��<#��<#�{<#�<#�<#��<#��<#��<#��<#��<#�<#ג<#�<#��<#�8<#��<#�c<#��<#��<#�u<#��<$l<$�<#�<#�o<#�C<#�<#׏<#� <#ۓ<#�<#�<#��<#�N<#��<#��<#ׇ<#׵<#�<#�&<#�~<#�}<#�J<#�m<#�<#�<#��<#�t<#�:<#�B<#�8<#��<#�<#�<#�<#޹<#�<#ڞ<#�A<#�<#�N<#ؽ<#�S<#��<#מ<#�<#׫<#ۣ<#�<#�<#��<#�<#�<#��<#��<#ە<#�g<#ۀ<#ވ<#�g<#�n<#�f<#�L<#�<#�<#ۓ<#�.<#�<#ۚ<#�<#��<#�v<#ۑ<#�<#ו<#�6<#��<#�^<#ؖ<#ׇ<#י<#�L<#�G<#כ<#�<#۸<#�h<#ؼ<#�}<#�<#��<#׎<#ו<#�<#�<#�<#׫<#׆<#�v<#�m<#�
<#�<#�I<#ר<#��<#�{<#�}<#�<#�*<#�<#�i<#�<#�<#ט<#ۆ<#�~<#�<#�;<#�W<#�<#�{<#�<#�v<#�<#�~<#�
<#�<#�<#״<#�i<#�
<#�<#�<#׆<#ׇ<#י<#�n<#�r<#�	<#מ<#�<#��<#ף<#ޢ<#�N<#ל<#�<#��<#�U<#ט<#ء<#��<#�<#�
<#�
<#�
<#מ<#�N<#ۛ<#޳<#��<#� <#�a<#�<#ו<#�8<#�	<#��<#ׄ<#�
<#׉<#��<#�<#��<#�k<#�]<#�
<#��<#�
<#�`<#��<#�e<#�<#�`<#�D<#ע<#י<#�g<#�d<#�<#�<#�F<#�\<#�<#�
<#�x<#�
<#�
<#�m<#ׁ<#�)<#�<#�}<#�A<#��<#�~<#�<#�<#ט<#�<#ׄ<#�<#�<#�w<#�<#؈<#��<#ט<#�3<#ح<#�W<#�S<#�s<#�
<#�p<#�w<#׵<#�u<#�j<#�
<#��<#�;<#ׇ<#ח<#ڜ<#�~<#�E<#��<#��<#׏<#ׇ<#ן<#�<#�W<#��<#�
<#�o<#��<#�n<#ׇ<#ך<#ۏ<#�A<#ۏ<#��<#�<#�n<#�i<#�U<#ת<#�<#�<#�<#�L<#��<#�z<#�.<#�j<#�
<#�
<#׆<#��<#ם<#�k<#�g<#�0<#�T<#�J<#��<#�<#ג<#۠<#�<#��<#�D<#�b<#�e<#۫<#�Z<#�<#ו<#�l<#�j<#��<#߇<#�A<#��<#׬<#׫<#�;<#מ<#��<#ؑ<#ٝ<#�<#�5<#�<#�p<#�,<#�n<#�+<#�<#�<#ב<#�)<#�<$q<#�;<#��<#��<#�<#ۋ<#��<#�t<#�2<#ޓ<#��<#�<#��<#�6<#�<#�<#�_<#�Z<#�k<#��<#�<#�<#׬<#י<#�i<#�8<#׮<#��<#�U<#��<#�a<#�V<#�<#ׇ<#א<#��<#�n<#ޅ<#�S<#�;<#۞<#��<#ׇ<#׈<#י<#�<#׎<#ׇ<#ׇ<#׈<#י<#�9<#צ<#�<#�<#��<#�<#ۑ<#�<#��<#�<#�<#�:<#�p<#�
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
<#�
<#�{<#�
<#�
<#ׅ<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151447201909291514472019092915145120190929151458202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@333@333@333G�O�G�O�G�O�G�O�D{� D{� D{� G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          