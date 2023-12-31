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
resolution        ?PbM���     �  G   TEMP         
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
resolution        ?PbM���     �  d�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  l�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  tl   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  vh   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  xd   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  z`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �D   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �(   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �l   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �l   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �l   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �l   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �0   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �T   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �p   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �4Argo profile    3.1 1.2 19500101000000  20210225042458  20210225042458  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125450                          2C  D   APEX                            6229                            120210                          846 @����� 1   @����� @Q�`A�7�4���`A�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @,��@�  @�  A   AffA@  A`  A�  A�  A�  A�  A�  A���A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`ffBg��Bo��Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^�C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dqy�Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|` B
&�B
&�B
&�B
&�B
&�B
&�B
%�B
&�B
%�B
&�B
&�B
&�B
%�B
&�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
,B
+B
>wB
�3B
�B�B;dBB�BL�BR�B]/BdZBgmBjBl�Bk�Bm�Bp�Bp�Bp�Bp�Bp�Bn�Bm�Bm�Bl�Bl�Bm�Bm�Bn�Bn�Bn�Bn�Bo�Bo�Bn�Bn�Bm�Bm�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bm�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bp�Bp�Bq�Bq�Br�Br�Bs�Br�Bs�Bs�Bt�Bt�Bv�Bw�Bx�Bx�By�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{�B{�B{�B{�B{�B{�B{�B|�B}�B}�B}�B}�B}�B~�B~�B~�B~�B~�B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�+B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�7B�7B�7B�7B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�DB�DB�DB�DB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�PB�PB�VB�VB�VB�VB�\B�\B�bB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@XA�@XA�@XA�@X1'@XA�@XA�@XA�@XA�@XA�@X1'@X �@X1'@XA�@X1'@XA�@XA�@XA�@X1'@X1'@X1'@X1'@X1'@X1'@X1'@X1'@X �@X �@W�w@WK�@Q�@(r�@l�@n�?��?�1?�t�?�V?��?�  ?�b?�J?��?� �?��P?���?��?�\)?�ƨ?�
=?�Ĝ?o\)?e�T?a��?Vȴ?O�?K�?H1'?DZ?@A�?;"�?5?}?2�!?2-?,1?(�9?�w?j?ȴ?��?�;?�D?1?C�?	7L?Z?%>��m>���>�E�>�9X>�&�>��>�1>�>�ff>�`B>���>��/>���>�5?>�b>�b>���>�hs>�hs>�bN>��;>���>�V>���>�I�>ě�>Ƨ�>�C�>���>�+>�|�>�j>�X>�K�>�>�>�ȴ>��F>� �>�1>��>��/>���>�x�>���>��/>��R>�"�>���>���>�b>��+>���>��>�n�>��>��>��>��>�hs>��>�hs>�bN>�bN>�V>�C�>�1'>�+>���>���>��\>��>{�m>vȴ>u>q��>o��>m�h>k�>fff>aG�>`A�>]/>\(�>Xb>T��>T��>S��>N�>I�^>G�>C��>A�7>=p�>8Q�>8Q�>49X>.{>)��>'�>%�T>"��>$�/> Ĝ>�R>�R>�>��>��>�u>�P>z�>bN>\)>O�>I�>C�>C�>1'>$�>�=��=���=��=�=�`B=�G�=�;d=�/=��=�
==���=��`=���=Ƨ�=Ƨ�=Ƨ�=��=�^5=�Q�=�-=���=��T=��T=��T=��T=��T=��T=���=��w=��-=��=�hs=�O�=�7L=��=�o=��=�o=�o=�o=��=�7L=�C�=�\)=�t�=��=�\)=�O�=�+=�C�=�O�=�C�=�O�=�\)=�\)=�O�=�7L=��=��=�o=�%=�o=�O�=�t�=�\)=��=y�#=y�#=�%=�o=�+=�\)=�O�=�O�=�O�=��=m�h=e`B=aG�=aG�=Y�=H�9=@�=D��=<j=8Q�=8Q�=8Q�=49X=8Q�=<j=D��=0 �=#�
=�w=��=#�
=,1=,1=,1='�='�='�='�=#�
=#�
=�w=�P=t�=+<��<�<�`B<���<�j<u<e`B<u<u<T��<T��<T��<D��<#�
<t�<o;ě�;�o;D��;o    ��o��o�ě��t��D����C�������`B���C��C��t��#�
�49X�@��P�`�T���Y��Y��aG��u�����t�����ě����`������/��G���G���G���l���xս�����ٽ����m���#���#�����#�������J�o���o�J�%�%���������+�
=q�I��n���+��+��+��+��P��u��u�������-��w� Ĝ�!���"��#�
�$�/�$�/�%�T�&�y�&�y�(�þ)��,1�/��0 ž2-�1&�2-�49X�49X�5?}�7KǾ?|�G��H�9�L�;O�;�Q녾T���W
=�]/�_;d�aG��cS��hr��o���t�j�w�پ{�m��  ���\�����������ƨ��I���\)��n���񪾓t����Ͼ����b��(���A����徥�T���þ��D�� ž��!��������X��X���۾����+��1'�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9��7L��7L��7L��7L��7L��7L��7L��7L��7L�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @,�e@�@��L@��LA`�A?�&A_�&A�&A��A��A��A��A���A��A��A��B��B��B��B��B'��B/��B7��B?��BG��BO��BW��B`d�Bg�$Bo�$Bw��B��B��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB�2xB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EB��EC��C��C��C��C	��C��C��C��C��C��C��C��C��C��C��C��C!��C#��C%��C'��C)��C+��C-��C/��C1��C3��C5��C7��C9��C;��C=��C?��CA��CC��CE��CG��CI��CK��CM��CO��CQ��CS��CU��CW��CY��C[��C^<C_��Ca��Cc��Ce��Cg��Ci��Ck��Cm��Co��Cq��Cs��Cu��Cw��Cy��C{��C}��C��C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C��C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV�DV��DW�DW��DX�DX��DY�DY��DZ�DZ��D[�D[��D\�D\��D]�D]��D^�D^��D_�D_��D`�D`��Da�Da��Db�Db��Dc�Dc��Dd�Dd��De�De��Df�Df��Dg�Dg��Dh�Dh��Di�Di��Dj�Dj��Dk�Dk��Dl�Dl��Dm�Dm��Dn�Dn��Do�Do��Dp�Dp��Dqy�Dq��Dr�Dr��Ds�Ds��Dt�Dt��Du�Du��Dv�Dv��Dw�Dw��Dx�Dx��Dy�Dy��Dz�Dz��D{�D|�D|_�B
&�B
&�B
&�B
&�B
&�B
&�B
%�B
&�B
%�B
&�B
&�B
&�B
%�B
&�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
&<B
,�B
0lB
Z+B
��B
�B!�B=~BE�BN�BX�Bb�BgRBi�BnNBojBn�Bo�Bq&BqBq�BrvBsBq�BoeBn]Bn�BnBnDBn7BoTBobBo�Bo�BpBo�Bo�BoGBo9Bn8Bn�Bn"BnWBn1Bn�Bn�Bn�Bn~Bo;Bo-Bn�Bn�Bn�Bo�Bo�Bo�Bp�BqBq�Bq`BsBr�Bt%BsFBs�BtBuBt�Bv�Bw�ByBx�By�BzBz�By�BztBz�B{:B{�B{.B{,B{B{Bz�B{�B|0B|?B|IB{�B|}B{�B|�B~B~YB~�B~KB}�B*B!B!B7BB�B�B�B��B�#B�B�B�B�+B�B�EB�`B�\B�-B�QB�B�JB�LB�XB�YB�,B�PB�9B�9B�:B�]B�]B�.B�CB�.B�QB�DB�!B�-B�]B�]B�:B�SB�>B�XB�bB�)B�WB�oB�WB�?B�?B�IB�B�VB�>B�*B�AB�9B�CB�;B�<B�WB�bB�@B�JB�>B�>B�3B�UB�KB�LB�oB�[B�ZB�]B�]B�RB�GB�JB�VB�JB�JB�WB�bB�VB�>B�?B�bB�bB�KB�cB�qB�[B�DB�DB�GB�HB�KB�dB�[B�^B��B�jB�jB�iB�iB�\B�FB�\B�TB�RB�IB�>B�IB�@B�BB�RB��B�qB��B�JB�ZB�sB�_B�bB�pB�|B��B��B�pB�B��B�hB�9B�VB��B��B��B�zB�cB�qB�eB�SB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�-B�CB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�#B� B��B��B��B��B��B��B�
B��B��B��B��B�B��B��B��B��B��B��B��B�B��B��B�
B�B��B��B��B��B��B�"B�"B��B�
B�
B�B�B��B�	B�B��B��B�FB�RB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@XA�@XA�@XA�@X1'@XA�@XA�@XA�@XA�@XA�@X1'@X �@X1'@XA�@X1'@XA�@XA�@XA�@X1'@X1'@X1'@X1'@X1'@X1'@X1'@X1'@X �@X �@W�w@WK�@Q�@(r�@l�@n�?��?�1?�t�?�V?��?�  ?�b?�J?��?� �?��P?���?��?�\)?�ƨ?�
=?�Ĝ?o\)?e�T?a��?Vȴ?O�?K�?H1'?DZ?@A�?;"�?5?}?2�!?2-?,1?(�9?�w?j?ȴ?��?�;?�D?1?C�?	7L?Z?%>��m>���>�E�>�9X>�&�>��>�1>�>�ff>�`B>���>��/>���>�5?>�b>�b>���>�hs>�hs>�bN>��;>���>�V>���>�I�>ě�>Ƨ�>�C�>���>�+>�|�>�j>�X>�K�>�>�>�ȴ>��F>� �>�1>��>��/>���>�x�>���>��/>��R>�"�>���>���>�b>��+>���>��>�n�>��>��>��>��>�hs>��>�hs>�bN>�bN>�V>�C�>�1'>�+>���>���>��\>��>{�m>vȴ>u>q��>o��>m�h>k�>fff>aG�>`A�>]/>\(�>Xb>T��>T��>S��>N�>I�^>G�>C��>A�7>=p�>8Q�>8Q�>49X>.{>)��>'�>%�T>"��>$�/> Ĝ>�R>�R>�>��>��>�u>�P>z�>bN>\)>O�>I�>C�>C�>1'>$�>�=��=���=��=�=�`B=�G�=�;d=�/=��=�
==���=��`=���=Ƨ�=Ƨ�=Ƨ�=��=�^5=�Q�=�-=���=��T=��T=��T=��T=��T=��T=���=��w=��-=��=�hs=�O�=�7L=��=�o=��=�o=�o=�o=��=�7L=�C�=�\)=�t�=��=�\)=�O�=�+=�C�=�O�=�C�=�O�=�\)=�\)=�O�=�7L=��=��=�o=�%=�o=�O�=�t�=�\)=��=y�#=y�#=�%=�o=�+=�\)=�O�=�O�=�O�=��=m�h=e`B=aG�=aG�=Y�=H�9=@�=D��=<j=8Q�=8Q�=8Q�=49X=8Q�=<j=D��=0 �=#�
=�w=��=#�
=,1=,1=,1='�='�='�='�=#�
=#�
=�w=�P=t�=+<��<�<�`B<���<�j<u<e`B<u<u<T��<T��<T��<D��<#�
<t�<o;ě�;�o;D��;o    ��o��o�ě��t��D����C�������`B���C��C��t��#�
�49X�@��P�`�T���Y��Y��aG��u�����t�����ě����`������/��G���G���G���l���xս�����ٽ����m���#���#�����#�������J�o���o�J�%�%���������+�
=q�I��n���+��+��+��+��P��u��u�������-��w� Ĝ�!���"��#�
�$�/�$�/�%�T�&�y�&�y�(�þ)��,1�/��0 ž2-�1&�2-�49X�49X�5?}�7KǾ?|�G��H�9�L�;O�;�Q녾T���W
=�]/�_;d�aG��cS��hr��o���t�j�w�پ{�m��  ���\�����������ƨ��I���\)��n���񪾓t����Ͼ����b��(���A����徥�T���þ��D�� ž��!��������X��X���۾����+��1'�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9��7L��7L��7L��7L��7L��7L��7L��7L��7L�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9�ȴ9��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'��1'1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
<#�
<#�e<#�s<#�<#�<#�<#�<#��<#�;<#�q<#�e<#��<#��<#�<#�
<#�j<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�k<#�<#�<$<9,�<�"<��A<��?<E�r<'>	<+{�<'b<;�/<9�]<*��<(`<.�B<*%;<+�,<'`U<$
#<$h<%C<&f�<(:�<,v9<&lR<$U5<'dl<%�=<$7�<$*n<$B_<$Rm<$�<$�z<$�<#ڷ<$�<$3<%��<$+z<$�B<$�<$M<$$�<#�<#��<#�t<$��<$'^<$�<#�U<#��<#�=<#�*<#�*<#�h<#�G<#��<#�h<#�<#��<#��<#�t<$�<#�<#�<#��<#�/<#�Y<#ף<#݉<#�V<#�Y<#�Y<$5v<#�z<#�p<#�f<#�<$<�<#��<#��<#��<#��<#�<#ز<#��<#�<#��<#�<$�<#�i<#ס<#��<#��<$�<#��<#�=<#��<#�n<#ۈ<#�H<#��<#ׇ<#�I<#�<#؈<#ؙ<#׮<#ע<#�~<#��<#�<#��<#�'<#�<#�@<#��<#�<#�W<#��<#��<#�{<#ׇ<#�G<#��<#��<#�'<#�<#�_<#׭<#��<#׮<#�f<#�<#�<#מ<#�<#�<#�9<#�<#��<#޸<#�)<#�<#�t<#�<#�q<#�
<#� <#��<#؀<#� <#��<#�<#؁<#ט<#ر<#��<#�c<#�\<#�C<#ף<#��<#ׇ<#�z<#�<#��<#�<#�&<#�<#�)<#ڵ<#�X<#�><#�><#�W<#׊<#��<#׆<#׋<#��<#�+<#��<#�
<#�<#�<#�'<#ל<#�D<#�4<#؛<#�
<#�
<#�<#�<#�<#� <#�[<#ם<#�-<#��<#��<#��<#��<#�u<#�_<#�%<#�<#�<#ה<#��<#ם<#ؒ<#�<#�Y<#�D<#ץ<#ڸ<#��<#ת<#�b<#�O<#י<#�
<#׃<#��<#��<#�
<#�X<#�p<#ן<#�'<#�D<#��<#�]<#�A<#�<#��<#�f<#ك<#ݼ<#�j<#�<#�<#޺<#�k<#�<#ׂ<#�<#�<#�R<#��<#�3<#�u<#�y<#�
<#�
<#�j<#ׄ<#ב<#؀<#�<#��<#ׁ<#�i<#��<#��<#�<#�<#�r<#�
<#�
<#�
<#�R<#�<#׃<#��<#ח<#�<#��<#ב<#�<#�1<#�n<#�<#׈<#�f<#�<#�<#�<#�
<#��<#�h<#׊<#׋<#��<#��<#׋<#׊<#��<#ך<#�<#�<#�P<#ۑ<#�g<#��<#�[<#�U<#�<#�<#��<#ޛ<#ޕ<#ۈ<#�.<#ה<#�z<#�<#�<#�<#��<#�<$�<$'p<#�<#�<#�C<#��<#�<#�<#�<#ג<#׌<#�	<#ޒ<#�z<#ף<#ׁ<#�
<#آ<#�j<#�S<#�<#ޔ<#�7<#ד<#�q<#�n<#ף<#ׁ<#�<#��<#�
<#�S<#�<#��<#�9<#�A<#�<<#�q<#�<#�
<#�<#�y<#�z<#�
<#ׇ<#��<#��<#��<#ׂ<#ׂ<#ׂ<#ׂ<#ב<#�<#�z<#�z<#�<#��<#ב<#��<#�<#ב<#ع<#�d<#�z<#��<#�
<#ׇ<#��<#�<#�<#׹<#�<<#�<#��<#�D<#��<#�<#��<#�#<#�*<#��<#��<#��<#�v<#ޥ<#��<#�<#۩<#�p<#�<#�N<#׹<#�`<#�d<#װ<#׃<#ו<#ۀ<#�Q<#�(<#�,<#�G<#�M<#�<#�<#�D<#�l<#�;<#�k<#ל<#�<$(<$$<#�t<#��<#�}<#�
<#�
<#�
<#�
<#�t<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�<#�|<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.0014286                                                                                                                                                                                                                                                   none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151405201909291514052019092915140920190929151416202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@,��@,��@,��G�O�G�O�G�O�G�O�D|` D|` D|` G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          