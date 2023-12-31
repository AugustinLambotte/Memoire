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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  E�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  Lx   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  S$   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V|   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X(   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l,   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  m�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  q0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  w�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ~�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �4   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �<   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �X   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �|   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �$   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �@   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �\Argo profile    3.1 1.2 19500101000000  20210225043239  20210225043239  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125494                          2C  D   APEX                            6229                            120210                          846 @�JW/hK�1   @�JW/hK�@P��hr��4��S���1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�  @�  @���A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C-�fC0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D��D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9fD9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DC��DD� DE  DEy�DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DSfDS�fDT  DT� DUfDU�3DU�3B
r�B
r�B
r�B
q�B
q�B
r�B
s�B
t�B
r�B
t�B
u�B
w�B
y�B
z�B
{�B
~�B
�B
�B
�=B
�bB
�hB
�oB
��B
�{B
��B
��B
��B
��B
��B
��B
��B
��B
�-B
��B
��B
�5B
�HB
�NB
�ZB
�ZB
�`B
�B
��B
��B
��B
��BB	7BPB\B\B\BhB{B�B�B�B�B"�B&�B&�B&�B&�B&�B&�B'�B,B/B2-B6FB8RB:^B<jB>wB@�BE�BI�BI�BI�BM�BS�BVBXB[#B^5BaHBcTBffBk�Bo�Bo�Bp�Br�Bu�Bw�Bx�Bz�B|�B�B�B�B�B�B�B�B�%B�%B�+B�1B�7B�7B�=B�=B�7B�=B�=B�=B�=B�=B�=B�=B�7B�1B�+B�%B�B�B� B� B� B� B~�B|�B{�Bz�Bz�By�By�Bx�Bx�Bx�Bx�Bv�Bu�Bu�Bt�Bt�Bs�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bv�Bw�Bw�Bw�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�Bz�Bz�Bz�Bz�B{�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B}�B}�B}�B~�B~�B~�B~�B~�B� B� B� B� B~�B� B� B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�+B�1B�1B�1B�1B�1B�1B�1B�7B�=B�=B�=B�JB�JB�PB�PB�VB�VB�VB�VB�VB�VB�VB�\B�bB�bB�bB�bB�bB�\B�\B�\B�bB�\B�\B�\B�\B�bB�bB�bB�bB�bB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?&�?hs?hs?�`?&�?�`?��?�?��?��?�D?	��?�?l�?l�?�?�?�?��?	��?	��?	��?	��?	��?	��?	��?	��?	��?	��?
=q?I�?�?hs?��?O�?r�?+?ff?�?�
? �>�/>�^5>��#>��H>��7>և+>�r�>�Q�?��>��>�dZ?   ?J?�?	��??�D?�`?&�?�?O�?O�?O�?��?{?�j?�u?dZ?!��?%�T?'�?*��?-�h?.�?333?1�?.�?-V?0bN?4��?7
=?8b?9��?;dZ?<�?=/?>v�?@�?BJ?B�\?Co?D��?G+?H��?Gl�?E��?D�??|�?>v�?<�?;��?;"�?;"�?:�?8��?7�P?6E�?4�j?2�?2-?0�`?/\)?-��?)7L?&$�?%�T?$Z?!��? Ĝ?��?�?b?�F?O�?�y?S�>�`B>���>�z�>���>�I�>�v�>�V>�"�>�t�>�I�>��9>���>~��>vȴ>p��>x��>G�>6E�>2->/�>#�
>��>�+>t�>I�>1'>J>J>   =��m=��#=��#>o=��>J>o>�>o>J>J>o>$�>+>+>o=�=��=�x�=�l�=�l�=�`B=�`B=�S�=�S�=�G�=�;d=�/=�
==���=��=��`=���=���=ě�=\=�v�=�E�=�9X=�9X=�9X=�E�=�j=�v�=��=\=\=\=ě�=ě�=ě�=ě�=\=\=��=�j=�^5=�Q�=�E�=�9X=��w=�hs=��=�hs=�\)=�hs=��P=��P=���=���=��P=�7L=��=q��=m�h=ix�=ix�=ix�=aG�=e`B=T��=H�9=8Q�=49X=8Q�=D��=@�=Y�=]/=H�9=8Q�=#�
=�P=\)=C�=+='�=0 �=,1=\)=<j=T��=]/=e`B=y�#=y�#=m�h=Y�=T��=D��=@�=e`B=aG�=aG�=aG�=aG�=]/=Y�=H�9=0 �=0 �='�=�w=��=t�=o<��<��<�`B<���<ě�<�1<�C�<e`B<T��<T��<D��<49X<49X<49X<49X<49X<49X<#�
<#�
<t�<t�;ě�;D����o�o�D���T���T����t���j������`B��h�������o���,1�49X�8Q�@��D���P�`�Y��aG��u��%��+��C���O߽�t�������w���
��1��1��9X��E���vɽȴ9������������/��l�������h���#�%�o�$ݾC��O߾n�����u����w�"��%�T�',1�333�<j�?|�A�7�G��J���R�Y��hr��t�j�~�۾�����1'���;�
=��b��;d��A���Ĝ��Ĝ��Ĝ��G���G���G���������������MӾ�MӾ�MӾ��徢�徢�徣S���S���S����
���
���
��Z���T��l����羬1���h�� ž�����-��33���j��?}��?}��?}���پ�X��p���  ��%�\��o��o1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�  @�  @���A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C-�fC0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D��D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9fD9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DC��DD� DE  DEy�DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DSfDS�fDT  DT� DUfDU�3DU�3B
r�B
r�B
r�B
q�B
q�B
r�B
s�B
tzB
sB
u!B
vDB
x4B
y�B
z�B
{�B
~�B
�B
��B
�B
�ZB
�kB
�vB
�xB
��B
��B
��B
��B
��B
��B
�>B
�@B
��B
�3B
�aB
��B
�gB
�kB
�B
�B
�B
��B
��B
��B
��B
�DB
�BdB�BGB�B�B�B BlB&BmB]B�B"�B'KB'3B&�B&�B&�B&�B&�B+JB.�B1B5B7�B9�B;�B>3B?�BE�BJ?BJBIBL�BS�BU�BW�BZ�B]�Ba6BcBfBk;BoyBo�Bp[Br5BuoBxBy/B{2B}�B�7B�VB�PB�2B�B�MB�]B�cB�eB�vB��B�]B�vB��B��B�B��B�NB��B��B�zB�bB��B�pB�B�cB�ZB��B�B��B�B�B��B�IB~�B}�B{�B{�Bz6BzBy�By<ByBx�ByBv�Bu�Bt�BuFBt4Bs�Bt�BuBt�BuBt�Bt�Bt�Bu�Bu�Bu�Bv�Bw�Bw�Bw�Bx�Bx�Bx�Bx�Bx�By�By�BzBz<By�B{Bz�Bz�Bz�B{�B|�B|�B|�B|�B|�B}B|�B}B|�B}B|�B}B|�B}B}B|�B}�B}�B}�B~�B~�B~�B~�B~�B�B�B� B�B B�	B�B�B�B�B�B�B�B�B�\B��B�%B�B� B��B�B��B�B�*B�cB�'B�TB�B�B�B�B�&B�B�AB�9B�EB�&B�B��B�-B��B�B�jB�dB�oB�XB�JB�>B�:B��B�B�IB��B��B��B�2B�5B�B�SB�|B��B�dB��B�`B��B�gB�bB�bB�cB�oB�pB��B��B�bB�yB�vB�jB�vB��B�nB�dB��B�|B�}B��B��B��B�wB�pB�|B�{B�pB�oB�oB�oB�pB�{B�pB�{B�tB��B��B��B��B��B��B�xB��B��B��B��B��B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B��B��B��B��B��B��B�RB�2B�B� B�B�B��B��B�HB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B��B��B��B��B��B��?&�?hs?hs?�`?&�?�`?��?�?��?��?�D?	��?�?l�?l�?�?�?�?��?	��?	��?	��?	��?	��?	��?	��?	��?	��?	��?
=q?I�?�?hs?��?O�?r�?+?ff?�?�
? �>�/>�^5>��#>��H>��7>և+>�r�>�Q�?��>��>�dZ?   ?J?�?	��??�D?�`?&�?�?O�?O�?O�?��?{?�j?�u?dZ?!��?%�T?'�?*��?-�h?.�?333?1�?.�?-V?0bN?4��?7
=?8b?9��?;dZ?<�?=/?>v�?@�?BJ?B�\?Co?D��?G+?H��?Gl�?E��?D�??|�?>v�?<�?;��?;"�?;"�?:�?8��?7�P?6E�?4�j?2�?2-?0�`?/\)?-��?)7L?&$�?%�T?$Z?!��? Ĝ?��?�?b?�F?O�?�y?S�>�`B>���>�z�>���>�I�>�v�>�V>�"�>�t�>�I�>��9>���>~��>vȴ>p��>x��>G�>6E�>2->/�>#�
>��>�+>t�>I�>1'>J>J>   =��m=��#=��#>o=��>J>o>�>o>J>J>o>$�>+>+>o=�=��=�x�=�l�=�l�=�`B=�`B=�S�=�S�=�G�=�;d=�/=�
==���=��=��`=���=���=ě�=\=�v�=�E�=�9X=�9X=�9X=�E�=�j=�v�=��=\=\=\=ě�=ě�=ě�=ě�=\=\=��=�j=�^5=�Q�=�E�=�9X=��w=�hs=��=�hs=�\)=�hs=��P=��P=���=���=��P=�7L=��=q��=m�h=ix�=ix�=ix�=aG�=e`B=T��=H�9=8Q�=49X=8Q�=D��=@�=Y�=]/=H�9=8Q�=#�
=�P=\)=C�=+='�=0 �=,1=\)=<j=T��=]/=e`B=y�#=y�#=m�h=Y�=T��=D��=@�=e`B=aG�=aG�=aG�=aG�=]/=Y�=H�9=0 �=0 �='�=�w=��=t�=o<��<��<�`B<���<ě�<�1<�C�<e`B<T��<T��<D��<49X<49X<49X<49X<49X<49X<#�
<#�
<t�<t�;ě�;D����o�o�D���T���T����t���j������`B��h�������o���,1�49X�8Q�@��D���P�`�Y��aG��u��%��+��C���O߽�t�������w���
��1��1��9X��E���vɽȴ9������������/��l�������h���#�%�o�$ݾC��O߾n�����u����w�"��%�T�',1�333�<j�?|�A�7�G��J���R�Y��hr��t�j�~�۾�����1'���;�
=��b��;d��A���Ĝ��Ĝ��Ĝ��G���G���G���������������MӾ�MӾ�MӾ��徢�徢�徣S���S���S����
���
���
��Z���T��l����羬1���h�� ž�����-��33���j��?}��?}��?}���پ�X��p���  ��%�\��o��o1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�z<#�<#�k<#�F<#ױ<#ׯ<#�}<#��<#��<#�(<$�<#��<#��<#�
<#�*<#�
<#�<#ށ<#�:<#�p<#�
<#�<#�t<#� <#�
<#�
<#�
<#�
<#��<#�:<$	�<#��<#�<$a�<$r?<#��<#�t<#�<#�<$Bj<,��<+��<#�U<#ْ<$0�<&�8<&�<%�5<$��<$p<#�<#��<#�<$�_<#�]<#�e<#�q<$Ny<#׹<#�<#�B<#�<#�<#�4<#ׇ<%�<$F�<$�<$�s<$QS<#�4<$�<$�<#�<$S�<#��<$�<#�<$!�<$f�<#��<#�z<#�,<#�#<#�<#�<#��<#�4<#�	<#�<<#��<#�z<$f<#��<#�<#�e<#�<$h�<#�0<#�<#�<#��<#�<#�+<#��<#�<#�<#�2<#�><#�Y<#�#<#�=<#�P<$nE<$D<#��<#� <$�<#�7<$��<#�<#��<$_+<$�0<$��<$PL<*��<%Ŏ<#׊<#�&<$?2<%+�<& /<&	M<$O<$/�<#�P<#�-<$1�<#��<#�\<#�<'��<$np<#ߊ<#��<$<$�<#��<#�4<#��<#�5<#��<#�<#��<#��<#�]<#�<#�o<#�_<#�/<#ו<#�G<#�W<#�}<#�
<#׏<#�v<#כ<#�<#��<#�a<#��<#�a<#א<#�<#׬<#�<#�f<#�<#�~<#ׇ<#ט<#�<#ת<#��<#�[<#��<#�G<#��<#ף<#�<#�2<#��<#�<#�
<#�}<#۟<#ט<#׆<#�~<#�<#�<#�v<#�
<#�<#�#<#�G<#�<#׏<#��<#׏<#׿<#�s<#��<$�<#�V<#�<#��<#�z<#׆<#��<#�<#؞<#�<#ؽ<#��<#�@<#��<#״<#�~<#�<#�<#�<#�(<#ݜ<#�~<#��<#׍<#��<#ڱ<#�8<#�<#׺<#�<#��<#�<#�~<#��<#ׇ<#�E<#��<#�)<#�|<#�<$ <#�p<#��<#�^<#�<#�<#�J<#�?<#ס<#�<#�Q<#��<#�i<#�
<#�
<#�
<#ׁ<#ף<#�c<#��<#�%<#ؖ<#�<#ט<#�
<#��<#�~<#�<#� <#�
<#�<<#�0<#ބ<#�I<#׶<#�
<#׀<#�}<#�
<#�
<#�
<#�
<#�
<#�x<#�
<#�y<#�<#ڔ<#��<#��<#�"<#�^<#ב<#�<#�<#��<#��<#��<#׉<#�~<#�<#ׂ<#׶<#�Z<#߱<#إ<#ג<#��<#ס<#�<#�<#�/<#�<<#��<#�K<#��<#ל<#�	<#׿<#�2<#�<#�<#�<#��<#׮<#��<#�<#��<#�M<#��<#�/<#�&<#�<#�<#ײ<#�/<#߷<#��<#ۛ<#�p<#�h<#�b<#�}<#ۂ<#ތ<#ۂ<#��<#�P<#�.<#��<#�<#�7<#�1<#ٗ<#�-<#�+<#��<#�<$>:<$7<$<<#��<#�I<$�<$�j<#��<$,�<#��<#׿<#�<#�<#�w<#�
<#�
<#׉<#�<#�<#׼<#�<#�
<#�w<#�
<#�<#�~<#�<#�<#�w<#�
<#�<#��<#�;<#ۡ<#�	<#�p<#�<#�<#�M<#ט<#��<#�<#׊<#�
<#�<#�<#ۧ<#�m<#�.<#�@<#��<#�r<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929163125201909291631252019092916312920190929163137202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�  @�  @�  G�O�G�O�G�O�G�O�DU�3DU�3DU�3G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          