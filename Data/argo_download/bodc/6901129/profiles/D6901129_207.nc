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
_FillValue                  p  �\Argo profile    3.1 1.2 19500101000000  20210225043752  20210225043752  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125493                          2C  D   APEX                            6229                            120210                          846 @�G�3�J`1   @�G�3�J`@P�C���4�-V1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�ff@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B��B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS�fDTfDT� DU  DU�fDU� B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
$�B
%�B
$�B
$�B
%�B
%�B
%�B
%�B
%�B
%�B
&�B
(�B
.B
/B
0!B
1'B
2-B
1'B
0!B
/B
.B
/B
1'B
2-B
33B
33B
33B
49B
5?B
5?B
5?B
6FB
6FB
8RB
9XB
:^B
:^B
8RB
6FB
33B
1'B
0!B
/B
.B
,B
%�B
 �B
"�B
�B
(�B
1'B
5?B
9XB
<jB
?}B
A�B
F�B
M�B
`BB
{�B
��B
��BB�B$�B1'B>wBL�BXBYB[#B\)B_;BcTBgmBiyBiyBjBjBjBk�Bm�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Br�Br�Br�Br�Bs�Bs�Br�Bs�Bt�Bt�Bt�Bt�Bu�Bu�Bv�Bv�Bu�Bu�Bu�Bu�Bv�Bu�Bu�Bu�Bu�Bu�Bt�Bt�Bt�Bt�Bt�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bs�Bt�Bs�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bv�Bv�Bw�Bw�Bx�Bx�Bx�Bx�By�By�By�By�Bz�Bz�Bz�B{�B{�B{�B|�B|�B|�B}�B}�B}�B}�B~�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�+B�+B�+B�1B�1B�7B�7B�7B�=B�=B�DB�DB�DB�DB�JB�PB�VB�VB�VB�VB�VB�VB�VB�\B�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@	x�@	��@	�^@	��@	��@	�^@	��@	��@	�^@	��@	�#@	�@ �@1'@��@��@��@��@��@�w@��@	x�@
n�@
�@@
��@	�@�u@�;@��@�;@  @�;@��@  @ �@b@b@�w@
=@E�@?}@�j@��@(�@t�@J?���?�V?�P?���?� �?�I�?��#?���?�hs?�5??+?
=?�P?ȴ?b?�u?�#?5??%�T?2�?CS�?L�D?z�?�1'?�n�?��y?��9?��9?�+?�M�?��w?��?��?�ȴ?�G�?�C�?��u?�+?��/?��!?���?z�H?m��?d�/?a%?[dZ?U�?P��?E�T?=p�?6?2-?,1?'l�?$Z?"M�?   ?dZ?�P?9X?bN?I�?l�?�?%>��>�K�>�!>��>�V>�r�>�l�>��
>�b>�
=>Ձ>���>���>�bN>�+>ě�>Ƨ�>��>���>��>��#>� �>��D>�~�>��>��>��>�>���>��;>���>��^>��9>��9>��9>���>�$�>��>�  >|�>s�F>k�>aG�>_;d>Z�>W
=>R�>R�>L��>H�9>D��>@�>@�><j>9X>333>1&�>,1>&�y>#�
>�w>�u>t�>�>�u>V>V>hs>
=q>1'=��=��=�=�/=���=ȴ9=ě�=��=�j=�Q�=� �=���=��
=��w=��-=��P=��=�hs=�hs=�\)=�C�=�C�=�C�=�C�=�7L=�O�=�\)=��P=���=��w=��-=��-=���=���=��=�t�=�\)=�O�=�O�=�\)=�\)=�\)=�\)=�C�=�7L=�7L=�7L=�+=��=��=�o=�o=�%=q��=u=q��=y�#=q��=aG�=]/=e`B=Y�=T��=L��=L��=L��=P�`=]/=aG�=e`B=e`B=aG�=P�`=H�9=L��=L��=L��=H�9=L��=D��=L��=L��=L��=L��=T��=P�`=P�`=D��=8Q�=��<��<ě�<ě�<�1<e`B;�`B;D��;o:�o��o���
��o�D���D���D���D����o�t��#�
�D���e`B��C����㼓t���o�D���e`B�u��t����㼴9X������`B��h�o�t����#�
�,1�0 Ž8Q�D���T���m�h��%�����7L��C���hs���w���
���1�� Ž� Ž� Ž�Q�\�������
=��
=��"ѽ�G���l���h�����F�����$ݾ1'�I��V������R��w��w�%�T�(�þ,1�/��5?}�@��D���H�9�I�^�Kƨ�Q녾\(��`A��o���vȴ�}󶾁%�������˾�1'���^�����녾�zᾕ���
=������5?��G����徤�/��ff���y���y��ff��ff��ff���T���T���T��ff��`B���T��ff���y��l����y���T���y��l���~������1�����1��1������D���D���D���D��1��1��1���D��1��1��1���D���D��1���D���D���D���D���D1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�z@ǮA�
A#�
AC�
Ac�
A��A��A��A��A��A��A��A��B ��B��B��B�]B ��B(��B0��B8��B@��BH��BP��BX��B`��Bh��Bp��Bx��B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�C =qC=qC=qC=qC=qC
=qC=qC=qC=qC=qC=qC=qC=qC=qC=qC=qC =qC"=qC$=qC&=qC(=qC*=qC,=qC.=qC0=qC2=qC4=qC6=qC8=qC:=qC<=qC>=qC@=qCB=qCD=qCF=qCH=qCJ=qCL=qCN=qCP=qCR=qCT=qCV=qCX=qCZ=qC\=qC^=qC`=qCb=qCd=qCf=qCh=qCj=qCl=qCn=qCp=qCr=qCt=qCv=qCx=qCz=qC|=qC~=qC��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D \D �\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D	\D	�\D
\D
�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D\D�\D \D �\D!\D!�\D"\D"�\D#\D#�\D$\D$�\D%\D%�\D&\D&�\D'\D'�\D(\D(�\D)\D)�\D*\D*�\D+\D+�\D,\D,�\D-\D-�\D.\D.�\D/\D/�\D0\D0�\D1\D1�\D2\D2�\D3\D3�\D4\D4�\D5\D5�\D6\D6�\D7\D7�\D8\D8�\D9\D9�\D:\D:�\D;\D;�\D<\D<�\D=\D=�\D>\D>�\D?\D?�\D@\D@�\DA\DA�\DB\DB�\DC\DC�\DD\DD�\DE\DE�\DF\DF�\DG\DG�\DH\DH�\DI\DI�\DJ\DJ�\DK\DK�\DL\DL�\DM\DM�\DN\DN�\DO\DO�\DP\DP�\DQ\DQ�\DR\DR�\DS\DS��DT�DT�\DU\DU��DU�\B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
%�B
$�B
%�B
%sB
%�B
%�B
&OB
%�B
%�B
%�B
%�B
&�B
( B
-�B
.qB
/�B
1#B
2wB
1�B
1B
/�B
.JB
.�B
1B
2IB
3BB
3B
3$B
4IB
5EB
5�B
5�B
6�B
7B
8�B
9XB
:�B
:�B
9vB
9�B
8%B
3UB
2;B
/�B
/}B
-9B
.�B
,�B
*�B
0B
-/B
1%B
5fB
9 B
<SB
?=B
@�B
E3B
K~B
]VB
y�B
��B
�CB
��BfB$LB1-B?BN�BX�BY�B\B]�BaRBe�BhpBjBjXBkMBj�Bl'Bm�BoEBqbBq�Bq�Bq�Br�BsDBsBssBs�Bs�BsGBtBt(Bs�BtqBu_BuxBu�Bu�BvfBv[Bw/BwgBv2Bu�BvBv3Bv�Bv"Bv�Bu�Bu�Bu�Bt�BuBu�Bt�Bt�BvKBt�Bt�BuBu�BuBt�Bt�BuBu�Bt`Bt�BtBt2Bt�Bt�Bt�Bt�BuBt�Bt�BuBt�Bu+Bu"Bu8Bu�Bv Bu�Bu�Bu�BvBu�Bu�Bu�Bu�Bu�Bu�BvBu�BvBvBu�Bu�BvBv Bu�Bv�BwCBv�Bv�BwBw�Bw9BwBw�Bx)By8Bx�Bx�Bx�By�By�BzBzB{Bz�Bz�B|	B{�B|B|�B|�B}	B}�B}�B}�B~B~�B�B�B��B��B�B�B�B�B�%B�B�*B�B�B�B�B�B�B�5B�,B� B� B�.B�/B�&B�1B�&B�3B�UB�B�3B�B�DB�\B�8B�B�RB�@B�NB�8B�:B�.B�B�5B�8B�DB�SB�{B�iB�IB�VB�WB�bB�LB�mB�@B�\B�`B�aB�MB�tB�jB��B��B��B��B��B�mB��B��B��B��B�vB�wB��B�vB�_B�aB�oB�qB�sB��B��B��B��B��B��B��B�oB�eB�PB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B��B��B��B��B��B�B��B�OB��B��B��B��B��B��B��B�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@	x�@	��@	�^@	��@	��@	�^@	��@	��@	�^@	��@	�#@	�@ �@1'@��@��@��@��@��@�w@��@	x�@
n�@
�@@
��@	�@�u@�;@��@�;@  @�;@��@  @ �@b@b@�w@
=@E�@?}@�j@��@(�@t�@J?���?�V?�P?���?� �?�I�?��#?���?�hs?�5??+?
=?�P?ȴ?b?�u?�#?5??%�T?2�?CS�?L�D?z�?�1'?�n�?��y?��9?��9?�+?�M�?��w?��?��?�ȴ?�G�?�C�?��u?�+?��/?��!?���?z�H?m��?d�/?a%?[dZ?U�?P��?E�T?=p�?6?2-?,1?'l�?$Z?"M�?   ?dZ?�P?9X?bN?I�?l�?�?%>��>�K�>�!>��>�V>�r�>�l�>��
>�b>�
=>Ձ>���>���>�bN>�+>ě�>Ƨ�>��>���>��>��#>� �>��D>�~�>��>��>��>�>���>��;>���>��^>��9>��9>��9>���>�$�>��>�  >|�>s�F>k�>aG�>_;d>Z�>W
=>R�>R�>L��>H�9>D��>@�>@�><j>9X>333>1&�>,1>&�y>#�
>�w>�u>t�>�>�u>V>V>hs>
=q>1'=��=��=�=�/=���=ȴ9=ě�=��=�j=�Q�=� �=���=��
=��w=��-=��P=��=�hs=�hs=�\)=�C�=�C�=�C�=�C�=�7L=�O�=�\)=��P=���=��w=��-=��-=���=���=��=�t�=�\)=�O�=�O�=�\)=�\)=�\)=�\)=�C�=�7L=�7L=�7L=�+=��=��=�o=�o=�%=q��=u=q��=y�#=q��=aG�=]/=e`B=Y�=T��=L��=L��=L��=P�`=]/=aG�=e`B=e`B=aG�=P�`=H�9=L��=L��=L��=H�9=L��=D��=L��=L��=L��=L��=T��=P�`=P�`=D��=8Q�=��<��<ě�<ě�<�1<e`B;�`B;D��;o:�o��o���
��o�D���D���D���D����o�t��#�
�D���e`B��C����㼓t���o�D���e`B�u��t����㼴9X������`B��h�o�t����#�
�,1�0 Ž8Q�D���T���m�h��%�����7L��C���hs���w���
���1�� Ž� Ž� Ž�Q�\�������
=��
=��"ѽ�G���l���h�����F�����$ݾ1'�I��V������R��w��w�%�T�(�þ,1�/��5?}�@��D���H�9�I�^�Kƨ�Q녾\(��`A��o���vȴ�}󶾁%�������˾�1'���^�����녾�zᾕ���
=������5?��G����徤�/��ff���y���y��ff��ff��ff���T���T���T��ff��`B���T��ff���y��l����y���T���y��l���~������1�����1��1������D���D���D���D��1��1��1���D��1��1��1���D���D��1���D���D���D���D���D1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�R<#ו<#ݹ<#�<#��<#��<#�<#�6<#�Q<#�<$;<$`�<#�/<$n<#�<#�<#�	<#�u<#�*<$<�<#�<$4<#��<#�J<#��<$0�<$�<$,�<#��<#�$<#�W<#�<#ݦ<#�<#��<#ݶ<#�N<#��<$/G<$9o<$w<$

<#�<$<$.�<%*<-�]<6��<'��<'��<$S]<%��<%.S<Z$�<|�s<Nf<���<2<#��<#�<#ؠ<#�Y<#�<$3f<%6"<'��<)��<&z�<We�<��<0v�</f\<#��<#�L<$JV<&�<$��<$a<$�[<&��<'�"<(/�<$�<$5<$�%<$��<$�<&[9<(�}<&|�<$tz<$��<%�<$�Z<'p0<&5�<%��<${<%$�<$�Z<$B�<${<$w<$�]<$p�<$P�<$qs<$��<$��<$P�<$Di<$�<$KO<$7<#�4<$�<$#<#�<$�<$��<#�<#��<#��<#��<$	,<$��<#��<#�<$24<#��<#�<$<$�<$
-<#��<#�<$4<$�<$Y�<#��<$<$&q<#�W<#��<#ۂ<#�<#��<#��<#��<$�<#�p<$�<$y<$%�<#�<#��<#�<#��<#�
<#�X<#�5<#�<#��<#��<#��<#�B<#�c<#�<#�j<#�s<#�<#�[<$9<#�<#�<#�<$$><#��<#�,<$�<#��<$R<#��<#�<$P<$P<#��<#�"<#�(<#�<#�x<#�<#�<#�Q<#��<#��<#�\<#�$<#�)<#�L<#��<#�<#�<#�a<#�w<#�<#פ<#ؒ<#�[<#�K<#ן<#�U<#�w<#�
<#��<#��<#�I<#�d<#޵<#�<#��<#�Z<#�[<#�C<#�<#ޢ<#�c<#�r<#�|<#��<#�p<#�|<#ۍ<#��<#�*<#�f<#�U<#��<#��<#��<#ޕ<#�<#��<#�J<#�.<#ۀ<#��<#؍<#�
<#؋<#��<#�V<#�d<#��<#��<#ؿ<#�><#�g<#�P<#�A<#��<#׹<#�9<#ڷ<#�<#��<#އ<#ۮ<#�<#�<$�<$�<$�<#�^<#�<$	<$<#�}<#�<#�.<#�5<#��<#م<#ث<#�I<#��<#��<#�<#�6<#��<#��<#�<#��<#�m<#��<#װ<#ׁ<#�a<#��<#�v<#�<#�O<#��<#�:<#�<#�	<#��<#�<#��<#��<#��<#�<#�F<#�<#��<#�><#�<#��<#�<#�<$H<#�~<#��<#��<#�<#�X<#۬<#��<#��<#��<#�<#�<#�
<#�a<#�<#��<#��<#��<#ޘ<#�<#�f<$a<#�<#�i<#�<$D<$�<#�<#ޏ<#��<#�Q<#�g<#�/<#��<#�h<$0\<#�e<#��<#�"<#�s<#��<$% <#�n<$g�<$
}<$	<#�,<#�&<#��<#�<#��<$#�<#��<#�P<#�<#�s<#�?<$<#��<#�<#� <#��<#ޟ<#�"<#��<#� <#�!<#��<#�0<#�P<#�<#ף<#ޏ<#ގ<#ދ<#�G<#��<#ל<#�"<#�]<#��<#�<#�_<#�<#�+<#�'<#�<#�<#�f<#�2<#�<#��<#�<#�I<#�4<#��<#�<#�I<#�H<#�<#�<#�<#�@<#�.<#�-<#�-<#�-;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.24                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929163046201909291630462019092916305020190929163058202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�DU� DU� DU� G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          