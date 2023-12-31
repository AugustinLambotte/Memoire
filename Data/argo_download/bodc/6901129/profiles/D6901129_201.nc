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
resolution        ?PbM���     �  L�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  S0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X4   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  lD   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  m�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  qH   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  w�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ~�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �X   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �(   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �D   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �`   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �|   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �<   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �,   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �H   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �d   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225044544  20210225044544  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125487                          2C  D   APEX                            6229                            120210                          846 @�9B��Z@1   @�9B��Z@@P�t�j~��57���+1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @9��@�33@�  A��A!��AA��A`  A�  A�33A�  A���A���A���A�  A�  B   B��B  BffB   B'��B0  B8  B@ffBH  BO��BX  B`  Bh  Bp  Bx  B�  B�  B�  B�33B�33B�  B���B���B�  B�  B�33B�  B���B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C�C  C
  C  C  C�C  C  C  C�fC�C  C  C �C"  C#�fC&  C(�C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>�C@  CB  CD  CF  CH  CI�fCL  CN  CP�CR  CT  CV�CX�CZ  C\  C^  C`  Ca�fCd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv�Cx  Cz  C|  C~  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C��C�  C�  C��C�  C��3C�  C�  C��3C�  C�  C��3C�  C�  C�  C�  C��3C�  C��C�  C��3C��3C��3C��3C�  C��C�  C�  C�  C�  C�  C��3C��3C�  C��C�  C�  C�  C��C��C�  C�  C��3C��3C��3C��3C�  C��C��C��C�  C�  C��3C��3C��3C��3C��3C��3C��3C��3C��3C��3C��3C��3C�  C��C��C��C�  C�  C�  C�  C�  C��3C�  C��3C��3C�  C��C��C�  C�  C��C��C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C��C��C�  C��3C��3C��3C��3C�  C�  C�  C��C��C��C�  C�  C��C��C��D fD � D ��Dy�D��D� D  D� DfD�fDfD�fDfD� D��Dy�D��Dy�D��D	y�D	��D
y�D
��Dy�D��Dy�D  D� DfD�fDfD�fDfD� D��D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� DfD�fDfD�fD fD �fD!  D!y�D!��D"y�D#  D#� D$fD$� D$��D%y�D%��D&y�D&��D'� D(  D(�fD)  D)y�D)��D*� D+  D+�fD,fD,� D-  D-�fD.  D.y�D/  D/�fD0  D0y�D1  D1� D2  D2� D3  D3�fD4  D4y�D5  D5� D5��D6� D7  D7y�D8fD8�fD9fD9� D:  D:� D:��D;� D<fD<� D=  D=� D>  D>� D?fD?� D@  D@� DA  DA�fDBfDB� DC  DCy�DDfDD� DE  DE� DF  DFy�DG  DG� DG��DH� DI  DIy�DJ  DJ� DKfDK� DL  DL� DMfDM� DN  DN�fDOfDO� DP  DP� DQ  DQ�fDRfDR� DS  DS� DT  DTy�DU  DU��DU�3B
+B
+B
+B
+B
+B
+B
)�B
(�B
)�B
(�B
(�B
(�B
'�B
&�B
'�B
&�B
&�B
&�B
$�B
"�B
#�B
&�B
&�B
(�B
/B
49B
5?B
5?B
6FB
7LB
8RB
9XB
9XB
9XB
:^B
=qB
B�B
C�B
E�B
I�B
J�B
K�B
M�B
N�B
P�B
P�B
Q�B
Q�B
R�B
R�B
R�B
R�B
R�B
S�B
S�B
S�B
S�B
S�B
T�B
T�B
T�B
T�B
T�B
T�B
VB
W
B
YB
[#B
aHB
y�B
�B
�{B
��B
�B
��B
=B!�B2-B>wBI�BN�BS�BXBZB[#BZBYBZB_;BaHBcTBe`BcTBe`BgmBhsBgmBhsBhsBiyBjBjBjBhsBl�Bl�Bm�Bm�Bp�Bo�Bo�Bn�Bn�Bo�Bo�Bo�Bo�Bp�Bo�Bp�Bq�Br�Br�Bq�Bq�Bq�Bp�Bp�Bq�Bp�Bo�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bo�Bp�Bp�Bp�Bp�Bo�Bp�Bp�Bq�Bq�Bq�Br�Br�Br�Br�Br�Bs�Bs�Bt�Bu�Bu�Bv�Bw�Bw�Bw�Bw�Bv�Bw�Bw�Bw�Bw�Bx�Bx�By�By�By�By�By�By�By�By�By�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{�B|�B{�B{�B|�B|�B}�B}�B}�B}�B}�B~�B~�B~�B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�+B�+B�+B�1B�1B�1B�7B�7B�7B�=B�=B�DB�DB�JB�JB�JB�PB�PB�PB�PB�PB�\B�\B�bB�bB�bB�hB�hB�hB�oB�oB�uB�oB�uB�uB�uB�{B�{B�{B�{B�{B�{B�uB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?��?���?��P?�K�?��P?�l�?��?��P?��?�z�?���?��j?��?���?�&�?��?� �?�R?�\)?�|�?�C�?���?�V?�O�?ݑh?�p�?���?���?�I�?��m?ۅ?�?���?��H?���?ٙ�?�+?�$�?�`B?��
?�-?�G�?��`?�%?�hs?щ7?��?��?�J?�J?�-?�-?�-?�M�?�M�?�hs?��`?� �?�Ĝ?��`?ϝ�?�V?Ͳ-?Ͳ-?͑h?�I�?�(�?�V?У�?�E�?׍P?�^5?�1'?��?Ƨ�?��w?���?��T?ȓu?�x�?��H?�?�/?��#?��?���?�ff?�t�?��?�o?��?��+?�E�?��?���?�?�S�?�&�?|j?wK�?st�?mO�?h1'?bJ?NV?C�
?>�R?;"�?(�9? Ĝ?��?��?�?�?��?`B?��>��m>�&�>�S�>ؓu>�z�>�n�>��`>���>�7L>�1'>Ƨ�>��>�?}>��>��>�bN>��`>��>��`>���>��^>�+>�$�>��>�  >vȴ>hr�>cS�>^5?>["�>I�^>8Q�>0 �>+>(��>#�
>�>z�>\)>%=�F=�=���>V>bN>hs>\)>\)>O�>O�>V>%=�=�=�h=�S�=��`=���=���=���=��=���=���=���=���=���=���=���=���=��=�Q�=�9X=��
=���=�\)=�O�=��=�o=aG�=49X=0 �=0 �='�=�w=C�=o<��<��<��<�`B<�h<�h<�<�<�h<��<��=C�=D��=L��=L��=H�9=H�9=D��=@�=<j=8Q�=8Q�=<j=8Q�=8Q�=0 �='�=t�=+<�`B<�j<�9X<�9X<�9X<�9X<�9X<�1<�1<�9X<�9X<�1<�j<ě�<ě�<ě�<�9X<���<�t�<�o<u<u<e`B<T��<D��<D��<49X<49X<#�
;ě�;�o;ě�;�`B<o;��
;D��;�o;ě�;�`B;�`B;�`B;�`B<o;�`B;��
;��
;�o:�o            ��o�D����`B�t��49X�D���T���e`B��C���C���C���t���1��1��9X��9X��9X��j�C��D���]/�q���u�}󶽙�����w������ Ž�9X��E���Q����\�ě����������"ѽ�h���ٽ��#�o�%���$ݾo���$ݾ1'�
=q�C��O߾����!���#�
�&�y�+�,1�1&�5?}�7KǾ=p��A�7�F��J���W
=�\(��cS��k��o����  ���𾇮�����������1'���9���9��7L��7L���^��C���C���ƨ��ƨ��ƨ��ƨ��ƨ��I����;��;��;��;��;��;�O߾�O߾�O߾�O߾�������V�����������������\)��\)��\)��\)���;��hs��񪾕�������㾜(����-��;d���w��A���A���Ĝ��G���G�����������������������MӾ�����MӾ�MӾ��徢�徢�徤Z��Z��Z���T����111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@6fg@���@�ffA ��A ��A@��A_33A33A���A���A�fgA�fgA�fgAߙ�AA���BfgB��B33B��B'fgB/��B7��B@33BG��BOfgBW��B_��Bg��Bo��Bw��B��B��fB��fB��B��B��fB��3B��3B��fB��fB��B��fB��3B��fB��fB��fB��3B��fB��fB��fB��fBӳ3B��fB��fB��fB��fB��fB��fB��fB��fB��fB��fB��fC�3C�3C�C�3C	�3C�3C�3C�C�3C�3C�3CٙC�C�3C�3C �C!�3C#ٙC%�3C(�C)�3C+�3C-�3C/�3C1�3C3�3C5�3C7�3C9�3C;�3C>�C?�3CA�3CC�3CE�3CG�3CIٙCK�3CM�3CP�CQ�3CS�3CV�CX�CY�3C[�3C]�3C_�3CaٙCc�3Ce�3Cg�3Ci�3Ck�3Cm�3Co�3Cq�3Cs�3Cv�Cw�3Cy�3C{�3C}�3C�3C�gC���C���C���C���C���C���C���C���C�gC���C���C�gC���C���C�gC���C���C���C���C���C���C���C���C���C���C���C���C���C���C�gC���C���C���C���C���C���C�gC���C���C���C���C���C���C���C���C�gC���C���C���C�gC�gC���C���C���C���C���C���C���C�gC�gC�gC���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C�gC�gC�gC���C���C���C���C���C���C���C���C���C���C�gC�gC���C���C�gC�gC���C���C���C���C���C���C���C���C���C���C���C���C�gC�gC�gC���C���C���C���C���C���C���C���C�gC�gC�gC���C���C�gC�gC�gD 3D |�D �gDvgD�gD|�D��D|�D3D�3D3D�3D3D|�D�gDvgD�gDvgD�gD	vgD	�gD
vgD
�gDvgD�gDvgD��D|�D3D�3D3D�3D3D|�D�gD|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D��D|�D3D�3D3D�3D 3D �3D ��D!vgD!�gD"vgD"��D#|�D$3D$|�D$�gD%vgD%�gD&vgD&�gD'|�D'��D(�3D(��D)vgD)�gD*|�D*��D+�3D,3D,|�D,��D-�3D-��D.vgD.��D/�3D/��D0vgD0��D1|�D1��D2|�D2��D3�3D3��D4vgD4��D5|�D5�gD6|�D6��D7vgD83D8�3D93D9|�D9��D:|�D:�gD;|�D<3D<|�D<��D=|�D=��D>|�D?3D?|�D?��D@|�D@��DA�3DB3DB|�DB��DCvgDD3DD|�DD��DE|�DE��DFvgDF��DG|�DG�gDH|�DH��DIvgDI��DJ|�DK3DK|�DK��DL|�DM3DM|�DM��DN�3DO3DO|�DO��DP|�DP��DQ�3DR3DR|�DR��DS|�DS��DTvgDT��DU��DU� B
*�B
+ B
+B
*�B
+B
*�B
*B
)�B
*YB
(�B
(�B
)1B
(B
'�B
(1B
'B
'eB
&�B
$�B
$�B
(B
'�B
&�B
(�B
/,B
4jB
5EB
5uB
6pB
7vB
8�B
9sB
9RB
9mB
:�B
>eB
B�B
C�B
F5B
JRB
KB
K�B
M�B
N�B
P�B
P�B
Q�B
Q�B
R�B
R�B
R�B
R�B
R�B
S�B
TOB
T.B
T@B
S�B
T�B
U|B
UzB
U?B
UB
UB
V}B
WB
X�B
Y�B
_AB
y:B
�@B
��B
�B
�B
��B	�B 3B12B>BI;BP�BWBYCBZ�B[5BZ�BZLB]�Ba�Bc�Be�Be�Bd�Bh�Bh�BiwBhIBi�BinBjBBk�Bk|Bk�BlMBn�Bm�BnZBq Br1Bq�Bp�Bo?BoABp:Bo�Bo�Bp�Bq�Bp�Bq�BrBr�Br�Bq�BrBq�Bp�BqoBr�Bq�BqXBqBp�Bp�Bp�BqBp�Bp�Bp�Bp�BqBpBqDBp�Bp�Bp�BptBq}Bq
Bq�Bq�Bq�BsBsBr�BsbBsBs�Bs�Bs�Bu�Bu�Bv�Bw�Bw�Bw�Bw�BwjBxBw�BxBxBy?Bx�By�By�By�By�By�By�By�By�By�By�By�BzB{Bz�B{CB{$B{Bz�B{Bz�B{SB|mB|�B{�B|B}	B}-B~B~B}�B}�B~B~�B~�B~�B�B�B�B�B��B�_B��B�B�B�B�B�B�B� B�B�B�"B�B�0B�1B�WB�@B�ZB�ZB�-B�!B�$B�%B�&B�2B�%B�B�+B�9B�B�%B�4B�5B�QB�]B�IB�YB�NB�EB�UB�WB�ZB�OB�\B�QB�]B��B�rB�FB�SB�WB��B�B�\B�RB�`B�qB�uB�sB�gB��B��B�{B��B��B��B�{B�yB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�1B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB��B��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�6B��B��B�B��B�ZB�@B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��?��?���?��P?�K�?��P?�l�?��?��P?��?�z�?���?��j?��?���?�&�?��?� �?�R?�\)?�|�?�C�?���?�V?�O�?ݑh?�p�?���?���?�I�?��m?ۅ?�?���?��H?���?ٙ�?�+?�$�?�`B?��
?�-?�G�?��`?�%?�hs?щ7?��?��?�J?�J?�-?�-?�-?�M�?�M�?�hs?��`?� �?�Ĝ?��`?ϝ�?�V?Ͳ-?Ͳ-?͑h?�I�?�(�?�V?У�?�E�?׍P?�^5?�1'?��?Ƨ�?��w?���?��T?ȓu?�x�?��H?�?�/?��#?��?���?�ff?�t�?��?�o?��?��+?�E�?��?���?�?�S�?�&�?|j?wK�?st�?mO�?h1'?bJ?NV?C�
?>�R?;"�?(�9? Ĝ?��?��?�?�?��?`B?��>��m>�&�>�S�>ؓu>�z�>�n�>��`>���>�7L>�1'>Ƨ�>��>�?}>��>��>�bN>��`>��>��`>���>��^>�+>�$�>��>�  >vȴ>hr�>cS�>^5?>["�>I�^>8Q�>0 �>+>(��>#�
>�>z�>\)>%=�F=�=���>V>bN>hs>\)>\)>O�>O�>V>%=�=�=�h=�S�=��`=���=���=���=��=���=���=���=���=���=���=���=���=��=�Q�=�9X=��
=���=�\)=�O�=��=�o=aG�=49X=0 �=0 �='�=�w=C�=o<��<��<��<�`B<�h<�h<�<�<�h<��<��=C�=D��=L��=L��=H�9=H�9=D��=@�=<j=8Q�=8Q�=<j=8Q�=8Q�=0 �='�=t�=+<�`B<�j<�9X<�9X<�9X<�9X<�9X<�1<�1<�9X<�9X<�1<�j<ě�<ě�<ě�<�9X<���<�t�<�o<u<u<e`B<T��<D��<D��<49X<49X<#�
;ě�;�o;ě�;�`B<o;��
;D��;�o;ě�;�`B;�`B;�`B;�`B<o;�`B;��
;��
;�o:�o            ��o�D����`B�t��49X�D���T���e`B��C���C���C���t���1��1��9X��9X��9X��j�C��D���]/�q���u�}󶽙�����w������ Ž�9X��E���Q����\�ě����������"ѽ�h���ٽ��#�o�%���$ݾo���$ݾ1'�
=q�C��O߾����!���#�
�&�y�+�,1�1&�5?}�7KǾ=p��A�7�F��J���W
=�\(��cS��k��o����  ���𾇮�����������1'���9���9��7L��7L���^��C���C���ƨ��ƨ��ƨ��ƨ��ƨ��I����;��;��;��;��;��;�O߾�O߾�O߾�O߾�������V�����������������\)��\)��\)��\)���;��hs��񪾕�������㾜(����-��;d���w��A���A���Ĝ��G���G�����������������������MӾ�����MӾ�MӾ��徢�徢�徤Z��Z��Z���T����111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�v<#��<#ר<#٤<#�<#�+<#�`<$F�<#�~<#��<#��<#ݜ<#��<$��<#�t<#�2<#��<#�<#�<&"�<0͸<$F�<#؅<#�+<#�<#��<#�,<#�8<#١<#٫<#��<#פ<#� <#�B<#�|<$zz<#�Q<#�#<$*<$S<#��<#��<#؏<#۫<#�)<#�J<#�U<#��<#�m<#�'<#�8<#�=<#��<#�6<#�a<#܉<#�<#�<#�<#��<#�x<#߼<#�<#�O<#�><#�*<#��<%:�<'(O<$1�<$u�<%�<+wN<*��<(�<$q<%�k<$�V<#�U<$:<&�C<*�~<$�|<$:�<#�R<$�<$��<.�;<'΍<(0�<'��<#�<%+5<+��<%_�<$�m<$_ <$��<$��<$H�<$�<$��<$�;<.��<&�8<$��<$H�<-I�<%�B<'}�<$�+<$$<$&<$�<#��<#�<$��<$�E<% <<$�b<#�L<#۬<#�5<#��<#��<#�<#�(<$J�<$c�<%W<&8<#��<#��<#��<#�<#�<#�<#��<#�5<#�;<#�c<#��<$|<#�L<#�<#�u<$X2<$[�<#�<#��<#�a<#�<#�H<#�<#��<$.�<#�O<#�><#ٟ<${:<#�j<#��<#��<#�3<#ױ<#�5<#��<$<#��<#�<#��<#�w<#�I<#�<#�v<#�D<#�7<#�<#�6<#�;<#�;<#�;<#�:<#�5<#�'<#�\<#�k<#�\<#�<#�f<#��<#�B<#��<#�Q<#�c<$m<#�<#�<#� <#�<#�@<#��<#�<#�8<#�.<#ٿ<#�^<#�A<#�<#�`<#�<#�	<#�R<#�Z<$4�<#��<#�=<#�<#�0<#�<#�<#�<#�<#�<#��<#�4<#�J<#؁<#��<#߿<#��<#��<#��<#�%<#�+<#�I<#�:<#�5<#�<#�@<#�c<#�:<#�&<#�J<#�G<#�"<#�X<#��<#��<#�<#�1<#�<#�4<#�<#�<#�8<#�K<#�<#�2<#�<#�j<#ם<#ٺ<#غ<#�<#��<#׻<#�O<#��<#ع<#�'<#�:<#�<#؂<#�<#��<#�><#� <#�y<#�<#�9<#�X<#�2<#�-<#��<#�<#��<#�$<#�"<#�<#�<#�+<#� <#�$<#�s<#�<#�<#�*<#�2<#�*<$r<$'*<#��<#��<#�<#��<$]<#��<#�g<#�-<#�<#�<#�(<#ܩ<#�<#�%<#�<#��<#�e<#�*<#�<#�N<#��<#�A<#��<#�<#ܱ<#��<#�*<#��<#��<#�%<#�6<#��<#�d<#�-<#׾<#ڣ<#�=<#�*<#ߠ<#��<#�1<#��<#��<#�H<#ޝ<$r<#�X<#�<#�3<#�K<$=�<$�<#�<#�2<#�9<#�N<#�g<#�<#�1<#�<#�<#�5<#�&<#�3<#�<#�\<#�<#�9<#�'<#�2<#�<#�4<#�9<#�9<#�9<#�4<#�<#�7<#�9<#�3<#�<#�8<#�4<#�<#�<#�9<#�H<#�<#א<#�2<#�9<#�#<#�3<#�(<#�l<#�	<#��<$z<#��<#�,<#�s<#�0<#�<#�<#�?<#�U<#�<#�1<#�<#�5<#�8<#�8<#�3<#�<#�<#�
<#�&<#�'<#�q<#�<#��<#�@<#�*<#ڞ<%9�<#ߥ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.05                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929162634201909291626342019092916263820190929162645202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@9��@9��@9��G�O�G�O�G�O�G�O�DU�3DU�3DU�3G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          