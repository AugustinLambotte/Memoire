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
resolution        ?PbM���     �  dh   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  l@   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  t   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  v   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  x   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  z    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �<   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �X   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �t   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �l   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �\   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �x   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225043448  20210225043448  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125476                          2C  D   APEX                            6229                            120210                          846 @�}� �1   @�}� �@P��G�{�4E����1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�ff@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DUy�DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Dby�Dc  Dc�fDd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� B
��B
�B
�B
�B
�B
�B
�`B
�BB
�B
��B
��B
��B
��B
��B
��B
��B
��BBoB�B�B�B{B$�B!�B!�B&�B'�B(�B%�B#�B"�B%�B!�B#�B$�B%�B'�B+B1'B49B49B7LB:^B<jB=qB@�B?}B>wBA�BC�BG�BK�BM�BN�BM�BM�BM�BR�BT�BVBVBVBW
BZB\)B[#B[#B\)B\)B\)B[#B[#B\)B]/B\)B[#B[#B[#B[#B[#B[#B[#B[#B[#B[#B[#B[#B[#B[#B\)B\)B\)B\)B\)B]/B]/B]/B^5B^5B^5B^5B^5B^5B^5B_;B_;B`BB`BB`BB`BB_;B_;B`BB`BB_;B_;B_;B`BB`BBaHBaHBaHBaHBaHBbNBbNBcTBcTBdZBcTBcTBdZBdZBe`Be`BffBffBffBffBffBgmBhsBhsBiyBiyBiyBiyBjBjBk�Bk�Bk�Bk�Bl�Bl�Bl�Bm�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bp�Bp�Bp�Bq�Bq�Br�Br�Br�Br�Br�Bs�Bt�Bv�Bv�Bv�Bv�Bw�Bw�Bx�By�Bz�B{�B|�B}�B}�B~�B� B� B� B~�B}�B|�B|�B|�B}�B~�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�1B�1B�7B�7B�7B�=B�=B�=B�DB�DB�JB�JB�JB�JB�PB�PB�PB�PB�PB�VB�VB�\B�\B�bB�hB�hB�hB�hB�hB�hB�bB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B��B�{B��B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��Az�A33A �+@�`B@�ff@��@��@�9X@�S�@���@��D@��@��-@���@�bN@�(�@���@|��@k"�@b�@[�@R��@GK�@0bN@,�@%/@K�@��@�?��m?���?���?Ұ!?�dZ?�"�?��w?�
=?�{?�
=?��-?��#?�b?���?��?�`B?�`B?�J?�&�?|(�?qhs?i7L?a�7?_|�?_;d?^�R?[��?MV?CS�?B�\?@Ĝ?>v�?8Q�?5�?5�?6�+?4�j?2�?-V?&ff?!�7?dZ?�?t�?1?`B?`B?`B>�v�>�r�>�b>��`>ɺ^>�%>�^5>�>�33>���>�r�>�G�>���>���>��>�hs>�V>���>���>�J>��>x��>r�!>k�>dZ>\(�>Y�>Xb>Y�>^5?>["�>Xb>T��>N�>@�>:^5>;dZ>?|�>8Q�>%�T>��>�w>��>�u>z�>\)>1'>�>�>	7L>
=q>C�>�>J=��=��>%>�>+=��=��#=��#=��#=��#=�=�=�=�=�=�=���>   >J>o>%>%>%=��#=��=�>   >o>$�>1'>J=�=���=�F=�G�=�/=�"�=�"�=�G�=�S�=�G�=�h=�F==�l�=�l�=�l�>   >C�>t�>n�>o>O�>\)>V>�P>�R>#�
>-V>0 �>49X>;dZ>=p�>=p�>>v�>?|�>333>z�>I�>J>�>I�>��>$�/>5?}>49X>&�y>�->�w>�->�>�>�->�->�R>�w>$�/>&�y>'�>+>-V>/�>0 �>0 �>0 �>1&�>1&�>1&�>0 �>/�>.{>.{>.{>-V>-V>-V>.{>/�>5?}>7K�>49X>7K�>C��>D��>A�7>@�>=p�>=p�>5?}>0 �>'�>#�
>"��>"��>!��> Ĝ> Ĝ>�R>��>��>�+>t�>bN>\)>O�>I�>
=q>
=q>+>$�>�>�>�>J=��=��#=�F==�h=�=�`B=�S�=�G�=��=��=�v�=�9X=�{=��=��
=��w=��P=��=�t�=�hs=�O�=�7L=��=}�=y�#=e`B=]/=T��=L��=D��=<j=0 �=#�
=�w=��=t�<�<�`B<�/<���<ě�<ě�<��
<�C�<�o<u<T��<D��<#�
<t�;�o;D���o�ě��t��49X�e`B�u��o���
�ě������o�o�C��'49X�8Q�@��P�`�]/�ixսy�#�}󶽇+��O߽��P���P�������-���w���
���1�� Ž�-��Q�ě�����������`���`��
=��"ѽ�/��G���S���l���xս����ٽ��#���%�J���+�1'�1'�
=q�I��\)�t�����+��P�������w�&�y�,1�/��2-�5?}�9X�;dZ�?|�A�7�E�˾J���O�;�S�ϾT���V�W
=�Xb�Y��Y��Z��^5?�dZ�gl��hr��l�D�n���p�׾t�j�vȴ�z�H�{�m�}󶾂J�����$ݾ�����9���^����������ƨ���;�\)��hs���Ͼ���b��������������㾞5?��Ĝ��S���`B��l����羬1���h����������׾��!��ȴ��X��^5���H���H���m���۾Ƨ�����bN����t����ϾՁ�׍P�ٙ��ٙ��ٙ��ٙ��ٙ��ٙ�����������������������ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��������������ؓu1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@�G�A ��A ��A@��A`��A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�B (�B(�B(�B(�B (�B((�B0(�B8(�B@(�BH(�BP(�BX(�B`(�Bh(�Bp(�Bx(�B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{C 
=C
=C
=C
=C
=C

=C
=C
=C
=C
=C
=C
=C
=C
=C
=C
=C 
=C"
=C$
=C&
=C(
=C*
=C,
=C.
=C0
=C2
=C4
=C6
=C8
=C:
=C<
=C>
=C@
=CB
=CD
=CF
=CH
=CJ
=CL
=CN
=CP
=CR
=CT
=CV
=CX
=CZ
=C\
=C^
=C`
=Cb
=Cd
=Cf
=Ch
=Cj
=Cl
=Cn
=Cp
=Cr
=Ct
=Cv
=Cx
=Cz
=C|
=C~
=C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU|)DV�DV��DW�DW��DX�DX��DY�DY��DZ�DZ��D[�D[��D\�D\��D]�D]��D^�D^��D_�D_��D`�D`��Da�Da��Db�Db|)Dc�Dc��Dd�Dd��De�De��Df�Df��Dg�Dg��Dh�Dh��Di�Di��Dj�Dj��Dk�Dk��Dl�Dl��Dm�Dm��Dn�Dn��Do�Do��Dp�Dp��Dq�Dq��Dr�Dr��Ds�Ds��Dt�Dt��Du�Du��Dv�Dv��Dw�Dw��Dx�Dx��Dy�Dy��Dz�Dz��D{�D{��B
�	B
�B
�*B
��B
��B
��B
�zB
�|B?B�B�B�BcB�BtB�B	�B�B�B^B�BB$�B'�B'{B+�B/�B.XB,�B*�B(B) B(�B'�B("B(/B)AB*�B.yB2�B4�B5$B7|B:GB<rB>�B@�B@�B@�BCBEBHBK�BM�BOzBP�BO�BNBSABUsBW*BV�BVBV�BZrB\�B\AB\gB]B]OB\�B\7B\�B]dB]2B\5B\QB]4B\�B[�B[�B[�B[�B[�B[fB[�B[�B[�B[�B[^B\�B\�B\vB\�B\�B]�B]XB]�B^�B^�B^�B^�B^\B^BB^*B_ B_aB`eB`jB`�B`�B_�B_4B`B`�B`B_�B^�B`rB`jBaxBa�Ba�BamBaUBbBbCBcJBc�Bd|Bc|BcWBdABd7Be=Be�BfBfgBfgBfgBf�BgnBhqBhvBiwBiyBilBiLBjcBjvBk�Bk�Bk�Bk�Bl�BlrBlQBmkBnqBn�Bo�Bo�Bo�Bo�BpBo�Bo�Bp�Bp�Bp�Bq�BqfBr�Br�Br�Br�Br�Bs Bt5BvbBv�BwyBvUBw�Bw�BxhBy�Bz�B{{B|�B}�B}�B~�B�B�B�B�B\B}_B}fB|�B}�B~RB��B�EB�B��B�zB��B�$B�B�B�	B�B�B�B��B�B�B�B�B�B�+B�7B�:B�/B�=B�?B�NB�TB�TB�JB�KB�YB�NB�PB�DB�BB�B�<B�yB�2B��B�SB��B�vB��B�mB��B��B��B��B�vB�iB�uB�tB�mB��B��B�rB��B��B��B�}B��B��B��B�wB��B��B��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�jB�tB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��Az�A33A �+@�`B@�ff@��@��@�9X@�S�@���@��D@��@��-@���@�bN@�(�@���@|��@k"�@b�@[�@R��@GK�@0bN@,�@%/@K�@��@�?��m?���?���?Ұ!?�dZ?�"�?��w?�
=?�{?�
=?��-?��#?�b?���?��?�`B?�`B?�J?�&�?|(�?qhs?i7L?a�7?_|�?_;d?^�R?[��?MV?CS�?B�\?@Ĝ?>v�?8Q�?5�?5�?6�+?4�j?2�?-V?&ff?!�7?dZ?�?t�?1?`B?`B?`B>�v�>�r�>�b>��`>ɺ^>�%>�^5>�>�33>���>�r�>�G�>���>���>��>�hs>�V>���>���>�J>��>x��>r�!>k�>dZ>\(�>Y�>Xb>Y�>^5?>["�>Xb>T��>N�>@�>:^5>;dZ>?|�>8Q�>%�T>��>�w>��>�u>z�>\)>1'>�>�>	7L>
=q>C�>�>J=��=��>%>�>+=��=��#=��#=��#=��#=�=�=�=�=�=�=���>   >J>o>%>%>%=��#=��=�>   >o>$�>1'>J=�=���=�F=�G�=�/=�"�=�"�=�G�=�S�=�G�=�h=�F==�l�=�l�=�l�>   >C�>t�>n�>o>O�>\)>V>�P>�R>#�
>-V>0 �>49X>;dZ>=p�>=p�>>v�>?|�>333>z�>I�>J>�>I�>��>$�/>5?}>49X>&�y>�->�w>�->�>�>�->�->�R>�w>$�/>&�y>'�>+>-V>/�>0 �>0 �>0 �>1&�>1&�>1&�>0 �>/�>.{>.{>.{>-V>-V>-V>.{>/�>5?}>7K�>49X>7K�>C��>D��>A�7>@�>=p�>=p�>5?}>0 �>'�>#�
>"��>"��>!��> Ĝ> Ĝ>�R>��>��>�+>t�>bN>\)>O�>I�>
=q>
=q>+>$�>�>�>�>J=��=��#=�F==�h=�=�`B=�S�=�G�=��=��=�v�=�9X=�{=��=��
=��w=��P=��=�t�=�hs=�O�=�7L=��=}�=y�#=e`B=]/=T��=L��=D��=<j=0 �=#�
=�w=��=t�<�<�`B<�/<���<ě�<ě�<��
<�C�<�o<u<T��<D��<#�
<t�;�o;D���o�ě��t��49X�e`B�u��o���
�ě������o�o�C��'49X�8Q�@��P�`�]/�ixսy�#�}󶽇+��O߽��P���P�������-���w���
���1�� Ž�-��Q�ě�����������`���`��
=��"ѽ�/��G���S���l���xս����ٽ��#���%�J���+�1'�1'�
=q�I��\)�t�����+��P�������w�&�y�,1�/��2-�5?}�9X�;dZ�?|�A�7�E�˾J���O�;�S�ϾT���V�W
=�Xb�Y��Y��Z��^5?�dZ�gl��hr��l�D�n���p�׾t�j�vȴ�z�H�{�m�}󶾂J�����$ݾ�����9���^����������ƨ���;�\)��hs���Ͼ���b��������������㾞5?��Ĝ��S���`B��l����羬1���h����������׾��!��ȴ��X��^5���H���H���m���۾Ƨ�����bN����t����ϾՁ�׍P�ٙ��ٙ��ٙ��ٙ��ٙ��ٙ�����������������������ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��ٙ��������������ؓu1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<1T�<���<kk<��,<$��<2\/<�<�M�<�Zy<[�><yy�<-�:<�k�<y�<w��<M�j<���<�4P<B�}<5}�<A\�<VL
<�Q}<*��<;�<f��<VLM<Ajd<.�<5	o<1c�<@��<*W�<>`�<1��<,W�<,��<)�	<-9<%��<$9�<$�w<#߶<#��<#ן<%3<#�Y<$�V<'�<%ǯ<%� <#�M<#��<#�_<$+�<)�{<&nl<#��<#��<$�<$��<$<#ׁ<#�(<#�I<#�
<$�w<%"<$�i<$�<$�<$�D<%q�<%L<#�H<#��<$��<'6Q<%�<$BO<$=�<$fI<$-�<$�<#�<$Q<$"Q<$:�<$�<#�<#�A<#�b<#��<#��<$|<#�<#ݗ<#��<#�<#�<#�}<#�o<#�B<#� <#� <#ߞ<#��<#�<<#�c<#�<$5�<#��<#�<#��<#�)<$v�<$�<#�O<#߿<#�|<#�<#��<#��<#ܺ<#�#<#�U<#�<#�<#��<#��<#�h<#�A<#�.<#ٙ<#ٲ<#�?<#��<#�-<#�*<#�3<#�.<#�*<#�<#�M<#�<#�%<#�5<#��<#؊<#�<#��<#�0<#�><#�1<#ߍ<#�)<#�J<#�:<#�Z<#��<#�}<#�<#�<#ړ<#��<#څ<#�|<#�<#��<#�<#�{<#��<#��<#�<#��<#�:<#�<$�<$	�<#��<#�l<$;�<#�J<#ط<#��<#�s<#�<#�<#�"<#ڢ<#ݗ<#�<#�A<#�<#�$<#�
<$& <%m�<$$<$b<#��<#�C<$'Z<$v<$A�<#��<$'�<$ 5<#׿<#��<#�o<#�<#�<#�%<#�<#�i<#�<#�H<#�)<#��<#؍<#��<#�%<#�&<#�F<#�><#�'<#�:<#��<#�z<#��<#�)<#�-<#�]<#�<#�$<#�)<#�C<#��<#�C<#�2<#��<$�<#�H<#�8<#�B<#�g<#�j<#�p<#�%<#��<#��<#�4<#�/<#�<#�<#�`<#ٴ<#�2<#�M<#�<#�<#ܲ<#�5<#�e<#��<#��<#�6<#�t<#�7<#�<#�<<#�<#�<#ܳ<#�/<#��<#�<#�+<#؍<#�<#�s<#׸<#�<#�<#��<#��<#�)<#��<#�<#�5<#��<#��<#�<#�+<#��<#�<#�<#ܪ<#�_<#�|<#ٽ<#�
<#�
<#�<#� <#��<#ܱ<#�(<#�-<#�K<#��<#��<#��<#��<#�<#�0<#�0<#�<#�8<#�*<#��<#�3<#��<#�_<#�<#�k<#��<#��<#�t<#�<#ܝ<#�9<#�F<#��<#�<#��<#�`<#�<#�2<#�Y<#��<#�<#�6<#�-<#�_<#�<#�
<#�Y<#�#<#�5<#�C<#�<#�6<#�<#��<#�3<#�
<#�
<#�	<#��<#�?<#�1<#�,<#��<#٢<#�<#�:<#�x<#�<#�/<#��<#�6<#��<#�j<#�P<#�K<#�*<#��<#��<#�=<#ܚ<#�<#�<#�1<#��<#� <#��<#�k<#�<#�)<#�7<#ܨ<#�#<#�U<#�<#�<#ܷ<#��<#�<#�Z<#�T<#��<#�V<#�<#�v<#�m<#�f<#�@<#�<#�<#�<#�<#�1<#�:<#��<#�<#�<#�`<#�/<#�<#�5<#�6<#�`<#��<#�<#�N<#�'<#�<#�<#��<#�<#ڎ<#٩<#�(<#�i<#�M<#�<#��<#�T<#��<#��<#�W<#�<#�6<#٘<#�<#�<#�w<#��<#��<#�3<#��<#݈<#ܹ<#�=<#�:<#��<#��<#�<#�<#�<#�6<#�@<#��<$K�<$P�<#��<#�<#�2<#�8<#��<#�<#�7<#�'<#�'<#�'<#�'<#�,<#�
<#�+<#�'<#�'<#�'<#�'<#�#<#�$<#�&<#�'<#�'<#�'<#�'<#�'<#�'<#�'<#�#<#�&<#�%<#�'<#�#<#�<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.04                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929160339201909291603392019092916034320190929160350202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�D{� D{� D{� G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          