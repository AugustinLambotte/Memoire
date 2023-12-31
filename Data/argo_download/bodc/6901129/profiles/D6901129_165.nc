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
_FillValue                 �  Xx   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Zl   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  \`   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  d,   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  k�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  u�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  w�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  y�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �l   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �8   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �d   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �d   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �d   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �d   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �(   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �L   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �h   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �,Argo profile    3.1 1.2 19500101000000  20210225042526  20210225042526  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125453                          2C  D   APEX                            6229                            120210                          846 @��n��� 1   @��n��� @Q$�/�4l�����1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            &A   A   A   @�  @���A   A   A@  A`  A�  A�  A�  A�  A�  A���A�  A�  B ffB  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#y�D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� DafDa�fDbfDb�fDc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw�fDx  Dx� Dy  Dy�fDy��B
�dB
�dB
�^B
�dB
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�^B
�XB
�^B
�XB
�XB
�XB
�^B
�^B
�^B
�^B
�dB
�dB
�dB
�dB
�jB
�jB
�qB
�qB
�wB
�}B
�}B
��B
ÖB
��B
��B
�B
�sB
��BB
=BuB�B�B1'B<jBC�BJ�BJ�BO�BO�B[#BffBhsBiyBiyBiyBjBjBjBk�Bk�BjBk�Bl�Bl�Bl�Bm�Bm�Bm�Bn�Bp�Bp�Bp�Br�Br�Br�Br�Bq�Bq�Bq�Bq�Bq�Br�Bs�Bt�Bv�Bw�By�Bz�B{�B|�B~�B� B�B�B�B�B�B�B�B�B�B�%B�+B�+B�%B�B�B� B|�B{�B|�B}�B}�B}�B|�B{�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{�B|�B}�B}�B}�B}�B}�B}�B}�B~�B~�B~�B~�B~�B~�B~�B� B� B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�1B�7B�1B�1B�1B�1B�1B�1B�1B�1B�1B�7B�7B�=B�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�PB�PB�PB�PB�PB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�bB�bB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�uB�uB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@Q�@Q�@Q�@Q�@Q�@bN@Q�@A�@A�@Q�@Q�@Q�@Q�@Q�@Q�@Q�@A�@Q�@�y@v�@@p�@��@��@33@�\@��@%@bN@\)@�T@�@j@1@"�@
M�@	�#@	&�@�@�j@��@ r�?�x�?�J?��?�?�G�?�I�?��T?˅?ě�?���?�Ĝ?�33?�;d?���?U?}?,��?!��?�?��?�j?hs?\)?	x�?�?��>��>�9X>�V>�1>��T>��
>���>���>��
>�~�>�`B>�;d>�A�>�Ĝ>�
=>��`>�O�>�I�>�C�>��m>��#>�dZ>�j>��>��>�I�>�t�>��>�A�>�Z>�>��>�9X>�F>��>�`B>�M�>��T>�x�>�h>�V>�V>�->��>�h>ݲ->�>��7>�r�>��>�\)>�b>���>�>��`>���>}�>y�#>j~�>\(�>Y�>Z�>["�>["�>["�>Xb>Xb>V>T��>Q�>O�;>M��>Kƨ>I�^>I�^>I�^>H�9>H�9>G�>F��>E��>E��>D��>C��>A�7>A�7>D��>F��>E��>D��>;dZ>8Q�>8Q�>?|�>F��>F��>F��>F��>J��>J��>J��>J��>I�^>D��>B�\>;dZ>5?}>2->1&�>,1>(��>(��>&�y>&�y>%�T>%�T>$�/>$�/>#�
>�R>+>(��>!��>�w>�P>t�>t�>\)>O�>C�>+>o>o>   =��#=�=���=��=�x�=�`B=�`B=�`B=�l�=�`B=�G�=��=��=ȴ9=Ƨ�=ě�=��=�^5=�-=���=�1=���=��
=��=�t�=�hs=�O�=�+=��=�%=y�#=y�#=u=q��=m�h=e`B=e`B=]/=T��=L��=L��=L��=L��=D��=@�=<j=8Q�=��=�P=�P=�P=�P=��=#�
=#�
=,1=0 �=,1=��=��=,1=@�=Y�=u=�%=�o=�o=�+=�+=ix�=T��=L��=H�9=��=�P=�P=�P=��=�w=�w=�P=�P=t�=+=+<��<�<�<�h<�h<�h<�`B<�`B<�/<�/<���<���<ě�<�j<�1<���<�t�<�C�<�o<e`B<t�;�`B;�`B;ě�;ě�;ě�;��
;��
;D���o��o��o���
��`B�o�t��#�
�49X�49X�D���D���e`B��o��o��o��C���1��9X��j��/���o�C��+�+�C��t����'49X�@��H�9�H�9�L�ͽT���Y��]/�}󶽇+��+��7L��C���O߽�\)��hs���P�������w�������T��1��-��E���Q콼j�ȴ9��
=��l�������#�%�$ݾC��bN�bN�n���R�#�
�%�T�(�þ+�,1�,1�-V�-V�0 ž1&�49X�6E��;dZ�?|�B�\�D���H�9�J���Kƨ�M��P�`�S�Ͼ["Ѿ\(��^5?�_;d�bMӾcS��fff�l�D�n���o���z�H�{�m�}󶾀����7�������9��ƨ��O߾�����;��녾������+���u��"Ѿ�/���R��Ĝ��S���ff��1�������׾�-���!��Q쾺^5��^5��^5���#���#���#��X���#��^5��푾�vɾ�  ��%�\�����+��1'�ȴ9��7L�ɺ^��=q��C���I������V��\)���;��bN���`���`����t���t���t���t���t���t���t����Ͼ�z������z��z���Ͼ�����z���Ձ��
=��
=��
=��
=�և+�׍P�ؓu�ۥ�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��A�\A"�\AB�\Ab�\A�G�A�G�A�G�A�G�A�G�A�{A�G�A�G�B
=B��B��B��B ��B(��B0��B8��B@��BH��BP��BX��B`��Bh��Bp��Bx��B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�B�Q�C (�C(�C(�C(�C(�C
(�C(�C(�C(�C(�C(�C(�C(�C(�C(�C(�C (�C"(�C$(�C&(�C((�C*(�C,(�C.(�C0(�C2(�C4(�C6(�C8(�C:(�C<(�C>(�C@(�CB(�CD(�CF(�CH(�CJ(�CL(�CN(�CP(�CR(�CT(�CV(�CX(�CZ(�C\(�C^(�C`(�Cb(�Cd(�Cf(�Ch(�Cj(�Cl(�Cn(�Cp(�Cr(�Ct(�Cv(�Cx(�Cz(�C|(�C~(�C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�!HC�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{C�{D 
=D �=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D	
=D	�=D

=D
�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D
=D�=D 
=D �=D!
=D!�=D"
=D"�=D#
=D#��D$
=D$�=D%
=D%�=D&
=D&�=D'
=D'�=D(
=D(�=D)
=D)�=D*
=D*�=D+
=D+�=D,
=D,�=D-
=D-�=D.
=D.�=D/
=D/�=D0
=D0�=D1
=D1�=D2
=D2�=D3
=D3�=D4
=D4�=D5
=D5�=D6
=D6�=D7
=D7�=D8
=D8�=D9
=D9�=D:
=D:�=D;
=D;�=D<
=D<�=D=
=D=�=D>
=D>�=D?
=D?�=D@
=D@�=DA
=DA�=DB
=DB�=DC
=DC�=DD
=DD�=DE
=DE�=DF
=DF�=DG
=DG�=DH
=DH�=DI
=DI�=DJ
=DJ�=DK
=DK�=DL
=DL�=DM
=DM�=DN
=DN�=DO
=DO�=DP
=DP�=DQ
=DQ�=DR
=DR�=DS
=DS�=DT
=DT�=DU
=DU�=DV
=DV�=DW
=DW�=DX
=DX�=DY
=DY�=DZ
=DZ�=D[
=D[�=D\
=D\�=D]
=D]�=D^
=D^�=D_
=D_�=D`
=D`�=Da�Da��Db�Db��Dc
=Dc�=Dd
=Dd�=De
=De�=Df
=Df�=Dg
=Dg�=Dh
=Dh�=Di
=Di�=Dj
=Dj�=Dk
=Dk�=Dl
=Dl�=Dm
=Dm�=Dn
=Dn�=Do
=Do�=Dp
=Dp�=Dq
=Dq�=Dr
=Dr�=Ds
=Ds�=Dt
=Dt�=Du
=Du�=Dv
=Dv�=Dw
=Dw��Dx
=Dx�=Dy
=Dy��Dz�B
�gB
�dB
�cB
�dB
�XB
�lB
�lB
�`B
�UB
�`B
�aB
�aB
�aB
�aB
�aB
�lB
�^B
�\B
��B
��B
��B
��B
��B
��B
��B
�B
��B
��B
�)B
�zB
� B
��B
��B
�B
�B
��B
�B
�|B
��B
ΰB
�OB
��B
�7B
�,B�B�BSBB#�B3�B?ZBHBO�BR3BU BX�Bb�BhBihBjmBjBjBj�Bk�BkXBl(Bl)BkWBl1Bl�BmBl�Bm�Bm�BmuBnBqBq4Bp�Br�Bs�BsHBsBq�Bq�BsBq�Bq�Br�BsTBtcBvBw!By7BzLB{B|OB~[B�B�B�WB�B�PB��B��B��B�%B�B��B�5B��B��B��B��B�cBB|-B|(B}�B~dB~wB~B|uB{B{�B{�B{
Bz�Bz�Bz�B{�B}B}�B~B~B~B~B~B~BB~�B~�BB~�BB
B�B�B�B�B�B�B��B��B�B�B�tB�0B�B��B��B�B�B�B��B�B�B�B�-B�]B�>B�xB�pB�KB�4B�bB�KB�'B�>B�&B�1B�&B�1B�'B�6B�aB��B�MB��B�RB��B�cB�4B�aB�KB�KB�dB�dB�8B�_B�aB�UB�4B�bB�nB�VB�CB�DB�9B�QB�^B�vB��B�lB�WB�WB�dB�pB�|B�zB�AB�WB�qB��B�YB�WB�bB�qB�XB�cB�fB�OB�]B�]B�]B�kB�VB�oB�rB�rB�]B�\B�]B�tB�iB�jB�pB��B�qB�cB�bB�fB�\B�SB�mB�WB�cB�B��B�sB�FB�9B�2B�+B�aB�zB��B�rB��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B�B�!B��B��B��B��B��B��B��B��B��B�NB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B�CB��B��B��B��B��B�;B�
B��B��B��B��B�B��B��B��B��B��B��B��B�B�DB��B��B��B��B�BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��@Q�@Q�@Q�@Q�@Q�@bN@Q�@A�@A�@Q�@Q�@Q�@Q�@Q�@Q�@Q�@A�@Q�@�y@v�@@p�@��@��@33@�\@��@%@bN@\)@�T@�@j@1@"�@
M�@	�#@	&�@�@�j@��@ r�?�x�?�J?��?�?�G�?�I�?��T?˅?ě�?���?�Ĝ?�33?�;d?���?U?}?,��?!��?�?��?�j?hs?\)?	x�?�?��>��>�9X>�V>�1>��T>��
>���>���>��
>�~�>�`B>�;d>�A�>�Ĝ>�
=>��`>�O�>�I�>�C�>��m>��#>�dZ>�j>��>��>�I�>�t�>��>�A�>�Z>�>��>�9X>�F>��>�`B>�M�>��T>�x�>�h>�V>�V>�->��>�h>ݲ->�>��7>�r�>��>�\)>�b>���>�>��`>���>}�>y�#>j~�>\(�>Y�>Z�>["�>["�>["�>Xb>Xb>V>T��>Q�>O�;>M��>Kƨ>I�^>I�^>I�^>H�9>H�9>G�>F��>E��>E��>D��>C��>A�7>A�7>D��>F��>E��>D��>;dZ>8Q�>8Q�>?|�>F��>F��>F��>F��>J��>J��>J��>J��>I�^>D��>B�\>;dZ>5?}>2->1&�>,1>(��>(��>&�y>&�y>%�T>%�T>$�/>$�/>#�
>�R>+>(��>!��>�w>�P>t�>t�>\)>O�>C�>+>o>o>   =��#=�=���=��=�x�=�`B=�`B=�`B=�l�=�`B=�G�=��=��=ȴ9=Ƨ�=ě�=��=�^5=�-=���=�1=���=��
=��=�t�=�hs=�O�=�+=��=�%=y�#=y�#=u=q��=m�h=e`B=e`B=]/=T��=L��=L��=L��=L��=D��=@�=<j=8Q�=��=�P=�P=�P=�P=��=#�
=#�
=,1=0 �=,1=��=��=,1=@�=Y�=u=�%=�o=�o=�+=�+=ix�=T��=L��=H�9=��=�P=�P=�P=��=�w=�w=�P=�P=t�=+=+<��<�<�<�h<�h<�h<�`B<�`B<�/<�/<���<���<ě�<�j<�1<���<�t�<�C�<�o<e`B<t�;�`B;�`B;ě�;ě�;ě�;��
;��
;D���o��o��o���
��`B�o�t��#�
�49X�49X�D���D���e`B��o��o��o��C���1��9X��j��/���o�C��+�+�C��t����'49X�@��H�9�H�9�L�ͽT���Y��]/�}󶽇+��+��7L��C���O߽�\)��hs���P�������w�������T��1��-��E���Q콼j�ȴ9��
=��l�������#�%�$ݾC��bN�bN�n���R�#�
�%�T�(�þ+�,1�,1�-V�-V�0 ž1&�49X�6E��;dZ�?|�B�\�D���H�9�J���Kƨ�M��P�`�S�Ͼ["Ѿ\(��^5?�_;d�bMӾcS��fff�l�D�n���o���z�H�{�m�}󶾀����7�������9��ƨ��O߾�����;��녾������+���u��"Ѿ�/���R��Ĝ��S���ff��1�������׾�-���!��Q쾺^5��^5��^5���#���#���#��X���#��^5��푾�vɾ�  ��%�\�����+��1'�ȴ9��7L�ɺ^��=q��C���I������V��\)���;��bN���`���`����t���t���t���t���t���t���t����Ͼ�z������z��z���Ͼ�����z���Ձ��
=��
=��
=��
=�և+�׍P�ؓu�ۥ�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#�n<#�)<#�o<#��<#��<#��<#��<#׈<#��<#��<#��<#��<#��<#��<#��<#�y<$ǟ<#�v<#��<$$<$V<$6�<$.-<$c<$Gh<$�<$�<$h�<$��<$T�<$?<#�<$H�<$9.<#�t<$*�<$��<(c�<$��<(�<<*F�<*�<(v�<%��<(��<&�!<(�C</�,<)_ <*�5<3��<8q�<Kp|<8�x<Xu�<L&l<'��<$�<<$�<$8<$@�<$�<$��<$��<$B�<$C�<$��<$M{<#��<$1�<#�#<#��<#��<#�<$K<$�<$.<#�N<#�x<$��<$6T<#�9<#ߙ<#�<%��<#�t<#�:<#�<#�W<#�<$3K<$<$�<$�<#�%<$<$�<#�?<#ڼ<#�<$�<#��<#�<#�	<#�\<#�B<#�j<#�?<#ڥ<$;<%�<$v1<'�<(��<'��<#�><$2�<#۔<$<$!<%n<$+�<#��<$a�<$M0<#�]<#ת<#ג<#�)<#�4<#�2<#��<#ޱ<#۴<#�<#��<#��<#�+<#��<#�!<#�<#�Q<#�<#�Q<#�
<#�<#�<#�T<#ۂ<#�l<#��<#ג<#�<#�+<#�@<$%<#�<#�Q<#�<#��<#ث<#��<#�<#��<#��<#�g<#ؠ<#ۼ<#�<#�s<#�p<#�a<#�<#��<#�h<#��<#�A<#ށ<#�"<#�1<#�<#�?<#�C<#�<#�<#��<#�q<#� <#�P<$-<#��<#�]<#�P<#�
<#�
<#�f<#�<#��<#�<#�{<#�k<#��<#�<#��<#�|<#ؿ<#��<#ן<#�g<#��<#�?<#��<#��<#�_<#�z<#��<#�<#�K<#�g<#��<#�c<#�<#�O<#��<#�v<#�h<#�k<#ێ<#ޜ<#�_<#��<#�\<#�]<#�t<#�(<#��<#ޠ<#ߋ<#ݛ<#�
<#��<#�<#�T<#�q<#ۥ<#۰<#��<#��<#��<#��<#ؔ<#ׂ<#�<#ؗ<#�<#׌<#�#<#�6<#ج<#؜<#�,<#�<#�}<#ן<#�z<#ؽ<#�<#ٕ<$�<#�<#�*<#�k<$!:<#݊<#ؑ<#��<#׀<#�d<#�x<#�m<#�=<#�5<#�7<#�E<#�Z<#�7<#�	<#�<#��<#��<#�p<#ب<#�<#�<#�N<#�N<#�O<#�h<#ޚ<#ކ<#�d<#��<#چ<#�<#�b<#��<#�<#�(<#��<#��<#�<#�<#�<#�?<#��<#�<#�Z<#��<#��<#�J<#�I<#�1<#�<#�7<#�<#ޅ<#�]<#�<#��<#ۃ<#�<#�
<#��<#��<#�f<#ު<#�9<#נ<#��<#�a<#މ<#��<#��<#��<#��<#�x<#�<#�M<#�{<#�]<#� <$�<#�<#�<#�.<#�D<#�D<#�F<#�w<#�D<#ۦ<#�t<#�s<#ޠ<#�<#�<#�y<#�s<#�<#��<#��<$x<#�<#��<#��<#�g<#�H<#�h<#�G<#�/<#��<$-�<#�_<#�<#�<#ގ<#�B<#�<#� <#�<#�&<#۪<#�<#�<#��<#�-<#��<#��<#�<#޳<#�u<#ޤ<#��<#�<#��<#��<#�j<#�}<#�I<#۫<#�)<#�<#��<#�g<$"�<#�|<#ޫ<#�<#��<#�!<$�<#��<#�g<#�<#��<#�m<#�e<#�<#�A<#�5<#�#<#�A<#�<<#�<#��<$#�<#�<#�!<#�<#�v<$!<#�<#��<#��<#׈<#��<#��<#׍<#�<#ۭ<#��<#�<#�<#��<#��<#�e<#��<#��<#�V<#�<<#�=<#�V<#�w<#޴<#�r<#�U<#�Y<#�T<#�:<#�$<#�(<#��<#�l<#��<#��<#��<#��<#��<#��<#�.<#�6<#�<#׊<#��<#ו<#��<#ײ<#��<#ק<#��<#�<#��<#��<#ז<#�(<#޳<#�:<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.16                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151607201909291516072019092915161120190929151617202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�  @�  @�  G�O�G�O�G�O�G�O�Dy��Dy��Dy��G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          