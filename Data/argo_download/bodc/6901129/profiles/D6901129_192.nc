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
_FillValue                 �  Vl   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X`   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ZT   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  \H   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  d   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  k�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  u�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  w|   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  yp   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �p   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �    HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225044323  20210225044323  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125479                          2C  D   APEX                            6229                            120210                          846 @�#M�[ 1   @�#M�[ @P�I�^�4.��+1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�ff@���A   A   A@  A`  A�  A�  A�  A�  A���A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C�fC  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CI�fCK�fCN  CP  CR  CT  CV  CX�CZ  C\  C^  C`  Cb�Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C��3C�  C��C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D �fD  D� D  D� D  D� D  D� D  D� D  D� D  D� DfD� D	  D	� D
  D
� D  Dy�D  D�fD  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D��D� D  D� D  D� D  D� DfD� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!y�D!��D"� D#  D#� D$  D$� D%  D%�fD&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/�fD0  D0y�D1  D1� D2  D2y�D3  D3� D4  D4� D4��D5� D6  D6� D7  D7� D7��D8� D9fD9�fD:  D:� D;  D;� D<  D<� D<��D=y�D>  D>�fD?fD?� D@  D@� DA  DAy�DB  DB� DCfDC� DC��DD� DE  DE�fDFfDF� DG  DG� DH  DH� DH��DI� DJ  DJ�fDK  DKy�DK��DL� DMfDM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DUy�DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]fD]� D^  D^� D_  D_�fD`  D`� Da  Day�Da��Db� DcfDc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� DkfDk�fDl  Dly�Dm  Dm� Dn  Dn�fDofDo� Do��Dp� Dq  Dq� DrfDr� Dr��Dsy�Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
��B
��B
��B
�B
�FB
�jB
�B{B!�B/B49BC�BiyBu�Bt�Bm�BbNB`BB^5B^5B_;B\)BVB@�B9XBF�BA�BD�BD�B>wB9XB5?B33B0!B.B-B-B0!B33B8RB8RB9XB:^B:^B7LB7LB8RB;dB:^B9XB8RB9XB<jB@�BB�BE�BI�BJ�BM�BO�BR�BT�BXB[#B^5B`BBdZBgmBgmBhsBiyBiyBiyBjBiyBhsBffBffBgmBgmBe`BdZBe`BffBffBhsBiyBjBjBk�Bk�BjBjBjBk�Bk�Bm�Bm�Bl�Bl�Bl�Bl�Bm�Bm�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bo�Bo�Bp�Bp�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bn�Bm�Bm�Bm�Bm�Bl�Bm�Bm�Bm�Bm�Bn�Bo�Bo�Bo�Bo�Bo�Bn�Bn�Bn�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bm�Bm�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Br�Br�Br�Br�Br�Br�Bs�Bs�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bv�Bw�Bx�Bx�Bx�Bx�By�By�By�Bz�B{�B{�B|�B|�B}�B}�B~�B~�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�+B�+B�1B�1B�1B�7B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�JB�PB�VB�bB�hB�oB�uB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�`B@�hs@�p�@�x�@�x�@�x�@�x�@�@�@�x�@�x�@�O�@��`@���@�p�@�=q@���@�S�@�Q�@��D@��@�-@�G�@��m@�1'@���@��T@���@�/@y��@w�;@so@i�@W|�@B�@G�@AX@4Z@/;d@$9X@��@A�@��@�-@Ĝ@��?�Q�?���?���?��?և+?Ұ!?�;d?ǍP?�S�?���?���?�I�?�^5?�$�?��?�33?��?�G�?�G�?�M�?�J?�A�?xQ�?t�j?lI�?e�T?c�
?b�\?`A�?]�-?X�u?]/?Z�?X��?X��?Q&�?PbN?Q&�?PbN?J��?D��?9��?6�+?9�#?%��?ȴ??�!?&�?��?	�^?��?r�?��?��?��?M�?%>�X>�F>>�{>�h>�>�ff>���>�"�>׍P>�>���>�n�>��>ɺ^>�7L>��>���>��>�Q�>�9X>��!>���>�V>���>�A�>�>�t�>��`>��>�O�>��\>}�>z�H>dZ>o��>dZ>Y�>`A�>dZ>^5?>W
=>O�;>F��>B�\>>v�>7K�>+>�=��#==�`B=���=��`=\=�j=��=��
=�7L=}�=u=m�h=e`B=}�=]/=Y�=Y�=T��=0 �='�='�=�P=��=�w=�w=��=�w=�P=\)=C�=+=o=+=+<��<���<�9X<�9X<�9X<�1<���<�C�<�C�<�o<�o<e`B<e`B<e`B<T��<D��;�`B;D��;�o;D��;o    ;o<49X<�j<���<�/<u<#�
<t�<D��<t�<���<��
<�C�<�1<���<�9X<���<�`B<�h=+=�P=��=�w=,1=,1='�=�P=,1=@�=T��=]/=e`B=aG�=@�=D��=Y�=]/=L��=P�`=@�=8Q�=C�<�`B<�h<�h<�=o=\)=t�=�P=,1=,1=8Q�=<j=<j=D��=L��=L��=L��=L��=L��=L��=L��=e`B=�o=��=�j=��`=�l�=�x�=�S�=��=�
==���=��=��=��=ȴ9=ě�=\=��=��=��=\=\=��=��=�v�=�j=�^5=�{=���=�O�=�C�=�+=��=�%=}�=ix�=P�`=@�=�w=t�=�P=\)<���<��
<u<e`B<#�
;�`B;��
;o��o�D����o�ě��t��49X�T���e`B�e`B��C����
��1��9X�ě���/��`B���o�+�+�C��\)�t���P�����#�
�#�
�#�
�0 Ž0 Ž���P���#�
�0 Ž<j�L�ͽaG��m�h�y�#�����C���\)��hs��t����T��1��{��9X��E���j��vɽ���ě��\�Ƨ���`����
=��������/��`B��l���l���xս�h������F���ٽ��m�%�o�+�O߾bN�z��+��P�������R�!���!���!���"��"��%�T�)��/��2-�333�6E��8Q�;dZ�>vɾ@��C���E�˾G��I�^�J���L�;M��O�;�Q녾V�Z��^5?�`A��aG��aG��bMӾbMӾcS��cS��fff�ixվm�h�q���t�j�x���{�m��%������������I����;�O߾�����;��hs��񪾖���+��
=��b����������������/��;d���w��A���Ĝ������S���Z��`B���T��ff���y���xվ�1��{��&龳�F��ȴ��KǾ�KǾ�KǾ�KǾ�KǾ�ȴ��ȴ��ȴ��ȴ��ȴ��ȴ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��
@�=qA�RA"�RAB�RAb�RA�\)A�\)A�\)A�\)A�(�A�\)A�\)A�\)B �B�B�B�B �B(�B0�B8�B@�BH�BP�BX�B`�Bh�Bp�Bx�B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
C +�C+�C+�C+�C+�C
+�C+�C+�C�C+�C+�C+�C+�C+�C+�C+�C +�C"+�C$+�C&+�C(+�C*+�C,+�C.+�C0+�C2+�C4+�C6+�C8+�C:+�C<+�C>+�C@+�CB+�CD+�CF+�CH+�CJ�CL�CN+�CP+�CR+�CT+�CV+�CXECZ+�C\+�C^+�C`+�CbECd+�Cf+�Ch+�Cj+�Cl+�Cn+�Cp+�Cr+�Ct+�Cv+�Cx+�Cz+�C|+�C~+�C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C�"�C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C�"�C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D 
�D �GD
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��DGD��D	
�D	��D

�D
��D
�D�{D
�D�GD
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D{D��D
�D��D
�D��D
�D��DGD��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D 
�D ��D!
�D!�{D"{D"��D#
�D#��D$
�D$��D%
�D%�GD&
�D&��D'
�D'��D(
�D(��D)
�D)��D*
�D*��D+
�D+��D,
�D,��D-
�D-��D.
�D.��D/
�D/�GD0
�D0�{D1
�D1��D2
�D2�{D3
�D3��D4
�D4��D5{D5��D6
�D6��D7
�D7��D8{D8��D9GD9�GD:
�D:��D;
�D;��D<
�D<��D={D=�{D>
�D>�GD?GD?��D@
�D@��DA
�DA�{DB
�DB��DCGDC��DD{DD��DE
�DE�GDFGDF��DG
�DG��DH
�DH��DI{DI��DJ
�DJ�GDK
�DK�{DL{DL��DMGDM��DN
�DN��DO
�DO��DP
�DP��DQ
�DQ��DR
�DR��DS
�DS��DT
�DT��DU
�DU�{DV
�DV��DW
�DW��DX
�DX��DY
�DY��DZ
�DZ��D[
�D[��D\
�D\��D]GD]��D^
�D^��D_
�D_�GD`
�D`��Da
�Da�{Db{Db��DcGDc��Dd
�Dd��De
�De��Df
�Df��Dg
�Dg��Dh
�Dh��Di
�Di��Dj
�Dj��DkGDk�GDl
�Dl�{Dm
�Dm��Dn
�Dn�GDoGDo��Dp{Dp��Dq
�Dq��DrGDr��Ds{Ds�{Dt
�Dt��Du
�Du��Dv
�Dv��Dw
�Dw��Dx
�Dx��Dy*�B
� B
��B
� B
�
B
�B
�B
��B
�B
�B
�B
�@B
��B
�CB
�7B
��B
�,B
�+B-B'/B2�B9|BEBknBz�B|�By�Bf�Bh2Bd2B_�Bb�Bb�BcBN�B7+BKBJ�BH�BLjBC\B=]B9�B6rB3�B2�B1�B.�B6B6�B9�B9�B:�B=_BA�B=B9nB9?B<#B<B<�B;�B9�B<�B@~BB7BE�BJkBLGBN�BQ�BTBUaBXRB[�B^�Ba,Bc�Bg�Bg�Bh�Bj�Bi�Bi]Bj�Bj�Bi�BhoBgBf�Bk>Bh&Bd�Be�Bf�BgBi2Bi�Bj�Bj�Bk�Bk�Bj�Bj�BkQBlBk�Bm�Bm�Bl�BmBl�BmMBm�Bm�Bm�Bm�Bm�BnBn�Bn�Bo,Bo�BpBqBp�Bo�Bo�BpBonBo�Bn�Bn�Bn�Bm�Bn�Bm�Bm�Bm�BmBnBnBmFBnjBo�Bo�Bo�BpBo�Bn�Bn�Bo/Bn�Bm�Bl�Bl�Bl�Bl�Bl�Bl�BmBl�Bl&Bl�Bl�Bl�Bl�BmLBm�Bn�Bn�Bn�BoBn�Bn�Bn�Bn�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Br�Br�Br�Br�Br�Br�Bs�Bs�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�BvXBw6Bx�Bx�ByfByBy�By�By�BzjB{�B|B|�B}B}�B}�B~�B~�B�B�B��B��B��B�B�B�:B��B��B��B�B�B�.B��B�B��B�B�UB�B�UB�CB��B�nB�B�$B�B�
B�B�!B� B��B�0B�B�-B�<B�%B�&B�AB�DB�DB�DB�DB�DB� B��B�zB��B��B��B�kB��B��B��B��B��B�}B��B��B��B��B��B��B��B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�[B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�`B@�hs@�p�@�x�@�x�@�x�@�x�@�@�@�x�@�x�@�O�@��`@���@�p�@�=q@���@�S�@�Q�@��D@��@�-@�G�@��m@�1'@���@��T@���@�/@y��@w�;@so@i�@W|�@B�@G�@AX@4Z@/;d@$9X@��@A�@��@�-@Ĝ@��?�Q�?���?���?��?և+?Ұ!?�;d?ǍP?�S�?���?���?�I�?�^5?�$�?��?�33?��?�G�?�G�?�M�?�J?�A�?xQ�?t�j?lI�?e�T?c�
?b�\?`A�?]�-?X�u?]/?Z�?X��?X��?Q&�?PbN?Q&�?PbN?J��?D��?9��?6�+?9�#?%��?ȴ??�!?&�?��?	�^?��?r�?��?��?��?M�?%>�X>�F>>�{>�h>�>�ff>���>�"�>׍P>�>���>�n�>��>ɺ^>�7L>��>���>��>�Q�>�9X>��!>���>�V>���>�A�>�>�t�>��`>��>�O�>��\>}�>z�H>dZ>o��>dZ>Y�>`A�>dZ>^5?>W
=>O�;>F��>B�\>>v�>7K�>+>�=��#==�`B=���=��`=\=�j=��=��
=�7L=}�=u=m�h=e`B=}�=]/=Y�=Y�=T��=0 �='�='�=�P=��=�w=�w=��=�w=�P=\)=C�=+=o=+=+<��<���<�9X<�9X<�9X<�1<���<�C�<�C�<�o<�o<e`B<e`B<e`B<T��<D��;�`B;D��;�o;D��;o    ;o<49X<�j<���<�/<u<#�
<t�<D��<t�<���<��
<�C�<�1<���<�9X<���<�`B<�h=+=�P=��=�w=,1=,1='�=�P=,1=@�=T��=]/=e`B=aG�=@�=D��=Y�=]/=L��=P�`=@�=8Q�=C�<�`B<�h<�h<�=o=\)=t�=�P=,1=,1=8Q�=<j=<j=D��=L��=L��=L��=L��=L��=L��=L��=e`B=�o=��=�j=��`=�l�=�x�=�S�=��=�
==���=��=��=��=ȴ9=ě�=\=��=��=��=\=\=��=��=�v�=�j=�^5=�{=���=�O�=�C�=�+=��=�%=}�=ix�=P�`=@�=�w=t�=�P=\)<���<��
<u<e`B<#�
;�`B;��
;o��o�D����o�ě��t��49X�T���e`B�e`B��C����
��1��9X�ě���/��`B���o�+�+�C��\)�t���P�����#�
�#�
�#�
�0 Ž0 Ž���P���#�
�0 Ž<j�L�ͽaG��m�h�y�#�����C���\)��hs��t����T��1��{��9X��E���j��vɽ���ě��\�Ƨ���`����
=��������/��`B��l���l���xս�h������F���ٽ��m�%�o�+�O߾bN�z��+��P�������R�!���!���!���"��"��%�T�)��/��2-�333�6E��8Q�;dZ�>vɾ@��C���E�˾G��I�^�J���L�;M��O�;�Q녾V�Z��^5?�`A��aG��aG��bMӾbMӾcS��cS��fff�ixվm�h�q���t�j�x���{�m��%������������I����;�O߾�����;��hs��񪾖���+��
=��b����������������/��;d���w��A���Ĝ������S���Z��`B���T��ff���y���xվ�1��{��&龳�F��ȴ��KǾ�KǾ�KǾ�KǾ�KǾ�ȴ��ȴ��ȴ��ȴ��ȴ��ȴ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#�C<#�s<#�d<#ع<#؃<#�z<#ؔ<#�g<#ؙ<#�.<$M<2.H<?f�<hߓ<ψP<��<b�/<9�<.��<8ر<%�<'�<8�<P�]<}��<3��<P9m<>�H<&|�<-�_<E�<�^�<�e�<'2N<2��<]�<0��<N�<6?&<0��<4>q<,V+<.�<4L<4��<&`1<>xT<.z�<%}�<%�<%��<+3<LO�<<u�<'��<$��<$e�<&<�<-�y<-�A<$0�<#��<#�S<#�<#��<$T%<%�><$b�<&4�<%
2<$�<#�<$�<$�<$��<$8�<$�<#�w<#�.<%�<#��<#�<#��<$��<%'.<'^�<$=0<#�</m(<*-<#�|<$<�<#��<$G)<$g0<#��<#�n<$<#��<#��<$z<#��<$��<$'�<$�<#�|<#�<#��<$<#�g<$j�<#�<#�<#�<#�<#��<$5<#��<$�<$3�<#�P<$�<$�<#�<#��<#��<$	<$�g<$�t<#�<#�U<#��<#�<$�f<$<<#�<$ъ<#��<$&6<$"a<#�`<#�@<#�]<$8<$ <$5<#��<#�l<$ 8<$6y<$�?<%�<#��<#�a<$�<#�`<#��<#�<$p<#�<$A�<#�.<#߽<#�H<#ޠ<#܅<$><#��<#�f<#�W<$b<#�<#٫<#��<#��<#�|<#ټ<#�<#��<#�V<#߈<#��<#��<#ے<#��<#�Q<#߇<#�b<#�<#�<#��<#��<#�:<#�<#�r<#�<#�`<#�&<#��<#��<#۽<#�4<#�<#�b<#�<#��<#��<#ދ<#�<#�&<$�<#ך<#�%<$/�<#�<#�)<#�/<#�<#�t<#׏<#��<#�d<#ޝ<#�l<#�K<#ה<#ס<#؞<#�k<#�[<#��<#�m<#�J<#��<#�%<#ۡ<#��<#�<#�<#�<#�=<$T<#׹<#��<#�z<#�<#��<#��<#��<$%<#��<#��<#�<#׭<#�
<#�q<#��<#�M<#�a<#��<#�P<#�a<#��<#�<#�<#ة<#�,<#�,<#�,<#�'<#�2<#�<#��<$IU<#��<#��<#��<#��<#�3<#�X<#�
<#ۨ<#ې<#�_<#��<#�}<#��<#ۺ<#ۍ<#�7<#�<#׫<#�)<#�s<#�H<#�y<#۪<#�3<#�s<#��<$E<#��<#އ<#��<#��<#۪<#�,<#��<#��<$�<#��<#ף<#��<$;<#��<#�-<#�[<#�<#�<#��<#�]<#�<#��<#�U<#�8<#�6<#�%<#ޢ<#��<#�p<#�+<#�/<#��<#ۿ<#�<#�<#��<#��<#�<#ۚ<#�5<#ۏ<#ۡ<#ۡ<#۠<#ۇ<#�U<#�?<#�e<#�_<#�3<#��<#��<#�}<#ۓ<#�<#�V<#�<#�E<#�><#��<#�<#�<#�f<#�h<#�J<#�D<$�<#��<#�<#��<#��<#�<#��<#��<#ލ<#��<#��<#� <#ۭ<#��<#ی<#�W<#�v<#��<#ۀ<#�;<#��<#�%<#ۋ<#ۘ<#ێ<#��<#ߜ<#�N<#�p<#��<#��<#�<#�p<#�!<#�<#��<#�s<#�(<#�=<#�J<#�0<#�t<#�d<#�"<#��<#�<#�Z<#�<#��<#�;<#�a<#�H<#�N<#�)<#�<#�
<#�><#ۦ<#��<#ۺ<#��<#�O<#�<#��<#�g<#�<#��<#�1<#�n<#�:<#ڻ<#ٱ<#��<#�<#�<#�<#�<#�<#�'<#��<#�)<#�F<#�q<#��<#��<#��<#�!<#�H<#�H<#��<#�<#��<#ۭ<#�<#�<#ۋ<#�;<#��<#�4<#�;<#��<#۷<#��<#ޭ<#�+<#�<#ފ<#ۙ<#۷<#��<#ީ<#�<#�<#�8<#��<#�j<#��<#�C<#�<#�<#�<#�<#נ<#�<#ؤ<#�h<#�<#ٞ<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.17                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929160552201909291605522019092916055620190929160603202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�Dy  Dy  Dy  G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          