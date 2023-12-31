CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  Y   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
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
resolution        ?�������     	d  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  H�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  Q�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 \  [L   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 \  ]�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 \  `   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     	d  b`   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  k�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  u(   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 \  ~�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 \  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 \  �D   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     	d  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  �   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     	d  �h   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �,   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �,   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �,   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �,   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �0   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225035259  20210225035259  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               _A   BO  125383                          2C  D   APEX                            6229                            120210                          846 @�6Ax	 1   @�6Ax	 @Q���+�0�����1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   B   A   @y��@�  A   A!��A>ffA`  A���A�  A�  A�  A�  A���A�  A�  B   B��B��B  B   B(ffB0  B8  B@  BH  BO��BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C�C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4�C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�3D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�3D�C3D�� D�� D�  D�@ D�� D�� D�  D�@ D��3BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BJ�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BK�BL�BL�BL�BL�BL�BL�BM�BM�BM�BM�BL�BL�BL�BL�BL�BL�BL�BL�BN�BP�BM�BM�BN�BO�BP�BR�BT�B\)B\)B_;B^5B^5B_;B_;B_;B_;B_;B`BB`BBaHBcTBffBiyBk�BiyBiyBhsBhsBhsBl�Bs�Bz�B�B�B�+B�JB�hB�{B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�'B�-B�-B�3B�3B�3B�3B�9B�9B�9B�?B�?B�FB�FB�FB�FB�FB�LB�FB�FB�FB�FB�FB�FB�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�3B�3B�3B�-B�3B�3B�3B�3B�3B�3B�3B�3B�-B�-B�'B�'B�!B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�'B�!B�!B�'B�'B�!B�'B�'B�!B�!B�'B�'B�!B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�!B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B
�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�ݽ��`�����/��;d��/�����������
=��"ѽ�/������ͽ��ͽ�����ͽ��ͽ��ͽ��ͽ��ͽ��ͽƧ�ȴ9�ȴ9�ȴ9�ě���Q콸Q콰 Ž��
���w������P��������������O߽����+�����C���O߽�O߽�O߽��������������hs�����C���o��O߽�hs��7L�T���ixս}�P�`��+������\)�����+��+���-��E���C��o�e`B�}�e`B�@���j����;D��=H�9=�-=�F=�==�x�>   =��#=���=�l�>$�>$�>   >$�>.{>hr�>gl�>o��>Z�>_;d>P�`>O�;>Z�>��>���>޸R>�9X>�j?r�?��?"��?"��?0 �?:��?;dZ?<(�?;"�?>��?T�j?Y�?R�!?F��?BJ?F��?G�?Co?D��?Fff?E��?E�?D�?C�
?CS�?B��?BJ?A�7?B�\?B��?@Ĝ??�w?@A�?A�7?A%?@�?@�?@  ??�w??;d?<�?8b?5�?4�j?4�j?3��?0bN?-V?+�?(1'?&��?&��?&��?&ff?&$�?%�?"�\?!�7?|�?j?��?�!?��?�?
=q?
~�?
��?
~�?	�^?�?l�?$�>�^5>�!>�D>�r�>ܬ>�z�>�hs>�bN>�hs>��>���>���>�O�>��>�1'>Õ�>�+>��m>��H>�X>�K�>��j>��j>�>�?}>�?}>��!>��!>�33>��F>�&�>���>��h>�1>�1>�1>�~�>��T>���>��w>�Ĝ>��w>�/>���>��+>��>���>��9>���>���>���>�o>�o>��>���>��>���>��>�$�>���>��>���>��>�=q>��>��\>vȴ>ix�>fff>fff>cS�>]/>Y�>V>V>W
=>S��>S��>S��>R�>Q�>L��>H�9>F��>G�>C��>5?}>49X>333>1&�>1&�>/�>-V>)��>"��>�w>�->��>O�>
=q>+>$�>%>%=��m=���=���=���=�=�=�F=�F=�=�F=�F=��==�=�=�G�=�/=�"�=���=��=���=���=���=ȴ9=��=�E�=�-=� �=�1=��=��=���=���=��=�t�=�o=}�=�%=q��=@�=0 �='�=�P=t�=\)=C�=C�=+=+=+=+=C�=C�=+=o<��=o<�h<ě�<�1<���<u<T��<49X;��
�D���ě��u��t������\)�,1�H�9�]/�e`B�e`B�q���u�q���u�}󶽃o��+��+��C���\)��\)���P���
��j���ͽ�S���F���پ   ���C���+��u����������R��R�&�y�'.{�/��.{�/��9X�>vɾKƨ�L�;P�`�S�ϾY��Z��["Ѿ\(��]/�]/�`A��`A��aG��bMӾbMӾhr��l�D�n���p�׾vȴ�x����%���\����������������$ݾ�+��1'��1'���9��=q����������ƨ��ƨ��ƨ��ƨ��I�����t�����(����/��l���1������ ž������D���D��V������~���~������D��{�������������%�����+�ȴ9��7L��7L������ƨ���;�I����;�����=q��=q��C���O߾����bN��t��׍P��b�ؓu�ؓu��"Ѿ�5?��Ĝ������/��`B�������������D��1������V�� ž� ž� ž� ž� ž� ž�� ž� ž�&��&�����F���j���j����E���E���E���E���E���ȴ��E���E�����?}������ȴ��E���Q������dZ��dZ��dZ��dZ��dZ���H���H��^5���H��^5��dZ���   �   ��|��|���۾��۾��۾�|��|��|�   ���۾�|���۾�|�   ��|�   ��|��|���۾�vɾ��۾��ۿ ��%�%� Ĝ� Ĝ�   � �� �� �� �� �� �� �� �� �� �� �� �� Ĝ� Ĝ� Ĝ�%� �� Ĝ�%�%�%� Ĝ� Ĝ� Ĝ� Ĝ� �� �� �� Ĝ� �� �� A�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @u�@�@�A z�A=G�A^�HA�=qA�p�A�p�A�p�A�p�A�=qA�p�A�p�A�p�BQ�BQ�B�RB�RB(�B/�RB7�RB?�RBG�RBOQ�BW�RB_�RBg�RBo�RBw�RB�RB��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)C�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C4�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C��
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
C��=C��=C��
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
D {�D ��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D	{�D	��D
{�D
��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D {�D ��D!{�D!��D"{�D"��D#{�D#��D${�D$��D%{�D%��D&{�D&��D'{�D'��D({�D(��D){�D)��D*{�D*��D+{�D+��D,{�D,��D-{�D-��D.{�D.��D/{�D/��D0{�D0��D1{�D1��D2{�D2��D3{�D3��D4{�D4��D5{�D5��D6{�D6��D7{�D7��D8{�D8��D9{�D9��D:{�D:��D;{�D;��D<{�D<��D={�D=��D>{�D>��D?{�D?��D@{�D@��DA{�DA��DB{�DB��DC{�DC��DD{�DD��DE{�DE��DF{�DF��DG{�DG��DH{�DH��DI{�DI��DJ{�DJ��DK{�DK��DL{�DL��DM{�DM��DN{�DN��DO{�DO��DP{�DP��DQ{�DQ��DR{�DR��DS{�DS��DT{�DT��DU{�DU��DV{�DV��DW{�DW��DX{�DX��DY{�DY��DZ{�DZ��D[{�D[��D\{�D\��D]{�D]��D^{�D^��D_{�D_��D`{�D`��Da{�Da��Db{�Db��Dc{�Dc��Dd{�Dd��De{�De��Df{�Df��Dg{�Dg��Dh{�Dh��Di{�Di��Dj{�Dj��Dk{�Dk��Dl{�Dl��Dm{�Dm��Dn{�Dn��Do{�Do��Dp{�Dp��Dq{�Dq��Dr{�Dr��Ds{�Ds��Dt{�Dt��Du{�Du��Dv{�Dv��Dw{�Dw��Dx{�Dx��Dy{�Dy��Dz{�Dz��D{{�D{��D|{�D|��D}{�D}��D~{�D~��D{�D��D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D� �D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D���D�=�D�}�D���D� �D�@�D�}�D���D���D�=�D�}�D���D���D�=�D���BK�BK�BK�BK�BK�BK�BK�BK�BK�BK�BKxBK�BK�BK�BK�BK�BJ�BK�BK�BK�BK�BK�BK�BK�BK�BKBK�BK�BK~BK�BK�BK�BLBK�BKtBL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BL�BK�BL�BL�BM
BL�BL�BLBNBNBM[BN�BMjBLaBL�BL�BL�BMYBMVBK�BM3BQ�BNBM�BNbBN�BP�BQ�BR�BZ7BZ�B_aB^#B^VB^�B_aB_LB_�B^mB`AB`�B`�BasBc�Bi�Bk(BjpBiJBiBh}Bg�BilBpWBv�B~�B�;B�CB��B��B�MB�$B��B��B��B��B�B��B�:B�@B�?B��B�.B��B��B��B��B�OB�JB�`B�AB�LB�OB�[B�PB�B�5B��B�tB�.B�B�]B�aB�IB�`B�SB�bB��B�!B��B�lB�IB�mB��B��B��B��B��B�?B�@B�KB�MB�rB��B�oB��B��B��B��B��B�!B�qB�(B�(B�@B�ZB��B�DB�zB��B��B��B��B�;B��B�kB�:B�B�ZB�MB�B�B��B��B��B��B�B�:B�;B�IB�NB�B��B� B�B�OB�B�B�
B�OB�%B�]B�8B�B�B�8B�~B�qB�?B��B�%B�FB�2B��B�wB�iB�~B�MB�B�B�B�B��B��B��B�B��B��B�B�B��B��B��B�IB��B��B��B�8B�B�6B�WB�AB�3B�B�B�1B�B�B�B�B�KB�@B�(B�	B�FB��B�!B�!B�,B�B�-B�.B�;B�iB�;B�.B�>B��B�BB�@B�*B�VB�B�@B�3B�B�B�'B�B�'B�B�B�'B�B�'B�(B�3B�B�YB�8B�/B�FB�FB�.B�-B�"B�/B�SB�^B�;B�.B�:B�;B�$B�IB�TB�?B�5B��B�?B�B�WB��B�YB�@B�TB�2B�4B�3B�(B�3B�(B�'B�'B�B�'B�3B�4B�3B�B�IB�dB�LB�@B�TB�;B�<B�mB��B�LB��B�IB��B��B�tB�rB�YB�4B�B�?B�'B�B�'B�4B�4B�3B�B�3B�4B�B�MB�dB��B�{B��B�uB�.B�CB�CB�hB��B�+B�B�B�B�2B�B�jB�B�PB�B��B�B��B�JB��B�B�:B�/B�BB�B�B�B�B�B�&B��B�
B�B�B�MB�4B�B�B�HB�B�lB�'B�B�	B�	B�	B�B�B�B��B�	B�!B�B��B�B��B��B��B�B�"B�|B�9B��B��B�9B�aB�FB�B��B��B��B��B�B��B��B��B�B�B�B�3B��B��B�9B�4B�B��B��B�B�B�B��B��B��B��B��B�B�B�B�B�4B�HB��B��B��B�!B�.B�"B�B�B��B�"B�-B��B��B��B��B�B��B��B��B�B�(B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B��B��B�B��B�B��B��B��B��B��B��B��B��B��B��B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�ܽ��`�����/��;d��/�����������
=��"ѽ�/������ͽ��ͽ�����ͽ��ͽ��ͽ��ͽ��ͽ��ͽƧ�ȴ9�ȴ9�ȴ9�ě���Q콸Q콰 Ž��
���w������P��������������O߽����+�����C���O߽�O߽�O߽��������������hs�����C���o��O߽�hs��7L�T���ixս}�P�`��+������\)�����+��+���-��E���C��o�e`B�}�e`B�@���j����;D��=H�9=�-=�F=�==�x�>   =��#=���=�l�>$�>$�>   >$�>.{>hr�>gl�>o��>Z�>_;d>P�`>O�;>Z�>��>���>޸R>�9X>�j?r�?��?"��?"��?0 �?:��?;dZ?<(�?;"�?>��?T�j?Y�?R�!?F��?BJ?F��?G�?Co?D��?Fff?E��?E�?D�?C�
?CS�?B��?BJ?A�7?B�\?B��?@Ĝ??�w?@A�?A�7?A%?@�?@�?@  ??�w??;d?<�?8b?5�?4�j?4�j?3��?0bN?-V?+�?(1'?&��?&��?&��?&ff?&$�?%�?"�\?!�7?|�?j?��?�!?��?�?
=q?
~�?
��?
~�?	�^?�?l�?$�>�^5>�!>�D>�r�>ܬ>�z�>�hs>�bN>�hs>��>���>���>�O�>��>�1'>Õ�>�+>��m>��H>�X>�K�>��j>��j>�>�?}>�?}>��!>��!>�33>��F>�&�>���>��h>�1>�1>�1>�~�>��T>���>��w>�Ĝ>��w>�/>���>��+>��>���>��9>���>���>���>�o>�o>��>���>��>���>��>�$�>���>��>���>��>�=q>��>��\>vȴ>ix�>fff>fff>cS�>]/>Y�>V>V>W
=>S��>S��>S��>R�>Q�>L��>H�9>F��>G�>C��>5?}>49X>333>1&�>1&�>/�>-V>)��>"��>�w>�->��>O�>
=q>+>$�>%>%=��m=���=���=���=�=�=�F=�F=�=�F=�F=��==�=�=�G�=�/=�"�=���=��=���=���=���=ȴ9=��=�E�=�-=� �=�1=��=��=���=���=��=�t�=�o=}�=�%=q��=@�=0 �='�=�P=t�=\)=C�=C�=+=+=+=+=C�=C�=+=o<��=o<�h<ě�<�1<���<u<T��<49X;��
�D���ě��u��t������\)�,1�H�9�]/�e`B�e`B�q���u�q���u�}󶽃o��+��+��C���\)��\)���P���
��j���ͽ�S���F���پ   ���C���+��u����������R��R�&�y�'.{�/��.{�/��9X�>vɾKƨ�L�;P�`�S�ϾY��Z��["Ѿ\(��]/�]/�`A��`A��aG��bMӾbMӾhr��l�D�n���p�׾vȴ�x����%���\����������������$ݾ�+��1'��1'���9��=q����������ƨ��ƨ��ƨ��ƨ��I�����t�����(����/��l���1������ ž������D���D��V������~���~������D��{�������������%�����+�ȴ9��7L��7L������ƨ���;�I����;�����=q��=q��C���O߾����bN��t��׍P��b�ؓu�ؓu��"Ѿ�5?��Ĝ������/��`B�������������D��1������V�� ž� ž� ž� ž� ž� ž�� ž� ž�&��&�����F���j���j����E���E���E���E���E���ȴ��E���E�����?}������ȴ��E���Q������dZ��dZ��dZ��dZ��dZ���H���H��^5���H��^5��dZ���   �   ��|��|���۾��۾��۾�|��|��|�   ���۾�|���۾�|�   ��|�   ��|��|���۾�vɾ��۾��ۿ ��%�%� Ĝ� Ĝ�   � �� �� �� �� �� �� �� �� �� �� �� �� Ĝ� Ĝ� Ĝ�%� �� Ĝ�%�%�%� Ĝ� Ĝ� Ĝ� Ĝ� �� �� �� Ĝ� �� �� A�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<#פ<#�
<#��<#��<#�p<#�<#�<#׏<#�<#�<#��<#�n<#�
<#ص<#מ<#�o<#ס<#�p<#׌<#��<#�
<#�h<#׃<#�<#�<#��<#�u<#�<#�p<#�m<#�
<#�|<#ؙ<#�<#�w<#�l<#��<#�
<#�w<#��<#�<#�j<#�Q<#�-<#ב<#گ<#�~<#�<#׹<#�<#�Y<#�]<#ׁ<#�<$Am<#�7<#�E<$<$)�<$<$�<#�<#�<#�&<$�<$�<$�o<&�<$�P<#�0<#��<$<$�x<#�z<%�&<'�A<&��<%��<#�<#ٯ<#�m<$<#�<#�<#�S<$f�<#ה<#�<#�s<&��<)��<#�<#�s<$�<#�e<$"@<#�<$�<+<,�&<1�3<'W�<$�#<&�Y<(oy<)	�<#�.<(d�<'�<#�.<#�f<#�,<$a�</��<$n	<$�(<'hP<$a9<$w!<#�|<$Ip<#��<#�<#ر<#�<#ڐ<#�<#ץ<#��<#؝<#�x<#�<#�g<#��<#��<#��<#�<#�z<#��<#׮<#׮<#�<#��<#��<$Z�<#��<#�V<#�=<#�c<$'�<$�<#��<$�<#�a<#�g<#�`<#�<#�<#��<#�<#��<#�<$�<%<�<#�X<#��<$|�<#��<#ؑ<#ؖ<#�<#�x<#��<#�&<#�!<%�(<$7�<$�<#�<$�&<$B<#��<#ם<#�<#�~<#ۣ<#�d<#��<#�<$�<#��<#�,<$��<#�=<#�|<#�j<#ݿ<#�T<#ڗ<#�<#�B<#��<#�Q<#ت<#�X<#�<#�#<#�<#�<#�X<#ׂ<#��<#��<#�<#ۂ<#�<#׀<#�<#��<#��<#�C<#�<#��<#�<#׏<#�<#�<#�q<#ڸ<#��<#ؙ<#�
<#ت<#ڞ<#�<#�<#؊<#�+<#�t<#��<#�<$'�<$3<#��<#�Q<#٨<#��<#��<#�<#�d<#�m<#��<#�X<#�_<#�<#�<#ި<#ۚ<#ׯ<#ؕ<#ۚ<$%w<#�<#�<#׎<#�L<#ף<#׶<#�y<#�<#�`<#׸<#��<$�<#�t<#�<#�<#�(<#�4<#�<#ח<#�\<#�_<#�<#�\<#�<#�g<#؈<#�<#�^<#�<#�<#ב<#�@<#��<#ׄ<#�<#�&<#�<#�<#�<#�W<#�<#��<#��<#׿<#�<#ת<#��<#׬<#ٓ<#��<#ג<#�<#�<#׌<#�c<#��<$y<#۹<#׬<#��<#�
<#�<#�<#�^<#�<#�^<#�f<#�m<#؞<#�d<#�<#�<#�<#�`<#د<#��<#�<#ש<#�<#��<#��<#�)<#�9<#�<#�<#٢<#�!<#�T<#�<#��<#��<#׮<#�D<#��<#�<#؈<#�<#ף<#׬<#ך<#�M<#י<#כ<#�D<#۠<#�<$�<#�.<$�<#�<#׻<#�<#�|<#��<$�<#�<#�<#�<#�<#��<#�E<#�
<#� <#�E<#�<#�{<#�%<#�H<#�E<$�<#�<#۠<#�}<#ݵ<#�<#�<#�<#�<#�J<#��<#׭<#ؔ<#�Z<#�9<#�<#��<#��<#כ<#��<#��<#�n<#�3<#׫<#�<#�<#�<#צ<#׮<#כ<#�P<#�<#�<#�<#�Q<#׉<#�W<#�d<#�\<#�
<#�O<$A<#�|<$�<$I<#�<#�	<#�W<#�%<#�[<#�<#�\<#�<#�P<#�O<#�_<#�<#�*<#�;<#��<#��<$+5<$@|<#�<#�n<#�<<#�<#�H<#�<#ׯ<#ד<#�i<#�<#�<#ش<#�[<#ױ<#ۢ<#�O<#�2<#�.<#�|<#�<#�<#�<<#��<#�<#�P<#��<#ۖ<#�<#�[<#�	<#�<#؄<#�h<#�O<#��<#�v<#ذ<#�c<#ۚ<#�<#�=<#�D<#�c<#�D<#�j<#؅<#�<#�R<#׍<#�M<#�<#ۭ<#׬<#�H<#ז<#�<#�`<#ד<#�9<#�A<#�bG�O�<#�<#��<#ؑ<#�<#�O<#ׄ<#�G<#ۘ<#�<#�<<#�E<#�b<#�b<#�i<#ؑ<#�n<#؁<#�<#�`<#׶<#�%<#ۄ<#�P<#ؑ<#�n<#؟<#�g<#�[<#�<#�_<#�[<#�
<#�h<#�
<#؏<#�<#�<#�u<#�<#ؒ<#�p<#؝<#؏<#�<#�B<#�H<#װ<#�c<#ؙ<#�<#�0<#׀<#�a<#�a<#�a<#�a<#�a<#�a<#�a<#�a<#�a<#�`<#�Z<#�<#�Z<#�Z<#�<#�i<#�<#�<#�a<#�g<#�|<#�e<#�a<#�h<#؎<#�g<#�Z<#�<#؇<#�o<#��<#�~;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.07                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXSCUTSCUTSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929142749201909291427492021022313492720210223134927202102231349272021022314543820210224114355  CV  CV  QCF$QCF$QCF$IP  IP                                  PSAL            TEMP            TEMP                                            G�O�G�O�Gc� Gc� Gc� G�O�G�O�G�O�G�O�Gc� Gc� Gc� G�O�G�O�G��G��G��G��G��G��G��                                131072          131072          131072                                          