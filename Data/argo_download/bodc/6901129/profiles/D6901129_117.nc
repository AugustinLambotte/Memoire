CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
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
resolution        ?�������     
L  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  Il   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  S�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ^   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  `�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  c,   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
L  e�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  p   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  zX   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �8   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
L  �`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
L  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �D   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �0   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �L   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �h   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �(   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �4   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �P   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �lArgo profile    3.1 1.2 19500101000000  20210225040831  20210225040831  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               uA   BO  125405                          2C  D   APEX                            6229                            120210                          846 @�l1�i�@1   @�l1�i�@@QS��$��3���
=q1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @�ff@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  BffB  B   B(  B0  B8  B@  BH  BP  BX  B`  BhffBp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�|�D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D���D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�|�D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�C3D�� D�� D��B��B��B��B��B��B�yB	+B	�B	�B	�B	%�B	:^B	aHB	|�B	��B	�B	��B	�B	��B
B
B
$�B
7LB
P�B
]/B
p�B
�=B
��B
��B
��B
�B
��B
�B
�`B
�BBJB�B+B33B<jBF�BO�BYBbNBdZBiyBt�By�B�B�\B��B��B��B��B�B�9B�LB�XB�dB�qBBƨBǮBǮBȴBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȴBĜB��B�}B�wB�qB�^B�XB�XB�RB�LB�FB�?B�9B�9B�9B�9B�FB�LB�FB�?B�9B�-B�'B�!B�B�B�!B�!B�!B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�!B�!B�'B�-B�-B�'B�-B�-B�-B�-B�-B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�uB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�-@?}@j@�F@�@��@�@��@&�@ bN?��9?�h?ٙ�?ݑh?�-?�C�?�\)?�  ?��u?��#?���?;��?J?��?j?,�D?<�?A��?D�?D�/?E�?J��?Y�?_|�?g�?}p�?��?��H?�b?�\)?�?�O�?�9X?�C�?�%?�M�?���?��m?θR?�ff?��?�o?��?���?��H?��;?���?�`B?�ff?��?�r�?�^5?�I�?��?��?�O�?��h?��?��R?���?��w?���?��?���?�V?��?�j?��D?��D?��D?�j?�(�?�?��j?��?�?� �?��?׮?ա�?У�?̋D?˅?ʟ�?��?�
=?�M�?�&�?���?�{?���?�|�?���?�p�?��H?��u?�E�?��!?���?�~�?�+?�7L?�x�?���?��?���?�{?�C�?�$�?��
?�S�?�  ?���?�C�?���?��^?��?�$�?��j?��?���?�%?}�-?z��?u�?m�h?l�D?p�`?wK�?z�?��?|�?{dZ?}�-?{��?s��?l�D?a�7?`Ĝ?f��?m�h?n{?pbN?q��?r�?st�?st�?st�?s33?s33?s33?r�!?r�?o��?l1?`A�?N��?KC�?Hr�?DZ?AG�?>��?;��?:�H?8�u?:�?>��?=�-?9��?4z�?2n�?/\)?+ƨ?*=q?'�?$��?!�7? �?   ? A�?   ?p�?�?�?�?p�?v�?v�?5??/?dZ?�?E�?��?��? �?�?	7L?�y?�?�\? �>�j>��m>�F>�~�>�r�>�r�>�>>�&�>�->�&�>� �>��>�&�>� �>��/>޸R>ؓu>�z�>ڟ�>�M�>�G�>�5?>��T>�V>>��y>�"�>�>���>�t�>�hs>�bN>��;>���>���>�=q>�J>�p�>�dZ>���>�`B>��+>�hs>�bN>��>�\)>�V>�C�>���>�  >t�j>r�!>p��>_;d>_;d>\(�>O�;>=p�>%�T>"��>$�/>/�>(��>$�/> Ĝ>!��>�u>�>hs=��=���>o>+>C�>1'>
=q>o=�`B>%>�=�F=�l�=��=��P=@�=C�<�/<�9X<D��;D��;o;D��;�`B<T��<u<�j<�=�P=t�=o=o<�/<D��<T��<�t�<�9X<�h='�=49X=8Q�=@�='�<��<�t�;�`B;�`B<���<��<�t�;ě�;�`B<#�
<�o=o=e`B=ix�=�%=�C�=�7L=�o=�o=��=�+=�O�=�hs=�C�=��=��=�7L=�+=y�#=y�#=}�=e`B=Y�=T��=P�`=<j=�w=t�=C�=+<�<�h<�`B<���<ě�<�9X<��
<�C�<e`B<#�
<t�<t�<49X<o:�o�o��o�T����o��C���C���t���9X��`B��P�,1�0 Ž49X�8Q�@��L�ͽ]/�]/�]/�m�h��%�����+��7L��O߽������㽓t���hs��hs���P���㽟�w���
��1��j�\�ȴ9���ͽ�����������
=��;d��S���`B����F��h���پ��1'�
=q�O߾n�����u��w�&�y�)��-V�.{�0 ž0 ž0 ž333�49X�333�49X�5?}�8Q�:^5�:^5�=p��>vɾ?|�B�\�F��G��J���L�;M��O�;�R�V�Xb�Xb�Z��\(��_;d�`A��aG��bMӾe`B�ixվj~��l�D�m�h�o���s�F�t�j�u�vȴ�w�پz�H������\��o��o��������������𾉺^��ƨ��O߾�����;���;��bN��hs��hs���������+��
=���u���㾝�-���w��Ĝ�������
���y���r���r����þ��þ��羪~���1��{�� ž�-���!���!��33���F��9X��?}��ȴ�������H��푾�p���󶾽󶾾vɾ��۾�|��  ������7��o�ě�����š˾�$ݾƧ��+�Ǯ��1'�ȴ9��=q��I���I����;�O߾�V��hs�����Ͼ�z��z�����׍P����ݲ-��A���MӾ��
���T���þ���h������׾�33��F��F��9X��?}��?}����ȴ���پ�X��p���|�   � A�� �� A�� ���7�o�o�o�S��S��Z���`B�$ݿr���ÿ	7L�	�^�	�^�	�^�	��	���1�I���ͿO߿������� ſ ſ�׿�`�-��9X��+�b�X�X�X�������������������������������������������u��u��u��u��u�Q�Q�b���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�G�@��H@��HAp�A=p�A]p�A}p�A��RA��RA��RA��RAθRA޸RA�RA��RB\)BB\)B\)B'\)B/\)B7\)B?\)BG\)BO\)BW\)B_\)BgBo\)Bw\)B\)B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BîBǮBˮBϮBӮB׮BۮB߮B�B�B�B�B�B��B��B��C�
C�
C�
C�
C	�
C�
C�
C�
C�
C�
C�
C�
C�
C�
C�
C�
C!�
C#�
C%�
C'�
C)�
C+�
C-�
C/�
C1�
C3�
C5�
C7�
C9�
C;�
C=�
C?�
CA�
CC�
CE�
CG�
CI�
CK�
CM�
CO�
CQ�
CS�
CU�
CW�
CY�
C[�
C]�
C_�
Ca�
Cc�
Ce�
Cg�
Ci�
Ck�
Cm�
Co�
Cq�
Cs�
Cu�
Cw�
Cy�
C{�
C}�
C�
C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D u�D ��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��D	u�D	��D
u�D
��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��Du�D��D u�D ��D!u�D!��D"u�D"��D#u�D#��D$u�D$��D%u�D%��D&u�D&��D'u�D'��D(u�D(��D)u�D)��D*u�D*��D+u�D+��D,u�D,��D-u�D-��D.u�D.��D/u�D/��D0u�D0��D1u�D1��D2u�D2��D3u�D3��D4u�D4��D5u�D5��D6u�D6��D7u�D7��D8u�D8��D9u�D9��D:u�D:��D;u�D;��D<u�D<��D=u�D=��D>u�D>��D?u�D?��D@u�D@��DAu�DA��DBu�DB��DCu�DC��DDu�DD��DEu�DE��DFu�DF��DGu�DG��DHu�DH��DIu�DI��DJu�DJ��DKu�DK��DLu�DL��DMu�DM��DNu�DN��DOu�DO��DPu�DP��DQu�DQ��DRu�DR��DSu�DS��DTu�DT��DUu�DU��DVu�DV��DWu�DW��DXu�DX��DYu�DY��DZu�DZ��D[u�D[��D\u�D\��D]u�D]��D^u�D^��D_u�D_��D`u�D`��Dau�Da��Dbu�Db��Dcu�Dc��Ddu�Dd��Deu�De��Dfu�Df��Dgu�Dg��Dhu�Dh��Diu�Di��Dju�Dj��Dku�Dk��Dlu�Dl��Dmu�Dm��Dnu�Dn��Dou�Do��Dpu�Dp��Dqu�Dq��Dru�Dr��Dsu�Ds��Dtu�Dt��Duu�Du��Dvu�Dv��Dwu�Dw��Dxu�Dx��Dyu�Dy��Dzu�Dz��D{u�D{��D|u�D|��D}u�D}��D~u�D~��Du�D��D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�w�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�w�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�:�D�z�D���D���D�>D�z�D���D�{B�IBσBπBнB�WB��B	
fB	�B	=B	�B	*%B	A�B	_�B	v$B	�B	��B	πB	�B
fB
�B
�B
/oB
3B
P�B
ZDB
m�B
�+B
�B
��B
��B
�B
��B
��B
�B
��B
��B	�B�B(B0�B9�BDBM>BV�Ba�BcOBf�Bs�BwB��B�B��B�tB��B�B��B��B��B��B�B��B��B�wBǔBǋBțBɔB�zB�zB��B�.B�B�WB�AB�B�B��B��B��B�B�B�pB�:BϭBʼB�hBÿB�pB�GB�@B��B��B��B��B�B�B��B��B��B�MB��B��B��B�9B�B�B��B��B��B�GB�jB�B�%B��B��B�nB�1B��B��B�7B�/B��B�B�/B�JB��B��B�yB��B�bB�*B��B�zB��B�HB�B�B��B�gB��B�jB��B��B�}B��B�xB��B�'B��B��B��B��B��B��B�B�*B�+B�7B�-B�.B�DB�(B��B��B�PB�>B��B�uB��B�vB�NB�vB�B�AB��B�B�4B��B��B�NB�sB��B�$B�PB�fB�dB�	B��B��B��B�OB��B��B��B��B��B��B��B�
B�.B�EB�ZB��B��B��B�B��B�2B�B�7B� B�%B��B�yB��B��B��B�iB�nB��B��B��B��B��B��B��B��B�eB�]B�&B�4B� B��B�B�'B�9B��B��B��B�XB�B��B�B��B��B�B��B�B��B�DB�
B��B��B�B�6B��B��B��B��B��B�B�HB�(B��B��B�dB��B��B�-B�uB��B��B�oB�B��B��B��B��B� B��B��B�lB��B�BB�hB�oB��B��B��B�[B�B�B�B��B�8B��B��B�'B��B��B��B��B��B�jB�FB�5B�eB�$B�5B�<B��B��B��B��B�BB��B�TB�bB�AB�B��B��B��B��B�%B�GB�)B��B�B�B�CB�;B��B��B�]B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B�B��B��B�+B�B��B��B�!B�:B�B�B��B�B��B��B�B��B�B�B�B�B�B��B��B��B�B�FB�B�B�\B�B��B��B��B�B�6B�XB�(B��B��B��B�B�B�B��B��B�B�#B�B��B��B�B�1B��B��B��B��B�B�B�B�B�)B�XB�B�B�B��B��B�&B�B�'B�B�B�TB��B��B�7B�`B�B�
B�B�.B�B�B�EB�@B�B�B��B�B��B��B�B��B��B��B��B�B�	B��B�B��B��B�B�!B��B�B�	B��B�	B�B�B�B��B�	B�
B�B��B��B��B�B�B��B�B��B�
B�!B��B��B��B��B�B�8B�B��B��B��B�B��B�B�3B�B�B�B�B��B��B�B��B�=B�B��B��B�B�3B�B�B�B�B�B�2B�B��B��B��B��B�B��B�B�B�B�B��B��B��B��B��B��B�	B�B�B�B��B��B��B��B��B��B��B��B��B�	B�B��B��B��B��B��B��B��B��B�	B�B��B��B��B��B�-B�
B��B��B��B��B�!B�#B�:B�"B�B�
B�B�.B�"B�B��B��B�
B� B��B��B��B��B��B��B��B��B�B�@B�B��B��B��B��B��B�B�.B��B��B��B��B�B�B��B�B�RB�B��B��B��B��B��B��B�B�B��B�B�B��B�B�B�B��B�B��B�$B�B�%B�UB�3B�&B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�-@?}@j@�F@�@��@�@��@&�@ bN?��9?�h?ٙ�?ݑh?�-?�C�?�\)?�  ?��u?��#?���?;��?J?��?j?,�D?<�?A��?D�?D�/?E�?J��?Y�?_|�?g�?}p�?��?��H?�b?�\)?�?�O�?�9X?�C�?�%?�M�?���?��m?θR?�ff?��?�o?��?���?��H?��;?���?�`B?�ff?��?�r�?�^5?�I�?��?��?�O�?��h?��?��R?���?��w?���?��?���?�V?��?�j?��D?��D?��D?�j?�(�?�?��j?��?�?� �?��?׮?ա�?У�?̋D?˅?ʟ�?��?�
=?�M�?�&�?���?�{?���?�|�?���?�p�?��H?��u?�E�?��!?���?�~�?�+?�7L?�x�?���?��?���?�{?�C�?�$�?��
?�S�?�  ?���?�C�?���?��^?��?�$�?��j?��?���?�%?}�-?z��?u�?m�h?l�D?p�`?wK�?z�?��?|�?{dZ?}�-?{��?s��?l�D?a�7?`Ĝ?f��?m�h?n{?pbN?q��?r�?st�?st�?st�?s33?s33?s33?r�!?r�?o��?l1?`A�?N��?KC�?Hr�?DZ?AG�?>��?;��?:�H?8�u?:�?>��?=�-?9��?4z�?2n�?/\)?+ƨ?*=q?'�?$��?!�7? �?   ? A�?   ?p�?�?�?�?p�?v�?v�?5??/?dZ?�?E�?��?��? �?�?	7L?�y?�?�\? �>�j>��m>�F>�~�>�r�>�r�>�>>�&�>�->�&�>� �>��>�&�>� �>��/>޸R>ؓu>�z�>ڟ�>�M�>�G�>�5?>��T>�V>>��y>�"�>�>���>�t�>�hs>�bN>��;>���>���>�=q>�J>�p�>�dZ>���>�`B>��+>�hs>�bN>��>�\)>�V>�C�>���>�  >t�j>r�!>p��>_;d>_;d>\(�>O�;>=p�>%�T>"��>$�/>/�>(��>$�/> Ĝ>!��>�u>�>hs=��=���>o>+>C�>1'>
=q>o=�`B>%>�=�F=�l�=��=��P=@�=C�<�/<�9X<D��;D��;o;D��;�`B<T��<u<�j<�=�P=t�=o=o<�/<D��<T��<�t�<�9X<�h='�=49X=8Q�=@�='�<��<�t�;�`B;�`B<���<��<�t�;ě�;�`B<#�
<�o=o=e`B=ix�=�%=�C�=�7L=�o=�o=��=�+=�O�=�hs=�C�=��=��=�7L=�+=y�#=y�#=}�=e`B=Y�=T��=P�`=<j=�w=t�=C�=+<�<�h<�`B<���<ě�<�9X<��
<�C�<e`B<#�
<t�<t�<49X<o:�o�o��o�T����o��C���C���t���9X��`B��P�,1�0 Ž49X�8Q�@��L�ͽ]/�]/�]/�m�h��%�����+��7L��O߽������㽓t���hs��hs���P���㽟�w���
��1��j�\�ȴ9���ͽ�����������
=��;d��S���`B����F��h���پ��1'�
=q�O߾n�����u��w�&�y�)��-V�.{�0 ž0 ž0 ž333�49X�333�49X�5?}�8Q�:^5�:^5�=p��>vɾ?|�B�\�F��G��J���L�;M��O�;�R�V�Xb�Xb�Z��\(��_;d�`A��aG��bMӾe`B�ixվj~��l�D�m�h�o���s�F�t�j�u�vȴ�w�پz�H������\��o��o��������������𾉺^��ƨ��O߾�����;���;��bN��hs��hs���������+��
=���u���㾝�-���w��Ĝ�������
���y���r���r����þ��þ��羪~���1��{�� ž�-���!���!��33���F��9X��?}��ȴ�������H��푾�p���󶾽󶾾vɾ��۾�|��  ������7��o�ě�����š˾�$ݾƧ��+�Ǯ��1'�ȴ9��=q��I���I����;�O߾�V��hs�����Ͼ�z��z�����׍P����ݲ-��A���MӾ��
���T���þ���h������׾�33��F��F��9X��?}��?}����ȴ���پ�X��p���|�   � A�� �� A�� ���7�o�o�o�S��S��Z���`B�$ݿr���ÿ	7L�	�^�	�^�	�^�	��	���1�I���ͿO߿������� ſ ſ�׿�`�-��9X��+�b�X�X�X�������������������������������������������u��u��u��u��u�Q�Q�b���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<$
�<$$<)�><7��<EX<+y<$|�<$�<*��<0�/<F�H<%��<E��<U�<��<:ڰ<(+P<9�<Dzq<���<jA<2%<#�r<*�J<*�@<$�<$h<#��<#��<$�5<)`�<%7<&O�</�#<'�c<)�<6<*��<)�<)��<)uq<)��<'�+<$?l<$�<)�<%6�<*
<(�;<(�<$�<%06<$o<&�<%ۤ<$3=<$�<$<#��<$PC<$XD<#�[<#��<#��<#ޔ<#��<#�<#�6<#�K<#�A<#�1<#�:<#�*<#ג<#�<#��<#��<#��<#�x<#�<#�w<'�p<)�}<&�P<)s�<'L�<$i<$>,<&�<%g<#�<#�<$�<$:`<%�3<#�N<$�<$@<#�<$0�<#߆<#��<$hM<$N�<$Qg<%R<%u�<%o�<$�c<$J�<#��<$<#��<#�r<'�<$�y<&SD<$L�<#��<$�<%��<#י<#�i<#�<$-�<$z<#��<$�<#�a<#ݑ<$>�<$<$m�<%?�<#ؗ<$}�<%�<$9�<%8�<#��<$'�<#�y<#�6<%=H<%<�<&�N<#��<$�<%KM<#��<$.<#�W<#��<#�n<#�q<#�K<#�H<#��<#��<#�<#ٺ<$6<$'�<'G$<+Z�<$"�<#�<$0�<$	?<#�<$t<#��<#�
<#��<$��<#ܒ<$4�<$c�<#�Z<$�<$H<#�q<#��<$B<$<#��<#�
<#��<#�F<#��<#ץ<#ػ<#� <#�<#�&<#�<#�m<#�<#��<#�<#�=<$I�<#��<#�
<#��<$�B<#�<#�<#��<#�M<#�<#�"<$3<$I<#ٔ<#�;<#�0<#�4<#��<#�Z<#�
<#�<#�^<#ہ<#�5<$�G<$
<$\<#��<$5�<$St<#�o<#��<$Y�<$J�<#�`<$!<$��<#�:<#�k<#�R<#س<#�
<#�V<#�x<#�:<#�:<$1�<#��<#�:<${�<$�W<%"�<#�<#�<#�G<#��<#�<#޽<#��<$<#�}<#�<#�T<$6<#ر<#��<$`<$R<$��<#׾<#�I<$�<#ܤ<#�<#ض<#��<#��<#ט<#��<$G!<#׺<#��<#�.<#��<#�:<#��<#�	<$&�<$;H<#�<#�<#��<$?<$�;<$�<$�<#��<#��<#��<#�D<#ב<#�<#�j<#�<#߷<$�<#�<#�^<#ה<#��<#�f<#��<$�<#ڐ<#�<#��<#�l<$/@<#�<#�}<#�l<#އ<#��<$�<#�L<#��<$.�<$){<$�<$�<#�d<#��<#�4<$n<%<#ݛ<#�><#�+<#מ<#כ<#�<#�9<#�t<#�Z<#�k<#א<#�v<#�i<#޵<#�e<#��<#إ<#�{<#��<#׀<#�n<#�Z<#�[<#��<#׀<#�<#�c<#�u<#�y<#�j<#�q<#�d<#�
<#�<#׋<#כ<#��<#ׂ<#�F<#ݠ<#׏<#�h<#ך<#�<#�<#ק<#ׁ<#��<#נ<#�z<#�E<#��<#�}<#�[<#�x<#�p<#�<#ד<#��<#��<#��<#�R<#��<#�<#�x<#�o<#�<#��<#ג<#�<#��<#غ<#פ<#�<#�<#�<#�><#�q<#׏<#��<#�<#�[<#ؚ<#��<#�<#��<#�<#�D<#��<#ډ<#�C<#��<#�#<#��<#�<#ץ<#�K<#ש<#׿<#�0<#��<#ׅ<#ׅ<#�b<#�
<#��<#ة<#�l<#�{<#�<#ׁ<#�f<#�~<#�
<#؜<#�o<#�f<#�f<#נ<#��<#�Q<#�<#�<#�j<#�<#׏<#׎<#�
<#ؠ<#�
<#�<#׃<#�o<#�t<#�e<#�m<#��<#ׂ<#�
<#�e<#�<#��<#�`<#�t<#�t<#�d<#�<#�?<#إ<#�s<#ؽ<#�s<#�<#�@<#׷<#އ<#�2<#ע<#א<#�
<#ا<#�t<#�
<#�g<#��<#�<#�l<#�c<#װ<#�x<#�F<#��<#�<#�<#�"<#�Q<#�<#�x<#ش<#׃<#ا<#�<#�^<#�]<#ل<#�#<#�<#��<#�N<#�x<#�r<#׽<#�<#מ<#�<#�$<#��<#�\<#�z<#��<#נ<#׭<#�)<#נ<#�V<#�<#ד<#׋<#�q<#�q<#�q<#�q<#�q<#�q<#�q<#�a<#י<#��<#ؒ<#�x<#�h<#�<#ކ<#׮<#�<#�s<#ػ<#�W<#�I<#��<#�"<#ۨ<#�)<#׳<#�C<#�<#ۥ<#�<#�<#�]<#ת<#�<#�Z<#ر<#�i<#�
<#؜<#�o<#�<#�<#��<#�w<#�#<#�o<#ט<#׾<#ٸ<#׈<#�l<#��<#��<#ػ<#�v<#��<#�L<#��<#ׁ<#��<#��<#�<#ׇ<#�<#��<#غ<#�}<#ج<#�B<#�g<#ׄ<#�<#�<#�<#�t<#�m<#�<#؜<#�
<#�{<#�8<#׸<#܆<#�{<#ޠ<#�@<#ت<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#��<#�R<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.16                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929144235201909291442352019092914423920190929144245202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�D��D��D��G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          