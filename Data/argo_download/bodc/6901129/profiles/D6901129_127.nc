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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  J�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  V8   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  a�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  d�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  jp   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  u�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �L   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �d   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �4   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �P   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �l   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ˈ   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �H   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ͬ   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �8   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �T   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �p   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ΌArgo profile    3.1 1.2 19500101000000  20210225040439  20210225040439  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               A   BO  125415                          2C  D   APEX                            6229                            120210                          846 @ׄ��#�1   @ׄ��#�@Q���l��7���vȴ1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @33@�33@�  A   A   A@  A`  A���A�  A�  A���A���A���A���A�  A�33B��B��B  B ffB(  B0  B8  B@  BHffBP  BW33B_��BhffBpffBx  B�  B�  B�33B�  B�  B�  B���B���B�  B�33B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  C   C  C  C�C  C
  C  C  C  C  C  C�fC  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D���D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D���D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�<�D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D���D��3B
p�B
p�B
p�B
q�B
q�B
q�B
q�B
q�B
q�B
q�B
p�B
p�B
q�B
q�B
q�B
q�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
q�B
t�B
r�B
q�B
q�B
r�B
s�B
s�B
t�B
v�B
w�B
w�B
u�B
t�B
s�B
w�B
w�B
z�B
z�B
y�B
x�B
z�B
z�B
{�B
|�B
}�B
}�B
�B
�7B
�1B
�7B
�DB
�DB
�PB
�VB
�bB
�oB
�uB
��B
��B
��B
��B
��B
�LB
�wB
�)B
�fB
�yB
��BhB(�BD�BS�B^5B`BB`BBaHBaHBdZBiyBn�Bq�Bt�Bu�Bu�Bu�Bv�Bv�Bv�Bw�Bx�By�By�By�By�By�Bz�B|�B~�B�B�B�B�B�B�B�%B�%B�%B�%B�+B�+B�+B�+B�1B�1B�1B�1B�1B�+B�+B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B� B�B� B}�B~�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�1B�=B�VB�bB�bB�hB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�'B�'B�'B�'B�'B�-B�-B�3B�?B�LB�RB�XB�^B�dB�^B�dB�jB�qB�wB��B��B��BBBBB��B��B��B��B�wB��B��B��B��B�}B�}B�}B�}B�wB�wB�qB�qB�jB�jB�dB�^B�XB�XB�RB�RB�RB�RB�RB�RB�LB�LB�LB�LB�FB�FB�FB�?B�?B�3B�-B�-B�-B�-B�'B�'B�'B�'B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>B�\>A�7>A�7>B�\>B�\>A�7>A�7>A�7>B�\>C��>C��>C��>A�7>Xb>P�`>F��>F��>E��>M��>T��>Xb>e`B>gl�>gl�>fff>\(�>R�>hr�>fff>s�F>}�>w��>m�h>w��>y�#>�%>��>�o>�%>���>�5?>�"�>�5?>��
>���>�l�>���>��h>�?}>��j>�p�>���>��;>��>�
=?��?
��?.V?:��?<(�?Rn�?e�T?�A�?�v�?��T?�ƨ?��?���?�V?�/?��?��?��T?���?�"�?��?�dZ?��?���?��m?�1?�1?��?���?��?��?��?��?��h?�V?���?�&�?��?�-?���?�9X?�z�?��/?���?��j?��j?��/?��?�?}?�`B?�?�?��T?��T?�`B?�?}?���?��!?�G�?� �?�A�?��7?�\)?�V?���?��9?�r�?�~�?��?���?���?��!?��\?���?���?���?�A�?��`?��7?�S�?�o?��7?�Ĝ?��;?��R?��?��?�;d?� �?��;?��;?��w?�A�?���?�hs?�o?��
?��/?���?�dZ?��h?���?���?��?�%?��?�9X?�?�l�?��P?�b?�K�?�
=?�ȴ?�E�?�ff?�ff?��+?���?���?��+?��+?���?�ȴ?�x�?�I�?�V?��?��w?��w?��?��`?���?°!?\?�33?Ł?Ǯ?�Q�?�ff?���?ļj?ļj?ļj?���?�?�
=?�
=?�K�?�
=?�?���?�S�?�J?�M�?�33?��?��`?���?�%?�Ĝ?� �?� �?� �?�;d?�{?��m?�7L?��y?���?�z�?��?���?�hs?��;?�v�?�{?��?�"�?�x�?�1'?�l�?�ff?�t�?���?���?��u?�K�?�ȴ?�E�?�?�?��?�?}?��?���?��!?�n�?�J?��`?��w?���?�b?���?��F?��F?��?��7?��?�A�?�w?;d?~5??z^5?t�j?n�?k�?j=q?i7L?g+?e��?co?bJ?_�w?\�?T�j?R-?PbN?O��?O\)?N�?N{?M��?MO�?L�D?KC�?J��?I7L?G�?G�?F�y?E�T?E��?B��?A��?>v�?=�?=/?<�?<j?:�H?9X?6E�?3��?0�`?/\)?.{?+?(�9?(r�?(1'?'l�?&$�?"�\?"M�?"J?"J?"J?!��?!�7?   ?�?��?�+?E�??��?�j?z�?��?t�?n�?�?&�?�`?bN? �?V?I�?	��?�?`B?Z?o?�\?�\?��?o?��?�\?M�?��?%? Ĝ?   >�|�>�v�>�p�>��H>�ȴ>�?}>��>�V>���>�ff>�`B>��
>���>�M�>���>���>�G�>�;d>޸R>ܬ>�(�>��>և+>���>���>�t�>�n�>��`>��`>��;>�O�>�I�>�C�>���>ȴ9>��>ě�>�o>�o>�o>\>�%>�|�>���>�v�>��>��>�dZ>�^5>�^5>���>�ȴ>�>��F>�33>���>���>� �>��D>�1>�1>��>�~�>�~�>���>�x�>�r�>��y>�S�>�M�>�G�>�Ĝ>�A�>��w>��w>�;d>��R>���>��>��u>�
=>��+>�>���>���>�>�>�t�>��`>�\)>��>���>�=q>��^>��^>���>���>�1'>���>�J>{�m>y�#>w��>n��>j~�>gl�>cS�>e`B>dZ>`A�>Z�>Y�>W
=>W
=>V>V>T��>R�>Q�>O�;>L��>G�>>v�>8Q�>6E�>49X>2->2->0 �>0 �>/�>'�>&�y>$�/>#�
>#�
>�R>�R> Ĝ>�R>��>�u>z�>�+>�+>�>�P>��>z�>t�>V>I�>V>V>n�>\)>o>   =�x�=�`B=�S�=�G�=�/=�
==���=\=�j=�Q�=�E�=� �=�1=��T=��
=��w=��w=��-=���=��P=��P=�t�=�t�=�hs=�O�=�C�=��=}�=q��=m�h=ix�=H�9=H�9=#�
=t�=C�<�h<�j<�9X<��
<�o<49X<D��<#�
<t�;�`B:�o��`B�u�u��o��t����㼴9X��9X���ͼ��#�
�,1�49X�P�`�ixսixս�%���w��1����Ƨ�Ƨ�ȴ9�������ͽ����F�   �C��V�bN�z����$�/�2-�=p��C���F��L�;N��N��P�`�S�ϾT���W
=�`A��gl��o���w�پy�#��  ����������������˾��˾��˾��𾇮���9���9���9��7L��=q��V������������������(���(������5?��;d���w���徦ff���y���þ���1���h��V������33���j����E����پ��پ�X���m��p����۾��7��  �����J��o��o�š˾�1'�������;�V��\)��\)���`�և+�������ڟ��ڟ��ۥ��/�߾w��A���Ĝ��Ĝ��G���������������MӾ�������S����
���
���
��Z��Z��Z��Z���/���T��l����xվ�1��1��{��33���j��E���ȴ��KǾ��#��푾�vɾ�|� A�� Ĝ� Ĝ�%����+�	7L�	x�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @\)@�G�@�{A
=A#
=AC
=Ac
=A�Q�A��A��A�Q�A�Q�A�Q�A�Q�A�B \)B\)B\)BB!(�B(B0B8B@BI(�BPBW��B`\)Bi(�Bq(�BxB�aHB�aHB��{B�aHB�aHB�aHB�.B�.B�aHB��{B�aHB�aHB�aHB�aHB�.B�aHB�aHB�aHB�aHB�aHB�aHB�aHB�aHB�aHB�.B�aHB�aHB�aHB�aHB�.B�aHB�aHC 0�C0�C0�CJ>C0�C
0�C0�C0�C0�C0�C0�C
C0�C0�C0�C0�C 0�C"0�C$0�C&0�C(0�C*0�C,0�C.0�C00�C20�C40�C60�C80�C:0�C<0�C>0�C@0�CB0�CD0�CF0�CH0�CJ0�CL0�CN0�CP0�CR0�CT0�CV0�CX0�CZ0�C\0�C^0�C`0�Cb0�Cd0�Cf0�Ch0�Cj0�Cl0�Cn0�Cp0�Cr0�Ct0�Cv0�Cx0�Cz0�C|0�C~0�C�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC��C�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RC�RD )D �)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D	)D	�)D
)D
�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D)D�)D )D �)D!)D!�)D")D"�)D#)D#�)D$)D$�)D%)D%�)D&)D&�)D')D'�)D()D(�)D))D)�)D*)D*�)D+)D+�)D,)D,�)D-)D-�)D.)D.�)D/)D/�)D0)D0�)D1)D1�)D2)D2�)D3)D3�)D4)D4�)D5)D5�)D6)D6�)D7)D7�)D8)D8�)D9)D9�)D:)D:�)D;)D;�)D<)D<�)D=)D=�)D>)D>�)D?)D?�)D@)D@�)DA)DA�)DB)DB�)DC)DC�)DD)DD�)DE)DE�)DF)DF�)DG)DG�)DH)DH�)DI)DI�)DJ)DJ�)DK)DK�)DL)DL�)DM)DM�)DN)DN�)DO)DO�)DP)DP�)DQ)DQ�)DR)DR�)DS)DS�)DT)DT�)DU)DU�)DV)DV�)DW)DW�)DX)DX�)DY)DY�)DZ)DZ�)D[)D[�)D\)D\�)D])D]�)D^)D^�)D_)D_�)D`)D`�)Da)Da�)Db)Db�)Dc)Dc�)Dd)Dd�)De)De�)Df)Df�)Dg)Dg�)Dh)Dh�)Di)Di�)Dj)Dj�)Dk)Dk�)Dl)Dl�)Dm)Dm�)Dn)Dn�)Do)Do�)Dp)Dp�)Dq)Dq�)Dr)Dr�)Ds)Ds�)Dt)Dt�)Du)Du�)Dv)Dv�)Dw)Dw�)Dx)Dx�)Dy)Dy�)Dz)Dz�)D{)D{�)D|)D|�)D})D}�)D~)D~�)D)D�)D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D��D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D��D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�B�D��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD��D��D�D�FD���D��GB
p�B
p�B
p�B
q�B
q�B
q�B
q�B
q�B
q�B
q�B
p�B
p�B
q�B
q�B
q�B
q�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
p�B
t�B
s#B
q�B
q�B
rSB
sbB
s�B
t"B
v�B
w�B
w�B
vCB
u!B
r�B
w�B
w7B
zeB
{'B
zVB
xkB
z�B
z�B
{�B
|�B
~B
}B
�XB
�bB
��B
��B
�gB
��B
�B
��B
��B
�mB
��B
�B
�oB
��B
��B
�!B
�B
�EB
ىB
��B
�CB
�B�B#�BA�BQ�B]�B`DB`-Ba:B`�Bb�Bg�Bm{Bp�Bt�Bu�Bu�Bu�Bv�Bv�Bv�Bw�Bx�By�By�By�By�By�Bz�B|tB~iB��B��B��B��B�B��B�B�<B�%B�B�B�B�B�B�.B�&B�3B�aB�BB��B��B��B��B�B��B��B��B��B��B�*B�eB�&B��B�vB��B�B��B�|B~HB~!B��B��B�sB�.B��B�_B�gB�yB�B��B��B��B�"B�B�B��B��B��B�|B��B��B��B�'B��B�=B�B��B�?B��B�B��B�B��B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�/B��B��B��B��B�HB��B�B��B�3B�OB��B��B��B�BB�+B�*B�B��B��B�EB�<B�tB��B�B��B��B�RB�$B��B��B��B�hB¤B��BBB��B��B�VB��B�ZB��B��B��B��B�B�B�B��B��B�B�B��B��B��B�B�B�.B��B��B��B��B�kB�kB�hB�gB�\B��B��B�eB�nB��B��B��B�QB�sB��B�3B�uB��B��B�CB�LB�?B�XB��B�)B�/B��B�RB�EB�oB�[B��B�BB�zB��B�xB��B�\B�.B�B�!B�,B�B�!B�.B�EB�$B�QB�EB�B�+B�5B�B��B�8B��B�B� B�	B�B�BB�CB��B�gB��B�=B�2B��B�^B� B��B�B�2B��B�B��B��B��B��B��B�=B��B��B�_B� B�
B�B��B��B�	B�
B� B�
B�B��B�B� B�FB�SB�ZB�\B�XB�!B�-B�	B��B��B��B��B��B��B�	B�B��B�B��B�B�B�(B�IB�B�AB�XB�LB�'B�B�B�B��B��B��B��B�B��B�B��B�B�7B�	B��B��B��B�
B��B�B�"B� B��B��B�B�7B��B�B��B��B��B�B�B��B��B��B��B�B��B��B�B��B��B�B��B�B��B��B�5B��B��B��B��B��B��B��B��B�B�1B��B��B��B��B��B��B��B��B�,B� B��B�B��B��B��B��B��B��B�(B�'B�B��B�B�&B��B��B��B��B�!B�B�QB�EB��B��B�NB�B�	B�B��B��B�B�,B��B��B��B��B��B��B�B��B��B�	B�#B�SB�2B�B�B�B��B�B��B��B�<B��B�B��B��B�%B��B��B�B�B�B�B��B��B��B��B��B�#B��B�&B�B��B��B��B�B�xB�B�lB�B��B��B�B�B�3B�B�B�B��B�B�B�B��B�B��B��B��B�B��B� B��B��B�B��B�B�B�B��B��B�GB��B�VB�B�B�'B�2B��B� B�B� B��B��B��B��B�/B�GB�DB��B��B��B��B�B��B�
B�!B�XB��B��B�3B�&B��B�*B��B�+B�UB�B��B��B��B��B�$B�wB�'B�[B��B��B�B�B�NB�rB�YB�B��B�B��B��B��B��B��B��B�=B�'B�/B�-B��B�B�,B��B��B��B��B��B��B��B��B��B��B��B��B��B�/B�hB��B�B��B��B��B��B��B��B��B��B�B� B��B��B��B��B��B��B�8B��B��B��B��B��B��B��B�	B��B��B�B��B��B��B��B��B�B�	B�
B��B��B��B��B��B�JB�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�@B��B��B��B��B�B�B��B��B��B��B��B��B�rB�EB�B��B��>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>A�7>B�\>A�7>A�7>B�\>B�\>A�7>A�7>A�7>B�\>C��>C��>C��>A�7>Xb>P�`>F��>F��>E��>M��>T��>Xb>e`B>gl�>gl�>fff>\(�>R�>hr�>fff>s�F>}�>w��>m�h>w��>y�#>�%>��>�o>�%>���>�5?>�"�>�5?>��
>���>�l�>���>��h>�?}>��j>�p�>���>��;>��>�
=?��?
��?.V?:��?<(�?Rn�?e�T?�A�?�v�?��T?�ƨ?��?���?�V?�/?��?��?��T?���?�"�?��?�dZ?��?���?��m?�1?�1?��?���?��?��?��?��?��h?�V?���?�&�?��?�-?���?�9X?�z�?��/?���?��j?��j?��/?��?�?}?�`B?�?�?��T?��T?�`B?�?}?���?��!?�G�?� �?�A�?��7?�\)?�V?���?��9?�r�?�~�?��?���?���?��!?��\?���?���?���?�A�?��`?��7?�S�?�o?��7?�Ĝ?��;?��R?��?��?�;d?� �?��;?��;?��w?�A�?���?�hs?�o?��
?��/?���?�dZ?��h?���?���?��?�%?��?�9X?�?�l�?��P?�b?�K�?�
=?�ȴ?�E�?�ff?�ff?��+?���?���?��+?��+?���?�ȴ?�x�?�I�?�V?��?��w?��w?��?��`?���?°!?\?�33?Ł?Ǯ?�Q�?�ff?���?ļj?ļj?ļj?���?�?�
=?�
=?�K�?�
=?�?���?�S�?�J?�M�?�33?��?��`?���?�%?�Ĝ?� �?� �?� �?�;d?�{?��m?�7L?��y?���?�z�?��?���?�hs?��;?�v�?�{?��?�"�?�x�?�1'?�l�?�ff?�t�?���?���?��u?�K�?�ȴ?�E�?�?�?��?�?}?��?���?��!?�n�?�J?��`?��w?���?�b?���?��F?��F?��?��7?��?�A�?�w?;d?~5??z^5?t�j?n�?k�?j=q?i7L?g+?e��?co?bJ?_�w?\�?T�j?R-?PbN?O��?O\)?N�?N{?M��?MO�?L�D?KC�?J��?I7L?G�?G�?F�y?E�T?E��?B��?A��?>v�?=�?=/?<�?<j?:�H?9X?6E�?3��?0�`?/\)?.{?+?(�9?(r�?(1'?'l�?&$�?"�\?"M�?"J?"J?"J?!��?!�7?   ?�?��?�+?E�??��?�j?z�?��?t�?n�?�?&�?�`?bN? �?V?I�?	��?�?`B?Z?o?�\?�\?��?o?��?�\?M�?��?%? Ĝ?   >�|�>�v�>�p�>��H>�ȴ>�?}>��>�V>���>�ff>�`B>��
>���>�M�>���>���>�G�>�;d>޸R>ܬ>�(�>��>և+>���>���>�t�>�n�>��`>��`>��;>�O�>�I�>�C�>���>ȴ9>��>ě�>�o>�o>�o>\>�%>�|�>���>�v�>��>��>�dZ>�^5>�^5>���>�ȴ>�>��F>�33>���>���>� �>��D>�1>�1>��>�~�>�~�>���>�x�>�r�>��y>�S�>�M�>�G�>�Ĝ>�A�>��w>��w>�;d>��R>���>��>��u>�
=>��+>�>���>���>�>�>�t�>��`>�\)>��>���>�=q>��^>��^>���>���>�1'>���>�J>{�m>y�#>w��>n��>j~�>gl�>cS�>e`B>dZ>`A�>Z�>Y�>W
=>W
=>V>V>T��>R�>Q�>O�;>L��>G�>>v�>8Q�>6E�>49X>2->2->0 �>0 �>/�>'�>&�y>$�/>#�
>#�
>�R>�R> Ĝ>�R>��>�u>z�>�+>�+>�>�P>��>z�>t�>V>I�>V>V>n�>\)>o>   =�x�=�`B=�S�=�G�=�/=�
==���=\=�j=�Q�=�E�=� �=�1=��T=��
=��w=��w=��-=���=��P=��P=�t�=�t�=�hs=�O�=�C�=��=}�=q��=m�h=ix�=H�9=H�9=#�
=t�=C�<�h<�j<�9X<��
<�o<49X<D��<#�
<t�;�`B:�o��`B�u�u��o��t����㼴9X��9X���ͼ��#�
�,1�49X�P�`�ixսixս�%���w��1����Ƨ�Ƨ�ȴ9�������ͽ����F�   �C��V�bN�z����$�/�2-�=p��C���F��L�;N��N��P�`�S�ϾT���W
=�`A��gl��o���w�پy�#��  ����������������˾��˾��˾��𾇮���9���9���9��7L��=q��V������������������(���(������5?��;d���w���徦ff���y���þ���1���h��V������33���j����E����پ��پ�X���m��p����۾��7��  �����J��o��o�š˾�1'�������;�V��\)��\)���`�և+�������ڟ��ڟ��ۥ��/�߾w��A���Ĝ��Ĝ��G���������������MӾ�������S����
���
���
��Z��Z��Z��Z���/���T��l����xվ�1��1��{��33���j��E���ȴ��KǾ��#��푾�vɾ�|� A�� Ĝ� Ĝ�%����+�	7L�	x�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#��<#�t<#�j<#��<#��<#��<#��<#��<#��<#�[<#�\<#��<#ܵ<#��<#��<#�o<#�V<#��<#��<#�<#�<#��<#��<#�i<$�4<#��<$�<#�!<#�F<#�Y<#��<#�I<$#<#�<#�}<#�!<$�<$�<$Cu<#�J<$�<#�<#��<$�<#�B<#�<#�4<#�)<#�a<#�<$LA<%�l<#�<#��<#��<#�<#��<#��<#�r<$�<#��<$6�<%Z�<#�p<#ښ<$ZV<3�2<$�<>�6<(��<#�d<0��<.¢<44<7	�<)�N<'-N<#��<#��<#�;<#׬<#�1<%1�<%L�<$�B<$B%<#��<#�<#��<#��<#�<#��<#٪<#��<#��<#�c<#ٷ<#�{<#�<#��<#��<#��<$ <#�<#�<#��<#�'<#�<#�<#ת<#�6<#ٚ<#��<#�<#��<#��<#�-<#�<#��<#��<#��<#�{<$?G<$ �<$)�<$%<#��<#�T<$v�<$�
<$��<$k�<#�h<$�<$�P<$q~<$<$W�<#�R<#��<%��<#��<$C<#�X<#ے<$
�<#�z<$1 <#�A<#�9<$<#�+<#�<#נ<#�<#�<#��<#��<#�[<#�1<#�v<$	�<#�\<#� <#�<&ۖ<$;u<#�E<#�<$�<#��<$a<#��<$�<#�c<#�m<#ׂ<#�B<#��<#�<#�V<#�$<#�l<#��<#��<#ْ<#��<#ى<#��<#�7<$i<${�<#��<$"<#�<#�H<#ۀ<#ו<$3<#��<#�v<#��<$E%<$4�<#��<$Wa<$$�<#�B<#�9<#��<#�<#��<#�z<#�p<#�*<#��<$�<$L�<#��<$�<#�<#�&<$U<$�<#��<#�<#�y<#��<#ٱ<#��<#��<$<$�M<$�i<$�<<$�<$�<#��<$�<$+�<$4�<$&0<#�<$)h<$2�<$D�<$K<#�<$4<$�*<&T<&ō<$/�<$h<#��<#�><#��<#߽<#��<#�N<#�<$<$�<#��<#�@<$	<$�<$i
<'��<%T<$�<#ک<#��<$(�<$A<#�<#��<#ߎ<#�2<$`�<$�<$�z<$B�<#�<#�!<$<#��<$�<#��<$M<$<<%�<$!�<#��<#�<#�}<#��<#��<#܎<#��<#�<#�]<#��<#�<<#�v<#�r<#�w<#�u<#ݰ<$%�<#�<$>�<#�<#��<#�X<#ߢ<#�t<#�Q<$5%<$�<$6_<#�<#�2<$6*<$;<#�!<#�u<#�<#�D<$I�<#ݷ<#�<#ّ<#٤<#�"<#��<#�'<$EV<$` <$�<#��<#��<#��<#�r<#�U<#��<#�<#�<#�?<#��<#�T<#ߜ<#�<#��<$�<$q<$�<$�<#��<#�<#��<#١<#��<#��<#�!<#�8<#�T<#��<#�<#ܕ<#�<#��<#�(<#�E<#�`<$�<#�<$ �<$�<$�<#�?<#�!<#�<#��<#�P<#�<#٢<#�u<#�<#��<#�<#�c<#�4<#��<#�<#߫<#�d<#��<#�<#٘<#�<#�<#�9<#߬<#܌<#�<#�U<#�<#�<#٢<#ٟ<#�N<#�<#��<#�W<#�/<#�<#��<#�<#��<#�<#�<#��<#�<#��<#ܴ<#�<#߹<#�D<$ m<#�<#ٛ<#�4<#�Z<#پ<#�"<#��<#�<#��<#��<#��<#ߟ<#�
<#ۈ<#�<#١<#�<#ܼ<#�g<#�L<#ܱ<#�<#�\<#�D<#�g<#�<#�?<#�q<#�9<#�<#�N<#ܳ<#�<#�)<#�l<#�o<#�<#��<#�<#�<$<$Y<#�<#��<$<#�O<#�p<#�j<#�5<#�3<#��<#�<#��<#�Q<#ٯ<#�	<#�<#�q<#�x<#�<#�8<#�<#�<$<#��<#�!<#߹<#�t<#ٵ<#�V<#ٺ<#ܰ<#��<#�><#�z<#�%<#��<#�<#�<#�&<#�O<#�%<#�g<#�_<#�.<#�l<#��<#�<#�*<#��<#��<#�i<#�<#�
<#ك<#�<#�U<$0<#�V<$$s<#ߏ<#�1<#�=<#�#<#��<#��<#��<#�9<#ߩ<#�T<#��<#��<#�<#�L<#�Z<#ٰ<#� <#�5<#�h<#ٺ<#� <#ٯ<#�)<#�}<#�k<#��<#�:<#�<#�<#��<$<#�w<$7<#�:<#�i<#��<#�,<#�<#��<#�<#�?<#��<#�<#�S<#�+<#�<$�<$w<#�<#�'<#�z<#�e<#�[<#��<#�<#��<$<#�<#�f<#��<#�/<#�y<#��<$_�<#�<$�<#�<#�l<#�U<#�<#ܯ<#��<$B�<#�Q<$%K<#�@<#�l<#�<#�I<$z<$D,<$(�<#��<#�J<#��<#ߠ<#٭<#ߐ<#�<#�a<#�<$m<$<$	<$^<#�<#�/<$S<#ܟ<#��<#ٔ<#��<#ن<#ٛ<#ߐ<#ߪ<#�b<#َ<#ه<#�#<#�r<$	,<$?x<#�M<#�J<#�,<#߿<#��<#پ<#�0<#�<#ߚ<#ܲ<#��<#��<#��<#�<#�V<#�}<#��<#�T<$�<#�*<#�<#ߠ<#�H<#�v<#��<#�<#�<#�]<#�<#��<#�!<#�H<#��<#߇<#ـ<#�<#�<#��<#�<#�"<#�a<#�p<#�U<$%p<#�<#�K<#�<#٪<#ߛ<#�<#�"<#��<#�S<#��<#�r<#�N<#�<#�)<#�<#�t<#�<#ܛ<#��<#�s<#��<#�%<#�x<#�k<#��<#�<#ߠ<#��<#�Z<#�.<#�<#�L<#��<$7<#�2<#��<#�!<#�0<#��<#�I<#�<#ߜ<#ߚ<#�T<#٨<#�.<$R<$!<#�<#�'<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.19                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929144932201909291449322019092914493520190929144941202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@33@33@33G�O�G�O�G�O�G�O�D��3D��3D��3G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          