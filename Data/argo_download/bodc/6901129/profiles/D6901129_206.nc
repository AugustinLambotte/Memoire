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
_FillValue                 �  Sl   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  U    TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  _L   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  f   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  n�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  p<   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  q�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  x�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  x   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �<   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �(   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �D   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �`   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �    HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �,   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �H   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �dArgo profile    3.1 1.2 19500101000000  20210225043725  20210225043725  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125492                          2C  D   APEX                            6229                            120210                          846 @�Ez&� 1   @�Ez&� @P�C���4�^5?|�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�33@�  @���A   A@  A`  A���A�  A�  A�  A�  A���A�  A�33B   B  B  B  B   B(  B0  B8  B@  BH  BPffBXffB`  Bh  Bo��Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D��D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DVfDV�fDW  DW� DX  DX� DY�B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
B
	7B
R�B
�?B
�/B8RB<jB@�BI�BL�BN�BVBZB`BBffBjBl�Bq�Bt�Bs�Bq�Bq�Br�Br�Bt�Bt�Bs�Bs�Bs�Bt�Bx�B� B�B�B�%B�+B�=B�7B�1B�1B�1B�1B�7B�=B�=B�=B�7B�1B�%B�B�B{�Bz�By�Bw�Bu�Bt�Bq�Bo�Bo�Bp�Br�Bu�Bw�Bw�Bv�Bu�Bt�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Bt�Bu�Bu�Bu�Bu�Bt�Bu�Bv�Bt�Br�Br�Br�Bs�Bt�Bv�Bx�Bz�Bz�Bz�By�Bx�Bx�Bw�Bw�Bv�Bv�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bt�Bt�Bs�Bs�Bs�Br�Br�Br�Br�Bs�Bs�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bx�Bx�Bx�By�Bz�Bz�Bz�B{�B{�B|�B|�B}�B|�B|�B}�B}�B}�B|�B{�B{�B|�B|�B}�B}�B~�B~�B~�B~�B� B�B�B�B�B�B�%B�%B�%B�%B�%B�+B�1B�1B�1B�1B�7B�7B�=B�=B�=B�DB�DB�DB�JB�JB�JB�JB�JB�PB�PB�PB�PB�VB�VB�\B�\B�\B�\B�bB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�oB�oB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?#�
?$Z?%�T?$�/?$Z?$Z?#��?#S�?#��?#��?#o?#S�?#o?#o?"��?#o?"��?#S�?$�?$Z?$��?$�/?%�?$�/?$�/?%�?%�?%�?%�?$�/?$�/?$�/?$��?$��?$�/?%�?%�?%`B?%�?$��?$��?#�
?#�
?#��?#��?#��?#��?#S�?#��?#��?_|�?�V?�K�?�1?��D?��?��?�n�?�j?�l�?��?���?�v�?��w?���?��+?���?��?�z�?�j?�&�?�bN?�/?��D?�=q?�t�?v�+?mV?f$�?z��?}/?���?��\?{�m?�n�?�A�?y�?vE�?t��?t9X?s�F?st�?p�`?o�?kC�?h1'?`�?S�F?N{?6�+?+C�?)x�?#o?�?��?�h>��m>�33>� �>�ȴ?�T?
��?�9?$�?��? �>�?}>� �>� �>� �>�!>�>�1>�V>��T>�S�>���>�Ĝ>�n�>�>ٙ�>Ձ>�Q�>��h>�V>��h>��>��D>��m>�1'>��7>��>�X>���>�M�>�"�>��+>�z�>�O�>��>y�#>o��>s�F>m�h>V>I�^>Kƨ>N�>P�`>T��>R�>A�7>'�>�->�+>	7L=��m==�h=�"�=�
==��`=��=��=��=���=ě�=\=\=\=�v�=�j=�j=�j=�v�=��=��=���=���=���=���=��=�
==�;d=�;d=�/=�;d=�/=�`B=�S�=ě�=\=��`=��`=�E�=�t�=�7L=�+=�7L=�C�=��w=�1=��
=�\)=�hs=��-=�{=�
====��=�h=�x�=�S�=�l�=�`B=�;d=�l�=�l�=�G�=�/=�S�====�h=�=�x�=�x�=�l�=�S�=�S�=�`B=�S�=�S�=�G�=�;d=�;d=�"�=�
==��=��=��=���=���=Ƨ�=ě�=\=��=�v�=�j=�E�=�-=� �=���=��T=��-=�t�=m�h=Y�=D��=Y�=]/=P�`=<j=8Q�=49X=0 �=��=�w=��=\)=\)=C�=+=+<��<�/<���<���<�j<�9X<�o<T��<T��<T��<T��<D��<o;�`B;ě�;ě�;o            ��o��o�o�D���D�����
��`B�t��T���u��o��o��t���1��j��j��j��j��9X��9X��1��9X��j���ͽC��49X�<j�@��P�`�Y��e`B�u��o��+��+��\)���P���w���
���1��-��Q����������\�\�ȴ9���`��
=��"ѽ�;d��G���G���;d��;d��G���`B��`B��xս�h���#�%���C��V�bN������-� Ĝ�$�/�)��-V�0 ž/��333�8Q�=p��Kƨ�bMӾhr��l�D�o���r�!�t�j�u�u�{�m�~�۾�������7L���;��`��񪾔zᾕ����+��
=��b������"Ѿ�"Ѿ�/��Ĝ���徢�徣�
���/���y��xվ���������p���C���\)��녾�z�����������և+��
=�և+1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�33@�  @���A   A@  A`  A���A�  A�  A�  A�  A���A�  A�33B   B  B  B  B   B(  B0  B8  B@  BH  BPffBXffB`  Bh  Bo��Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D��D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DVfDV�fDW  DW� DX  DX� DY�B
�B
�B
?B
+B
B
1B
B
 B
B
#B
B
B
B
B
B
B
 �B
 �B
 �B
 B
 B
B
B
B
B
B
B
B
B
B
B
B
B
 B
 B
B
B
B
$B
B
/B
B
B
B
B
B
B
 �B
�B	��B
B�B
��B
��B7
B<nBBNBK�BOBP�BU�BYB^�Be�Bk�Bn�Br�Bt�Bs�Bt�Bu�BsBs�BuBu�BvMBv�Bu�Bu�Bu%BiB��B�2B��B��B��B��B��B�pB�WB�MB�HB��B��B��B��B��B��B�RB�OB~B{LB{Bx�BwGBv1Bt�BpxBo�Bp Bp�Bt�Bx(BxJBwjBv;Bu�Bt7Bs�Bs�Bs�BtZBs�Bs�BueBvBu�Bu�Bw	BteBusBw4BwdBs�Br�Br�Bs�Bt�BuTBw�B{lB{=B{[B{aBy�By�BxBBxBwyBw�Bv�Bu6Bt�BuBu�BuPBt�Bu�Bu�Bu�Bu�Bu�Bu�Bt=BtBtYBs;Br�Br�BsBs�Bs�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bx�Bx�Bx�By�Bz�Bz�Bz�B{�B{�B|�B|�B}�B}�B|�B}�B}�B~�B}�B|(B{�B|�B|�B}xB}�B(BqB~�B~�B�B�B�~B�B�B�:B�<B�IB�B�2B�JB��B�3B�VB�IB�B��B�7B�;B�IB�LB�PB�EB�TB�`B�JB�?B�VB�NB�[B�\B�RB�lB�mB�qB�sB�^B�hB�kB�{B�oB�oB�oB�pB�uB��B��B�vB��B��B��B��B�B��B��B�5B�aB��B��B�~B�|B�~B��B�iB��B��B�vB��B��B�wB��B��B��B��B��B��B��B��B�|B�{B�|B��B��B��B��B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�MB��B��B��B��B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B��B��B��B��B��B��B��B��B��?#�
?$Z?%�T?$�/?$Z?$Z?#��?#S�?#��?#��?#o?#S�?#o?#o?"��?#o?"��?#S�?$�?$Z?$��?$�/?%�?$�/?$�/?%�?%�?%�?%�?$�/?$�/?$�/?$��?$��?$�/?%�?%�?%`B?%�?$��?$��?#�
?#�
?#��?#��?#��?#��?#S�?#��?#��?_|�?�V?�K�?�1?��D?��?��?�n�?�j?�l�?��?���?�v�?��w?���?��+?���?��?�z�?�j?�&�?�bN?�/?��D?�=q?�t�?v�+?mV?f$�?z��?}/?���?��\?{�m?�n�?�A�?y�?vE�?t��?t9X?s�F?st�?p�`?o�?kC�?h1'?`�?S�F?N{?6�+?+C�?)x�?#o?�?��?�h>��m>�33>� �>�ȴ?�T?
��?�9?$�?��? �>�?}>� �>� �>� �>�!>�>�1>�V>��T>�S�>���>�Ĝ>�n�>�>ٙ�>Ձ>�Q�>��h>�V>��h>��>��D>��m>�1'>��7>��>�X>���>�M�>�"�>��+>�z�>�O�>��>y�#>o��>s�F>m�h>V>I�^>Kƨ>N�>P�`>T��>R�>A�7>'�>�->�+>	7L=��m==�h=�"�=�
==��`=��=��=��=���=ě�=\=\=\=�v�=�j=�j=�j=�v�=��=��=���=���=���=���=��=�
==�;d=�;d=�/=�;d=�/=�`B=�S�=ě�=\=��`=��`=�E�=�t�=�7L=�+=�7L=�C�=��w=�1=��
=�\)=�hs=��-=�{=�
====��=�h=�x�=�S�=�l�=�`B=�;d=�l�=�l�=�G�=�/=�S�====�h=�=�x�=�x�=�l�=�S�=�S�=�`B=�S�=�S�=�G�=�;d=�;d=�"�=�
==��=��=��=���=���=Ƨ�=ě�=\=��=�v�=�j=�E�=�-=� �=���=��T=��-=�t�=m�h=Y�=D��=Y�=]/=P�`=<j=8Q�=49X=0 �=��=�w=��=\)=\)=C�=+=+<��<�/<���<���<�j<�9X<�o<T��<T��<T��<T��<D��<o;�`B;ě�;ě�;o            ��o��o�o�D���D�����
��`B�t��T���u��o��o��t���1��j��j��j��j��9X��9X��1��9X��j���ͽC��49X�<j�@��P�`�Y��e`B�u��o��+��+��\)���P���w���
���1��-��Q����������\�\�ȴ9���`��
=��"ѽ�;d��G���G���;d��;d��G���`B��`B��xս�h���#�%���C��V�bN������-� Ĝ�$�/�)��-V�0 ž/��333�8Q�=p��Kƨ�bMӾhr��l�D�o���r�!�t�j�u�u�{�m�~�۾�������7L���;��`��񪾔zᾕ����+��
=��b������"Ѿ�"Ѿ�/��Ĝ���徢�徣�
���/���y��xվ���������p���C���\)��녾�z�����������և+��
=�և+1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�x<#��<#�-<#��<#�<#� <#�{<#�s<#�
<#؜<#�[<#�~<#�<#�d<#�f<#�d<#��<#��<#ה<#ׂ<#ׂ<#ך<#�j<#�
<#�z<#�
<#�
<#�<#נ<#�
<#�
<#�r<#�
<#�<#�z<#�<#�><#�><#��<#�<#ڱ<#�<#�w<#�
<#�
<#�
<#�j<#��<#�I<i3�<�6<2q�<H��<%�<#�<&Q�<'.r<'�<&�.<#�:<$��<%`�<$<$�<'�h<$��<#צ<#�Q<*�<0�<#��<$�<#�<$|�<(�.<*��<&p�<$�`<.$<$[<%�i<#�<%�Q<%ڐ<$Hj<%F4<$<#��<#�F<#�]<#��<$�<#�$<$;M<$$B<%}<(6h<$�^<1��<'��<#��<$��<$�<%��<%{�<*=�<$g�<#�^<$)3<&� <$r�<#�0<$<$%A<$�<$�q<$	r<#�<#�<#�<$(<#�)<#׻<$-�<#��<#ۇ<#�$<%w<#�n<#�<#��<)BP<$�<#׿<#�1<#�4<#ڃ<%~�<$��<$�<#�<$<%�.<$3<$/k<#�<#�5<$5\<$M�<$N�<$X<#ܸ<#�<$�@<$�<#�
<#۳<#�+<#ݫ<#��<$e<$�<$�<#�<$'�<$�<#�~<#�	<#��<#�<#۵<#�Z<#�
<#�<#��<#��<#�]<#�
<#�<#�2<#�S<#�
<#�
<#׀<#�N<#�<#�<#�<#�o<#ה<#��<#��<#�=<#�<#�F<#�\<#�H<#��<#�d<$:<#�K<#��<#�<$"�<$Z�<#�<#��<#׵<#װ<$�<#�<#�h<$�<#׎<#�<#�2<$�<$�<#�
<#׎<#�;<#إ<#��<#إ<#�<#�(<#ݛ<#�<#�<#س<#ڪ<#�<#�
<#�<#�~<#ױ<#�~<#�<#��<#�k<#�
<#�r<#�o<#�<#�\<#�~<#�<#�M<#ؗ<#�'<#ج<#�<#�~<#׻<#��<#׌<#׆<#׆<#ן<#ׁ<#�<#��<#מ<#�L<#ق<#��<#�y<$1<#�<#�Z<#�<#נ<#ۗ<#�:<#׫<#ׇ<#ת<#�N<#�z<#׆<#��<#�
<#׀<#�~<#�<#�<#�E<#ן<#��<#ؕ<#׵<#�$<#�b<#�<#�
<#�
<#י<#��<#ט<#�<#�<#��<#ؔ<#�
<#�
<#�v<#�<#�?<#�~<#�<#��<#��<#�<#�f<#�<#׆<#�<#�<#�7<#��<#�<#�<#�<#�v<#�<#�p<#׆<#ב<#�k<#��<$k<#��<#ת<#�!<#�%<#ۀ<#ޯ<#ށ<#��<#�<#��<#�?<#�l<#�<#��<#�<#��<#��<#�2<#�<#�
<#�v<#�<#ۀ<#ގ<#�e<#�	<#��<#ׁ<#�
<#�G<#�<#׏<#�<#�<#ټ<#��<#�F<#��<#��<#�.<#�O<#�?<#��<#�<#�[<#�5<#��<#�<#ۈ<#�<#�K<#�r<#�<#��<$7�<$�p<#��<#��<#�i<#�N<#��<#ׇ<#�<#�~<#�s<$N<#ߌ<#�<#�<#��<#�<#�<#�<#��<#י<#�<#�i<#� <#�<#޿<#��<#�r<#�<#��<#�"<#��<#�<#�p<#��<#��<%LA<%-`<#��<#�W<#�<#ם<#�<#��<#׊<#�w<#�S<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929163001201909291630012019092916300620190929163014202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�33@�33@�33G�O�G�O�G�O�G�O�DY�DY�DY�G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          