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
_FillValue                 �  ST   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  U   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  Xd   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  _    TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  nH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  q�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  xd   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �      	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �<   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �<   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �<   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �<   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �    HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �$   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �@   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210225044033  20210225044033  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125490                          2C  D   APEX                            6229                            120210                          846 @�@�P�� 1   @�@�P�� @P�j~��#�4��
=p�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @@  @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C)�fC,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C��C��C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<y�D<��D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS�fDT  DT� DU  DU� DV  DV� DWfDWs3B	l�B	m�B	m�B	n�B	p�B	s�B	t�B	x�B	}�B	�B	�B	�B	�B	�B	�B	�+B	�1B	�=B	�=B	�=B	�DB	�JB	�PB	�bB	�uB	��B	��B	��B	��B	��B	��B	��B	��B	��B	�;B	��B
B

=B
�B
,B
?}B
N�B
XB
_;B
cTB
gmB
iyB
m�B
o�B
p�B
r�B
u�B
y�B
|�B
�B
�1B
�VB
��B
��B
��B
�B
�B
�?B
��B
�#B
�B
�B
��BB�B�B�B!�B#�B(�B0!B@�BP�BZB]/BgmBjBk�Bk�Bl�Bl�Bl�Bk�BjBhsBl�Bl�Bl�Bm�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bn�Bo�Bo�Bo�Bo�Bo�Bp�Bo�Bo�Bn�Bo�Bp�Bp�Bp�Bq�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bq�Bq�Bp�Bp�Bp�Bq�Bs�Bs�Br�Br�Br�Br�Br�Br�Br�Br�Bq�Bq�Br�Br�Bs�Bs�Bs�Bt�Bu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bw�Bw�Bw�Bx�Bx�Bx�Bx�By�By�Bz�Bz�Bz�Bz�Bz�B{�B{�B{�B|�B|�B}�B}�B~�B~�B~�B~�B� B� B� B� B� B�B�B�B�B�B� B� B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�+B�+B�1B�1B�7B�7B�=B�=B�=B�DB�=B�=B�DB�DB�DB�JB�JB�PB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�\B�VB�\B�\B�\B�\B�bB�\B�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�uB�uB�uB�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?t��?u?}?t��?u�?w�P?z�H?{�m?�w?�Q�?���?�%?�&�?���?���?���?���?�5??�  ?�Ĝ?��`?��`?���?�t�?��+?���?���?�{?��?���?�o?�t�?��F?�V?��7?�p�?�;d?��?�t�?���?�b?�{?�M�?�v�?��?���?�j?��D?���?�/?�p�?�p�?�p�?��?�5??�;d?� �?���?�&�?���?���?��7?�&�?�|�?�;d?���?�K�?��9?�7L?��?��-?�v�?���?�\)?��R?�"�?�+?���?���?��/?� �?���?x�u?pbN?kC�?e�?`Ĝ?\�?W�P?PbN?C�
?$�?��?�j?hs?I�?��?1'?ff?��>�j>�X>�ȴ>�>��j>�->��>�V>�r�>�Z>�;d>�
=>�bN>ȴ9>��>��>�o>�p�>�j>��H>��F>���>�{>��D>��>���>�x�>��T>���>��+>��>�>��>�t�>�
=>�z�>��>�\)>�I�>��9>��\>y�#>u>j~�>dZ>cS�>bM�>Xb>Q�>N�>G�>49X>�R>�+>�+>��>&�y>'�>'�>�>n�>�>\)>
=q>J=�F=�h=�=�S�=��=�
==���=�"�=�`B=�x�=�S�=�S�=�G�=��=���=�Q�=�^5=�j=�9X=�{=���=��=��=��T=��T=��
=���=���=���=��P=�t�=q��=m�h=y�#=�o=}�=e`B=ix�=y�#=y�#=�%=�o=�%=�7L=�+=�+=�+=�7L=�7L=�+=��=y�#=P�`=<j=49X=<j=8Q�=,1=T��=Y�=H�9=@�=49X=,1='�=,1='�=�P=�w=\)=t�=�P=#�
=�P=�P=49X=49X=#�
=0 �='�=C�<�h<�`B<�=o=+=��=�w=�w=��=�P=��=�w=��=�P=C�=o<�<�<�<�h<�`B<�`B<�/<�j<�1<��
<�1<�1<�t�<u<D��<t�;ě�;ě�;�o;D��;�o;�o    ��o��`B�#�
�49X�T���e`B�u�u�e`B�u��C���t���t����㼣�
�ě���/��/����P�8Q�<j�<j�<j�@��@��P�`�aG��ixսu�y�#��o��C���t����P�������-���
�� Ž�E���j��vɽ���ȴ9������`������;d��xս���F���#���%���$ݾ$ݾ1'�C��O߾V�V�V�V�\)�\)�hs�t�����+��u������ Ĝ�$�/�&�y�(�þ-V�0 ž49X�7KǾ8Q�:^5�;dZ�>vɾA�7�B�\�C���D���I�^�Q녾Xb�]/�dZ�l�D�s�F�w�پx���y�#�z�H�z�H�{�m��  ���������=q��ƨ��I���O߾��;��񪾕������(����w��G����徤Z��l���~�������h����������j����Q쾹X��^5��dZ���m��j��j��j��j��푾�푾�j��j��j��j��j��j��j��j���m11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @@  @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C)�fC,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C��C��C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<y�D<��D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS�fDT  DT� DU  DU� DV  DV� DWfDWs3B	l�B	m�B	m�B	nCB	pB	s�B	s�B	u�B	|BB	��B	�B	��B	�YB	��B	��B	�BB	��B	��B	�8B	�DB	��B	��B	�<B	� B	�B	�fB	�1B	��B	�?B	��B	��B	��B	��B	�UB	ޫB	�B
 mB
	�B
hB
)�B
=�B
P'B
X�B
_IB
czB
g`B
ifB
mnB
o�B
p�B
r�B
u�B
y�B
|�B
��B
�B
�,B
�aB
��B
��B
�4B
��B
��B
��B
�*B
��B
�tB
�BZB9BcB�B"B%3B*�B2�BB�BSB\BbBi�Bl
Bl�Bl�BmVBmUBm�Bl�Bl�BnJBnkBm�Bm4Bn�Bo7Bn�Bn�Bo%Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�BpBpBpBpeBpBBoSBn�Bo�Bo�Bp"Bo�Bo�BoDBo�Bo�Bo�Bo�Bo�Bp�Bo�Bp�Bn�BogBp�BqBp�BqXBr�Br�Br�Br�BsBsEBs8Br�Bs9Br�Br�Br�Bs+Br�Br�Bs	Br�Br�BqBp�Bp~BqBs�Bs�Bs8Bs)Br�Br�Br�BsBsBr�Bq�Bq�Br�Br�Bs�Bs�Bs|Bt�Bu�Bu�Bv�Bv�Bw"Bw6Bv�Bw�BxBw�Bw�Bw�Bx�Bx�Bx�Bx�By�BzBz�Bz�Bz�B{|Bz�B{�B{�B|B}1B|�B}�B}�B~�B~�BB~�B�	B� B� B�B�B�B�B�:B�B�CB�B�B�B�(B��B�B�GB�4B�>B�3B�*B�B�/B�RB�B�XB�!B�!B�B�VB�3B��B�;B�mB�B�XB��B�|B�OB�,B�.B�:B�B�AB�PB�\B�]B�IB�KB�bB�dB�|B�pB�oB�ZB�ZB�eB�fB�]B�jB��B�xB�nB�TB�cB��B��B��B��B��B�dB�zB�nB�XB�iB��B��B��B��B�vB��B�uB�tB�iB�aB�{B��B��B�vB��B��B��B��B�xB��B��B��B��B�yB�vB��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?t��?u?}?t��?u�?w�P?z�H?{�m?�w?�Q�?���?�%?�&�?���?���?���?���?�5??�  ?�Ĝ?��`?��`?���?�t�?��+?���?���?�{?��?���?�o?�t�?��F?�V?��7?�p�?�;d?��?�t�?���?�b?�{?�M�?�v�?��?���?�j?��D?���?�/?�p�?�p�?�p�?��?�5??�;d?� �?���?�&�?���?���?��7?�&�?�|�?�;d?���?�K�?��9?�7L?��?��-?�v�?���?�\)?��R?�"�?�+?���?���?��/?� �?���?x�u?pbN?kC�?e�?`Ĝ?\�?W�P?PbN?C�
?$�?��?�j?hs?I�?��?1'?ff?��>�j>�X>�ȴ>�>��j>�->��>�V>�r�>�Z>�;d>�
=>�bN>ȴ9>��>��>�o>�p�>�j>��H>��F>���>�{>��D>��>���>�x�>��T>���>��+>��>�>��>�t�>�
=>�z�>��>�\)>�I�>��9>��\>y�#>u>j~�>dZ>cS�>bM�>Xb>Q�>N�>G�>49X>�R>�+>�+>��>&�y>'�>'�>�>n�>�>\)>
=q>J=�F=�h=�=�S�=��=�
==���=�"�=�`B=�x�=�S�=�S�=�G�=��=���=�Q�=�^5=�j=�9X=�{=���=��=��=��T=��T=��
=���=���=���=��P=�t�=q��=m�h=y�#=�o=}�=e`B=ix�=y�#=y�#=�%=�o=�%=�7L=�+=�+=�+=�7L=�7L=�+=��=y�#=P�`=<j=49X=<j=8Q�=,1=T��=Y�=H�9=@�=49X=,1='�=,1='�=�P=�w=\)=t�=�P=#�
=�P=�P=49X=49X=#�
=0 �='�=C�<�h<�`B<�=o=+=��=�w=�w=��=�P=��=�w=��=�P=C�=o<�<�<�<�h<�`B<�`B<�/<�j<�1<��
<�1<�1<�t�<u<D��<t�;ě�;ě�;�o;D��;�o;�o    ��o��`B�#�
�49X�T���e`B�u�u�e`B�u��C���t���t����㼣�
�ě���/��/����P�8Q�<j�<j�<j�@��@��P�`�aG��ixսu�y�#��o��C���t����P�������-���
�� Ž�E���j��vɽ���ȴ9������`������;d��xս���F���#���%���$ݾ$ݾ1'�C��O߾V�V�V�V�\)�\)�hs�t�����+��u������ Ĝ�$�/�&�y�(�þ-V�0 ž49X�7KǾ8Q�:^5�;dZ�>vɾA�7�B�\�C���D���I�^�Q녾Xb�]/�dZ�l�D�s�F�w�پx���y�#�z�H�z�H�{�m��  ���������=q��ƨ��I���O߾��;��񪾕������(����w��G����徤Z��l���~�������h����������j����Q쾹X��^5��dZ���m��j��j��j��j��푾�푾�j��j��j��j��j��j��j��j���m11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#׃<#�=<#�<#��<$�<#�<$X�<*�h<&,S<%�t<#ב<$\$<$��<#�<#�><$�:<$,�<#��<#נ<#�<#�^<$&�<$� <%#t<%VB<#��<#�<$�5<$}<#۱<#�<'P�<5/�<%|x<$�<$�0<$$<$Q<$�<'�:<%��<%!�<$9<#�]<#ڡ<#��<#ؗ<#ۗ<#�9<#�<#�
<#�4<#��<#��<#�<#�><#�h<#��<#�<#�<#�7<$'8<#�v<$�<&��<$�<#�s<$J�<$O�<#�Z<#�w<#��<#�m<%D<%�8<(��<'�<'p<&�<5�s<(xq<%�n<$�o<$��<$S�<$R�<$�Q<%_�<(�<<��<&��<$��<$-�<$�<$#�<#܌<#�<$�<$�(<#�<#�<#�<#�+<#�=<#��<#�	<#�9<#��<$�<$O
<$(X<$A <#�h<#�<#߰<$/<#ٸ<#��<$1<#�S<#�<#ۮ<#�<#�2<#�{<#�x<$�J<#�V<#�|<#�<#�h<#ڏ<#��<#�<#�<#�H<#�<#��<$8<$A<#�r<$�<#�<#׸<#��<$�<#�,<#�'<#�)<$��<$�m<#�w<#�<#ۆ<$ �<#��<#�<$x<$^<#ڄ<#�<#��<#�D<#��<#�$<#��<#�k<#�<#ׅ<#�F<#�4<#�^<#غ<#ڨ<#�<#�c<#��<#�.<#��<#� <#׎<#�d<#�V<#��<#��<#�<#�~<#�<#׽<#�f<#�b<#�w<#ד<#ٴ<$ X<#�f<#۔<#��<#١<#�h<#�s<#޳<#�
<#�!<#�w<#�R<#ܶ<#�D<#�
<#�
<#�v<#�<#�]<#פ<#�U<$�<#�_<#��<#؃<#�Q<#ژ<$|<#�^<#�|<#�/<#�5<#�<#�d<#�g<#��<#�.<#ؙ<#�<#�`<#�e<#��<#�<#�<#�<#�<#��<#�<#�1<#�<#�<#�e<#��<#؁<#��<#�v<#��<#�
<#�~<#׋<#ז<#�u<#�v<#ט<#�L<#�<#��<#�<#�<#״<#�Q<#�<#כ<#ނ<#�k<#�<#ע<#�<#�H<#�f<#�f<#�c<#�<#�<#��<#�}<#�b<#�
<#�P<#ާ<#�v<#�3<#ל<#��<#׎<#�}<#�<#ע<#�~<#�-<#�~<#�<#�<#ע<#�{<#�%<#�<#ۅ<#�N<#� <#�\<#�<#�
<#ג<#�<#�P<#ކ<#�$<#�4<#�|<#�O<#ޫ<#ބ<#�&<#�f<#��<#ۣ<#��<#ۜ<#�3<#ח<#ע<#��<#ڷ<#�j<#�|<#�'<#�<#�M<#�<#��<#�D<#�<#�<#�:<#��<#�<#�b<#�6<#�<#�l<#�
<#�
<#�
<#�~<#�<#��<#��<#��<#ד<#��<#ף<#�_<#��<#��<#�<#�<#ޑ<#ۍ<#ދ<#�N<#ל<#��<#ן<#�4<#�4<#׎<#׈<#׮<#�<#�<#�<#�<#�<#�!<#�<#ޑ<#ן<#׆<#�~<#�
<#ם<#��<#�W<#�S<$Q<#��<#ף<#�'<#�_<#�=<#�<#��<#��<#�<#��<#�o<#۾<#�<#��<#�k<#�x<#ۯ<#�=<#�<#�}<#�^<#��<#��<#��<#׏<#ז<#�<#�
<#�
<#�V<#�
<#ח<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#׃<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929162837201909291628372019092916284120190929162848202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@@  @@  @@  G�O�G�O�G�O�G�O�DWs3DWs3DWs3G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          