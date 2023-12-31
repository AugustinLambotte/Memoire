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
resolution        ?PbM���     �  L8   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  R�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Th   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  W�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^<   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  d�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  kT   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  n�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  p@   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  v�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  }X   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �D   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �D   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �D   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �D   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �,   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �H   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210225044611  20210225044611  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125488                          2C  D   APEX                            6229                            120210                          846 @�;�¿��1   @�;�¿��@P�n��O��5(�\)1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�33@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP�fDQfDQ�3DQ� B
T�B
S�B
T�B
T�B
T�B
T�B
S�B
Q�B
P�B
O�B
O�B
N�B
M�B
M�B
L�B
L�B
L�B
L�B
L�B
K�B
K�B
J�B
J�B
H�B
E�B
C�B
@�B
?}B
;dB
9XB
8RB
8RB
8RB
9XB
;dB
B�B
J�B
L�B
K�B
H�B
H�B
K�B
J�B
G�B
G�B
I�B
J�B
M�B
M�B
M�B
N�B
P�B
T�B
W
B
YB
]/B
dZB
k�B
|�B
��B
�B
�!B
�'B
�wB
�ZB
��B�B)�B<jBF�BH�BH�BH�BH�BH�BG�BE�BD�BG�BK�BO�BO�BS�BS�BVBT�BW
BW
BVBYBVBYBZB[#B[#B`BBcTBcTBe`BiyBk�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bn�Bn�Bn�Bn�Bn�Bn�Bm�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bl�Bl�Bl�Bl�Bm�Bm�Bm�Bn�Bo�Bp�Bp�Bp�Bp�Bq�Br�Bs�Bs�Bs�Bt�Bt�Bs�Bs�Bt�Bt�Bt�Bu�Bu�Bv�Bv�Bv�Bu�Bu�Bu�Bv�Bv�Bv�Bw�Bx�Bx�Bx�By�By�Bz�Bz�B{�B{�B{�B{�B{�Bz�B{�B{�B{�B{�B{�B|�B|�B|�B|�B}�B~�B~�B~�B~�B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�+B�+B�1B�1B�1B�1B�7B�=B�DB�DB�DB�DB�DB�DB�DB�JB�PB�PB�VB�VB�VB�VB�VB�\B�\B�bB�bB�\B�\B�\B�\B�bB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�{B�uB�uB�uB�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�j@�/@��@O�@�-@�@O�@�
@dZ@  @l�@��@��@�@��@��@
��@	�@b@  @��@5?@v�@�@(�?�;d?��m?�;d?�M�?�\?���?֧�?��?���?щ7?�l�?��;?�&�?�;d?�  ?͑h?�Ĝ?�hs?���?�z�?��?ě�?öF?�M�?�Ĝ?�j?�Q�?�x�?�dZ?��?���?�-?�33?��
?�b?�hs?���?��T?�V?��m?��#?���?�;d?��H?�Q�?�
=?�ȴ?�ff?�$�?��T?�?}?�`B?�hs?��/?�C�?��/?°!?�/?�dZ?���?�l�?���?�33?�bN?��?���?�  ?��?��/?��?pbN?d�?a%?Y�#?Y�#?X��?T9X?Qhs?O��?K?E��?BM�?:�H?7�P?7K�?3��?/�;?,�D?)x�?%��?"��?|�?�?�m?�m?�H?X?K�??}?33? �?1?��?�?��?%>��>��#>�E�>�~�>�G�>�t�>���>��>��7>�dZ>���>�9X>���>� �>��D>�r�>�A�>�5?>�/>��->��w>�M�>�Z>���>��>�bN>�=q>�O�>���>��>��>|�>q��>cS�>^5?>T��>O�;>O�;>Q�>O�;>L��>F��>C��><j>6E�>7K�>9X>5?}>!��>�P>z�>O�>hs>O�>1'>$�>$�>�>   =���=��=�l�=���=ȴ9=\=�v�=�Q�=�Q�=�E�=�E�=��=��w=���=��=��
=���=���=��
=��T=��
=��w=��-=�O�=y�#=m�h=P�`=<j='�=,1=,1=,1='�='�=�w=��=�P=t�=t�=#�
='�=��=�P=t�=C�=C�<��<��<���<���<�/<�h<�`B<�h<���<�j<���<���<�/<�`B<�<�<��=o<��<�h<�`B<�/<���<���<���<�/<�`B<�/<���<�9X<�t�<�o<#�
;�o;o;D��;o:�o��o��o��o��o���
��`B�t��T����C����
��9X��j���ͼ��ͼ�/��`B��h���C���P���'49X�49X�<j�D���D���L�ͽ]/�}�y�#�y�#��o��C���O߽�hs���P���P���P���P���w��1��Q���ͽ�������S���l���l���S���G���l������#���#���J���
=q�O߾O߾\)�hs��P����R� Ĝ�#�
�#�
�$�/�%�T�(�þ+�1&�5?}�8Q�:^5�<j�@��F��I�^�L�;L�;M��O�;�Q녾^5?�fff�j~��n���t�j�vȴ�y�#�|푾�$ݾ�7L��ƨ��O߾�V��\)��녾��Ͼ�zᾕ�������b���������(����R��;d���w��A���Ĝ��G����徣S����
��Z���/���T���y����羬1���h��{�����������׾�&龱����&龱&龲-��33���j����E�����E���ȴ���پ���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�$@��_@��_AZ�A>Z�A^Z�A~Z�A�-XA�-XA�-XA�-XA�-XA�-XA�-XA�-XB��B��B��B��B'��B/��B7��B?��BG��BO��BW��B_��Bg��Bo��Bw��B��B��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VB��VC�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C3�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���D ykD �kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kD	ykD	�kD
ykD
�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kDykD�kD ykD �kD!ykD!�kD"ykD"�kD#ykD#�kD$ykD$�kD%ykD%�kD&ykD&�kD'ykD'�kD(ykD(�kD)ykD)�kD*ykD*�kD+ykD+�kD,ykD,�kD-ykD-�kD.ykD.�kD/ykD/�kD0ykD0�kD1ykD1�kD2ykD2�kD3ykD3�kD4ykD4�kD5ykD5�kD6ykD6�kD7ykD7�kD8ykD8�kD9ykD9�kD:ykD:�kD;ykD;�kD<ykD<�kD=ykD=�kD>ykD>�kD?ykD?�kD@ykD@�kDAykDA�kDBykDB�kDCykDC�kDDykDD�kDEykDE�kDFykDF�kDGykDG�kDHykDH�kDIykDI�kDJykDJ�kDKykDK�kDLykDL�kDMykDM�kDNykDN�kDOykDO�kDP�DP��DQ��DQ�kB
T�B
T	B
T�B
T�B
T�B
UyB
UB
R_B
SWB
P[B
O�B
P>B
N;B
M�B
M�B
MxB
NB
M�B
L�B
L�B
LUB
J�B
K,B
JB
H�B
D�B
EB
>�B
@�B
<�B
9�B
9KB
9.B
9jB
9HB
?�B
J<B
M�B
QB
I�B
G�B
K�B
NB
H�B
GyB
I�B
K!B
NUB
NtB
OsB
PWB
P�B
TJB
V�B
X�B
[�B
c�B
egB
u�B
�!B
�XB
��B
��B
��B
�kB
�(BLB+�B=LBGBH�BH�BH�BH�BH�BG�BG8BI5BKIBN(BP�BQ�BT�BUGBV1BU�BW�BX(BY7BZ�B]KB\6B[WB\WB^�Bb�Bc�Bd�BecBi�BlZBmBl�BmkBm�Bm8Bm�Bm)Bl�Bm'BlPBm-Bm#BmDBmBm,BmBm�Bn�Bn�Bn�Bn�Bn�Bn�Bn,Bo_Bp0Bo�BpkBpBpBo�Bo�Bo�BnwBn�BmfBm-Bl�BmBl�Bl�Bm�Bl�Bl�Bl�BmOBm�Bm�Bm�BngBoaBpvBp�BquBqzBr:BrmBs�BtBtbBt�BuDBtdBs�Bu*Bt�Bt�Bu�Bu�Bv�BwBv�BvBvBu�Bv�Bw Bw�BxJBx�By'Bx�BzBzBz�Bz�B{�B|&B|B|B|%B{SB|1B|B|B|B{�B|�B|�B}BB}B}�BBBB~�B�B��B�B�%B�B�lB�lB�/B�^B�FB�CB��B�B�B�B�B�'B�B�B�B�B��B�B�=B�'B�&B�2B�!B�FB�)B�bB�+B�!B�B�<B�'B�cB�PB�&B�6B�8B�7B�,B�CB�8B�;B�VB�fB�^B�bB�bB�WB�VB�JB�OB�kB�nB��B��B�xB��B��B�zB�[B�wB�zB��B��B�pB�pB�|B��B��B��B��B��B��B��B��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�yB��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�0B��B��B��B��B��B��B��B�TB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�j@�/@��@O�@�-@�@O�@�
@dZ@  @l�@��@��@�@��@��@
��@	�@b@  @��@5?@v�@�@(�?�;d?��m?�;d?�M�?�\?���?֧�?��?���?щ7?�l�?��;?�&�?�;d?�  ?͑h?�Ĝ?�hs?���?�z�?��?ě�?öF?�M�?�Ĝ?�j?�Q�?�x�?�dZ?��?���?�-?�33?��
?�b?�hs?���?��T?�V?��m?��#?���?�;d?��H?�Q�?�
=?�ȴ?�ff?�$�?��T?�?}?�`B?�hs?��/?�C�?��/?°!?�/?�dZ?���?�l�?���?�33?�bN?��?���?�  ?��?��/?��?pbN?d�?a%?Y�#?Y�#?X��?T9X?Qhs?O��?K?E��?BM�?:�H?7�P?7K�?3��?/�;?,�D?)x�?%��?"��?|�?�?�m?�m?�H?X?K�??}?33? �?1?��?�?��?%>��>��#>�E�>�~�>�G�>�t�>���>��>��7>�dZ>���>�9X>���>� �>��D>�r�>�A�>�5?>�/>��->��w>�M�>�Z>���>��>�bN>�=q>�O�>���>��>��>|�>q��>cS�>^5?>T��>O�;>O�;>Q�>O�;>L��>F��>C��><j>6E�>7K�>9X>5?}>!��>�P>z�>O�>hs>O�>1'>$�>$�>�>   =���=��=�l�=���=ȴ9=\=�v�=�Q�=�Q�=�E�=�E�=��=��w=���=��=��
=���=���=��
=��T=��
=��w=��-=�O�=y�#=m�h=P�`=<j='�=,1=,1=,1='�='�=�w=��=�P=t�=t�=#�
='�=��=�P=t�=C�=C�<��<��<���<���<�/<�h<�`B<�h<���<�j<���<���<�/<�`B<�<�<��=o<��<�h<�`B<�/<���<���<���<�/<�`B<�/<���<�9X<�t�<�o<#�
;�o;o;D��;o:�o��o��o��o��o���
��`B�t��T����C����
��9X��j���ͼ��ͼ�/��`B��h���C���P���'49X�49X�<j�D���D���L�ͽ]/�}�y�#�y�#��o��C���O߽�hs���P���P���P���P���w��1��Q���ͽ�������S���l���l���S���G���l������#���#���J���
=q�O߾O߾\)�hs��P����R� Ĝ�#�
�#�
�$�/�%�T�(�þ+�1&�5?}�8Q�:^5�<j�@��F��I�^�L�;L�;M��O�;�Q녾^5?�fff�j~��n���t�j�vȴ�y�#�|푾�$ݾ�7L��ƨ��O߾�V��\)��녾��Ͼ�zᾕ�������b���������(����R��;d���w��A���Ĝ��G����徣S����
��Z���/���T���y����羬1���h��{�����������׾�&龱����&龱&龲-��33���j����E�����E���ȴ���پ���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#�<#�I<#�?<#�C<#�U<$�<#�<((�<#��<#��<%/�<#�<#��<$��<$�<$ت<$1C<#�!<$2N<$x<#�(<#��<%%�<+�9<%8�<2�'<$b�<9�k<+k�<$�K<$u�<$Pk<#�<'��<+:<$<$G8<81 <$g�<$�M<#�9<,�<$�<#�C<#��<#��<#�<$:<%��<%c�<#��<$P<$E<$�<%�n<$6�<?xW<J:<,�<'��<%�Q<)��<k�W<#�
<&2g<;�S<%λ<$\6<#�<#�7<#�X<#�L<#�Z<#��<#׎<%��<3.b<-]�<'�T<$V�<&�<$(D<%z<#ٙ<$f�<$#<$��<+Z[<%��<Hd <*� <$�<$��<-xZ<'�!<$�<%<#פ<#ڵ<$LR<$�<#�<$Z�<$�y<$!�<%)�<$�<#�!<$`<$A<$d<$�<$-�<$�<$�<$c<#ؗ<#��<#�m<#��<#�`<#�<#�<$�<$<B<$
�<#ܸ<$Cm<#�"<#�6<#�<#�]<$�<$a�<%&<$S�<$.<#��<$^<#�<#�<#܀<#س<#�v<#�G<$8�<#ړ<#�;<#�'<#�(<#�<#��<#�w<$G�<$M�<$�<#�Q<#�<#��<$ Z<#�J<$b<$"�<#�!<#�<#�B<#׃<#�Q<#�]<#أ<#��<#�|<#�
<#�<#�8<#�<#�~<$e�<#��<#س<#�<#��<#�-<#ݜ<#�-<#�x<#�<#ݘ<#ڗ<#ب<#�D<#�<#�<#ؑ<#�V<#�C<#��<#�)<#׆<#�<#�)<#�%<#�=<#�><#�<#��<#ـ<#�g<#�+<#�C<#�
<#�3<#�J<#��<#�<#��<#� <#س<#��<#��<#�<#׵<#�Z<#�<#�<#�<#��<#�><#ي<#�,<#�<#�<#�@<#ק<#��<#׉<#�5<#��<#�<#܉<#�<#�-<#ڌ<#�?<#۾<#��<#ٍ<#ٟ<#��<#��<#ي<#�<#�<#�$<#�<#�<#�<#��<#��<#ى<#پ<#�<#�<#�L<#ڷ<#�v<#�:<#��<#�;<#ٔ<#�<#�#<#�P<#�K<#׶<#��<#�<#�M<#�Z<#�<#�	<#�t<#�Z<#�:<#�f<#׮<#�F<#�<#�<#�<#�t<#�g<#�<#�<#�S<#׼<#�D<#�A<#ש<#�V<#��<#�7<#��<#׸<#؎<#��<#�<#�M<#�N<#��<#��<#צ<#ڛ<#��<#�<#��<#��<#׈<#�/<#�Q<#��<#��<#�6<#ؓ<#�<#�9<#ײ<#�O<#�u<#ؠ<#�b<#�V<#צ<#�C<#�k<#�<#��<#�;<#�T<#�L<#��<#�<#�<#�\<#�n<#�_<#��<#؂<#�S<#�]<#ھ<#�*<#ث<#�}<#��<#�<#�K<#ן<$�<#��<#��<#�U<#�6<#�a<#�}<#��<$/	<#��<#ݯ<#ؓ<#�V<#�d<#ݒ<#ڂ<#�<#�D<#�<#ן<#�@<#�v<#�<#ڮ<#�@<#�<#�4<#�<#�<#�<#�N<#�<#�<#�%<#�
<#�K<#�P<#ׂ<#��<#ڝ<#�n<#�<#�<#ל<#�1<#�<#�<#�I<#��<#�K<#�W<#�j<#�O<#�<#�\<#�<#�<#�<#��<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.10286                                                                                                                                                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929162717201909291627172019092916272220190929162729202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�33@�33@�33G�O�G�O�G�O�G�O�DQ� DQ� DQ� G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          