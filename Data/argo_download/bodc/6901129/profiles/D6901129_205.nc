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
_FillValue                 �  Sx   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  U,   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  _\   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  f$   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  n�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  pT   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  r   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  x�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �`   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �0   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �L   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �h   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �D   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �4   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �P   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �l   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225044637  20210225044637  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125491                          2C  D   APEX                            6229                            120210                          846 @�C
F�� 1   @�C
F�� @P�p��
=�4��;dZ1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @333@�  @�33A   A   A@  A`  A���A�  A���A�  A�33A�  A�  A���B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bg��Bp  Bw��B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Ca�fCc�fCf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7y�D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV�fDWfDW� DX  DX� DY  B�B�B�B�B�B�B�B�B�B�B�B�B�B�B"�B-B;dBq�B�\B�fB	,B	ffB	~�B	�1B	�hB	�bB	�{B	��B	��B	��B	�/B
  B
�B
?}B
dZB
q�B
{�B
x�B
�B
�1B
�B
� B
y�B
u�B
o�B
ffB
`BB
\)B
W
B
Q�B
^5B
cTB
hsB
q�B
v�B
z�B
�B
}�B
x�B
x�B
y�B
y�B
z�B
�B
�DB
��B
��B
�-B
�RB
��B
��B
��B
��B
��B
�/B
�HB
�ZB
��B+B�B�B"�B;dBI�BQ�BS�BVBW
BZB\)B\)B\)B]/B^5B`BBbNBcTBcTBdZBffBiyBjBjBk�BjBl�Bn�Bn�Bn�Bp�Bq�Br�Bt�Bt�Bt�Bt�Bu�Bv�Bv�Bv�Bu�Bv�Bv�Bv�Bw�Bw�Bv�Bv�Bt�Bt�Bt�Bu�Bt�Bt�Bt�Bt�Bs�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bv�Bv�Bv�Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bw�Bw�Bw�Bw�Bw�Bw�Bv�Bv�Bv�Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�By�By�By�By�By�By�Bz�By�By�Bz�Bz�B{�B{�B{�B{�B{�B|�B}�B}�B}�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B~�B~�B~�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�B�B�B�B�B�%B�%B�+B�1B�7B�7B�7B�=B�=B�=B�DB�DB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�PB�PB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�bB�bB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?Q��?R-?T9X?T�j?St�?R�?Qhs?Q�?Rn�?R�?R�?R�!?St�?Vȴ?[dZ?`A�?g�?��P?�bN?��\?�@	�@ �@�@!�#@�R@9X@   @+@&ȴ@0Ĝ@B�H@I�7@TZ@b��@f�+@j�H@fȴ@e�@f{@a��@[�@Sƨ@Lj@H��@?�w@9x�@6�@2�H@.ff@#33@!��@�@��@�@@ff@v�?���?���?�  ?�dZ?�7L?�b?ԛ�?�"�?�z�?�  ?�?��?�=q?���?��?�x�?��m?��?�!?��?�r�?�
=?�33?�o@x�@��?�(�?���?�$�?��?��?��?�Ĝ?�|�?�X?ҏ\?�x�?�M�?���?��D?��?�hs?�C�?��u?�
=?��?�9X?�|�?���?��H?��9?�Z?�G�?���?�r�?���?���?��!?}�?w��?s33?q��?o�?h��?a%?U�?Qhs?P�`?P �?O��?MV?A��?2n�?)��?'�?%��?!�7?�#?��?	�^?J? �>���>��>�j>�X>��j>>�ff>߾w>ܬ>ؓu>�>��`>���>Ƨ�>�o>\>�%>�p�>�j>��m>��m>���>���>��T>�/>��>�>�bN>���>}�>j~�>dZ>aG�>^5?>["�>Xb>V>O�;>G�>@�>?|�>;dZ>5?}>.{>)��>"��>�->��>I�>	7L>�>%=��=��=���=��=��==�l�=��=��=Ƨ�=�j=� �=�1=��
=���=�\)=�7L=�7L=�\)=�t�=�hs=�\)=�7L=�o=q��=ix�=e`B=aG�=aG�=ix�=e`B=aG�=L��=49X=�w=C�='�=L��=@�=@�=P�`=Y�=H�9=L��=D��=,1=+<���<�9X<��=C�=o<ě�<�o<o<�o<���<�o<e`B;�`B;�o    �D���D���D����o�o��o��o��o��o;o;�o;�o;o;�o;ě�;ě�<t�<t�;ě�;��
;D��;�o�D���49X��o�o�49X�u��1�ě��ě���9X�ě���`B��`B��`B��`B��`B�����o�t����0 Ž49X�8Q�<j�D���H�9�P�`�P�`�P�`�P�`�T���T���T���T���T���Y��Y��Y��]/�aG��aG��e`B�e`B�e`B�e`B�ixսixսe`B�ixսq���q���}󶽇+��+��+�����o��%��o��+��hs�������-���w���
���罬1��1�� Ž�E��\�\�\�\������
=��S���S���S���l���������+�\)��� Ĝ�.{�<j�?|�>vɾ?|�F��H�9�I�^�Kƨ�Kƨ�J���Kƨ�W
=�]/�_;d�cS��ixվr�!�vȴ�{�m��  ��J��o����������1'���^����������ƨ��O߾�\)��n���t����Ͼ���������G���xվ�~���V������������������������ ž� ž��׾�&龱&龱&龱��������������������&龱&龱&龱��������������������������11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @333@�  @�33A   A   A@  A`  A���A�  A���A�  A�33A�  A�  A���B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bg��Bp  Bw��B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Ca�fCc�fCf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7y�D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV�fDWfDW� DX  DX� DY  BnBCB�B�B�B�B�B�B�B�B�BBB�B!�B+�B4�Bm�B�dB�B	kB	Z�B	{ B	��B	��B	�&B	�,B	�B	��B	��B	�B	��B
�B
5kB
a+B
n�B
~�B
y�B
��B
�B
��B
��B
 B
x�B
vB
j�B
b^B
_B
Z|B
Y�B
_rB
e;B
j�B
s�B
y,B
��B
��B
�@B
~%B
|�B
{�B
z�B
{qB
�8B
��B
�B
��B
��B
��B
��B
��B
�*B
��B
�)B
ޕB
�B
��B
�~B�B�B�BWB9.BN[BR�BUeBW�BZ6B\�B\�B\�B^yB_�Ba�Bb�Bc�Bd$BeBf�Bh�BjyBkBk8Bk�BlDBm�Bn�BonBp8Bq�BsYBtdBumBuYBu�Bv!Bv�Bw�BwBwIBv�BxLBx�Bw�Bw�Bw�Bv�BwUBv�Bw�Bv[Bv*Bu1Bu�Bv.BvJBu:Bv(BuBt�Bt�Bt�BuBu,Bu=Bu�Bu`BvBv(Bw
BwGBw`Bx5By)Bx�Bx�By)Bx�Bx�Bx�By Bx�Bx�Bx�Bx6BxBxYBw�BwqBw�BxBx�Bx�Bx�Bx�Bx�By!By9By+Bx�ByBy!By,By
By+ByByBytBzBzBzBy�By�Bz Bz By�By�BzBz4B{ BzBzB{(Bz�B|B|B|1B|B{�B|�B}�B~ B~BB!B4B~B~B~B}�B}�B~B~B~2B~?B~6B4B~�B~�B�!B�B��B��B�7B�B�'B�XB�~B��B��B��B��B�0B�|B�~B�tB��B��B�EB�?B�vB�GB�QB�DB�"B�#B�3B�B�%B�4B�7B�8B�B�%B�@B�XB�-B�1B�KB�*B�QB�uB�^B�hB�LB�B�DB��B�=B��B��B��B�yB�TB�?B�pB��B�[B�\B�\B�]B�tB�iB�jB��B�xB��B�nB�oB�sB��B�vB��B�iB�hB�lB�yB�oB�oB�oB�pB�{B�oB�pB�{B�{B�pB�B�vB�uB�vB��B�xB�mB��B��B�}B��B��B��B��B�xB�xB�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�(B��B�4B�AB��B��B��B��B��B��B��B��B��B��B�!B��B��B��B��B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�dB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?Q��?R-?T9X?T�j?St�?R�?Qhs?Q�?Rn�?R�?R�?R�!?St�?Vȴ?[dZ?`A�?g�?��P?�bN?��\?�@	�@ �@�@!�#@�R@9X@   @+@&ȴ@0Ĝ@B�H@I�7@TZ@b��@f�+@j�H@fȴ@e�@f{@a��@[�@Sƨ@Lj@H��@?�w@9x�@6�@2�H@.ff@#33@!��@�@��@�@@ff@v�?���?���?�  ?�dZ?�7L?�b?ԛ�?�"�?�z�?�  ?�?��?�=q?���?��?�x�?��m?��?�!?��?�r�?�
=?�33?�o@x�@��?�(�?���?�$�?��?��?��?�Ĝ?�|�?�X?ҏ\?�x�?�M�?���?��D?��?�hs?�C�?��u?�
=?��?�9X?�|�?���?��H?��9?�Z?�G�?���?�r�?���?���?��!?}�?w��?s33?q��?o�?h��?a%?U�?Qhs?P�`?P �?O��?MV?A��?2n�?)��?'�?%��?!�7?�#?��?	�^?J? �>���>��>�j>�X>��j>>�ff>߾w>ܬ>ؓu>�>��`>���>Ƨ�>�o>\>�%>�p�>�j>��m>��m>���>���>��T>�/>��>�>�bN>���>}�>j~�>dZ>aG�>^5?>["�>Xb>V>O�;>G�>@�>?|�>;dZ>5?}>.{>)��>"��>�->��>I�>	7L>�>%=��=��=���=��=��==�l�=��=��=Ƨ�=�j=� �=�1=��
=���=�\)=�7L=�7L=�\)=�t�=�hs=�\)=�7L=�o=q��=ix�=e`B=aG�=aG�=ix�=e`B=aG�=L��=49X=�w=C�='�=L��=@�=@�=P�`=Y�=H�9=L��=D��=,1=+<���<�9X<��=C�=o<ě�<�o<o<�o<���<�o<e`B;�`B;�o    �D���D���D����o�o��o��o��o��o;o;�o;�o;o;�o;ě�;ě�<t�<t�;ě�;��
;D��;�o�D���49X��o�o�49X�u��1�ě��ě���9X�ě���`B��`B��`B��`B��`B�����o�t����0 Ž49X�8Q�<j�D���H�9�P�`�P�`�P�`�P�`�T���T���T���T���T���Y��Y��Y��]/�aG��aG��e`B�e`B�e`B�e`B�ixսixսe`B�ixսq���q���}󶽇+��+��+�����o��%��o��+��hs�������-���w���
���罬1��1�� Ž�E��\�\�\�\������
=��S���S���S���l���������+�\)��� Ĝ�.{�<j�?|�>vɾ?|�F��H�9�I�^�Kƨ�Kƨ�J���Kƨ�W
=�]/�_;d�cS��ixվr�!�vȴ�{�m��  ��J��o����������1'���^����������ƨ��O߾�\)��n���t����Ͼ���������G���xվ�~���V������������������������ ž� ž��׾�&龱&龱&龱��������������������&龱&龱&龱��������������������������11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<#�1<#ٖ<#�<#ٿ<#�<#�<#�y<#�=<#�<#�c<#�1<$�<$\�<$�&<%�<B��<0�<|�0<�<�<���<vue</��<+�.<'d�<&)�<((<#ޑ<7�<M�i<~<<7��<N̩<e��<+�l<*y�<)(F<$q�<#�<)�<3��<:��<8*�<)��<A�~<2 <':p<*=H<,��<O��<% 2<&��<(�<'7`<(,�<:�<<��<I�<8|2</�!<&�<$b�<$<$�<(��<-#<0��<%��<%��<$~�<#�<#��<#��<$\�<%W�<&S<%�J<#�k<$)<%$�<#��<9��<'��<3��<$|�<%j�<&!i<+�/<)<$<$<'��<(�M<,sJ<)Nu<%V]<$Y�<&/�<(j�<'ك<$��<$�<$><<#�<&CV<%e<#�B<$a�<%�<$�<&	<<&�<$5�<$!�<$t�<%X�<$�k<$h><#��<$�<$��<%��<'WT<$P<#�u<#ڄ<#�<$�<'e`<*�<%�<#�<$ �<$UO<%v�<%��<%��<%gj<#��<#��<#�R<#��<#�<#�U<$	g<$n�<$(s<#��<#�<#�<$�<$(<#��<#�I<#י<#۠<#�~<#�2<#ׂ<#�<#�,<$R�<$�'<$\�<#�<#��<$#<$��<$,]<$|�<#�x<#ۉ<#�<#�_<#�L<#�T<#�<#�f<#�<#��<#��<#�<#��<#ߚ<#�<#�$<#�<$#o<#��<#�<#�I<#��<#�<#�<#�<#�<#כ<#��<#��<#��<#߬<#�E<#�V<#�m<#�W<#�<#�r<#�J<#�<#�<#�L<#�q<#ײ<#�<#ۚ<#�A<#�b<#׆<#�}<#�
<#ت<#�}<#׭<#�<#�<#��<#�(<#�<<#�_<#�-<#�<#�H<#ش<#�7<#�u<#�F<#�<#��<$+�<#�<#�v<#��<#ء<#��<#��<#�3<#�<#��<#�X<#��<#�/<#��<#ލ<#�<#�<#�<#י<#��<#ׁ<#�<#�
<#�<#�K<#��<#�<#�B<#س<#�<#�
<#�s<#�
<#�3<#כ<#��<#�<$F<#׀<#��<#�<#��<#�<#�<#�<#�<#ب<#� <#��<#�<#�
<#�
<#�<#��<#׏<#ע<#�`<#�O<#�<#�w<#ׇ<#�i<#��<#ח<#��<#�
<#�
<#�<#�U<#�
<#�
<#�
<#�
<#�z<#�
<#�
<#�~<#�~<#�
<#�U<#�
<#�
<#�
<#�z<#�<#ק<#׃<#��<#�<#��<#�a<#�<#�
<#�R<#��<#�u<#א<#��<#�y<#ޡ<#��<#ד<#��<#�<#׶<#�<#��<#ۘ<#��<#�<#�<#�<#�s<#�)<#�C<#�
<#�<#��<#�<#۰<#�^<#��<#��<$<#�|<$%�<$3,<#ۊ<#�l<#׮<#�G<#�'<#׏<#��<#�
<#�l<#��<$d<#�<#�k<#��<#� <#�=<#ߥ<#��<#ނ<#ތ<#�<#��<#�<#�k<#�j<#�O<#��<#�<#�<#�q<#��<#�s<#��<#׫<#��<#�r<#�*<$ <$K�<#ٿ<#�<#�<#�<#�
<#�<#�w<#�h<#د<#�<#ׁ<#�~<#�
<#�
<#�w<#�
<#�
<#�
<#�
<#�x<#�
<#�
<#�v<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929162924201909291629242019092916292820190929162935202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@333@333@333G�O�G�O�G�O�G�O�DY  DY  DY  G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          