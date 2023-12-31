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
resolution        ?PbM���     �  Lx   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  S$   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V|   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X(   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l,   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  m�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  q0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  w�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ~�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �4   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �<   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �X   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �|   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �$   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �@   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �\Argo profile    3.1 1.2 19500101000000  20210225044003  20210225044003  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125486                          2C  D   APEX                            6229                            120210                          846 @�6�y�@1   @�6�y�@@P�n��P�58Q��1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�ff@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  A�33B  B  B  B   B(  B0  B8  B@  BH  BPffBX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D3��D4y�D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DR��DS� DT  DT� DU  DU�fDU�3B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�dB	��B	��B	�mB
�B
1'B
I�B
]/B
l�B
�B
��B
�^B
�B
��BuB �B'�B/B:^B8RB<jB?}BD�BI�BG�BJ�BO�BT�BW
BYBYB^5BaHBbNBe`BffBffBffBffBffBhsBl�Bo�Bp�Bp�Bp�Bo�Bp�Bq�Bp�Bp�Bq�Bp�Bm�BjBiyBhsBiyBjBiyBiyBk�Bq�Bq�Bp�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bo�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bm�Bm�Bn�Bn�Bm�Bn�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bo�Bp�Bp�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Bq�Bq�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bs�Bs�Bs�Bs�Bs�Bs�Bt�Bt�Bt�Bu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�Bz�Bz�Bz�B{�B{�B{�B{�B|�B|�B|�B}�B}�B~�B~�B~�B� B� B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�1B�1B�1B�1B�1B�7B�DB�DB�JB�JB�PB�PB�VB�\B�\B�\B�\B�\B�bB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�oB�uB�oB�hB�hB�hB�oB�oB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@  �@  �@ b@  �@ 1'@ bN@ �@ Q�@  �@   @ b@  �@  �@  �@ b@ A�@ r�@ 1'@ 1'@ 1'@ A�@ �`@ ��@ �9@ �9@ �9@ �9@ �9@ ��@ �9@ Ĝ@ Ĝ@ Ĝ@ Ĝ@ A�@�@r�@Q�@�@?}@bN@�@!%@$�@&��@&�+@*�@4Z@4��@)&�@"��@ �u@9X@�@Ĝ@	�@?�1'?�R?�O�?�Z?�ƨ?�33?θR?ɺ^?�?}?�O�?��y?��
?�A�?�I�?���?���?��?��?�&�?��`?��?�p�?�"�?���?�l�?� �?�~�?�r�?���?���?}p�?mV?T9X?F$�??;d?9X?49X?2-?-V?(��?6E�?6?/\)?$�?��?|�?;d?v�?(�?n�?�;?�h?��>�>�~�>�x�>���>�x�>�ff>ݲ->��>��>���>���>�^5>�K�>�->���>�x�>��y>��w>�"�>�
=>���>�n�>�V>�$�>�J>��>x��>ix�>Q�>L��>G�>D��>@�>8Q�>7K�>6E�>,1>��>\)>C�>
=q>	7L>�=��m==�;d=�/=��=���=��`=���=ě�=ȴ9=���=ě�=��=�-=�1=��T=��
=���=���=��w=��w=��w=���=��=�hs=�\)=�C�=�o=y�#=m�h=e`B=]/=D��=8Q�=49X=0 �=0 �=,1='�=#�
=�w=o<�<�h<�h<�h<�h<�h<�`B<�`B<�`B<�/<�/<�/<ě�<ě�<�j<�j<�j<�j<�j<�9X<��
<�1<�9X<ě�<�`B<�`B<�<�=o=\)=\)=@�=D��=H�9=L��=Y�=T��=@�=8Q�=0 �=<j=@�=H�9=H�9=@�=0 �=,1=�w=t�=+=o<�`B<���<���<���<���<�9X<�o<T��;ě�:�o    �ě���o;D��<o;ě����
�o<o<#�
<#�
<e`B<�1<���=+=+<��<�<�<�h<�<��<��<��<�h<�h<�`B<�9X<���<�t�<u;��
<t�<o�o��9X��/�ě���9X��1��t��T���49X�49X��o��C���t����
�ě���/��h�o�+�+�C���w�#�
�,1�<j�P�`�aG��e`B�m�h�m�h�q���y�#��o�y�#�y�#�}󶽁%��7L��C���hs���P���㽝�-���
���T��1��1�� Ž�9X��j���������������
=��������������/��l���F���ٽ��m���m�o�$ݾ	7L�	7L�	7L�	7L�I��V�hs�hs�t���P��-� Ĝ�'-V�1&�2-�333�333�5?}�8Q�8Q�=p��@��A�7�C���E�˾G��E�˾E�˾E�˾I�^�L�;M��R�Xb�aG��e`B�fff�j~��l�D�l�D�l�D�m�h�p�׾vȴ�}󶾀  ���\���^���Ͼ�����w��S����xվ��羬1��{�������������������������������������������1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @xQ�@�@�A�HA:�HAZ�HAz�HA�p�A�p�A�p�A�p�A�p�A�p�A�p�A���B�RB�RB�RB�RB&�RB.�RB6�RB>�RBF�RBO�BV�RB^�RBf�RBn�RBv�RB~�RB�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)B�\)C�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C3�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C��
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
C��=C��
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
D k�D �Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�D	k�D	�D
k�D
�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�Dk�D�D k�D �D!k�D!�D"k�D"�D#k�D#�D$k�D$�D%k�D%�D&k�D&�D'k�D'�D(k�D(�D)k�D)�D*k�D*�D+k�D+�D,k�D,�D-k�D-�D.k�D.�D/k�D/�D0k�D0�D1k�D1�D2k�D2�D3k�D3�D4eD4�D5k�D5�D6k�D6�D7k�D7�D8k�D8�D9k�D9�D:k�D:�D;k�D;�D<k�D<�D=k�D=�D>k�D>�D?k�D?�D@k�D@�DAk�DA�DBk�DB�DCk�DC�DDk�DD�DEk�DE�DFk�DF�DGk�DG�DHk�DH�DIk�DI�DJk�DJ�DKk�DK�DLk�DL�DMk�DM�DNk�DN�DOk�DO�DPk�DP�DQk�DQ�DRk�DR�DSk�DS�DTk�DT�DUq�DU��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	��B	��B	��B	�pB	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�lB	��B	�TB	�B	լB	�lB
$B
.B
F�B
ZTB
k B
�JB
�B
�YB
��B,BB"�B+3B4�B<�B=�B@�BD�BHBJ[BK BM�BS	BV�BX�BZ�B\B`�BbsBc�Bf�BhBgvBg*Bf�Bf�Bh�Bl�Bp�Bq�BqvBqMBrVBr�Br}BrBq�Br�Bs�Br*Bm*Bj�Bi�BjoBj�BjqBj)Bi=Bq�Br�Br�Bp�Bo�Bp�Bp�Bq BrgBq*BpBpzBp?Bo�Bn�Bn�Bn�Bn�BngBn�Bn�Bn�Bn�BoBo�BpBp�BqBp�BqQBqBpBp�Bp�BpBpaBpBo�BpBpYBo�Bn�Bn�Bn�Bn�Bo�Bo�Bo�BoBoiBo"Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Bq�Bq�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bs�Bs�Bs�Bs�Bs�Bs�Bt�Bt�BuBu�Bu�Bv�Bv�Bv�Bv�Bv�Bw�Bx%Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�Bz�Bz�Bz�B{�B{�B{�B{�B|�B|�B|�B}�B}�B~�B~�B~�B�BsB��B� B�B��B�$B�UB�4B�0B��B�B�
B�$B�?B�VB�4B�JB�KB�JB�4B�VB�JB�%B�B�2B�LB�oB�MB�{B�cB�5B�lB��B��B��B�MB��B�B��B�%B�BB�B��B�B��B�VB�rB�iB�]B�hB�RB�SB�bB�cB�{B�gB�wB��B��B�}B��B��B�CB��B�2B�B��B�CB�TB�dB�MB�9B�cB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�kB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B�SB��B�>B�0B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@  �@  �@ b@  �@ 1'@ bN@ �@ Q�@  �@   @ b@  �@  �@  �@ b@ A�@ r�@ 1'@ 1'@ 1'@ A�@ �`@ ��@ �9@ �9@ �9@ �9@ �9@ ��@ �9@ Ĝ@ Ĝ@ Ĝ@ Ĝ@ A�@�@r�@Q�@�@?}@bN@�@!%@$�@&��@&�+@*�@4Z@4��@)&�@"��@ �u@9X@�@Ĝ@	�@?�1'?�R?�O�?�Z?�ƨ?�33?θR?ɺ^?�?}?�O�?��y?��
?�A�?�I�?���?���?��?��?�&�?��`?��?�p�?�"�?���?�l�?� �?�~�?�r�?���?���?}p�?mV?T9X?F$�??;d?9X?49X?2-?-V?(��?6E�?6?/\)?$�?��?|�?;d?v�?(�?n�?�;?�h?��>�>�~�>�x�>���>�x�>�ff>ݲ->��>��>���>���>�^5>�K�>�->���>�x�>��y>��w>�"�>�
=>���>�n�>�V>�$�>�J>��>x��>ix�>Q�>L��>G�>D��>@�>8Q�>7K�>6E�>,1>��>\)>C�>
=q>	7L>�=��m==�;d=�/=��=���=��`=���=ě�=ȴ9=���=ě�=��=�-=�1=��T=��
=���=���=��w=��w=��w=���=��=�hs=�\)=�C�=�o=y�#=m�h=e`B=]/=D��=8Q�=49X=0 �=0 �=,1='�=#�
=�w=o<�<�h<�h<�h<�h<�h<�`B<�`B<�`B<�/<�/<�/<ě�<ě�<�j<�j<�j<�j<�j<�9X<��
<�1<�9X<ě�<�`B<�`B<�<�=o=\)=\)=@�=D��=H�9=L��=Y�=T��=@�=8Q�=0 �=<j=@�=H�9=H�9=@�=0 �=,1=�w=t�=+=o<�`B<���<���<���<���<�9X<�o<T��;ě�:�o    �ě���o;D��<o;ě����
�o<o<#�
<#�
<e`B<�1<���=+=+<��<�<�<�h<�<��<��<��<�h<�h<�`B<�9X<���<�t�<u;��
<t�<o�o��9X��/�ě���9X��1��t��T���49X�49X��o��C���t����
�ě���/��h�o�+�+�C���w�#�
�,1�<j�P�`�aG��e`B�m�h�m�h�q���y�#��o�y�#�y�#�}󶽁%��7L��C���hs���P���㽝�-���
���T��1��1�� Ž�9X��j���������������
=��������������/��l���F���ٽ��m���m�o�$ݾ	7L�	7L�	7L�	7L�I��V�hs�hs�t���P��-� Ĝ�'-V�1&�2-�333�333�5?}�8Q�8Q�=p��@��A�7�C���E�˾G��E�˾E�˾E�˾I�^�L�;M��R�Xb�aG��e`B�fff�j~��l�D�l�D�l�D�m�h�p�׾vȴ�}󶾀  ���\���^���Ͼ�����w��S����xվ��羬1��{�������������������������������������������1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#��<#ڔ<#��<#��<#�"<#��<#�g<#�~<#��<#�<#�<#�<#��<#�7<#��<#��<#�<#��<#�<#�<$$<#ۇ<#��<#��<#��<#��<#��<#�b<#�<#߯<#�<#��<#��<#�T<%�F<.O <#��<(��<#�<)�
<,V�<,�<*�\<%�0<#��<*�<Ho2<$"<O��<2,�<%��<+B<;/�<(8�<9��<24�<5��<,�<#�.<+M�<*�<*��<%��<&d<%�P<)�<'z�<$��<$�.<%q<%~B<$m^<$<#�0<#ٽ<#�<#�<$~�<$1�<$%<$�<(�*<&�&<$&_<%	�<$��<$��<*1<2Q�<(��<$֯<$��<$K�<#�"<$M/<$<(��<#��<$�g<&��<$Sv<#�<#��<#�R<#�]<%�"<#�,<#�<&�<%{J<$l�<#�;<#ۘ<#�<#��<$'�<$Zm<#�N<#ן<$�<#��<#�5<#�<#�@<#��<#��<$,<#�<#�<#�
<#נ<#�,<$�<#��<#�Q<#�v<$�<$x�<#�<#�r<#�u<#�<#ܥ<#�!<#��<#�e<$"�<#��<#�<#۱<#�<#�o<#��<#�<#ݗ<#ګ<#�<#��<#�A<#ד<#׹<#��<#�'<#��<#؜<#ڧ<#׉<#��<#�;<#ۆ<#ޢ<#ې<#ް<#ވ<#�<#א<#�U<#�3<#��<#�<#׃<#�_<#�_<#ظ<#��<#ױ<#ڔ<#�"<#ޚ<#�o<#�A<#ڻ<#۞<#��<#؆<#�<#ީ<#��<#��<#ޤ<#ۗ<#ެ<#ޢ<#�<#߇<#�]<#ת<#�w<#��<#߈<#��<#޾<#�<<#�^<#�:<#�u<#�<#�<#�<#�
<#��<#��<#�)<#�u<#��<$G`<#��<#�<#�<#�<#۬<#�V<#��<#�G<#�<#�<#�e<#�4<#��<#�<#��<#׎<#�<#׏<#��<#�<#ד<#��<#�<#�G<#�t<#��<#�V<#��<#�<#ڔ<#�[<#��<#��<#�p<#؇<#��<#��<$+[<#��<#�2<#�<$�<$ �<$5<#޳<#ك<#�B<#�a<#�v<#��<#��<#޲<#�R<#��<#��<#��<#��<#׊<#��<#�V<#��<#��<#�y<$<$�<#׃<#�&<#�<#�K<#��<#��<#��<#��<#�O<#��<#�<#��<#�
<#�x<#��<#׌<#�#<#�B<#��<#�t<#�I<#��<#�
<#׊<#�
<#�O<#و<#�5<#�<#��<#׶<#�^<#މ<#�9<#��<#�<#��<#�u<#�}<#��<#��<#׋<#��<#מ<#��<#��<#ش<#�
<#؄<#�z<#�
<#�<#�-<#ڣ<#�<#�x<#�2<#ؿ<#׌<#��<#ع<#��<#ݒ<#�m<#�s<#׌<#�5<#�h<#� <#׃<#ز<#ד<#��<#ح<#�<#�I<#�k<#�M<#כ<#�
<#��<#�<<#�,<#��<#ב<#�u<#ג<#ש<#�Q<#��<#��<#�<#�<#ߖ<#��<#�
<#�{<#ڐ<#�<#��<#��<#�<#ڎ<#�<#��<#�N<#�E<#��<#�a<#�A<#�%<#�J<#׫<$
\<$P�<#��<#��<#�<<#�<#׬<#�<#�
<#�<#�D<#�7<#�I<#�j<#�h<#�C<#�i<#�K<#�g<#�f<#�f<#�%<#�e;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.32                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929162549201909291625492019092916255320190929162601202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�DU�3DU�3DU�3G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          