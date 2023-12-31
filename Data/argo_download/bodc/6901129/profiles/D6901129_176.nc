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
_FillValue                 �  V0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Z   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  [�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  c�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  kT   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  s   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  t�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  v�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  x�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �x   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �(   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �8   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �8   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �8   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �8   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �    HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �<   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  � Argo profile    3.1 1.2 19500101000000  20210225043100  20210225043100  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125463                          2C  D   APEX                            6229                            120210                          846 @��K�s�1   @��K�s�@P��1&��4D�t�j1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @�  @�  @���A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  C   C  C�fC  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D-��D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm�fDn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du�fDv,�Dv@ BffBffBe`Be`BffBffBgmBhsBe`BdZB_;B]/B`BB^5B\)B\)B\)B\)B\)B\)B\)B\)B\)B\)B[#B[#BZBYBVBR�BQ�BP�BP�BP�BP�BP�BP�BO�BQ�BQ�BT�BXBW
BXBZBZBYBXBW
BVBS�BR�BR�BR�BVBW
BYBYBXBW
BT�BO�BO�BP�BP�BP�BQ�BR�BT�BVBW
BYBZB[#B\)B]/B]/B\)B^5B`BB`BB^5BaHBcTBe`BhsBiyBjBk�Bl�Bm�Bn�Bq�Bq�Bq�Bq�Bq�Br�Br�Bt�Bw�By�By�Bz�B{�B|�B}�B~�B�B�B�%B�1B�7B�7B�7B�7B�7B�7B�=B�=B�=B�=B�DB�DB�DB�DB�JB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�DB�DB�DB�JB�JB�JB�PB�PB�VB�\B�\B�\B�\B�\B�\B�VB�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�PB�PB�PB�VB�VB�VB�VB�VB�\B�\B�VB�\B�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@/
=@/;d@/K�@.E�@+ƨ@*�!@*�\@+o@/
=@/\)@0  @#ƨ@@=q@��@G�@��@�u@%@�`@Ĝ@��@�u@��@|�@
=@�y@5?@
=q@{@~�@�?���?��?�(�?��H?�X?�ff?�V?��?ٙ�?���?���?��?�"�?��y?°!?�v�?�?�J?��h?��^?�$�?��?�{?�x�?��+?�z�?��H?��\?��?a��?I��?C��??;d?=�-?;�m?;��?<j?;dZ?9�#?9��?9�?8Q�?8Q�?6E�?4�j?2�?&��?�+?
��>��>ٙ�>��>�G�>�ff>���>�>>�->�9X>��j>���>�Q�>��j>�!>�{>�`B>�Z>���>��T>�~�>�`B>�S�>��/>��/>�A�>߾w>�S�>��T>>�!>�33>�F>�F>�33>�33>�!>�&�>�&�>�1>��/>�Z>���>߾w>޸R>��>�b>Ձ>�ƨ>�$�>��H>�->��D>��>�~�>���>�r�>�l�>��y>�`B>�A�>��->�/>�(�>�(�>��>��>�hs>�7L>��9>��>�$�>�$�>��>��>���>|�>{�m>|�>�  >��>��7>�J>��>��>���>�O�>��^>��\>}�>|�>y�#>s�F>n��>k�>k�>l�D>gl�>fff>cS�>Z�>O�;>M��>I�^>B�\>>v�>8Q�>6E�>5?}>49X>2->&�y>"��>�->�P>�>
=q>�>J=��=��=�l�=�G�=�/=�"�=��=���=��=���=\=�9X=�{=�1=��=��=��=���=��w=��-=���=�hs=�O�=�O�=�C�=�7L=�7L=�+=�+=�+=�7L=�7L=�7L=�7L=�+=�7L=�+=�+=�%=�%=}�=m�h=ix�=]/=D��=0 �=,1='�=#�
=�w=��=��=�P=��=t�=\)=\)=\)=+=o=o<�`B<�j<�9X<�1<�1<��
<���<�C�<�o<u<u<u<T��<D��<#�
<t�<o;��
;o�o�ě��t��D���e`B�u���㼴9X��j��j�ě��ě����ͼ�`B���C���w�''',1�0 Ž8Q�<j�D���D���D���H�9�L�ͽT���]/�aG��e`B�u�}󶽃o�����+��7L��C���\)��������������-��1��9X��^5��vɽ���ě��������`�����G���l�����h������������#�   �o�1'�
=q�O߾n��t�����P�����w� Ĝ�!���"��$�/�%�T�+�-V�0 ž5?}�6E��7KǾ:^5�<j�?|�?|�A�7�D���H�9�J���N��P�`�R�S�ϾV�["Ѿ^5?�dZ�gl��ixվk��l�D�p�׾p�׾q���vȴ�{�m�}�~�۾�%���\������������$ݾ��𾆧�+��1'���9������I���\)���;���`��t����������P������-��5?���R���R���w��A�������S���l����h��33��9X��?}��E���KǾ�KǾ�Q쾹�#���H���m��p����۾�|��%��J��J�����$ݾ�+��1'��1'�ȴ9�ɺ^���;�\)��t���
=����ٙ��ٙ�����ڟ�����ٙ���"Ѿ޸R��;d��;d��A���A���Ĝ��G���������������MӾ�����MӾ�����������������������������G���G���G���G���Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��A���;d111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��\@\@�\)A!G�AAG�AaG�A���A���A���A���A���AУ�A��A��B Q�BQ�BQ�BQ�B Q�B(Q�B0Q�B8Q�B@Q�BHQ�BPQ�BXQ�B`Q�BhQ�BpQ�BxQ�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B�(�B���B�(�C {C{C��C{C{C
{C{C{C{C{C{C{C{C{C{C{C {C"{C${C&{C({C*{C,{C.{C0{C2{C4{C6{C8{C:{C<{C>{C@{CB{CD{CF{CH{CJ{CL{CN{CP{CR{CT{CV{CX{CZ{C\{C^{C`{Cb{Cd{Cf{Ch{Cj{Cl{Cn{Cp{Cr{Ct{Cv{Cx{Cz{C|{C~{C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
C�
C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=C�
=D D �DD�DD�DD�DD�DD�DD�DD�DD�D	D	�D
D
�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�D D �D!D!�D"D"�D#D#�D$D$�D%D%�D&D&�D'D'�D(D(�D)D)�D*D*�D+D+�D,D,�D-D-�D-��D.�D/D/�D0D0�D1D1�D2D2�D3D3�D4D4�D5D5�D6D6�D7D7�D8D8�D9D9�D:D:�D;D;�D<D<�D=D=�D>D>�D?D?�D@D@�DADA�DBDB�DCDC�DDDD�DEDE�DFDF�DGDG�DHDH�DIDI�DJDJ�DKDK�DLDL�DMDM�DNDN�DODO�DPDP�DQDQ�DRDR�DSDS�DTDT�DUDU�DVDV�DWDW�DXDX�DYDY�DZDZ�D[D[�D\D\�D]D]�D^D^�D_D_�D`D`�DaDa�DbDb�DcDc�DdDd�DeDe�DfDf�DgDg�DhDh�DiDi�DjDj�DkDk�DlDl�DmDm��DnDn�DoDo�DpDp�DqDq�DrDr�DsDs�DtDt�DuDu��Dv1�DvEBfCBf`BfBg(Bg0Bf}Bf�Be�BeBd8Bg�Bc�BaBb�B^%B\uB\mB[�B\BB\BB\BB\;B\�B\gB[yB[@BZ�B\BYBU�BSBQ�BQ�BQ�BQcBQ}BRBR�BU1BVZBW~BX}BX(BYxB[�B[�BZ�B[!BX�BW�BUmBTWBUBS�BW�BX+BY�B\�B[BW5B[EBThBQ+BQ�BQ0BQ=BQ�BR�BU/BVMBWBY/BZBB[(B\�B]zB]�B^�Ba<Bb�Bb�BawBaIBb�Bd�Bh2Bi:Bj4Bk0BlZBm�Bn>Bq�Bq�Bq�BrBryBr�Br�BtXBweBzHBzBz�B{�B}XB~B~�B��B�GB��B�B�*B�8B�CB�8B�DB�\B�AB��B��B�OB��B�vB�bB��B�eB��B�)B��B�OB�B��B�ZB�WB�bB�KB�VB�KB�eB��B�yB�IB�UB�?B�nB�]B��B� B�MB�VB�aB�?B�VB�>B�8B��B�QB�7B�#B�:B�2B�;B�B�B�B�"B��B�B��B�lB��B��B��B��B�\B�SB��B�lB��B��B��B�tB��B��B��B��B�nB�cB�cB�sB��B��B��B��B�vB��B��B�qB�B��B��B�}B�tB�iB�iB�vB��B�yB��B��B��B�vB��B�iB�jB��B�vB�xB��B��B��B�qB�{B�{B�pB�}B�sB�uB�iB�uB�uB�yB��B�pB��B�}B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�/B��B��B��B��B��B��B��B��B�B�GB�;B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B��B��B��B��B��B��B�B��B�B�B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�nB��@/
=@/;d@/K�@.E�@+ƨ@*�!@*�\@+o@/
=@/\)@0  @#ƨ@@=q@��@G�@��@�u@%@�`@Ĝ@��@�u@��@|�@
=@�y@5?@
=q@{@~�@�?���?��?�(�?��H?�X?�ff?�V?��?ٙ�?���?���?��?�"�?��y?°!?�v�?�?�J?��h?��^?�$�?��?�{?�x�?��+?�z�?��H?��\?��?a��?I��?C��??;d?=�-?;�m?;��?<j?;dZ?9�#?9��?9�?8Q�?8Q�?6E�?4�j?2�?&��?�+?
��>��>ٙ�>��>�G�>�ff>���>�>>�->�9X>��j>���>�Q�>��j>�!>�{>�`B>�Z>���>��T>�~�>�`B>�S�>��/>��/>�A�>߾w>�S�>��T>>�!>�33>�F>�F>�33>�33>�!>�&�>�&�>�1>��/>�Z>���>߾w>޸R>��>�b>Ձ>�ƨ>�$�>��H>�->��D>��>�~�>���>�r�>�l�>��y>�`B>�A�>��->�/>�(�>�(�>��>��>�hs>�7L>��9>��>�$�>�$�>��>��>���>|�>{�m>|�>�  >��>��7>�J>��>��>���>�O�>��^>��\>}�>|�>y�#>s�F>n��>k�>k�>l�D>gl�>fff>cS�>Z�>O�;>M��>I�^>B�\>>v�>8Q�>6E�>5?}>49X>2->&�y>"��>�->�P>�>
=q>�>J=��=��=�l�=�G�=�/=�"�=��=���=��=���=\=�9X=�{=�1=��=��=��=���=��w=��-=���=�hs=�O�=�O�=�C�=�7L=�7L=�+=�+=�+=�7L=�7L=�7L=�7L=�+=�7L=�+=�+=�%=�%=}�=m�h=ix�=]/=D��=0 �=,1='�=#�
=�w=��=��=�P=��=t�=\)=\)=\)=+=o=o<�`B<�j<�9X<�1<�1<��
<���<�C�<�o<u<u<u<T��<D��<#�
<t�<o;��
;o�o�ě��t��D���e`B�u���㼴9X��j��j�ě��ě����ͼ�`B���C���w�''',1�0 Ž8Q�<j�D���D���D���H�9�L�ͽT���]/�aG��e`B�u�}󶽃o�����+��7L��C���\)��������������-��1��9X��^5��vɽ���ě��������`�����G���l�����h������������#�   �o�1'�
=q�O߾n��t�����P�����w� Ĝ�!���"��$�/�%�T�+�-V�0 ž5?}�6E��7KǾ:^5�<j�?|�?|�A�7�D���H�9�J���N��P�`�R�S�ϾV�["Ѿ^5?�dZ�gl��ixվk��l�D�p�׾p�׾q���vȴ�{�m�}�~�۾�%���\������������$ݾ��𾆧�+��1'���9������I���\)���;���`��t����������P������-��5?���R���R���w��A�������S���l����h��33��9X��?}��E���KǾ�KǾ�Q쾹�#���H���m��p����۾�|��%��J��J�����$ݾ�+��1'��1'�ȴ9�ɺ^���;�\)��t���
=����ٙ��ٙ�����ڟ�����ٙ���"Ѿ޸R��;d��;d��A���A���Ĝ��G���������������MӾ�����MӾ�����������������������������G���G���G���G���Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��Ĝ��A���;d111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#ح<#� <$R�<&mK<$b<#ڹ<#��<)�V<#�<#ئ<W��<A�U<$v.<3DY<'[<#�B<#�4<#�<#�<#�<#�
<#٫<$�<#�<#��<#�<$*�<*��<+K�<)n�<$�+<$�@<$1�<$;E<$�<$(�<$��<*�<,68<2ʾ<(ѯ<$?<$�U<%y�<%ݛ<%�<&�<+?l<%ձ<&#�<%��<%t<'E<$��<&C><$�<$}�<-x�<+/J<#�U<@��<3|�<%=�<$w�<#��<#�N<#�<#�?<#�b<#��<#��<#��<#�Q<#��<#��<#�<$ �<(J�<+x<'��<(S�<,*�<#ט<$ �<#��<#߷<#��<#�<#�v<#�D<#�<#�R<#׿<#�<#�5<$		<$h�<#܊<#�-<#�<#��<$�<#��<#ؐ<#��<$�<#��<#�<#�g<$F'<#��<#� <#�<#׋<#��<#י<#�<#�a<#��<$�<$:�<#��<#�N<#��<#܋<$<#ۿ<#��<$��<$!Y<$Ϻ<$p=<$�<#�V<#�|<#ހ<#�#<#�2<#�2<#�@<$j<#�<#��<#�!<#׬<#�d<#�<$L�<$Y�<#�z<#�e<#�><#׬<#�/<#׎<#�1<$9c<#��<#�
<#�_<#�<#�y<#�<#�&<#�&<#�<#݊<#��<$<�<#�&<#�k<#�-<#�(<#��<#ޖ<#׉<#�<#�<#�i<#��<$�<$<<#�v<#�S<#��<#�<#�<#�<#�<#�<#�D<$�<#�X<#�<#��<#�u<$�<#��<#ۭ<#߼<#��<#�6<#�<#�<#�	<#�<#�k<#ަ<#�.<#�2<#��<#�^<#�<#�<#ט<#נ<#ހ<#�<#ه<#��<#�<#�C<#מ<#��<#��<#׍<#�+<#�Y<#�|<#�
<#ׁ<#ׅ<#��<#؁<#�<#��<#ש<#�<#י<#�$<#�)<#��<#��<#�r<#�<#��<#��<#��<#��<#��<#ו<#��<#�<#ۇ<#؊<#׌<#ז<#�0<#��<#�w<#�a<#�<#��<#��<#ו<#��<#�<#�-<#�<#�5<#ׄ<#�n<#�y<#ذ<#�9<#��<#�<#�|<#��<#�<#��<#��<#ޚ<#ە<#��<#�<#ݽ<#��<#׌<#��<#ו<#�<#ޣ<#��<#�<#�c<#�]<#׃<#׌<#��<#�<#�9<#�<#۫<#�^<#׋<#��<#�<#�A<#�A<#��<#�(<#�L<#ێ<#�A<#��<#��<#��<#�<#�j<#ވ<#�><#ר<#ۺ<#�<#�?<#ߣ<#��<#�<#�k<#ޞ<#��<#��<#��<#ޮ<#�Y<#�<#��<#��<#ל<#�(<#�t<#޽<#��<#�i<#ܾ<#ݸ<#�j<#�C<#�A<#�v<#��<#�M<#�7<#��<#�<#�8<#�H<#�<#۬<#��<#�M<#�8<#�<#�}<#�{<#�5<#װ<#�P<#��<#�<#��<#�p<#�l<#�?<#�<#ۍ<#�<#�!<#�z<#��<#��<#�A<#�+<#��<#װ<#��<#�_<#�<#��<#�$<#ޝ<#ޕ<#�U<#�<#�E<#�=<#��<#ג<#��<#�6<#��<#�<#�[<#�<#�<#�<#�<#�
<#ۃ<#�?<$�<#�^<#��<#��<#ט<#�<#�#<#�J<#��<#�<$&�<$5<#�X<#�X<#�V<#�%<#ך<#�]<#�x<#�r<#�r<#޷<#�u<#�6<#�W<#�2<#��<#�<#ۿ<#�U<#�#<#ו<#��<#۪<#�<#�<$ �<#�<#��<#��<#׎<#��<#��<#�
<#�<#��<#��<#��<#ב<#� <#ך<#��<#��<#��<#ׇ<#ׇ<#��<#�
<#ص<#�<#�w<#�<#�<#�<#�~<#�v<#�
<#�y<#�~<#�v<#�
<#�{<#�~<#�~<#�~<#�~<#�~<#�x<#�<#�R<#�n;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.08                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929155423201909291554232019092915542720190929155434202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�  @�  @�  G�O�G�O�G�O�G�O�Dv@ Dv@ Dv@ G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          