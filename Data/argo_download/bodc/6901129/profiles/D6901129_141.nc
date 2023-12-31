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
resolution        ?PbM���     �  G   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  O   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    W   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    Y   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    [   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  ]   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  m   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                    u   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                    w   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                    y   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  {   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �    	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �\   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �\   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �\   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �\   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �    HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �D   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �`   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �$Argo profile    3.1 1.2 19500101000000  20210225040608  20210225040608  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125429                          2C  D   APEX                            6229                            120210                          846 @צ�3�a�1   @צ�3�a�@Q���l��4�1&�x�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� DqfDq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� DfD�fBl�BjBgmBjBk�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bn�Bo�Bn�Bo�Bo�Bp�Bo�Bp�Bp�Bq�Bp�Bq�Bq�Bq�Bp�Bo�Bm�Bl�Bk�Bk�BjBjBiyBiyBiyBhsBhsBgmBgmBgmBffBe`Be`BdZBcTBcTBcTBcTBdZBdZBcTBdZBdZBdZBdZBe`Be`Be`BffBffBffBffBgmBhsBjBk�Bl�Bm�Bn�Bp�Bp�Bq�Bq�Bq�Bq�Bq�Br�Br�Bs�Bt�Bu�Bu�Bu�Bv�Bx�By�Bz�Bz�B{�B{�B|�B~�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�1B�7B�7B�=B�DB�DB�JB�JB�PB�\B�\B�VB�\B�\B�\B�\B�bB�bB�bB�hB�hB�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?>�R??;d?8Q�?X?�?�
>�|�>�v�>�j>��H>��H>��#>���>�K�>�ȴ>�E�>�E�>�>��j>��j>�9X>�F>�F>�F>�9X>�33>�33>�F>��j>�9X>�33>�F>�F>�9X>�?}>�?}>�>�>�E�>�ȴ>�ȴ>�ȴ>�K�>�K�>���>���>���>���>��#>�dZ>��m>��m>�dZ>��m>��m>��>��>�v�>��>���? Ĝ?M�?J?�7?��?Z?��?��?��?�T?`B?��>�!>�Ĝ>�;d>Ձ>���>�=q>�J>���>���>��/>��->��>�I�>y�#>k�>R�>Kƨ>B�\>-V>hs>J>o>%>   =��=��m=��m=��#=��#=��#=��m=��>   >%>J>�>+>I�>z�>��>�R>+>2->;dZ>C��>C��>E��>F��>F��>F��>I�^>Kƨ>N�>S��>Y�>Y�>Z�>_;d>ix�>n��>r�!>vȴ>w��>{�m>}�>���>�7L>�C�>�C�>�ƨ>�O�>�C�>�=q>��^>�7L>�1'>��>���>�$�>���>�o>�o>��>��>���>���>�+>���>�J>�  >~��>|�>y�#>x��>w��>s�F>|�>�%>���>��>���>��^>���>���>��>��P>��>��;>���>�I�>���>��\>���>��9>��9>��9>��^>�hs>�>�
=>��>��`>�z�>��u>��u>�b>�b>�>�t�>�n�>�hs>�bN>���>�I�>�I�>�7L>�+>�$�>���>��\>��>z�H>u>s�F>q��>m�h>j~�>gl�>bM�>_;d>Z�>T��>R�>Q�>P�`>N�>J��>G�>F��>E��>D��>=p�>49X>1&�>0 �>-V>)��>'�>'�>&�y>"��> Ĝ> Ĝ> Ĝ>�w>�R>�->��>�P>�+>�>hs>\)>O�>I�>I�>C�>
=q>	7L>+>�>   =��#=��=�h=�x�=�`B=�;d=��=���=���=��=���=��=�9X=�-=�{=�t�=�O�=�O�=�O�=�7L=}�=m�h=aG�=Y�=T��=P�`=,1=�w=��=t�=t�<�<���<ě�<ě�<�9X<��
<T��<o;ě�:�o��o�o�D���D�����
�ě���`B�o�T����o��C���t���9X���������C��t��,1�8Q�<j�<j�@��H�9�Y��aG��e`B�q����C���\)������㽛�㽝�-���T����置{�� Ž�-��E���vɽ�����ͽ�
=��"ѽ�;d��S���l�������#���#�   �J�J�J���O߾hs�t���+��R��w� Ĝ�!���!���(�þ-V�333�?|�@��B�\�F��G��F��G��H�9�I�^�Kƨ�O�;�R�T���W
=�["Ѿ_;d�bMӾe`B�hr��ixվl�D�m�h�n���r�!�u�x���}󶾁�7���7���7���\���\��o��o��o��������+��I����;���`��녾�񪾔zᾘ�u��������;d��A���Ĝ��MӾ�`B��ff���y���þ��þ��羬�D��{�������׾�-���!���F���j��KǾ��#���m��j��푾�p���%�\��o�Ǯ���;��;���Ձ�և+��
=�ٙ�����ܬ�ݲ-��5?�޸R��;d��Ĝ��S���Z���/���/���T���y���þ�������~�����D��D��{���׾�-��E����#��j��푾�vɿ A�� �� Ĝ� Ĝ�%�G��G��G��G��G��G��G��G��G��%�%�%� Ĝ� Ĝ� �� �� �� �� �� ���|�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @5�@{�@�@�A�HA>�HA^�HA~�HA�p�A�p�A�p�A�p�A�p�A�p�A�p�A�p�B�RB�RB�RB�RB'�RB/�RB7�RB?�RBG�RBO�RBW�RB_�RBg�RBo�RBw�RB�RB��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)B��)C�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C3�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C��
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
D {�D ��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D	{�D	��D
{�D
��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D{�D��D {�D ��D!{�D!��D"{�D"��D#{�D#��D${�D$��D%{�D%��D&{�D&��D'{�D'��D({�D(��D){�D)��D*{�D*��D+{�D+��D,{�D,��D-{�D-��D.{�D.��D/{�D/��D0{�D0��D1{�D1��D2{�D2��D3{�D3��D4{�D4��D5{�D5��D6{�D6��D7{�D7��D8{�D8��D9{�D9��D:{�D:��D;{�D;��D<{�D<��D={�D=��D>{�D>��D?{�D?��D@{�D@��DA{�DA��DB{�DB��DC{�DC��DD{�DD��DE{�DE��DF{�DF��DG{�DG��DH{�DH��DI{�DI��DJ{�DJ��DK{�DK��DL{�DL��DM{�DM��DN{�DN��DO{�DO��DP{�DP��DQ{�DQ��DR{�DR��DS{�DS��DT{�DT��DU{�DU��DV{�DV��DW{�DW��DX{�DX��DY{�DY��DZ{�DZ��D[{�D[��D\{�D\��D]{�D]��D^{�D^��D_{�D_��D`{�D`��Da{�Da��Db{�Db��Dc{�Dc��Dd{�Dd��De{�De��Df{�Df��Dg{�Dg��Dh{�Dh��Di{�Di��Dj{�Dj��Dk{�Dk��Dl{�Dl��Dm{�Dm��Dn{�Dn��Do{�Do��Dp{�Dq�Dq{�Dq��Dr{�Dr��Ds{�Ds��Dt{�Dt��Du{�Du��Dv{�Dv��Dw{�Dw��Dx{�Dx��Dy{�Dy��Dz{�Dz��D{{�D{��D|{�D|��D}{�D}��D~{�D�D��Bl�Bk�Bm4Bl�Bm�BmWBl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�BlBltBl�Bl�Bl�Bl�BlBluBm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�BmyBmxBnrBn�Bn�Bn�Bn�Bn�BnmBo�Bn�Bo�Bn�Bo^BoXBp�Bo�BpDBp�Bq�Bp�BqyBq�Bq�Bq
Bq�Bo@Bl�BlmBk�BkiBkIBjuBj]BjcBi&Bh�Bh�Bh�Bh&Bg�Be�Be�BeaBd�BdBcIBcnBdeBdeBccBdZBdfBd[Bd\BeTBeTBeWBfWBfYBfABfNBg.BhBj)BkZBk�Bm6Bn)BpABp�Bq�Bq�Bq�Bq�Bq�Br�Br�BsxBt�Bu�Bu�Bu�BvQBx�By�Bz�Bz�B{�B{�B|jB~�B��B�B��B��B�AB�.B�#B�&B�2B�(B�HB�B�CB�DB�"B�B�$B�B��B� B�hB�jB�\B�:B�DB�PB�8B�9B�YB��B��B��B�.B�2B��B�B�qB�^B�%B��B��B��B��B�B��B�B�0B�cB�fB�KB��B�B�eB��B��B�0B�)B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B�`B��B��B��B��B�B��B��B��B��B��B�3B��B��B��B��B�B��B��B��B��B��B�!B�B��B�B��B��B��B��B��B��B��B��B�B��B��B��B�B�B�B��B��B��B�B��B��B��B��B��B�B��B��B��B�=B��B��B��B��B��B�B��B��B��B��B��B��B�B��B�B�B��B��B��B��B�B�B��B��B��B��B��B��B�3B�B��B��B�1B��B��B��B��B�&B�B�B�aB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�B��B��B��B��B��B��B��B��B��B��B�BB�B��B��B��B��B�&B��B�B�B��B��B��B�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�0B�;B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�#B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?>�R??;d?8Q�?X?�?�
>�|�>�v�>�j>��H>��H>��#>���>�K�>�ȴ>�E�>�E�>�>��j>��j>�9X>�F>�F>�F>�9X>�33>�33>�F>��j>�9X>�33>�F>�F>�9X>�?}>�?}>�>�>�E�>�ȴ>�ȴ>�ȴ>�K�>�K�>���>���>���>���>��#>�dZ>��m>��m>�dZ>��m>��m>��>��>�v�>��>���? Ĝ?M�?J?�7?��?Z?��?��?��?�T?`B?��>�!>�Ĝ>�;d>Ձ>���>�=q>�J>���>���>��/>��->��>�I�>y�#>k�>R�>Kƨ>B�\>-V>hs>J>o>%>   =��=��m=��m=��#=��#=��#=��m=��>   >%>J>�>+>I�>z�>��>�R>+>2->;dZ>C��>C��>E��>F��>F��>F��>I�^>Kƨ>N�>S��>Y�>Y�>Z�>_;d>ix�>n��>r�!>vȴ>w��>{�m>}�>���>�7L>�C�>�C�>�ƨ>�O�>�C�>�=q>��^>�7L>�1'>��>���>�$�>���>�o>�o>��>��>���>���>�+>���>�J>�  >~��>|�>y�#>x��>w��>s�F>|�>�%>���>��>���>��^>���>���>��>��P>��>��;>���>�I�>���>��\>���>��9>��9>��9>��^>�hs>�>�
=>��>��`>�z�>��u>��u>�b>�b>�>�t�>�n�>�hs>�bN>���>�I�>�I�>�7L>�+>�$�>���>��\>��>z�H>u>s�F>q��>m�h>j~�>gl�>bM�>_;d>Z�>T��>R�>Q�>P�`>N�>J��>G�>F��>E��>D��>=p�>49X>1&�>0 �>-V>)��>'�>'�>&�y>"��> Ĝ> Ĝ> Ĝ>�w>�R>�->��>�P>�+>�>hs>\)>O�>I�>I�>C�>
=q>	7L>+>�>   =��#=��=�h=�x�=�`B=�;d=��=���=���=��=���=��=�9X=�-=�{=�t�=�O�=�O�=�O�=�7L=}�=m�h=aG�=Y�=T��=P�`=,1=�w=��=t�=t�<�<���<ě�<ě�<�9X<��
<T��<o;ě�:�o��o�o�D���D�����
�ě���`B�o�T����o��C���t���9X���������C��t��,1�8Q�<j�<j�@��H�9�Y��aG��e`B�q����C���\)������㽛�㽝�-���T����置{�� Ž�-��E���vɽ�����ͽ�
=��"ѽ�;d��S���l�������#���#�   �J�J�J���O߾hs�t���+��R��w� Ĝ�!���!���(�þ-V�333�?|�@��B�\�F��G��F��G��H�9�I�^�Kƨ�O�;�R�T���W
=�["Ѿ_;d�bMӾe`B�hr��ixվl�D�m�h�n���r�!�u�x���}󶾁�7���7���7���\���\��o��o��o��������+��I����;���`��녾�񪾔zᾘ�u��������;d��A���Ĝ��MӾ�`B��ff���y���þ��þ��羬�D��{�������׾�-���!���F���j��KǾ��#���m��j��푾�p���%�\��o�Ǯ���;��;���Ձ�և+��
=�ٙ�����ܬ�ݲ-��5?�޸R��;d��Ĝ��S���Z���/���/���T���y���þ�������~�����D��D��{���׾�-��E����#��j��푾�vɿ A�� �� Ĝ� Ĝ�%�G��G��G��G��G��G��G��G��G��%�%�%� Ĝ� Ĝ� �� �� �� �� �� ���|�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<%_�<;��<&�<'�<$G@<#��<#ۀ<#�	<#�I<#ר<#�c<#�<#�<#�<#�a<#�<#׈<#�V<#�<#�<#�g<#�s<#؆<#ׅ<#�n<#��<#ڬ<#�<#ׇ<#ؒ<#�|<#��<#ڕ<#א<#ؤ<#�{<#ض<#ر<#�q<#�s<#أ<#�v<#ت<#�s<#�{<#��<#�$<#ގ<#��<#�j<#�
<#إ<#ד<#�><#��<#ء<#�<#��<#�<#�<#�<#׭<#�~<#��<#ر<#�J<#�j<#ش<#׼<#�<&{<%�	<#ܼ<$kI<#�<$m�<$Eg<$�G<$e<$m�<$,T<#��<%!�<%]�<$2�<$�%<#�<#�<$�|<%<$3�<#ؔ<#״<#�<#�<#�<#�v<#�<#�d<#�R<#��<#��<#�^<#�+<#��<#�$<#��<#��<#��<#�0<#��<$!�<#��<$�<#��<#ׄ<#�k<#ؼ<#�k<#׃<#�<#��<#�4<#�<#��<#��<#��<#��<$�<#��<#�<#�<#��<#��<#�~<$<$n<#�<#�^<#�5<#��<#��<#��<#�<#�<#ם<#�<#��<#��<#��<#� <#�><#�L<#�}<#��<#�g<#؋<#އ<#�<#�{<#�<#צ<#�<#�<#�<#��<$<#��<#�z<#�<#ٷ<#��<#�<#�R<$��<#�<$?<#ܤ<#�<#�><$-�<#ڒ<#��<#�w<#�W<#ד<#�<$F�<$k<#�Z<#��<#ׅ<#�<#��<#ר<#�<#�F<#�&<#�!<#ר<#��<#�T<#�l<#�E<#�H<#�y<#��<#��<#�<#�h<#�~<#�=<#��<#��<#ׄ<#�v<#�I<#�`<#��<#و<#�<#�7<#ב<#�<#�<#׶<#۩<#�#<#�<#�<#�<#�<#��<#ه<#�<#�<<#��<#ן<#�S<#�<#�q<#צ<#�Y<#�`<#�<#�<#�<#�<#�<#�<#�<#�E<#׷<#ן<#�<#׍<#�<#�<#�<#׫<#�><#۝<#�P<#ہ<#׾<#ת<#״<#�7<#ٝ<#�g<#�U<#�<#�_<#�A<#��<#�<#�"<$6<#��<#�W<#�W<#׻<#��<#��<#�4<#פ<#�<#� <#�'<#�<#�<#א<#�5<#�<#۫<#ש<#�K<#ן<#��<#�#<#�#<#ר<#޿<#׭<#�<#�<#�S<#ס<#�<#�<#�<#��<#�+<#�<#�<#ۙ<#��<#ہ<#�<#�<#��<#�	<#��<#�<#�|<#�<#׹<#�f<#״<#�<#٦<#��<#�<#�5<#�<#�O<#�<#��<#�<#�<#�h<#�=<#�
<#��<#��<#�!<#��<#��<#׸<#׭<#׮<#��<#�4<#�}<#�9<#�<#ם<#�_<#�O<#�l<#��<#ۡ<#��<#ِ<#�<#�<#�<#�<#�3<#��<#�<#��<$�<#�<#ף<#۶<#�<#�w<#�<#�<#�<#׸<#ۍ<#�;<#ײ<#��<#ۣ<#۾<#�H<#�F<#�#<#�<#�<#�<#�<#�p<#�V<#�o<#�/<#޲<#�F<#�W<#ד<#�Q<#�<#�_<#�V<#ע<#��<#ۯ<#��<#�<#��<#ׯ<#׻<#٩<#��<#��<#�'<#��<#ױ<#�<#�\<#��<#צ<#�<#�'<#�?<#��<#��<#�l<#�8<#��<#�&<#�<#ק<#ה<#ߢ<#�F<#ۥ<#�(<#�
<#�<#�<#��<#�<#��<#�s<#�<#��<#��<#ק<#�'<#�D<#�
<#�*<#ײ<#�<#�<#�<#�r<#��<#��<#�<#�S<#ש<#��<#�v<#ת<#�\<#�<#�<#��<#�G<#�7<#�<<#��<#٣<#�<#�}<#�<#�<#�Q<#ۉ<#�<#�<#�W<#�<#�<#�\<#�c<#�c<#�c<#�c<#�c<#�c<#�k<#ؘ<#�k<#�j<#؝<#�m<#ؔ<#�j<#�c<#�c<#�c<#�v<#ݕ<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.07                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929145904201909291459042019092914590820190929145914202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@9��@9��@9��G�O�G�O�G�O�G�O�D�fD�fD�fG�O�G�O�G���G���G���G���G���G���G���                                6389758         6389758         131072                                          