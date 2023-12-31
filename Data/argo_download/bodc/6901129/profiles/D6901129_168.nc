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
resolution        ?PbM���     �  NP   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  U�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  W�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Y�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  [�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  c8   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  j�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  rh   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  tP   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  v8   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  x    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �P   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �H   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �H   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �H   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �H   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �0   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �L   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210225042301  20210225042301  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125456                          2C  D   APEX                            6229                            120210                          846 @��s�`1   @��s�`@P���-V�4�$�/�1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @333@�  @�  A��AffA@  Aa��A�  A�  A�  A���A�  A�  A�  A�  B   B  B��B��B   B(  B0  B8  B@  BH  BP  BW��B_��Bg��Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
�C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cm�fCp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  Dy�D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr�fDr�fB1B1B1B1B1B+B+B1B+B+B1B+B+B+B1B1B1B1B1B1B+B+B+B1B1B+B1B	7B	7B	7B1B1B	7B
=B
=B	7B	7B
=B	7B	7B	7B	7B	7B	7B
=B	7BJBbBhBoBuB{B{B{B�B�B�B�B�B�B�B�B�B�B�B �B5?B=qBC�BH�BR�Be`Bv�By�Bz�B~�B~�B}�B}�B|�B}�B~�B~�B}�B~�B~�B�B�B�B�B�%B�+B�+B�+B�1B�1B�7B�7B�7B�7B�7B�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�1B�1B�1B�+B�1B�1B�1B�1B�1B�7B�7B�=B�=B�DB�JB�JB�JB�PB�PB�\B�bB�hB�uB�{B�{B�{B�{B�{B�{B�{B�{B�uB�oB�hB�\B�PB�JB�JB�JB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�PB�VB�VB�PB�PB�PB�VB�VB�VB�VB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�bB�bB�bB�bB�bB�bB�hB�bB�hB�hB�hB�hB�hB�hB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?���?���?�?�?}?���?��y?�ȴ?���?�`B?�`B?�`B?���?�?�`B?��/?���?�z�?���?�9X?��?�?���?�`B?��?�9X?���?���?��?��?��?�33?�33?�n�?�-?�-?��\?��\?�M�?��!?���?���?��\?���?��\?�J?�J?�Ĝ?�\)?���?��?�v�?�V?�V?�V?�V?���?��-?��-?���?�O�?�V?��?�1?��m?�ƨ?���?�`B?���?���?� �?�I�?�Z?wK�?q�?lI�?c�
?bJ?a��?a��?]p�?R-?Kƨ?E��?>��?2-?$Z?|�?p�?�?�-?��?|�? �? Ĝ?!G�? Ĝ? Ĝ? A�?�R?��?��?�?��?
��?�9?�?�?$�?�
?��?G�>�p�>�{>�r�>�G�>�(�>ڟ�>�b>�>Ձ>��>�bN>��>��>��;>��>���>Ձ>�>�>�5?>�M�>�`B>�h>�33>�F>�F>�33>�!>�->�&�>��>�->�1>�M�>��>��>�9X>�V>�>��y>��
>���>�M�>�Ĝ>�;d>�5?>�(�>�"�>���>���>��>��>���>��>��`>�\)>��>���>�7L>�$�>��>���>��>��>���>��\>�J>}�>r�!>gl�>gl�>gl�>hr�>gl�>fff>fff>fff>gl�>gl�>gl�>gl�>fff>fff>e`B>_;d>["�>V>T��>S��>S��>O�;>H�9>?|�>;dZ>5?}>333>2->/�>)��>"��>�R>��>�>�->��>��>�u>�P>�P>�+>z�>n�>\)>
=q>+>J=���=��=��=���=���=���=ě�=ě�=\=\=\=��=��=�1=��
=���=��=��=��=�t�=�t�=�t�=�t�=�\)=�O�=�O�=�\)=�C�=�+=�+=�o=}�=y�#=q��=e`B=L��=D��=@�=8Q�=0 �=�w=t�=o=o=o=\)=\)=+=t�=t�=\)<�<�<�<�<�/<�o<T��<T��<T��<D��<49X<49X<t�<o;ě�;ě�;D��    ��o��o�o�D���D����o�ě��#�
�D���e`B��C����
��9X��1��9X��9X��j���ͼ�/��`B��`B��h���C���P���0 ŽH�9�P�`�T���]/�q���}󶽇+��t����P���w���T���T���T���罰 Ž�E���vɽ������ě����ͽ�
=��S���xս�F�%�J�o���	7L�
=q�I��I��V�z������P��u����������R��w��w��w��w� Ĝ� Ĝ�!���"��$�/�&�y�(�þ+�0 ž2-�333�8Q�8Q�8Q�:^5�>vɾ@��?|�B�\�C���C���E�˾G��M��Q녾T���Xb�\(��]/�]/�_;d�cS��dZ�fff�hr��hr��ixվixվj~��q���x���z�H��  ���7���\��o�����������˾���7L��=q��O߾�V���;������b�������㾟;d���徤�/��`B��ff���y��r���xվ�1�����������!��9X��?}������^5��j��p���󶾾�۾����J����Ƨ��7L���;��녾����և+�և+��
=��
=��
=��
=��
=��
=��
=�և+��
=��
=��
=��
=��
=��
=��
=��
=��
=��
=��
=�և+�և+�և+�և+�և+�և+�և+�և+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @8z�@���@£�A�A�RAAQ�Ab�A���A���A���A�u�A���AШ�A��A��B T{BT{B�B�B T{B(T{B0T{B8T{B@T{BHT{BPT{BW�B_�Bg�BpT{BxT{B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=B�*=C CCCCC
.�CCCCCCCCCCC C"C$C&C(C*C,C.C0C2C4C6C8C:C<C>C@CBCDCFCHCJCLCNCPCRCTCVCXCZC\C^C`CbCdCfChCjClCm��CpCrCtCvCxCzC|C~C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�C�
�D HD �HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HD	HD	�HD
HD
�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD~�DHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HDHD�HD HD �HD!HD!�HD"HD"�HD#HD#�HD$HD$�HD%HD%�HD&HD&�HD'HD'�HD(HD(�HD)HD)�HD*HD*�HD+HD+�HD,HD,�HD-HD-�HD.HD.�HD/HD/�HD0HD0�HD1HD1�HD2HD2�HD3HD3�HD4HD4�HD5HD5�HD6HD6�HD7HD7�HD8HD8�HD9HD9�HD:HD:�HD;HD;�HD<HD<�HD=HD=�HD>HD>�HD?HD?�HD@HD@�HDAHDA�HDBHDB�HDCHDC�HDDHDD�HDEHDE�HDFHDF�HDGHDG�HDHHDH�HDIHDI�HDJHDJ�HDKHDK�HDLHDL�HDMHDM�HDNHDN�HDOHDO�HDPHDP�HDQHDQ�HDRHDR�HDSHDS�HDTHDT�HDUHDU�HDVHDV�HDWHDW�HDXHDX�HDYHDY�HDZHDZ�HD[HD[�HD\HD\�HD]HD]�HD^HD^�HD_HD_�HD`HD`�HDaHDa�HDbHDb�HDcHDc�HDdHDd�HDeHDe�HDfHDf�HDgHDg�HDhHDh�HDiHDi�HDjHDj�HDkHDk�HDlHDl�HDmHDm�HDnHDn�HDoHDo�HDpHDp�HDqHDq�HDrHDr��Dr�B1B(B^B	B�B<B�BB.B.BBBgB`BIB>B(BQB�BBQBIB�B'BBqB�B	9B	8B	 B7B{B	QB
=B
B	;B	QB
B	-B	9B	OB	#B	RB	jB
?B	�B�B�BwB�B�B|B|B�B�B�B�B�B�B�B�B�B�B�BkB"�B5�B>BDIBJ7BU�Bh�Bw�Bz�B|sBOBB}�B~�BB4B�%B�HB�eB��B�B�lB�5B��B��B�
B��B�B�B�EB�5B�RB��B��B��B��B� B�FB��B�RB�:B�vB��B�?B��B��B��B��B��B��B�XB�oB�eB�@B�qB�sB�_B�;B�"B��B�B�=B�?B�HB��B��B�B��B��B�hB�{B��B��B��B��B�pB�uB�B�`B��B��B��B��B��B��B��B�^B�RB�lB�nB�dB�{B�dB�XB�cB�WB�NB��B�}B�|B�rB�\B��B��B��B�kB�]B�\B�QB�]B�iB�`B��B��B��B�UB�PB�EB�\B�\B�QB�PB�DB�SB�TB�VB�_B�UB�bB��B��B��B�bB�bB�XB��B��B��B��B��B�pB�dB�|B��B��B��B�{B�KB�KB�oB�oB�cB�bB�WB�cB�oB�pB�}B��B��B��B��B�0B�uB�cB�WB�dB�|B�[B�hB�\B�]B�hB�aB��B��B��B�yB�`B�cB�nB�cB�bB�eB��B�sB�eB�]B��B��B�jB��B��B�|B��B��B��B��B��B��B��B��B��B��B�zB�}B�\B��B��B�_B��B��B��B��B��B��B��B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�3B�B��B��B�B�B��B��B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B�B��B��B�XB��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?���?���?�?�?}?���?��y?�ȴ?���?�`B?�`B?�`B?���?�?�`B?��/?���?�z�?���?�9X?��?�?���?�`B?��?�9X?���?���?��?��?��?�33?�33?�n�?�-?�-?��\?��\?�M�?��!?���?���?��\?���?��\?�J?�J?�Ĝ?�\)?���?��?�v�?�V?�V?�V?�V?���?��-?��-?���?�O�?�V?��?�1?��m?�ƨ?���?�`B?���?���?� �?�I�?�Z?wK�?q�?lI�?c�
?bJ?a��?a��?]p�?R-?Kƨ?E��?>��?2-?$Z?|�?p�?�?�-?��?|�? �? Ĝ?!G�? Ĝ? Ĝ? A�?�R?��?��?�?��?
��?�9?�?�?$�?�
?��?G�>�p�>�{>�r�>�G�>�(�>ڟ�>�b>�>Ձ>��>�bN>��>��>��;>��>���>Ձ>�>�>�5?>�M�>�`B>�h>�33>�F>�F>�33>�!>�->�&�>��>�->�1>�M�>��>��>�9X>�V>�>��y>��
>���>�M�>�Ĝ>�;d>�5?>�(�>�"�>���>���>��>��>���>��>��`>�\)>��>���>�7L>�$�>��>���>��>��>���>��\>�J>}�>r�!>gl�>gl�>gl�>hr�>gl�>fff>fff>fff>gl�>gl�>gl�>gl�>fff>fff>e`B>_;d>["�>V>T��>S��>S��>O�;>H�9>?|�>;dZ>5?}>333>2->/�>)��>"��>�R>��>�>�->��>��>�u>�P>�P>�+>z�>n�>\)>
=q>+>J=���=��=��=���=���=���=ě�=ě�=\=\=\=��=��=�1=��
=���=��=��=��=�t�=�t�=�t�=�t�=�\)=�O�=�O�=�\)=�C�=�+=�+=�o=}�=y�#=q��=e`B=L��=D��=@�=8Q�=0 �=�w=t�=o=o=o=\)=\)=+=t�=t�=\)<�<�<�<�<�/<�o<T��<T��<T��<D��<49X<49X<t�<o;ě�;ě�;D��    ��o��o�o�D���D����o�ě��#�
�D���e`B��C����
��9X��1��9X��9X��j���ͼ�/��`B��`B��h���C���P���0 ŽH�9�P�`�T���]/�q���}󶽇+��t����P���w���T���T���T���罰 Ž�E���vɽ������ě����ͽ�
=��S���xս�F�%�J�o���	7L�
=q�I��I��V�z������P��u����������R��w��w��w��w� Ĝ� Ĝ�!���"��$�/�&�y�(�þ+�0 ž2-�333�8Q�8Q�8Q�:^5�>vɾ@��?|�B�\�C���C���E�˾G��M��Q녾T���Xb�\(��]/�]/�_;d�cS��dZ�fff�hr��hr��ixվixվj~��q���x���z�H��  ���7���\��o�����������˾���7L��=q��O߾�V���;������b�������㾟;d���徤�/��`B��ff���y��r���xվ�1�����������!��9X��?}������^5��j��p���󶾾�۾����J����Ƨ��7L���;��녾����և+�և+��
=��
=��
=��
=��
=��
=��
=�և+��
=��
=��
=��
=��
=��
=��
=��
=��
=��
=��
=�և+�և+�և+�և+�և+�և+�և+�և+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�o<#�<#�<#ٙ<#�<#�}<$>�<#�~<#׮<#װ<#י<#�.<#�w<#��<#��<#��<#�<#ܧ<#��<#��<#ވ<#�><$a<#�<#�T<#��<#��<#ה<#�}<#�v<#��<#�Z<#�f<#�r<#��<#��<#�C<#��<#�<#ד<#��<#�I<#ۆ<#��<#ב<$�<$-<#ݿ<#�;<#�<#��<#׊<#׋<#��<#��<#��<#�p<#�<#��<#�<#��<#��<#�k<#�<$L�<&d?<$'�<$1A<$F2<%�G<+"g<+�<$�<$��<%��<#��<#��<#ذ<$��<'e`<%'�<$�j<%D�<(�T<)�<$��<#��<#��<#��<#�	<#ץ<#�!<#�
<#�4<#�D<#��<#��<#�o<$W/<#��<$5�<$��<$�<$#<#�<#�f<#�<$�<#�7<$�<$�<%w�<$%�<$><$�<#�-<#��<#�<#�Y<#��<#��<#ߣ<#�g<#ס<#�<#��<#�<#�<#�<$2<#�w<#��<$?�<$<#�*<#׎<#��<#�<#�<#�<#�<#�*<$'�<$��<%��<%Tl<%�x<$I`<#�<#�}<#��<#۩<#�?<#�t<#�f<#۸<#�<#ۓ<#�%<#�]<#�<#��<#�<#�<#�<#ߏ<#��<#�j<#��<#�<#��<#�<#��<#כ<#�	<#�I<#�u<#�<$K<$^<#��<#ׅ<#�<#��<#��<#ד<#ׄ<#�<#��<#�^<#�<#�N<#�m<#��<#�<#�d<#�<#��<#��<#׹<#�<#��<$p<#�^<#��<#ۧ<#�><#��<#��<#�<#�<#�u<#�<#�<#�A<#�\<#�<#��<#י<#�<#�^<#ۓ<#�$<#��<#އ<#��<#�<$xs<#��<#�<#ל<#�(<#��<#�u<#��<#א<#ו<#��<#� <$<#�h<#��<#�D<#�j<#ו<#��<#ה<#׌<#��<#�Y<#٬<#�X<#�<#�U<#�B<#ׯ<#��<#��<#�<#�
<#�!<#�N<#ۄ<#�<#�_<#۪<#��<#��<#�6<#�t<#ק<#��<#ב<#��<#�q<#׃<#�J<#� <#�<#�\<#ת<#�j<$|<#ޜ<#מ<#ה<#��<#�<#��<#�<#�<#�)<#׻<#ސ<#ޒ<#ػ<#ל<#��<#��<#ח<#�<#ۋ<#�<#�r<#�<#޵<#ޱ<#�.<#�<#��<#ל<#�<#�o<#�U<#�{<#�i<#��<#�7<#�<#޷<#��<#�<#��<#۾<#�$<#۵<#��<#�m<#�<#�<#ۥ<#�<#�x<#�F<#׶<#�v<#޹<#��<#�Z<#��<#��<#�r<#�<#�4<#�t<#߁<#�6<#��<#�<#�<#ۉ<#�j<#�&<#�2<#׮<#ۣ<#�<#��<#ע<#�"<#�~<#��<#�<#�V<#�S<#��<#׈<#׵<#�\<#��<#ו<#��<#�<#�<#��<#�p<#�B<#�<#�<#�<#�B<#�P<#��<#�|<#��<#�D<#�<#��<#ؤ<#׸<#�B<#��<#��<#�<#��<#��<#�<#�<#מ<#�p<#�\<#�5<#�j<#�8<#ף<#��<#ז<#�f<#��<#��<#�<#�<#��<#�j<#�<#��<#ױ<#�G<#��<#�<#��<#�V<#ڤ<#��<$I<#��<#�\<#��<#�t<#�O<#��<#��<#ڪ<#ً<#��<#�}<#��<#�<#��<#�U<#ޕ<#�j<#��<#��<#�<#ۆ<#�<#ۇ<#��<#�=<#�<#��<#�<$1�<#��<#��<#�m<#ף<#��<#׎<#ׇ<#ׇ<#׆<#׆<#�~<#�<#��<#׏<#׆<#׆<#׆<#׆<#׆<#׆<#׆<#׆<#�}<#�<#ׂ<#׆<#׆<#׆<#׆<#׆<#׆<#׆;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.0825                                                                                                                                                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151816201909291518162019092915182020190929151826202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@333@333@333G�O�G�O�G�O�G�O�Dr�fDr�fDr�fG�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          