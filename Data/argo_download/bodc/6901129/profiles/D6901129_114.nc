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
�  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  J   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  U    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  b�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  eh   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
�  h$   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  s   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  ~   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �l   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
�  �(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  �   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  �   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �X   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �X   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �X   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �X   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    Ĭ   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �@   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �\   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  � Argo profile    3.1 1.2 19500101000000  20210225040027  20210225040027  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               rA   BO  125402                          2C  D   APEX                            6229                            120210                          846 @�d�{��1   @�d�{��@Q��C���2x���F1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @���@���@���AffA>ffA`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D��D�#3B�B5?B`BB'�B[#BhsB�\B��B�B�B�XB��B��B�#B	B	>wB	L�B	W
B	t�B	�uB	ŢB	�B
B
�B
8RB
VB
y�B
��B
��B
�XB
��B
�B
�BB
�mB
�B
�BBuB�B$�B-B2-B<jBF�BT�BcTBl�Bt�B}�B�B�hB��B��B��B��B��B��B��B�B�B�B�B�B�^BĜBĜB�^B�dBŢB��B��B��B��B��B��BȴBƨBÖBB��B�jB�RB�9B�B�B�B�B�B��B��B��B��B��B��B��B��B��B�B�B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�oB�oB�bB�\B�VB�VB�VB�VB�\B�bB�hB�oB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�hB�bB�bB�bB�bB�hB�oB�oB�uB�uB�uB�uB�{B�{B��B��B�{B�{B�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B��B��B��B��B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�{B�uB�uB�uB�uB�oB�oB�hB�hB�bB�bB�bB�\B�bB�bB�oB�uB�{B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@q��@lZ@b-@X��@W|�@E�@H�u@H1'@DI�@>�R@<�@@Q�@6��@'�P@	G�?�o?� �?"�=�h>1&�>�`B?�!?1'>��y>��H??!��?4z�?BM�?O�?Y��?f��?l��?t��?xQ�?|(�?��?��`?�ȴ?�=q?�v�?���?���?�\)?��9?���?š�?�ȴ?Ł?�33?�Ĝ?��T?�S�?�?�
=?�O�?���?���?� �?�?Ձ?�  ?�J?ݲ-?�O�?��?ߝ�?��?�X?��
?�?�!?�|�?�j?���?�^?��T?޸R?�p�?�^5?�S�?�/?Ǯ?��?���?��?���?�b?�Ĝ?�ȴ?�33?�Z?���?�?�x�?�~�?��?��7?���?���?� �?��?���?���?�?��9?��/?��F?�J?���?�V?��D?�"�?�~�?��/?�?�A�?u�?e��?Q��?KC�?I�^?A��??;d?=�-?=�?<j?:�H?:�H?;�m?=�-?=�-?>v�?@Ĝ?B�\?A%?>�R?>v�?@�?H1'?I�^?J=q?J~�?Hr�?F�y?F$�?BM�??;d?;dZ?7�P?5?4�j?4z�?1hs?(��?#��?#�
? Ĝ?X?hs?�?��?�?7
=?8Q�?7
=?6ȴ?6ȴ?8Q�?3��?9�?9��?6?.�?'l�?"M�?;d?"J?"��? Ĝ?�-?��?�+?33?��?1?r�?��?��>�v�>�X>�{>��
>ܬ>�(�>޸R>߾w>ݲ->�(�>�b>��;>�o>��m>�ȴ>�?}>���>���>���>��>��>w��>aG�>Z�>L��>E��>Kƨ>M��>L��>M��>M��>J��>N�>R�>S��>O�;>V>I�^>@�><j>49X>/�>&�y> Ĝ>�->�+>z�>�P>�u>�+>n�>\)>
=q>+=��#=��=�=�/=��`=\=� �=��T=���=�-=��
=��P=��-=��=m�h=<j=t�<�<�t�<D��<49X��o�49X�t�;ě�<���<�/<�h<��=o=�w=\)<�/<���<�j<�o;��
��o�ě��49X�D���49X;�o<D��<49X;�`B;D��    �D���o���
��o�o;o<o<T��<��
<���<���=t�=�w=0 �=,1=#�
=�P=t�=#�
=�w=#�
=49X=H�9=H�9=D��=<j='�='�=#�
=#�
=�w=�P=�P=t�=t�=��=�w=#�
='�='�='�='�='�=#�
=�P=C�<��<�`B<�9X<���<���<���<�/<�/<���<�1<���<��
<�o<D��<o;ě�;��
    �t��D���49X��o���
��9X���C��o�C����t��\)�,1�,1�'C��+�t��'8Q�@��L�ͽixսy�#��%��%��+��O߽�hs����������
���-��E���E���j����\�Ƨ�������
=��/��G���`B���F���#�%�$ݾ
=q�V�hs�n��z��u������ Ĝ�"��%�T�,1�.{�/��1&�8Q�;dZ�>vɾ@��A�7�D���F��J���M��P�`�R�V�Z��]/�bMӾfff�j~��l�D�m�h�q���r�!�vȴ�vȴ�w�پ{�m�~�۾�%��J����������˾�+���9��=q������C���ƨ�������hs��񪾓�Ͼ�
=��b���u������(����-���w���徥�T��l�����þ�xվ��羫���V���h�����&龲-��33��?}��KǾ�Q쾹X���#���H��dZ���m��󶾾vɾ�|���7�Õ��ě���$ݾǮ��7L��=q��O߾�V��\)���`����z�����Ձ�����և+�և+�׍P�ؓu�ٙ���"Ѿۥ��(��ܬ��/�ݲ-�޸R��Ĝ��G���MӾ�����
��`B��ff��l�����þ�xվ�������~�����1��D��V��{��{���׾�׾�&��-��!��F��9X����E���ȴ��ȴ��KǾ��پ�Q��X���#���H��dZ��푾�p���p�����vɿ   � ��%�G�����Mӿ����o�S�������Z�����/�������`B��T�ff�ff�ff������y�+�l��1'��9��ÿ	7L�	xտ	�^�	��
~��
���C��ƨ�1��D�V��h�{�V��������� ſbN��`�&�hs����녿n���!��33�9X��j����?}��E���+�
=��P�Q���X�X���X�X��#����"ѿdZ��m���p���-��-�p���-��-��-�;d�|�;d�vɿ5?�|� A�� Ĝ�!G��!�7�!�7�!�7�!���"J�!���!���!�7�!�7�!���"J�"Mӿ"�\�"Mӿ"�\�#S��$��$���$�/�%��%�˿%�˿%�˿%�˿%�T�%�T�%�˿%�˿%�T�%�˿%�T�%�˿%��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�
>@�=q@�=qA�A=�A^�RA~�RA�\)A�\)A�\)A�\)A�\)A�\)A�\)A�\)B�B�B�B�B'�B/�B7�B?�BG�BO�BW�B_�Bg�Bo�Bw�B�B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
B��
C�C�C�C�C	�C�C�C�C�C�C�C�C�C�C�C�C!�C#�C%�C'�C)�C+�C-�C/�C1�C3�C5�C7�C9�C;�C=�C?�CA�CC�CE�CG�CI�CK�CM�CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc�Ce�Cg�Ci�Ck�Cm�Co�Cq�Cs�Cu�Cw�Cy�C{�C}�C�C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���D z�D ��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��D	z�D	��D
z�D
��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��Dz�D��D z�D ��D!z�D!��D"z�D"��D#z�D#��D$z�D$��D%z�D%��D&z�D&��D'z�D'��D(z�D(��D)z�D)��D*z�D*��D+z�D+��D,z�D,��D-z�D-��D.z�D.��D/z�D/��D0z�D0��D1z�D1��D2z�D2��D3z�D3��D4z�D4��D5z�D5��D6z�D6��D7z�D7��D8z�D8��D9z�D9��D:z�D:��D;z�D;��D<z�D<��D=z�D=��D>z�D>��D?z�D?��D@z�D@��DAz�DA��DBz�DB��DCz�DC��DDz�DD��DEz�DE��DFz�DF��DGz�DG��DHz�DH��DIz�DI��DJz�DJ��DKz�DK��DLz�DL��DMz�DM��DNz�DN��DOz�DO��DPz�DP��DQz�DQ��DRz�DR��DSz�DS��DTz�DT��DUz�DU��DVz�DV��DWz�DW��DXz�DX��DYz�DY��DZz�DZ��D[z�D[��D\z�D\��D]z�D]��D^z�D^��D_z�D_��D`z�D`��Daz�Da��Dbz�Db��Dcz�Dc��Ddz�Dd��Dez�De��Dfz�Df��Dgz�Dg��Dhz�Dh��Diz�Di��Djz�Dj��Dkz�Dk��Dlz�Dl��Dmz�Dm��Dnz�Dn��Doz�Do��Dpz�Dp��Dqz�Dq��Drz�Dr��Dsz�Ds��Dtz�Dt��Duz�Du��Dvz�Dv��Dwz�Dw��Dxz�Dx��Dyz�Dy��Dzz�Dz��D{z�D{��D|z�D|��D}z�D}��D~z�D~��Dz�D��D�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD��qD�=qD�}qD��qD�D� �B!�B<jBf�B(<BfpBf�B��B��B��B�B��BԘB�sB�bB	�B	Q�B	ezB	nB	ryB	��B	��B	�5B
�B
�B
3�B
S�B
v]B
��B
��B
�WB
�TB
��B
��B
��B
�B
�.B
��B:B\B#QB,B/�B9�BC*BP�Bb�Bl#Bu7B~�BeB�rB�mB��B�AB�B�-B�XB��B��B�B��B�5B��B��B�+B� B�xB�B��B�B�B�B��BʮB��B�-B�7B�*B��B�B��B�rB�}B�UB�wB�B�B��B��B�FB��B�EB� B��B�{B�B�)B��B�B��B�7B�mB��B��B��B�gB�nB��B��B�dB��B�hB�1B��B��B��B��B��B�5B��B��B��B��B��B�OB��B��B�[B�0B�B�kB�LB�B�)B��B��B��B�B�3B�MB��B��B�B��B��B�jB�IB�fB�aB��B��B��B�;B�)B��B��B�#B��B��B�cB�|B�pB��B��B�#B��B��B��B��B�B��B��B�RB�PB��B�oB�`B��B�HB�|B��B�B�vB�VB��B�yB�B�WB��B�>B��B��B�XB��B�sB��B��B��B�B�{B��B�RB�B��B�<B��B�VB��B�0B��B�yB��B�B��B�B�OB�wB�eB�uB��B�FB�FB�qB��B�?B�B��B��B��B��B��B��B��B��B��B�VB�pB��B��B��B��B��B��B��B��B��B��B��B��B��B�]B�LB��B��B�XB�B��B�B��B��B��B��B�vB�B��B�DB��B��B�B�_B�hB�wB�6B��B��B��B��B��B�B��B��B��B��B�lB��B�'B��B��B��B��B��B��B�iB�~B�sB�aB�MB�]B�IB�iB��B�8B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�2B�B��B��B��B��B�B�'B�B��B�(B�(B�(B�B�B�6B�`B�B��B�5B�(B�B�_B�B��B�	B�B��B��B�IB��B��B��B��B�!B�:B�.B�B�#B�QB�/B�B��B�!B�!B�B�B�DB�B�B�9B�B��B� B�B�	B�B�-B�B�B�!B�B�B�9B�B�"B�/B�:B�.B�*B�B�
B�B�-B�
B�B�9B�B�B�?B�B�B�B�KB�B�B�B�B�B�B�(B�B�B�B�B�(B�B�4B�)B�(B�B�B�&B�B�&B��B� B�!B�B�B�B� B��B�
B�B�B�B��B��B��B� B�B�.B�B�B�8B�
B��B�B�B�B�B�4B�3B�B��B�B��B��B�B�B��B�B�B�B�B�B�B�B�B��B��B��B��B�B��B��B�B�B�B�B�B�B��B�,B��B��B�	B�B�	B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B�B��B�B��B��B��B��B��B�B��B�B�B��B�B�B�	B�	B��B��B��B��B�%B��B�B��B��B��B��B�B��B��B��B�B�B��B��B�B��B��B�	B�
B�B�B��B��B��B��B��B�	B�B�B��B�B�B�B��B��B��B��B��B��B�0B��B��B��B��B�+B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�	B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��@q��@lZ@b-@X��@W|�@E�@H�u@H1'@DI�@>�R@<�@@Q�@6��@'�P@	G�?�o?� �?"�=�h>1&�>�`B?�!?1'>��y>��H??!��?4z�?BM�?O�?Y��?f��?l��?t��?xQ�?|(�?��?��`?�ȴ?�=q?�v�?���?���?�\)?��9?���?š�?�ȴ?Ł?�33?�Ĝ?��T?�S�?�?�
=?�O�?���?���?� �?�?Ձ?�  ?�J?ݲ-?�O�?��?ߝ�?��?�X?��
?�?�!?�|�?�j?���?�^?��T?޸R?�p�?�^5?�S�?�/?Ǯ?��?���?��?���?�b?�Ĝ?�ȴ?�33?�Z?���?�?�x�?�~�?��?��7?���?���?� �?��?���?���?�?��9?��/?��F?�J?���?�V?��D?�"�?�~�?��/?�?�A�?u�?e��?Q��?KC�?I�^?A��??;d?=�-?=�?<j?:�H?:�H?;�m?=�-?=�-?>v�?@Ĝ?B�\?A%?>�R?>v�?@�?H1'?I�^?J=q?J~�?Hr�?F�y?F$�?BM�??;d?;dZ?7�P?5?4�j?4z�?1hs?(��?#��?#�
? Ĝ?X?hs?�?��?�?7
=?8Q�?7
=?6ȴ?6ȴ?8Q�?3��?9�?9��?6?.�?'l�?"M�?;d?"J?"��? Ĝ?�-?��?�+?33?��?1?r�?��?��>�v�>�X>�{>��
>ܬ>�(�>޸R>߾w>ݲ->�(�>�b>��;>�o>��m>�ȴ>�?}>���>���>���>��>��>w��>aG�>Z�>L��>E��>Kƨ>M��>L��>M��>M��>J��>N�>R�>S��>O�;>V>I�^>@�><j>49X>/�>&�y> Ĝ>�->�+>z�>�P>�u>�+>n�>\)>
=q>+=��#=��=�=�/=��`=\=� �=��T=���=�-=��
=��P=��-=��=m�h=<j=t�<�<�t�<D��<49X��o�49X�t�;ě�<���<�/<�h<��=o=�w=\)<�/<���<�j<�o;��
��o�ě��49X�D���49X;�o<D��<49X;�`B;D��    �D���o���
��o�o;o<o<T��<��
<���<���=t�=�w=0 �=,1=#�
=�P=t�=#�
=�w=#�
=49X=H�9=H�9=D��=<j='�='�=#�
=#�
=�w=�P=�P=t�=t�=��=�w=#�
='�='�='�='�='�=#�
=�P=C�<��<�`B<�9X<���<���<���<�/<�/<���<�1<���<��
<�o<D��<o;ě�;��
    �t��D���49X��o���
��9X���C��o�C����t��\)�,1�,1�'C��+�t��'8Q�@��L�ͽixսy�#��%��%��+��O߽�hs����������
���-��E���E���j����\�Ƨ�������
=��/��G���`B���F���#�%�$ݾ
=q�V�hs�n��z��u������ Ĝ�"��%�T�,1�.{�/��1&�8Q�;dZ�>vɾ@��A�7�D���F��J���M��P�`�R�V�Z��]/�bMӾfff�j~��l�D�m�h�q���r�!�vȴ�vȴ�w�پ{�m�~�۾�%��J����������˾�+���9��=q������C���ƨ�������hs��񪾓�Ͼ�
=��b���u������(����-���w���徥�T��l�����þ�xվ��羫���V���h�����&龲-��33��?}��KǾ�Q쾹X���#���H��dZ���m��󶾾vɾ�|���7�Õ��ě���$ݾǮ��7L��=q��O߾�V��\)���`����z�����Ձ�����և+�և+�׍P�ؓu�ٙ���"Ѿۥ��(��ܬ��/�ݲ-�޸R��Ĝ��G���MӾ�����
��`B��ff��l�����þ�xվ�������~�����1��D��V��{��{���׾�׾�&��-��!��F��9X����E���ȴ��ȴ��KǾ��پ�Q��X���#���H��dZ��푾�p���p�����vɿ   � ��%�G�����Mӿ����o�S�������Z�����/�������`B��T�ff�ff�ff������y�+�l��1'��9��ÿ	7L�	xտ	�^�	��
~��
���C��ƨ�1��D�V��h�{�V��������� ſbN��`�&�hs����녿n���!��33�9X��j����?}��E���+�
=��P�Q���X�X���X�X��#����"ѿdZ��m���p���-��-�p���-��-��-�;d�|�;d�vɿ5?�|� A�� Ĝ�!G��!�7�!�7�!�7�!���"J�!���!���!�7�!�7�!���"J�"Mӿ"�\�"Mӿ"�\�#S��$��$���$�/�%��%�˿%�˿%�˿%�˿%�T�%�T�%�˿%�˿%�T�%�˿%�T�%�˿%��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<@C<G�:<B�<#�<r}�<&�E<#�<)M�<.�'<&n�<)$�<Cn�<u�<�%3<�_<��<�{<��'<(<��<=��<%��<.-�<&�r<3�6<(9�<-m:<*;<(s�<' V<(�d<%Nw<%�r<$5;<$o5<-9�<,��<'޽<%W�<%ӫ<$�<)W+<*b<-_n<0�<$`
<$ �<#�R<$"<4�6<&��<$f�<$q�<#�z<,��<%�~<(1/<=��<%p�<&�<&��<$g<1��<;�5<%��<@� <$��<29[</W;<$S�<$��<$�<$��<#�
<$ú<%k�<(�,<$�<$��<(�G<'Ϛ<'+p<+��<'�!<#��<#�<#�w<)��<-<%H<#�:<#�e<*�)<$!�<#�Z<$q�<&m,<#��<#�/<$z<'\*<#��<#��<#�B<$X�<%L�<#�N<$�<$f�<$ D<$$�<$�<#��<&��<#�5<&��<'"�<*��<-�l<$�<#�`<%Z�<$ �<#��<#�:<#��<#��<#ה<#�$<#��<#��<#��<$�<#�U<#�T<#�<#�<$J<%k�<#��<#�r<#�`<#�q<#��<#�<$0�<$*<$2�<$3�<#�<#�<#�<$|<%�d<$�#<#ס<$�<%7 <%i�<#�1<#��<'t[<4�e<#�<#�P<#�<#�p<#�<$=�<$��<#��<$"w<%,�<%Dj<$�K<$�<$�<#�4<#�<$K<$��<#�G<$q<$�<$V�<$&<#�`<$<<$G�<$�<$�<$��<$$�<#�<#��<#ڿ<#�v<#��<#�<$C$<$�<$*�<#�(<#��<$!,<$ݎ<$JQ<$��<$-�<$�<$�<#�j<$G<#��<#�<#�u<#�<#ت<#׆<#�`<#��<#��<#ئ<#ژ<#��<$�<#��<#��<#�<#ߟ<#��<#�<#�j<#�<#׻<#ފ<#��<#׻<#ڀ<#�-<#ޑ<#�@<#�P<#��<#�N<#�u<#�<#�<#�<#��<#�<#��<#�(<#�@<#�l<$<#�7<$R<#�*<#�<$
�<#�\<#�2<$.T<#�<#ܐ<$H<$3H<#�X<#��<#�H<#ٗ<#�<#��<#�<#ײ<#׼<#�<$�<#��<#��<#�I<#�<#�X<$I�<#�\<#�<#ڄ<#��<#�	<#�!<#ݎ<#ܫ<#�u<#ے<#�.<#�<#��<#�<#��<#�c<$	G<#�<#�2<#�<#��<#�p<#�<#�x<#� <#�<#��<#�N<#��<#�
<#ם<#�<#�W<#�<#�q<#�
<#�z<#�M<#�<#׋<#��<#�%<#��<#�u<#ײ<#ׁ<#�]<#״<#�<#��<#�<#�<#�8<#�P<#�r<#�<#��<#�E<#�r<#ז<#�<#׃<#�{<#�t<#�r<#�?<#ה<#�<#�-<#�<#٧<#�#<#޷<#�Y<#ס<#�f<#��<#ڻ<#�}<#ڭ<#�_<#�p<#�<#�U<#�G<#��<#��<#��<#�v<#�U<#ח<#�B<#�.<#ۆ<#ׂ<#�a<#��<#��<#ה<#ײ<#��<#�<#נ<#�O<#׆<#�[<#��<#׊<#�<#י<#�<#מ<#ז<#��<#ז<#ש<#�7<#ח<#�<#�{<#ޣ<#�@<#��<#ؔ<#�<#י<#�.<#�<#נ<#�&<#׋<#�h<#�*<#ש<#�<#׾<#�&<#�P<#��<#׉<#�<#��<#ת<#�@<#�	<#��<#ס<#�<#�\<#�6<#ދ<#ۗ<#�I<#׆<#�<#��<#�<#��<#�~<#�<#�<#�<#��<#ק<#�<#�<#ו<#�<#�
<#��<#�<#�
<#�<#�<#׽<#�<#�<#׾<#� <#נ<#�<#۬<#�<#� <#ۯ<#�2<#��<#�<#�<#׆<#�<#�<#��<#��<#�<#��<#�6<#ם<#ף<#��<#ܾ<#�m<#ׇ<#�
<#׷<#�<#�<#�7<#�<#ך<#��<#�<#׉<#�m<#ي<#�k<#׶<#�<<#׷<#ך<#�<#�a<#��<#�<#�
<#�
<#�p<#�
<#�f<#׊<#ב<#ך<#��<#�<#�
<#�
<#�
<#�<#ן<#�&<#�<#׀<#�<#ך<#��<#ך<#׈<#�<#ׅ<#�<#�
<#�n<#�<#ׂ<#�<#�
<#�<#�}<#�X<#��<#׆<#�t<#�<#׀<#�<#׆<#�<#��<#�<#�
<#�t<#�
<#�<#�<#׆<#�<#ׂ<#�<#��<#�
<#�n<#�
<#�<#��<#ז<#׊<#�<#׌<#ג<#ׁ<#�d<#�
<#�<#�<#ׂ<#�<#�<#�<#�
<#�s<#�r<#�<#׏<#ׁ<#�k<#�s<#�
<#�m<#�
<#�<#�<#ׇ<#�<#ׂ<#�<#�<#�<#�<#�<#ׇ<#�<#׋<#��<#�<#ה<#�m<#ד<#׋<#�<#�
<#ם<#�<#�<#�<#ׅ<#�<#�<#�<#�<#׃<#�<#�<#�<#�'<#ם<#�<#�<#��<#�<#�<#׏<#מ<#�_<#��<#�<#�~<#ؐ<#�<#�b<#ב<#�<#�k<#�<#ז<#�<#��<#�<#�y<#ش<#�
<#�r<#�M<#�<#�<#��<#�R<#�<#ݞ<#�<#כ<#׌<#�<#�u<#�q<#�
<#�
<#ع<#ׁ<#��<#�I<#�<#�<#�<#�
<#؜<#�<#�<#�<#ז<#�<#�<#׀<#�j<#�x<#�p<#�
<#�z<#ؼ<#�w<#�
<#ذ<#�<#�D<#�}<#�x;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=0.08                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929144034201909291440342019092914403820190929144044202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@���@���@���G�O�G�O�G�O�G�O�D�#3D�#3D�#3G�O�G�O�G�� G�� G�� G�� G�� G�� G��                                 6389758         6389758         131072                                          