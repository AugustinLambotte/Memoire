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
�  I�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  T�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  bD   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  d�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
�  g�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  r|   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  }L   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     
�  �8   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  �   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     
�  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �\   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �x   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    Ô   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ð   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  Ō   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �|   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        Ƙ   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ƴ   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225041100  20210225041100  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               {A   BO  125411                          2C  D   APEX                            6229                            120210                          846 @�z�+<@1   @�z�+<@@Q8ě��T�6"��`B1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @@  @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@ffBH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*�C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DAfDA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�fB��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B	B	+B	1B	
=B	JB	VB	bB	uB	�B	�B	.B	6FB	G�B	[#B	k�B	��B	�-B	�^B	��B	�B	��B
.B
O�B
l�B
�DB
��B
�LB
�B�B(�B-B.B9XBD�BB�B7LB;dBE�BJ�BVB\)B[#BVBN�BQ�BW
BaHBcTBdZBe`BffBe`BffBgmBgmBgmBhsBiyBm�Bp�Bq�Bo�Bm�Bm�Bm�Bn�Bo�Bo�Bp�Bq�Bs�Bv�Bu�Bu�Bu�Bt�Bv�Bx�By�Bz�Bz�By�Bx�Bx�Bz�B|�B}�B|�B~�B�B�B�%B�DB�\B�\B�\B�hB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B����m�������#��xս�`B��l���`B��`B��S���`B��S��������j��j�\��������vɽȴ9��/����xս�/��j�q����/=�w=���>�->p��>�V?��?$�>�j?n�?��?>��?_|�?t9X?�K�?�$�?��?�S�?�`B?��?�?ۅ?�1?߾w?�|�?�j?և+?��
?��?~��?o�;?k�?c��?Tz�?+>��>�?�9?
��?��?	7L?1'?`B?�7? Ĝ>�j>�K�>�&�>�D>��?1'?
=q>���>�7L>�K�>�{>�>�x�>���>�1>�V>�9X>�>���>�"�>�z�>�bN>�n�>���>�
=>�b>�b>��9>o��>o��>k�>��>���>�o>z�H>�z�>��>�z�>��>Ǯ>Ǯ>Õ�>�p�>և+?�?33?��?dZ?dZ?dZ?dZ?dZ?dZ?��?�H?dZ?"�?"�?��?�H??}? �?�?Q�?�P?ȴ??}?z�?t�?�`?\)?	��?��?��?�?�/?J? �>��>�^5>�^5>���?Z?��? A�>��>�j>�p�>�|�? Ĝ? �>��>�Q�>�E�>� �>�D>�~�>���>��>��>�l�>��/>��
>�S�>���>ۥ�>�"�>��>��>׍P>�
=>���>���>�t�>��>��;>�\)>��>�V>���>�I�>�7L>ě�>��7>��>��m>��H>�X>���>�K�>�>�->���>�V>��>��T>��T>��/>�G�>�/>�"�>���>���>��u>��P>��P>�>�z�>�z�>���>��>��>��`>��>�O�>�C�>��^>��>t�j>t�j>t�j>r�!>q��>q��>p��>o��>m�h>l�D>k�>gl�>bM�>\(�>\(�>S��>G�><j>9X>6E�>333>0 �>/�>.{>(��>$�/>!��>�->��>�P>�>n�>O�>
=q>%=��#=�F==�l�=�l�=�`B=�/=�
==��=�"�=�/=�>o>�>+>$�>$�>$�>$�>�>1'>1'>+>+>+>1'>1'>O�>bN>bN>I�>$�>�>+>$�>J>%=��m=��#=��#==�x�=�`B=���=��
=��=�C�=��=�+=�+=y�#=L��=L��=H�9=H�9=D��='�=�w=#�
=#�
=#�
=#�
=#�
=#�
=#�
=#�
=#�
=t�=o=o<��<�<�h<�h<�/<���<ě�<�1<���<T��;�`B;o��o�o�o�o��o�ě��ě���`B�D���o�o�o:�o;��
<t�;ě�;o;��
;��
;D��;o��o�ě��#�
�T���e`B���㼴9X��j�ě�������h�+�\)���'49X�@��L�ͽe`B�q����o��+���P��1��{��9X��Q����\�Ƨ�������`������/��l������ٽ��m�   �+�	7L�C��O߾1'�$ݾ
=q�\)�z�z��+�����-��R� Ĝ�$�/�')��/��:^5�;dZ�<j�A�7�F��G��I�^�M��O�;�R�S�ϾT���V�V�\(��_;d�_;d�^5?�_;d�]/�^5?�_;d�`A��_;d�bMӾdZ�fff�n���o���p�׾p�׾p�׾r�!�s�F�t�j�{�m�~�۾�  ��%��%���7���7��o�����$ݾ�$ݾ���+������9��=q��ƨ��I���I���I����;�I���I���I���V��t�����������+���u��"Ѿ����/���-��Ĝ���徣�
��`B��l���1����������9X��?}��ȴ��Q쾹�#���m��󶾾�۾�  ��  ��%��%��%��o��+��=q���;����n�����z���ٙ��ڟ���"Ѿۥ��(��ܬ�ݲ-��;d��A���S���Z���T���y���y��l����y��r����������D��{��{��{���{��{����׾���Q�������#���#��^5��^5��^5��^5���H��dZ���m��p���푾��m��j��푾���vɾ�vɾ��۾�|� Ĝ��7����Mӿ��S��S��������S��Z��/��/���`B��˿$ݿ$ݿff�ff�ff����y�+�1'��ÿ�9��ÿ�ÿ	xտ	�^�
~��
~��
~��
~��
~��
~��
~���1�I��I��I��I��I���Ϳ�ͿO߿V���������bN�bN�bN�bN�&�&�&����-�n��-����33�t��z��j���������������E��E��ȴ�ȴ��+��+��+�E��E�����������������������������111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@@  @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@ffBH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*�C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DAfDA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�fB��B��B��B�pB��B��B��B��B��B��B��B�uB�rB��B��B�B��B�B	B	uB	�B	
�B	B	B	�B	�B	B	�B	,5B	2�B	C�B	X�B	aB	�bB	�gB	��B	�B	�9B	��B
)�B
J�B
g B
�8B
��B
�B
��B�B&�B,�B,�B9rBE�BEBD(BI�BJ�BM�BV�B]�B^!B]�BX;BR�BS�B`�Bc�BdSBe�Bf�BfBf�Bg�Bg�BhBh�Bg�Bk�BpJBs�Bt�BofBnuBm�Bn�Bo�Bo[Bp�Bq Bs�BwqBw�BviBv"Bt�Bv&By	By�Bz�B|LB{kBx�Bx�ByuB{�B�B}aB|�B�nB��B��B��B�CB��B��B��B�JB�~B��B��B��B��B��B��B��B�B��B��B��B��B��B�B��B��B�B�1B�B�B�AB�!B�.B�rB�HB��B�$B��B�B�yB�kB�1B�/B�8B��B�xB�B�B��B�4B�B��B��B��B��B�<B�vB�'B��B�LB�)B�B�B��B�B�3B�B�B�B��B�B�B�B�B�B�'B�B�B�B�'B�B�B�B�B�B�>B�gB�:B�GB�"B�B�B�B��B�B�FB�FB�B�iB�#B��B�
B�BB�KB�B��B�B�B�B��B�B�B��B��B�B�B�B�B�B�B�B�kB��B��B��B��B��B��B��B��B��B��B��B�B�B�#B��B�<B�lB�^B��B��B��B��B��B��B�B�B��B�B�B��B��B��B�	B��B�9B��B��B��B��B��B��B�B��B��B��B��B�IB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�/B�EB�
B��B�
B�-B�B� B�	B��B�9B� B�B��B��B�HB�(B�B��B��B�"B�fB��B��B��B��B�9B�B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B�	B�	B��B�B�B�;B�;B�/B�B�B��B��B�B�B��B� B��B��B��B��B��B��B��B�%B�,B��B��B�B�B�(B�@B�4B�'B�B�?B�(B�B�B�B�(B�3B�B�(B�(B�+B�*B�'B�OB�*B�?B�B�gB�{B�B�%B�B�3B�B�B�B�#B�B�.B�:B�LB�	B�B�B�KB�B�B�B��B��B�-B�8B�6B��B�B�'B�B�B�B�'B�B�B�6B�vB�B��B�.B�-B��B�
B�!B�
B�B��B��B��B��B�9B�B��B��B��B��B��B��B��B��B�B�B�B�VB�B��B��B��B�
B��B��B�FB�B��B�B��B��B��B�B�!B�	B��B��B��B��B�
B�B�B��B��B��B��B��B��B��B�$B�iB�$B��B��B��B�"B�.B�B��B��B�5B�B�B�B�B�WB�)B�4B�'B�B�B�B�B�B�B�B�B��B�B��B��B�B�EB�/B�"B�B�8B��B�	B�B�8B��B��B��B��B��B��B�B��B�,B��B�B��B��B��B��B�B�	B��B��B��B�B��B��B��B��B��B�B��B��B�YB�B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B�B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B�B��B��B�B�B��B��B��B�B��B��B�B�B��B��B�B��B��B��B��B�B��B��B�B��B��B��B��B��B��B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B����m�������#��xս�`B��l���`B��`B��S���`B��S��������j��j�\��������vɽȴ9��/����xս�/��j�q����/=�w=���>�->p��>�V?��?$�>�j?n�?��?>��?_|�?t9X?�K�?�$�?��?�S�?�`B?��?�?ۅ?�1?߾w?�|�?�j?և+?��
?��?~��?o�;?k�?c��?Tz�?+>��>�?�9?
��?��?	7L?1'?`B?�7? Ĝ>�j>�K�>�&�>�D>��?1'?
=q>���>�7L>�K�>�{>�>�x�>���>�1>�V>�9X>�>���>�"�>�z�>�bN>�n�>���>�
=>�b>�b>��9>o��>o��>k�>��>���>�o>z�H>�z�>��>�z�>��>Ǯ>Ǯ>Õ�>�p�>և+?�?33?��?dZ?dZ?dZ?dZ?dZ?dZ?��?�H?dZ?"�?"�?��?�H??}? �?�?Q�?�P?ȴ??}?z�?t�?�`?\)?	��?��?��?�?�/?J? �>��>�^5>�^5>���?Z?��? A�>��>�j>�p�>�|�? Ĝ? �>��>�Q�>�E�>� �>�D>�~�>���>��>��>�l�>��/>��
>�S�>���>ۥ�>�"�>��>��>׍P>�
=>���>���>�t�>��>��;>�\)>��>�V>���>�I�>�7L>ě�>��7>��>��m>��H>�X>���>�K�>�>�->���>�V>��>��T>��T>��/>�G�>�/>�"�>���>���>��u>��P>��P>�>�z�>�z�>���>��>��>��`>��>�O�>�C�>��^>��>t�j>t�j>t�j>r�!>q��>q��>p��>o��>m�h>l�D>k�>gl�>bM�>\(�>\(�>S��>G�><j>9X>6E�>333>0 �>/�>.{>(��>$�/>!��>�->��>�P>�>n�>O�>
=q>%=��#=�F==�l�=�l�=�`B=�/=�
==��=�"�=�/=�>o>�>+>$�>$�>$�>$�>�>1'>1'>+>+>+>1'>1'>O�>bN>bN>I�>$�>�>+>$�>J>%=��m=��#=��#==�x�=�`B=���=��
=��=�C�=��=�+=�+=y�#=L��=L��=H�9=H�9=D��='�=�w=#�
=#�
=#�
=#�
=#�
=#�
=#�
=#�
=#�
=t�=o=o<��<�<�h<�h<�/<���<ě�<�1<���<T��;�`B;o��o�o�o�o��o�ě��ě���`B�D���o�o�o:�o;��
<t�;ě�;o;��
;��
;D��;o��o�ě��#�
�T���e`B���㼴9X��j�ě�������h�+�\)���'49X�@��L�ͽe`B�q����o��+���P��1��{��9X��Q����\�Ƨ�������`������/��l������ٽ��m�   �+�	7L�C��O߾1'�$ݾ
=q�\)�z�z��+�����-��R� Ĝ�$�/�')��/��:^5�;dZ�<j�A�7�F��G��I�^�M��O�;�R�S�ϾT���V�V�\(��_;d�_;d�^5?�_;d�]/�^5?�_;d�`A��_;d�bMӾdZ�fff�n���o���p�׾p�׾p�׾r�!�s�F�t�j�{�m�~�۾�  ��%��%���7���7��o�����$ݾ�$ݾ���+������9��=q��ƨ��I���I���I����;�I���I���I���V��t�����������+���u��"Ѿ����/���-��Ĝ���徣�
��`B��l���1����������9X��?}��ȴ��Q쾹�#���m��󶾾�۾�  ��  ��%��%��%��o��+��=q���;����n�����z���ٙ��ڟ���"Ѿۥ��(��ܬ�ݲ-��;d��A���S���Z���T���y���y��l����y��r����������D��{��{��{���{��{����׾���Q�������#���#��^5��^5��^5��^5���H��dZ���m��p���푾��m��j��푾���vɾ�vɾ��۾�|� Ĝ��7����Mӿ��S��S��������S��Z��/��/���`B��˿$ݿ$ݿff�ff�ff����y�+�1'��ÿ�9��ÿ�ÿ	xտ	�^�
~��
~��
~��
~��
~��
~��
~���1�I��I��I��I��I���Ϳ�ͿO߿V���������bN�bN�bN�bN�&�&�&����-�n��-����33�t��z��j���������������E��E��ȴ�ȴ��+��+��+�E��E�����������������������������111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#ײ<#�<#�<#�<#�-<#�d<#�|<#�<#�l<#�d<#��<#�3<#�`<#��<#�<#��<#�f<#�<<#װ<#�<$8<$�<#�<#�<$Y�<%��<%�|<*�:<&��<.[�</K<(��<k"i<#��<$��<,��<&G<D٫<>*�<0֖<6a$<9z�<' Y<*�O<@9<;�<BwW<'�n<#�<$�<#�<$��<(��<�yB<��<7,P<*A<$lx<%�<*��<O[:<^+�<$5�<,K/<#��<#��<#�1<#�
<$G<$<�<#�<$�<$�<$<#�<%�y<&_M<#��<'tC<5�n<&o&<$t�<#�A<#��<#�7<#�<#�<$0<#��<$,�<&W�<$*�<#�|<#�m<$(�<#�O<#��<#�g<%fO<%��<#�><#�'<%i�<%hI<%�k<#��<'&<#�@<$w<(��<'ڧ<#�<#�F<$�<(�(<8^T<'�<%�<#�<#�<#�
<#�
<#�
<#�<#��<#�l<#ش<#�r<#�<#�{<#��<$�u<$�$<#��<&#w<#�&<#��<#��<#�%<#�~<$�<#�^<$��<#�<#�<#�?<$\<$<#�S<#�x<#��<#�
<#�<$p`<#�:<$!�<#��<#�
<#��<#ޱ<#�<#׈<#�e<$�<#�;<$�<#�q<#��<#�A<#��<#�<#ס<#�F<#��<#ף<#�<$V<#׎<#��<#�<#�<#ר<#�L<#�<#ס<#�R<#�S<#ט<#׆<#ח<#�<#׿<#��<$�<#�<#��<#޾<#�,<#�F<#�,<#נ<#ۜ<#��<#��<#�<${<#� <#�<#�<#��<#�~<#��<#ע<#��<#��<#��<#�<#�?<#�<#�<#׆<#��<#��<#�<#�[<#��<#��<#�Q<$F<$w�<#�$<#�<#��<#׆<#�<#׃<#׎<#��<#׏<#ע<#��<#�X<#�[<#�&<#�$<$9<$<#��<#��<#�a<#�.<#י<#׫<#�<#ޯ<#ە<#ޠ<#�~<#��<#�a<#�}<#�<#�<#��<#�<<#�d<#�'<#�<<#�
<#ט<#�D<#�<#�y<#�[<#�O<$J<#�S<#��<#�9<#�h<#�
<#�
<#�<#�6<#��<#�
<#�w<#�
<#�<#ג<#�<#�<#ڟ<#�<#�<#�<#ח<#ؙ<#ז<#�4<#ױ<#��<#׊<#�<#�0<#��<#�R<$�<$��<#�z<#�<#�<#�|<#�<#�B<$}<#�<#�z<#�<#׵<#�<#��<#�m<#�
<#�
<#�<#�<#�<#�
<#�
<#�<#�<#�d<#�<#׃<#׆<#�~<#�<#��<#��<#מ<#�<#�S<#�<#�7<#��<#�<#�d<#�
<#�<#��<#��<#�<#�\<#��<#ח<#�
<#�<#��<#߄<#�<#گ<#�"<#��<#�<#��<#י<#�t<#�<#޾<#�<<#׳<#�.<#�a<#ג<#׏<#��<#�c<#�c<#�<#�Z<#�d<#�<#ڀ<#�E<#��<#��<#�\<#ف<#��<$�<#׮<#� <#�<#�2<#׫<#��<#��<#ۇ<#�/<#ު<#�<#�V<#׋<#؂<#��<#�<#�<#��<#��<#�<#ؕ<#�^<#��<#�L<#�<#�<#�g<#�K<#ם<#�<#�v<#�h<#�E<#�|<$�<#��<#׬<#�<#�z<#ת<#�<#�J<#�<#�(<#ז<#׆<#�<#�<#�<#�c<#�<#�s<#�h<#؜<#�n<#׆<#�v<#�<#��<#�<#�f<#��<#�<#�~<#�<#�<#�"<#ב<#��<#�D<#�\<#ן<#��<#�<#�<#�<#�4<#�e<#��<#�<#�<#ׇ<#א<#�<#�S<#�2<#׍<#�
<#�
<#�n<#�n<#�
<#�<#�;<$u<#�<#׈<#�<#ך<#އ<#�<#�3<#ה<#�y<#�a<#ޣ<#�<#�p<#�[<#�<#�W<#�$<#��<#�=<#�_<#�j<#ۄ<#ޫ<#މ<#��<#��<#�<#ػ<#�<#�<#��<#�<#�<#�<#�N<#�<#��<#�Y<#��<#��<#�B<#׍<#ׇ<#ׇ<#א<#��<#�I<#�K<#�<#�V<#�=<#��<#�<#�o<#�Y<#�A<#�O<#��<#י<#� <#�<#�<#�
<#�o<#�w<#�<#��<#�<#� <$�<#��<#�x<#��<#�<#�{<#�
<#�
<#�
<#ׇ<#ׇ<#�a<#�<#י<#��<#�q<#׏<#��<#�~<#�<#�N<#��<#ރ<#ۆ<#�D<#��<#خ<#��<#�<#�G<#� <#ׅ<#ޟ<#�d<#�<#�]<#ׇ<#�i<#��<#�<#�{<#�
<#�<#�<#�c<#ע<#�<#��<#׊<#�x<#�<#��<#ס<#��<#�
<#�
<#�
<#�
<#�
<#�
<#��<#�O<#א<#�
<#�
<#�
<#�<#��<#�<#��<#�%<#ך<#י<#�C<#��<#�<#�
<#�
<#�<#�<#�<#��<#��<#�<#�Z<#��<#�
<#�
<#�<#ע<#�g<#ך<#ס<#��<#�
<#�<#�<#ة<#�f<#�}<#�<#��<#�<#��<#�
<#�w<#�
<#�
<#�~<#�<#�v<#�
<#�
<#�
<#�
<#�
<#�x<#�
<#�x<#�
<#�<#�~<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�v<#�
<#�
<#�
<#�
<#�
<#�M<#�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929144645201909291446452019092914464920190929144655202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@@  @@  @@  G�O�G�O�G�O�G�O�D�fD�fD�fG�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          