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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  J�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  VP   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  a�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  d�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  j�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  v8   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �h   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �P   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �8   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �P   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �H   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �H   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �H   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �H   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ˜   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ˸   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �0   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �L   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        μ   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210225041001  20210225041001  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               xA   BO  125408                          2C  D   APEX                            6229                            120210                          846 @�s�6�`1   @�s�6�`@Q7�vȴ9�4j=p��
1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @�ff@�  A   A   A@  A^ffA�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�|�D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�C3D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D���By�Bz�B{�B|�B�B�VB��B��B	�B	E�B	u�B	��B	�3B	��B	��B	�oB	��B	��B	�'B	�
B	�BB	�B	�)B	��B
0!B
XB
t�B
�JB
�LB
�B
��B  BoB(�B8RBG�BT�B^5BcTBl�Bx�B�B�%B�JB��B��B��B��B��B�B�B�-B�'B�'B�-B�9B�?B�3B�?B�RB�RB�jB��BŢBĜB��B�dB�RB�RB�RB�XB�dB�jB�^B�XB�LB�?B�3B�-B�!B�B�B�B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�uB�oB�oB�oB�hB�hB�hB�hB�\B�PB�PB�PB�PB�PB�PB�\B�bB�oB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�oB�uB�uB�uB�uB�uB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�1?�j?�/?�p�?���?��\?��9?���?�=q?θR?�o@��@��@S�?�?�{?��?˥�?��?��;@
��@�^?�E�?��7?�  ?y��?y�?}�?�33?�x�?�p�?�hs?��?���?��?���?�-?�b?�"�?�;d?�z�?�dZ?�I�?�;d?��y?���?�~�?��?�5??��?�  ?�33?���?���?�1'?��#?�^?�-?�&�?�?�J?�E�?�X?�R?���?陚?ޗ�?�ȴ?�33?�o?�A�?�J?ӶF?���?Η�?�I�?��y?�n�?�v�?�dZ?��y?��!?���?���?�b?�b?��P?��?�A�?��?�Ĝ?�Ĝ?���?��?�V?�5??�1?�ȴ?�z�?�9X?��F?�-?��?���?��?��j?�n�?��`?���?��?|�?z�H?z�H?z��?x�u?t��?v?u?r�!?o�?m��?mV?l��?l��?nV?l��?hr�?fff?e`B?d�/?d�?d�?c��?b��?a��?`A�?_�w?^�R?]�-?\j?Z��?WK�?T9X?Rn�?Q&�?N��?MO�?LI�?K?I��?H1'?G+?E`B?C�
?BJ?@  ?<(�?7�P?5?2n�?0��?/��?.��?,��?&��? Ĝ?/?^5?��?Q�?Q�?��?ȴ?�+?�j?�F?�F?�?n�?��?&�?��?{?�D?
~�?��?��?�7>��>�Q�>�ȴ>�>�->>�{>�h>�~�>��y>�l�>�ff>�G�>޸R>ݲ->�"�>Ձ>��>���>�C�>��>���>�"�>���>�t�>�O�>��>~��>�%>}�>|�>x��>vȴ>u>t�j>t�j>vȴ>vȴ>t�j>n��>l�D>j~�>hr�>fff>`A�>Xb>Q�>J��><j>A�7>I�^>O�;>W
=>Y�>Y�>S��>L��>J��>J��>I�^>?|�>+>��>�>�u>\)=��m=�;d=�"�=���=��`=��=��
=��w=���=�O�=�C�=�t�=H�9<�C�;�o:�o�D���D���T����`B;�o=\)<��=�w=e`B=���=��=�;d=�"�=�"�=��`=��T=�t�=�+=H�9='�=�P<�=P�`=�7L=�7L=�+=m�h=D��<��
<�`B=8Q�;o;��
<��
;o�49X����P��P�#�
�<j�Y��,1��/���
��h�o����`B���ͼ�j��9X��1��1����D���o���
�o    ;�o;��
;��
;��
;ě�;�`B<o;ě�;�o;D��;�o;ě�;ě�<#�
<D��<T��<e`B<u<�t�<���<��
<�1<�9X<ě�<ě�<���<���<���<�/<�`B<�h<�h<�<�<�<��<��<��<��<�<�<�h<�`B<�/<���<���<���<ě�<�9X<��
<�C�<�o<u<u<u<e`B<49X<t�;��
;�o;D��    ��o�D����o���
���
��`B��`B�t��#�
�D���e`B�u��o��t����
���ͼ��+�C���P�#�
�0 Ž<j�@��D���P�`�Y��aG��ixսq���u�}󶽇+��\)���������P���-���
���罬1��1�� Ž�j�Ƨ�ȴ9�����
=��;d��G���S���l��������ٽ��#�   �J���+�
=q�O߾V�V�V�z��u��������-��R��R��R�#�
�$�/�&�y�(�þ/��333�6E��9X�<j�=p��?|�B�\�D���I�^�M��P�`�Q녾T���W
=�Xb�Z��\(��\(��^5?�`A��_;d�_;d�aG��cS��fff�hr��j~��l�D�n���p�׾r�!�t�j�t�j�vȴ�z�H��  ���7����������𾇮���9��=q��ƨ��O߾����\)���;��bN��hs��n����Ͼ��+���P�����������������㾟;d��MӾ�S����/����þ�xվ�xվ�xվ�xվ�xվ�~���V��{����� ž��׾�&龱&龲-���j��Q쾸����X��X���#��^5��dZ��푾�  ��%��%��%���7�\����š˾�+�ȴ9�ɺ^��=q��C���O߾����V��V��\)���;���`��hs��t����Ͼ�����
=�ؓu�����(��ݲ-��5?�߾w��G���MӾ�Z���T��ff���y���r���r���xվ����~���~�����1������ ž�!��F���j����ȴ���پ��پ��پ�Q��^5��p�������|� A�� Ĝ�Mӿo�����
�Z����˿ff�l����1'��9�	xտ	�^�
=q�
���ƨ�I���D��h�{�������\)�����׿�`�t��z�z��j��j�?}�E���P��P��u�������"ѿdZ����m�j�푿�-��R��ۿ�ۿ;d�|��w�   �!%�"�\�#o�#o�#o�#o�#o�"��"��"��"��"�\�"�\�"�\�"�\�"��"��"�\�"�\�"�\�"�\�"�\�"�\�"�\�"�\�"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�!���!���!���!���!���!���!�7�!�7�!�7�!�71111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@�G�A ��A ��A@��A_
=A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�B (�B(�B(�B(�B (�B((�B0(�B8(�B@(�BH(�BP(�BX(�B`(�Bh(�Bp(�Bx(�B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{C 
=C
=C
=C
=C
=C

=C
=C
=C
=C
=C
=C
=C
=C
=C
=C
=C 
=C"
=C$
=C&
=C(
=C*
=C,
=C.
=C0
=C2
=C4
=C6
=C8
=C:
=C<
=C>
=C@
=CB
=CD
=CF
=CH
=CJ
=CL
=CN
=CP
=CR
=CT
=CV
=CX
=CZ
=C\
=C^
=C`
=Cb
=Cd
=Cf
=Ch
=Cj
=Cl
=Cn
=Cp
=Cr
=Ct
=Cv
=Cx
=Cz
=C|
=C~
=C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C��C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�C�D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV�DV��DW�DW��DX�DX��DY�DY��DZ�DZ��D[�D[��D\�D\��D]�D]��D^�D^��D_�D_��D`�D`��Da�Da��Db�Db��Dc�Dc��Dd�Dd��De�De��Df�Df��Dg�Dg��Dh�Dh��Di�Di��Dj�Dj��Dk�Dk��Dl�Dl��Dm�Dm��Dn�Dn��Do�Do��Dp�Dp��Dq�Dq��Dr�Dr��Ds�Ds��Dt�Dt��Du�Du��Dv�Dv��Dw�Dw��Dx�Dx��Dy�Dy��Dz�Dz��D{�D{��D|�D|��D}�D}��D~�D~��D�D��D�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD�~D��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�D{D��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��HD��HD�HD�AHD��{D��By�Bz�B{�B|�B��B�(B��B�B	�B	8:B	h�B	��B	��B	��B	�MB	��B	��B	��B	��B	�NB	�#B	�B	�B	�?B
6�B
XBB
s�B
��B
�B
�B
�5B
��B�B%�B4�BC�BR�B\�Ba�Bj�BvRB��B�B��B��B�B�DB��B��B��B�B�]B��B�B��B�XB��B��B��B��B��B�8B��BŝBƠBĄB�SB��B�lB�LB��B��B��B��B�DB�FB��B��B�[B��B��B��B�9B��B�B�4B��B��B��B��B��B�B��B�3B�B��B��B��B�B�B��B��B��B��B�+B��B�qB��B��B�+B��B��B��B�MB��B��B��B��B��B�1B�B��B��B��B�DB��B�XB�,B�B�"B�B�B�*B�8B�QB�"B�;B�>B�JB�cB��B��B�aB�IB��B�NB�DB�RB�HB�jB�IB�gB�\B�fB�sB��B��B�_B��B�XB�2B�0B�[B�B�B��B�jB�(B��B��B��B�B��B�1B�B��B�B��B�B��B�"B�"B� B�6B� B�sB��B�HB�KB��B��B�B�B��B��B�B�B��B��B�5B��B��B��B�BB�SB��B��B�B��B�pB��B�B�B�DB�B�]B��B��B��B��B��B��B��B�pB��B��B��B��B��B��B��B��B��B��B��B�'B�LB�%B�@B�9B�wB��B��B��B��B��B��B�B��B�MB�yB��B��B�MB�$B��B��B��B��B�B��B��B��B�uB�DB��B��B��B�{B��B��B�\B�B��B��B��B�B��B�yB�PB�zB��B��B��B��B�B��B�ZB��B��B��B��B��B��B��B��B�B��B�.B��B�yB�bB��B�dB��B��B��B�xB��B��B��B��B��B�0B��B��B�{B�eB�eB�tB��B��B��B�xB�<B�gB�wB�{B��B�vB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B��B��B�B�B�B�'B�B�B�B�B�B�B�B��B�B��B�B�
B�B�B�	B�	B�B�B�:B�:B�"B�
B�"B�"B�"B�!B�	B�
B�!B�B�B�B�B�
B�B�.B�.B�!B��B�	B�"B�"B�!B�	B��B�B�EB�:B�
B�"B�.B�-B�
B�	B�B�B�8B�	B�
B�!B�B�!B�B�"B�!B�	B��B��B�AB�)B�)B��B��B�B�B��B��B�4B�B�B�B�@B�)B�B�B�B�B�B�B�B�4B�(B�B�B�B�B�B�B�B��B�B�B��B��B�B�B�B�B�B�B�B�B�B�B��B�
B�#B�.B�B�"B�B�B�
B�
B�B�B�B�B��B��B��B�	B�
B�B�'B�B�B��B��B��B�B�?B�3B�B�B�2B�B��B��B��B��B��B�B�&B�B�B�B��B��B��B�B�)B�>B��B��B��B��B��B��B�
B�-B��B��B��B��B��B�	B�	B�	B�	B��B��B��B�B��B��B��B��B��B��B��B�B��B��B�B�	B�	B�B�B��B�B�	B��B�B�	B��B��B��B��B��B��B��B��B��B��B��B�+B��B��B� B��B��B��B��B��B��B��B��B�B�,B��B��B�B��B��B�,B�
B��B��B��B�B��B�	B�B��B��B��B��B�B��B��B�B�B�B��B�B�B�B��B��B�B��B�B��B�_B�B��B��B��B�B�B�&B��B�B�B�B�B��B��B��B�B�B�B�B��B��B��B��B��B��B�B�3B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�1?�j?�/?�p�?���?��\?��9?���?�=q?θR?�o@��@��@S�?�?�{?��?˥�?��?��;@
��@�^?�E�?��7?�  ?y��?y�?}�?�33?�x�?�p�?�hs?��?���?��?���?�-?�b?�"�?�;d?�z�?�dZ?�I�?�;d?��y?���?�~�?��?�5??��?�  ?�33?���?���?�1'?��#?�^?�-?�&�?�?�J?�E�?�X?�R?���?陚?ޗ�?�ȴ?�33?�o?�A�?�J?ӶF?���?Η�?�I�?��y?�n�?�v�?�dZ?��y?��!?���?���?�b?�b?��P?��?�A�?��?�Ĝ?�Ĝ?���?��?�V?�5??�1?�ȴ?�z�?�9X?��F?�-?��?���?��?��j?�n�?��`?���?��?|�?z�H?z�H?z��?x�u?t��?v?u?r�!?o�?m��?mV?l��?l��?nV?l��?hr�?fff?e`B?d�/?d�?d�?c��?b��?a��?`A�?_�w?^�R?]�-?\j?Z��?WK�?T9X?Rn�?Q&�?N��?MO�?LI�?K?I��?H1'?G+?E`B?C�
?BJ?@  ?<(�?7�P?5?2n�?0��?/��?.��?,��?&��? Ĝ?/?^5?��?Q�?Q�?��?ȴ?�+?�j?�F?�F?�?n�?��?&�?��?{?�D?
~�?��?��?�7>��>�Q�>�ȴ>�>�->>�{>�h>�~�>��y>�l�>�ff>�G�>޸R>ݲ->�"�>Ձ>��>���>�C�>��>���>�"�>���>�t�>�O�>��>~��>�%>}�>|�>x��>vȴ>u>t�j>t�j>vȴ>vȴ>t�j>n��>l�D>j~�>hr�>fff>`A�>Xb>Q�>J��><j>A�7>I�^>O�;>W
=>Y�>Y�>S��>L��>J��>J��>I�^>?|�>+>��>�>�u>\)=��m=�;d=�"�=���=��`=��=��
=��w=���=�O�=�C�=�t�=H�9<�C�;�o:�o�D���D���T����`B;�o=\)<��=�w=e`B=���=��=�;d=�"�=�"�=��`=��T=�t�=�+=H�9='�=�P<�=P�`=�7L=�7L=�+=m�h=D��<��
<�`B=8Q�;o;��
<��
;o�49X����P��P�#�
�<j�Y��,1��/���
��h�o����`B���ͼ�j��9X��1��1����D���o���
�o    ;�o;��
;��
;��
;ě�;�`B<o;ě�;�o;D��;�o;ě�;ě�<#�
<D��<T��<e`B<u<�t�<���<��
<�1<�9X<ě�<ě�<���<���<���<�/<�`B<�h<�h<�<�<�<��<��<��<��<�<�<�h<�`B<�/<���<���<���<ě�<�9X<��
<�C�<�o<u<u<u<e`B<49X<t�;��
;�o;D��    ��o�D����o���
���
��`B��`B�t��#�
�D���e`B�u��o��t����
���ͼ��+�C���P�#�
�0 Ž<j�@��D���P�`�Y��aG��ixսq���u�}󶽇+��\)���������P���-���
���罬1��1�� Ž�j�Ƨ�ȴ9�����
=��;d��G���S���l��������ٽ��#�   �J���+�
=q�O߾V�V�V�z��u��������-��R��R��R�#�
�$�/�&�y�(�þ/��333�6E��9X�<j�=p��?|�B�\�D���I�^�M��P�`�Q녾T���W
=�Xb�Z��\(��\(��^5?�`A��_;d�_;d�aG��cS��fff�hr��j~��l�D�n���p�׾r�!�t�j�t�j�vȴ�z�H��  ���7����������𾇮���9��=q��ƨ��O߾����\)���;��bN��hs��n����Ͼ��+���P�����������������㾟;d��MӾ�S����/����þ�xվ�xվ�xվ�xվ�xվ�~���V��{����� ž��׾�&龱&龲-���j��Q쾸����X��X���#��^5��dZ��푾�  ��%��%��%���7�\����š˾�+�ȴ9�ɺ^��=q��C���O߾����V��V��\)���;���`��hs��t����Ͼ�����
=�ؓu�����(��ݲ-��5?�߾w��G���MӾ�Z���T��ff���y���r���r���xվ����~���~�����1������ ž�!��F���j����ȴ���پ��پ��پ�Q��^5��p�������|� A�� Ĝ�Mӿo�����
�Z����˿ff�l����1'��9�	xտ	�^�
=q�
���ƨ�I���D��h�{�������\)�����׿�`�t��z�z��j��j�?}�E���P��P��u�������"ѿdZ����m�j�푿�-��R��ۿ�ۿ;d�|��w�   �!%�"�\�#o�#o�#o�#o�#o�"��"��"��"��"�\�"�\�"�\�"�\�"��"��"�\�"�\�"�\�"�\�"�\�"�\�"�\�"�\�"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"Mӿ"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�"J�!���!���!���!���!���!���!�7�!�7�!�7�!�71111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<#��<#ر<#�<%��<'��<1�Q<Kh<OTU<�I<��<0�a<a��<r��<.l<2�<'��<*Y�<��<M^M<ET+=
]<��~<3�<A��<#�\<$�9<%�<'�p<&��<%��<+�~<3��<+��<,n(<.��<'m!<$�Y<%��<&��<(�<#�J<$��<)?�<$j�<$�<#��<$k<#�~<#��<$��<'�<$l8<$�R<$�<#�6<)/�<$ +<%��<%�<%eJ<$��<&��<#�<'�<0�<*�6<%Ht<#�<$��<$l<$*<#��<%��<$��<&�p<&=<%��<$��<&)<%�&<&u;<%�<$��<#��<#��<$��<&I=<#�<#��<#�/<#؀<$+�<#�v<#�1<$m�<&�N<$�<#�/<#��<$'z<' <&�|<$t4<#� <$j�<$�<#ڃ<#�e<#��<$f�<#��<#�z<#�6<$/i<#��<#�D<$!�<$5d<#�<#��<#�8<#�<#�U<#�r<$[�<#�<#�<#�><#�<#�+<#�<#ݕ<#��<#��<#�,<#�#<#��<#�<#�<$+�<$Y<#�}<#�<$	�<#��<#��<#�<#�*<#�<#�<#�/<#�<#�<#��<$Jr<$m�<#�!<$)�<#��<#�2<#�><#�T<$�<$�b<$@�<$�<#�N<#�'<#�<<#�<#�G<#؎<#�<#��<#�?<#ܢ<#�.<#ܫ<#�<#�6<#�;<#�8<#�><#�<$,�<$TL<$�<$�<#��<#�<#�W<#�q<#�?<#�c<#�w<#�<#�<#ڼ<$�<#�&<#ڌ<#�v<$<$&<#�d<#�<%Yv<%3�<&��<#��<$�<$$�<$Q�<$�<#��<#�W<#�$<#�	<#ڡ<#��<#�<#�$<#��<#�,<#�b<#��<#�b<#�<#�<#�g<#�	<#��<#��<#��<$0�<#ݼ<#�	<#�v<#�<#ؖ<#�8<#�M<#��<#�*<#�8<#�B<$f<$�p<$UW<#�A<#�I<$(<$eR<$;<#�<#��<#ڥ<#�U<$6�<#��<#�J<#��<#�)<#��<$�S<%��<$7�<#�V<#�w<$ <#�<#�0<$"<%i�<#ܗ<#�C<$d�<'	,<#�<#�7<#�<#�P<#��<$�><$'<#�(<$^C<#��<#�<#�A<$�<$B`<#�
<#��<#��<$�<%;<#��<$:S<&��<#��<$+�<$s:<'^�<$�<#�$<#�P<#��<#��<#��<$	�<$9�<#�,<#�@<#�v<#�'<#�f<#ڄ<#�<#�,<#�<#�<#�V<#��<#� <#��<#�<#� <#܋<#�U<#�$<#�$<#�*<#�*<#�<#�Q<#ٺ<#�<#�<#�W<#�<#��<#�i<#�<#�K<#�4<#٤<#�V<#�+<#�+<#�/<#�><#�<#�'<#�&<#�-<#�E<#�*<#�&<#�$<#�9<#�#<#�$<#�<#�<#�(<#�,<#�C<#�,<#�<#�<#�m<#��<#�.<#�<#�)<#�^<#�<#ܖ<#�4<#�<#�)<#�,<#�*<#ܐ<#�?<#�J<#�7<#�5<#�w<#�2<#��<#�&<#�h<#� <#��<#�7<#��<#�4<#��<#��<#�'<#�(<#��<#�H<#�(<#�*<#��<#�@<#��<#��<#��<#ܪ<#�<#�5<#܎<#�<#�<#�<#��<#�6<#�"<#�<#��<#ܩ<#�9<#�0<#��<#��<#ܨ<#� <#�1<#�?<#�;<#��<#�K<#��<#��<#�J<#�6<#�(<#��<#�E<#�V<#�0<#�5<#ܔ<#�0<#ܴ<#�<#��<#ܨ<#�$<#�*<#�7<#�<#�<#��<#�'<#�b<#��<#��<#�E<#�&<#�<#�]<#��<#�]<#�K<#��<#�<#��<#ܪ<#�?<#�<#ܶ<#�G<#��<#�<#��<#�?<#ܯ<#�<#�4<#��<#��<#�7<#��<#��<#�"<#�.<#��<#�<#ܳ<#�<#�<#�<#�<#�<#��<#�[<#�:<#�<#��<#�<#�@<#�v<#��<#��<#�<#�<#��<#��<#��<#ܩ<#�1<#�<#�'<#��<#�<#܁<#�<#ڱ<#��<#�0<#��<#�M<#�`<#�<#�<#�<#�*<#��<#�:<#�<#�(<#�'<#�'<#�0<#�D<#�<#�D<#�<#��<#�<#�<#�5<#�)<#��<#�><#�i<#�<#�-<#�<#��<#ڛ<#�)<#�t<#�<#�/<#�+<#�!<#�<#��<#��<#��<#��<#��<#�4<#�<#�,<#�'<#�<#�4<#��<#�%<#��<#�M<#��<#�6<#�&<#�T<#�<#�<#��<#ܿ<#�N<#ܯ<#��<#�7<#�R<#��<#�+<#�%<#��<#�<#�3<#��<#�<#�<#�/<#�<#�X<#�<#�]<#�N<#�<#�<#�<#�<#�<#��<#�,<#�+<#�6<#��<#�<#�G<#�9<#܉<#�<#�Z<#��<#�<#��<#�+<#�<#ܪ<#�(<#�<#��<#�,<#�i<#��<#�H<#܎<#�?<#�Z<#�-<#ߐ<#�S<#��<#�E<#�'<#�p<#�4<#�/<#��<#�<#��<#�<$
�<#�<#�<<#��<#�3<#�/<#�<#�~<#�T<#��<#�<#�p<#�<#� <#�<#�%<#��<#�<#��<#�:<#�-<#�/<#�<#�<#�<#�A<#�<#�O<#��<#�N<#�
<#�<<#�"<#�&<#�"<#�&<#�-<#�A<#�><#�&<#�*<#�	<#�$<#�%<#�"<#�&<#�&<#�&<#�G<#�<#�"<#�&<#�"<#�&<#�&<#�&<#�&<#�%<#�%<#�%<#�<#�F<#�><#�%<#�%<#�%<#�%<#�%<#�%<#�%<#�G<#�<#�?<#�<#�C<#�%<#�;<#�%<#�!<#�%<#�"<#�%<#�%<#�%<#�!<#�&<#�%<#�%<#�%<#�%;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.04                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929144438201909291444382019092914444220190929144449202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�D���D���D���G�O�G�O�G�� G�� G�� G�� G�� G�� G��                                 6389758         6389758         131072                                          