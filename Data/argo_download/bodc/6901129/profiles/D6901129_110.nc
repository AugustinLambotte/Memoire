CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  -   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  K�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  X�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  e<   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  hl   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  k�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  n�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  {�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �4   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 0  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 0  �   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 0  �H   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  �x   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �,   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  Ĕ   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �H   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �d   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ׀   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ל   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ׸   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �x   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �h   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ڄ   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ڠ   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ڼArgo profile    3.1 1.2 19500101000000  20210225034743  20210225034743  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               nA   BO  125398                          2C  D   APEX                            6229                            120210                          846 @�[7��`1   @�[7��`@QÕ�$��2l1&�y1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @9��@�33@�ffA   A  A>ffA`  A�  A���A�ffA�33A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  Dy�D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D D�� D�  D�@ DÀ D�� D�  D�@ DĀ D�� D�  D�@ Dŀ D�� D�  D�@ Dƀ D�� D�  D�@ Dǀ D�� D�  D�@ DȀ D�� D�  D�@ Dɀ D�� D�  D�@ Dʀ D�� D��D�  BR�BT�Bq�B�B�`B��B�HB/BcTB�%B	B	m�B	�^B	��B
$�B
\)B
�1B
��B
�FB
��B
�`B
�yB
�B
��B
��B
��BBB+B	7BDBJB\BbBhBoBuB{B�B�B�B�B�B�B �B"�B&�B(�B)�B,B-B-B/B33B49B49B6FB;dB>wB@�BA�BD�BF�BF�BH�BL�BM�BQ�BZBaHBl�Bu�Bz�B~�B~�B|�Bw�Br�Bn�BjBffBdZBcTBdZBdZBcTBcTBbNBbNBbNBaHB^5B]/B]/B]/B^5B_;B`BBbNBdZBe`BffBhsBjBk�Bk�Bk�Bl�Bm�Bo�Bq�Br�Bs�Bu�Bv�Bx�B{�B~�B� B�B�B�B�B�+B�7B�DB�PB�\B�\B�bB�bB�hB�hB�hB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�{B�uB�uB�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�bB�bB�bB�bB�bB�\B�\B�\B�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�\B�\B�bB�bB�hB�oB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�9?�7L?��?�7L?�j?�ƨ?���?n{?\j?\�?;��? Ĝ?l�? Ĝ>�z�>�K�>���>�A�>��
>��P>t�j>gl�>Xb>Q�>O�;>J��>(��=�x�=��=���=�\)=]/=�P<��<�1<o��o�ě������`B�ě������J���"��n��I��333�^5?��J���˾�1'��I���S����r���xվ�����������������m�h��$ݾ���r�!�49X�}�$�/�ȴ9<���>$�/>["�>�>�/>��R>��+>!��=���=0 ŽP�`����V�I��	7L�z��+�&�y�(�þ)��.{�R�vȴ��I����;��\)��I���+�|푾hr��e`B�cS��Q녾-V�#�
�#�
�#�
�'�R�C���xս��`�\��1��hs��+��`B<t�=+=�P=\)=\)<�/=]/=ȴ9=�;d>�>&�y>333>5?}>9X>;dZ><j>:^5>;dZ>;dZ>:^5>:^5>:^5>9X>5?}>.{>)��>+>'�>%�T>!��>!��>"��> Ĝ>��>�>bN>C�>%=�====�F=�=�h=�`B=�l�=�/=ȴ9=�j=�E�=�-=�-=�-=�1=��w=���=���=�C�=}�=m�h=e`B=aG�=]/=H�9=\)<�h<��
<D��<49X;ě�:�o�o���
��/���C��C��+�\)�<j�e`B�q���q���m�h�ixսe`B�aG��]/�T���H�9�H�9�49X��w����j�D��;o<���<��=,1=aG�=}�=��=�hs=��
=� �=�Q�=��=ȴ9=���=�
==�
==ě�=�^5=Ƨ�=ȴ9=���=�j=�1=��-=���=��T=�-=�E�=�-=��w=�hs=�C�=��=u=m�h=q��=aG�=�7L=��=�C�=aG�=P�`=D��=@�=�+=��-=���=��P=q��=u=q��=q��=u=u=y�#=q��=L��=@�=H�9=T��=aG�=y�#=Y�=P�`=D��=<j=<j=49X=,1=��=\)<��<��<��<��<�<�/<�9X<���<u<e`B<D��<49X<t�<o<o;�`B;�`B;��
;�o;D��;o;o;��
;�`B;ě�;��
;�o;D���D����`B��`B�ě��ě��ě��o�u��t���C����
��j���ͼ�/��`B��`B��`B��`B���+�C��t���P����w�#�
�,1�0 Ž<j�H�9�D���T���]/�u�u�y�#�u�y�#��%��O߽�\)���P���-������1��{�� Ž�Q������`������
=��
=��
=�����"ѽ�;d��S����\)�n��bN�bN�V�C��I��C��I��\)�t���������u��������R� Ĝ�#�
�$�/�%�T�%�T�(�þ)��+�,1�,1�,1�-V�.{�1&�333�9X�<j�?|�@��C���G��I�^�J���L�;N��S�ϾY��\(��]/�_;d�cS��fff�j~��q���t�j�x���y�#�y�#�y�#�|푾~�۾�%���\������˾�$ݾ���+�����������=q��C���I���O߾�����;��n���t���zᾔ������������zᾔ�����������+�����������㾜���5?��A��������徢�徣S����
��Z��`B��ff���y���r���r����þ�xվ�xվ���1��V��V��{����� ž��׾�&龲-���F��9X���j����E���ȴ��ȴ��KǾ������#��^5��^5���H��dZ��dZ��dZ���m��j��푾�vɾ�  ��  ���7�ě���+��1'��1'��7L��=q������C���I����;�V������;��bN��녾�n���z�Ձ����
=��b�ؓu����������ڟ���"Ѿ�"Ѿۥ��(���(��ܬ��5?��A�����������MӾ�S����
���/��`B��`B��ff���y��l���r���r���r����þ���V��V��{������� ž�&��������-��!��F��9X��9X��?}����ȴ��KǾ�KǾ��پ�Q��Q�����������#��j��푾�p�������vɾ��۾�vɾ��۾�|� ��%�G��%�G������7�J�J�J�Mӿ��o�����
�Z��/�����`B��˿�T�$ݿ����l����1'��9�	7L�	xտ	��	��	��	��
=q�
~��
���
����C��C��C��ƨ�I���D��Ϳ�ͿV�O߿O߿�h�����{������\)�����;� ſbN��׿�`��`�&�&�&�hs����녿n��n��n���33�t��t���F��Ͽ9X�z��j��������+�ȴ�ȴ�
=�
=�
=�Kǿ�P�b��u��u��u��������������#��#�����������^5����"ѿ"ѿdZ����m����m��m��m�(��j�����푿/�p��5?�vɿ�R��ۿ�ۿ�ۿ;d�;d�|��w��w��w�   �   �   � A�� A�� A�� A�� A�� A�� A�� A�� �� �� �� �� �� �� �� Ĝ� Ĝ�!%�!G��!�7�"J�"Mӿ"��"��#o�#o�#���$��$�/�%��%`B�%`B�%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%`B�%`B�%`B�%`B�%`B�%`B�%`B�%`B�%��%��%��%��%��%��%��%��%��%��%��$�/�$�/�$�/�$�/�$�/111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @9��@�33@�ffA   A  A>ffA`  A�  A���A�ffA�33A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  Dy�D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D D�� D�  D�@ DÀ D�� D�  D�@ DĀ D�� D�  D�@ Dŀ D�� D�  D�@ Dƀ D�� D�  D�@ Dǀ D�� D�  D�@ DȀ D�� D�  D�@ Dɀ D�� D�  D�@ Dʀ D�� D��D�  BR�BS�Bs�B{B��B�!B�B2-Bc�B��B	�B	r(B	�QB	ԷB
'nB
^B
�TB
�{B
�}B
τB
��B
�=B
��B
��B
�GB �BvB�B$B	�BB"B�B�BB�B�BB�B�B>B B�BB B"�B(�B+B+�B,kB-OB-�B1@B3�B4GB4GB4�B:YB?zB@xB?�BCBG�BGpBF�BJBP�BM�BV�B[!BfBr�Bw"B~?B~�B}�B~1BviBp�BoBi�BeoBcABd:Bd�Bc{BdBbhBb]Bb�Bb�B_�B^�B]�B]'B]�B^�B_nBaZBd+Be@Be�Bf�BjBk�Bk�Bk�BlBl�Bn�BqBrTBs/Bu Bv{BwBzB}�B�B�B�B�^B��B�B��B�FB��B��B�@B�1B�LB�[B��B�_B�lB�{B�qB�uB��B��B��B��B�lB��B��B��B�wB�mB��B��B��B��B��B��B��B��B�vB�tB�]B�kB��B��B�mB��B��B��B��B��B�oB�qB��B��B��B�B��B��B��B��B�|B�~B��B�B��B��B��B�sB��B��B��B��B��B�}B��B�XB�KB�tB��B��B�xB�PB�DB�DB�DB�FB�IB�=B�5B�VB�B� B��B�B��B��B��B��B�B��B�IB�}B�aB�CB�mB��B��B��B�vB��B��B�-B��B�uB��B��B�B�)B�B��B��B�B��B��B�6B�B��B��B�B��B��B��B�:B��B�	B�jB��B��B��B��B�LB��B��B��B��B��B��B��B��B��B��B�IB�B��B��B��B��B�GB�B�B�B��B�B�B�B�B�B��B��B��B��B�B�'B�B�B��B�	B��B�B��B��B��B��B�B��B��B��B��B��B��B�B�	B�	B�B�GB�-B��B��B��B��B�B�PB�!B��B�!B�"B�B�B�B��B��B��B�B�!B�
B�B�	B�	B�	B�	B�B�
B�!B� B��B�2B�B�IB�B�B��B�B�B�JB�B�4B�(B�B�>B�B�B�0B�3B�_B�B�B�B�B�B�B�B�B��B��B�B��B��B��B��B�B��B�	B�"B�-B�B��B��B�!B�!B�	B��B�B�B�!B�	B�B��B� B�
B�	B�B��B��B�	B�
B�!B�B�EB�"B�!B�
B�"B�-B�B�B�B�B�5B�4B�B�B�B�"B�B�#B�CB�B�B��B��B��B�B�
B�B�B�B�B��B��B��B��B��B��B�,B�
B�	B�
B�B�B�.B�
B�	B��B��B��B��B��B��B��B�B�8B�B�B�B�B�B�B�B��B��B��B��B�B�B��B�B��B��B��B��B��B�B�B�B��B�B�B�B��B��B�B�B��B��B�B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B�B�B��B�B�0B�!B��B��B�B��B��B��B��B��B�B��B��B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B�B�B�B�B�B�B��B�B��B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�9?�7L?��?�7L?�j?�ƨ?���?n{?\j?\�?;��? Ĝ?l�? Ĝ>�z�>�K�>���>�A�>��
>��P>t�j>gl�>Xb>Q�>O�;>J��>(��=�x�=��=���=�\)=]/=�P<��<�1<o��o�ě������`B�ě������J���"��n��I��333�^5?��J���˾�1'��I���S����r���xվ�����������������m�h��$ݾ���r�!�49X�}�$�/�ȴ9<���>$�/>["�>�>�/>��R>��+>!��=���=0 ŽP�`����V�I��	7L�z��+�&�y�(�þ)��.{�R�vȴ��I����;��\)��I���+�|푾hr��e`B�cS��Q녾-V�#�
�#�
�#�
�'�R�C���xս��`�\��1��hs��+��`B<t�=+=�P=\)=\)<�/=]/=ȴ9=�;d>�>&�y>333>5?}>9X>;dZ><j>:^5>;dZ>;dZ>:^5>:^5>:^5>9X>5?}>.{>)��>+>'�>%�T>!��>!��>"��> Ĝ>��>�>bN>C�>%=�====�F=�=�h=�`B=�l�=�/=ȴ9=�j=�E�=�-=�-=�-=�1=��w=���=���=�C�=}�=m�h=e`B=aG�=]/=H�9=\)<�h<��
<D��<49X;ě�:�o�o���
��/���C��C��+�\)�<j�e`B�q���q���m�h�ixսe`B�aG��]/�T���H�9�H�9�49X��w����j�D��;o<���<��=,1=aG�=}�=��=�hs=��
=� �=�Q�=��=ȴ9=���=�
==�
==ě�=�^5=Ƨ�=ȴ9=���=�j=�1=��-=���=��T=�-=�E�=�-=��w=�hs=�C�=��=u=m�h=q��=aG�=�7L=��=�C�=aG�=P�`=D��=@�=�+=��-=���=��P=q��=u=q��=q��=u=u=y�#=q��=L��=@�=H�9=T��=aG�=y�#=Y�=P�`=D��=<j=<j=49X=,1=��=\)<��<��<��<��<�<�/<�9X<���<u<e`B<D��<49X<t�<o<o;�`B;�`B;��
;�o;D��;o;o;��
;�`B;ě�;��
;�o;D���D����`B��`B�ě��ě��ě��o�u��t���C����
��j���ͼ�/��`B��`B��`B��`B���+�C��t���P����w�#�
�,1�0 Ž<j�H�9�D���T���]/�u�u�y�#�u�y�#��%��O߽�\)���P���-������1��{�� Ž�Q������`������
=��
=��
=�����"ѽ�;d��S����\)�n��bN�bN�V�C��I��C��I��\)�t���������u��������R� Ĝ�#�
�$�/�%�T�%�T�(�þ)��+�,1�,1�,1�-V�.{�1&�333�9X�<j�?|�@��C���G��I�^�J���L�;N��S�ϾY��\(��]/�_;d�cS��fff�j~��q���t�j�x���y�#�y�#�y�#�|푾~�۾�%���\������˾�$ݾ���+�����������=q��C���I���O߾�����;��n���t���zᾔ������������zᾔ�����������+�����������㾜���5?��A��������徢�徣S����
��Z��`B��ff���y���r���r����þ�xվ�xվ���1��V��V��{����� ž��׾�&龲-���F��9X���j����E���ȴ��ȴ��KǾ������#��^5��^5���H��dZ��dZ��dZ���m��j��푾�vɾ�  ��  ���7�ě���+��1'��1'��7L��=q������C���I����;�V������;��bN��녾�n���z�Ձ����
=��b�ؓu����������ڟ���"Ѿ�"Ѿۥ��(���(��ܬ��5?��A�����������MӾ�S����
���/��`B��`B��ff���y��l���r���r���r����þ���V��V��{������� ž�&��������-��!��F��9X��9X��?}����ȴ��KǾ�KǾ��پ�Q��Q�����������#��j��푾�p�������vɾ��۾�vɾ��۾�|� ��%�G��%�G������7�J�J�J�Mӿ��o�����
�Z��/�����`B��˿�T�$ݿ����l����1'��9�	7L�	xտ	��	��	��	��
=q�
~��
���
����C��C��C��ƨ�I���D��Ϳ�ͿV�O߿O߿�h�����{������\)�����;� ſbN��׿�`��`�&�&�&�hs����녿n��n��n���33�t��t���F��Ͽ9X�z��j��������+�ȴ�ȴ�
=�
=�
=�Kǿ�P�b��u��u��u��������������#��#�����������^5����"ѿ"ѿdZ����m����m��m��m�(��j�����푿/�p��5?�vɿ�R��ۿ�ۿ�ۿ;d�;d�|��w��w��w�   �   �   � A�� A�� A�� A�� A�� A�� A�� A�� �� �� �� �� �� �� �� Ĝ� Ĝ�!%�!G��!�7�"J�"Mӿ"��"��#o�#o�#���$��$�/�%��%`B�%`B�%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%�˿%`B�%`B�%`B�%`B�%`B�%`B�%`B�%`B�%��%��%��%��%��%��%��%��%��%��%��$�/�$�/�$�/�$�/�$�/111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<%w�<'�<OpV<�-<:�<ax<*��<#�w<B�<4I<3J�<$�@<5�<(��<&��<#�<#��<$�<)uZ<$G<$J<#�o<#؇<#�<%�r<(P|<$ �<$��<$+�<$O <$c�<#�<$�<$&�<#��<#�G<$'N<$ �<#��<)P<(3<$��<#�<$B/<#�)<&|�<'!-<&Y<#�|<#�<$`<'a.<$]<#ג<#ג<%rJ<$�c<$��<#�u<&l�<%�#<%7�<$P<' <)��<+?�</�<+��<?$)<BKr<*&<.vV<$A�<#�<$��<A
�<.R�<'%�<3%�<,$b<$��<#�,<#�S<$�<#ۆ<$ER<#��<#׷<#�j<&^<&0<%�w<#��<#�E<#��<$
�<$a<$��<#��<#�L<$j�<&�<$x<#�<#�<#�Y<#�<$�<$�<$,<#�!<$�<$(<#�<&n�<&^+<$�H<#�<#؎<#�+<#�Q<%%j<'K�<$�<$�W<%ݟ<${<#�k<#�[<#ؙ<#ל<#��<#�L<#�<#�v<#�<#�<#כ<#��<#�<#ݑ<#�M<#��<#ٞ<#�{<#�<#�><#ؠ<#�<#�,<#�n<#�8<$�<#�2<#�8<#�<#�<#��<#�^<#�S<#�(<#�;<#�<$�<#��<#ۚ<#��<#�
<#�<#�i<#�,<#�	<#��<#�X<#�M<#޹<#��<#׎<#װ<#�<$0�<#�d<#�;<#�<#��<#�;<#��<#�"<$�<#�<#ۮ<#�*<#�<#�f<#٪<$�<$�<#ۿ<#�
<#ׄ<#ׇ<#ׇ<#�b<#א<#�<#�k<#�'<#�x<#�`<#�<#�_<$
<$7�<$_.<$!P<$�<$%�<#�n<#�<#�<#��<#�<#�]<#ߺ<#�<#��<#פ<#�<#�><#��<#�u<#��<#�I<#�X<#��<#�<#�|<#�A<#�3<#ز<#�I<#��<#�<#ے<#ۚ<#�U<#��<#�T<#��<$�<#�Y<#��<$#|<#� <#�6<#�<$t�<$�<#�$<#�'<$5�<#�I<#�w<#�<#ײ<#�<#�k<#�S<#�T<#�j<#ؾ<#ۤ<#�<#�h<#��<#�<#�8<#��<#�<#��<#�<#�k<#ۗ<#�(<#�
<#�
<#�
<#ב<#ۅ<#�<#�S<#�j<#ע<#��<#ו<#��<#א<#�<#�{<#�<#��<#׌<#׆<#צ<#�<#�<#آ<#�g<#׆<#׈<#׶<#�<#�l<#�<#�v<#�
<#�<#�3<#�<#�H<#�e<#�<#�L<#�	<#��<#ׁ<#�
<#�
<#�<#�<#�!<#ס<#��<#׊<#׆<#ׇ<#׏<#��<#י<#�?<#�<#�(<#��<#�<#��<#�
<#�p<#�u<#ׄ<#�/<#��<#׵<#�~<#�f<#�A<#�"<#ׇ<#ע<#ݒ<#�/<#��<#��<#ׅ<#�
<#�
<#ׇ<#�i<#�l<#�:<$j<$�<#۞<#؋<#�<#�
<#�e<#�X<#�g<#ג<#�l<#�f<#��<#�<#�<#�1<#�2<#׎<#�<#��<#�<#�<#ג<#�<#�<#��<#ז<#׆<#�~<#�
<#�
<#ׁ<#ט<#�-<#�U<#�-<#ۅ<#�5<#ת<#�N<#�n<#��<#�o<#�d<#�2<#��<#��<#�b<#ף<#��<#ެ<#ۘ<#�<#��<#��<#�W<#��<#�<#�<#�O<#�<#�O<#�i<#�g<#�4<#ז<#ׇ<#׆<#�~<#�
<#�<#��<#�<#��<#�<#�8<#�=<#�<#� <#��<#ׅ<#�
<#�
<#�n<#�v<#�~<#�<#�)<#��<#�;<#��<#�><#�l<#މ<#�_<#��<#�<#�~<#ׇ<#א<#��<#��<#ט<#��<#�~<#�<#�~<#�<#�<#�1<#�
<#��<#�<#��<#��<#��<#׏<#א<#�<#�<#י<#׏<#��<#׌<#�~<#�<#ו<#�:<#��<#ׇ<#�<#׆<#�~<#�
<#�
<#׃<#ׇ<#י<#�4<#�<#�<#�<#�
<#��<#�<#�<#ٺ<#�#<#׹<#׏<#�P<#�;<#�z<#�<#�<#ס<#�<#װ<#މ<#�<#י<#��<#��<#׊<#�~<#�<#��<#׆<#�~<#�<#�<#�~<#�<#ב<#�P<#ޑ<#�<#�<#ׇ<#��<#ד<#��<#׆<#�<#��<#ׇ<#׏<#��<#�<#�
<#ר<#�<<#�O<#�<#��<#��<#ׇ<#�<#׌<#��<#׆<#�<#ׁ<#׏<#��<#ׇ<#�<#��<#ד<#��<#ׇ<#�<#׆<#�~<#�<#�w<#�<#�"<#�<#׫<#׆<#�~<#�<#�~<#�v<#�o<#׀<#י<#�<#��<#�v<#�j<#�~<#ػ<#�e<#��<#�<#�
<#׉<#��<#ז<#��<#ט<#��<#��<#׆<#�<#ׁ<#ׇ<#ׇ<#׏<#��<#�<#�<#ז<#ׇ<#א<#��<#��<#ט<#��<#�
<#�
<#�
<#�~<#׆<#�~<#�<#�~<#�~<#�
<#�<#��<#��<#׌<#�~<#�<#�~<#�~<#�<#׃<#�~<#�<#׉<#��<#��<#׏<#׆<#׆<#׆<#׆<#��<#�=<#�
<#�V<#�
<#�
<#׃<#�q<#��<#��<#�<#�<#��<#׏<#�~<#�<#׃<#׆<#׆<#׆<#ׇ<#׏<#�a<#ؙ<#��<#ׁ<#�<#�v<#�
<#�
<#ׁ<#׏<#��<#��<#�
<#�
<#�~<#�~<#�
<#�<#��<#ׅ<#�<#�z<#�
<#�
<#�
<#�
<#ׁ<#׏<#��<#�<#׃<#׆<#ײ<#׿<#�~<#�
<#�<#�R<#׆<#�~<#�<#�~<#ׇ<#ט<#�<#ט<#׆<#�~<#�
<#�
<#�v<#�<#׆<#�~<#�
<#�
<#�z<#�
<#�
<#�u<#�
<#�
<#�<#�
<#�
<#�
<#�
<#�v<#�
<#�
<#�
<#�
<#�
<#�
<#�{<#�<#�<#׆<#׏<#��<#ו<#��<#�<#�{<#�
<#�F<#�<#� <#׆<#�~<#�<#�v<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�x<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�u<#�
<#�
<#�<#�<#�<#�<#�<#� <#�%<#�<#�}<#�
<#�
<#�<#�<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929143752201909291437522019092914375620190929143802202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@9��@9��@9��G�O�G�O�G�O�G�O�D�  D�  D�  G�O�G�O�G�� G�� G�� G�� G�� G�� G��                                 6389758         6389758         131072                                          