CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       06-Aug-2021 03:59:42Zcreation      
references        (http://www.argodatamgt.org/Documentation   comment       bThis netCDF file is generated using BODC's argoReader and netCDF writer software (argo@bodc.ac.uk)     user_manual_version       3.4    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^X   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  m�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  q�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  u�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  y�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �H   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �h   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  �P   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  �8   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �p   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �l   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �l   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �l   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �l   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �0   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �T   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �p   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �4Argo profile    3.1 1.2 19500101000000  20210806035957  20210806035957  6901128 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               ~A   BO  124373                          2C  D   APEX                            6232                            120210                          846 @ׁ�d� 1   @ׁ�d� @R���v�?ش9XbN1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @�ff@���A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B ffB  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR�CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0fD0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq�fDr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D���D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�C3D�� D�� D�  D�@ D�� D�� D�3D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�|�D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D D�� D�  D�@ DÃ3D�� D�  D�@ DĀ D�� D�  D�@ Dŀ D�� D�  D�@ Dƀ D�� D�  D�@ Dǀ D�� D�  D�@ DȀ D�� D�  D�@ Dɀ D�� D�  D�@ Dʀ D�� D�  D�@ Dˀ D�� D�  D�@ D̀ D�� D�  D�@ D̀ D�� D�  D�@ D΀ D�� D�  D�@ Dπ D�� D�  D�@ DЀ D�� D�  D�@ Dр D�� D�  D�@ DҀ D�� D�  D�@ DӀ D�� D�  D�@ DԀ D�� D�  D�@ DՀ D�� D�  D�@ Dր D�� D�  D�@ D׀ D�� D�  D�@ D؀ D�� D�  D�@ Dـ D�� D�  D�@ Dڀ D�� D�  D�@ Dۀ D�� D�  D�@ D܀ D�� D�  D�@ D݀ D�� D�  D�@ Dހ D�� D�  D�@ D߀ D�� D�  D�@ D�� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D�� D�� D�  D�@ D� D�� D�  D�@ D� D��D�  D�@ D� D�� D�  D�@ D� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D��3D���D��BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BI�BI�BI�BJ�BH�BH�BH�BH�BH�BH�BH�BI�BJ�BL�BM�BN�BO�BR�BS�BS�BXBZBZBZBZBZBZBZB[#B]/B`BBn�Bm�Bw�B� B�B�B�%B�%B�+B�1B�7B�DB�JB�PB�PB�VB�\B�\B�oB�{B�{B�uB�oB�hB�bB�bB�bB�hB�oB�oB�hB�bB�bB�hB�bB�bB�bB�hB�hB�hB�hB�oB�uB�uB�oB�oB�oB�{B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?F$�?F$�?F$�?F$�?F$�?F$�?F$�?Fff?Fff?Fff?Fff?F��?F�y?F��?G+?H1'?Fff?F��?F��?F��?F��?F��?F�y?F�y?Hr�?I��?K?LI�?L��?O\)?P�`?P�`?T�j?W
=?WK�?WK�?W�P?W�P?W��?Xb?W�P?P�`?@�?��?G�>��R>�=���=�w<�1;ě�<�/�o;o<�1<e`B<�9X<�o<�C�<�`B<�h=L��=�1=���=�+=ix�=+�o��`B��o<t�<e`B=t�<o��1�'<j�Y���%�������㽅�������T���-��hs�������\��xս� ŽD���#�
��P��P����w��O߽����/��h���o�����o�%�����m��F����������l���l���l������ȴ9�Ƨ�ě��\�\��������vɽ�vɽ�����ͽ��ͽ��������ȴ9�\�����vɽ�j��E��� Ž��罧���罬1��������w���㽗�P���P��C����-��7L�����O߽����hs��hs��\)��C���O߽�t����P��O߽�+��%�y�#�y�#�q���ixսaG��m�h��%��+��o�aG��e`B��%���P���㽼j�ȴ9�������������"ѽ�"ѽ����;d��"ѽ��`���������"ѽ��������   �o�$ݾ+�����   ���m���ٽ��#���m���%�$ݾ1'�
=q�	7L�1'�1'�
=q�I��V�\)�bN�bN�z����P������R��R��R��w� Ĝ�"��%�T�''(�þ+�.{�1&�333�49X�5?}�6E��8Q�8Q�8Q�9X�9X�:^5�<j�=p��>vɾ>vɾ?|�?|�A�7�B�\�C���E�˾G��I�^�Kƨ�L�;N��P�`�Q녾Q녾S�ϾV�Xb�Z��\(��]/�]/�^5?�`A��aG��dZ�e`B�hr��ixվk��k��k��l�D�l�D�l�D�l�D�l�D�k��j~��gl��gl��ixվl�D�k��l�D�o���p�׾o���m�h�j~��k��k��l�D�l�D�n���q���p�׾q���r�!�s�F�t�j�t�j�s�F�t�j�t�j�vȴ�w�پx���x���y�#�z�H�z�H�{�m�y�#�vȴ�t�j�s�F�u�vȴ�vȴ�x���vȴ�|푾y�#�y�#�{�m�}󶾀  ��J��o���\�~�۾~�۾{�m�~�۾���������������������������������˾�$ݾ�$ݾ���+�����������������1'��1'���9��7L���^���^���^���^��=q������C���C���C���ƨ��I����;��;�O߾�O߾���V��V��V������;���;���;��bN��bN���`��hs��hs��녾�n���t���t���t���zᾔ����������������+��
=���P��b��b��b��b���������������������㾜(���(���������/��/��/���-��5?��;d��;d��A���Ĝ���������������徣S���Z���/���/��`B��`B��`B���T���T��ff��l�����r���r����þ�xվ��義�羪~������1���D��V��V���h��{��{��{���������������� ž��׾�&龲-���!���!���!��33��33��33���F���F��E���ȴ��KǾ�KǾ�Q쾹X��^5���m��j��푾�푾�󶾾vɾ��۾�|��  ��%���7��J��J��o��������$ݾƧ��1'��1'�ȴ9�ȴ9��7L�ɺ^������ƨ��O߾�����;���`��n���t����Ͼ�z�����Ձ�և+��
=��b�ؓu����ٙ��ٙ��ٙ�����ڟ��ۥ�ܬ��/�ݲ-��5?�޸R�޸R��;d�߾w�߾w��A���A���A���A���A���Ĝ��G���G�������MӾ����S���S���Z��`B��`B��`B��`B���T��ff���y��l����r���r����þ�xվ����~���~�������1��1��D��V��V��h��h��{������������ ž� ž� ž� ž�׾�&��-��-��!��33��F��9X��9X���j��?}��?}��?}����E���E���E���ȴ��ȴ��ȴ��ȴ��KǾ�KǾ�KǾ��پ��پ��پ�Q��Q��Q��Q��Q������X���#��^5���H��dZ���m��j��j��j��푾�푾�푾�푾�푾�푾�푾�푾�p�����vɾ�vɾ��۾�|��|���۾�|�   �   � A�� A�� A�� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ�%�G���7�������J�J�MӿMӿMӿ�\��\��\���o�o�����������������
�����Z��������/�����`B��˿�T��T��T�$ݿ$ݿ$ݿ$ݿ$ݿff�ff���+�l����1'�1'�1'�r���ÿ�ÿ	7L�	xտ	xտ	�^�	�^�	��	��
=q�
~��
~��
���
�������C��C������ƨ�ƨ�1�1�I���D��Ϳ�Ϳ�ͿV�V�V�V�V�O߿O߿O߿�h��h��h���{�{�{�V�����������\)�����������;��;� ſ ſ ſ ſbN��׿�`�&�&�&�hs����녿녿녿녿녿-�n��n��n��n��n���!����33�t���F��Ͽ�Ͽ�Ͽ9X�z�z�z�z��j��j��j����?}��������������E���+�ȴ�ȴ�
=�
=�
=�
=�KǿKǿ�P��P��P��ٿ�P��P��P��P��ٿ�ٿb�b�b�Q�Q�Q�Q��u��u��u��u��u�����������������X�X�����������#��#���^5�^5�������������"ѿdZ������m��m��m�(��j�j�푿/�p��p��p���-������5?�vɿvɿ�R��R��ۿ�ۿ;d�;d�;d�|�|��w��w��w��w��w�   �   � A�� A�� A�� �� Ĝ� Ĝ�!%�!G��!G��!G��!�7�!�7�"J�"Mӿ"�\�"�\�"�\�"��#o�#o�#o�#o�#S��#S��#S��#S��#���#�
�$��$��$��#�
�$��$��#�
�$��$��$Z�$���$���$���$���$���$���$�/�%��%��%`B�%`B�%`B�%`B�%�˿%�˿%�T�%�T�%�T�&$ݿ&ff�&ff�&ff�&��&��&��&��&�y�'+�'l��'(1'�(r��(�9�)7L�)xտ)xտ)7L�)xտ)J�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@ָRA�\A.�\AN�\An�\A�G�A�G�A�G�A�G�A�G�A�G�A�G�A�G�B
=B��B��B��B#��B+��B3��B;��BC��BK��BS��B[��Bc��Bk��Bs��B{��B���B���B���B���B���B���B���B���B���B���B���B���B���B�B�B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���C ��C��C��C��C��C
��C��C��C��C��C��C��C��C��C��C��C ��C"��C$��C&��C(��C*��C,��C.��C0��C2��C4��C6��C8��C:��C<��C>��C@��CB��CD��CF��CH��CJ��CL��CN��CP��CS�CT��CV��CX��CZ��C\��C^��C`��Cb��Cd��Cf��Ch��Cj��Cl��Cn��Cp��Cr��Ct��Cv��Cx��Cz��C|��C~��C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�g�C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{C�t{D :=D �=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D	:=D	�=D
:=D
�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D:=D�=D :=D �=D!:=D!�=D":=D"�=D#:=D#�=D$:=D$�=D%:=D%�=D&:=D&�=D':=D'�=D(:=D(�=D):=D)�=D*:=D*�=D+:=D+�=D,:=D,�=D-:=D-�=D.:=D.�=D/:=D/�=D0@�D0�=D1:=D1�=D2:=D2�=D3:=D3�=D4:=D4�=D5:=D5�=D6:=D6�=D7:=D7�=D8:=D8�=D9:=D9�=D::=D:�=D;:=D;�=D<:=D<�=D=:=D=�=D>:=D>�=D?:=D?�=D@:=D@�=DA:=DA�=DB:=DB�=DC:=DC�=DD:=DD�=DE:=DE�=DF:=DF�=DG:=DG�=DH:=DH�=DI:=DI�=DJ:=DJ�=DK:=DK�=DL:=DL�=DM:=DM�=DN:=DN�=DO:=DO�=DP:=DP�=DQ:=DQ�=DR:=DR�=DS:=DS�=DT:=DT�=DU:=DU�=DV:=DV�=DW:=DW�=DX:=DX�=DY:=DY�=DZ:=DZ�=D[:=D[�=D\:=D\�=D]:=D]�=D^:=D^�=D_:=D_�=D`:=D`�=Da:=Da�=Db:=Db�=Dc:=Dc�=Dd:=Dd�=De:=De�=Df:=Df�=Dg:=Dg�=Dh:=Dh�=Di:=Di�=Dj:=Dj�=Dk:=Dk�=Dl:=Dl�=Dm:=Dm�=Dn:=Dn�=Do:=Do�=Dp:=Dp�=Dq:=Dq��Dr:=Dr�=Ds:=Ds�=Dt:=Dt�=Du:=Du�=Dv:=Dv�=Dw:=Dw�=Dx:=Dx�=Dy:=Dy�=Dz:=Dz�=D{:=D{�=D|:=D|�=D}:=D}�=D~:=D~�=D:=D�=D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D��D�]D��D��D�D�]D��D��D�D�]D��D��D�D�`RD��D��D�D�]D��D��D� RD�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D���D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]DD��D�D�]DàRD��D�D�]DĝD��D�D�]DŝD��D�D�]DƝD��D�D�]DǝD��D�D�]DȝD��D�D�]DɝD��D�D�]DʝD��D�D�]D˝D��D�D�]D̝D��D�D�]D͝D��D�D�]DΝD��D�D�]DϝD��D�D�]DНD��D�D�]DѝD��D�D�]DҝD��D�D�]DӝD��D�D�]DԝD��D�D�]D՝D��D�D�]D֝D��D�D�]DםD��D�D�]D؝D��D�D�]DٝD��D�D�]DڝD��D�D�]D۝D��D�D�]DܝD��D�D�]DݝD��D�D�]DޝD��D�D�]DߝD��D�D�]D��D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D��D�D�]D�D���D�D�]D�D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��D��D�D�]D��RD���D��BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BH�BI�BI�BI�BK	BH�BH�BH�BH�BH�BH�BH�BIrBJzBL�BM�BN�BOnBR�BS�BSMBW�BZ	BZBZBZBZBZBZBB\vB`tBh�Bq�Bw~B�vB��B�B��B��B�DB�xB��B�VB��B��B��B�?B��B�?B�]B��B��B�FB��B��B�B��B�;B��B�!B�qB��B��B�^B��B��B��B��B�mB��B��B��B�:B�#B�/B��B�_B�KB�B��B�B�gB��B��B��B��B��B�B��B��B��B��B�rB�qB�kB�mB�yB�YB�qB��B��B��B�zB�{B��B��B�0B�UB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�8B�ZB��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�'B�B��B��B��B�EB�qB�	B��B�1B��B��B�	B�B��B��B�B��B��B��B�B�B�MB�TB��B��B��B�B�B��B��B��B��B��B��B��B��B��B�B�-B�B�B��B��B��B�	B�	B�	B��B��B��B�!B��B�
B�!B��B�	B��B��B��B��B�
B�B�	B��B��B�
B�B�B�	B��B��B��B�
B��B��B��B��B��B�	B��B��B��B��B��B�B��B��B�	B�	B�	B�	B��B�	B�	B��B��B�	B�	B�	B�	B�	B��B��B��B�	B��B�B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B�
B�B��B��B�B��B��B��B��B�B��B�B��B�B�B��B�B�B�B�B��B��B�B��B�B�B�B��B�B�B��B�B��B��B��B��B�B�	B��B�B��B�DB��B��B�B�B�B�-B�B��B��B��B��B�%B�ZB�B��B��B��B�B��B��B�	B�B�	B��B�B�	B�B��B��B��B��B�B��B�B�	B�B��B��B��B�	B�	B�B��B��B�	B�	B�B��B�B��B�B�B��B��B�	B�B��B��B�B��B�B�B��B�	B�	B�B��B��B�B�	B��B�	B�B��B�B�	B�	B�B��B��B��B�B�	B�	B�B��B�B�	B��B�B��B�B��B��B�	B�	B�B��B�B�B�B��B��B�B�B�B�B��B�B��B��B�B��B�B�B�B��B�B��B�B�B�B��B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B��B��B��B��B��B��B��B��B�,B��B��B��B�B�B�B�B��B��B��B�B��B��B��B��B�B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B�B�B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?F$�?F$�?F$�?F$�?F$�?F$�?F$�?Fff?Fff?Fff?Fff?F��?F�y?F��?G+?H1'?Fff?F��?F��?F��?F��?F��?F�y?F�y?Hr�?I��?K?LI�?L��?O\)?P�`?P�`?T�j?W
=?WK�?WK�?W�P?W�P?W��?Xb?W�P?P�`?@�?��?G�>��R>�=���=�w<�1;ě�<�/�o;o<�1<e`B<�9X<�o<�C�<�`B<�h=L��=�1=���=�+=ix�=+�o��`B��o<t�<e`B=t�<o��1�'<j�Y���%�������㽅�������T���-��hs�������\��xս� ŽD���#�
��P��P����w��O߽����/��h���o�����o�%�����m��F����������l���l���l������ȴ9�Ƨ�ě��\�\��������vɽ�vɽ�����ͽ��ͽ��������ȴ9�\�����vɽ�j��E��� Ž��罧���罬1��������w���㽗�P���P��C����-��7L�����O߽����hs��hs��\)��C���O߽�t����P��O߽�+��%�y�#�y�#�q���ixսaG��m�h��%��+��o�aG��e`B��%���P���㽼j�ȴ9�������������"ѽ�"ѽ����;d��"ѽ��`���������"ѽ��������   �o�$ݾ+�����   ���m���ٽ��#���m���%�$ݾ1'�
=q�	7L�1'�1'�
=q�I��V�\)�bN�bN�z����P������R��R��R��w� Ĝ�"��%�T�''(�þ+�.{�1&�333�49X�5?}�6E��8Q�8Q�8Q�9X�9X�:^5�<j�=p��>vɾ>vɾ?|�?|�A�7�B�\�C���E�˾G��I�^�Kƨ�L�;N��P�`�Q녾Q녾S�ϾV�Xb�Z��\(��]/�]/�^5?�`A��aG��dZ�e`B�hr��ixվk��k��k��l�D�l�D�l�D�l�D�l�D�k��j~��gl��gl��ixվl�D�k��l�D�o���p�׾o���m�h�j~��k��k��l�D�l�D�n���q���p�׾q���r�!�s�F�t�j�t�j�s�F�t�j�t�j�vȴ�w�پx���x���y�#�z�H�z�H�{�m�y�#�vȴ�t�j�s�F�u�vȴ�vȴ�x���vȴ�|푾y�#�y�#�{�m�}󶾀  ��J��o���\�~�۾~�۾{�m�~�۾���������������������������������˾�$ݾ�$ݾ���+�����������������1'��1'���9��7L���^���^���^���^��=q������C���C���C���ƨ��I����;��;�O߾�O߾���V��V��V������;���;���;��bN��bN���`��hs��hs��녾�n���t���t���t���zᾔ����������������+��
=���P��b��b��b��b���������������������㾜(���(���������/��/��/���-��5?��;d��;d��A���Ĝ���������������徣S���Z���/���/��`B��`B��`B���T���T��ff��l�����r���r����þ�xվ��義�羪~������1���D��V��V���h��{��{��{���������������� ž��׾�&龲-���!���!���!��33��33��33���F���F��E���ȴ��KǾ�KǾ�Q쾹X��^5���m��j��푾�푾�󶾾vɾ��۾�|��  ��%���7��J��J��o��������$ݾƧ��1'��1'�ȴ9�ȴ9��7L�ɺ^������ƨ��O߾�����;���`��n���t����Ͼ�z�����Ձ�և+��
=��b�ؓu����ٙ��ٙ��ٙ�����ڟ��ۥ�ܬ��/�ݲ-��5?�޸R�޸R��;d�߾w�߾w��A���A���A���A���A���Ĝ��G���G�������MӾ����S���S���Z��`B��`B��`B��`B���T��ff���y��l����r���r����þ�xվ����~���~�������1��1��D��V��V��h��h��{������������ ž� ž� ž� ž�׾�&��-��-��!��33��F��9X��9X���j��?}��?}��?}����E���E���E���ȴ��ȴ��ȴ��ȴ��KǾ�KǾ�KǾ��پ��پ��پ�Q��Q��Q��Q��Q������X���#��^5���H��dZ���m��j��j��j��푾�푾�푾�푾�푾�푾�푾�푾�p�����vɾ�vɾ��۾�|��|���۾�|�   �   � A�� A�� A�� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ�%�G���7�������J�J�MӿMӿMӿ�\��\��\���o�o�����������������
�����Z��������/�����`B��˿�T��T��T�$ݿ$ݿ$ݿ$ݿ$ݿff�ff���+�l����1'�1'�1'�r���ÿ�ÿ	7L�	xտ	xտ	�^�	�^�	��	��
=q�
~��
~��
���
�������C��C������ƨ�ƨ�1�1�I���D��Ϳ�Ϳ�ͿV�V�V�V�V�O߿O߿O߿�h��h��h���{�{�{�V�����������\)�����������;��;� ſ ſ ſ ſbN��׿�`�&�&�&�hs����녿녿녿녿녿-�n��n��n��n��n���!����33�t���F��Ͽ�Ͽ�Ͽ9X�z�z�z�z��j��j��j����?}��������������E���+�ȴ�ȴ�
=�
=�
=�
=�KǿKǿ�P��P��P��ٿ�P��P��P��P��ٿ�ٿb�b�b�Q�Q�Q�Q��u��u��u��u��u�����������������X�X�����������#��#���^5�^5�������������"ѿdZ������m��m��m�(��j�j�푿/�p��p��p���-������5?�vɿvɿ�R��R��ۿ�ۿ;d�;d�;d�|�|��w��w��w��w��w�   �   � A�� A�� A�� �� Ĝ� Ĝ�!%�!G��!G��!G��!�7�!�7�"J�"Mӿ"�\�"�\�"�\�"��#o�#o�#o�#o�#S��#S��#S��#S��#���#�
�$��$��$��#�
�$��$��#�
�$��$��$Z�$���$���$���$���$���$���$�/�%��%��%`B�%`B�%`B�%`B�%�˿%�˿%�T�%�T�%�T�&$ݿ&ff�&ff�&ff�&��&��&��&��&�y�'+�'l��'(1'�(r��(�9�)7L�)xտ)xտ)7L�)xտ)J�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
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
<#�
<#�
<#�
<#�
<#�
<(^�<#�
<:��<(��<#�
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
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.91                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 2/0.5.MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2021v01 and ARGO2020v03 ref. data.                                                                                                                    202108041430502021080511580420210804143050202108041430502021080511580420210805115804BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190828172144201908281721442019082817214920190828172156202006251118432021080414305020210805115804  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�D��D��D��G�O�G�O�G��G��G��G��G��G��G��                                6389758         6389758         131072                                          