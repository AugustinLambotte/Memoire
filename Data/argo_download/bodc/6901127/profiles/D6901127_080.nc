CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  H   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       07-Aug-2021 06:45:22Zcreation      
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
resolution        ?�������        ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        L@   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        Y`   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  f�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  i�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  m   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������        pX   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        }x   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        ��   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 H  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 H  �    TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 H  �H   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������        ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���        ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �P   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �P   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �P   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �P   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ۤ   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �8   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �T   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210807064533  20210807064533  6901127 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               PA   BO  123218                          2C  D   APEX                            6231                            120210                          846 @���*�@1   @���*�@@R�1&�x�?�O�;dZ1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @�ff@�  A   A   A@  A`  A�  A���A�ffA�ffA�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� DfD� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� Dz  Dz� D{  D{� D|  D|� D}  D}� D~  D~� D  D� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D��3D�  D�@ D�� D�� D�  D�@ D�� D�� D���D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D�� D�� D�  D�@ D D�� D�  D�@ DÀ D�� D�  D�@ DĀ D�� D�  D�@ Dŀ D�� D�  D�@ Dƀ D�� D�  D�@ Dǀ D�� D�  D�@ DȀ D�� D�  D�@ Dɀ D�� D�  D�@ Dʀ D�� D�  D�@ Dˀ D�� D�  D�@ D̀ D�� D�  D�@ D̀ D�� D�  D�@ D΀ D�� D�  D�@ Dπ D�� D�  D�@ DЀ D�� D�  D�C3Dр D�� D�fD�)�B+B+B+B+B+B+B%B%B
��B
��BVBP�BQ�BL�BL�BM�BO�BP�BR�BVBXB\)BbNBffBhsBiyBl�Bo�Bs�Bv�Bw�Bz�B|�B�B�B�B�%B�+B�=B�\B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@ɑh@ɲ-@�@�@ɩ�@�p�@�hs@�`B@�1'@��
?�V?���?��P?��>O߽�u��ƨ����/���������{�m�N��e`B�`A��R�J���6E��#�
��+�
=q��l���j�q����w�8Q콶E���+�T��=+=��=ȴ9>�>t�=�=�7L=��
=�Q�=�v�=�{=��w=��P=��
=�t�=y�#=T��=�w:�o�8Q�aG������o��o�q����o�������w��{��E��\��vɽ��㽝�-�����9X���-���P�� Ž�j����������㽉7L��+��\)��hs��\)��\)��\)��hs���罶E���vɽ\������G�������#���#�������J�!���/��'+�49X�9X�C���G��L�;O�;�P�`�R�S�ϾQ녾L�;N��Xb�W
=�R�T���Xb�["ѾY��Q녾F��B�\�C���B�\�L�;Y��V�Z��bMӾo���k��j~��p�׾t�j�z�H�}�|푾vȴ�z�H�}󶾀  ��  ���7��J��o������o�~�۾~�۾{�m�l�D�hr��cS��e`B�gl��gl��hr��hr��hr��gl��e`B�j~��k��r�!�t�j�u�t�j�t�j�w�پx���y�#�{�m�}󶾀  �����%���7�����+��1'��1'��C���ƨ��C�����������=q��=q���;�bN���`��hs��녾�n���n���n���n����Ͼ��Ͼ�����
=��b��b���u���u���u��������������������"Ѿ��㾛�㾛�㾛�㾛"Ѿ��㾛�㾛�㾛�㾛�㾜(���(������������(���(���/�����/��/���-���-���-���-���-���-��/���-��5?��;d��A���A���A���A����w���w���w��Ĝ��G���A���A���A���A���A���Ĝ��Ĝ��Ĝ��Ĝ����������������������MӾ�MӾ�MӾ���������Ĝ����������G���G���G�������S���S���S����御G���Ĝ��Ĝ��Ĝ���w���w��G���MӾ�MӾ�MӾ��徢MӾ��御��������A����R��;d���R���R��;d���w��Ĝ��G���G���MӾ�Ĝ��Ĝ��Ĝ��A���A����w��A����w��A���A���A���A���A���A���Ĝ��A���A���A���A����w���w���w��A���Ĝ��Ĝ��Ĝ��A���A���Ĝ��Ĝ��Ĝ��G�����������MӾ��徣S���S����
���
��Z��Z���/���T���T���T��ff��`B���
���
���
���
���/��`B��ff��l����þ��羫����D��V���h��{������������� ž� ž��׾��׾�&龱��������-���!���F��9X���j��?}����E���KǾ��پ�Q쾸����X��X��X��^5���H��dZ��j��푾�󶾾vɾ�|��  ��  ��  �����%���7�\�\�Õ��Õ���������$ݾ�+�Ǯ�Ǯ�ȴ9�ȴ9�ɺ^������ƨ��I���O߾�O߾�������\)��bN��hs��녾�n������Ͼ�z��z�����Ձ�և+��
=��b����ٙ��ٙ��������ڟ��ڟ��ڟ��ۥ�ۥ��(���(��ܬ��/�ݲ-��5?�޸R�޸R��;d�߾w��A���Ĝ��G����������S����
���
��Z��Z���/��`B��`B��`B���T���T��ff���y��l���l����r����þ��þ�xվ����~���~�����������1��D��V��h��h��{��{��{����������� ž�׾�׾�&��&�����-��-��-��!��!��!��!��33��33��F��F��9X��9X��9X���j���j���j��?}��?}������E���ȴ��ȴ��ȴ��KǾ�KǾ�KǾ�KǾ��پ��پ��پ�Q��Q��Q����������X��X���#���#���#���#��^5��^5��^5���H���H��dZ��dZ��dZ��dZ���m���m���m���m��j��j��푾�푾�푾�p���p���p�������������vɾ�vɾ�vɾ��۾��۾��۾��۾�|��|��|��|��|��|�   �   �   � A�� A�� A�� A�� A�� �� �� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ�%� Ĝ� Ĝ�%�%�%�%�%�%�G��G��G���7��7��7��7��7��7��7��7��7��7��7�������������J�J�J�J�MӿMӿMӿMӿMӿMӿ�\��\��\��\�����o�o�o�S��S��S������������������������
��
��
��
��
�������Z�Z�����������/��/�������`B�`B�`B��˿�˿�T��T�$ݿ$ݿ$ݿff�ff�ff�ff����������y��y�+�+�l��l��l������������1'�1'�1'�r��r��r���9��9��ÿ�ÿ	7L�	xտ	xտ	�^�	�^�	�^�	��	��
=q�
~��
���
���
����C��C����ƨ�ƨ�ƨ�1�1�1�I��I���D��Ϳ�ͿV�V�O߿O߿O߿�h��h��h��h�����{�{�V�V������������\)�����������;��;��;� ſ ſ ſbN�bN��׿�׿�`�&�&�hs�hs�������녿녿-�n��n���!��!��!���33�33�t��t���F��F��F��F��Ͽ�Ͽ�Ͽ9X�9X�z��j��j��j�������������������?}111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��
@�p�A�RA"�RAB�RAb�RA�\)A�(�A�A�A�\)A�\)A�\)A�\)B �B�B�B�B �B(�B0�B8�B@�BH�BP�BX�B`�Bh�Bp�Bx�B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
C +�C+�C+�C+�C+�C
+�C+�C+�C+�C+�C+�C+�C+�C+�C+�C+�C +�C"+�C$+�C&+�C(+�C*+�C,+�C.+�C0+�C2+�C4+�C6+�C8+�C:+�C<+�C>+�C@+�CB+�CD+�CF+�CH+�CJ+�CL+�CN+�CP+�CR+�CT+�CV+�CX+�CZ+�C\+�C^+�C`+�Cb+�Cd+�Cf+�Ch+�Cj+�Cl+�Cn+�Cp+�Cr+�Ct+�Cv+�Cx+�Cz+�C|+�C~+�C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C�"�C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D 
�D ��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D	
�D	��D

�D
��D
�D��D
�D��D
�D��D
�D��DGD��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D
�D��D 
�D ��D!
�D!��D"
�D"��D#
�D#��D$
�D$��D%
�D%��D&
�D&��D'
�D'��D(
�D(��D)
�D)��D*
�D*��D+
�D+��D,
�D,��D-
�D-��D.
�D.��D/
�D/��D0
�D0��D1
�D1��D2
�D2��D3
�D3��D4
�D4��D5
�D5��D6
�D6��D7
�D7��D8
�D8��D9
�D9��D:
�D:��D;
�D;��D<
�D<��D=
�D=��D>
�D>��D?
�D?��D@
�D@��DA
�DA��DB
�DB��DC
�DC��DD
�DD��DE
�DE��DF
�DF��DG
�DG��DH
�DH��DI
�DI��DJ
�DJ��DK
�DK��DL
�DL��DM
�DM��DN
�DN��DO
�DO��DP
�DP��DQ
�DQ��DR
�DR��DS
�DS��DT
�DT��DU
�DU��DV
�DV��DW
�DW��DX
�DX��DY
�DY��DZ
�DZ��D[
�D[��D\
�D\��D]
�D]��D^
�D^��D_
�D_��D`
�D`��Da
�Da��Db
�Db��Dc
�Dc��Dd
�Dd��De
�De��Df
�Df��Dg
�Dg��Dh
�Dh��Di
�Di��Dj
�Dj��Dk
�Dk��Dl
�Dl��Dm
�Dm��Dn
�Dn��Do
�Do��Dp
�Dp��Dq
�Dq��Dr
�Dr��Ds
�Ds��Dt
�Dt��Du
�Du��Dv
�Dv��Dw
�Dw��Dx
�Dx��Dy
�Dy��Dz
�Dz��D{
�D{��D|
�D|��D}
�D}��D~
�D~��D
�D��D�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD�ȤD�qD�EqD��qD��qD�qD�EqD��qD��qD�>D�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqD��qD��qD�qD�EqDqD��qD�qD�EqDÅqD��qD�qD�EqDąqD��qD�qD�EqDŅqD��qD�qD�EqDƅqD��qD�qD�EqDǅqD��qD�qD�EqDȅqD��qD�qD�EqDɅqD��qD�qD�EqDʅqD��qD�qD�EqD˅qD��qD�qD�EqD̅qD��qD�qD�EqDͅqD��qD�qD�EqD΅qD��qD�qD�EqDυqD��qD�qD�EqDЅqD��qD�qD�H�DхqD��qD��D�/B�BB-BPByB8BRB
BV|B`bBoBf�Bh�Bb�BY�BS�BQ�BQ�BS�BU�BW�BY�B`8BgOBh8Bh�BlBn�Br�BvBw2By�B{�BlB�#B�wB�&B�B��B�B��B��B�B��B�B��B�&B�6B��B�%B� B��B��B�4B�]B�FB��B��B��B�PB�;B��B��B��B�B�SB��B�'B�B�B��B�B��B��B��B�aB��B�xB�<B�uB��B��B�B��B��B�4B�B��B�B�	B�B��B�QB�4B�B�pB�PB�eB�4B�B��B�B�~B�tB��B��B�B�aB�3B�jB�#B�(B�B��B�B��B��B��B�
B�\B��B��B�B�B�B��B��B�sB��B�B��B�zB��B��B�.B�_B��B��B��B�=B�(B�:B�B��B��B�!B�B�B��B�B�B�B�B��B��B��B��B�JB��B��B�#B�!B�	B�B�B�B��B��B�CB�B�^B�#B�B��B�	B�,B�B�B�B�B�B�B�B�B�AB�IB�B��B�EB�
B��B��B��B��B��B�:B�RB�
B�	B�	B�B��B��B��B�B��B�B�'B�B��B�B��B��B�B�B��B�B�B�B�B��B��B��B��B�B��B��B��B��B�B��B�B��B��B��B��B�B��B�B��B�B��B��B��B��B��B��B�B�B�B�B��B��B��B��B��B��B�B�B��B��B��B��B��B�B��B��B��B�B��B��B��B��B�B��B��B��B��B��B�B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B�B�B��B��B�B��B�B��B��B��B��B�B��B��B�	B�	B�B�B��B�B��B��B��B��B��B��B�B��B�B��B��B��B��B��B�B��B��B��B��B��B��B��B�	B�B��B��B��B��B�B��B��B�	B�B��B�	B�	B�B��B�B��B�B��B�B�B��B��B�B��B��B��B��B��B�B�	B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�	B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B�B��B�B��B�B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@ɑh@ɲ-@�@�@ɩ�@�p�@�hs@�`B@�1'@��
?�V?���?��P?��>O߽�u��ƨ����/���������{�m�N��e`B�`A��R�J���6E��#�
��+�
=q��l���j�q����w�8Q콶E���+�T��=+=��=ȴ9>�>t�=�=�7L=��
=�Q�=�v�=�{=��w=��P=��
=�t�=y�#=T��=�w:�o�8Q�aG������o��o�q����o�������w��{��E��\��vɽ��㽝�-�����9X���-���P�� Ž�j����������㽉7L��+��\)��hs��\)��\)��\)��hs���罶E���vɽ\������G�������#���#�������J�!���/��'+�49X�9X�C���G��L�;O�;�P�`�R�S�ϾQ녾L�;N��Xb�W
=�R�T���Xb�["ѾY��Q녾F��B�\�C���B�\�L�;Y��V�Z��bMӾo���k��j~��p�׾t�j�z�H�}�|푾vȴ�z�H�}󶾀  ��  ���7��J��o������o�~�۾~�۾{�m�l�D�hr��cS��e`B�gl��gl��hr��hr��hr��gl��e`B�j~��k��r�!�t�j�u�t�j�t�j�w�پx���y�#�{�m�}󶾀  �����%���7�����+��1'��1'��C���ƨ��C�����������=q��=q���;�bN���`��hs��녾�n���n���n���n����Ͼ��Ͼ�����
=��b��b���u���u���u��������������������"Ѿ��㾛�㾛�㾛�㾛"Ѿ��㾛�㾛�㾛�㾛�㾜(���(������������(���(���/�����/��/���-���-���-���-���-���-��/���-��5?��;d��A���A���A���A����w���w���w��Ĝ��G���A���A���A���A���A���Ĝ��Ĝ��Ĝ��Ĝ����������������������MӾ�MӾ�MӾ���������Ĝ����������G���G���G�������S���S���S����御G���Ĝ��Ĝ��Ĝ���w���w��G���MӾ�MӾ�MӾ��徢MӾ��御��������A����R��;d���R���R��;d���w��Ĝ��G���G���MӾ�Ĝ��Ĝ��Ĝ��A���A����w��A����w��A���A���A���A���A���A���Ĝ��A���A���A���A����w���w���w��A���Ĝ��Ĝ��Ĝ��A���A���Ĝ��Ĝ��Ĝ��G�����������MӾ��徣S���S����
���
��Z��Z���/���T���T���T��ff��`B���
���
���
���
���/��`B��ff��l����þ��羫����D��V���h��{������������� ž� ž��׾��׾�&龱��������-���!���F��9X���j��?}����E���KǾ��پ�Q쾸����X��X��X��^5���H��dZ��j��푾�󶾾vɾ�|��  ��  ��  �����%���7�\�\�Õ��Õ���������$ݾ�+�Ǯ�Ǯ�ȴ9�ȴ9�ɺ^������ƨ��I���O߾�O߾�������\)��bN��hs��녾�n������Ͼ�z��z�����Ձ�և+��
=��b����ٙ��ٙ��������ڟ��ڟ��ڟ��ۥ�ۥ��(���(��ܬ��/�ݲ-��5?�޸R�޸R��;d�߾w��A���Ĝ��G����������S����
���
��Z��Z���/��`B��`B��`B���T���T��ff���y��l���l����r����þ��þ�xվ����~���~�����������1��D��V��h��h��{��{��{����������� ž�׾�׾�&��&�����-��-��-��!��!��!��!��33��33��F��F��9X��9X��9X���j���j���j��?}��?}������E���ȴ��ȴ��ȴ��KǾ�KǾ�KǾ�KǾ��پ��پ��پ�Q��Q��Q����������X��X���#���#���#���#��^5��^5��^5���H���H��dZ��dZ��dZ��dZ���m���m���m���m��j��j��푾�푾�푾�p���p���p�������������vɾ�vɾ�vɾ��۾��۾��۾��۾�|��|��|��|��|��|�   �   �   � A�� A�� A�� A�� A�� �� �� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ�%� Ĝ� Ĝ�%�%�%�%�%�%�G��G��G���7��7��7��7��7��7��7��7��7��7��7�������������J�J�J�J�MӿMӿMӿMӿMӿMӿ�\��\��\��\�����o�o�o�S��S��S������������������������
��
��
��
��
�������Z�Z�����������/��/�������`B�`B�`B��˿�˿�T��T�$ݿ$ݿ$ݿff�ff�ff�ff����������y��y�+�+�l��l��l������������1'�1'�1'�r��r��r���9��9��ÿ�ÿ	7L�	xտ	xտ	�^�	�^�	�^�	��	��
=q�
~��
���
���
����C��C����ƨ�ƨ�ƨ�1�1�1�I��I���D��Ϳ�ͿV�V�O߿O߿O߿�h��h��h��h�����{�{�V�V������������\)�����������;��;��;� ſ ſ ſbN�bN��׿�׿�`�&�&�hs�hs�������녿녿-�n��n���!��!��!���33�33�t��t���F��F��F��F��Ͽ�Ͽ�Ͽ9X�9X�z��j��j��j�������������������?}111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�<#�
<#ؽ<#�{<#��<#ڹ<#�<0#�=���=�pu<ه�<��o<�^!<��n<�J<@><&M_<$�#<$R�<#��<#��<'S><&�n<$�s<#�<$<#�<$h6<$T$<$N<$
<$�4<$��<%�<$ht<#��<'Gs<$��<({=<'f�<%�<%�
<%a�<$:<%�<<(C<#��<#�<#�#<$<$g<#��<#��<$�<$,a<$�<$S�<&�U<'��<$2<$}<#��<#��<#�;<#��<$(�<#��<#��<#�<#��<#�<$5�<#��<$n3<#�_<#��<#�*<$/�<#�b<$$�<#��<#��<$E<#�<#��<#�y<#ݓ<#ה<#�L<#�m<#�<$.:<#��<#�F<#��<$k<#�<$�<#�<#�J<#�:<#��<$�<%�Z<$M�<#ݞ<#�<$�<#�Z<$.<#��<#�<#�<#�C<#�?<#ۣ<#�<#ڤ<#߃<$�<#ז<#�><#��<#��<#�<#�
<#�<#�0<#�p<#ۼ<#�I<$ <$41<#�<#�g<$ <$<�<#��<#��<#��<#��<#�S<#�<#ב<#��<#�<#� <#�9<#٣<#�4<#�5<#�H<#۴<#׈<#ߢ<#�r<#��<$!�<#�<#�]<#��<#�0<#�<#۫<#�V<#�6<#׬<#� <#��<#ܹ<#��<#��<#��<#��<#ق<#�J<#�<#�\<#��<#�p<#�L<#��<#��<#�R<#�<#�_<#��<#��<#�|<#� <#��<#׵<#�4<#׾<#ٗ<#�g<#��<#�<#��<#��<#ۿ<#�E<#�F<#�}<#�<#ٹ<#ߖ<#�<#�i<#�q<#۽<#�Q<#�V<#ۿ<#۽<#�]<#ۿ<#��<#��<#ۺ<#�O<#�A<#�0<#��<#ۤ<#�Q<#�A<#�A<#�S<#ۦ<#�e<#۟<#�F<#�-<#׼<#�Q<#��<#׼<#ۑ<#�d<#۝<#�E<#�><#�><#�=<#�,<#׻<#۪<#��<#�8<#�<#�<#ٮ<#٨<#׋<#�=<#�a<#��<#ے<#�<#��<#�+<#��<#�L<#ۮ<#�G<#�;<#�_<#��<#�H<#�9<#�9<#�K<#۬<#�D<#�&<#׷<#�<#�<#މ<#�I<#׵<#�7<#�I<#�<#��<#�K<#�"<#י<#�R<#ן<#�4<#�<#�<#�L<#�e<#�<#�V<#�E<#�x<#׵<#�K<#�<#��<#�S<#�e<#�+<#��<#�8<#۾<#��<#�<#�`<#�[<#��<#�e<#��<#�<#ת<#�<#׾<#�V<#׼<#ۈ<#�?<#�/<#�/<#�/<#�@<#�g<#׻<#�&<#�-<#�<#׳<#�!<#�@<#ۡ<#ۙ<#�<<#�<#׭<#�.<#�}<#�9<#�><#۟<#ۘ<#�N<#ۧ<#۰<#ۗ<#�O<#ۄ<#�<<#ۃ<#��<#�W<#��<#�8<#�:<#�C<#�<#�]<#�<#�)<#�M<#��<#��<#�A<#�B<#�Z<#�N<#�V<#�)<#��<#۫<#�f<#�O<#۪<#۩<#ې<#�D<#�u<#�H<#۝<#ێ<#�H<#ۏ<#��<#�<#ۺ<#ۧ<#ۧ<#ۨ<#�{<#߶<#��<#ۦ<#ۥ<#ۊ<#�6<#�H<#��<#۰<#ۿ<#��<#��<#��<#��<#��<#�<#�Z<#�4<#ۋ<#ۣ<#ۼ<#��<#�_<#ޔ<#�y<#��<#�\<#��<#��<#ۙ<#�C<#ޥ<#��<#߰<#�<#��<#��<#ް<#�K<#۹<#��<#��<#��<#��<#ۨ<#۞<#�<#�}<#ۂ<#�<<#ۃ<#۶<#��<#��<#��<#��<#ۓ<#�-<#�n<#�?<#�g<#�#<#�?<#ލ<#�Q<#�l<#�-<#ۆ<#ۘ<#��<#��<#��<#�v<#�1<#ۖ<#ۖ<#ۖ<#ۗ<#۰<#��<#۟<#�{<#�<<#�h<#�*<#ہ<#�y<#�(<#�)<#�o<#�9<#�x<#ۑ<#�x<#�5<#�x<#ې<#�w<#�/<#�y<#��<#��<#�'<#�}<#�t<#�%<#�&<#ۄ<#ێ<#ۍ<#�t<#�%<#�a<#�%<#�%<#ہ<#�s<#�4<#�X<#�1<#ۋ<#�q<#�5<#�Z<#�+<#�<#�o<#�#<#�"<#�n<#�<#�<#�"<#�d<#�)<#�W<#�3<#�[<#�<#�!<#�R<#�<#� <#�b<#�0<#�R<#�-<#ۄ<#�j<#� <#�<#�<#�D<#�<#�<#�V<#�<#�<#�N<#�<#�<#�]<#�-<#�M<#�)<#�f<#�<#�<#�<#�Z<#�<#�<#�K<#�'<#�c<#�<#�
<#�<#�X<#�<#�	<#�<#�P<#�<#�O<#�<#�<#�V<#�<#�<#�F<#�<#�<#�<#�<#�G<#�<#�<#�F<#�<#�<#�<#�J<#�<#�<#�<#�<#�<#�A<#�<#�<#�H<#�<#�<#�<#�<#�B<#�<#�<#�X<#�<#�<#�<#�<#�3<#פ<#�<#�D<#�<#� <#� <#�<#�<#�J<#�<#�<#�B<#��<#��<#��<#��<#��<#��<#��<#��<#��<#�<#�?<#�<#��<#�<#�F<#�<#��<#�<#�O<#�<#��<#��<#��<#�<#�7<#�<#��<#�<#�<<#�<#�5<#�<#��<#�y<#�<<#�`<#�K<#�p<#�s<#��<#��<#��<#�
<#�/<#�<#��<#��<#�	<#�1<#��<#�<#�1<#�<#�5<#��<#�<#�-<#�<#�;<#��<#�<#�+<#�<#�<#�3<#�<#�-<#�<#�9<#�<#�<#�0<#��<#��<#�<#�*<#��<#��<#�<#�(<#�<#�@<#�<#�&<#�<#�<#�4<#؞<#�j<#��<#�<#�,<#�<#�&<#��<#�<#�%<#�<#�<#�<<#�<#�#<#�<#�I<#�;<#�<#�!<#��<#��<#�(<#�<#�<<#�Q<#�8<#��<#��<#�?<#�8<#ؼ<#��<#�6<#��<#��<#�<#��<#��<#�6<#�	<#�5<#�5<#�<#�+<#�<#�<#��<#��<#�!<#��<#��<#��<#�)<#�
<#�<#�<#�2<#�<#�0<#�1<#� <#�'<#�	<#ڎ<#ۨ<#��<#��<#�<#��<#��<#�$<#��<#��<#�<#��<#�<#�	<#�.<#�-<#��<#�<#�<#�<#��<#�<#�<#�,<#�+<#��<#�<#��<#��<#�)<#� <#�<#�<#�<#��<#�<#؈<#�5<#��<#�<#��<#��<#�<#��<#�'<#�%<#��<#��<#�<#��<#��<#��<#��<#��<#ܸ<#�b;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.17                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 2/0.5.MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2021v01 and ARGO2020v03 ref. data.                                                                                                                    202108041426132021080614305220210804142613202108041426132021080614305220210806143052BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190824141056201908241410562019082414110120190824141107202006251447052021080414261320210806143052  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@�ff@�ff@�ffG�O�G�O�G�O�G�O�D�)�D�)�D�)�G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          