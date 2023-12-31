CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  1   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:11Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    6�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7T   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    7�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    7�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     7�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     88   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           8\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            8h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  B`   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  MX   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  XP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  cH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  n@   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �$   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �d   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �t   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �x   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090443  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�5Ѝ��z1   @�5�Xf��@OcEP�-��@|գ@�32   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AQ��A`  Ap  A~ffA�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A���A���A���A�  B   B  B  B  B  B  B  B��B��B$  B(ffB,  B0  B4ffB8  B<  B@  BD  BH  BL  BPffBT  BX  B\ffB_��Bc��Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  C  C� C  C	� C  C��C  C� C  C� C  C� C   C"� C%  C'� C*  C,��C/�C1� C4  C6� C9  C;ffC>  C@��CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY��C\  C^ffCa  Cc� Cf�Ch� Ck  Cm� Cp  Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�33C�s3C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�L�Cʀ C�� C��C�L�Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C�  C�@ C� C�� C�  C�@ C��C�� C�  C�@ C� C�� C�  C�@ C� C���C��C�@ C�� C�� C��C�� C�  D � D��D@ D�fD� D��D@ D	� D
�fD  D@ D� D� D  D@ D�fD� D  D@ D� D� D  D@ D� D� D   D!FfD"� D#��D%  D&@ D'�fD(� D*  D+9�D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;y�D<� D>fD?@ D@� DA��DB��DD9�DEy�DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db9�Dcy�Dd��Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}��D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D��3D�C3D�� D��3D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�` D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�3D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�  D���D�\�D�  Dà D�C3D�� Dŀ D�  D�� D�\�D���DȠ D�C3D�� Dʀ D�  D�� D�` D�  D͠ D�@ D���Dπ D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  D�� D�c3D�  Dנ D�@ D���Dـ D�#3D��3D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D���D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�3D��3D�@ D��3D��3D�#3D��3D�` D�  D��3D�C3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AQ��A`  Ap  A~ffA�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A���A���A���A�  B   B  B  B  B  B  B  B��B��B$  B(ffB,  B0  B4ffB8  B<  B@  BD  BH  BL  BPffBT  BX  B\ffB_��Bc��Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  C  C� C  C	� C  C��C  C� C  C� C  C� C   C"� C%  C'� C*  C,��C/�C1� C4  C6� C9  C;ffC>  C@��CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY��C\  C^ffCa  Cc� Cf�Ch� Ck  Cm� Cp  Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�33C�s3C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�L�Cʀ C�� C��C�L�Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C�  C�@ C� C�� C�  C�@ C��C�� C�  C�@ C� C�� C�  C�@ C� C���C��C�@ C�� C�� C��C�� C�  D � D��D@ D�fD� D��D@ D	� D
�fD  D@ D� D� D  D@ D�fD� D  D@ D� D� D  D@ D� D� D   D!FfD"� D#��D%  D&@ D'�fD(� D*  D+9�D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;y�D<� D>fD?@ D@� DA��DB��DD9�DEy�DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db9�Dcy�Dd��Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}��D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D��3D�C3D�� D��3D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�` D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�3D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�  D���D�\�D�  Dà D�C3D�� Dŀ D�  D�� D�\�D���DȠ D�C3D�� Dʀ D�  D�� D�` D�  D͠ D�@ D���Dπ D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  D�� D�c3D�  Dנ D�@ D���Dـ D�#3D��3D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D���D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�3D��3D�@ D��3D��3D�#3D��3D�` D�  D��3D�C3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��9@��@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��/@��/@��`@��`@��`@��`@��`@��`@��`@��`@��@���@���@���@��j@���@��j@��9@���@��/@��`@��@��`@��`@���@���@���@���@���@��@���@���@���@���@���@���@���@���@��`@���@���@���@���@���@���@��@��@��@���@���@��u@��u@��D@��D@��9@��@�z�@�r�@�z�@�z�@�z�@��@��D@��D@��@�z�@�z�@��D@��D@��u@��D@�z�@�z�@�r�@�r�@�Z@��D@���@��@��D@�j@�Q�@�I�@�1'@�1'@�1'@�bN@�r�@�r�@�z�@�ƨ@��P@���@��P@��@��@�dZ@�C�@�K�@�S�@�dZ@���@�l�@�|�@�t�@�o@��@�ȴ@���@��!@��\@��\@���@���@���@���@���@��\@��\@��+@�v�@�v�@�v�@�n�@�v�@�n�@�-@��@��T@���@�p�@�O�@�`B@�7L@��`@�Z@��;@�+@�5?@���@���@�X@�7L@�X@���@��7@�7L@�bN@�|�@�l�@�"�@���@�v�@�=q@���@���@���@���@�Ĝ@��@�9X@� �@��;@��F@���@�dZ@��@���@��+@�n�@�E�@��#@��@�X@�G�@�G�@��@��/@��D@�z�@�Q�@�9X@�  @��@��+@�v�@�V@�=q@�J@�@���@��-@���@�hs@���@�Z@�9X@�  @��;@���@�ƨ@���@��@�\)@�S�@�;d@�"�@�
=@��y@��!@�ff@�M�@�E�@�-@��@�J@��@��@���@���@��7@�x�@�?}@�&�@��@�V@�V@�%@�%@��@���@��9@���@���@�1'@��;@��@��@�S�@���@���@��\@�ff@�=q@��@���@��#@�@���@��h@�O�@���@���@���@��j@��u@�I�@�A�@�1'@�w@��@l�@~ȴ@}�-@|�/@|Z@|�@{�
@{��@{o@z��@z=q@y�@y�#@y�@y�#@y��@y��@y�7@yx�@y�^@y�^@x�`@x �@w�w@w|�@w|�@w|�@w;d@v�y@v��@vv�@vV@v{@u�@up�@up�@u`B@u`B@up�@up�@u�@t�@sC�@r��@q��@q�#@qhs@p�9@pA�@p  @o;d@o
=@nȴ@n��@n��@oK�@p��@p��@p�@pr�@pA�@pQ�@pQ�@pQ�@pA�@p �@o�@o�@o��@p1'@qhs@q��@q��@r=q@rn�@qx�@q�@p��@pQ�@pQ�@p��@p��@p��@p�@p�u@p�u@p�u@p��@pĜ@p��@qG�@qX@qG�@qx�@q��@q��@rJ@r�@r-@r-@r-@r-@r-@r-@r�@q��@q�@q��@q��@q��@q��@qx�@qx�@qhs@q7L@q%@q%@p��@p��@p��@p�9@p�@pr�@pr�@pA�@p  @o�w@o��@o�P@ol�@oK�@o;d@o�@n�@nȴ@n�R@n��@nV@n@m�@m�@m��@m�-@m�-@m�h@m?}@mV@l�/@l��@lI�@k��@kƨ@k�F@kt�@ko@k33@ko@j�@j�!@j�!@j�!@j�!@j��@j�\@j�\@j�\@jn�@jM�@j=q@j-@jJ@i�#@i��@i�7@iX@h��@hbN@h1'@h1'@hb@g�@g�@g\)@g
=@f�R@fff@f$�@e�@e@e@e@e��@e�@e`B@eO�@e/@d�@d�j@dI�@d(�@c��@c�F@c�F@c��@c33@b��@bM�@a��@a�#@a�^@a��@a��@aX@a�@`�`@`��@a%@a&�@`�9@`��@`�`@`�`@`Ĝ@`�@`1'@_�@_�w@_+@^��@^5?@]p�@\��@\I�@[�
@Y��@Y��@Y�7@YG�@X��@X �@W��@W;d@V��@V�+@Vff@V5?@U�@U�h@T��@T�j@T��@Tz�@Tj@T�@T1@S��@S�m@S�F@SS�@S"�@S@R��@Rn�@Q��@Q��@QG�@Q%@PĜ@PA�@P  @O�;@O�P@OK�@O;d@O+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��9@��@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��/@��/@��`@��`@��`@��`@��`@��`@��`@��`@��@���@���@���@��j@���@��j@��9@���@��/@��`@��@��`@��`@���@���@���@���@���@��@���@���@���@���@���@���@���@���@��`@���@���@���@���@���@���@��@��@��@���@���@��u@��u@��D@��D@��9@��@�z�@�r�@�z�@�z�@�z�@��@��D@��D@��@�z�@�z�@��D@��D@��u@��D@�z�@�z�@�r�@�r�@�Z@��D@���@��@��D@�j@�Q�@�I�@�1'@�1'@�1'@�bN@�r�@�r�@�z�@�ƨ@��P@���@��P@��@��@�dZ@�C�@�K�@�S�@�dZ@���@�l�@�|�@�t�@�o@��@�ȴ@���@��!@��\@��\@���@���@���@���@���@��\@��\@��+@�v�@�v�@�v�@�n�@�v�@�n�@�-@��@��T@���@�p�@�O�@�`B@�7L@��`@�Z@��;@�+@�5?@���@���@�X@�7L@�X@���@��7@�7L@�bN@�|�@�l�@�"�@���@�v�@�=q@���@���@���@���@�Ĝ@��@�9X@� �@��;@��F@���@�dZ@��@���@��+@�n�@�E�@��#@��@�X@�G�@�G�@��@��/@��D@�z�@�Q�@�9X@�  @��@��+@�v�@�V@�=q@�J@�@���@��-@���@�hs@���@�Z@�9X@�  @��;@���@�ƨ@���@��@�\)@�S�@�;d@�"�@�
=@��y@��!@�ff@�M�@�E�@�-@��@�J@��@��@���@���@��7@�x�@�?}@�&�@��@�V@�V@�%@�%@��@���@��9@���@���@�1'@��;@��@��@�S�@���@���@��\@�ff@�=q@��@���@��#@�@���@��h@�O�@���@���@���@��j@��u@�I�@�A�@�1'@�w@��@l�@~ȴ@}�-@|�/@|Z@|�@{�
@{��@{o@z��@z=q@y�@y�#@y�@y�#@y��@y��@y�7@yx�@y�^@y�^@x�`@x �@w�w@w|�@w|�@w|�@w;d@v�y@v��@vv�@vV@v{@u�@up�@up�@u`B@u`B@up�@up�@u�@t�@sC�@r��@q��@q�#@qhs@p�9@pA�@p  @o;d@o
=@nȴ@n��@n��@oK�@p��@p��@p�@pr�@pA�@pQ�@pQ�@pQ�@pA�@p �@o�@o�@o��@p1'@qhs@q��@q��@r=q@rn�@qx�@q�@p��@pQ�@pQ�@p��@p��@p��@p�@p�u@p�u@p�u@p��@pĜ@p��@qG�@qX@qG�@qx�@q��@q��@rJ@r�@r-@r-@r-@r-@r-@r-@r�@q��@q�@q��@q��@q��@q��@qx�@qx�@qhs@q7L@q%@q%@p��@p��@p��@p�9@p�@pr�@pr�@pA�@p  @o�w@o��@o�P@ol�@oK�@o;d@o�@n�@nȴ@n�R@n��@nV@n@m�@m�@m��@m�-@m�-@m�h@m?}@mV@l�/@l��@lI�@k��@kƨ@k�F@kt�@ko@k33@ko@j�@j�!@j�!@j�!@j�!@j��@j�\@j�\@j�\@jn�@jM�@j=q@j-@jJ@i�#@i��@i�7@iX@h��@hbN@h1'@h1'@hb@g�@g�@g\)@g
=@f�R@fff@f$�@e�@e@e@e@e��@e�@e`B@eO�@e/@d�@d�j@dI�@d(�@c��@c�F@c�F@c��@c33@b��@bM�@a��@a�#@a�^@a��@a��@aX@a�@`�`@`��@a%@a&�@`�9@`��@`�`@`�`@`Ĝ@`�@`1'@_�@_�w@_+@^��@^5?@]p�@\��@\I�@[�
@Y��@Y��@Y�7@YG�@X��@X �@W��@W;d@V��@V�+@Vff@V5?@U�@U�h@T��@T�j@T��@Tz�@Tj@T�@T1@S��@S�m@S�F@SS�@S"�@S@R��@Rn�@Q��@Q��@QG�@Q%@PĜ@PA�@P  @O�;@O�P@OK�@O;d@O+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�sB�sB�sB�sB�sB�sB�yB�yB�yB�yB�yB�yB�yB�sB�sB�sB�sB�sB�mB�fB�fB�`B�`B�`B�ZB�NB�BB�5B�#B�B��B��B��B��B��B��B��BɺBǮBƨBŢBŢBĜBĜBĜBÖBÖBBBBBB��BBBBB��B��B��B�}B�}B�wB�jB�jB�jB�jB�dB�dB�^B�^B�^B�XB�XB�LB�FB�FB�FB�?B�?B�9B�9B�3B�3B�-B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�uB�oB�\B�VB�PB�PB�JB�DB�=B�7B�1B�1B�1B�+B�1B�=B�VB�VB�VB�\B�\B�\B�bB�bB�bB�bB�bB�bB�bB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�!B�!B�'B�-B�3B�3B�9B�9B�9B�9B�?B�?B�?B�?B�?B�?B�?B�FB�FB�FB�FB�LB�LB�LB�LB�LB�LB�RB�RB�RB�RB�XB�XB�XB�XB�XB�XB�XB�RB�RB�RB�RB�RB�XB�XB�^B�^B�dB�dB�dB�dB�jB�qB�qB�wB�wB�wB�wB�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�wB�wB�qB�qB�qB�wB�wB�wB�}B�}B�}B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BBÖBÖBĜBƨBƨBƨBƨBƨBƨBƨBƨBŢBŢBĜBĜBÖBB��B�}B��B��B��B�}B�}B�wB�wB�}B�}B�}B�}B�}B�}B�}B��B��B��B��B��B��B��B��BBBB��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�?B�9B�9B�9B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�3B�9B�9B�9B�9B�9B�3B�3B�3B�3B�9B�3B�3B�3B�3B�3B�9B�3B�9B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�,B�,B�,B�,B�,B�,B�,B�,B�,B�&B� B� B� B� B� B� B�B� B� B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��BֿBӬBΎB�iB�]B�]B�]B�]B�]B�VB�JB�EB�?B�?B�9B�9B�9B�3B�3B�,B�,B�,B�,B�,B�&B�,B�,B�,B�,B�&B� B� B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�uB�uB�uB�uB�uB�uB�oB�oB�jB�jB�jB�dB�dB�]B�]B�WB�QB�QB�QB�KB�?B�9B�9B�2B�2B�2B�,B�,B�&B�&B�&B�&B�&B�&B�&B� B�&B�,B�,B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�,B�2B�?B�EB�EB�9B�9B�9B�2B�9B�?B�?B�EB�KB�KB�KB�KB�QB�QB�]B�dB�jB�jB�vB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�!B�'B�'B�'B�!B�!B�!B�!B�!B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�-B�4B�4B�:B�FB�FB�FB�FB�FB�FB�FB�FB�@B�@B�:B�:B�4B�-B�!B�B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�'B�'B�'B�'B�'B�-B�-B�-B�'B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0042832                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904442018110709044420181107090444  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172711  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172711  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090444  IP  PSAL            @ffD�C3G�O�                