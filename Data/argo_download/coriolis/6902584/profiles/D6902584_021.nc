CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:10Z last update (coriolis COQC software)   
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
_FillValue                 4  Bh   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  Mh   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Xh   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a4   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ch   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  nh   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w4   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �    PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �4   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �4   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �\   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �`   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �d   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �h   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �l   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �    SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �0   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �0   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �0   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090442  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�.QI�1   @�.Q��@O9y�3��AO�����1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�33@�  @�33A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B��B  BffB  B  B  B  B   B$  B(ffB,  B0  B4  B8  B<  B?��BC��BH  BL  BP  BS��BW��B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B���B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B���B�  B�33B�33B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B���C  C� C  C	� C  C� C�C� C  C� C  C� C   C"� C%  C'� C)�fC,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR�CT� CW  CY� C\�C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�@ C���C���C��C�L�C�� C��3C��3C�33C�� C�� C�  C�@ C�s3C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C��C�@ C�s3C�� C�  C�@ Cʀ C���C��C�@ Cπ C�� C��C�@ CԀ C�� C�  C�L�Cـ C�� C�  C�33Cހ C�� C��3C�@ C��C���C�  C�@ C� C�3C�  C�@ C� C�� C�  C�@ C� C���C�  C�@ C���C�� C�  C���C��D �fDfDFfD�fD� D  D@ D	� D
� DfD@ D� D�fDfD@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D59�D6� D7� D9fD:@ D;� D<� D>  D?FfD@� DA� DC  DD@ DEy�DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS9�DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm�fDn�fDpfDq@ Dr� Ds� Du  Dv@ Dw�fDx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D��3D�` D�  D��3D�@ D�� D�� D�  D���D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�3D��3D�@ D�� D��3D�#3D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�` D�  D���D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D���D�|�D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�c3D�  DҠ D�@ D�� DԀ D�  D�� D�` D���Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�<�D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D��D�� D�` D�  D��D�<�D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D�� D�` D�  D��3D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�33@�  @�33A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B��B  BffB  B  B  B  B   B$  B(ffB,  B0  B4  B8  B<  B?��BC��BH  BL  BP  BS��BW��B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B���B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B���B�  B�33B�33B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B���C  C� C  C	� C  C� C�C� C  C� C  C� C   C"� C%  C'� C)�fC,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR�CT� CW  CY� C\�C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�@ C���C���C��C�L�C�� C��3C��3C�33C�� C�� C�  C�@ C�s3C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C��C�@ C�s3C�� C�  C�@ Cʀ C���C��C�@ Cπ C�� C��C�@ CԀ C�� C�  C�L�Cـ C�� C�  C�33Cހ C�� C��3C�@ C��C���C�  C�@ C� C�3C�  C�@ C� C�� C�  C�@ C� C���C�  C�@ C���C�� C�  C���C��D �fDfDFfD�fD� D  D@ D	� D
� DfD@ D� D�fDfD@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D59�D6� D7� D9fD:@ D;� D<� D>  D?FfD@� DA� DC  DD@ DEy�DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS9�DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm�fDn�fDpfDq@ Dr� Ds� Du  Dv@ Dw�fDx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D��3D�` D�  D��3D�@ D�� D�� D�  D���D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�3D��3D�@ D�� D��3D�#3D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�` D�  D���D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D���D�|�D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�c3D�  DҠ D�@ D�� DԀ D�  D�� D�` D���Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�<�D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D��D�� D�` D�  D��D�<�D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D�� D�` D�  D��3D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��\@�~�@��\@���@���@���@���@���@���@�-@�J@�{@��@��@�J@�$�@��@�@��T@��#@��#@��#@��@��@��@�{@�@�$�@�$�@�5?@�5?@�J@��@��@��@�J@�-@�=q@�-@�5?@�=q@�=q@��@��@��@��T@��T@��#@��@�@�@�@�@��@��@��@��@��@��@���@���@��@��@��@��@��T@��T@��#@��#@���@��#@���@���@���@���@���@���@���@���@��#@��#@��#@��T@��T@���@���@���@��#@�@��-@��^@�J@��@�@��#@��7@��7@�p�@�p�@���@��#@���@��h@�7L@�hs@�?}@��@��@�&�@��@�&�@�&�@���@���@���@��9@�bN@�I�@��@�  @���@�ƨ@��y@�E�@�$�@�-@�$�@��-@��^@���@��@�$�@�=q@�^5@�v�@���@�
=@�
=@�C�@��
@�  @�b@��@��@� �@� �@�(�@�1'@� �@� �@� �@�(�@� �@�b@��m@��w@��F@��F@��F@��F@��F@��w@��F@��F@��F@��F@��F@��F@��F@��w@��w@�ƨ@�ƨ@���@���@���@�dZ@�K�@�C�@��@��y@�ȴ@��!@��+@�v�@�5?@�@�p�@��@� �@��@�J@��h@���@���@�9X@��@��@���@�v�@�5?@���@�%@��@�j@� �@���@�@��@���@���@��+@�v�@�V@���@���@��-@���@�`B@�G�@��@���@���@�Q�@��w@�t�@�dZ@��@���@�E�@�E�@�-@�{@��@��^@���@��h@�hs@�G�@�%@��9@��u@�r�@�1'@��@��@��@�  @��m@�|�@���@��H@�o@�+@�+@�;d@�n�@�J@��@�/@���@�9X@�(�@���@��F@���@���@���@��F@��@���@�t�@�dZ@�;d@�o@��y@��!@�^5@�5?@�@���@��-@��@�hs@��@��/@���@�r�@�I�@� �@�b@�  @�@��@�@+@~�@}��@}�@}��@}�-@}��@}V@|9X@{S�@z��@z^5@zJ@y�#@y�^@y�#@z��@{33@z�H@z�\@y�#@zJ@y��@y�@x�@v��@uV@t�D@tj@t(�@s�
@sS�@r�!@p�@p �@o�@n�R@n�@nff@m�-@m�@m?}@l�@k�F@k�F@m�@m��@m�-@l1@ko@k��@kS�@l�/@j~�@j�H@ko@k�@l9X@j�@j�@j�H@n�+@oK�@ko@j^5@j�\@j��@k��@lI�@lI�@lj@l��@l�/@m`B@m�-@o+@pA�@r�\@r^5@rJ@u/@u��@u�@u��@u/@uV@up�@uV@uV@u�@t��@uV@t�/@t�@t�@t�D@tj@tj@t��@t�@t�@uV@u/@t��@t�/@t�@t��@t�@uO�@up�@uV@s�
@sƨ@tZ@tz�@tZ@t�@sƨ@s�@s�@st�@s"�@r�@r��@r^5@rM�@r-@rJ@q�7@pr�@o��@o�w@o�P@pQ�@o�w@n��@n�y@n��@nv�@nV@nE�@m��@m�@lI�@l�@k�F@kdZ@kC�@kdZ@kC�@ko@j�!@j�!@jJ@jJ@jJ@i�@i��@iX@ihs@ihs@iX@i&�@h�`@h�u@h �@g|�@g�@g
=@f�y@fȴ@fv�@f{@e�@f$�@e�@e�-@e?}@e`B@eV@d�/@d�/@d�j@dz�@d9X@d9X@dI�@dI�@c��@cƨ@ct�@c��@cC�@b�@b��@b~�@a�@a�@a%@`��@aG�@a&�@`Ĝ@_��@^ȴ@^ff@]�@^5?@^V@^ff@^v�@^v�@^{@]��@]p�@]p�@]�@\��@\��@]/@]`B@]�@\�j@\I�@[�@[o@Z�@Z�\@Y��@Yhs@Y&�@Y%@X��@X�9@Xr�@X1'@W��@W�@W|�@W+@V�y@V��@V��@V�+@Vff@V5?@U@U`B@U�@T�@T�j@Tz�@T�@S�m@S�F@SC�@R��@R�\@R^5@R-@Q��@Q��@Qx�@QX@QG�@Q&�@P�`@P�911111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��\@�~�@��\@���@���@���@���@���@���@�-@�J@�{@��@��@�J@�$�@��@�@��T@��#@��#@��#@��@��@��@�{@�@�$�@�$�@�5?@�5?@�J@��@��@��@�J@�-@�=q@�-@�5?@�=q@�=q@��@��@��@��T@��T@��#@��@�@�@�@�@��@��@��@��@��@��@���@���@��@��@��@��@��T@��T@��#@��#@���@��#@���@���@���@���@���@���@���@���@��#@��#@��#@��T@��T@���@���@���@��#@�@��-@��^@�J@��@�@��#@��7@��7@�p�@�p�@���@��#@���@��h@�7L@�hs@�?}@��@��@�&�@��@�&�@�&�@���@���@���@��9@�bN@�I�@��@�  @���@�ƨ@��y@�E�@�$�@�-@�$�@��-@��^@���@��@�$�@�=q@�^5@�v�@���@�
=@�
=@�C�@��
@�  @�b@��@��@� �@� �@�(�@�1'@� �@� �@� �@�(�@� �@�b@��m@��w@��F@��F@��F@��F@��F@��w@��F@��F@��F@��F@��F@��F@��F@��w@��w@�ƨ@�ƨ@���@���@���@�dZ@�K�@�C�@��@��y@�ȴ@��!@��+@�v�@�5?@�@�p�@��@� �@��@�J@��h@���@���@�9X@��@��@���@�v�@�5?@���@�%@��@�j@� �@���@�@��@���@���@��+@�v�@�V@���@���@��-@���@�`B@�G�@��@���@���@�Q�@��w@�t�@�dZ@��@���@�E�@�E�@�-@�{@��@��^@���@��h@�hs@�G�@�%@��9@��u@�r�@�1'@��@��@��@�  @��m@�|�@���@��H@�o@�+@�+@�;d@�n�@�J@��@�/@���@�9X@�(�@���@��F@���@���@���@��F@��@���@�t�@�dZ@�;d@�o@��y@��!@�^5@�5?@�@���@��-@��@�hs@��@��/@���@�r�@�I�@� �@�b@�  @�@��@�@+@~�@}��@}�@}��@}�-@}��@}V@|9X@{S�@z��@z^5@zJ@y�#@y�^@y�#@z��@{33@z�H@z�\@y�#@zJ@y��@y�@x�@v��@uV@t�D@tj@t(�@s�
@sS�@r�!@p�@p �@o�@n�R@n�@nff@m�-@m�@m?}@l�@k�F@k�F@m�@m��@m�-@l1@ko@k��@kS�@l�/@j~�@j�H@ko@k�@l9X@j�@j�@j�H@n�+@oK�@ko@j^5@j�\@j��@k��@lI�@lI�@lj@l��@l�/@m`B@m�-@o+@pA�@r�\@r^5@rJ@u/@u��@u�@u��@u/@uV@up�@uV@uV@u�@t��@uV@t�/@t�@t�@t�D@tj@tj@t��@t�@t�@uV@u/@t��@t�/@t�@t��@t�@uO�@up�@uV@s�
@sƨ@tZ@tz�@tZ@t�@sƨ@s�@s�@st�@s"�@r�@r��@r^5@rM�@r-@rJ@q�7@pr�@o��@o�w@o�P@pQ�@o�w@n��@n�y@n��@nv�@nV@nE�@m��@m�@lI�@l�@k�F@kdZ@kC�@kdZ@kC�@ko@j�!@j�!@jJ@jJ@jJ@i�@i��@iX@ihs@ihs@iX@i&�@h�`@h�u@h �@g|�@g�@g
=@f�y@fȴ@fv�@f{@e�@f$�@e�@e�-@e?}@e`B@eV@d�/@d�/@d�j@dz�@d9X@d9X@dI�@dI�@c��@cƨ@ct�@c��@cC�@b�@b��@b~�@a�@a�@a%@`��@aG�@a&�@`Ĝ@_��@^ȴ@^ff@]�@^5?@^V@^ff@^v�@^v�@^{@]��@]p�@]p�@]�@\��@\��@]/@]`B@]�@\�j@\I�@[�@[o@Z�@Z�\@Y��@Yhs@Y&�@Y%@X��@X�9@Xr�@X1'@W��@W�@W|�@W+@V�y@V��@V��@V�+@Vff@V5?@U@U`B@U�@T�@T�j@Tz�@T�@S�m@S�F@SC�@R��@R�\@R^5@R-@Q��@Q��@Qx�@QX@QG�@Q&�@P�`@P�911111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBBBBBBBBBB��B��B��B��BB��B��B��BB��B��B��B��BB��B��BB��BB��B��B��B��B��BBB��BB��BB��B��B��BBB��B��B��BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B��B��B��B��BBĜBĜBŢBƨBǮBȴBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBǮBŢBĜBŢBĜBBB��B�}B�wB�qB�qB�jB�dB�dB�^B�XB�XB�XB�RB�XB�XB�RB�XB�XB�RB�RB�RB�RB�RB�RB�RB�LB�FB�FB�9B�9B�9B�3B�-B�'B�'B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B�{B�hB�VB�PB�PB�PB�PB�JB�=B�%B�%B�B�B�B�B�B�B�B~�B}�B~�B�B�B�B�B~�B�B� B�B~�B� B�B�B�B�B�B�B�JB�PB�B�B�B�B�+B�7B�=B�DB�JB�PB�VB�bB��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�'B�-B�3B�3B�9B�?B�?B�FB�FB�RB�^B�dB�^B�RB�RB�dB�jB�jB�jB�jB�qB�qB�qB�qB�qB�qB�qB�qB�wB�wB�qB�jB�dB�dB�jB�}B�}B�qB�wB�wB�wB�wB�wB�qB�jB�^B�dB�dB�dB�dB�dB�jB�jB�jB�qB�dB�jB�qB�qB�qB�qB�wB�wB�wB�qB�qB�wB�qB�jB�jB�qB�jB�jB�qB�qB�qB�wB�wB�qB�qB�wB�wB�wB�}B�}B�}B�}B�}B��B��B��B��B��B��B��B��B��B��B�wB�qB�qB�qB�}B�wB�qB�dB�^B�^B�^B�dB�dB�jB�jB�qB�qB�jB�dB�jB�jB�jB�jB�wB�}B�}B�}B�}B�wB�qB�qB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�qB�qB�qB�qB�qB�qB�qB�jB�dB�dB�dB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�j11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�B�B�B�B�$B�*B�0B�6B�=B�CB�IB�OB�aB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�aB�hB�hB�hB�hB�hB�hB�aB�aB�aB�aB�[B�UB�UB�UB�OB�OB�IB�=B�6B�*B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�yB�yB�B��B��B�B�yB�mB�mB�[B�OB�OB�OB�IB�IB�IB�IB�IB�OB�OB�UB�OB�OB�OB�OB�IB�IB�CB�<B�<B�<B�<B�6B�6B�0B�0B�*B�*B�$B�$B�$B�$B�$B�$B�B�$B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B~�B}�B}�B{yBzsB{yB�B��B��B}�B{yB}�B|B��B{yB|B}�B�B��B~�B}�B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B�B�6B�6B�<B�nB�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�B�B�B�B�B�B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0034291                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904422018110709044220181107090442  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172710  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172710  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090442  IP  PSAL            @ffD�� G�O�                