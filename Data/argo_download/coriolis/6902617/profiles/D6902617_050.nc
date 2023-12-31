CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:20Z creation; 2018-10-26T19:00:45Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181026190020  20181119104027  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               2A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @���!�1   @����Sa@S
���h@�xh�`1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A~ffA�  A�  A�  A�  A�  A�  A���A���A���A�  A�33A�33A�33A�  A�  B   B  B  BffBffBffBffBffB ffB$  B(  B,  B0  B4  B7��B<  B@  BDffBH  BK��BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  BtffBx  B|  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B���B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  C  C� C  C	� C�C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C.�fC1� C4  C6� C8�fC;� C>  C@� CC  CE� CH�CJ� CL�fCO� CR�CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp�Cr� Cu  Cw� Cz  C|� C~�fC��3C�  C�L�C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�s3C��3C��3C�@ C���C�� C��3C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cٌ�C�� C�  C�@ Cހ C߳3C�  C�L�C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D��D  D@ Dy�D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2�fD4  D5@ D6� D7� D9  D:FfD;� D<��D>  D?@ D@� DA� DCfDDFfDE�fDF� DHfDI@ DJ� DK� DM  DN@ DO� DP�fDR  DS@ DT� DU� DW  DX@ DY� DZ�fD\fD]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di�fDk  Dl@ Dm� Dn� Do��Dq@ Dr�fDs�fDufDv@ Dw� Dx� Dy��D{9�D|y�D}� D  D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�|�D�  D��3D�c3D�3D��3D�@ D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  Dã3D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�3Dͣ3D�@ D�� Dπ D�  D�� D�` D�  DҜ�D�<�D�� DԀ D�  D��3D�` D�  Dנ D�<�D�� Dـ D�  D�� D�` D�3Dܠ D�@ D��3Dހ D�  D߼�D�\�D�  D� D�@ D�� D� D�  D�� D�` D���D� D�@ D�� D� D�  D�� D�c3D�3D� D�<�D���D�|�D��D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A~ffA�  A�  A�  A�  A�  A�  A���A���A���A�  A�33A�33A�33A�  A�  B   B  B  BffBffBffBffBffB ffB$  B(  B,  B0  B4  B7��B<  B@  BDffBH  BK��BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  BtffBx  B|  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B���B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  C  C� C  C	� C�C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C.�fC1� C4  C6� C8�fC;� C>  C@� CC  CE� CH�CJ� CL�fCO� CR�CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp�Cr� Cu  Cw� Cz  C|� C~�fC��3C�  C�L�C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�s3C��3C��3C�@ C���C�� C��3C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cٌ�C�� C�  C�@ Cހ C߳3C�  C�L�C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D��D  D@ Dy�D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2�fD4  D5@ D6� D7� D9  D:FfD;� D<��D>  D?@ D@� DA� DCfDDFfDE�fDF� DHfDI@ DJ� DK� DM  DN@ DO� DP�fDR  DS@ DT� DU� DW  DX@ DY� DZ�fD\fD]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di�fDk  Dl@ Dm� Dn� Do��Dq@ Dr�fDs�fDufDv@ Dw� Dx� Dy��D{9�D|y�D}� D  D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�|�D�  D��3D�c3D�3D��3D�@ D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  Dã3D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�3Dͣ3D�@ D�� Dπ D�  D�� D�` D�  DҜ�D�<�D�� DԀ D�  D��3D�` D�  Dנ D�<�D�� Dـ D�  D�� D�` D�3Dܠ D�@ D��3Dހ D�  D߼�D�\�D�  D� D�@ D�� D� D�  D�� D�` D���D� D�@ D�� D� D�  D�� D�c3D�3D� D�<�D���D�|�D��D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@���@��@��9@��j@��j@�Ĝ@��j@��j@��9@��9@��j@���@���@��u@��D@��@���@��9@��j@��j@��j@��j@�Ĝ@�Ĝ@��j@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@��9@��9@��9@��9@��j@��j@��j@��j@��j@��j@��@��@��@�b@�(�@�b@�@�;@�;@�  @�1@�@�  @�  @�  @�  @�;@��@�w@��@��@��@��@��@��@��@��@�w@�w@�w@��@�@��@�@�P@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@�P@�P@�P@�P@�P@�P@�P@�P@|�@|�@|�@|�@|�@|�@|�@�P@�P@�P@�P@|�@�P@\)@;d@�@~�R@~E�@~$�@}�@}�h@}/@|��@}/@}?}@}/@|�/@|�/@|I�@{�
@{��@{��@{t�@{o@z��@z-@x��@x  @w�;@wl�@v�@v��@vE�@u/@s��@sƨ@s��@s��@s��@s��@st�@r�!@o�@m�h@i��@f@b�H@^��@Z��@Up�@Q��@P �@Nȴ@M�@Ihs@F�R@B��@A%@@1'@>v�@=/@<9X@;dZ@;"�@:��@9��@8�`@7�@6�y@5V@1��@1��@1��@1�^@/l�@0  @)�@(��@)�^@)X@'�@&�@%p�@�
@��@��@�P@�@�P@M�?��^?��
?�Q�?�|�?�V?�+?̋D?���?�1?�`B?�p�?��y?��-?�?���?��?~��?z��?x��?xb?vȴ?m�h?X��?NV?Gl�?BJ?>��?2�!?%��?�?bN?V?I�?�?
~�?	x�?`B>�v�>�ȴ>�>޸R>׍P>�t�>�bN>�C�>��>���>�ȴ>�1>��T>�/>�bN>���>fff>]/>V>0 �>�->\)>	7L=��=���=�E�=� �=��T=��=Y�='�<�`B<�9X<T���o�u��t��ě��C���P�H�9�ixս�+��C���O߽�hs���w��{��Q���������h���پo�1'�I��\)�hs����R�%�T�)��0 ž333�7KǾ;dZ�=p��A�7�D���E�˾Kƨ�P�`�Xb�["Ѿ\(��^5?�bMӾe`B�ixվk��o���u�x���{�m�����o���˾�7L������I����;��녾��Ͼ��������P���P���u�������㾜����-���w��G�������MӾ�S����
��`B��l���xվ�~���V�� ž�����-���!���j����KǾ�Q쾺�H����o�š˾�1'���;�\)���ϾՁ��
=�ۥ�ݲ-�ݲ-��5?��;d������Z���T��r����þ�xվ�D��V��h���F���j����KǾ�X���m��푾����ۿ ��%��\�o��
������˿l���9�	xտ	xտ	��
~��C����I���ͿV��������������;� ſ�׿hs�녿-�n��n��t��9X�z��j��j����?}���
=�Kǿb��u���X�X�����#����"ѿ�m�푿�-��-��-�p���-��R��ۿ;d�   � A�� A�� Ĝ�!G��!G��!G��!�7�!���"Mӿ"�\�#o�"��"��"��#o�#S��#�
�$���%��%�˿%�˿&$ݿ&��'+�'��'(r��(�9�(�ÿ)7L�)��*���,�D�-V�,�Ϳ,�Ϳ-O߿-�h�-��.V�.���.��.��.��.���.���.��.��.��/��/��/\)�/\)�/���/\)�/�;�0bN�0�`�1녿2n��2�!�2�333�333�3�Ͽ4z�4�j�5��6�6�+�6�+�7
=�7�P�7�ٿ7�ٿ8Q�8�u�9X�9��9��9��9��9��9���:���;�m�<(��<(��<j�<��<��<��<푿=/�=/�=/�<푿<�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @���@���@��@��9@��j@��j@�Ĝ@��j@��j@��9@��9@��j@���@���@��u@��D@��@���@��9@��j@��j@��j@��j@�Ĝ@�Ĝ@��j@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@��9@��9@��9@��9@��j@��j@��j@��j@��j@��j@��@��@��@�b@�(�@�b@�@�;@�;@�  @�1@�@�  @�  @�  @�  @�;@��@�w@��@��@��@��@��@��@��@��@�w@�w@�w@��@�@��@�@�P@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@l�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@|�@�P@�P@�P@�P@�P@�P@�P@�P@|�@|�@|�@|�@|�@|�@|�@�P@�P@�P@�P@|�@�P@\)@;d@�@~�R@~E�@~$�@}�@}�h@}/@|��@}/@}?}@}/@|�/@|�/@|I�@{�
@{��@{��@{t�@{o@z��@z-@x��@x  @w�;@wl�@v�@v��@vE�@u/@s��@sƨ@s��@s��@s��@s��@st�@r�!@o�@m�h@i��@f@b�H@^��@Z��@Up�@Q��@P �@Nȴ@M�@Ihs@F�R@B��@A%@@1'@>v�@=/@<9X@;dZ@;"�@:��@9��@8�`@7�@6�y@5V@1��@1��@1��@1�^@/l�@0  @)�@(��@)�^@)X@'�@&�@%p�@�
@��@��@�P@�@�P@M�?��^?��
?�Q�?�|�?�V?�+?̋D?���?�1?�`B?�p�?��y?��-?�?���?��?~��?z��?x��?xb?vȴ?m�h?X��?NV?Gl�?BJ?>��?2�!?%��?�?bN?V?I�?�?
~�?	x�?`B>�v�>�ȴ>�>޸R>׍P>�t�>�bN>�C�>��>���>�ȴ>�1>��T>�/>�bN>���>fff>]/>V>0 �>�->\)>	7L=��=���=�E�=� �=��T=��=Y�='�<�`B<�9X<T���o�u��t��ě��C���P�H�9�ixս�+��C���O߽�hs���w��{��Q���������h���پo�1'�I��\)�hs����R�%�T�)��0 ž333�7KǾ;dZ�=p��A�7�D���E�˾Kƨ�P�`�Xb�["Ѿ\(��^5?�bMӾe`B�ixվk��o���u�x���{�m�����o���˾�7L������I����;��녾��Ͼ��������P���P���u�������㾜����-���w��G�������MӾ�S����
��`B��l���xվ�~���V�� ž�����-���!���j����KǾ�Q쾺�H����o�š˾�1'���;�\)���ϾՁ��
=�ۥ�ݲ-�ݲ-��5?��;d������Z���T��r����þ�xվ�D��V��h���F���j����KǾ�X���m��푾����ۿ ��%��\�o��
������˿l���9�	xտ	xտ	��
~��C����I���ͿV��������������;� ſ�׿hs�녿-�n��n��t��9X�z��j��j����?}���
=�Kǿb��u���X�X�����#����"ѿ�m�푿�-��-��-�p���-��R��ۿ;d�   � A�� A�� Ĝ�!G��!G��!G��!�7�!���"Mӿ"�\�#o�"��"��"��#o�#S��#�
�$���%��%�˿%�˿&$ݿ&��'+�'��'(r��(�9�(�ÿ)7L�)��*���,�D�-V�,�Ϳ,�Ϳ-O߿-�h�-��.V�.���.��.��.��.���.���.��.��.��/��/��/\)�/\)�/���/\)�/�;�0bN�0�`�1녿2n��2�!�2�333�333�3�Ͽ4z�4�j�5��6�6�+�6�+�7
=�7�P�7�ٿ7�ٿ8Q�8�u�9X�9��9��9��9��9��9���:���;�m�<(��<(��<j�<��<��<��<푿=/�=/�=/�<푿<�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B`BB_;B_;B_;B_;B`BBaHB`BB`BB`BB`BB_;B_;B_;B_;B_;B_;B_;B^5B^5B]/B\)B[#B[#B[#B\)B\)B]/B]/B^5B^5B_;B^5B^5B]/B[#BW
BQ�BL�BG�BB�B<jB5?B2-B-B+B(�B$�B"�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B{BoBoBhBbB\BDB1B+B	7B
=B+B%B  B��B�B�B�B�B�B�B�sB�`B�ZB�HB�;B�#B�B�B�B��B��B��B��BɺBǮBƨBŢBŢBĜBĜB��B�}B�wB�qB�jB�jB�XB�XB�LB�FB�FB�FB�FB�?B�?B�9B�3B�3B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 BaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBbNBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B`BB_;B_;B_;B_;B`BBaHB`BB`BB`BB`BB_;B_;B_;B_;B_;B_;B_;B^5B^5B]/B\)B[#B[#B[#B\)B\)B]/B]/B^5B^5B_;B^5B^5B]/B[#BW
BQ�BL�BG�BB�B<jB5?B2-B-B+B(�B$�B"�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B{BoBoBhBbB\BDB1B+B	7B
=B+B%B  B��B�B�B�B�B�B�B�sB�`B�ZB�HB�;B�#B�B�B�B��B��B��B��BɺBǮBƨBŢBŢBĜBĜB��B�}B�wB�qB�jB�jB�XB�XB�LB�FB�FB�FB�FB�?B�?B�9B�3B�3B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040272018111910402720181119104027  IF  ARFMCODA024c                                                                20181026190020                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190045  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181026190045  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104028  IP  PSAL            @ffD�� G�O�                