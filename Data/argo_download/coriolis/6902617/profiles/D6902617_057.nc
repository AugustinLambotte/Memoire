CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:21Z creation; 2018-10-26T19:00:48Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181026190021  20181119104031  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               9A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @� ��{�1   @� �"R��@SBV��(@!X?�<�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @   @Fff@�  @�33@�  @�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  B   B  B  B  B  BffB  B  B��B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BL  BP  BT  BW��B\  B`ffBd  Bh  Bl  BpffBt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�33B�33B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  C�fC� C�C	��C  C� C  C��C  C� C  C� C   C"� C%  C'� C)�fC,� C/�C1� C3�fC6ffC9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^��Ca  Cc� Cf�Ch��Ck�Cm� Cp  Cr� Ct�fCw� Cz�C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C�  C�@ Cʀ C�� C�  C�33Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C��C���C�  C�@ C��C�� C��3C�@ C� C�� C�  C�33C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D  DFfD� D� D  D@ D� D�fDfD@ D� D� D  D@ D� D� D��D!@ D"� D#� D%  D&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2� D4  D5@ D6y�D7� D9  D:@ D;� D<� D>  D?FfD@�fDA� DB��DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DOy�DP��DQ��DS9�DTy�DU��DW  DX@ DY� DZ� D\  D]@ D^� D_� DafDb@ Dc� Dd� Df  Dg@ Dh� Di� Dk  DlFfDm� Dn� Dp  Dq@ Dr� Ds� Du  DvFfDw�fDx� Dz  D{@ D|�fD}� D  D�#3D�� D�` D�3D�� D�@ D�� D�� D�#3D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D��3D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�3DȠ D�@ D��3Dʃ3D�#3D��3D�c3D�3Dͣ3D�C3D��3Dσ3D�  D��3D�c3D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D�  D�� D�c3D�  Dܜ�D�@ D�� Dހ D��D�� D�` D�  D� D�C3D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D�  D��D�<�D���D� D�  D��3D�` D�  D� D�@ D�� D� D�#3D��3D�` D�  D�� D�<�D�� D��3D�#3D��3D�` D�  D��3D�C3D��fD�Y�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @   @Fff@�  @�33@�  @�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  B   B  B  B  B  BffB  B  B��B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BL  BP  BT  BW��B\  B`ffBd  Bh  Bl  BpffBt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�33B�33B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  C�fC� C�C	��C  C� C  C��C  C� C  C� C   C"� C%  C'� C)�fC,� C/�C1� C3�fC6ffC9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^��Ca  Cc� Cf�Ch��Ck�Cm� Cp  Cr� Ct�fCw� Cz�C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C�  C�@ Cʀ C�� C�  C�33Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C��C���C�  C�@ C��C�� C��3C�@ C� C�� C�  C�33C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D  DFfD� D� D  D@ D� D�fDfD@ D� D� D  D@ D� D� D��D!@ D"� D#� D%  D&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2� D4  D5@ D6y�D7� D9  D:@ D;� D<� D>  D?FfD@�fDA� DB��DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DOy�DP��DQ��DS9�DTy�DU��DW  DX@ DY� DZ� D\  D]@ D^� D_� DafDb@ Dc� Dd� Df  Dg@ Dh� Di� Dk  DlFfDm� Dn� Dp  Dq@ Dr� Ds� Du  DvFfDw�fDx� Dz  D{@ D|�fD}� D  D�#3D�� D�` D�3D�� D�@ D�� D�� D�#3D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D��3D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�3DȠ D�@ D��3Dʃ3D�#3D��3D�c3D�3Dͣ3D�C3D��3Dσ3D�  D��3D�c3D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D�  D�� D�c3D�  Dܜ�D�@ D�� Dހ D��D�� D�` D�  D� D�C3D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D�  D��D�<�D���D� D�  D��3D�` D�  D� D�@ D�� D� D�#3D��3D�` D�  D�� D�<�D�� D��3D�#3D��3D�` D�  D��3D�C3D��fD�Y�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@n�y@n��@oK�@o��@o��@o��@o�P@ol�@ol�@ol�@oK�@oK�@oK�@oK�@oK�@o\)@o|�@o\)@oK�@o\)@o\)@o;d@o;d@o+@o+@oK�@ol�@oK�@o;d@oK�@oK�@o;d@o;d@o;d@o;d@o\)@ol�@ol�@o\)@ol�@o|�@ol�@o|�@o|�@o�P@o��@o��@o��@o��@o�P@o�P@o��@o��@o��@o��@o�@o��@o�;@o�w@o�;@o�;@o�;@o�;@o��@o��@o�;@o��@o�@o�@o�@o�w@o�w@o�w@o�w@o�w@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o��@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�@o�;@o�@o�@o�@o�@p  @o�@p  @p  @pb@pb@pb@p �@p �@p �@p1'@p1'@p1'@p1'@p1'@p1'@pA�@pA�@pQ�@pQ�@pA�@pA�@pA�@pA�@pA�@pA�@pA�@pA�@pQ�@pQ�@pQ�@pQ�@pA�@pA�@pA�@pA�@pQ�@pQ�@pA�@pbN@pbN@pr�@p�@p�u@p�u@o��@o\)@o\)@o+@n�y@m�@k��@ix�@d�/@d�D@dI�@d�@c�
@c�m@cƨ@c�m@d�D@d(�@dI�@c�F@cC�@b^5@aG�@_��@_\)@^��@]�@]��@^@]?}@\Z@[t�@Z^5@X�u@W�P@V�y@V{@T�@Sƨ@R~�@O|�@N$�@M�@K��@HA�@E�T@C�F@A�7@>��@9hs@6�y@1��@-@+t�@(bN@#�m@��@%@$�@�m@;d@`B@�@
�!@	G�@Ĝ@Q�@�P@$�@�
@n�@ b?��+?�o?��#?��?�t�?�G�?�  ?�;d?ޗ�?ۅ?���?׍P?�E�?ա�?�S�?���?�x�?�?���?�ff?��?�(�?�^5?�r�?��y?���?�ȴ?���?��?�-?�bN?�x�?���?�9X?���?�&�?��?��/?��/?��;?���?�`B?��?�G�?�`B?�Z?���?���?� �?�o?�M�?�n�?���?z�?pbN?co?T��?P�`?BM�?6�+?)x�?��?33>�X>�dZ>���>��^>hr�>N�>Kƨ>H�9>G�>)��>�>bN>�=�l�=�{=T��='�<�<�C�;��
�D�����
�����e`B����������1��9X�ȴ9��������G����#�	7L�I��bN��+����"��-V�333�9X�E�˾L�;N��N��S�ϾY��]/�aG��dZ�m�h�w�پ{�m��  ��J���˾��9��O߾��`����������-���R��Ĝ��MӾ�`B��l���xվ�����h�����33��?}��ȴ��^5��p���|���7�ě��ɺ^������`��n������ϾՁ�և+�׍P�ܬ�޸R��Ĝ��MӾ�S����T��r���{��33����X��dZ���m���ۿ Ĝ�%����J�S�����/��/���`B�`B��T�1'��ÿ
�����1�V����������;�bN��׿�`��`�-�t������+��P������H�����푿p����5?���vɿ;d�|��w� ��!G��!�7�!���"�\�#���$Z�$Z�$���%�T�'+�'��'��(1'�)7L�)xտ*���+C��,�D�,�Ϳ-V�.{�.V�.V�.��/\)�/���0 ſ0bN�0bN�0�`�1hs�1녿2�!�333�3t��3�F�3�Ͽ49X�49X�49X�4�j�4���5��6E��6�+�6ȴ�7Kǿ7�P�7�ٿ7�ٿ8Q�8���9��9���9���9�#�:��:��:^5�:^5�:���:���:�H�;"ѿ;��;�m�;�m�;�m�<(��<j�<��<j�<��<푿<푿<��=/�=p��>5?�>5?�>vɿ>�R�>�R�>�ۿ?;d�?|�?�w�?�w�?�w�@  �@  �@A��@��@Ĝ�@Ĝ�A%�A%�AG��A%�AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��A�7�A�7�A�7�A�7�A�7�A�711111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @n�y@n��@oK�@o��@o��@o��@o�P@ol�@ol�@ol�@oK�@oK�@oK�@oK�@oK�@o\)@o|�@o\)@oK�@o\)@o\)@o;d@o;d@o+@o+@oK�@ol�@oK�@o;d@oK�@oK�@o;d@o;d@o;d@o;d@o\)@ol�@ol�@o\)@ol�@o|�@ol�@o|�@o|�@o�P@o��@o��@o��@o��@o�P@o�P@o��@o��@o��@o��@o�@o��@o�;@o�w@o�;@o�;@o�;@o�;@o��@o��@o�;@o��@o�@o�@o�@o�w@o�w@o�w@o�w@o�w@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o��@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o��@o�;@o�;@o�;@o�;@o�;@o�;@o�;@o�@o�;@o�@o�@o�@o�@p  @o�@p  @p  @pb@pb@pb@p �@p �@p �@p1'@p1'@p1'@p1'@p1'@p1'@pA�@pA�@pQ�@pQ�@pA�@pA�@pA�@pA�@pA�@pA�@pA�@pA�@pQ�@pQ�@pQ�@pQ�@pA�@pA�@pA�@pA�@pQ�@pQ�@pA�@pbN@pbN@pr�@p�@p�u@p�u@o��@o\)@o\)@o+@n�y@m�@k��@ix�@d�/@d�D@dI�@d�@c�
@c�m@cƨ@c�m@d�D@d(�@dI�@c�F@cC�@b^5@aG�@_��@_\)@^��@]�@]��@^@]?}@\Z@[t�@Z^5@X�u@W�P@V�y@V{@T�@Sƨ@R~�@O|�@N$�@M�@K��@HA�@E�T@C�F@A�7@>��@9hs@6�y@1��@-@+t�@(bN@#�m@��@%@$�@�m@;d@`B@�@
�!@	G�@Ĝ@Q�@�P@$�@�
@n�@ b?��+?�o?��#?��?�t�?�G�?�  ?�;d?ޗ�?ۅ?���?׍P?�E�?ա�?�S�?���?�x�?�?���?�ff?��?�(�?�^5?�r�?��y?���?�ȴ?���?��?�-?�bN?�x�?���?�9X?���?�&�?��?��/?��/?��;?���?�`B?��?�G�?�`B?�Z?���?���?� �?�o?�M�?�n�?���?z�?pbN?co?T��?P�`?BM�?6�+?)x�?��?33>�X>�dZ>���>��^>hr�>N�>Kƨ>H�9>G�>)��>�>bN>�=�l�=�{=T��='�<�<�C�;��
�D�����
�����e`B����������1��9X�ȴ9��������G����#�	7L�I��bN��+����"��-V�333�9X�E�˾L�;N��N��S�ϾY��]/�aG��dZ�m�h�w�پ{�m��  ��J���˾��9��O߾��`����������-���R��Ĝ��MӾ�`B��l���xվ�����h�����33��?}��ȴ��^5��p���|���7�ě��ɺ^������`��n������ϾՁ�և+�׍P�ܬ�޸R��Ĝ��MӾ�S����T��r���{��33����X��dZ���m���ۿ Ĝ�%����J�S�����/��/���`B�`B��T�1'��ÿ
�����1�V����������;�bN��׿�`��`�-�t������+��P������H�����푿p����5?���vɿ;d�|��w� ��!G��!�7�!���"�\�#���$Z�$Z�$���%�T�'+�'��'��(1'�)7L�)xտ*���+C��,�D�,�Ϳ-V�.{�.V�.V�.��/\)�/���0 ſ0bN�0bN�0�`�1hs�1녿2�!�333�3t��3�F�3�Ͽ49X�49X�49X�4�j�4���5��6E��6�+�6ȴ�7Kǿ7�P�7�ٿ7�ٿ8Q�8���9��9���9���9�#�:��:��:^5�:^5�:���:���:�H�;"ѿ;��;�m�;�m�;�m�<(��<j�<��<j�<��<푿<푿<��=/�=p��>5?�>5?�>vɿ>�R�>�R�>�ۿ?;d�?|�?�w�?�w�?�w�@  �@  �@A��@��@Ĝ�@Ĝ�A%�A%�AG��A%�AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��AG��A�7�A�7�A�7�A�7�A�7�A�711111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBl�Bm�Bk�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bk�Bk�Bl�Bl�Bk�Bk�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bl�Bk�Bl�Bk�Bk�Bl�Bk�Bl�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�BjBjBiyBhsBhsBgmBe`BbNB_;B\)B\)B[#B[#B\)B\)B\)B]/B_;B`BB`BB_;B_;B]/B[#BZBZBYBXBXBYBYBXBW
BT�BR�BQ�BP�BO�BN�BM�BK�BH�BG�BE�BB�B?}B;dB8RB49B0!B)�B%�B!�B�B�B�B{BbBJB1B%BBBB%BBBBB  B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B�B�B�B�sB�ZB�TB�NB�HB�HB�BB�BB�BB�BB�;B�5B�)B�)B�#B�B�#B�HB�fB�yB�B�yB�yB�`B�ZB�mB�B�B�B�B�B�B�B�B�B�B�sB�ZB�HB�BB�#B�
B��B��B��BĜB�qB�^B�FB�?B�?B�?B�?B�?B�9B�9B�3B�-B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 Bl�Bm�Bk�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bk�Bl�Bl�Bl�Bl�Bk�Bk�Bl�Bl�Bk�Bk�Bl�Bl�Bl�Bk�Bl�Bl�Bl�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bl�Bk�Bl�Bk�Bk�Bl�Bk�Bl�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�Bk�BjBjBiyBhsBhsBgmBe`BbNB_;B\)B\)B[#B[#B\)B\)B\)B]/B_;B`BB`BB_;B_;B]/B[#BZBZBYBXBXBYBYBXBW
BT�BR�BQ�BP�BO�BN�BM�BK�BH�BG�BE�BB�B?}B;dB8RB49B0!B)�B%�B!�B�B�B�B{BbBJB1B%BBBB%BBBBB  B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B�B�B�B�sB�ZB�TB�NB�HB�HB�BB�BB�BB�BB�;B�5B�)B�)B�#B�B�#B�HB�fB�yB�B�yB�yB�`B�ZB�mB�B�B�B�B�B�B�B�B�B�B�sB�ZB�HB�BB�#B�
B��B��B��BĜB�qB�^B�FB�?B�?B�?B�?B�?B�9B�9B�3B�-B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040312018111910403120181119104031  IF  ARFMCODA024c                                                                20181026190021                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190048  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181026190048  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104031  IP  PSAL            @   D�Y�G�O�                