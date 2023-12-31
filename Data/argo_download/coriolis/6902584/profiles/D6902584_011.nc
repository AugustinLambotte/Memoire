CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:03Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090436  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�QK�1   @�Q���@N�+�u��@��!���1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   AffA   A0  AA��AP  A`  Ap  A�  A�  A�  A���A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  BffB  B  B  B  BffB   B$  B(  B,  B0  B4  B8  B<  B@ffBD  BH  BLffBP  BT  BXffB\  B`  BdffBh  Bl  Bp  Bt  BxffB|ffB�  B���B���B���B���B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�33B�33B�  B�33C  C� C  C	� C  C� C  C� C  C� C  C� C �C"� C%  C'� C*  C,��C/  C1� C4  C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C���C�  C�L�C���C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C��C�L�C���C���C�  C�@ Cŀ C�� C��C�@ Cʀ C�� C��C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��C�@ Cހ C�� C�  C�@ C� C�� C��C�@ C� C�� C�  C�33C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C���C��D � D  D@ D�fD� D  D@ D	� D
� D  D@ D� D� D  D@ Dy�D��D��D@ D� D� D��D9�D� D� D   D!@ D"y�D#� D%  D&@ D'� D(��D*  D+FfD,� D-� D/fD0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>  D?9�D@� DA� DC  DD@ DE� DF� DH  DI9�DJ� DK�fDM  DN@ DO� DP� DR  DS@ DT� DU� DWfDX@ DY� DZ�fD\  D]@ D^� D_�fDa  Db@ Dc� Dd� De��Dg9�Dh� Di� Dk  Dl@ Dmy�Dn��Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D��3D�� D�  D�� D�\�D���D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D��3D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D��D���D�\�D���D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�C3D�� D�� D�  D�� D�c3D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� D�|�D�  D�� D�\�D�  D͠ D�@ D�� Dπ D�#3D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�3D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D�3D��D�� D�` D���D��D�@ D��3D�3D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   AffA   A0  AA��AP  A`  Ap  A�  A�  A�  A���A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  BffB  B  B  B  BffB   B$  B(  B,  B0  B4  B8  B<  B@ffBD  BH  BLffBP  BT  BXffB\  B`  BdffBh  Bl  Bp  Bt  BxffB|ffB�  B���B���B���B���B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�33B�33B�  B�33C  C� C  C	� C  C� C  C� C  C� C  C� C �C"� C%  C'� C*  C,��C/  C1� C4  C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C���C�  C�L�C���C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C��C�L�C���C���C�  C�@ Cŀ C�� C��C�@ Cʀ C�� C��C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��C�@ Cހ C�� C�  C�@ C� C�� C��C�@ C� C�� C�  C�33C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C���C��D � D  D@ D�fD� D  D@ D	� D
� D  D@ D� D� D  D@ Dy�D��D��D@ D� D� D��D9�D� D� D   D!@ D"y�D#� D%  D&@ D'� D(��D*  D+FfD,� D-� D/fD0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>  D?9�D@� DA� DC  DD@ DE� DF� DH  DI9�DJ� DK�fDM  DN@ DO� DP� DR  DS@ DT� DU� DWfDX@ DY� DZ�fD\  D]@ D^� D_�fDa  Db@ Dc� Dd� De��Dg9�Dh� Di� Dk  Dl@ Dmy�Dn��Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D��3D�� D�  D�� D�\�D���D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D��3D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D��D���D�\�D���D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�C3D�� D�� D�  D�� D�c3D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� D�|�D�  D�� D�\�D�  D͠ D�@ D�� Dπ D�#3D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�3D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D�3D��D�� D�` D���D��D�@ D��3D�3D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�TA�#A�TA�A�A�A�A�A�A�A�mA�mA�A�A�A�A�A�A�A�A�A�A�A�A�A�A�AƨA�^A�^A�wAAƨA�A��AdZA�!An�A-AVA��A�uAA	�;A�A�-@��#@�J@�ƨ@�E�@�1'@���@ȣ�@ǍP@�7L@ȣ�@�1'@�  @�K�@��@˅@�&�@Гu@�
=@ͩ�@�O�@��@�7L@���@Ɂ@�33@��@���@�1'@�5?@ə�@���@ɡ�@ɺ^@�j@�^5@�ff@ȴ9@���@ȋD@��#@��#@�ff@���@���@�-@�hs@��@���@�Ĝ@�A�@�  @���@�ƨ@���@�K�@��@�1@�I�@��F@�t�@��@�$�@���@�&�@��@�Q�@���@�ƨ@��P@�K�@�@�^5@���@��@��`@��D@���@�l�@�@�n�@��@�@�G�@�%@�z�@�1@��F@�@��H@��@��R@�^5@���@��^@�G�@�V@�z�@�1@��@��@��@���@���@�-@���@�G�@���@�I�@���@��@�\)@�dZ@�t�@�\)@�C�@�33@���@���@��@��@��9@�I�@�t�@�ȴ@���@��h@�?}@��`@�Q�@�t�@�^5@���@�&�@�bN@��@���@�;d@���@��R@�M�@��h@��@���@���@��@���@�dZ@�"�@���@�t�@��P@���@��P@�C�@��@���@�"�@��R@�n�@�M�@��@�p�@�/@���@��j@�1@��w@�C�@�"�@��@���@�ff@�@��@��@��j@�bN@��@�\)@��@���@�V@���@�`B@��@�bN@�1@�dZ@���@�v�@�=q@��#@��h@�p�@�?}@��@���@�z�@�I�@�1'@�b@���@�33@�@��\@�n�@�-@��@���@��-@�p�@�?}@���@��D@�A�@��;@�ƨ@�|�@�;d@��y@���@�v�@�=q@�J@��^@���@�hs@�&�@��/@��@�r�@�A�@��;@��m@�t�@��F@���@�t�@�@�n�@�V@�V@�E�@�V@��!@��@��@���@�ȴ@���@��\@�v�@�^5@�5?@�J@��@���@��^@��7@�hs@�O�@�/@�%@��`@���@��u@�z�@�Z@�A�@�1@�@K�@�@~ȴ@~v�@~V@~E�@}@}`B@}V@|�/@|9X@{ƨ@{dZ@{"�@{o@z~�@y��@yhs@y&�@x��@xĜ@x�9@x��@x1'@xb@w��@w��@wK�@w�@x��@y&�@y�@x�9@x1'@w��@wl�@w\)@w�@v�@vV@v{@u`B@t��@tI�@t(�@s�
@s��@st�@s�F@uV@vȴ@vȴ@vV@vE�@tz�@t9X@t�@t��@t��@t��@t��@t��@tI�@u/@t��@s�m@s33@r��@r��@st�@st�@sS�@s33@s33@sC�@sS�@s��@s�
@t(�@tI�@tZ@t��@t�@w��@x�9@y7L@y��@xbN@vff@u�@v@v@u��@u��@v@v��@v��@v�y@vff@u�h@u`B@u��@u�@t��@t�@uV@t�j@tj@t(�@t�@t1@s�F@s�@s��@s�@st�@st�@s�
@s�F@s��@s�@sC�@r�H@so@s@r�!@r=q@r�@q��@r�\@r��@q�^@q%@pbN@p �@p  @o�w@o�w@o��@oK�@o;d@oK�@o�@n�+@nV@n��@n��@n�R@m�T@m��@m@m`B@m�@m?}@l�@k��@k��@kƨ@k�
@k��@j�@j��@jn�@j~�@j~�@i��@i��@iX@iX@i�@hĜ@h�u@hr�@hA�@hb@g��@gK�@gK�@g+@fȴ@f��@fff@f@e��@e�@e?}@d�/@d�j@d�D@d9X@c��@c�F@ct�@c"�@b��@c@b��@bn�@bM�@a�#@a�7@a7L@a%@`�@`bN@`1'@`  @_��@_�P@_K�@_
=@^�@^ȴ@^�+@^V@^5?@^V@^{@\�/@\�@\z�@\9X@\�@[��@[�
@[��@[dZ@[33@[o@[o@Z�!@Z�\@ZM�@Zn�@[33@[S�@[o@Z�@Z�!@YG�@XĜ@Y7L@Y%@Y%@X�@XbN11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A�TA�#A�TA�A�A�A�A�A�A�A�mA�mA�A�A�A�A�A�A�A�A�A�A�A�A�A�A�AƨA�^A�^A�wAAƨA�A��AdZA�!An�A-AVA��A�uAA	�;A�A�-@��#@�J@�ƨ@�E�@�1'@���@ȣ�@ǍP@�7L@ȣ�@�1'@�  @�K�@��@˅@�&�@Гu@�
=@ͩ�@�O�@��@�7L@���@Ɂ@�33@��@���@�1'@�5?@ə�@���@ɡ�@ɺ^@�j@�^5@�ff@ȴ9@���@ȋD@��#@��#@�ff@���@���@�-@�hs@��@���@�Ĝ@�A�@�  @���@�ƨ@���@�K�@��@�1@�I�@��F@�t�@��@�$�@���@�&�@��@�Q�@���@�ƨ@��P@�K�@�@�^5@���@��@��`@��D@���@�l�@�@�n�@��@�@�G�@�%@�z�@�1@��F@�@��H@��@��R@�^5@���@��^@�G�@�V@�z�@�1@��@��@��@���@���@�-@���@�G�@���@�I�@���@��@�\)@�dZ@�t�@�\)@�C�@�33@���@���@��@��@��9@�I�@�t�@�ȴ@���@��h@�?}@��`@�Q�@�t�@�^5@���@�&�@�bN@��@���@�;d@���@��R@�M�@��h@��@���@���@��@���@�dZ@�"�@���@�t�@��P@���@��P@�C�@��@���@�"�@��R@�n�@�M�@��@�p�@�/@���@��j@�1@��w@�C�@�"�@��@���@�ff@�@��@��@��j@�bN@��@�\)@��@���@�V@���@�`B@��@�bN@�1@�dZ@���@�v�@�=q@��#@��h@�p�@�?}@��@���@�z�@�I�@�1'@�b@���@�33@�@��\@�n�@�-@��@���@��-@�p�@�?}@���@��D@�A�@��;@�ƨ@�|�@�;d@��y@���@�v�@�=q@�J@��^@���@�hs@�&�@��/@��@�r�@�A�@��;@��m@�t�@��F@���@�t�@�@�n�@�V@�V@�E�@�V@��!@��@��@���@�ȴ@���@��\@�v�@�^5@�5?@�J@��@���@��^@��7@�hs@�O�@�/@�%@��`@���@��u@�z�@�Z@�A�@�1@�@K�@�@~ȴ@~v�@~V@~E�@}@}`B@}V@|�/@|9X@{ƨ@{dZ@{"�@{o@z~�@y��@yhs@y&�@x��@xĜ@x�9@x��@x1'@xb@w��@w��@wK�@w�@x��@y&�@y�@x�9@x1'@w��@wl�@w\)@w�@v�@vV@v{@u`B@t��@tI�@t(�@s�
@s��@st�@s�F@uV@vȴ@vȴ@vV@vE�@tz�@t9X@t�@t��@t��@t��@t��@t��@tI�@u/@t��@s�m@s33@r��@r��@st�@st�@sS�@s33@s33@sC�@sS�@s��@s�
@t(�@tI�@tZ@t��@t�@w��@x�9@y7L@y��@xbN@vff@u�@v@v@u��@u��@v@v��@v��@v�y@vff@u�h@u`B@u��@u�@t��@t�@uV@t�j@tj@t(�@t�@t1@s�F@s�@s��@s�@st�@st�@s�
@s�F@s��@s�@sC�@r�H@so@s@r�!@r=q@r�@q��@r�\@r��@q�^@q%@pbN@p �@p  @o�w@o�w@o��@oK�@o;d@oK�@o�@n�+@nV@n��@n��@n�R@m�T@m��@m@m`B@m�@m?}@l�@k��@k��@kƨ@k�
@k��@j�@j��@jn�@j~�@j~�@i��@i��@iX@iX@i�@hĜ@h�u@hr�@hA�@hb@g��@gK�@gK�@g+@fȴ@f��@fff@f@e��@e�@e?}@d�/@d�j@d�D@d9X@c��@c�F@ct�@c"�@b��@c@b��@bn�@bM�@a�#@a�7@a7L@a%@`�@`bN@`1'@`  @_��@_�P@_K�@_
=@^�@^ȴ@^�+@^V@^5?@^V@^{@\�/@\�@\z�@\9X@\�@[��@[�
@[��@[dZ@[33@[o@[o@Z�!@Z�\@ZM�@Zn�@[33@[S�@[o@Z�@Z�!@YG�@XĜ@Y7L@Y%@Y%@X�@XbN11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB>wB>wB<jB=qB=qB=qB=qB=qB=qB=qB>wB>wB=qB=qB=qB=qB=qB=qB=qB=qB=qB=qB=qB=qB=qB<jB<jB=qB>wB=qB<jB;dB:^B:^B9XB;dBB�BA�B>wB6FB1'BR�Br�Br�B��B��B�B��BƨB�FB��B��B�LB�RB��BB�}BÖBŢBŢB�B��B��B��B�B�B�B�B�B�B�HB�`B�B�sB�B��B�B��B��B��B�B�B��B��B��B��B��B��B  BB��B��B��B��B��B��B�B�`B�sB�B�B�sB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�B�B�B�B�B�B�B�yB�sB�sB�mB�mB�mB�mB�mB�fB�fB�`B�TB�ZB�ZB�`B�fB�sB�yB�B�B�B�yB�sB�fB�`B�ZB�HB�;B�)B�B�B�B�
B��B��B��B��B��B��B��B��BɺBɺBɺBǮBĜBÖBÖBĜBŢBÖBĜBĜBǮBɺB��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBǮBƨBŢBĜBĜBÖBÖBB��B��B��B��B�}B�qB�jB�jB�dB�^B�XB�RB�LB�FB�9B�3B�-B�-B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�oB�oB�hB�hB�bB�\B�\B�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�\B�oB�uB�uB�oB�hB�hB�bB�bB�bB�bB�\B�VB�PB�JB�DB�DB�DB�DB�=B�JB�bB�{B�{B�{B�{B�bB�bB�bB�oB�uB�uB�uB�uB�uB��B��B�{B�oB�oB�oB�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�'B�'B�'B�-B�3B�9B�9B�9B�9B�9B�9B�?B�?B�9B�9B�?B�LB�RB�FB�?B�?B�?B�?B�?B�FB�FB�FB�LB�RB�LB�LB�LB�XB�^B�dB�XB�^B�dB�dB�jB�jB�dB�^B�^B�dB�jB�jB�^B�^B�^B�dB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�^B�^B�^B�^B�^B�XB�XB�XB�XB�^B�^B�^B�^B�^B�^B�XB�XB�^B�^B�XB�^B�XB�XB�XB�XB�RB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�^B�^B�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�^B�dB�wB�}B�}B�}B�wB�dB�dB�qB�qB�qB�qB�q11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B=�B=�B;�B<�B<�B<�B<�B<�B<�B<�B=�B=�B<�B<�B<�B<�B<�B<�B<�B<�B<�B<�B<�B<�B<�B;�B;�B<�B=�B<�B;�B:�B9�B9�B8�B:�BA�B@�B=�B5�B0�BRXBrBrB��B�aB�yB�VB�B��B�WB��B��B��B��B��B��B��B�B�B�uB�@B�SB�.B�B�B�B�B�B��B�B��B�B��B�B�4B�B�4B�.B�.B�	B�	B�YB�_B�SB�LB�4B�YB�eB kB�_B�SB�4B�_B�FB�.B�	B��B��B��B��B��B�@B�FB�FB�FB�LB�FB�FB�@B�@B�@B�@B�@B�@B�@B�:B�4B�.B�(B�(B�!B�B�B�	B�	B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B��B��B��B��B��B��B��B�B�BޠBۏBكBكB�vB�pB�^B�KB�?B�9B�3B�'B�'B�'B� B� B� B�B�B��B��B�B�B��B�B�B�B� B�'B�3B�3B�-B�-B�9B�3B�-B�-B�-B�'B� B� B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB��B��B�tB�tB�nB�hB�bB�bB�\B�\B�\B�\B�VB�VB�OB�OB�IB�OB�IB�CB�=B�=B�=B�=B�=B�7B�7B�1B�1B�+B�+B�$B�$B�B�B�B�B�B�B�B� B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�B�B�B�B�B�OB�hB�tB�zB�hB�OB�IB�OB�OB�VB�VB�\B�hB�tB�tB�tB�hB�hB�tB�nB�tB�tB�zB�zB�tB�tB�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.00058773                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904362018110709043620181107090436  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172703  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172703  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090436  IP  PSAL            @ffD��3G�O�                