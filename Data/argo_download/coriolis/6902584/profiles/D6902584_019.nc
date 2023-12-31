CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  1   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:08Z last update (coriolis COQC software)   
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
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090441  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�)QI�^51   @�)R~�@N�(E�p��B�o�%�1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A.ffA@  AP  A`  Ap  A�  A�33A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�  A�33A�  B ffB  B  B  B  B  B  B��B   B$  B(  B,  B/��B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B[��B`  BdffBhffBl  Bp  Bt  Bx  B|  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�  B�  C  C� C  C	� C�C� C�fC� C  C� C�C� C   C"� C%  C'ffC*  C,� C/  C1ffC4  C6� C9  C;ffC>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Co�fCrffCu  Cw��Cz�C|��C�C�� C��3C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33Cŀ C�� C�  C�@ C�s3C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-�fD/  D0@ D1�fD2� D4  D5@ D6� D7� D9fD:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF� DH  DI@ DJy�DK� DMfDN@ DO� DP� DR  DS@ DT� DU��DV��DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di�fDk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv9�Dw� Dx� Dz  D{@ D|� D}� DfD�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D��3D�� D��D�� D�` D���D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�<�D���D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�#3D��3D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�C3D��3D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D��3D�c3D�  Dà D�@ D��3Dŀ D�#3D�� D�\�D�  Dȣ3D�@ D�� D�|�D�  D�� D�` D�  D͠ D�@ D�� Dσ3D�  D�� D�` D�  DҠ D�C3D�� DԀ D�  D�� D�` D�  Dל�D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�C3D��3D� D�  D�� D�` D�  D� D�C3D��3D�3D�  D�� D�c3D�  D� D�<�D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�` 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @ff@@  @�  @�  @�  @�  A   A  A   A.ffA@  AP  A`  Ap  A�  A�33A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�  A�33A�  B ffB  B  B  B  B  B  B��B   B$  B(  B,  B/��B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B[��B`  BdffBhffBl  Bp  Bt  Bx  B|  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�  B�  C  C� C  C	� C�C� C�fC� C  C� C�C� C   C"� C%  C'ffC*  C,� C/  C1ffC4  C6� C9  C;ffC>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Co�fCrffCu  Cw��Cz�C|��C�C�� C��3C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33Cŀ C�� C�  C�@ C�s3C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-�fD/  D0@ D1�fD2� D4  D5@ D6� D7� D9fD:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF� DH  DI@ DJy�DK� DMfDN@ DO� DP� DR  DS@ DT� DU��DV��DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di�fDk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv9�Dw� Dx� Dz  D{@ D|� D}� DfD�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D��3D�� D��D�� D�` D���D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�<�D���D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�#3D��3D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�C3D��3D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D��3D�c3D�  Dà D�@ D��3Dŀ D�#3D�� D�\�D�  Dȣ3D�@ D�� D�|�D�  D�� D�` D�  D͠ D�@ D�� Dσ3D�  D�� D�` D�  DҠ D�C3D�� DԀ D�  D�� D�` D�  Dל�D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�C3D��3D� D�  D�� D�` D�  D� D�C3D��3D�3D�  D�� D�c3D�  D� D�<�D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�` 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�X@�hs@�x�@�x�@��@��@��@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��h@��h@��h@��h@��h@��h@��h@���@���@���@���@���@���@���@���@��h@��h@��h@��7@��h@��7@���@���@���@���@���@���@��h@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��-@��^@��^@��-@��-@���@��-@��^@��^@��^@��^@��^@�@�@�@�@�@�@�@���@���@���@���@���@���@�@�@�@�@���@�@���@���@���@���@���@���@���@���@���@���@�@���@��#@���@���@���@��#@��T@��#@��^@�@�@�@��^@���@���@��7@�hs@�Ĝ@�  @�C�@�  @�{@���@�@��y@��!@�$�@�@��T@��^@�O�@��@���@��u@�1'@�b@��F@�S�@�o@���@�V@��#@�?}@��`@��@���@�Z@�(�@��w@��H@��+@�^5@��T@�@���@��^@�p�@�%@���@���@�Q�@��@���@���@�|�@�"�@�@��@���@���@�v�@�=q@�J@�@�G�@�/@��@���@��@�j@�bN@�Z@�1'@�b@�  @��m@���@��@�l�@�K�@��y@��@���@���@���@��+@�ff@�=q@�$�@��@��h@�X@�&�@�&�@��@���@�z�@�bN@�Z@�Z@�9X@��m@���@���@���@��@��@�|�@�t�@�l�@�S�@�;d@�33@���@��\@�v�@�ff@�E�@�-@�$�@��@��T@��T@���@��^@���@��@��@�p�@�hs@�`B@�?}@�/@���@��D@��@�@l�@~ȴ@~v�@~ff@~E�@~�+@~�+@~@}�-@}?}@|�@|Z@|9X@{�m@{33@z��@zJ@yG�@w�w@v��@vȴ@v��@v�R@u@up�@u?}@t��@up�@t�@t�/@t�@t�@sƨ@s�F@sS�@rn�@q��@qhs@q%@q&�@q��@q��@q�7@qx�@q&�@p�@o��@o|�@o;d@n��@o
=@o
=@o�@n��@l��@l(�@k��@l1@l�D@l1@k�@k33@j=q@i�@i%@i%@h��@h �@g�;@g�w@g|�@gl�@f�y@f��@fV@e�T@e�h@e�h@e`B@e�@d�j@dj@d9X@d�@c�m@cƨ@c�F@c�F@c��@cdZ@cS�@cS�@ct�@c�m@d�@d(�@dI�@dj@e�-@fȴ@hb@h�9@i7L@i�^@i�@jn�@j�H@ko@k��@k��@l�@l��@n@oK�@p  @qhs@rJ@rJ@q�#@q��@q�7@qx�@qx�@qhs@qhs@qx�@qhs@qX@qG�@qG�@p��@p�9@p��@p��@p�u@pr�@pbN@p1'@pb@o�@pr�@pĜ@p�`@p��@p�`@p��@pbN@p �@p  @ol�@o�@n�y@n�R@n�+@nV@nff@n5?@m@mO�@m/@m�@m/@l�@lI�@k�
@k�
@k��@k�
@k��@j�@j�@kC�@kt�@j��@k33@k33@j��@j~�@j�!@j�!@j�@i��@ix�@i7L@h��@h��@h�`@hĜ@h�@hA�@h1'@g�w@g�;@g��@g�P@gK�@f��@f�@f�+@fȴ@fff@f$�@e�@e��@e�@e/@e�@e�@d��@d��@d�j@dI�@d�@c�F@cdZ@b~�@bJ@a�7@aG�@`��@`��@`bN@` �@`  @_�@_;d@_�@^�R@^E�@]�-@]p�@]O�@]V@\�@\��@\j@[ƨ@[��@[�@Z�@Z�\@Zn�@Z^5@Zn�@Y��@Y�#@Y��@Yhs@Y�@X�`@X�@X1'@W|�@W�@W�;@W
=@V�+@Vff@Vv�@V@Up�@Up�@Up�@U`B@U�@U��@U�@Tz�@S�@S"�@R�H@R~�@Q��@Q��@Q��@QX@Q&�@Q%@P��@Pr�@P �@O�@O��@O��@O|�@O|�@Ol�@O
=@N�R@NE�@N$�@N{@M�@M�-@M/@M/@L��@L�D@LI�@L�D111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�X@�hs@�x�@�x�@��@��@��@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��7@��h@��h@��h@��h@��h@��h@��h@���@���@���@���@���@���@���@���@��h@��h@��h@��7@��h@��7@���@���@���@���@���@���@��h@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��-@��^@��^@��-@��-@���@��-@��^@��^@��^@��^@��^@�@�@�@�@�@�@�@���@���@���@���@���@���@�@�@�@�@���@�@���@���@���@���@���@���@���@���@���@���@�@���@��#@���@���@���@��#@��T@��#@��^@�@�@�@��^@���@���@��7@�hs@�Ĝ@�  @�C�@�  @�{@���@�@��y@��!@�$�@�@��T@��^@�O�@��@���@��u@�1'@�b@��F@�S�@�o@���@�V@��#@�?}@��`@��@���@�Z@�(�@��w@��H@��+@�^5@��T@�@���@��^@�p�@�%@���@���@�Q�@��@���@���@�|�@�"�@�@��@���@���@�v�@�=q@�J@�@�G�@�/@��@���@��@�j@�bN@�Z@�1'@�b@�  @��m@���@��@�l�@�K�@��y@��@���@���@���@��+@�ff@�=q@�$�@��@��h@�X@�&�@�&�@��@���@�z�@�bN@�Z@�Z@�9X@��m@���@���@���@��@��@�|�@�t�@�l�@�S�@�;d@�33@���@��\@�v�@�ff@�E�@�-@�$�@��@��T@��T@���@��^@���@��@��@�p�@�hs@�`B@�?}@�/@���@��D@��@�@l�@~ȴ@~v�@~ff@~E�@~�+@~�+@~@}�-@}?}@|�@|Z@|9X@{�m@{33@z��@zJ@yG�@w�w@v��@vȴ@v��@v�R@u@up�@u?}@t��@up�@t�@t�/@t�@t�@sƨ@s�F@sS�@rn�@q��@qhs@q%@q&�@q��@q��@q�7@qx�@q&�@p�@o��@o|�@o;d@n��@o
=@o
=@o�@n��@l��@l(�@k��@l1@l�D@l1@k�@k33@j=q@i�@i%@i%@h��@h �@g�;@g�w@g|�@gl�@f�y@f��@fV@e�T@e�h@e�h@e`B@e�@d�j@dj@d9X@d�@c�m@cƨ@c�F@c�F@c��@cdZ@cS�@cS�@ct�@c�m@d�@d(�@dI�@dj@e�-@fȴ@hb@h�9@i7L@i�^@i�@jn�@j�H@ko@k��@k��@l�@l��@n@oK�@p  @qhs@rJ@rJ@q�#@q��@q�7@qx�@qx�@qhs@qhs@qx�@qhs@qX@qG�@qG�@p��@p�9@p��@p��@p�u@pr�@pbN@p1'@pb@o�@pr�@pĜ@p�`@p��@p�`@p��@pbN@p �@p  @ol�@o�@n�y@n�R@n�+@nV@nff@n5?@m@mO�@m/@m�@m/@l�@lI�@k�
@k�
@k��@k�
@k��@j�@j�@kC�@kt�@j��@k33@k33@j��@j~�@j�!@j�!@j�@i��@ix�@i7L@h��@h��@h�`@hĜ@h�@hA�@h1'@g�w@g�;@g��@g�P@gK�@f��@f�@f�+@fȴ@fff@f$�@e�@e��@e�@e/@e�@e�@d��@d��@d�j@dI�@d�@c�F@cdZ@b~�@bJ@a�7@aG�@`��@`��@`bN@` �@`  @_�@_;d@_�@^�R@^E�@]�-@]p�@]O�@]V@\�@\��@\j@[ƨ@[��@[�@Z�@Z�\@Zn�@Z^5@Zn�@Y��@Y�#@Y��@Yhs@Y�@X�`@X�@X1'@W|�@W�@W�;@W
=@V�+@Vff@Vv�@V@Up�@Up�@Up�@U`B@U�@U��@U�@Tz�@S�@S"�@R�H@R~�@Q��@Q��@Q��@QX@Q&�@Q%@P��@Pr�@P �@O�@O��@O��@O|�@O|�@Ol�@O
=@N�R@NE�@N$�@N{@M�@M�-@M/@M/@L��@L�D@LI�@L�D111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��BBBÖBÖBÖBÖBÖBBÖBBÖBBÖBÖBÖBÖBBBBBBBBBBBBBBBBBBÖBBBÖBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBÖBBBBBBÖBÖBBÖBÖBÖBBBBBÖBBÖBÖBBÖBBBBBÖBBBÖBBBÖBÖBBBBÖBBBBBBBBBB��B��B��B��B�}B�qB��B��BĜBŢBŢBĜBĜBÖBÖBÖBB��B��B��B��B�}B�}B�}B�wB�qB�jB�dB�XB�RB�XB�XB�RB�LB�FB�?B�9B�3B�3B�-B�-B�-B�'B�'B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�uB�oB�oB�bB�VB�JB�JB�JB�JB�DB�=B�=B�=B�DB�=B�DB�DB�DB�7B�7B�1B�+B�%B�B�B�B�+B�+B�+B�+B�%B�B�B�B�B�B�B�B�B�B~�B}�B|�B|�B~�B}�B|�B{�By�Bw�Bw�Bw�Bw�Bv�Bu�Bu�Bu�Bu�Bt�Bs�Bs�Br�Br�Bq�Bq�Bq�Bp�Bp�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bp�Bq�Br�Br�Br�Bs�Bw�Bz�B~�B�B�B�B�%B�+B�7B�=B�JB�PB�VB�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�!B�!B�!B�'B�'B�-B�9B�FB�LB�RB�XB�XB�RB�RB�LB�FB�FB�FB�FB�FB�FB�LB�LB�LB�LB�LB�RB�RB�RB�LB�LB�RB�RB�XB�XB�RB�XB�^B�dB�dB�jB�jB�jB�jB�qB�wB�qB�qB�qB�qB�qB�wB�wB�wB�wB�wB�wB�wB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�wB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�wB�wB�wB�wB�wB�wB�wB�wB�}B�wB�qB�qB�wB�wB�qB�qB�qB�wB�wB�wB�}B�}B�wB�wB�wB�wB�qB�wB�wB�qB�jB�qB�qB�qB�jB�jB�qB�qB�wB�}B�wB�qB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�jB�jB�jB�jB�qB�qB�qB�qB�qB�jB�qB�qB�qB�jB�jB�qB�qB�jB�jB�w111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B³B³B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�vB�jB�dB�jB�jB�dB�^B�XB�QB�KB�EB�EB�?B�?B�?B�9B�9B�9B�3B�3B�-B�-B�-B�-B�&B�&B�&B�&B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�iB�]B�]B�]B�]B�WB�PB�PB�PB�WB�PB�WB�WB�WB�JB�KB�EB�?B�9B�3B�3B�3B�?B�?B�?B�?B�9B�3B�&B B B B B B B B|B{BzBzB|B{BzBx�Bv�Bt�Bt�Bt�Bt�Bs�Br�Br�Br�Br�Bq�Bp�Bp�Bo�Bo�Bn�Bn�Bn�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bn�Bo�Bo�Bo�Bp�Bt�Bw�B|B~B�&B�3B�9B�?B�KB�QB�^B�cB�iB�uB��B��B��B��B��B��B��B��B��B��B��B�B�	B�B�B�B�B�!B�!B�'B�.B�.B�4B�4B�4B�:B�:B�?B�KB�XB�^B�dB�jB�jB�dB�dB�^B�XB�XB�XB�XB�XB�XB�^B�^B�^B�^B�^B�dB�dB�dB�^B�^B�dB�dB�jB�jB�dB�jB�pB�vB�vB�|B�|B�|B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B��B��B��B�}B�}B��B��B��B��B��B��B�}B�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�}B�}B�}B�}B��B��B��B��B��B�}B��B��B��B�}B�}B��B��B�}B�}B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0028598                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904412018110709044120181107090441  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172708  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172708  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090441  IP  PSAL            @ffD�` G�O�                