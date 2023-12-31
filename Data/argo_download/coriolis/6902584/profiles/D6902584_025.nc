CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:12Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090444  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�8P�V��1   @�8Q��U�@Oru�_��@�U:�=�2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @�  @�  @���@���A   A��A   A1��A@  AP  A`  Ap  A~ffA�  A�33A�  A�  A�33A�33A�33A�  A�  A�  A���A���A�  A�  A�  A�33B��B��B  BffB  B��B  B   B$  B(  B,  B/��B4  B8  B<  B@ffBD  BH  BL  BPffBTffBX  B\  B`  BdffBhffBl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  C  C� C  C	��C  C� C  C� C�C� C  C� C   C"� C%  C'� C*  C,� C/  C1ffC3�fC6� C9  C;� C>  C@� CC  CE��CH  CJffCL�fCO� CR  CT� CW  CYffC\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp�Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C��3C�  C�L�C���C���C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�33C�s3C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C��3C�33C� C�� C�  C�@ C��C�� C��3C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;�fD<� D>  D?@ D@� DA� DC  DDFfDE� DF� DG��DI@ DJ� DK� DM  DN@ DO� DP� DQ��DS@ DT�fDU� DW  DX@ DYy�DZ� D\fD]@ D^y�D_��Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn��Dp  DqFfDr� Ds�fDu  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D�  D�� D�` D�  D��3D�C3D�� D�|�D��D�� D�` D�  D�� D�@ D�� D��3D�  D��3D�` D�  D��3D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D���D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�|�D�  D�� D�` D�3D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D��3D�  D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�` D�  Dà D�@ D�� Dŀ D��DƼ�D�\�D�  DȠ D�@ D�� Dʃ3D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�3Dң3D�C3D�� DԀ D�  Dռ�D�` D�  Dף3D�@ D�� Dـ D�  D��3D�` D���Dܠ D�C3D�� Dހ D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D���D� D�@ D�� D�3D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D��fD�` 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@@  @�  @�  @���@���A   A��A   A1��A@  AP  A`  Ap  A~ffA�  A�33A�  A�  A�33A�33A�33A�  A�  A�  A���A���A�  A�  A�  A�33B��B��B  BffB  B��B  B   B$  B(  B,  B/��B4  B8  B<  B@ffBD  BH  BL  BPffBTffBX  B\  B`  BdffBhffBl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  C  C� C  C	��C  C� C  C� C�C� C  C� C   C"� C%  C'� C*  C,� C/  C1ffC3�fC6� C9  C;� C>  C@� CC  CE��CH  CJffCL�fCO� CR  CT� CW  CYffC\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp�Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C��3C�  C�L�C���C���C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�33C�s3C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C��3C�33C� C�� C�  C�@ C��C�� C��3C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;�fD<� D>  D?@ D@� DA� DC  DDFfDE� DF� DG��DI@ DJ� DK� DM  DN@ DO� DP� DQ��DS@ DT�fDU� DW  DX@ DYy�DZ� D\fD]@ D^y�D_��Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn��Dp  DqFfDr� Ds�fDu  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D�  D�� D�` D�  D��3D�C3D�� D�|�D��D�� D�` D�  D�� D�@ D�� D��3D�  D��3D�` D�  D��3D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D���D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�|�D�  D�� D�` D�3D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D��3D�  D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�` D�  Dà D�@ D�� Dŀ D��DƼ�D�\�D�  DȠ D�@ D�� Dʃ3D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�3Dң3D�C3D�� DԀ D�  Dռ�D�` D�  Dף3D�@ D�� Dـ D�  D��3D�` D���Dܠ D�C3D�� Dހ D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D���D� D�@ D�� D�3D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D��fD�` 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�\)@�C�@�;d@�C�@�S�@��@��P@��@�|�@�|�@�|�@�|�@�|�@�|�@�|�@��@��@��@��@��P@���@���@���@��F@���@���@���@���@���@���@��@���@���@���@���@��@��@��@��@���@��@��@��@���@���@���@���@��F@��w@��F@��F@��
@���@��;@��m@��m@��;@��
@��
@���@��
@��;@��;@��m@��m@��m@��m@��m@��@��;@�ƨ@��m@��m@��@��m@��m@��@��@��@��
@���@��w@��F@��F@��w@��w@��w@��w@��w@��w@�ƨ@��;@��;@��
@�ƨ@�ƨ@���@��@��@��;@��@��;@��@�1@� �@�1'@�9X@�9X@�1'@�A�@�I�@�I�@�A�@�(�@�1@�  @��;@��@�b@�Q�@�r�@�bN@�Q�@�9X@�I�@�bN@�bN@�I�@��@��@�  @� �@��@���@��P@��@���@���@��@�|�@�|�@�|�@��@�|�@�t�@�l�@�l�@�\)@��R@���@���@���@�V@�=q@��@��#@��^@��^@�`B@���@��u@�j@�r�@�Q�@�1'@��@��F@�|�@�C�@��@�@��@�V@���@��`@���@���@��j@���@��@��@��@��j@��9@��9@���@���@��u@�z�@��9@��9@��j@��j@��9@��j@��9@��9@��u@���@���@���@�j@�1@�1@�  @���@��@��m@��;@��
@��
@���@���@���@��w@��@���@�\)@�;d@�33@�"�@��@�o@�@���@�5?@���@�X@�%@��@�z�@�1'@� �@���@��@�|�@�dZ@�\)@�S�@�\)@�S�@�;d@�;d@�+@�o@�@��@��y@��y@��H@��@��R@�ff@��@��@��-@�G�@��@��`@���@��@��;@�K�@�@��\@�@�&�@��9@�Ĝ@��`@���@���@��@���@��j@�&�@��@���@���@��/@��j@���@���@���@�o@�~�@�@���@�K�@���@�M�@��#@�&�@���@�r�@�1@��
@��P@�S�@�^5@�-@��@��7@���@�A�@��@~��@~�+@~{@}@}/@|�@|z�@{��@{S�@z�H@z��@z~�@yx�@yX@x��@x�@xQ�@w�;@v�@v$�@u�T@u�-@u�h@up�@u/@t�/@tZ@sC�@r�!@r��@r��@rn�@r-@q��@p�`@q%@p�`@p�`@pĜ@pA�@o�@ol�@oK�@ol�@o�P@o;d@o
=@n�@n�@n{@m��@m�@mp�@mp�@m@m�T@m�T@m@n�R@o+@o\)@o�w@p  @p  @pbN@p��@q��@q��@rM�@r��@r�H@r�@s@so@s33@sC�@s��@sƨ@sƨ@s�F@s�
@s�m@s��@t�@t(�@tI�@tz�@tj@t�D@t�D@tz�@t(�@t1@t�@s�F@s�@st�@sS�@s33@s"�@sC�@s�F@s�@sdZ@sS�@sS�@s33@s33@s33@s"�@s"�@s"�@s"�@so@r��@r�\@r�@q&�@qX@qhs@qx�@q��@q7L@p�u@pQ�@p�u@p��@p �@o�;@o��@o|�@ol�@o�@n�y@n��@n�+@n��@n��@n��@n�R@n�R@o+@o|�@o�P@o\)@o�@nV@n$�@nV@nv�@n��@n�y@o+@o;d@o
=@nV@m��@m�h@m`B@l�@l�@l�@kƨ@kC�@k33@j�@j�\@i��@i&�@i&�@i&�@h��@hA�@g��@gK�@f�+@fff@fV@f$�@e�h@eV@d�/@dz�@d�@c�m@c��@c�@co@b^5@a��@a�#@a�7@aX@`�9@`bN@`1'@_�;@_�w@_��@_l�@_
=@^��@^�y@^��@^�+@^$�@]@]�@]O�@]V@\�/@\z�@\Z@\�@[��@[33@Z�H@Z��@Zn�@Z=q@Y�^@Yx�@Y&�@Y%@X��@X�9@Xr�@Xb@W�@W�P@W
=@Vff@V{@U�T@U��@U�@U?}@T�@T��@Tz�@Tj@T9X@T1@S�F@S�@SS�@So@R�@R�\@R=q@R�@Q�@Q��@Q�7@Qhs@QG�@Q&�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�\)@�C�@�;d@�C�@�S�@��@��P@��@�|�@�|�@�|�@�|�@�|�@�|�@�|�@��@��@��@��@��P@���@���@���@��F@���@���@���@���@���@���@��@���@���@���@���@��@��@��@��@���@��@��@��@���@���@���@���@��F@��w@��F@��F@��
@���@��;@��m@��m@��;@��
@��
@���@��
@��;@��;@��m@��m@��m@��m@��m@��@��;@�ƨ@��m@��m@��@��m@��m@��@��@��@��
@���@��w@��F@��F@��w@��w@��w@��w@��w@��w@�ƨ@��;@��;@��
@�ƨ@�ƨ@���@��@��@��;@��@��;@��@�1@� �@�1'@�9X@�9X@�1'@�A�@�I�@�I�@�A�@�(�@�1@�  @��;@��@�b@�Q�@�r�@�bN@�Q�@�9X@�I�@�bN@�bN@�I�@��@��@�  @� �@��@���@��P@��@���@���@��@�|�@�|�@�|�@��@�|�@�t�@�l�@�l�@�\)@��R@���@���@���@�V@�=q@��@��#@��^@��^@�`B@���@��u@�j@�r�@�Q�@�1'@��@��F@�|�@�C�@��@�@��@�V@���@��`@���@���@��j@���@��@��@��@��j@��9@��9@���@���@��u@�z�@��9@��9@��j@��j@��9@��j@��9@��9@��u@���@���@���@�j@�1@�1@�  @���@��@��m@��;@��
@��
@���@���@���@��w@��@���@�\)@�;d@�33@�"�@��@�o@�@���@�5?@���@�X@�%@��@�z�@�1'@� �@���@��@�|�@�dZ@�\)@�S�@�\)@�S�@�;d@�;d@�+@�o@�@��@��y@��y@��H@��@��R@�ff@��@��@��-@�G�@��@��`@���@��@��;@�K�@�@��\@�@�&�@��9@�Ĝ@��`@���@���@��@���@��j@�&�@��@���@���@��/@��j@���@���@���@�o@�~�@�@���@�K�@���@�M�@��#@�&�@���@�r�@�1@��
@��P@�S�@�^5@�-@��@��7@���@�A�@��@~��@~�+@~{@}@}/@|�@|z�@{��@{S�@z�H@z��@z~�@yx�@yX@x��@x�@xQ�@w�;@v�@v$�@u�T@u�-@u�h@up�@u/@t�/@tZ@sC�@r�!@r��@r��@rn�@r-@q��@p�`@q%@p�`@p�`@pĜ@pA�@o�@ol�@oK�@ol�@o�P@o;d@o
=@n�@n�@n{@m��@m�@mp�@mp�@m@m�T@m�T@m@n�R@o+@o\)@o�w@p  @p  @pbN@p��@q��@q��@rM�@r��@r�H@r�@s@so@s33@sC�@s��@sƨ@sƨ@s�F@s�
@s�m@s��@t�@t(�@tI�@tz�@tj@t�D@t�D@tz�@t(�@t1@t�@s�F@s�@st�@sS�@s33@s"�@sC�@s�F@s�@sdZ@sS�@sS�@s33@s33@s33@s"�@s"�@s"�@s"�@so@r��@r�\@r�@q&�@qX@qhs@qx�@q��@q7L@p�u@pQ�@p�u@p��@p �@o�;@o��@o|�@ol�@o�@n�y@n��@n�+@n��@n��@n��@n�R@n�R@o+@o|�@o�P@o\)@o�@nV@n$�@nV@nv�@n��@n�y@o+@o;d@o
=@nV@m��@m�h@m`B@l�@l�@l�@kƨ@kC�@k33@j�@j�\@i��@i&�@i&�@i&�@h��@hA�@g��@gK�@f�+@fff@fV@f$�@e�h@eV@d�/@dz�@d�@c�m@c��@c�@co@b^5@a��@a�#@a�7@aX@`�9@`bN@`1'@_�;@_�w@_��@_l�@_
=@^��@^�y@^��@^�+@^$�@]@]�@]O�@]V@\�/@\z�@\Z@\�@[��@[33@Z�H@Z��@Zn�@Z=q@Y�^@Yx�@Y&�@Y%@X��@X�9@Xr�@Xb@W�@W�P@W
=@Vff@V{@U�T@U��@U�@U?}@T�@T��@Tz�@Tj@T9X@T1@S�F@S�@SS�@So@R�@R�\@R=q@R�@Q�@Q��@Q�7@Qhs@QG�@Q&�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�sB�sB�fB�fB�fB�fB�fB�fB�fB�`B�`B�ZB�ZB�NB�NB�HB�HB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�;B�;B�;B�BB�BB�BB�BB�BB�BB�;B�;B�;B�;B�;B�;B�;B�5B�5B�5B�5B�5B�5B�/B�5B�5B�/B�/B�/B�/B�/B�/B�)B�)B�)B�)B�#B�#B�#B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBƨBĜBĜBB��B�wB�qB�qB�wB�wB�wB�wB�wB�wB��B��B��B��B��B�}B�}B�wB�jB�XB�LB�9B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�oB�oB�hB�hB�hB�bB�\B�VB�PB�VB�PB�PB�PB�JB�DB�=B�7B�=B�=B�=B�=B�7B�+B�1B�1B�1B�1B�+B�%B�%B�%B�+B�+B�1B�1B�1B�1B�+B�%B�%B�%B�+B�1B�7B�7B�=B�PB�VB�bB�hB�hB�oB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�!B�!B�!B�'B�-B�3B�9B�FB�LB�RB�XB�XB�dB�jB�qB�wB�wB�wB�}B��B�}B��B�}B�qB�wB�}B��B��B��B��B��BBBBBÖBÖBĜBĜBĜBĜBŢBŢBƨBƨBȴBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȴBǮBǮBǮBǮBƨBŢBĜBÖBÖBÖBBBBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B��B��B�}B�}B�}B�}B��B�}B�}B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�B�B�B�B�B� B� B� B� B� B� B� B� B� B� B� B�B�B�B�B� B�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B� B� B� B�B�B� B�B�B�B�B�B�B�B�B� B�B�B�B�B� B�B� B� B�B� B�B� B� B�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B� B�B�B�B�B�B� B�B�B� B� B� B� B� B��B��B��B��B� B� B�B�B� B� B� B�B�B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B߫B߫BݟBݟBܙBܙBۓBۓBۓBۓBۓBۓBۓBۓBۓBۓBۓBۓBڌBڌBڌBۓBۓBۓBۓBۓBۓBڌBڌBڌBڌBڌBڌBڌBنBنBنBنBنBنB؀BنBنB؀B؀B؀B؀B؀B؀B�zB�zB�zB�zB�tB�tB�tB�nB�hB�aB�UB�OB�OB�IB�IB�IB�CB�CB�DB�DB�DB�DB�DB�DB�DB�DB�DB�>B�>B�>B�>B�>B�>B�7B�7B�+B�%B�%B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�tB�[B�PB�JB�DB�7B�1B�1B�+B�+B�+B�+B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B��B��B��B��B��B��B��B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�%B�%B�+B�1B�1B�7B�=B�DB�DB�PB�PB�VB�\B�bB�\B�bB�oB�uB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B� B�&B� B� B�B�B�B�B� B�&B�,B�2B�,B�&B� B� B�&B�&B�&B�&B� B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0045691                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904442018110709044420181107090444  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172712  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172712  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090444  IP  PSAL            @��D�` G�O�                