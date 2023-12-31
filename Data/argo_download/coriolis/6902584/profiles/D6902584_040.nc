CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  '   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:20Z last update (coriolis COQC software)   
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
_FillValue                 (  B8   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D`   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 (  L�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O$   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  W�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 (  `\   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 (  k    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  mH   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  u�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 (  ~�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 (  �D   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �l   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �d   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �h   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �l   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �p   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �t   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �8   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �8   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �8   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �8             ,  �8Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090452  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               (A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�]�}�T�1   @�]�`V@PRF�#C�Ag��qy1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�33@�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A���A�  A�33A�33A�  A�  A�  A�  A�  A�  A�  A�33B   B  B  B  B  B  B  B  B ffB$ffB(ffB,ffB0ffB4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  BtffBxffB|ffB�33B�33B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C�fCffC�fC	� C  C� C  C� C  CffC  C� C   C"� C%  C'��C*  C,� C/�C1� C4  C6��C9  C;� C>  C@ffCC  CE� CG�fCJ� CM  CO��CR  CT� CW  CY� C[�fC^� Ca  Cc� Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C��C�@ C�s3C�� C��C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C��3C�@ CԌ�C���C�  C�@ Cٌ�C���C��C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D� DfDFfD� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2��D3��D5@ D6� D7��D9  D:@ D;� D<� D=��D?9�D@� DA� DC  DD@ DE� DF��DG��DI@ DJ� DK� DMfDNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr�fDs� Dt��Dv@ Dw� Dx� Dz  D{@ D|� D}� DfD�  D�� D�` D�  D�� D�C3D�� D�|�D��D���D�\�D���D�� D�C3D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�c3D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D��3D�C3D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�@ D��3D�� D�  D���D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D��3D�` D�  D͜�D�<�D�� Dπ D�  D�� D�` D���DҠ D�@ D�� DԀ D�  D��3D�` D�3Dף3D�<�D�� Dـ D�#3D��3D�` D�  Dܠ D�C3D�� Dހ D��D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�<�D�� D� D�#3D�� D�` D�  D� D�C3D��3D�3D�  D�� D�` D�  D� D�C3D��3D�3D�#3D�� D�i�D�#311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�33@�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A���A�  A�33A�33A�  A�  A�  A�  A�  A�  A�  A�33B   B  B  B  B  B  B  B  B ffB$ffB(ffB,ffB0ffB4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  BtffBxffB|ffB�33B�33B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C�fCffC�fC	� C  C� C  C� C  CffC  C� C   C"� C%  C'��C*  C,� C/�C1� C4  C6��C9  C;� C>  C@ffCC  CE� CG�fCJ� CM  CO��CR  CT� CW  CY� C[�fC^� Ca  Cc� Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C��C�@ C�s3C�� C��C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C��3C�@ CԌ�C���C�  C�@ Cٌ�C���C��C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D� DfDFfD� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2��D3��D5@ D6� D7��D9  D:@ D;� D<� D=��D?9�D@� DA� DC  DD@ DE� DF��DG��DI@ DJ� DK� DMfDNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr�fDs� Dt��Dv@ Dw� Dx� Dz  D{@ D|� D}� DfD�  D�� D�` D�  D�� D�C3D�� D�|�D��D���D�\�D���D�� D�C3D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�c3D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D��3D�C3D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�@ D��3D�� D�  D���D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D��3D�` D�  D͜�D�<�D�� Dπ D�  D�� D�` D���DҠ D�@ D�� DԀ D�  D��3D�` D�3Dף3D�<�D�� Dـ D�#3D��3D�` D�  Dܠ D�C3D�� Dހ D��D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�<�D�� D� D�#3D�� D�` D�  D� D�C3D��3D�3D�  D�� D�` D�  D� D�C3D��3D�3D�#3D�� D�i�D�#311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�n�@�v�@�n�@�v�@�v�@�+@旍@�\@�n�@�+@�\@�n�@�^5@�^5@�M�@�E�@�5?@�5?@�V@�M�@�V@�E�@�5?@�{@��#@噚@���@�X@�w@��m@�|�@�K�@�ȴ@�`B@�j@��
@�l�@��H@�ff@��@ݙ�@��`@�b@۶F@�;d@���@�V@؋D@�ƨ@ו�@֟�@�7L@��`@�j@��
@��@�@��@���@ҧ�@җ�@җ�@җ�@җ�@�~�@�n�@�^5@�^5@�=q@��#@ёh@�V@Ь@Л�@�z�@� �@��@�ƨ@�ƨ@Ϯ@υ@�dZ@�K�@�33@�+@��H@�~�@�-@�@�`B@�%@���@���@̣�@̣�@̓u@̃@�bN@�A�@�b@˥�@��y@�G�@Ǯ@�+@���@�v�@Ų-@�`B@�@��T@���@�ff@ư!@�v�@�=q@�{@��T@��#@ũ�@��/@ċD@�9X@�1@��
@þw@�;d@�E�@���@�X@�%@�Ĝ@�j@�I�@���@��m@��F@�@���@�^5@�E�@�{@���@�X@�Ĝ@��u@�bN@�(�@�b@�ƨ@�l�@�S�@�@���@���@�`B@��/@��@���@��P@�l�@�33@�
=@�ȴ@���@�J@��#@�@���@���@��-@���@���@��h@�hs@�&�@���@���@��@�bN@�Q�@� �@�1@��
@��w@��P@�"�@���@��+@�=q@�$�@�{@��@�@��7@�x�@�G�@��@��@�z�@�1'@���@��;@��w@�
=@�V@�@��@��h@�7L@�V@�V@�%@���@��@��/@��D@�z�@�Z@�1'@�b@��m@���@�1@��@���@��F@��P@�C�@�o@��H@��!@���@��@�Z@�1'@�1@���@��w@�t�@���@�M�@���@��@��`@��@��D@�j@�1@�+@���@�5?@���@�p�@��@��@�z�@��P@�
=@�~�@���@���@��m@�;d@�~�@�?}@�bN@�ƨ@�C�@��y@���@�-@��T@���@��7@���@��u@�r�@���@�S�@�"�@�V@���@��@��@�1'@���@�S�@��y@�5?@�{@��@��-@�X@��@� �@��w@�|�@��y@��+@�5?@�$�@��@�{@�@�@��@���@���@�`B@�/@�%@���@��9@�I�@���@��
@�ƨ@���@�t�@��@��H@��!@���@�~�@�{@��@���@�`B@��@�%@���@��/@�Q�@�33@��H@�M�@�@�X@�X@�X@�%@��D@�Z@�@�;@;d@;d@~5?@}O�@|9X@{�@z�@z-@z�@z-@z-@y��@y�@y�^@x��@xQ�@xb@xĜ@zn�@}O�@~@~V@~5?@}�-@|�@|�@|Z@|�@}/@}?}@}p�@}p�@}`B@}p�@}p�@}/@|�/@|Z@z�!@z�H@z~�@zn�@z=q@y��@y��@z=q@z�@{33@{o@{@{o@z��@y�@yhs@yG�@y&�@yG�@yhs@yx�@yG�@yhs@y��@yG�@xĜ@x��@x�`@x��@xbN@x �@x  @w�P@w�P@w+@v��@v5?@v$�@v$�@u�h@u�@t�/@t�D@tI�@t(�@s�m@s33@r^5@q�^@qX@q7L@p�@o��@o��@o��@o|�@o�@m��@m`B@m`B@m�@l��@l�/@l��@l�j@lI�@k��@k�
@k�F@kC�@k@j�@j��@j~�@jJ@i�^@i�^@i��@i��@i�^@i�^@i��@hbN@g�P@g;d@f��@e�@d�j@d1@cƨ@c��@bM�@a�7@ahs@aG�@`��@`��@`A�@`1'@_��@^�R@^E�@]@\��@[�F@[dZ@Z�H@X��@U�T@O��@Mp�@K�m@Jn�@HbN@G�P@GK�@F�y@E�h@E/@C�F@@�`@>��@=�@<�j@8��@3��@1�#@+t�@*��@(1'@%�h@$��@$�@$I�@#�F@!��@��@dZ@b@�H@��@$�@ �`?�|�?��?�I�?��9?���?�M�?�j?��?���?�A�?�7?�bN?�  ?�bN?��?���?���?�I�?�(�?�j?ܬ?���?�I�?�j11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�n�@�v�@�n�@�v�@�v�@�+@旍@�\@�n�@�+@�\@�n�@�^5@�^5@�M�@�E�@�5?@�5?@�V@�M�@�V@�E�@�5?@�{@��#@噚@���@�X@�w@��m@�|�@�K�@�ȴ@�`B@�j@��
@�l�@��H@�ff@��@ݙ�@��`@�b@۶F@�;d@���@�V@؋D@�ƨ@ו�@֟�@�7L@��`@�j@��
@��@�@��@���@ҧ�@җ�@җ�@җ�@җ�@�~�@�n�@�^5@�^5@�=q@��#@ёh@�V@Ь@Л�@�z�@� �@��@�ƨ@�ƨ@Ϯ@υ@�dZ@�K�@�33@�+@��H@�~�@�-@�@�`B@�%@���@���@̣�@̣�@̓u@̃@�bN@�A�@�b@˥�@��y@�G�@Ǯ@�+@���@�v�@Ų-@�`B@�@��T@���@�ff@ư!@�v�@�=q@�{@��T@��#@ũ�@��/@ċD@�9X@�1@��
@þw@�;d@�E�@���@�X@�%@�Ĝ@�j@�I�@���@��m@��F@�@���@�^5@�E�@�{@���@�X@�Ĝ@��u@�bN@�(�@�b@�ƨ@�l�@�S�@�@���@���@�`B@��/@��@���@��P@�l�@�33@�
=@�ȴ@���@�J@��#@�@���@���@��-@���@���@��h@�hs@�&�@���@���@��@�bN@�Q�@� �@�1@��
@��w@��P@�"�@���@��+@�=q@�$�@�{@��@�@��7@�x�@�G�@��@��@�z�@�1'@���@��;@��w@�
=@�V@�@��@��h@�7L@�V@�V@�%@���@��@��/@��D@�z�@�Z@�1'@�b@��m@���@�1@��@���@��F@��P@�C�@�o@��H@��!@���@��@�Z@�1'@�1@���@��w@�t�@���@�M�@���@��@��`@��@��D@�j@�1@�+@���@�5?@���@�p�@��@��@�z�@��P@�
=@�~�@���@���@��m@�;d@�~�@�?}@�bN@�ƨ@�C�@��y@���@�-@��T@���@��7@���@��u@�r�@���@�S�@�"�@�V@���@��@��@�1'@���@�S�@��y@�5?@�{@��@��-@�X@��@� �@��w@�|�@��y@��+@�5?@�$�@��@�{@�@�@��@���@���@�`B@�/@�%@���@��9@�I�@���@��
@�ƨ@���@�t�@��@��H@��!@���@�~�@�{@��@���@�`B@��@�%@���@��/@�Q�@�33@��H@�M�@�@�X@�X@�X@�%@��D@�Z@�@�;@;d@;d@~5?@}O�@|9X@{�@z�@z-@z�@z-@z-@y��@y�@y�^@x��@xQ�@xb@xĜ@zn�@}O�@~@~V@~5?@}�-@|�@|�@|Z@|�@}/@}?}@}p�@}p�@}`B@}p�@}p�@}/@|�/@|Z@z�!@z�H@z~�@zn�@z=q@y��@y��@z=q@z�@{33@{o@{@{o@z��@y�@yhs@yG�@y&�@yG�@yhs@yx�@yG�@yhs@y��@yG�@xĜ@x��@x�`@x��@xbN@x �@x  @w�P@w�P@w+@v��@v5?@v$�@v$�@u�h@u�@t�/@t�D@tI�@t(�@s�m@s33@r^5@q�^@qX@q7L@p�@o��@o��@o��@o|�@o�@m��@m`B@m`B@m�@l��@l�/@l��@l�j@lI�@k��@k�
@k�F@kC�@k@j�@j��@j~�@jJ@i�^@i�^@i��@i��@i�^@i�^@i��@hbN@g�P@g;d@f��@e�@d�j@d1@cƨ@c��@bM�@a�7@ahs@aG�@`��@`��@`A�@`1'@_��@^�R@^E�@]@\��@[�F@[dZ@Z�H@X��@U�T@O��@Mp�@K�m@Jn�@HbN@G�P@GK�@F�y@E�h@E/@C�F@@�`@>��@=�@<�j@8��@3��@1�#@+t�@*��@(1'@%�h@$��@$�@$I�@#�F@!��@��@dZ@b@�H@��@$�@ �`?�|�?��?�I�?��9?���?�M�?�j?��?���?�A�?�7?�bN?�  ?�bN?��?���?���?�I�?�(�?�j?ܬ?���?�I�?�j11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB:^B:^B:^B:^B:^B:^B9XB:^B;dB9XB9XB:^B:^B:^B;dB;dB<jB;dB:^B:^B:^B:^B:^B;dB<jB=qB<jB>wBJ�BG�BK�BL�BR�BbNBe`Be`BffBgmBgmBgmBgmBhsBiyBhsBhsBk�Bk�Bk�Bl�Bk�Bn�Bp�Bp�Bq�Bs�Bu�Bv�Bv�Bw�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�Bz�B{�B{�B|�B|�B|�B|�B}�B}�B|�B}�B|�B}�B|�B}�B|�B}�B}�B~�B~�B~�B~�B� B~�B~�B}�B}�B}�B|�B{�Bz�Bx�Bv�Br�Bm�Bl�Bk�BjBjBk�Bm�Bo�Bp�Br�Bu�Bu�Bu�Bt�Bt�Bt�Bs�Bq�Bp�Bo�Bo�Bn�Bn�Bm�Bk�BiyBhsBhsBgmBffBe`BdZBdZBcTBbNBaHB`BB_;B_;B^5B\)B[#BZBZBYBXBXBW
BVBT�BS�BQ�BO�BN�BK�BK�BJ�BJ�BI�BI�BH�BG�BF�BF�BF�BE�BE�BE�BE�BE�BE�BD�BC�BC�BB�BB�BA�BA�BA�B@�B@�B?}B?}B>wB=qB<jB<jB;dB;dB;dB:^B:^B9XB9XB8RB7LB6FB6FB5?B5?B49B2-B1'B0!B0!B/B.B.B.B.B.B.B.B.B-B-B-B-B,B-B-B-B,B+B+B)�B(�B'�B&�B#�B!�B �B �B�B�B�B�B�B�B�B�B�B�B�B�B{BoBbB\BVBPBDB
=B
=B1B%BBB  B��B��B��B��B�B�B�B�B�B�B�B�B�B�yB�sB�mB�`B�TB�NB�BB�5B�/B�#B�B�B�
B�B��B��B��B��B��B��B��B��B��BɺBȴBǮBǮBǮBǮBǮBǮBƨBƨBƨBŢBĜBĜBĜBÖBÖBBB��B��B��B��B�}B�}B�wB�wB�qB�qB�qB�jB�jB�jB�jB�dB�^B�LB�FB�9B�3B�-B�-B�-B�'B�!B�!B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�LB�XB�dB�dB�^B�RB�LB�XB�^B�qB�wB�}B�}B�}B�}B�}B��B��B��B�wB�wB�wB�wB��B��B��BBĜBƨBǮBǮBǮBǮBȴBȴBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBȴBȴBȴBȴBǮBǮBǮBǮBƨBƨBŢBŢBĜBÖBÖBB��B�wB�dB�XB�XB�LB�FB�FB�?B�?B�9B�9B�-B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B1DB1DB1DB1DB1DB1DB0>B1DB2JB0>B0>B1DB1DB1DB2JB2JB3PB2JB1DB1DB1DB1DB1DB2JB3PB4WB3PB5]BA�B>�BB�BC�BI�BY2B\DB\DB]JB^QB^QB^QB^QB_WB`]B_WB_WBbiBbiBbiBcoBbiBe|Bg�Bg�Bh�Bj�Bl�Bm�Bm�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bp�Bp�Bp�Bp�Bq�Br�Br�Bs�Bs�Bs�Bs�Bt�Bt�Bs�Bt�Bs�Bt�Bs�Bt�Bs�Bt�Bt�Bu�Bu�Bu�Bu�Bv�Bu�Bu�Bt�Bt�Bt�Bs�Br�Bq�Bo�Bm�Bi�BdvBcpBbjBadBadBbjBdvBf�Bg�Bi�Bl�Bl�Bl�Bk�Bk�Bk�Bj�Bh�Bg�Bf�Bf�Be}Be}BdvBbkB`_B_YB_YB^SB]LB\FB[@B[@BZ:BY4BX.BW)BV"BV"BUBSBR
BQBQBO�BN�BN�BM�BL�BK�BJ�BH�BF�BE�BB�BB�BA�BA�B@�B@�B?�B>�B=�B=�B=�B<�B<�B<�B<�B<�B<�B;�B:B:B9xB9xB8rB8rB8rB7lB7lB6fB6fB5`B4ZB3TB3TB2NB2NB2NB1HB1HB0BB0BB/<B.6B-0B-0B,)B,)B+#B)B(B'B'B&B$�B$�B$�B$�B$�B$�B$�B$�B#�B#�B#�B#�B"�B#�B#�B#�B"�B!�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BzBtBtBnBhB	\BPBJBDB>B2B+B+B�B�B�B��B��B��B��B��B�B�B�B�B�B�B�|B�vB�vB�pB�jB�dB�^B�RB�FB�@B�4B�'B�!B�B�B�B��B��B��B��B��B��B��B��B��B»B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�~B�xB�rB�rB�lB�lB�fB�fB�fB�_B�_B�_B�_B�ZB�TB�BB�<B�/B�)B�#B�#B�#B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�BB�NB�ZB�ZB�TB�HB�BB�NB�TB�gB�mB�sB�sB�sB�sB�sB�yB�yB�yB�mB�mB�mB�mB�B�B�B��B��B��B��B��B��B��B��B��B��B��B¼B¼B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B½B½B½B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�nB�\B�PB�PB�DB�>B�>B�8B�8B�2B�2B�&B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0088532                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904522018110709045320181107090452  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172720  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172720  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090453  IP  PSAL            @ffD�#3G�O�                