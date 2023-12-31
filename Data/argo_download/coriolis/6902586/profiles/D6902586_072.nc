CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T18:20:46Z creation; 2018-11-05T18:21:32Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105182046  20181107125410  6902586 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               HA   IF                                  2C  D   NOVA                            SN145                           n/a                             865 @��%] L1   @��%�l�@N9�i/n��C����0�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A���A�  A�  A�  A���A���A���A�  A�  A�  A�  A�  A�33A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BL  BP  BTffBX  B\  B`  Bd  BhffBl  Bp  Bt  Bx  B|ffB�  B���B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�33B�  B���B�  B�  B���B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C�fC� C  CffC�fC"� C$�fC'� C*  C,� C/  C1��C4�C6��C9  C;� C>  C@� CC  CE��CH�CJ��CM  COffCR  CT� CW  CY� C\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp�Cr� Cu  Cw� Cz  C|� C  C�� C�  C�33C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C��3C�33C�s3C��3C��3C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C��C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ C�s3C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�3C�  C�@ C�� C���C�  C��fC���C�  D � D  D@ D� D� D  D9�D	y�D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D� D  D@ D�fD� D   D!@ D"�fD#�fD%  D&9�D'� D(� D*fD+@ D,y�D-� D/  D09�D1y�D2��D4  D5@ D6�fD7�fD9  D:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY� DZ� D\  D]FfD^�fD_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq9�Dry�Ds��Du  Dv@ Dw� Dx� Dz  D{@ D|� D}��D  D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�c3D�  D���D�<�D���D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D��3D�` D���D�� D�@ D�� D�� D�#3D��3D�` D�3D�� D�<�D�� D��3D�  D�� D�` D�3D��3D�@ D�� D�� D��D���D�\�D�  D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�C3D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D��D���D�` D�  D��3D�@ D�� D�� D��D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�\�D���DȠ D�@ D�� Dʀ D�  D�� D�c3D�  D͠ D�@ D���Dπ D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  D�� D�` D�3Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D��3D�` D�  D�3D�C3D�� D� D�#3D�� D�` D���D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�3D�  D�� D�` D�  D� D�C3D�� D�|�D��D�� D�\�D�  D�� D�<�D�� D�� D��D�� D�c3D�3D��3D�C311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Aq��A�  A�  A���A�  A�  A�  A���A���A���A�  A�  A�  A�  A�  A�33A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BL  BP  BTffBX  B\  B`  Bd  BhffBl  Bp  Bt  Bx  B|ffB�  B���B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�33B�  B���B�  B�  B���B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C�fC� C  CffC�fC"� C$�fC'� C*  C,� C/  C1��C4�C6��C9  C;� C>  C@� CC  CE��CH�CJ��CM  COffCR  CT� CW  CY� C\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp�Cr� Cu  Cw� Cz  C|� C  C�� C�  C�33C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C��3C�33C�s3C��3C��3C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C��C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ C�s3C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�3C�  C�@ C�� C���C�  C��fC���C�  D � D  D@ D� D� D  D9�D	y�D
� D  D@ D� D� D  D@ D� D� D  D@ Dy�D� D  D@ D�fD� D   D!@ D"�fD#�fD%  D&9�D'� D(� D*fD+@ D,y�D-� D/  D09�D1y�D2��D4  D5@ D6�fD7�fD9  D:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY� DZ� D\  D]FfD^�fD_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq9�Dry�Ds��Du  Dv@ Dw� Dx� Dz  D{@ D|� D}��D  D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�c3D�  D���D�<�D���D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D��3D�` D���D�� D�@ D�� D�� D�#3D��3D�` D�3D�� D�<�D�� D��3D�  D�� D�` D�3D��3D�@ D�� D�� D��D���D�\�D�  D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�C3D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D��D���D�` D�  D��3D�@ D�� D�� D��D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�\�D���DȠ D�@ D�� Dʀ D�  D�� D�c3D�  D͠ D�@ D���Dπ D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  D�� D�` D�3Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D��3D�` D�  D�3D�C3D�� D� D�#3D�� D�` D���D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�3D�  D�� D�` D�  D� D�C3D�� D�|�D��D�� D�\�D�  D�� D�<�D�� D�� D��D�� D�c3D�3D��3D�C311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A-"�A-"�A-&�A-&�A-+A-&�A-&�A-&�A-&�A-&�A-&�A-&�A-�A-�A*VA'�A"�A!�
A!C�A�wA��A �AĜAG�A�FA�PAS�A �A��A�FA	�AA��@���@�
=@�I�@�J@��@�ff@�
=@�I�@�x�@�l�@��@��u@���@�z�@�K�@��w@���@��@��\@���@�t�@�  @��T@�I�@��@���@���@�dZ@�n�@�%@���@�bN@��@��h@�7L@�?}@���@���@�`B@���@�ƨ@���@��\@�-@���@��7@��@���@���@�?}@���@���@�
=@�
=@�33@�bN@��@�33@���@��
@��P@�
=@���@�@�@���@��h@�p�@��`@���@��@��7@���@�S�@��R@�Ĝ@���@�~�@�bN@��@�G�@��9@�1@��
@���@�V@��@���@��@���@��#@���@�I�@�Z@��@���@��P@�ȴ@�^5@�-@���@�`B@��j@��@���@�|�@��m@��m@�C�@�
=@���@�^5@��@�%@�P@~ff@|��@{�@zJ@x�`@xQ�@xr�@xbN@x�u@w�;@xb@x1'@xbN@xbN@xĜ@y��@z=q@y��@x �@v��@w
=@y%@xĜ@xb@x  @xbN@xQ�@w�@wl�@v�@u��@uV@t��@t�D@tj@s��@s�
@st�@s��@t9X@t�@tj@tZ@t�@sC�@sdZ@sS�@r�@rJ@p�@pb@p  @pb@p �@pbN@pr�@pbN@n�+@n�+@o+@n��@n��@n��@n��@n�R@n�R@n�@n�+@m@m��@m�@mp�@m`B@m�h@m��@m?}@mV@mV@m/@m/@m�@m/@l��@l9X@k�m@k�m@m@m��@m@l��@k�m@k�m@k��@j�@j��@j~�@j�\@j��@k"�@k�@k�@k�@j��@j��@j�\@j��@jM�@i��@iX@i�@h�`@h�9@h�u@h�@hr�@hA�@hb@h  @g�@h  @g�@g;d@g�@f��@g
=@g\)@g+@g��@h  @h  @g�@g�;@f�@fE�@e��@e�T@fE�@fE�@fV@e�T@e�@e�@e�T@e@e�-@e�-@e��@e��@e`B@d��@d��@d��@dI�@d�@d1@c�m@d9X@dZ@d9X@c�m@c�
@cdZ@c33@cC�@c�@ct�@c33@c@b��@b��@b^5@a�@a��@a��@b^5@b^5@b^5@b-@b=q@b~�@bM�@a�@a�#@a�#@a7L@a7L@b�@b-@a��@a7L@`bN@_l�@^�y@^v�@^v�@^5?@^{@^E�@^��@^�R@^�R@]��@]�-@]@]p�@]?}@]��@^@^5?@_;d@`1'@ax�@b��@b�H@a�7@^V@^E�@^V@]�T@^$�@^�+@^�@_;d@_
=@^�+@]p�@]p�@]��@]�-@]�-@]V@\�/@\�D@\�D@\Z@\I�@\(�@\9X@\9X@\(�@\9X@\I�@\j@\z�@\z�@\Z@\�D@\�/@\�@]V@\�/@\z�@\1@\�@[��@[�m@\1@[�F@Z^5@Z�\@[�
@]/@\�@\z�@\j@\��@]`B@\�@\z�@\Z@[�
@\�@Z�\@[dZ@[��@[�
@[�@[��@\�D@^�+@_�@`1'@c�@dZ@eV@e��@d�j@e`B@f�R@e�@e�-@f�+@gl�@hQ�@g�@g�@iG�@j�@j^5@jJ@i��@iX@h�`@i�^@j�!@j�!@j^5@i�@iX@iG�@h��@hQ�@gK�@f��@g�;@g�;@f�@fff@gl�@g|�@fȴ@f@e�h@e/@d��@d�@d��@d�D@dI�@c��@co@b~�@bM�@b-@bJ@a�^@a��@ahs@a7L@a&�@a�7@aG�@`��@`�u@`Q�@`1'@` �@_�w@_��@_\)@^��@^��@^V@^$�@]�T@]��@]�@\�/@\�D@\I�@[�@[33@[o@Z��@ZM�@Z�@Y��@Y�7@Y7L@Y�@Y%@X��@XA�@X �@W�;@W�;@W�P@W\)@W;d@W�@V��@V��@V5?@V@U�-@U��@U?}@U`B@U�@U�@Vff@V��@Vv�@V{@UO�@Sƨ@R�\@RM�@Q��@Q��@Q��@QX@PĜ@P�@Pr�@P  @O�@Ol�@O
=@N��@Nv�@N$�@M��@M�-@M�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A-"�A-"�A-&�A-&�A-+A-&�A-&�A-&�A-&�A-&�A-&�A-&�A-�A-�A*VA'�A"�A!�
A!C�A�wA��A �AĜAG�A�FA�PAS�A �A��A�FA	�AA��@���@�
=@�I�@�J@��@�ff@�
=@�I�@�x�@�l�@��@��u@���@�z�@�K�@��w@���@��@��\@���@�t�@�  @��T@�I�@��@���@���@�dZ@�n�@�%@���@�bN@��@��h@�7L@�?}@���@���@�`B@���@�ƨ@���@��\@�-@���@��7@��@���@���@�?}@���@���@�
=@�
=@�33@�bN@��@�33@���@��
@��P@�
=@���@�@�@���@��h@�p�@��`@���@��@��7@���@�S�@��R@�Ĝ@���@�~�@�bN@��@�G�@��9@�1@��
@���@�V@��@���@��@���@��#@���@�I�@�Z@��@���@��P@�ȴ@�^5@�-@���@�`B@��j@��@���@�|�@��m@��m@�C�@�
=@���@�^5@��@�%@�P@~ff@|��@{�@zJ@x�`@xQ�@xr�@xbN@x�u@w�;@xb@x1'@xbN@xbN@xĜ@y��@z=q@y��@x �@v��@w
=@y%@xĜ@xb@x  @xbN@xQ�@w�@wl�@v�@u��@uV@t��@t�D@tj@s��@s�
@st�@s��@t9X@t�@tj@tZ@t�@sC�@sdZ@sS�@r�@rJ@p�@pb@p  @pb@p �@pbN@pr�@pbN@n�+@n�+@o+@n��@n��@n��@n��@n�R@n�R@n�@n�+@m@m��@m�@mp�@m`B@m�h@m��@m?}@mV@mV@m/@m/@m�@m/@l��@l9X@k�m@k�m@m@m��@m@l��@k�m@k�m@k��@j�@j��@j~�@j�\@j��@k"�@k�@k�@k�@j��@j��@j�\@j��@jM�@i��@iX@i�@h�`@h�9@h�u@h�@hr�@hA�@hb@h  @g�@h  @g�@g;d@g�@f��@g
=@g\)@g+@g��@h  @h  @g�@g�;@f�@fE�@e��@e�T@fE�@fE�@fV@e�T@e�@e�@e�T@e@e�-@e�-@e��@e��@e`B@d��@d��@d��@dI�@d�@d1@c�m@d9X@dZ@d9X@c�m@c�
@cdZ@c33@cC�@c�@ct�@c33@c@b��@b��@b^5@a�@a��@a��@b^5@b^5@b^5@b-@b=q@b~�@bM�@a�@a�#@a�#@a7L@a7L@b�@b-@a��@a7L@`bN@_l�@^�y@^v�@^v�@^5?@^{@^E�@^��@^�R@^�R@]��@]�-@]@]p�@]?}@]��@^@^5?@_;d@`1'@ax�@b��@b�H@a�7@^V@^E�@^V@]�T@^$�@^�+@^�@_;d@_
=@^�+@]p�@]p�@]��@]�-@]�-@]V@\�/@\�D@\�D@\Z@\I�@\(�@\9X@\9X@\(�@\9X@\I�@\j@\z�@\z�@\Z@\�D@\�/@\�@]V@\�/@\z�@\1@\�@[��@[�m@\1@[�F@Z^5@Z�\@[�
@]/@\�@\z�@\j@\��@]`B@\�@\z�@\Z@[�
@\�@Z�\@[dZ@[��@[�
@[�@[��@\�D@^�+@_�@`1'@c�@dZ@eV@e��@d�j@e`B@f�R@e�@e�-@f�+@gl�@hQ�@g�@g�@iG�@j�@j^5@jJ@i��@iX@h�`@i�^@j�!@j�!@j^5@i�@iX@iG�@h��@hQ�@gK�@f��@g�;@g�;@f�@fff@gl�@g|�@fȴ@f@e�h@e/@d��@d�@d��@d�D@dI�@c��@co@b~�@bM�@b-@bJ@a�^@a��@ahs@a7L@a&�@a�7@aG�@`��@`�u@`Q�@`1'@` �@_�w@_��@_\)@^��@^��@^V@^$�@]�T@]��@]�@\�/@\�D@\I�@[�@[33@[o@Z��@ZM�@Z�@Y��@Y�7@Y7L@Y�@Y%@X��@XA�@X �@W�;@W�;@W�P@W\)@W;d@W�@V��@V��@V5?@V@U�-@U��@U?}@U`B@U�@U�@Vff@V��@Vv�@V{@UO�@Sƨ@R�\@RM�@Q��@Q��@Q��@QX@PĜ@P�@Pr�@P  @O�@Ol�@O
=@N��@Nv�@N$�@M��@M�-@M�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBVBVBVBVBPBPBPBJBDB
=B1B%BB
��B
��B
��BBBB%BB
��B
��BbB�B�B\B\B �B#�B"�B�B%�B=qBI�BR�BVBH�B[#B]/BjBjBjBm�BiyBe`BhsBp�Bn�BjBo�By�B~�B�=B��B�jB��B�)B�HB�/B�B�B�B��B��B��B��BɺB��BȴBĜBB��B��B�}B�}B��B��B��BÖBŢBƨBȴB��B�B�TB�TB�TB�mB�yB�`B�mB�yB�yB�mB�fB�ZB�ZB�TB�TB�NB�BB�)B�B�
B�B�5B�#B�B��B��BǮBÖB�wB�wB�wB��BĜB��B��BȴBƨB�wB�}B�}B�dB�dB�dB�dB�dB�^B�^B�XB�RB�LB�?B�-B�'B�9B�FB�FB�9B�3B�-B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�uB�{B��B��B��B�{B�{B�{B�{B�uB�oB�oB�hB�hB�hB�hB�bB�bB�bB�bB�bB�\B�bB�bB�\B�VB�VB�\B�bB�bB�bB�hB�hB�hB�hB�\B�VB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�JB�JB�DB�DB�=B�=B�DB�DB�DB�=B�=B�7B�7B�7B�=B�=B�7B�7B�7B�7B�+B�+B�+B�+B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�%B�%B�1B�1B�1B�%B�B�B�B� B� B~�B~�B� B�B�B� B~�B~�B~�B}�B}�B~�B� B�B�B�%B�7B�DB�DB�+B�B�B�B�B�B�B�B�B�B�B� B� B�B�B� B~�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B}�B}�B}�B}�B}�B}�B~�B~�B� B~�B}�B}�B}�B|�B|�B|�B|�By�Bz�B}�B�B� B~�B� B� B�B�B� B� B~�B~�B|�B}�B~�B~�B~�B~�B�B�1B�JB�VB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�!B�'B�!B�!B�!B�!B�3B�?B�FB�FB�FB�?B�FB�?B�9B�3B�-B�?B�FB�9B�9B�RB�XB�LB�FB�FB�?B�?B�?B�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�?B�FB�FB�LB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�XB�XB�^B�dB�jB�}B��B��B�}B�wB�^B�RB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�XB�RB�RB�RB�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 BVBVBVBVBPBPBPBJBDB
=B1B%BB
��B
��B
��BBBB%BB
��B
��BbB�B�B\B\B �B#�B"�B�B%�B=qBI�BR�BVBH�B[#B]/BjBjBjBm�BiyBe`BhsBp�Bn�BjBo�By�B~�B�=B��B�jB��B�)B�HB�/B�B�B�B��B��B��B��BɺB��BȴBĜBB��B��B�}B�}B��B��B��BÖBŢBƨBȴB��B�B�TB�TB�TB�mB�yB�`B�mB�yB�yB�mB�fB�ZB�ZB�TB�TB�NB�BB�)B�B�
B�B�5B�#B�B��B��BǮBÖB�wB�wB�wB��BĜB��B��BȴBƨB�wB�}B�}B�dB�dB�dB�dB�dB�^B�^B�XB�RB�LB�?B�-B�'B�9B�FB�FB�9B�3B�-B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�uB�{B��B��B��B�{B�{B�{B�{B�uB�oB�oB�hB�hB�hB�hB�bB�bB�bB�bB�bB�\B�bB�bB�\B�VB�VB�\B�bB�bB�bB�hB�hB�hB�hB�\B�VB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�JB�JB�DB�DB�=B�=B�DB�DB�DB�=B�=B�7B�7B�7B�=B�=B�7B�7B�7B�7B�+B�+B�+B�+B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�%B�%B�1B�1B�1B�%B�B�B�B� B� B~�B~�B� B�B�B� B~�B~�B~�B}�B}�B~�B� B�B�B�%B�7B�DB�DB�+B�B�B�B�B�B�B�B�B�B�B� B� B�B�B� B~�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B}�B}�B}�B}�B}�B}�B~�B~�B� B~�B}�B}�B}�B|�B|�B|�B|�By�Bz�B}�B�B� B~�B� B� B�B�B� B� B~�B~�B|�B}�B~�B~�B~�B~�B�B�1B�JB�VB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�!B�'B�!B�!B�!B�!B�3B�?B�FB�FB�FB�?B�FB�?B�9B�3B�-B�?B�FB�9B�9B�RB�XB�LB�FB�FB�?B�?B�?B�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�?B�FB�FB�LB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�XB�XB�^B�dB�jB�}B��B��B�}B�wB�^B�RB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�XB�RB�RB�RB�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811071254102018110712541020181107125410  IF  ARFMCODA024c                                                                20181105182046                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105182132  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105182132  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107125410  IP  PSAL            @ffD�C3G�O�                