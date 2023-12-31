CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:05Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090437  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�P���1   @�Q2@yX@N��&]a��B����1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A���A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�33A�  A�  B   BffB  B  BffB  B  B  B   B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BLffBP  BS��BX  B\  B`  Bd  Bh  Bl  Bp  BtffBxffB|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B���B�  B�33C�C��C  C	� C  C� C  C� C  C��C�C��C �C"��C%  C'� C*  C,��C/�C1� C4  C6� C9  C;� C>  C@� CC  CE� CG�fCJffCL�fCO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck�Cm� Cp  Cr� Cu�Cw� Cz  C|� C  C�� C��3C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C��C�L�Cʌ�C�� C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C�  C�@ C� C�� C�  C�33C� C�� C�  C�@ C� C�� C�  C�33C� C�� C�  C�33C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	�fD
� D  DFfD�fD� D  D@ D� D� D  D@ Dy�D� D  D@ D� D��D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:FfD;� D<� D>  D?@ D@� DA�fDC  DD@ DE� DF� DH  DI@ DJy�DK� DM  DN@ DOy�DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx�fDz  D{@ D|� D}� DfD�#3D�� D�c3D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�c3D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʃ3D�#3D��3D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  D��3D�` D�3Dנ D�<�D�� Dـ D�  D�� D�` D�  Dܠ D�<�D���Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D���D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�c3D�3D��3D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A���A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�33A�  A�  B   BffB  B  BffB  B  B  B   B$  B(  B,  B0  B4  B8  B<ffB@  BD  BH  BLffBP  BS��BX  B\  B`  Bd  Bh  Bl  Bp  BtffBxffB|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B���B�  B�33C�C��C  C	� C  C� C  C� C  C��C�C��C �C"��C%  C'� C*  C,��C/�C1� C4  C6� C9  C;� C>  C@� CC  CE� CG�fCJffCL�fCO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck�Cm� Cp  Cr� Cu�Cw� Cz  C|� C  C�� C��3C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C��C�L�Cʌ�C�� C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cތ�C�� C�  C�@ C� C�� C�  C�33C� C�� C�  C�@ C� C�� C�  C�33C� C�� C�  C�33C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	�fD
� D  DFfD�fD� D  D@ D� D� D  D@ Dy�D� D  D@ D� D��D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:FfD;� D<� D>  D?@ D@� DA�fDC  DD@ DE� DF� DH  DI@ DJy�DK� DM  DN@ DOy�DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx�fDz  D{@ D|� D}� DfD�#3D�� D�c3D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�c3D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D���D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʃ3D�#3D��3D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  D��3D�` D�3Dנ D�<�D�� Dـ D�  D�� D�` D�  Dܠ D�<�D���Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D���D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D��3D�c3D�3D��3D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��P@��P@���@���@���@���@���@���@��P@��@��P@�S�@�+@��@���@�E�@���@�?}@�7L@��@��`@�A�@��@�
=@���@��@���@�@��h@�x�@�p�@�p�@�hs@�hs@�hs@�hs@�hs@�p�@�p�@�hs@�`B@�hs@�hs@�hs@�hs@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�hs@�hs@�`B@�X@�O�@�7L@�Ĝ@�Q�@�F@�o@��@�p�@��y@���@�z�@��@�(�@㝲@��T@�@�=q@��@�?}@ӝ�@�t�@ϥ�@�1@ϥ�@ϥ�@��/@Ѳ-@Ѳ-@�hs@�X@�&�@���@мj@�9X@���@͑h@��@�G�@́@́@�X@��@�bN@���@���@��@�;d@Ɨ�@š�@�G�@���@�1@��y@�ff@��T@��@�ƨ@��@��7@�@��@�9X@��m@�t�@���@�|�@���@�|�@�5?@�X@�&�@��w@�dZ@�\)@�"�@�ff@�M�@�ff@���@�Ĝ@��w@��@���@�v�@�ff@�-@��@�-@��T@���@�1'@�ƨ@�\)@�;d@�C�@�{@�p�@�%@���@� �@�l�@�o@���@��7@�;d@�M�@�=q@��@��-@���@�Z@�1'@��@��w@��@�K�@��@���@��!@���@�$�@�O�@�  @�M�@�J@�G�@�Ĝ@���@�z�@�A�@�1'@��F@�C�@��!@�J@�G�@��@���@���@��/@���@��9@�Q�@�  @�ƨ@�ƨ@���@�;d@�\)@�+@��y@��R@���@���@�n�@�$�@�{@�J@��@�J@���@���@��@��u@�Q�@�1@���@��m@��F@�"�@���@��\@��+@�M�@�J@���@��@��j@�r�@��m@��@�K�@���@�ff@�{@�@�@�p�@���@���@�r�@�1'@��D@��u@�j@�Q�@��@�dZ@�;d@���@�ff@�5?@���@��@�z�@�C�@��H@�n�@���@�G�@���@�z�@�b@��;@��@�dZ@�S�@��@���@�ȴ@�n�@��#@��h@��@�G�@�/@�7L@�V@���@��D@�(�@�;@��@��@\)@~ȴ@~��@~v�@~ff@}��@}�@}V@|�/@{ƨ@{"�@z~�@z^5@z�@y��@y�@x��@xbN@x1'@wK�@w
=@v�R@vV@u��@u�@t��@t1@t1@s��@s�@sC�@r�@s"�@sdZ@s"�@r�H@r��@r~�@rJ@q��@q�^@q�7@qX@q7L@p��@o��@n��@n��@n�y@o�@o�w@o��@o��@o�@pr�@p�u@pA�@o�@p��@q7L@qG�@qhs@r-@r��@s"�@s33@s"�@s"�@sS�@st�@sdZ@s"�@s33@s��@t9X@t�@u�@u��@v$�@u�-@u`B@uO�@u�@u�@uO�@u?}@u?}@u/@u�@u�@u?}@uO�@u`B@v@vff@v�+@vff@v�R@vȴ@vȴ@v�+@v��@v�@v�R@v�@v�@v��@u�@uO�@t�D@s�F@r�@r��@rn�@rM�@rJ@q��@qx�@qX@qG�@q&�@q�@q7L@q7L@q�@q�@q%@q&�@q�@rM�@qX@r�@s�m@s�m@st�@r��@r~�@r-@rJ@q��@qG�@p��@o�@o;d@o�@n�y@n��@o�w@ol�@o|�@n�y@n�R@n��@nȴ@n5?@m�@l�@lI�@l1@lZ@m�@m�@l��@l��@l��@k�
@j~�@i�#@i��@i��@i�7@i&�@h�9@g�w@f�@fE�@e�@e��@ep�@eO�@e`B@ep�@eV@d��@d�D@d1@c�
@cƨ@c��@c��@cS�@c33@co@b��@bn�@b-@bJ@a��@aX@a�@`r�@`  @_�w@_�P@^�+@]��@]@]��@]�@\�D@\�@[��@[ƨ@[dZ@["�@Z�H@Z��@Z~�@Z=q@Z�@Y�@Y�@Y��@Y��@Y�7@YX@Y7L@Y%@Y�@Y&�@Y%@Y%@XĜ@XbN@X  @W�@Wl�@W+@V��@V�y@Vv�@VE�@V@U��@U@Up�@UV@T�@TZ@TI�@T�D@Tj@T1@Sƨ@SdZ@R�\@Qhs@P�`@Pr�@P1'@O�@P �@O�;@O�w@Ol�@O;d@O
=@N��@NV@O
=@O|�@Ol�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��P@��P@���@���@���@���@���@���@��P@��@��P@�S�@�+@��@���@�E�@���@�?}@�7L@��@��`@�A�@��@�
=@���@��@���@�@��h@�x�@�p�@�p�@�hs@�hs@�hs@�hs@�hs@�p�@�p�@�hs@�`B@�hs@�hs@�hs@�hs@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�hs@�hs@�`B@�X@�O�@�7L@�Ĝ@�Q�@�F@�o@��@�p�@��y@���@�z�@��@�(�@㝲@��T@�@�=q@��@�?}@ӝ�@�t�@ϥ�@�1@ϥ�@ϥ�@��/@Ѳ-@Ѳ-@�hs@�X@�&�@���@мj@�9X@���@͑h@��@�G�@́@́@�X@��@�bN@���@���@��@�;d@Ɨ�@š�@�G�@���@�1@��y@�ff@��T@��@�ƨ@��@��7@�@��@�9X@��m@�t�@���@�|�@���@�|�@�5?@�X@�&�@��w@�dZ@�\)@�"�@�ff@�M�@�ff@���@�Ĝ@��w@��@���@�v�@�ff@�-@��@�-@��T@���@�1'@�ƨ@�\)@�;d@�C�@�{@�p�@�%@���@� �@�l�@�o@���@��7@�;d@�M�@�=q@��@��-@���@�Z@�1'@��@��w@��@�K�@��@���@��!@���@�$�@�O�@�  @�M�@�J@�G�@�Ĝ@���@�z�@�A�@�1'@��F@�C�@��!@�J@�G�@��@���@���@��/@���@��9@�Q�@�  @�ƨ@�ƨ@���@�;d@�\)@�+@��y@��R@���@���@�n�@�$�@�{@�J@��@�J@���@���@��@��u@�Q�@�1@���@��m@��F@�"�@���@��\@��+@�M�@�J@���@��@��j@�r�@��m@��@�K�@���@�ff@�{@�@�@�p�@���@���@�r�@�1'@��D@��u@�j@�Q�@��@�dZ@�;d@���@�ff@�5?@���@��@�z�@�C�@��H@�n�@���@�G�@���@�z�@�b@��;@��@�dZ@�S�@��@���@�ȴ@�n�@��#@��h@��@�G�@�/@�7L@�V@���@��D@�(�@�;@��@��@\)@~ȴ@~��@~v�@~ff@}��@}�@}V@|�/@{ƨ@{"�@z~�@z^5@z�@y��@y�@x��@xbN@x1'@wK�@w
=@v�R@vV@u��@u�@t��@t1@t1@s��@s�@sC�@r�@s"�@sdZ@s"�@r�H@r��@r~�@rJ@q��@q�^@q�7@qX@q7L@p��@o��@n��@n��@n�y@o�@o�w@o��@o��@o�@pr�@p�u@pA�@o�@p��@q7L@qG�@qhs@r-@r��@s"�@s33@s"�@s"�@sS�@st�@sdZ@s"�@s33@s��@t9X@t�@u�@u��@v$�@u�-@u`B@uO�@u�@u�@uO�@u?}@u?}@u/@u�@u�@u?}@uO�@u`B@v@vff@v�+@vff@v�R@vȴ@vȴ@v�+@v��@v�@v�R@v�@v�@v��@u�@uO�@t�D@s�F@r�@r��@rn�@rM�@rJ@q��@qx�@qX@qG�@q&�@q�@q7L@q7L@q�@q�@q%@q&�@q�@rM�@qX@r�@s�m@s�m@st�@r��@r~�@r-@rJ@q��@qG�@p��@o�@o;d@o�@n�y@n��@o�w@ol�@o|�@n�y@n�R@n��@nȴ@n5?@m�@l�@lI�@l1@lZ@m�@m�@l��@l��@l��@k�
@j~�@i�#@i��@i��@i�7@i&�@h�9@g�w@f�@fE�@e�@e��@ep�@eO�@e`B@ep�@eV@d��@d�D@d1@c�
@cƨ@c��@c��@cS�@c33@co@b��@bn�@b-@bJ@a��@aX@a�@`r�@`  @_�w@_�P@^�+@]��@]@]��@]�@\�D@\�@[��@[ƨ@[dZ@["�@Z�H@Z��@Z~�@Z=q@Z�@Y�@Y�@Y��@Y��@Y�7@YX@Y7L@Y%@Y�@Y&�@Y%@Y%@XĜ@XbN@X  @W�@Wl�@W+@V��@V�y@Vv�@VE�@V@U��@U@Up�@UV@T�@TZ@TI�@T�D@Tj@T1@Sƨ@SdZ@R�\@Qhs@P�`@Pr�@P1'@O�@P �@O�;@O�w@Ol�@O;d@O
=@N��@NV@O
=@O|�@Ol�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�'B�'B�'B�'B�'B�'B�!B�!B�!B�!B�!B�B�!B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�!B�RB�B�NB�`B�TB�sB�yB�B�B��BB��B��BB+B	7BoB�B�B�B�B�B�B�BuB\BJB
=BPB\B\B\BVBDB+BBBBBBBB��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B�B�B�B�yB�yB�B�B�B�fB�mB�ZB�ZB�`B�sB�sB�yB�B�mB�fB�fB�fB�fB�mB�sB�fB�`B�`B�`B�TB�TB�BB�/B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BǮBǮBŢBĜBĜBĜBÖBÖBBBB��B��B��B��B��B��B��B��B��B��B��B��B��BÖBƨBƨBƨBŢBŢBŢBĜBŢBƨBŢBŢBǮBǮBȴBȴBȴBȴBȴBȴBȴBȴBǮBƨBƨBŢBŢBŢBĜBÖB��B��B��B�}B�wB�jB�jB�^B�dB�^B�XB�LB�LB�FB�FB�RB�XB�XB�RB�LB�FB�FB�?B�9B�9B�3B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�oB�oB�hB�hB�hB�bB�bB�\B�\B�\B�VB�VB�PB�PB�DB�DB�DB�=B�=B�=B�=B�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�7B�+B�%B�%B�%B�+B�7B�7B�=B�DB�PB�PB�PB�PB�\B�hB�hB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�!B�'B�3B�9B�?B�FB�LB�LB�XB�RB�XB�XB�^B�^B�XB�RB�FB�?B�9B�9B�?B�?B�9B�9B�?B�FB�FB�FB�LB�LB�RB�XB�dB�jB�wB��BBBǮB��B��BɺBǮBǮBǮBǮBƨBƨBŢBÖBBBBBƨBƨBƨBƨBƨBƨBǮBŢBĜBÖBB��BÖBŢBŢBŢBŢBŢBĜB��B��B��BBBB��B�}B�wB�wB�wB�wB�wB�}B�}B�}B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�wB�qB�jB�jB�jB�jB�dB�dB�dB�^B�^B�^B�^B�dB�dB�dB�jB�dB�jB�jB�qB�wB�wB�wB�wB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BBBBBB��B�}B�wB�qB�qB�qB�wB�wB�wB�wB�wB�wB�qB�wB��BBÖ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�"B��B�B�0B�$B�CB�IB�B�B��B��B��B��B��B�BB>BiBoBiBbBbBVBPBDB+BB	BB+B+B+B%B
B�B�B�B �B�B��B��B��B��B��B��B��B�B�nB�OB�[B�tB�tB�hB�hB�hB�nB�zB�B�B�tB�hB�hB�UB�UB�UB�UB�IB�IB�UB�nB�OB�6B�=B�*B�*B�0B�CB�CB�IB�OB�=B�6B�6B�6B�6B�=B�CB�6B�0B�0B�0B�$B�$B�B��B��B��B��B��B��B��BмBмB��B��B��B��B��B��B��BмBϵBͩB˝B�B�B�sB�mB�mB�mB�gB�gB�`B�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�TB�TB�ZB�gB�yB�yB�yB�sB�sB�sB�mB�sB�yB�sB�sB�B�BǅBǅBǅBǅBǅBǅBǅBǅB�B�yB�yB�sB�sB�sB�mB�gB�ZB�TB�TB�NB�HB�;B�;B�/B�5B�/B�)B�B�B�B�B�#B�)B�)B�#B�B�B�B�B�
B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�wB�wB�wB�wB�wB�wB�qB�kB�kB�kB�kB�dB�dB�dB�^B�^B�XB�XB�XB�XB�LB�FB�FB�AB�AB�:B�:B�:B�4B�4B�.B�.B�.B�(B�(B�"B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B��B��B��B��B��B�	B�	B�B�B�"B�"B�"B�"B�.B�:B�:B�AB�SB�^B�dB�dB�kB�qB�wB�wB�}B�}B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�
B�B�B�B�B�)B�#B�)B�)B�/B�/B�)B�#B�B�B�
B�
B�B�B�
B�
B�B�B�B�B�B�B�#B�)B�5B�;B�HB�ZB�`B�`B�BɒBɒBȋB�B�B�B�B�yB�yB�sB�gB�`B�`B�`B�`B�yB�yB�yB�yB�yB�yB�B�sB�mB�gB�`B�ZB�gB�sB�sB�sB�sB�sB�mB�ZB�TB�ZB�`B�`B�`B�ZB�NB�HB�HB�HB�HB�HB�NB�NB�NB�NB�TB�TB�TB�TB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�TB�NB�NB�NB�HB�BB�;B�;B�;B�;B�5B�5B�5B�/B�/B�/B�/B�5B�5B�5B�;B�5B�;B�;B�BB�HB�HB�HB�HB�NB�TB�TB�ZB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�ZB�`B�`B�`B�`B�`B�TB�NB�HB�BB�BB�BB�HB�HB�HB�HB�HB�HB�BB�HB�ZB�`B�g11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0011569                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904382018110709043820181107090438  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172705  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172705  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090438  IP  PSAL            @ffD��3G�O�                