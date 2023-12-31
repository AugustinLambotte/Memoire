CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:18Z creation; 2018-11-05T17:27:22Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z        X  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        X  D   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   Ld   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     X  N|   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  V�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   _,   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  aD   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   i�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  k�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  t   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   |d   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  ~|   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �    HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �    	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �D   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �t   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �t   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �t   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �t             ,  �tArgo profile    3.1 1.2 19500101000000  20181105172618  20181107090454  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               ,A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�g�	U��1   @�g�U�lx@O0�.Ǔ��C���&�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @y��@�  @�  @�  @���A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A���A�  B   B  B  B  BffB  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%�C'� C*  C,� C/  C1� C4  C6ffC9  C;� C>  C@� CB�fCE� CH  CJ� CM�CO��CR�CT��CW�CY��C\�C^� C`�fCc� Cf  Ch� Ck�Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C��C�@ C�s3C�� C��C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�s3CƳ3C�  C�@ Cʀ C˳3C�  C�L�Cπ C�� C�  C�@ CԌ�C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C��C�� C�  C�L�C� C�� C�  C�@ C��C�� C�  C�@ C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D��D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!FfD"�fD#�fD%fD&@ D'�fD(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5FfD6� D7��D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DEy�DF��DH  DI@ DJy�DK� DMfDN@ DO� DP�fDR  DS@ DT� DU� DW  DX@ DY� DZ�fD\  D]@ D^� D_�fDafDb@ Dc� Dd� DffDg@ Dh� Di� Dj��Dl9�Dmy�Dn� Dp  Dq@ Dry�Ds��Dt��Dv@ Dw� Dx� Dy��D{@ D|�fD}�fDfD�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�<�D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D��D�� D�` D���D���D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�C3D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  Dã3D�C3D�� Dŀ D�#3D��3D�c3D�  DȠ D�@ D�� Dʀ D�  D��3D�c3D�  D͠ D�@ D�� Dσ3D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  Dռ�D�\�D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dހ D��D�� D�c3D�3D� D�@ D�� D� D�  D�� D�` D�  D��D�@ D�� D�3D�&fD���D�  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@@  @y��@�  @�  @�  @���A  A   A0  A@  AP  A`  Aq��A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A���A�  B   B  B  B  BffB  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%�C'� C*  C,� C/  C1� C4  C6ffC9  C;� C>  C@� CB�fCE� CH  CJ� CM�CO��CR�CT��CW�CY��C\�C^� C`�fCc� Cf  Ch� Ck�Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C��C�@ C�s3C�� C��C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�s3CƳ3C�  C�@ Cʀ C˳3C�  C�L�Cπ C�� C�  C�@ CԌ�C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C��C�� C�  C�L�C� C�� C�  C�@ C��C�� C�  C�@ C�� C�� C�  C���C�  D � D  D@ D� D� D  D@ D	� D
� D��D@ D� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!FfD"�fD#�fD%fD&@ D'�fD(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5FfD6� D7��D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DEy�DF��DH  DI@ DJy�DK� DMfDN@ DO� DP�fDR  DS@ DT� DU� DW  DX@ DY� DZ�fD\  D]@ D^� D_�fDafDb@ Dc� Dd� DffDg@ Dh� Di� Dj��Dl9�Dmy�Dn� Dp  Dq@ Dry�Ds��Dt��Dv@ Dw� Dx� Dy��D{@ D|�fD}�fDfD�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�<�D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D��D�� D�` D���D���D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�C3D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  Dã3D�C3D�� Dŀ D�#3D��3D�c3D�  DȠ D�@ D�� Dʀ D�  D��3D�c3D�  D͠ D�@ D�� Dσ3D�  D�� D�` D�3DҠ D�@ D�� DԀ D�  Dռ�D�\�D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dހ D��D�� D�c3D�3D� D�@ D�� D� D�  D�� D�` D�  D��D�@ D�� D�3D�&fD���D�  111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�A%AoA�A�A�A�A�A�A�AoAVAoA
=A%AA��AȴAI�A1A�TAhsA �yA �@�~�@��@��-@�7L@�7L@�7L@�?}@�7L@�&�@��@��@�&�@��@�`B@��7@��7@�`B@�?}@�7L@�G�@�G�@�7L@��@��@�%@�/@�?}@�G�@�G�@�G�@�?}@�7L@�G�@��`@�Ĝ@�Ĝ@�A�@�S�@�ȴ@���@�n�@�$�@�@��T@��#@�@��-@��h@���@���@�b@�ƨ@�l�@�33@��@���@��\@�V@���@���@�@�D@�D@���@��@���@��@�
=@���@�x�@��@�@�z�@�J@�/@��@��@�b@��@ݙ�@��@ۮ@�5?@�?}@��@�hs@ԣ�@�dZ@�ff@�5?@�-@���@��@��@���@θR@���@̬@˥�@�?}@��@Ɵ�@�ff@�$�@ź^@ģ�@ÍP@î@ċD@�r�@�33@��@��@���@�|�@��H@�p�@�1'@���@���@��P@�+@���@��`@��@��
@�ƨ@��w@��
@��@�b@�(�@���@���@��7@�V@��@��@�Q�@�  @��@�"�@�v�@��@��-@���@�j@�b@�ȴ@���@�hs@��@�I�@�(�@��m@�S�@�"�@�ȴ@���@�-@���@��-@��@�O�@�V@��@��@�ƨ@�\)@��y@�~�@�5?@�5?@�G�@��@���@��u@��@�\)@�
=@��!@���@��T@�O�@�7L@��@��@���@��@�z�@��F@��y@���@���@�ff@���@��7@�/@��@��P@�t�@��P@�
=@�=q@�@��7@�&�@�Ĝ@��@�I�@���@��;@��P@��+@���@�p�@�&�@��`@�1@��P@�33@��@��y@�ȴ@�ȴ@��y@���@���@�n�@�ff@�^5@�5?@�@��#@��h@��7@��@��/@��@�1@�b@��m@���@�C�@��@���@��@��@��j@��@�l�@�"�@�@�ȴ@��@�`B@��u@�Z@�9X@��m@��P@�K�@�@�ff@�E�@��-@�&�@���@��`@�Ĝ@���@�z�@��@���@��w@�t�@�"�@���@�v�@�M�@�5?@�@��h@�7L@�%@��`@��/@���@��j@���@���@��D@�A�@� �@�1@|�@K�@�@~��@~@}�-@}�h@}p�@}?}@|��@|��@|9X@|1@{�F@{dZ@{o@z��@z�!@z�!@z��@z��@z=q@y�#@y��@y��@yX@x�@xr�@xQ�@xbN@xQ�@x  @w�P@w;d@v�y@v��@vV@vE�@v5?@v{@v@v{@v@v@v@u�@u�T@u�-@u�-@u`B@t��@t�@t�D@t�@s33@s@r��@r~�@qhs@p��@o�@o\)@o�@n��@n��@m�T@m�@l�j@kdZ@j�\@i�@i��@i�@h1'@g��@g�w@gK�@f�+@e�@e`B@d�@d�@cdZ@c@b��@b~�@a��@`�9@`Q�@_�@_��@_|�@_\)@^�R@]@]@\�D@[�m@Z�\@Y�#@XQ�@X1'@X�u@W�@W�@V5?@V{@V5?@V5?@VE�@V5?@VE�@U@T��@St�@R��@Rn�@Q�^@QX@P1'@OK�@N��@N�+@L��@L�@L��@L��@L�@K33@Ko@I��@H�`@H��@HbN@H �@G�@Gl�@G+@F�@Fȴ@Fv�@E�@D�D@D�@Ct�@B�H@BM�@A��@Ahs@@�9@@Q�@?\)@>�@>��@<�D@;�m@;ƨ@;�@9��@9G�@9%@7��@7
=@6ff@5�@49X@3�@2~�@1x�@0��@/�@.��@.�@.E�@-��@,Z@*�@)�^@)�@&��@%V@$z�@#"�@!�@!&�@ bN@  �@�@?}@"�@�u@b@�w@��@�@�@��@��@Z@I�@9X@9X@(�@�m@ƨ@ƨ@ƨ@ƨ@��@t�@C�@"�@"�@t�@�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�A%AoA�A�A�A�A�A�A�AoAVAoA
=A%AA��AȴAI�A1A�TAhsA �yA �@�~�@��@��-@�7L@�7L@�7L@�?}@�7L@�&�@��@��@�&�@��@�`B@��7@��7@�`B@�?}@�7L@�G�@�G�@�7L@��@��@�%@�/@�?}@�G�@�G�@�G�@�?}@�7L@�G�@��`@�Ĝ@�Ĝ@�A�@�S�@�ȴ@���@�n�@�$�@�@��T@��#@�@��-@��h@���@���@�b@�ƨ@�l�@�33@��@���@��\@�V@���@���@�@�D@�D@���@��@���@��@�
=@���@�x�@��@�@�z�@�J@�/@��@��@�b@��@ݙ�@��@ۮ@�5?@�?}@��@�hs@ԣ�@�dZ@�ff@�5?@�-@���@��@��@���@θR@���@̬@˥�@�?}@��@Ɵ�@�ff@�$�@ź^@ģ�@ÍP@î@ċD@�r�@�33@��@��@���@�|�@��H@�p�@�1'@���@���@��P@�+@���@��`@��@��
@�ƨ@��w@��
@��@�b@�(�@���@���@��7@�V@��@��@�Q�@�  @��@�"�@�v�@��@��-@���@�j@�b@�ȴ@���@�hs@��@�I�@�(�@��m@�S�@�"�@�ȴ@���@�-@���@��-@��@�O�@�V@��@��@�ƨ@�\)@��y@�~�@�5?@�5?@�G�@��@���@��u@��@�\)@�
=@��!@���@��T@�O�@�7L@��@��@���@��@�z�@��F@��y@���@���@�ff@���@��7@�/@��@��P@�t�@��P@�
=@�=q@�@��7@�&�@�Ĝ@��@�I�@���@��;@��P@��+@���@�p�@�&�@��`@�1@��P@�33@��@��y@�ȴ@�ȴ@��y@���@���@�n�@�ff@�^5@�5?@�@��#@��h@��7@��@��/@��@�1@�b@��m@���@�C�@��@���@��@��@��j@��@�l�@�"�@�@�ȴ@��@�`B@��u@�Z@�9X@��m@��P@�K�@�@�ff@�E�@��-@�&�@���@��`@�Ĝ@���@�z�@��@���@��w@�t�@�"�@���@�v�@�M�@�5?@�@��h@�7L@�%@��`@��/@���@��j@���@���@��D@�A�@� �@�1@|�@K�@�@~��@~@}�-@}�h@}p�@}?}@|��@|��@|9X@|1@{�F@{dZ@{o@z��@z�!@z�!@z��@z��@z=q@y�#@y��@y��@yX@x�@xr�@xQ�@xbN@xQ�@x  @w�P@w;d@v�y@v��@vV@vE�@v5?@v{@v@v{@v@v@v@u�@u�T@u�-@u�-@u`B@t��@t�@t�D@t�@s33@s@r��@r~�@qhs@p��@o�@o\)@o�@n��@n��@m�T@m�@l�j@kdZ@j�\@i�@i��@i�@h1'@g��@g�w@gK�@f�+@e�@e`B@d�@d�@cdZ@c@b��@b~�@a��@`�9@`Q�@_�@_��@_|�@_\)@^�R@]@]@\�D@[�m@Z�\@Y�#@XQ�@X1'@X�u@W�@W�@V5?@V{@V5?@V5?@VE�@V5?@VE�@U@T��@St�@R��@Rn�@Q�^@QX@P1'@OK�@N��@N�+@L��@L�@L��@L��@L�@K33@Ko@I��@H�`@H��@HbN@H �@G�@Gl�@G+@F�@Fȴ@Fv�@E�@D�D@D�@Ct�@B�H@BM�@A��@Ahs@@�9@@Q�@?\)@>�@>��@<�D@;�m@;ƨ@;�@9��@9G�@9%@7��@7
=@6ff@5�@49X@3�@2~�@1x�@0��@/�@.��@.�@.E�@-��@,Z@*�@)�^@)�@&��@%V@$z�@#"�@!�@!&�@ bN@  �@�@?}@"�@�u@b@�w@��@�@�@��@��@Z@I�@9X@9X@(�@�m@ƨ@ƨ@ƨ@ƨ@��@t�@C�@"�@"�@t�@�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB7LB7LB7LB7LB7LB7LB7LB6FB6FB7LB6FB6FB5?B5?B5?B49B33B49B5?B5?B49B6FB5?B33B33B/B/B/B0!B0!B1'B1'B1'B1'B1'B2-B1'B49B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B7LB:^B;dB>wB>wB=qB>wB>wB>wB>wB>wB=qB?}B>wB>wB?}B?}B?}BA�BB�BC�BC�BD�BD�BE�BE�BE�BE�BE�BF�BE�BF�BE�BF�BD�BE�BH�BH�BH�BG�BF�BE�BK�BK�BN�BL�BJ�BO�BR�BW
BVBS�BVBXBaHB`BB`BBcTBffBffBgmBgmBhsBjBiyBhsBgmBffBhsBk�BjBl�Bm�Bm�BjBcTBaHBcTBcTBbNBaHB_;B_;BaHBgmBffBdZBbNBaHB^5B]/B[#BXBT�BR�BR�BQ�BO�BK�BH�BG�BG�BG�BG�BG�BG�BG�BG�BG�BE�BC�BA�BA�B@�B@�B?}B>wB=qB;dB9XB8RB6FB49B2-B.B,B)�B(�B'�B&�B&�B%�B%�B$�B#�B"�B"�B!�B!�B �B �B�B�B�B�B�B�B�B�B�B�B�B{BuBhBbBbB\BVBJBJBJBDBDB
=B	7B1B%B%BBBBBB��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�yB�sB�fB�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�NB�NB�HB�HB�BB�;B�;B�5B�5B�5B�5B�5B�5B�/B�)B�#B�B�B�B��B��B��B��B��B��B��B��B��B��B��BɺBȴBƨBƨBĜBÖBBB��B��B��B��B�}B�}B�}B�wB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�qB�wB�wB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BBBBBBBBBBBBÖBBBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBBBBBBBBBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBBBBBB��B��B��B��B�}B�wB�qB�qB�qB�qB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�XB�RB�LB�LB�FB�FB�?B�?B�9B�9B�9B�9B�3B�3B�3B�3B�-B�-B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B-B-B-B-B-B-B-B+�B+�B-B+�B+�B*�B*�B*�B)�B(�B)�B*�B*�B)�B+�B*�B(�B(�B$�B$�B$�B%�B%�B&�B&�B&�B&�B&�B'�B&�B)�B*�B*�B*�B*�B*�B*�B*�B*�B*�B*�B-B0B1B40B40B3*B40B40B40B40B40B3*B56B40B40B56B56B56B7BB8GB9NB9NB:TB:TB;ZB;ZB;ZB;ZB;ZB<`B;ZB<`B;ZB<`B:UB;[B>lB>lB>lB=fB<aB;[BABABD�BB�B@yBE�BH�BL�BK�BI�BK�BM�BW BU�BU�BYB\B\B]%B]%B^+B`7B_1B^,B]&B\B^,Ba=B`8BbDBcIBcJB`8BYBWBYBYBXBWBT�BT�BWB]'B\ BZBXBWBS�BR�BP�BM�BJ�BH�BH�BG�BE�BA�B>qB=kB=kB=kB=kB=kB=kB=kB=kB=kB;_B9SB7GB7GB6AB6AB5;B45B3/B1"B/B.B,B)�B'�B#�B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BByBsBmBgBbB\BUBIBCBCB
=B	7B*B%B%BBBBBBBB  B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�|B�vB�vB�jB�^B�QB�KB�FB�@B�:B�-B�'B�'B�!B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȻBƮBũBÝBBB��B��B��B�~B�sB�sB�gB�aB�ZB�ZB�TB�TB�NB�NB�HB�HB�HB�BB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�6B�6B�6B�6B�6B�6B�6B�6B�6B�6B�6B�6B�6B�6B�=B�CB�CB�IB�OB�UB�UB�UB�UB�UB�UB�UB�UB�UB�UB�UB�UB�UB�UB�[B�[B�[B�[B�[B�[B�[B�[B�[B�[B�[B�bB�[B�[B�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�[B�[B�[B�[B�[B�[B�[B�[B�[B�bB�bB�bB�bB�bB�bB�bB�bB�cB�cB�cB�cB�cB�\B�\B�\B�\B�\B�\B�VB�VB�VB�PB�JB�DB�>B�>B�>B�>B�8B�8B�2B�2B�2B�2B�2B�2B�2B�2B�&B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0099889                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904552018110709045520181107090455  IF  ARFMCODA024c                                                                20181105172618                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172722  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172722  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090455  IP  PSAL            @��D�  G�O�                