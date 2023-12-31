CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:20Z creation; 2018-10-26T19:00:47Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181026190020  20181119104029  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               6A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @���hs1   @���hs@SLf�Cf@��F��8   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @���@���@�33A��A  A   A0  A@  AP  A^ffAnffA�  A�  A�  A���A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A���B   B  B  B  BffB  B  B  B   B$  B(  B,ffB0  B4  B8ffB<  B@  BD  BG��BL  BP  BS��BX  B\  B_��Bd  Bh  Bl  Bp  BtffBx  B|ffB�33B�  B�  B�  B���B�  B�  B�  B�  B�33B�  B�  B���B�  B�  B�  B���B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	��C�C� C  CffC  C��C  C� C   C"� C$�fC'� C*  C,ffC/  C1��C4  C6� C9  C;� C=�fC@� CC  CE� CH  CJffCM  CO� CQ�fCT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck�Cm� Cp  Cr��Cu  Cw� Cz�C|� C~�fC�� C��C�@ C�� C���C�  C�L�C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C��3C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C߳3C��3C�@ C��C�� C�  C�@ C� C�� C�  C�@ C�s3C�� C��C�@ C�s3C�� C��C�@ C�� C�� C�  C�� C�  D �fD  D@ D� D� DfD@ D	y�D
��D  DFfD� D��D  D@ D� D��D  D@ D� D� D  D@ Dy�D� D   D!FfD"� D#� D%  D&@ D'� D(��D*  D+@ D,� D-� D/  D0FfD1� D2� D4  D5@ D6� D7��D9  D:@ D;�fD<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DL��DN@ DO�fDP� DR  DS@ DT� DU�fDWfDX@ DY� DZ��D\  D]@ D^� D_� D`��Db@ Dc� Dd� Df  Dg@ Dh�fDi�fDkfDlFfDm�fDn�fDp  DqFfDr� Ds��Du  DvFfDw� Dx� Dz  D{@ D|�fD}� D  D�  D�� D�` D�  D���D�<�D���D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�\�D�  D��3D�@ D���D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D���D�@ D�� D�� D�  D�� D�` D�3Dà D�@ D���D�|�D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�3D͠ D�@ D�� Dσ3D�#3D�� D�` D�  DҠ D�@ D���DԀ D��Dռ�D�` D�  Dנ D�@ D��3Dـ D�  Dڼ�D�\�D�  Dܠ D�@ D��3Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D��D�` D�  D� D�@ D��3D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D� D�@ D�� D� D�  D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @���@���@�33A��A  A   A0  A@  AP  A^ffAnffA�  A�  A�  A���A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A���B   B  B  B  BffB  B  B  B   B$  B(  B,ffB0  B4  B8ffB<  B@  BD  BG��BL  BP  BS��BX  B\  B_��Bd  Bh  Bl  Bp  BtffBx  B|ffB�33B�  B�  B�  B���B�  B�  B�  B�  B�33B�  B�  B���B�  B�  B�  B���B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	��C�C� C  CffC  C��C  C� C   C"� C$�fC'� C*  C,ffC/  C1��C4  C6� C9  C;� C=�fC@� CC  CE� CH  CJffCM  CO� CQ�fCT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck�Cm� Cp  Cr��Cu  Cw� Cz�C|� C~�fC�� C��C�@ C�� C���C�  C�L�C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C�� C��3C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C߳3C��3C�@ C��C�� C�  C�@ C� C�� C�  C�@ C�s3C�� C��C�@ C�s3C�� C��C�@ C�� C�� C�  C�� C�  D �fD  D@ D� D� DfD@ D	y�D
��D  DFfD� D��D  D@ D� D��D  D@ D� D� D  D@ Dy�D� D   D!FfD"� D#� D%  D&@ D'� D(��D*  D+@ D,� D-� D/  D0FfD1� D2� D4  D5@ D6� D7��D9  D:@ D;�fD<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DL��DN@ DO�fDP� DR  DS@ DT� DU�fDWfDX@ DY� DZ��D\  D]@ D^� D_� D`��Db@ Dc� Dd� Df  Dg@ Dh�fDi�fDkfDlFfDm�fDn�fDp  DqFfDr� Ds��Du  DvFfDw� Dx� Dz  D{@ D|�fD}� D  D�  D�� D�` D�  D���D�<�D���D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�\�D�  D��3D�@ D���D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D���D�@ D�� D�� D�  D�� D�` D�3Dà D�@ D���D�|�D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�3D͠ D�@ D�� Dσ3D�#3D�� D�` D�  DҠ D�@ D���DԀ D��Dռ�D�` D�  Dנ D�@ D��3Dـ D�  Dڼ�D�\�D�  Dܠ D�@ D��3Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D��D�` D�  D� D�@ D��3D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D� D�@ D�� D� D�  D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�@���@���@���@���@���@�@�@�@�@�@�@�J@�@�@�@�@�@�@�@�@���@��@��@���@���@�@�@�@���@�@�J@�@���@�@��@��@��@��@��@���@��@��@���@���@��@��T@��#@��#@���@���@��@��@��@��T@��T@��@��@��@��@��@�@�@���@��@��@�@�J@�J@�@�J@�J@�{@�{@�{@�J@�@�J@��@��@��@�{@��@�{@�{@��@��@�$�@�$�@�$�@�$�@�$�@�$�@�$�@�$�@�-@�$�@�-@�$�@�-@�-@�-@�-@�-@�$�@�$�@�-@�$�@�$�@�5?@�5?@�5?@�5?@�5?@�=q@�E�@�=q@�M�@�V@�V@�=q@�E�@�-@�M�@�^5@�V@�^5@�^5@�ff@�ff@�ff@�n�@�n�@�v�@�v�@�^5@�ff@�^5@�V@�V@�E�@�^5@�^5@�^5@�M�@�ff@�M�@��T@�O�@���@�`B@�7L@�&�@��9@�Z@�j@�Q�@~ȴ@{�F@{��@{dZ@z�!@z��@{o@{S�@{�@{�
@{�
@{�m@{��@{�m@{�m@{�
@{��@{dZ@z�\@x��@x �@w+@w�@v��@u�@s�F@so@q��@pĜ@p1'@o
=@n�y@nV@nv�@n{@mO�@l�j@l�@l�@l�j@l�@j�H@j�H@j�!@k@j-@h�`@h��@hbN@h��@ix�@iG�@h��@f��@eV@d�@c�@a��@`1'@^�+@^ff@\�@\I�@[�@Z��@Z^5@Yx�@X  @R��@M��@L��@L�@F��@C"�@B�!@?�w@?�w@@b@@��@A�^@B�!@Co@A��@@r�@:��@5V@/��@';d@!�^@{@�@+@@��@
�H@v�@�@~�@�@�@=q@ r�?���?��^?��?�ȴ?��?�7?� �?��;?߾w?��?㕁?���?�A�?�=q?�&�?�~�?�V?��j?�C�?�;d?��9?�-?��m?�S�?q��?_|�?Q��?C��?+C�?;d?��?�D?�9?��?o>�j>���>�=q>�%>�X>���>��/>���>�t�>��>]/>I�^>0 �>bN>o=�S�=��=\=���=��-=���=Y�=C�=o<�`B<�j<o;ě���o��t����
�+�,1�P�`��7L���㽩�罸Q�\�����m���   �o�	7L�
=q�C��hs�������w�&�y�/��333�6E��<j�?|�B�\�E�˾H�9�N��T���Xb�["Ѿ_;d�fff�o���s�F�vȴ�w�پ~�۾��\��$ݾ��9��C���I���녾�zᾗ
=���������"Ѿ�/���w��G���MӾ�MӾ�MӾ�`B���xվ�~���1�����33����ȴ������^5���m���۾�J�\�Õ��š˾�1'�ȴ9���;��;���;��n��Ձ�և+��b�ٙ���(���/�߾w��A��������T��r���xվ������{����������پ��m�   �Mӿ�7�J��
�o������T�$ݿ$ݿ$ݿ$ݿff�r��r��r���9�	7L�	�^�
=q��ƨ�1��D��D��D��ͿO߿�h��h��h����������� ſbN�hs�n���F�z�?}���E��Kǿ�ٿ�ٿ���#����dZ��m�������/��-�5?�vɿ5?��ۿ|�|�|�|�|�   �   � A�� �� ��!G��!���"J�!���!�7�"�\�#o�#S��#���#���#���#S��#���$��$�/�%��%�˿'+�'''(1'�(1'�(1'�(�9�(�9�)7L�)7L�)7L�)xտ)xտ)�^�)�^�*=q�*~��+�+��+ƨ�-V�-�h�-�h�-O߿-��-��-��/��/\)�/���/�;�/�;�0bN�0�׿0�`�0�`�0�`�1&�1녿1녿2-�2-�2-�2-�2n��2n�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @���@�@���@���@���@���@���@�@�@�@�@�@�@�J@�@�@�@�@�@�@�@�@���@��@��@���@���@�@�@�@���@�@�J@�@���@�@��@��@��@��@��@���@��@��@���@���@��@��T@��#@��#@���@���@��@��@��@��T@��T@��@��@��@��@��@�@�@���@��@��@�@�J@�J@�@�J@�J@�{@�{@�{@�J@�@�J@��@��@��@�{@��@�{@�{@��@��@�$�@�$�@�$�@�$�@�$�@�$�@�$�@�$�@�-@�$�@�-@�$�@�-@�-@�-@�-@�-@�$�@�$�@�-@�$�@�$�@�5?@�5?@�5?@�5?@�5?@�=q@�E�@�=q@�M�@�V@�V@�=q@�E�@�-@�M�@�^5@�V@�^5@�^5@�ff@�ff@�ff@�n�@�n�@�v�@�v�@�^5@�ff@�^5@�V@�V@�E�@�^5@�^5@�^5@�M�@�ff@�M�@��T@�O�@���@�`B@�7L@�&�@��9@�Z@�j@�Q�@~ȴ@{�F@{��@{dZ@z�!@z��@{o@{S�@{�@{�
@{�
@{�m@{��@{�m@{�m@{�
@{��@{dZ@z�\@x��@x �@w+@w�@v��@u�@s�F@so@q��@pĜ@p1'@o
=@n�y@nV@nv�@n{@mO�@l�j@l�@l�@l�j@l�@j�H@j�H@j�!@k@j-@h�`@h��@hbN@h��@ix�@iG�@h��@f��@eV@d�@c�@a��@`1'@^�+@^ff@\�@\I�@[�@Z��@Z^5@Yx�@X  @R��@M��@L��@L�@F��@C"�@B�!@?�w@?�w@@b@@��@A�^@B�!@Co@A��@@r�@:��@5V@/��@';d@!�^@{@�@+@@��@
�H@v�@�@~�@�@�@=q@ r�?���?��^?��?�ȴ?��?�7?� �?��;?߾w?��?㕁?���?�A�?�=q?�&�?�~�?�V?��j?�C�?�;d?��9?�-?��m?�S�?q��?_|�?Q��?C��?+C�?;d?��?�D?�9?��?o>�j>���>�=q>�%>�X>���>��/>���>�t�>��>]/>I�^>0 �>bN>o=�S�=��=\=���=��-=���=Y�=C�=o<�`B<�j<o;ě���o��t����
�+�,1�P�`��7L���㽩�罸Q�\�����m���   �o�	7L�
=q�C��hs�������w�&�y�/��333�6E��<j�?|�B�\�E�˾H�9�N��T���Xb�["Ѿ_;d�fff�o���s�F�vȴ�w�پ~�۾��\��$ݾ��9��C���I���녾�zᾗ
=���������"Ѿ�/���w��G���MӾ�MӾ�MӾ�`B���xվ�~���1�����33����ȴ������^5���m���۾�J�\�Õ��š˾�1'�ȴ9���;��;���;��n��Ձ�և+��b�ٙ���(���/�߾w��A��������T��r���xվ������{����������پ��m�   �Mӿ�7�J��
�o������T�$ݿ$ݿ$ݿ$ݿff�r��r��r���9�	7L�	�^�
=q��ƨ�1��D��D��D��ͿO߿�h��h��h����������� ſbN�hs�n���F�z�?}���E��Kǿ�ٿ�ٿ���#����dZ��m�������/��-�5?�vɿ5?��ۿ|�|�|�|�|�   �   � A�� �� ��!G��!���"J�!���!�7�"�\�#o�#S��#���#���#���#S��#���$��$�/�%��%�˿'+�'''(1'�(1'�(1'�(�9�(�9�)7L�)7L�)7L�)xտ)xտ)�^�)�^�*=q�*~��+�+��+ƨ�-V�-�h�-�h�-O߿-��-��-��/��/\)�/���/�;�/�;�0bN�0�׿0�`�0�`�0�`�1&�1녿1녿2-�2-�2-�2-�2n��2n�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bl�Bm�Bm�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bm�Bl�Bm�Bl�Bm�Bl�Bl�Bm�Bm�Bm�Bl�Bm�Bl�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bl�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bk�Bk�BjBjBiyBiyBiyBhsBgmBffBdZBcTBbNBbNBbNBcTBffBgmBhsBhsBhsBhsBhsBhsBhsBgmBgmBgmBffBe`BdZBcTBcTBbNBaHBaHBaHB`BB_;B_;B^5B^5B]/B]/B]/B\)B[#B[#B[#B[#B[#BZBZBZBZBYBXBXBW
BXBXBXBVBS�BQ�BO�BO�BN�BM�BK�BJ�BI�BH�BG�BE�BD�BA�B>wB:^B5?B33B0!B-B(�B'�B'�B+B-B0!B2-B33B5?B33B/B'�B�B�BVB
=BBBB  B��B��B��B��B  B  BB%BBB  B��B��B��B��B��B��B��BBBBB��B��B�B�sB�NB�/B�
B��B��B��B��BȴBŢBB�}B�jB�dB�^B�XB�XB�XB�RB�FB�?B�9B�9B�3B�3B�-B�-B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bl�Bm�Bm�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bm�Bl�Bm�Bl�Bm�Bl�Bl�Bm�Bm�Bm�Bl�Bm�Bl�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bl�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bl�Bm�Bm�Bl�Bl�Bk�Bk�BjBjBiyBiyBiyBhsBgmBffBdZBcTBbNBbNBbNBcTBffBgmBhsBhsBhsBhsBhsBhsBhsBgmBgmBgmBffBe`BdZBcTBcTBbNBaHBaHBaHB`BB_;B_;B^5B^5B]/B]/B]/B\)B[#B[#B[#B[#B[#BZBZBZBZBYBXBXBW
BXBXBXBVBS�BQ�BO�BO�BN�BM�BK�BJ�BI�BH�BG�BE�BD�BA�B>wB:^B5?B33B0!B-B(�B'�B'�B+B-B0!B2-B33B5?B33B/B'�B�B�BVB
=BBBB  B��B��B��B��B  B  BB%BBB  B��B��B��B��B��B��B��BBBBB��B��B�B�sB�NB�/B�
B��B��B��B��BȴBŢBB�}B�jB�dB�^B�XB�XB�XB�RB�FB�?B�9B�9B�3B�3B�-B�-B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040292018111910402920181119104029  IF  ARFMCODA024c                                                                20181026190020                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190047  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181026190047  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104030  IP  PSAL            @ffD��3G�O�                