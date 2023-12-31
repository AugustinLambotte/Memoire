CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T18:20:46Z creation; 2018-11-05T18:21:33Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105182046  20181107125410  6902586 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               IA   IF                                  2C  D   NOVA                            SN145                           n/a                             865 @�ɥ:��1   @�ɥOՅ�@N�h�  �C��*�p1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @���@�  @�  @���A  A   A0  A@  AP  A`  Ap  A�  A���A�  A�  A�33A�33A�  A���A���A���A���A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,ffB0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BXffB\  B`  Bd  Bh  Bk��Bp  Bt  Bx  B|  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B���B�  B�33B�33B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B���B�  B�  C  C� C  C	� C  CffC�fC� C�C� C  C� C   C"ffC$�fC'� C*  C,� C/  C1� C4  C6� C8�fC;� C>�C@� CC  CE� CH  CJffCM  CO��CR  CTffCW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|��C  C�� C��C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�s3C��3C�  C�@ C�s3C��3C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�33C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C���C��C�@ Cπ C���C��C�@ C�s3C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�33C�s3C��3C�  C��fC�� C�  D � D  D9�D� D�fD  D@ D	� D
� D  D@ D� D� DfD@ D� D� D  D@ D� D� DfD@ D� D� D   D!@ D"� D#�fD%fD&FfD'� D(��D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA� DCfDDFfDE�fDF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DSFfDT�fDU� DW  DX@ DY�fDZ�fD\  D]@ D^� D_� DafDbFfDc�fDd� Df  Dg@ Dh� Di� Dk  DlFfDm�fDn�fDp  Dq@ Dry�Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D���D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D�� D��3D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D���D�� D�C3D��3D�� D�  D��3D�c3D�  D�� D�@ D��3D�� D�  D��3D�` D�  Dà D�@ D�� Dŀ D�  D�� D�c3D�  DȠ D�@ D�� Dʀ D�  D�� D�c3D�  D͠ D�@ D�� Dπ D�  D�� D�\�D�  DҠ D�@ D�� DԀ D�  Dռ�D�\�D�  Dנ D�@ D�� Dـ D�#3D�� D�` D�  Dܜ�D�@ D�� Dހ D�  D�� D�` D�  D�3D�C3D�� D�|�D�  D�� D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�3D�  D��D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D��3D�@ 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @���@�  @�  @���A  A   A0  A@  AP  A`  Ap  A�  A���A�  A�  A�33A�33A�  A���A���A���A���A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,ffB0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BXffB\  B`  Bd  Bh  Bk��Bp  Bt  Bx  B|  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B���B�  B�33B�33B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B���B�  B�  C  C� C  C	� C  CffC�fC� C�C� C  C� C   C"ffC$�fC'� C*  C,� C/  C1� C4  C6� C8�fC;� C>�C@� CC  CE� CH  CJffCM  CO��CR  CTffCW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|��C  C�� C��C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�s3C��3C�  C�@ C�s3C��3C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�33C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C���C��C�@ Cπ C���C��C�@ C�s3C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�33C�s3C��3C�  C��fC�� C�  D � D  D9�D� D�fD  D@ D	� D
� D  D@ D� D� DfD@ D� D� D  D@ D� D� DfD@ D� D� D   D!@ D"� D#�fD%fD&FfD'� D(��D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA� DCfDDFfDE�fDF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DSFfDT�fDU� DW  DX@ DY�fDZ�fD\  D]@ D^� D_� DafDbFfDc�fDd� Df  Dg@ Dh� Di� Dk  DlFfDm�fDn�fDp  Dq@ Dry�Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D���D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D�� D��3D�  D���D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D���D�� D�C3D��3D�� D�  D��3D�c3D�  D�� D�@ D��3D�� D�  D��3D�` D�  Dà D�@ D�� Dŀ D�  D�� D�c3D�  DȠ D�@ D�� Dʀ D�  D�� D�c3D�  D͠ D�@ D�� Dπ D�  D�� D�\�D�  DҠ D�@ D�� DԀ D�  Dռ�D�\�D�  Dנ D�@ D�� Dـ D�#3D�� D�` D�  Dܜ�D�@ D�� Dހ D�  D�� D�` D�  D�3D�C3D�� D�|�D�  D�� D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�3D�  D��D�\�D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D��3D�@ 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A��A�uA�uA�uA�DA�\A�\A�+A�A~�A~�A~�A~�A~�A~�A~�A~�A~�Az�Av�Av�AffAQ�A��A�^A�!A	�A�+@���@�9X@ҟ�@�-@���@�Q�@��T@��@��P@���@�C�@�?}@�ff@��m@��@��H@�t�@�V@�l�@���@�j@�n�@���@��@�J@��-@���@��R@��@��/@��@��D@�Q�@��F@�j@��@��+@��@��@�X@���@�-@��\@�+@�\)@��@�  @�bN@� �@��
@�C�@��H@���@�J@�x�@��@���@���@���@���@�+@� �@��+@�1'@��u@�%@�&�@��@���@��9@�bN@�Z@��/@���@�z�@�A�@��@���@��j@��@�bN@��@���@�^5@�J@�@��@��D@�\)@��@�5?@��@��@��h@��@�Z@�1@���@��P@�o@�n�@��T@�J@�?}@���@��@�dZ@�5?@��;@��@�&�@�1'@�j@�b@~�y@}@}@�w@�@l�@~�@~ȴ@~��@~��@~ȴ@~��@}�T@}`B@|�@|�/@|�@|�j@|I�@{ƨ@{��@{o@y��@y�7@x�@w��@v��@vV@u��@u�h@u/@t9X@s�F@s�F@t(�@t�@s�F@s33@r��@r��@r�\@rM�@q��@qx�@qX@q7L@p��@p�@p1'@p �@pQ�@p��@qX@qhs@q�@pr�@p1'@p1'@o�@ol�@n�+@nff@n�+@n�+@n�+@n5?@m�@mO�@lI�@lz�@m�@n{@n{@n$�@n@m@m�@l��@l�j@l��@lj@l(�@j�H@jJ@jM�@i�#@i�^@hA�@hb@h  @g�@g�@hQ�@j�\@j��@j��@j�!@jM�@i�@i�7@iX@i%@gl�@f$�@fV@g
=@g�P@g\)@fȴ@g+@g\)@g;d@g
=@f�@f�R@f��@fV@f{@f@f@e�@e��@e�-@e�@e�@e@e��@e�T@e�h@e�@ep�@e`B@e�@e�T@e?}@ep�@e@e?}@e�@eO�@eO�@e/@dI�@cdZ@b�H@b�@b��@b��@b~�@cC�@b��@b�\@b��@c33@b�@b~�@b^5@a��@bJ@b-@b^5@b~�@bn�@b�!@a�#@`�`@`�9@a%@`Q�@^�+@^V@_l�@_�@_
=@^��@\��@\��@]��@\Z@[�
@[C�@[C�@[33@Z�@Z��@Z��@Z��@Z=q@Y��@Y��@Yx�@Yhs@Yx�@Yhs@Yx�@Y�#@Z�!@[��@[�
@[�@[S�@Z�@Z�H@["�@[��@\�/@\�@\��@]V@]�@]�h@]/@\�@\�@\�j@\�D@\(�@[�m@[��@[��@[�m@\j@\j@\9X@[ƨ@Z�H@Z~�@ZM�@Z-@ZJ@Z�\@\Z@\z�@\�D@\9X@Z�!@Zn�@Z�\@Z��@Z�!@Z^5@Y�^@YG�@Y7L@Y��@Y��@Yhs@Y&�@X�u@X  @W��@W�P@W;d@W�P@W��@W��@Wl�@W��@VE�@U��@Vff@V�y@W+@Wl�@W�w@W�@X  @X �@Z�@Y��@Yx�@Yx�@Z�\@[C�@[��@\��@^��@_l�@_��@_�@_�@_�;@_�P@^v�@^��@_�w@a�@c��@cƨ@g�w@i%@iX@ix�@i��@i��@hbN@gK�@f��@f�@f5?@f5?@fE�@fE�@fff@f�y@g+@g�w@h1'@h�@h �@g�@g�;@g�@gK�@f��@fff@f{@e@e��@e`B@e/@d��@d�D@dI�@cdZ@bJ@`Ĝ@`r�@`A�@_�@_�w@_�P@_K�@_l�@_\)@^�y@^�R@^$�@]�@\�@\�j@\z�@\9X@[��@[�@[S�@[@Z~�@Z-@Z�@Y�@Y�^@YG�@X�u@XbN@W�@W��@V�y@VV@V{@U@U?}@UV@T(�@S��@SdZ@R��@R^5@RJ@Q�#@Q��@Q�^@Q�7@QX@P�9@PA�@P  @Ol�@N�R@N��@Nv�@N@M�@M�@L�@L��@Kƨ@K�@K33@Ko@J�@J��@JM�@I�#@Ix�@I�@H�`@I�@I&�@H��@H�9@H��@H��@H�u@H �@G�P@G�@F��@F��@F�R@Fv�@E��@E?}@D��@D��@Dz�@DZ@D9X@C��@Cƨ@Ct�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A��A�uA�uA�uA�DA�\A�\A�+A�A~�A~�A~�A~�A~�A~�A~�A~�A~�Az�Av�Av�AffAQ�A��A�^A�!A	�A�+@���@�9X@ҟ�@�-@���@�Q�@��T@��@��P@���@�C�@�?}@�ff@��m@��@��H@�t�@�V@�l�@���@�j@�n�@���@��@�J@��-@���@��R@��@��/@��@��D@�Q�@��F@�j@��@��+@��@��@�X@���@�-@��\@�+@�\)@��@�  @�bN@� �@��
@�C�@��H@���@�J@�x�@��@���@���@���@���@�+@� �@��+@�1'@��u@�%@�&�@��@���@��9@�bN@�Z@��/@���@�z�@�A�@��@���@��j@��@�bN@��@���@�^5@�J@�@��@��D@�\)@��@�5?@��@��@��h@��@�Z@�1@���@��P@�o@�n�@��T@�J@�?}@���@��@�dZ@�5?@��;@��@�&�@�1'@�j@�b@~�y@}@}@�w@�@l�@~�@~ȴ@~��@~��@~ȴ@~��@}�T@}`B@|�@|�/@|�@|�j@|I�@{ƨ@{��@{o@y��@y�7@x�@w��@v��@vV@u��@u�h@u/@t9X@s�F@s�F@t(�@t�@s�F@s33@r��@r��@r�\@rM�@q��@qx�@qX@q7L@p��@p�@p1'@p �@pQ�@p��@qX@qhs@q�@pr�@p1'@p1'@o�@ol�@n�+@nff@n�+@n�+@n�+@n5?@m�@mO�@lI�@lz�@m�@n{@n{@n$�@n@m@m�@l��@l�j@l��@lj@l(�@j�H@jJ@jM�@i�#@i�^@hA�@hb@h  @g�@g�@hQ�@j�\@j��@j��@j�!@jM�@i�@i�7@iX@i%@gl�@f$�@fV@g
=@g�P@g\)@fȴ@g+@g\)@g;d@g
=@f�@f�R@f��@fV@f{@f@f@e�@e��@e�-@e�@e�@e@e��@e�T@e�h@e�@ep�@e`B@e�@e�T@e?}@ep�@e@e?}@e�@eO�@eO�@e/@dI�@cdZ@b�H@b�@b��@b��@b~�@cC�@b��@b�\@b��@c33@b�@b~�@b^5@a��@bJ@b-@b^5@b~�@bn�@b�!@a�#@`�`@`�9@a%@`Q�@^�+@^V@_l�@_�@_
=@^��@\��@\��@]��@\Z@[�
@[C�@[C�@[33@Z�@Z��@Z��@Z��@Z=q@Y��@Y��@Yx�@Yhs@Yx�@Yhs@Yx�@Y�#@Z�!@[��@[�
@[�@[S�@Z�@Z�H@["�@[��@\�/@\�@\��@]V@]�@]�h@]/@\�@\�@\�j@\�D@\(�@[�m@[��@[��@[�m@\j@\j@\9X@[ƨ@Z�H@Z~�@ZM�@Z-@ZJ@Z�\@\Z@\z�@\�D@\9X@Z�!@Zn�@Z�\@Z��@Z�!@Z^5@Y�^@YG�@Y7L@Y��@Y��@Yhs@Y&�@X�u@X  @W��@W�P@W;d@W�P@W��@W��@Wl�@W��@VE�@U��@Vff@V�y@W+@Wl�@W�w@W�@X  @X �@Z�@Y��@Yx�@Yx�@Z�\@[C�@[��@\��@^��@_l�@_��@_�@_�@_�;@_�P@^v�@^��@_�w@a�@c��@cƨ@g�w@i%@iX@ix�@i��@i��@hbN@gK�@f��@f�@f5?@f5?@fE�@fE�@fff@f�y@g+@g�w@h1'@h�@h �@g�@g�;@g�@gK�@f��@fff@f{@e@e��@e`B@e/@d��@d�D@dI�@cdZ@bJ@`Ĝ@`r�@`A�@_�@_�w@_�P@_K�@_l�@_\)@^�y@^�R@^$�@]�@\�@\�j@\z�@\9X@[��@[�@[S�@[@Z~�@Z-@Z�@Y�@Y�^@YG�@X�u@XbN@W�@W��@V�y@VV@V{@U@U?}@UV@T(�@S��@SdZ@R��@R^5@RJ@Q�#@Q��@Q�^@Q�7@QX@P�9@PA�@P  @Ol�@N�R@N��@Nv�@N@M�@M�@L�@L��@Kƨ@K�@K33@Ko@J�@J��@JM�@I�#@Ix�@I�@H�`@I�@I&�@H��@H�9@H��@H��@H�u@H �@G�P@G�@F��@F��@F�R@Fv�@E��@E?}@D��@D��@Dz�@DZ@D9X@C��@Cƨ@Ct�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB\B\B\BVBVBVBVBVBVBVBPBPBJBJBDB
=B1B%BB  B
��B
�B
�sB
�NBBBPB�BhB.BC�B[#BffBl�Bu�Bv�Bu�Bq�Bv�By�Bz�B{�B� B�B�%B�%B�1B�7B�JB�VB�DB�bB�oB�hB�hB�hB�hB�oB�oB�uB�uB��B��B��B��B��B��B��B��B�B�B�-B�-B�9B�9B�?B�9B�9B�-B�-B�'B�!B�B�B�B�B�B�3B�FB�XB��B��B��B��B��B��B��B��B��B��B��B�B�B�
B�B�#B�#B�#B�B�B�B�B�B��B��B��B��B��BɺBȴBǮB��B��BȴBȴBǮBƨBŢBÖBĜBƨBĜBB��B�qB�RB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�hB�oB�hB�bB�\B�VB�VB�VB�\B�bB�{B��B��B�{B�{B�uB�oB�oB�bB�VB�JB�PB�VB�\B�\B�VB�\B�\B�\B�\B�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�JB�PB�PB�VB�PB�PB�PB�PB�PB�\B�VB�PB�VB�\B�VB�\B�VB�VB�PB�JB�=B�7B�7B�7B�7B�7B�=B�7B�7B�=B�=B�=B�7B�7B�1B�1B�7B�7B�7B�7B�=B�1B�%B�%B�%B�B�B�B�B�B�B�B}�B}�B~�B|�B{�Bz�Bz�By�By�By�By�Bx�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bv�Bw�By�B{�B|�B{�B{�Bz�Bz�Bz�B|�B~�B� B~�B� B�B�B� B~�B~�B~�B~�B}�B}�B|�B|�B}�B~�B~�B~�B}�B{�Bz�Bz�By�By�B{�B� B� B� B~�B{�Bz�B{�B{�B{�Bz�By�Bx�Bx�By�By�Bx�Bx�Bv�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bq�Bp�Br�Bs�Bs�Bt�Bu�Bv�Bv�Bw�B{�Bz�Bz�Bz�B}�B� B�B�%B�=B�PB�PB�VB�VB�VB�PB�JB�PB�hB��B��B��B��B�B�B�B�'B�'B�B�B�B�B�B�B�B�B�B�!B�-B�9B�FB�RB�LB�RB�RB�XB�XB�RB�RB�RB�RB�RB�RB�RB�RB�XB�RB�FB�3B�-B�-B�-B�-B�3B�3B�3B�9B�9B�9B�9B�3B�-B�-B�'B�'B�'B�'B�'B�'B�'B�'B�-B�3B�3B�3B�3B�-B�-B�'B�'B�'B�-B�-B�-B�-B�'B�'B�!B�!B�!B�!B�!B�'B�-B�-B�-B�-B�-B�'B�'B�'B�!B�'B�'B�!B�!B�'B�'B�'B�'B�'B�'B�'B�-B�-B�-B�-B�-B�3B�9B�?B�FB�FB�LB�LB�LB�LB�FB�FB�FB�LB�LB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�RB�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B\B\B\BVBVBVBVBVBVBVBPBPBJBJBDB
=B1B%BB  B
��B
�B
�sB
�NBBBPB�BhB.BC�B[#BffBl�Bu�Bv�Bu�Bq�Bv�By�Bz�B{�B� B�B�%B�%B�1B�7B�JB�VB�DB�bB�oB�hB�hB�hB�hB�oB�oB�uB�uB��B��B��B��B��B��B��B��B�B�B�-B�-B�9B�9B�?B�9B�9B�-B�-B�'B�!B�B�B�B�B�B�3B�FB�XB��B��B��B��B��B��B��B��B��B��B��B�B�B�
B�B�#B�#B�#B�B�B�B�B�B��B��B��B��B��BɺBȴBǮB��B��BȴBȴBǮBƨBŢBÖBĜBƨBĜBB��B�qB�RB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�hB�oB�hB�bB�\B�VB�VB�VB�\B�bB�{B��B��B�{B�{B�uB�oB�oB�bB�VB�JB�PB�VB�\B�\B�VB�\B�\B�\B�\B�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�JB�PB�PB�VB�PB�PB�PB�PB�PB�\B�VB�PB�VB�\B�VB�\B�VB�VB�PB�JB�=B�7B�7B�7B�7B�7B�=B�7B�7B�=B�=B�=B�7B�7B�1B�1B�7B�7B�7B�7B�=B�1B�%B�%B�%B�B�B�B�B�B�B�B}�B}�B~�B|�B{�Bz�Bz�By�By�By�By�Bx�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bv�Bw�By�B{�B|�B{�B{�Bz�Bz�Bz�B|�B~�B� B~�B� B�B�B� B~�B~�B~�B~�B}�B}�B|�B|�B}�B~�B~�B~�B}�B{�Bz�Bz�By�By�B{�B� B� B� B~�B{�Bz�B{�B{�B{�Bz�By�Bx�Bx�By�By�Bx�Bx�Bv�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bq�Bp�Br�Bs�Bs�Bt�Bu�Bv�Bv�Bw�B{�Bz�Bz�Bz�B}�B� B�B�%B�=B�PB�PB�VB�VB�VB�PB�JB�PB�hB��B��B��B��B�B�B�B�'B�'B�B�B�B�B�B�B�B�B�B�!B�-B�9B�FB�RB�LB�RB�RB�XB�XB�RB�RB�RB�RB�RB�RB�RB�RB�XB�RB�FB�3B�-B�-B�-B�-B�3B�3B�3B�9B�9B�9B�9B�3B�-B�-B�'B�'B�'B�'B�'B�'B�'B�'B�-B�3B�3B�3B�3B�-B�-B�'B�'B�'B�-B�-B�-B�-B�'B�'B�!B�!B�!B�!B�!B�'B�-B�-B�-B�-B�-B�'B�'B�'B�!B�'B�'B�!B�!B�'B�'B�'B�'B�'B�'B�'B�-B�-B�-B�-B�-B�3B�9B�?B�FB�FB�LB�LB�LB�LB�FB�FB�FB�LB�LB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�RB�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811071254112018110712541120181107125411  IF  ARFMCODA024c                                                                20181105182046                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105182133  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105182133  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107125411  IP  PSAL            @ffD�@ G�O�                