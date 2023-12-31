CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090453  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               )A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�`PZ�.;1   @�`P�N �@Oߺ���`�A��QvX1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @   @@  @�33@�  @�  @���A   A��A!��A0  A@  AP  A`  Ap  A���A���A���A���A�  A�  A�33A�  A���A�  A�  A�  A�33A�  A���A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B;��B?��BC��BH  BL  BP  BT  BW��B[��B_��Bc��Bh  Bl  Bp  Bs��Bw��B{��B�  B�33B�33B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	ffC  C��C�C� C  C� C  C� C �C"� C%  C'��C*  C,� C/  C1� C4  C6ffC9  C;��C>  C@��CC  CE� CH  CJffCL�fCO� CR  CT� CW  CY��C\  C^� Ca  Cc� Cf  Ch��Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� DfDFfD	� D
� DfD@ D� D� D  D@ D� D� DfD@ D� D�fDfD@ D� D�fD   D!@ D"� D#� D%  D&@ D'y�D(� D*fD+FfD,� D-� D.��D09�D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF� DH  DI@ DJ� DK� DM  DN9�DO� DP� DR  DSFfDT�fDU� DV��DX@ DY� DZ� D\  D]@ D^y�D_� DafDbFfDc�fDd� De��Dg@ Dh� Di� DkfDlFfDm�fDn� Dp  Dq9�Dry�Ds� DufDv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�  D�� D�\�D���D���D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D��3D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  Dã3D�@ D�� Dŀ D��D�� D�` D�  DȜ�D�@ D��3Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dσ3D�#3D�� D�` D���DҠ D�@ D�� DԀ D�  D�� D�c3D�  Dנ D�@ D�� Dك3D�#3D�� D�\�D�  Dܣ3D�C3D��3Dހ D��D�� D�` D�  D� D�C3D��3D� D��D�� D�` D�  D� D�@ D�� D�3D�  D��3D�` D�  D�3D�C3D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�  D��D�` D�3D�� D�C3D�� D�|�D�  D��3D�` D�3D��fD�C3D��3D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @   @@  @�33@�  @�  @���A   A��A!��A0  A@  AP  A`  Ap  A���A���A���A���A�  A�  A�33A�  A���A�  A�  A�  A�33A�  A���A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B;��B?��BC��BH  BL  BP  BT  BW��B[��B_��Bc��Bh  Bl  Bp  Bs��Bw��B{��B�  B�33B�33B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	ffC  C��C�C� C  C� C  C� C �C"� C%  C'��C*  C,� C/  C1� C4  C6ffC9  C;��C>  C@��CC  CE� CH  CJffCL�fCO� CR  CT� CW  CY��C\  C^� Ca  Cc� Cf  Ch��Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� DfDFfD	� D
� DfD@ D� D� D  D@ D� D� DfD@ D� D�fDfD@ D� D�fD   D!@ D"� D#� D%  D&@ D'y�D(� D*fD+FfD,� D-� D.��D09�D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@� DA� DC  DDFfDE� DF� DH  DI@ DJ� DK� DM  DN9�DO� DP� DR  DSFfDT�fDU� DV��DX@ DY� DZ� D\  D]@ D^y�D_� DafDbFfDc�fDd� De��Dg@ Dh� Di� DkfDlFfDm�fDn� Dp  Dq9�Dry�Ds� DufDv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�  D�� D�\�D���D���D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D��3D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  Dã3D�@ D�� Dŀ D��D�� D�` D�  DȜ�D�@ D��3Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dσ3D�#3D�� D�` D���DҠ D�@ D�� DԀ D�  D�� D�c3D�  Dנ D�@ D�� Dك3D�#3D�� D�\�D�  Dܣ3D�C3D��3Dހ D��D�� D�` D�  D� D�C3D��3D� D��D�� D�` D�  D� D�@ D�� D�3D�  D��3D�` D�  D�3D�C3D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�  D��D�` D�3D�� D�C3D�� D�|�D�  D��3D�` D�3D��fD�C3D��3D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�/@�/@�/@�/@�/@�/@�7L@�/@�7L@�7L@��@��/@��`@��/@���@���@���@���@��/@���@��9@�@�C�@��#@��`@�j@�1@�@�S�@�@�\@�V@�5?@�$�@웦@���@�O�@�@⟾@��T@�j@�;d@�1@�{@���@�Q�@�l�@�=q@�j@�`B@�/@̼j@�bN@�Ĝ@���@���@�v�@�~�@Ώ\@���@�
=@�o@��y@�ff@ͩ�@̛�@�1'@�(�@���@��@�9X@�Q�@�I�@�r�@ȣ�@ȼj@���@���@ȼj@�Q�@�9X@�1@��m@Ǿw@ǝ�@Ǖ�@�dZ@�dZ@�\)@�;d@��H@ư!@�ff@���@ũ�@ř�@Ł@�p�@�G�@�/@��@��@ģ�@�j@�j@�j@��m@�\)@�;d@�@��@��@¸R@�V@�@���@�/@��9@�(�@��w@�dZ@�+@��@�^5@��h@��@��w@��@�l�@�S�@�K�@�C�@�C�@�+@�@�ȴ@�v�@�5?@��7@�7L@���@��@��@��F@���@��@�t�@�K�@�C�@�C�@�+@�;d@��@���@�v�@�V@�M�@�E�@�-@�{@��@���@��7@�O�@��j@�ƨ@�
=@��y@�ȴ@���@�ff@�=q@���@���@��^@���@�x�@�hs@�O�@�&�@��`@��@�bN@�bN@�I�@�1'@�b@�  @��m@��w@�l�@�;d@�@��R@�E�@���@��@��@���@�  @���@�ȴ@�~�@��#@��h@�x�@�?}@�V@���@��u@�r�@�bN@���@���@��\@�^5@�E�@�{@���@��T@���@�@��^@���@���@��h@�p�@�G�@��`@��9@�A�@��@�t�@�C�@�33@�"�@��@�E�@�J@��@���@��h@�G�@��`@���@�r�@�1'@��
@���@��@��!@�^5@��T@�O�@��j@�I�@�b@��m@��F@��@�v�@�x�@�I�@�1@��@��F@�K�@���@���@��@��m@���@�S�@��H@�M�@�7L@��@��D@���@�n�@�-@���@��@���@��j@���@��u@�(�@�\)@�@��H@�^5@��@�@���@�X@��@�%@��/@�Q�@�b@��@��;@�ƨ@��@��@��!@��@��#@��^@���@�&�@���@��@���@�r�@�Q�@�Z@�Q�@�A�@� �@�  @��@��@���@��+@�ff@�=q@�{@�@�@�G�@�%@��`@��@�Z@�A�@�(�@�  @�;@�@|�@;d@;d@�@~��@~�R@~��@~��@~5?@}�h@}V@|�@|z�@{�m@{dZ@{33@z�@z�@z��@z^5@z�@y��@y�7@yx�@yhs@y&�@y&�@yx�@y��@y��@y�^@y�^@y��@yG�@x��@xĜ@w�@w��@w�@w|�@w;d@w;d@w+@w�@w�@w�@w+@w;d@wK�@w\)@w|�@w�;@xb@w�@w�P@w�P@w|�@w|�@w|�@w|�@wl�@wK�@w+@w�@w
=@v��@v�R@v��@vff@v@u�-@u?}@t��@tj@t1@s�m@s�
@s�m@t1@s��@st�@s33@r��@r^5@r-@rJ@q��@q�@q�^@qhs@q7L@q�@p�@n�y@n�+@m�T@mp�@l��@lZ@lI�@l(�@l1@l1@l1@jJ@h��@g��@f�@f�+@f@e�@d(�@b�H@bJ@a�^@`bN@`  @]��@\�@\�j@\��@\z�@\Z@[�@[S�@[�@[C�@Z�H@Yhs@W��@W�@W�;@W�@W��@WK�@Vv�@V5?@U@U�h@UO�@UO�@U/@U��@T��@T��@S��@S�m@SC�@Qhs@Q&�@Q7L@Q&�@P�u@P1'@O�w@O\)@L��@EO�@D9X@CS�@Ahs@?K�@=�h@<�@<Z@;o@:J@9��@:-@:�\@:n�@9�@9��@9�^@9G�@7�@6�R@5�@4��@4(�@4�@3S�@3@2�@2�H@2�H@2�H@2�@2^5@0�`@.@,j@+��@+S�@+C�@&��@#�F@#t�@l�@�R@@`B@�D@�
@�H@�\@^5@?}@
-@�9@��@n�?���?���?��m?�7L?�`B?ޗ�?Լj?�M�?��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @���@�/@�/@�/@�/@�/@�/@�7L@�/@�7L@�7L@��@��/@��`@��/@���@���@���@���@��/@���@��9@�@�C�@��#@��`@�j@�1@�@�S�@�@�\@�V@�5?@�$�@웦@���@�O�@�@⟾@��T@�j@�;d@�1@�{@���@�Q�@�l�@�=q@�j@�`B@�/@̼j@�bN@�Ĝ@���@���@�v�@�~�@Ώ\@���@�
=@�o@��y@�ff@ͩ�@̛�@�1'@�(�@���@��@�9X@�Q�@�I�@�r�@ȣ�@ȼj@���@���@ȼj@�Q�@�9X@�1@��m@Ǿw@ǝ�@Ǖ�@�dZ@�dZ@�\)@�;d@��H@ư!@�ff@���@ũ�@ř�@Ł@�p�@�G�@�/@��@��@ģ�@�j@�j@�j@��m@�\)@�;d@�@��@��@¸R@�V@�@���@�/@��9@�(�@��w@�dZ@�+@��@�^5@��h@��@��w@��@�l�@�S�@�K�@�C�@�C�@�+@�@�ȴ@�v�@�5?@��7@�7L@���@��@��@��F@���@��@�t�@�K�@�C�@�C�@�+@�;d@��@���@�v�@�V@�M�@�E�@�-@�{@��@���@��7@�O�@��j@�ƨ@�
=@��y@�ȴ@���@�ff@�=q@���@���@��^@���@�x�@�hs@�O�@�&�@��`@��@�bN@�bN@�I�@�1'@�b@�  @��m@��w@�l�@�;d@�@��R@�E�@���@��@��@���@�  @���@�ȴ@�~�@��#@��h@�x�@�?}@�V@���@��u@�r�@�bN@���@���@��\@�^5@�E�@�{@���@��T@���@�@��^@���@���@��h@�p�@�G�@��`@��9@�A�@��@�t�@�C�@�33@�"�@��@�E�@�J@��@���@��h@�G�@��`@���@�r�@�1'@��
@���@��@��!@�^5@��T@�O�@��j@�I�@�b@��m@��F@��@�v�@�x�@�I�@�1@��@��F@�K�@���@���@��@��m@���@�S�@��H@�M�@�7L@��@��D@���@�n�@�-@���@��@���@��j@���@��u@�(�@�\)@�@��H@�^5@��@�@���@�X@��@�%@��/@�Q�@�b@��@��;@�ƨ@��@��@��!@��@��#@��^@���@�&�@���@��@���@�r�@�Q�@�Z@�Q�@�A�@� �@�  @��@��@���@��+@�ff@�=q@�{@�@�@�G�@�%@��`@��@�Z@�A�@�(�@�  @�;@�@|�@;d@;d@�@~��@~�R@~��@~��@~5?@}�h@}V@|�@|z�@{�m@{dZ@{33@z�@z�@z��@z^5@z�@y��@y�7@yx�@yhs@y&�@y&�@yx�@y��@y��@y�^@y�^@y��@yG�@x��@xĜ@w�@w��@w�@w|�@w;d@w;d@w+@w�@w�@w�@w+@w;d@wK�@w\)@w|�@w�;@xb@w�@w�P@w�P@w|�@w|�@w|�@w|�@wl�@wK�@w+@w�@w
=@v��@v�R@v��@vff@v@u�-@u?}@t��@tj@t1@s�m@s�
@s�m@t1@s��@st�@s33@r��@r^5@r-@rJ@q��@q�@q�^@qhs@q7L@q�@p�@n�y@n�+@m�T@mp�@l��@lZ@lI�@l(�@l1@l1@l1@jJ@h��@g��@f�@f�+@f@e�@d(�@b�H@bJ@a�^@`bN@`  @]��@\�@\�j@\��@\z�@\Z@[�@[S�@[�@[C�@Z�H@Yhs@W��@W�@W�;@W�@W��@WK�@Vv�@V5?@U@U�h@UO�@UO�@U/@U��@T��@T��@S��@S�m@SC�@Qhs@Q&�@Q7L@Q&�@P�u@P1'@O�w@O\)@L��@EO�@D9X@CS�@Ahs@?K�@=�h@<�@<Z@;o@:J@9��@:-@:�\@:n�@9�@9��@9�^@9G�@7�@6�R@5�@4��@4(�@4�@3S�@3@2�@2�H@2�H@2�H@2�@2^5@0�`@.@,j@+��@+S�@+C�@&��@#�F@#t�@l�@�R@@`B@�D@�
@�H@�\@^5@?}@
-@�9@��@n�?���?���?��m?�7L?�`B?ޗ�?Լj?�M�?��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB33B33B33B33B33B33B33B33B33B2-B2-B33B2-B2-B2-B2-B2-B1'B1'B1'B0!B/B-B1'B0!B0!B1'B0!B1'B1'B1'B2-B1'B/B-B0!B7LB;dB;dB;dB:^B9XB:^B?}B>wB<jB<jB:^B:^B9XB?}B@�BB�BD�BN�BT�B[#BbNBbNBbNBcTBe`BgmBgmBgmBdZBbNB^5B\)B]/BYB^5BaHBbNBe`BgmBgmBgmBhsBhsBhsBiyBiyBiyBiyBiyBiyBiyBhsBhsBiyBjBiyBjBk�Bk�Bk�Bk�Bk�Bl�Bl�Bl�Bl�Bl�Bl�Bn�Bp�Bp�Bo�Bo�Bo�Bn�Bn�Bm�Bl�Bk�BjBiyBhsBgmBffBe`BdZBcTBaHB^5B\)BZBYBYBYBXBXBXBXBW
BVBT�BS�BR�BQ�BP�BN�BN�BM�BM�BM�BL�BL�BL�BL�BL�BK�BK�BJ�BJ�BI�BI�BI�BI�BH�BH�BH�BG�BF�BD�BC�BA�B@�B@�B@�B?}B?}B>wB>wB=qB=qB=qB<jB<jB<jB;dB:^B:^B9XB9XB9XB9XB8RB8RB7LB7LB6FB5?B49B33B2-B1'B0!B.B,B+B'�B&�B%�B$�B#�B#�B"�B!�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBoBoBoBhB\B\BVBVBPBPBJBDB
=B
=B1B1B%BBBBB  B��B��B��B��B��B��B��B��B�B�B�B�B�B�yB�sB�mB�fB�`B�TB�HB�5B�/B�#B�B�B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBǮBǮBƨBƨBĜBĜBÖBÖBÖBÖBB��B�wB�qB�qB�jB�dB�dB�^B�^B�XB�XB�XB�RB�RB�RB�LB�FB�?B�9B�9B�3B�3B�-B�-B�'B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�'B�'B�-B�-B�3B�9B�9B�?B�?B�FB�LB�XB�dB�jB�qB�wB�wB�}B�}B�}B��B��B��B��BBBBBBBBBBBÖBÖBĜBƨBǮBȴBɺBɺBɺBɺB��B��B��B��B��B��BɺBɺBɺBǮBƨBƨBǮBȴBȴBȴBȴBȴBȴBȴBȴBǮBǮBƨBƨBƨBƨBŢBŢBŢBĜBĜBÖBÖBÖBÖBÖBÖBBBBB��B��B�wB�jB�qB�jB�jB�jB�^B�XB�RB�LB�FB�FB�FB�FB�FB�?B�9B�-B�-B�!B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B)�B)�B)�B)�B)�B)�B)�B)�B)�B(�B(�B)�B(�B(�B(�B(�B(�B'�B'�B'�B&�B%�B#�B'�B&�B&�B'�B&�B'�B'�B'�B(�B'�B%�B#�B&�B-�B1�B1�B2 B0�B/�B0�B6B5B3B3B0�B0�B/�B6B7 B9+B;8BEuBK�BQ�BX�BX�BX�BY�B[�B^B^B^BZ�BX�BT�BR�BS�BO�BT�BW�BX�B[�B^B^B^B_B_B_B`B`B`B`B`B`B`B_B_B`BaB`BaBbBbBbBbBbBc%Bc%Bc%Bc%Bc%Bc%Be2Bg>Bg>Bf8Bf8Bf8Be2Be2Bd+Bc%Bb BaB`B_B^B]B[�BZ�BY�BW�BT�BR�BP�BO�BO�BO�BN�BN�BN�BN�BM�BL�BK�BJ�BI�BH�BG�BEvBEvBDpBDpBDpBCjBCjBCjBCjBCjBBdBBdBA^BA^B@XB@XB@XB@XB?RB?RB?RB>LB=FB;:B:4B8'B7!B7!B7!B6B6B5B5B4B4B4B3	B3	B3	B2B0�B0�B/�B/�B/�B/�B.�B.�B-�B-�B,�B+�B*�B)�B(�B'�B&�B$�B"�B!�B�B�B�B~BxBxBrBlBlBfB_BYBNBHBHBBBBBBBBB<B<B<B<B<B5B5B/B)B)B#BB
B	B	B	BB�B�B�B�B�B�B�B�B �B �B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�nB�aB�\B�VB�VB�JB�>B� B�B�B�B�B��B��B��B��B��BϿB̬B˦BʡBȕBȕBǎBǎBǎBłB�|B�vB�pB�kB�dB�^B�^B�XB�XB�RB�RB�FB�FB�@B�@B�@B�@B�9B�.B�"B�B�B�B�B�B�	B�	B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�#B�#B�(B�(B�(B�.B�.B�4B�4B�:B�:B�:B�:B�:B�:B�:B�:B�:B�:B�AB�AB�GB�SB�YB�_B�eB�eB�eB�eB�lB�rB�rB�rB�lB�lB�eB�eB�eB�YB�SB�SB�YB�_B�_B�_B�_B�_B�_B�_B�_B�YB�ZB�TB�TB�TB�TB�NB�NB�NB�HB�HB�BB�BB�BB�BB�BB�BB�;B�;B�;B�;B�5B�0B�$B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�hB�bB�bB�\B�\B�VB�VB�PB�KB�KB�KB�\B�hB�uB�nB�nB�nB�nB�bB�bB�cB�]B�WB�QB�QB�QB�QB�QB�QB�QB�WB�]B�cB�cB�iB�iB�iB�cB�vB��B��B��B��B��B��B��B��B��B��B�}B�pB�jB�dB�^B�YB�YB�eB�eB�eB�eB�kB�~B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0091338                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904532018110709045320181107090453  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172720  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172720  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090453  IP  PSAL            @   D�� G�O�                