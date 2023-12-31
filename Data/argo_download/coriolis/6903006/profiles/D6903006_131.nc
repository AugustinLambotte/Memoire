CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  N   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:24Z creation; 2023-01-04T00:23:11Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_054b      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    :�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    :�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  ;   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  ;D   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  ;�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        ;�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     <   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     <8   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    <X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         <\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    <d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            <h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           <p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           <x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    <�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    <�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    <�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    <�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        =�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     
p  =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  H   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        8  I`   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  N�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        8  O�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  U    PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     8  Vp   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  [�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  `�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  b0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  gh   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  h�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  m�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  s(   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  tx   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  y�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     8  {    HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �    HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �$   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �(   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �,   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �0   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  �8   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �x   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �x   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �x   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  �xArgo profile    3.1 1.2 19500101000000  20230104002024  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @��>��1   @��>��@R�B;��h�&f�u�8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 2000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    A.�~    �H���   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �s333@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����7   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���5�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���O��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���H(  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��#Eh  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����0  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��GL�X  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��X^h  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���5��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��[�ޠ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���A;�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���t��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��	{B^  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��l�l  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���˩�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��`T�>  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��Յ��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��""""  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��N���  9990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990  A&ffA1��A;33AK33A[33Al��A~ffA�  A���A���A���A�  A�  A�  A���A�33A�33A���A�33A陚A�  A���B ffBffB��B��B33B��BffB33B   B$  B'��B+33B0ffB4  B8  B<  B@  BD  BHffBTffBhffB}33B�  B���B���B���B�ffB���B���B�  B���BᙚB뙚B���B���C  C
33C�fC�3C�fCL�C#33C(33C-�C2  C6�fC<  C@��CE�3CJ�3CO��CUL�CZL�C^�fCc� Ch� Cm��Cr�3Cw��C|��C��fC�s3C�  C���C��C���C�&fC�� C�ٚC�s3C��C�� C��fC�Y�C�ٚC�s3C�  C���C�ٚC�s3C�  C���C��C�ffC�  C���C��fCĀ C��Cɀ C��3C�s3C��fC�ffC��fC�Y�C�ٚC�@ C��fC� C��3C� C��fC�L�C��3C�C��C���C��C���C��D L�D�3D��D�fD,�Dy�D�fD��D
&fDy�D��D��D@ Dl�D��D3DS3D�3D�3D�DY�D� D�fD��D33Dy�D � D!�3D#&fD$` D%�fD&�fD(,�D)l�D*��D+��D-33D.� D/�3D0�fD29�D3y�D4�3D5�3D733D8y�D9� D;�D<9�D=l�D>� D@3DAS3DB�3DC��DE�DFFfDG� DH��DI��DK33DLs3DM�3DN�3DP33DQy�DR� DT  DU33DV� DW��DYfDZ@ D[� D\� D^fD_33D`ffDa��Db�3Dd,�De�fDf�fDg��Di33Djy�Dk� Dm�Dn@ Doy�Dp��Dq� Ds9�Dt�3Du��Dw3DxY�Dy��D{s3D}ٚD�0 D�ffD���D��D�&fD�s3D��fD�� D�0 D�p D�� D��3D�,�D�p D��3D���D�,�D�l�D��fD��D�)�D�l�D��fD�� D�,�D�l�D��3D���D�)�D�vfD��fD���D�0 D�l�D���D���D�0 D�l�D���D��fD�)�D�s3D���D�� D�,�D�l�D��fD���D�&fD�ffD���D�� D�,�D�l�Dì�D��3D�,�D�ffDȩ�D��fD�&fD�ffDͰ D��D�)�D�vfDҰ D��3D�  D�p D׳3D�� D�0 D�vfDܳ3D��3D�6fD�l�D�fD��3D�#3D�p D�fD���D�9�D�s3D�fD���D�0 D�c3D�3D��3D�)�D�c3D�� D�� D�#3D�i�D�� D��3D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A&ffA1��A;33AK33A[33Al��A~ffA�  A���A���A���A�  A�  A�  A���A�33A�33A���A�33A陚A�  A���B ffBffB��B��B33B��BffB33B   B$  B'��B+33B0ffB4  B8  B<  B@  BD  BHffBTffBhffB}33B�  B���B���B���B�ffB���B���B�  B���BᙚB뙚B���B���C  C
33C�fC�3C�fCL�C#33C(33C-�C2  C6�fC<  C@��CE�3CJ�3CO��CUL�CZL�C^�fCc� Ch� Cm��Cr�3Cw��C|��C��fC�s3C�  C���C��C���C�&fC�� C�ٚC�s3C��C�� C��fC�Y�C�ٚC�s3C�  C���C�ٚC�s3C�  C���C��C�ffC�  C���C��fCĀ C��Cɀ C��3C�s3C��fC�ffC��fC�Y�C�ٚC�@ C��fC� C��3C� C��fC�L�C��3C�C��C���C��C���C��D L�D�3D��D�fD,�Dy�D�fD��D
&fDy�D��D��D@ Dl�D��D3DS3D�3D�3D�DY�D� D�fD��D33Dy�D � D!�3D#&fD$` D%�fD&�fD(,�D)l�D*��D+��D-33D.� D/�3D0�fD29�D3y�D4�3D5�3D733D8y�D9� D;�D<9�D=l�D>� D@3DAS3DB�3DC��DE�DFFfDG� DH��DI��DK33DLs3DM�3DN�3DP33DQy�DR� DT  DU33DV� DW��DYfDZ@ D[� D\� D^fD_33D`ffDa��Db�3Dd,�De�fDf�fDg��Di33Djy�Dk� Dm�Dn@ Doy�Dp��Dq� Ds9�Dt�3Du��Dw3DxY�Dy��D{s3D}ٚD�0 D�ffD���D��D�&fD�s3D��fD�� D�0 D�p D�� D��3D�,�D�p D��3D���D�,�D�l�D��fD��D�)�D�l�D��fD�� D�,�D�l�D��3D���D�)�D�vfD��fD���D�0 D�l�D���D���D�0 D�l�D���D��fD�)�D�s3D���D�� D�,�D�l�D��fD���D�&fD�ffD���D�� D�,�D�l�Dì�D��3D�,�D�ffDȩ�D��fD�&fD�ffDͰ D��D�)�D�vfDҰ D��3D�  D�p D׳3D�� D�0 D�vfDܳ3D��3D�6fD�l�D�fD��3D�#3D�p D�fD���D�9�D�s3D�fD���D�0 D�c3D�3D��3D�)�D�c3D�� D�� D�#3D�i�D�� D��3D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@������P��l���+��+��+��
=���y��
=��
=��+��
=��+���y��$ݿ�$ݿ��Ͽ�t���S���S���S���t���t���S���F�㕁��F���
���Ͽ�����������Ͽ�����/��`B��`B��?}������Z��9X��n���%�����1�ա˿�G����Ϳ�/���y��$ݿƇ+�š˿����hr��+C���P�z���������=q���/�'+���þ��������P�hs�-V������
=\)>��?V?F$�?O�?PbN?P �?7�P?"M�?�F?)x�?)�^?-O�?/\)?+C�?!%? �?|�?5??^5?��?�?
=q?`B>�`B>�ff>>���>�"�>��>š�>\>���>Õ�>�v�>��>��>�O�>��R>�l�>�1>��>�5?>�n�>�C�>�I�>���>��>��>u>fff>S��>%�T>�+=ě�=�9X=�x�>J>J>I�>   =�h=�x�=�/=�E�=�j=��`=�G�>%=>�>1'>�>�>�=��=���=�l�=�/=��=�{=���=�C�=P�`=,1<��<�1<��
<�9X<�C�<u<o:�o    ��o�ě��o�#�
�������w�#�
��P��w�,1�,1�0 Ž@��H�9�]/�]/�ixս����\)�������罰 Ž�9X��^5���ͽ����;d��S���h��F�����پ   �+�C��I��t������-��w�#�
�(�þ,1�1&�6E��;dZ�@��F��H�9�Kƨ�L�;Kƨ�M��N��Q녾T���Xb�\(��^5?�bMӾgl��hr��hr��ixվl�D�p�׾w�پ{�m���7������˾�$ݾ��;�\)��n������
=������5?��;d��A���Z���y��~���{�� ž������!��?}���#���m��|�����1'��C������hs��z����ۥ��;d������`B��r���h��&��9X���پ�푾�vɿ ��G��%����o�Z����T���1'�
~��1���bN�&�녿33����?}�b���dZ����#�"ѿ(��p��vɿ;d�!�7�#o�$Z�%��&�y�(1'�)xտ*=q�,�Ϳ.V�/��0bN�1&�2�!�5��6ȴ�7Kǿ7�ٿ8Q�8�u�8���9X�9�#�9�#�:^5�;"ѿ;�m�<j�<�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ����P��l���+��+��+��
=���y��
=��
=��+��
=��+���y��$ݿ�$ݿ��Ͽ�t���S���S���S���t���t���S���F�㕁��F���
���Ͽ�����������Ͽ�����/��`B��`B��?}������Z��9X��n���%�����1�ա˿�G����Ϳ�/���y��$ݿƇ+�š˿����hr��+C���P�z���������=q���/�'+���þ��������P�hs�-V������
=\)>��?V?F$�?O�?PbN?P �?7�P?"M�?�F?)x�?)�^?-O�?/\)?+C�?!%? �?|�?5??^5?��?�?
=q?`B>�`B>�ff>>���>�"�>��>š�>\>���>Õ�>�v�>��>��>�O�>��R>�l�>�1>��>�5?>�n�>�C�>�I�>���>��>��>u>fff>S��>%�T>�+=ě�=�9X=�x�>J>J>I�>   =�h=�x�=�/=�E�=�j=��`=�G�>%=>�>1'>�>�>�=��=���=�l�=�/=��=�{=���=�C�=P�`=,1<��<�1<��
<�9X<�C�<u<o:�o    ��o�ě��o�#�
�������w�#�
��P��w�,1�,1�0 Ž@��H�9�]/�]/�ixս����\)�������罰 Ž�9X��^5���ͽ����;d��S���h��F�����پ   �+�C��I��t������-��w�#�
�(�þ,1�1&�6E��;dZ�@��F��H�9�Kƨ�L�;Kƨ�M��N��Q녾T���Xb�\(��^5?�bMӾgl��hr��hr��ixվl�D�p�׾w�پ{�m���7������˾�$ݾ��;�\)��n������
=������5?��;d��A���Z���y��~���{�� ž������!��?}���#���m��|�����1'��C������hs��z����ۥ��;d������`B��r���h��&��9X���پ�푾�vɿ ��G��%����o�Z����T���1'�
~��1���bN�&�녿33����?}�b���dZ����#�"ѿ(��p��vɿ;d�!�7�#o�$Z�%��&�y�(1'�)xտ*=q�,�Ϳ.V�/��0bN�1&�2�!�5��6ȴ�7Kǿ7�ٿ8Q�8�u�8���9X�9�#�9�#�:^5�;"ѿ;�m�<j�<�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�B	�B	�B	�B	 �B	 �B	 �B	!�B	 �B	�B	 �B	 �B	 �B	"�B	-B	(�B	A�B	E�B	C�B	F�B	G�B	I�B	I�B	J�B	M�B	L�B	L�B	N�B	Q�B	S�B	S�B	S�B	T�B	VB	YB	\)B	_;B	`BB	bNB	ffB	n�B	�\B	��B	�!B	�jB	��B	�B	�`B	��B
DB
oB
�B
!�B
-B
S�B
iyB
v�B
gmB
{�B
z�B
� B
�{B
�9B
��B
�B
�HB
��B
��BB1BbB �B2-BF�Bm�Bs�Bt�Bs�Bp�Bn�Br�Bz�B|�B�B�B�B�B�B�+B�+B�1B�7B�7B�7B�=B�%B�+B�7B�DB�=B�DB�DB�DB�PB�bB�hB�bB�=B�VB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�oB�VB�VB�hB�uB�uB�uB�{B�oB�uB�uB�oB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B	�B	�B	�B	�B	 �B	 �B	 �B	!�B	 �B	�B	 �B	 �B	 �B	"�B	-B	(�B	A�B	E�B	C�B	F�B	G�B	I�B	I�B	J�B	M�B	L�B	L�B	N�B	Q�B	S�B	S�B	S�B	T�B	VB	YB	\)B	_;B	`BB	bNB	ffB	n�B	�\B	��B	�!B	�jB	��B	�B	�`B	��B
DB
oB
�B
!�B
-B
S�B
iyB
v�B
gmB
{�B
z�B
� B
�{B
�9B
��B
�B
�HB
��B
��BB1BbB �B2-BF�Bm�Bs�Bt�Bs�Bp�Bn�Br�Bz�B|�B�B�B�B�B�B�+B�+B�1B�7B�7B�7B�=B�%B�+B�7B�DB�=B�DB�DB�DB�PB�bB�hB�bB�=B�VB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�oB�VB�VB�hB�uB�uB�uB�{B�oB�uB�uB�oB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002024                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002311  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002311  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A&ffD�� G�O�                