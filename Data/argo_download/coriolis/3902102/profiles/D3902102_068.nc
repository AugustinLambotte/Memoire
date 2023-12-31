CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  @   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:24Z creation; 2020-11-17T12:19:03Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8(   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8h   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9,   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         90   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    98   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9<   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9D   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9L   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9T   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9X   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9`   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9d   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9h   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9l   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :l   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	   :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  Cp   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	   E�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  N�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	   P�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   Y�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  b�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   e0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  n0   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   pp   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   yp   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  �p   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 @  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	   ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �L   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �P   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �T   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �X   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �\   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �    SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �    SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �    SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  � Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170406  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               DA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @�W�8�1   @�W�8�@R��_I��(� 4u<8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A0  A@  AQ��A`  Ap  A~ffA���A���A���A���A�  A�  A���A���A�33A�  A�  A�33A�ffA�33A���B   BffBffB��BffBffBffB��B��B#33B'��B+��B0  B4  B8ffB<  B?33BC��BH  BK��BO33BT  BW��B[33B`ffBd  Bh  Bl  Bo��Bs��Bw��B{��B��B�  B�  B�33B���B���B�  B���B���B���B�ffB���B�  B���B�  B�ffB�ffB�  B���B�ffB�33B���B���B�ffB�33B�  B���B���B�ffB�ffB�33B�  B���B�Bę�B�ffB�ffB�33B�  B�  B�  B�  B�33B�33B�33B�33B�ffB�B���B���B���C  C�C�C33C	L�CffC� C�3C�3C��CffC  C�C33C33CffC!ffC#ffC%ffC'ffC)ffC+� C-� C/� C1�3C3L�C5ffC7� C9� C;� C=��C?��CA33CCL�CEffCGffCI� CK��CM�3COL�CQffCS� CU��CW�3CYL�C[ffC]��C_L�Ca� Cc�3CeffCg33CiffCk��CmL�Co� Cq�3Cs� Cu33CwffCy��C{L�C}ffC��C��fC�� C�ٚC��3C���C�� C�ٚC�� C��fC���C��3C���C��3C���C�� C�ٚC�� C��fC���C�� C�ٚC���C�� C��fC���C���C���C�� C�� C�� C�� C��3C��3C��fC��fC��fC��fC��fC��fC�ٚC�ٚC���C���C���C���C���C���C���C���C���C���C���C���C���C�ٚC�ٚC��fC��3C��3C�� C���C���C���C��fC��3C���C���C�� C�� Cų3CƳ3Cǳ3Cȳ3Cɳ3Cʳ3C�� C�� C���C���CϦfCЦfC�� C�� C�ٚCԳ3C�� C�ٚC׳3C�� C�ٚCڳ3Cۙ�Cܳ3C���C޳3Cߌ�C�3C�� C�fC��C�3C�ٚC�� C�fC��C�3C�ٚC�� C�fC��C�� C��fC�ٚC�� C�3C�C��C�� C��3C��fC�ٚC���C�s3C�ٚD &fD� D��D  D9�Dl�D�fD� D
  D` D� D� D  D` D� D� D  D` D�fD��D,�Ds3D� D�3D9�D� D ��D"  D#9�D$y�D%�3D'fD(FfD)� D*� D,  D-FfD.l�D/��D1fD29�D3l�D4�fD5��D7S3D8�3D9ٚD;  D<,�D=y�D>��D?��DA9�DB�fDC��DD��DF,�DGffDH�fDI�fDK&fDLl�DM��DN�3DP9�DQ� DR��DT  DU33DV�fDW� DX�3DZ33D[s3D\�3D]��D_@ D`l�Da� DcfDd33De� Df�fDg�3Di@ Dj��Dk� Dl��Dn33Dol�Dp�fDq� Ds  Dt` Du�fDv��Dx,�Dys3Dz��D|�D}FfD~y�D�3D�|�D�#3D���D�VfD�  D��fD�C3D��3D��3D�fD��fD�Y�D���D�� D�9�D�� D�y�D�#3D���D�Y�D��fD��3D�0 D���D���D�,�D���D�\�D�� D��fD�6fD�� D��fD�#3D�� D�\�D���D��fD�9�D���D��3D��D��3D�VfD���D��3D�@ D���D�|�D�fD���D�\�D�  D���D�@ D��D���D�  D���D�VfD�  D���D�C3D���D�s3D��D��fD�` D���D��fD�33D��3D�s3D�3D�� D�P D��3D���D�<�D�ٚD�vfD�fD���D�\�D�  D��3D�<�D��fD�� D��D���D�VfD���D���D�@ D��3D�|�D�3D���D�c3D���D��fD�0 D���D��fD�#3D�� D�` D���Dę�D�@ D�ٚD�vfD�  D�ɚD�c3D�  Dɣ3D�FfD�� D�y�D��D̹�D�\�D��3DΙ�D�@ D���D�y�D�fDѶfD�VfD���Dә�D�9�D���D�|�D�  D��3D�Y�D�� Dؙ�D�@ D��fD�|�D�fD��3D�` D���DݖfD�33D�� D�|�D��D๚D�VfD�fD�3D�FfD��3D� D�#3D�3D�Y�D�  D�3D�<�D��fD�|�D�3D� D�\�D��D��D�L�D�� D��D�0 D�� D�\�D�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A0  A@  AQ��A`  Ap  A~ffA���A���A���A���A�  A�  A���A���A�33A�  A�  A�33A�ffA�33A���B   BffBffB��BffBffBffB��B��B#33B'��B+��B0  B4  B8ffB<  B?33BC��BH  BK��BO33BT  BW��B[33B`ffBd  Bh  Bl  Bo��Bs��Bw��B{��B��B�  B�  B�33B���B���B�  B���B���B���B�ffB���B�  B���B�  B�ffB�ffB�  B���B�ffB�33B���B���B�ffB�33B�  B���B���B�ffB�ffB�33B�  B���B�Bę�B�ffB�ffB�33B�  B�  B�  B�  B�33B�33B�33B�33B�ffB�B���B���B���C  C�C�C33C	L�CffC� C�3C�3C��CffC  C�C33C33CffC!ffC#ffC%ffC'ffC)ffC+� C-� C/� C1�3C3L�C5ffC7� C9� C;� C=��C?��CA33CCL�CEffCGffCI� CK��CM�3COL�CQffCS� CU��CW�3CYL�C[ffC]��C_L�Ca� Cc�3CeffCg33CiffCk��CmL�Co� Cq�3Cs� Cu33CwffCy��C{L�C}ffC��C��fC�� C�ٚC��3C���C�� C�ٚC�� C��fC���C��3C���C��3C���C�� C�ٚC�� C��fC���C�� C�ٚC���C�� C��fC���C���C���C�� C�� C�� C�� C��3C��3C��fC��fC��fC��fC��fC��fC�ٚC�ٚC���C���C���C���C���C���C���C���C���C���C���C���C���C�ٚC�ٚC��fC��3C��3C�� C���C���C���C��fC��3C���C���C�� C�� Cų3CƳ3Cǳ3Cȳ3Cɳ3Cʳ3C�� C�� C���C���CϦfCЦfC�� C�� C�ٚCԳ3C�� C�ٚC׳3C�� C�ٚCڳ3Cۙ�Cܳ3C���C޳3Cߌ�C�3C�� C�fC��C�3C�ٚC�� C�fC��C�3C�ٚC�� C�fC��C�� C��fC�ٚC�� C�3C�C��C�� C��3C��fC�ٚC���C�s3C�ٚD &fD� D��D  D9�Dl�D�fD� D
  D` D� D� D  D` D� D� D  D` D�fD��D,�Ds3D� D�3D9�D� D ��D"  D#9�D$y�D%�3D'fD(FfD)� D*� D,  D-FfD.l�D/��D1fD29�D3l�D4�fD5��D7S3D8�3D9ٚD;  D<,�D=y�D>��D?��DA9�DB�fDC��DD��DF,�DGffDH�fDI�fDK&fDLl�DM��DN�3DP9�DQ� DR��DT  DU33DV�fDW� DX�3DZ33D[s3D\�3D]��D_@ D`l�Da� DcfDd33De� Df�fDg�3Di@ Dj��Dk� Dl��Dn33Dol�Dp�fDq� Ds  Dt` Du�fDv��Dx,�Dys3Dz��D|�D}FfD~y�D�3D�|�D�#3D���D�VfD�  D��fD�C3D��3D��3D�fD��fD�Y�D���D�� D�9�D�� D�y�D�#3D���D�Y�D��fD��3D�0 D���D���D�,�D���D�\�D�� D��fD�6fD�� D��fD�#3D�� D�\�D���D��fD�9�D���D��3D��D��3D�VfD���D��3D�@ D���D�|�D�fD���D�\�D�  D���D�@ D��D���D�  D���D�VfD�  D���D�C3D���D�s3D��D��fD�` D���D��fD�33D��3D�s3D�3D�� D�P D��3D���D�<�D�ٚD�vfD�fD���D�\�D�  D��3D�<�D��fD�� D��D���D�VfD���D���D�@ D��3D�|�D�3D���D�c3D���D��fD�0 D���D��fD�#3D�� D�` D���Dę�D�@ D�ٚD�vfD�  D�ɚD�c3D�  Dɣ3D�FfD�� D�y�D��D̹�D�\�D��3DΙ�D�@ D���D�y�D�fDѶfD�VfD���Dә�D�9�D���D�|�D�  D��3D�Y�D�� Dؙ�D�@ D��fD�|�D�fD��3D�` D���DݖfD�33D�� D�|�D��D๚D�VfD�fD�3D�FfD��3D� D�#3D�3D�Y�D�  D�3D�<�D��fD�|�D�3D� D�\�D��D��D�L�D�� D��D�0 D�� D�\�D�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����`��%��%��G�������n������Mӿ�7��&��&��&��%��&��%��&��bN��|��|��|��|�޸R�ۥ���m�ߝ���/�ۥ�ٺ^��r��ؓu������1��(����m���ش9��l����Ͽ�G���;d��dZ���H��ff��o���7��r����7��S����7�������P��%������G���p��n������
=q�X��Z��"Ѿ���������#������D��Z=8Q�>���>���>��?^5?F$�?c�
?�b?�v�?��u?�?�Z?�E�?��u?�t�?�  ?�`B?̬?��y?���?�A�?���?���?ӶF?�V?�1'?�%?��H?�1'?���?���?�j?��?�^5?�5??�~�?�?��?���?���?�-?�A�?���?�b?��?�t�?�33?���?��?��-?� �?��?���?��
?��?�V?�j?�7L?�ȴ?��j?��!?���?���?�bN?�G�?�33?�33?�S�?���?�
=?�?}?�E�?���?��#?�p�?�/?���?�G�?�p�?�(�?�V?�/?�~�?��?��?��R?�p�?�dZ?���?��#?�1'?�K�?�E�?��?�%?�{?��?�1?�X?�K�?��9?���?�7L?�1'?���?��F?�n�?�M�?��!?���?��!?��!?��!?�n�?�M�?��`?��`?��?��?���?�&�?|�?{�m?y�#?y�?x��?xb?w�P?w
=?v?t�j?st�?st�?s��?t9X?t�j?s�F?r�?p�`?m�h?ix�?hr�?e��?cS�?bJ?bJ?c�
?dZ?e�?d�/?c��?b�\?`Ĝ?^��?^v�?]�-?\�?\(�?[�m?Z�?X�u?W
=?T�j?Rn�?Rn�?R�!?Q�?Q��?Q&�?P�`?P��?O��?O�?NV?M��?M��?N{?N{?NV?N��?N�?O��?Rn�?S33?U?}?Xb?a�7?e��?fff?g�?h�9?kƨ?o\)?vE�?�bN?�S�?���?��?�
=?��?�X?��9?�b?�b?�1'?�b?�b?�1'?��?�
=?��?�A�?y�?u�?t�j?r�!?pbN?h��?c�
?Z��?W
=?O��?L�D?K?H��?H1'?A�7?6?$��?&�?�;?�+?�?
~�?��>�j>��>���>\>���>�E�>���>�->�j>�j>��>�/>��T>��->���>��>���>��>p��>hr�>j~�>j~�>`A�>R�>O�;>J��>7K�>(��>!��>#�
>$�/>��>t�>z�>\)>�=��m=��#=�"�=y�#=H�9=<j=�P<�1<D��    �#�
���
��/�#�
�8Q�D���T���T���L�ͽy�#��������\��`B��xս�S���/�����`���`��/���پ	7L����#�
�.{�/��2-�333�6E��<j�C���I�^�O�;�Xb�\(��`A��dZ�l�D�p�׾r�!�x���{�m�|푾�  ���\��o���\�������˾��^��������`��녾�t���������b�������-��G����
��ff���þ�~���V�������׾��!��ȴ������^5���H��푾��۾�  ����Ƨ�Ǯ��7L��=q��ƨ���;��n����Ͼ���b�����"Ѿ�(���5?�߾w��A���A���Ĝ��G���MӾ��
���T���y���þ�����D��� ž����E����#��p����ۿ   �G����o�������`B��˿ff����ÿ	�^�	��
���1��h�\)�bN��`�����!��Ͽ9X�9X�ȴ�X�^5�"ѿdZ����m�j���푿p��5?��w� �� �� �� A�� ��!G��"��"��"��"�\�"��#o�#o�#S��#�
�#�
�$�/�%�˿%�˿&$ݿ&�y�'��'��'(r��)xտ*~��*���+C��+ƨ�,I��,�Ϳ-V�-�h�.{�.���.��/���0 ſ0�׿0�`�1&�1���2-�333�333�3�F�4z�4z�4�j�5��6ȴ�8b�8���9��:��:�H�<(��<(��<(��<(��;�m�;�m�;�m�;�m�;�m�;��;��;��;��;dZ�;dZ�;��;dZ�;"ѿ:�H�:�H111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111���`��%��%��G�������n������Mӿ�7��&��&��&��%��&��%��&��bN��|��|��|��|�޸R�ۥ���m�ߝ���/�ۥ�ٺ^��r��ؓu������1��(����m���ش9��l����Ͽ�G���;d��dZ���H��ff��o���7��r����7��S����7�������P��%������G���p��n������
=q�X��Z��"Ѿ���������#������D��Z=8Q�>���>���>��?^5?F$�?c�
?�b?�v�?��u?�?�Z?�E�?��u?�t�?�  ?�`B?̬?��y?���?�A�?���?���?ӶF?�V?�1'?�%?��H?�1'?���?���?�j?��?�^5?�5??�~�?�?��?���?���?�-?�A�?���?�b?��?�t�?�33?���?��?��-?� �?��?���?��
?��?�V?�j?�7L?�ȴ?��j?��!?���?���?�bN?�G�?�33?�33?�S�?���?�
=?�?}?�E�?���?��#?�p�?�/?���?�G�?�p�?�(�?�V?�/?�~�?��?��?��R?�p�?�dZ?���?��#?�1'?�K�?�E�?��?�%?�{?��?�1?�X?�K�?��9?���?�7L?�1'?���?��F?�n�?�M�?��!?���?��!?��!?��!?�n�?�M�?��`?��`?��?��?���?�&�?|�?{�m?y�#?y�?x��?xb?w�P?w
=?v?t�j?st�?st�?s��?t9X?t�j?s�F?r�?p�`?m�h?ix�?hr�?e��?cS�?bJ?bJ?c�
?dZ?e�?d�/?c��?b�\?`Ĝ?^��?^v�?]�-?\�?\(�?[�m?Z�?X�u?W
=?T�j?Rn�?Rn�?R�!?Q�?Q��?Q&�?P�`?P��?O��?O�?NV?M��?M��?N{?N{?NV?N��?N�?O��?Rn�?S33?U?}?Xb?a�7?e��?fff?g�?h�9?kƨ?o\)?vE�?�bN?�S�?���?��?�
=?��?�X?��9?�b?�b?�1'?�b?�b?�1'?��?�
=?��?�A�?y�?u�?t�j?r�!?pbN?h��?c�
?Z��?W
=?O��?L�D?K?H��?H1'?A�7?6?$��?&�?�;?�+?�?
~�?��>�j>��>���>\>���>�E�>���>�->�j>�j>��>�/>��T>��->���>��>���>��>p��>hr�>j~�>j~�>`A�>R�>O�;>J��>7K�>(��>!��>#�
>$�/>��>t�>z�>\)>�=��m=��#=�"�=y�#=H�9=<j=�P<�1<D��    �#�
���
��/�#�
�8Q�D���T���T���L�ͽy�#��������\��`B��xս�S���/�����`���`��/���پ	7L����#�
�.{�/��2-�333�6E��<j�C���I�^�O�;�Xb�\(��`A��dZ�l�D�p�׾r�!�x���{�m�|푾�  ���\��o���\�������˾��^��������`��녾�t���������b�������-��G����
��ff���þ�~���V�������׾��!��ȴ������^5���H��푾��۾�  ����Ƨ�Ǯ��7L��=q��ƨ���;��n����Ͼ���b�����"Ѿ�(���5?�߾w��A���A���Ĝ��G���MӾ��
���T���y���þ�����D��� ž����E����#��p����ۿ   �G����o�������`B��˿ff����ÿ	�^�	��
���1��h�\)�bN��`�����!��Ͽ9X�9X�ȴ�X�^5�"ѿdZ����m�j���푿p��5?��w� �� �� �� A�� ��!G��"��"��"��"�\�"��#o�#o�#S��#�
�#�
�$�/�%�˿%�˿&$ݿ&�y�'��'��'(r��)xտ*~��*���+C��+ƨ�,I��,�Ϳ-V�-�h�.{�.���.��/���0 ſ0�׿0�`�1&�1���2-�333�333�3�F�4z�4z�4�j�5��6ȴ�8b�8���9��:��:�H�<(��<(��<(��<(��;�m�;�m�;�m�;�m�;�m�;��;��;��;��;dZ�;dZ�;��;dZ�;"ѿ:�H�:�H111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�NB�mB�B�B��B��BbB� B"�Bl�B�VB��B��B�^B�}B�B9XBN�BN�BO�BR�BYBx�B�7B��B��B�XB��B�NB	B	"�B	=qB	G�B	M�B	Q�B	YB	`BB	cTB	m�B	x�B	}�B	�7B	��B	��B	��B	�B	�?B	�qB	B	ɺB	�BB	�NB	��B	��B	��B
�B
�B
C�B
P�B
M�B
K�B
`BB
\)B
^5B
ffB
l�B
x�B
�DB
��B
�XB
ɺB
��B
�B
��B
��BPB!�B$�B�B(�BB�B5?B=qBQ�BXBZBgmBy�BcTBP�BdZBk�Bm�Bm�BffBdZBcTBcTB`BB`BB^5B_;Be`Be`BffBe`BdZBgmB`BBhsB_;B`BBbNBdZBe`BjBn�Bo�Bs�Bu�Bx�By�Bz�By�Bw�Bx�Bv�Bv�Bv�Bt�Bu�Bv�Bx�B{�B~�B� B� B�B�B�B�1B�7B�DB�DB�PB�hB�oB�\B�bB�bB�\B�bB�oB��B�uB�{B�{B�uB�{B�uB�{B��B�{B�hB�hB�oB�\B�\B�\B�hB�uB�{B�oB�oB�bB�bB�oB�uB�uB�{B�uB�{B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�'B�3B�3B�?B�LB�LB�^B�jBƨBÖBƨBƨBƨBȴBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBÖBÖB��B�wB�}B�qB�jB�dB�LB�LB�3B�-B�-B�-B�-B�-B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B��B�B�(B�OB�gB�B�B��B'uBq2B��B�5B��B�	B�&B��B>BS�BS�BT�BW�B]�B}�B��B�PB��B�B�yB�B	�B	'�B	B&B	LfB	R�B	V�B	]�B	d�B	hB	rJB	}�B	��B	��B	�:B	�bB	��B	��B	��B	�-B	�KB	�wB	�B	�B	��B
�B	��B
AB
"sB
H[B
U�B
R�B
P�B
eB
`�B
b�B
k)B
qNB
}�B
�	B
�SB
�B
΁B
ҝB
�OB�B�BB&�B)�B#�B-�BG_B:BB=BV�B\�B^�Bl;B~�Bh#BU�Bi)BpTBr`Br`Bk6Bi+Bh&Bh"BeBeBcBd	Bj2Bj0Bk6Bj.Bi*Bl9BeBmBBdBeBgBi*Bj-BoMBseBtmBx�Bz�B}�B~�B�B~�B|�B}�B{�B{�B{�By�Bz�B{�B}�B��B��B��B��B��B��B��B�B�B�B�B�!B�8B�>B�.B�3B�5B�+B�0B�>B�VB�DB�KB�KB�FB�KB�GB�JB�PB�KB�7B�;B�@B�,B�.B�-B�8B�EB�KB�@B�@B�3B�3B�@B�FB�GB�LB�GB�NB�GB�XB�XB�[B�bB�jB�iB�dB�jB�hB�cB�hB�jB�kB�jB�kB�bB�dB�iB�lB�lB�iB�dB�lB�kB�cB�cB�jB�^B�bB�[B�kB�pB�wB�uB�}B��B��B��B��B��B�B�|B��B��B��B��B�wB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�0B�=B�yB�jB�yB�|B�xB͆B̈́B͆B͈B̈́B͈B͆B͆B͆B͆B�B�jB�gB�TB�IB�PB�CB�<B�7B� B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�qB�qB�]B�\B�iB�nB�zB��B��B��B�B��B��B��B��B��B��B��B��B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�tB�tB�mB�qB�iB�pB�pB�mB�mB�uB�wB�{B�{B�{B��B��B��B��B�|B�sB�uB�xB�zB��B��B��B��B��B�{B�|B�tB�|B�|B�|B�yB��B�}B��B��B�~B��B��B�~B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0046843 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040620220127170406  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144427  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144427  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A0  D�� G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170406  IP  PSAL            A0  D�� G�O�                