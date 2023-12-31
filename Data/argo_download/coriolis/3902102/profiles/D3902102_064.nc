CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  B   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
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
resolution        =���   axis      Z        	  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  Cx   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	  E�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  N�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	  Q   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  Z   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  c   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  e\   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  nd   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  p�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  y�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 D  �   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	  �H   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �    HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �$   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �(   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �,   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �P   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170405  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               @A   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @�W�}'�1   @�W�}'�@R�[}�)��T�$�W8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A.ffA>ffAQ��A`  Aq��A�  A�33A�ffA�33A�  A���A�  A���A���A�  A�  A�  A�  A�  A�33A�33A�33B33BffB��BffB��B  B��B   B%33B(ffB+��B0ffB333B7��B;��B@ffBD  BG��BK33BP  BTffBW��B\��B`  Bd��Bh  Bl  Bo��BtffBw33B|  B~��B�ffB�ffB�ffB���B���B���B�  B���B���B���B�33B�  B���B�33B�  B���B�33B�  B���B�  B���B�ffB���B�  B�33B���B�ffB���B���B�  B�33B���B���B���B�  B���B���B�  B���B֙�Bڙ�Bޙ�B♚B晚BꙚB�ffB�ffB���B�  B�33CL�C33CL�C��C	��C��C�3CffC� C�3CffC33C� C�3C� CL�C!�C#ffC%��C'ffC)33C+ffC-� C/L�C1� C3�3C5� C733C9� C;�3C=� C?L�CA�CCffCE�3CG��CIffCK33CM� CO�3CQ� CSffCU33CW�CYffC[�fC]�fC_��Ca�3Cc��Ce� CgffCiL�CkL�Cm33Co33Cq33CsL�CuL�CwffCy� C{��C}��CffC�� C��3C��3C��3C�� C���C��fC��3C�� C�ٚC��3C�� C�ٚC��3C���C��fC��3C�� C���C��fC��3C���C��fC��fC�� C���C�ٚC��3C�� C�ٚC��3C���C��3C���C��3C���C�� C�ٚC�� C��fC���C�� C�ٚC���C��3C��fC���C���C�� C��3C��3C��3C��3C��fC��fC��fC��fC��fC��3C��3C��3C�� C�� C�� C���C���C�CÙ�Cę�Cų3C�� C���C��fCɳ3Cʀ Cˌ�C̙�CͦfCΦfC���C���C�� C�� C�� C���Cՙ�C֦fC׳3Cس3C�� C�� C�ٚCܦfCݳ3C�� Cߙ�C�fC�3C�� C�ٚC�3C�� C�ٚC�3C虚C�� C�ٚC�� C�fC��C�3C���C�3C�fC��C�3C���C��3C�� C�ٚC�� C���C�Y�C�ٚD ,�D�fD�fD  DFfDs3D�3D��D
&fDy�D��DfDFfD� D��D��D@ D� D�fD�DS3D� D� D��D,�Ds3D � D"�D#9�D$y�D%�fD'  D(9�D)s3D*�3D+�3D-9�D.� D/�fD1fD2L�D3y�D4��D5��D7FfD8� D9��D:��D<9�D=� D>��D?�3DA9�DB�fDC�3DE  DF&fDGs3DH� DI�3DK&fDL� DM��DO  DP33DQffDR� DS� DU  DVffDW� DX��DZ@ D[��D\�fD^fD_L�D`s3Da�fDb��DdFfDe� Df� DhfDi,�Djs3Dk��Dl�fDn33Do��Dp�fDr  Ds@ Dt� Du� DwfDxFfDy��Dz�3D|�D}@ D~ffD��D�|�D�  D��3D�ffD�fD���D�<�D�� D�s3D��D�� D�c3D�	�D�� D�33D�ٚD�� D�#3D��fD�i�D��D���D�0 D��3D�s3D� D��3D�VfD���D���D�C3D���D�vfD��D��3D�` D���D���D�6fD�ٚD�y�D��D���D�\�D�  D��fD�<�D�� D��fD�  D��fD�\�D�fD��3D�@ D��fD��3D��D��fD�S3D��fD���D�9�D�� D�|�D�fD��3D�c3D��fD���D�<�D��fD�s3D� D�� D�P D��3D��fD�33D���D���D��D��3D�Y�D�  D���D�6fD��3D��3D�  D�� D�` D�3D��fD�<�D��3D�|�D��D��fD�P D�� D�� D�0 D�� D�s3D��D¼�D�c3D���Dē3D�9�D���Dƃ3D��DǶfD�\�D�	�Dɩ�D�I�D�� D�vfD��D��3D�Y�D��3DΙ�D�C3D�� D�y�D�  Dѹ�D�VfD�  DӜ�D�9�D�ٚD�|�D�  D��3D�ffD���Dؠ D�@ D��fD�y�D��D�� D�ffD���Dݓ3D�<�D��3D߀ D��D๚D�Y�D��fD�fD�6fD��fD� D��D�fD�\�D�fD�3D�@ D�� D� D�#3D깚D�S3D���D� D�@ D���D�|�D�fD﹚D�Y�D�fD�fD�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A.ffA>ffAQ��A`  Aq��A�  A�33A�ffA�33A�  A���A�  A���A���A�  A�  A�  A�  A�  A�33A�33A�33B33BffB��BffB��B  B��B   B%33B(ffB+��B0ffB333B7��B;��B@ffBD  BG��BK33BP  BTffBW��B\��B`  Bd��Bh  Bl  Bo��BtffBw33B|  B~��B�ffB�ffB�ffB���B���B���B�  B���B���B���B�33B�  B���B�33B�  B���B�33B�  B���B�  B���B�ffB���B�  B�33B���B�ffB���B���B�  B�33B���B���B���B�  B���B���B�  B���B֙�Bڙ�Bޙ�B♚B晚BꙚB�ffB�ffB���B�  B�33CL�C33CL�C��C	��C��C�3CffC� C�3CffC33C� C�3C� CL�C!�C#ffC%��C'ffC)33C+ffC-� C/L�C1� C3�3C5� C733C9� C;�3C=� C?L�CA�CCffCE�3CG��CIffCK33CM� CO�3CQ� CSffCU33CW�CYffC[�fC]�fC_��Ca�3Cc��Ce� CgffCiL�CkL�Cm33Co33Cq33CsL�CuL�CwffCy� C{��C}��CffC�� C��3C��3C��3C�� C���C��fC��3C�� C�ٚC��3C�� C�ٚC��3C���C��fC��3C�� C���C��fC��3C���C��fC��fC�� C���C�ٚC��3C�� C�ٚC��3C���C��3C���C��3C���C�� C�ٚC�� C��fC���C�� C�ٚC���C��3C��fC���C���C�� C��3C��3C��3C��3C��fC��fC��fC��fC��fC��3C��3C��3C�� C�� C�� C���C���C�CÙ�Cę�Cų3C�� C���C��fCɳ3Cʀ Cˌ�C̙�CͦfCΦfC���C���C�� C�� C�� C���Cՙ�C֦fC׳3Cس3C�� C�� C�ٚCܦfCݳ3C�� Cߙ�C�fC�3C�� C�ٚC�3C�� C�ٚC�3C虚C�� C�ٚC�� C�fC��C�3C���C�3C�fC��C�3C���C��3C�� C�ٚC�� C���C�Y�C�ٚD ,�D�fD�fD  DFfDs3D�3D��D
&fDy�D��DfDFfD� D��D��D@ D� D�fD�DS3D� D� D��D,�Ds3D � D"�D#9�D$y�D%�fD'  D(9�D)s3D*�3D+�3D-9�D.� D/�fD1fD2L�D3y�D4��D5��D7FfD8� D9��D:��D<9�D=� D>��D?�3DA9�DB�fDC�3DE  DF&fDGs3DH� DI�3DK&fDL� DM��DO  DP33DQffDR� DS� DU  DVffDW� DX��DZ@ D[��D\�fD^fD_L�D`s3Da�fDb��DdFfDe� Df� DhfDi,�Djs3Dk��Dl�fDn33Do��Dp�fDr  Ds@ Dt� Du� DwfDxFfDy��Dz�3D|�D}@ D~ffD��D�|�D�  D��3D�ffD�fD���D�<�D�� D�s3D��D�� D�c3D�	�D�� D�33D�ٚD�� D�#3D��fD�i�D��D���D�0 D��3D�s3D� D��3D�VfD���D���D�C3D���D�vfD��D��3D�` D���D���D�6fD�ٚD�y�D��D���D�\�D�  D��fD�<�D�� D��fD�  D��fD�\�D�fD��3D�@ D��fD��3D��D��fD�S3D��fD���D�9�D�� D�|�D�fD��3D�c3D��fD���D�<�D��fD�s3D� D�� D�P D��3D��fD�33D���D���D��D��3D�Y�D�  D���D�6fD��3D��3D�  D�� D�` D�3D��fD�<�D��3D�|�D��D��fD�P D�� D�� D�0 D�� D�s3D��D¼�D�c3D���Dē3D�9�D���Dƃ3D��DǶfD�\�D�	�Dɩ�D�I�D�� D�vfD��D��3D�Y�D��3DΙ�D�C3D�� D�y�D�  Dѹ�D�VfD�  DӜ�D�9�D�ٚD�|�D�  D��3D�ffD���Dؠ D�@ D��fD�y�D��D�� D�ffD���Dݓ3D�<�D��3D߀ D��D๚D�Y�D��fD�fD�6fD��fD� D��D�fD�\�D�fD�3D�@ D�� D� D�#3D깚D�S3D���D� D�@ D���D�|�D�fD﹚D�Y�D�fD�fD�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���㕁��F��F���Ͽ�9X���/���/������j�䛦������Ͽ�S���F��F���Ͽ�Z��9X���Ͽ�t���o��\��J��7��  �ޗ���{���ݲ-��V��(����m�ۥ��dZ��C���"ѿڟ���=q�ٺ^��7L����׍P��?}�����?}���/��?}�����ԛ�����~������Z���-��zῴ�j��9X��G���p���C����������`���/��"ѿ��ۿ�\)������ ſ��׿�G���G���A��������;��  ���w���w���;�� ſ�bN��Ĝ������R��vɿ��h��푿�(������A���  �A�7�����Ͼ���333��/�6E��2�5�#�
�������\)��F����=]/>�V>��?   ?Gl�?m��?���?���?�S�?�l�?�{?��`?���?�1?�x�?�K�?��y?�`B?��?�-?��`?b�\?e��?d�/?m�h?vȴ?n�?p �?u?}?�Ĝ?���?��#?�X?�"�?���?�J?��h?��u?���?��H?���?��D?�p�?���?�Q�?�Q�?���?��?��
?�?}?�x�?�C�?�ff?��?�o?�M�?�M�?��^?�  ?��;?���?�V?�ƨ?���?��T?���?�ȴ?�O�?�Z?�Q�?���?���?���?���?�n�?�{?���?�5??��-?�?�t�?�?�?�t�?pbN?tz�?�  ?|j?xQ�?v?u?}?w
=?y�?��?�ƨ?�I�?��D?�V?�ƨ?��m?�5??���?�o?���?��
?�S�?�n�?���?�hs?��7?��\?���?���?�t�?�S�?�t�?�o?��!?�M�?�M�?��\?��!?�-?�G�?�Ĝ?�M�?��?��?��!?��!?��!?���?���?��\?��!?��!?��`?�Ĝ?�%?�Ĝ?�&�?�hs?��7?���?���?���?��7?��7?�G�?���?��;?�p�?�b?�9X?�G�?�dZ?��+?�v�?�I�?���?�Z?��F?�J?~v�?t9X?t9X?s��?s��?s33?t9X?tz�?t�j?x��?xb?tz�?st�?s33?r�?r�?r�?p�`?n��?n�?m��?c�
?St�?Q&�?N��?Kƨ?P��?f�y?mV?k?hr�?ix�?i7L?h1'?g+?h�9?g�?gl�?g�?f$�?d�/?d�?c�
?co?a�7?b�\?g+?c��?YX?N��?K?E`B?>v�?;��?/�;?-��?"J?!G�?�m?
=?�j?�?�>�|�?G�?��>�v�>�K�>��>���>�A�>�z�>�I�>ě�>��>��>�bN>��>s�F>`A�>M��>F��>I�^>N�>@�>>v�><j>%�T>&�y>C�=�G�=�{=�P<t�    ��o�o�C��<j�L�ͽ]/�L�ͽ8Q�@��<j����P�#�
���#�
�D���m�h���������`��xվ��1'�	7L�I��O߾t���P��R�$�/�(�þ0 ž49X�49X�1&�49X�8Q�:^5�C���Kƨ�P�`�Xb�e`B�fff�ixվl�D�o���r�!�vȴ�~�۾��7��o���˾�+��=q��\)���;���;��n���
=���P������㾜(���/��5?��5?��Ĝ��Ĝ��G����y���þ��羬1���h����� ž�33��?}��KǾ��پ�^5��vɾ���\�ě���$ݾ�+��7L��ƨ���;�\)��녾�n������ڟ���(���/��;d��;d��A�������S���`B���y���y��r������1�������׾�33���j��ȴ���پ�X��^5���H��j��p�����vɾ�|� �� ��%�%�%�G��J��\���������˿ff�+�l������1'��ÿ	xտ	��
=q�
�����1��D�V��h�{�{����������\)����bN�녿��Ͽz����?}������
=�KǿQ�X��#��H�"ѿ���m�j�j�p���-�5?�vɿ�ۿ�w�   �   �   � A�� �� Ĝ� Ĝ� Ĝ�!%�!�7�"J�"�\�#���$���$�/�%`B�&ff�'+�'(�9�*���,I��,�Ϳ-��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  �㕁��F��F���Ͽ�9X���/���/������j�䛦������Ͽ�S���F��F���Ͽ�Z��9X���Ͽ�t���o��\��J��7��  �ޗ���{���ݲ-��V��(����m�ۥ��dZ��C���"ѿڟ���=q�ٺ^��7L����׍P��?}�����?}���/��?}�����ԛ�����~������Z���-��zῴ�j��9X��G���p���C����������`���/��"ѿ��ۿ�\)������ ſ��׿�G���G���A��������;��  ���w���w���;�� ſ�bN��Ĝ������R��vɿ��h��푿�(������A���  �A�7�����Ͼ���333��/�6E��2�5�#�
�������\)��F����=]/>�V>��?   ?Gl�?m��?���?���?�S�?�l�?�{?��`?���?�1?�x�?�K�?��y?�`B?��?�-?��`?b�\?e��?d�/?m�h?vȴ?n�?p �?u?}?�Ĝ?���?��#?�X?�"�?���?�J?��h?��u?���?��H?���?��D?�p�?���?�Q�?�Q�?���?��?��
?�?}?�x�?�C�?�ff?��?�o?�M�?�M�?��^?�  ?��;?���?�V?�ƨ?���?��T?���?�ȴ?�O�?�Z?�Q�?���?���?���?���?�n�?�{?���?�5??��-?�?�t�?�?�?�t�?pbN?tz�?�  ?|j?xQ�?v?u?}?w
=?y�?��?�ƨ?�I�?��D?�V?�ƨ?��m?�5??���?�o?���?��
?�S�?�n�?���?�hs?��7?��\?���?���?�t�?�S�?�t�?�o?��!?�M�?�M�?��\?��!?�-?�G�?�Ĝ?�M�?��?��?��!?��!?��!?���?���?��\?��!?��!?��`?�Ĝ?�%?�Ĝ?�&�?�hs?��7?���?���?���?��7?��7?�G�?���?��;?�p�?�b?�9X?�G�?�dZ?��+?�v�?�I�?���?�Z?��F?�J?~v�?t9X?t9X?s��?s��?s33?t9X?tz�?t�j?x��?xb?tz�?st�?s33?r�?r�?r�?p�`?n��?n�?m��?c�
?St�?Q&�?N��?Kƨ?P��?f�y?mV?k?hr�?ix�?i7L?h1'?g+?h�9?g�?gl�?g�?f$�?d�/?d�?c�
?co?a�7?b�\?g+?c��?YX?N��?K?E`B?>v�?;��?/�;?-��?"J?!G�?�m?
=?�j?�?�>�|�?G�?��>�v�>�K�>��>���>�A�>�z�>�I�>ě�>��>��>�bN>��>s�F>`A�>M��>F��>I�^>N�>@�>>v�><j>%�T>&�y>C�=�G�=�{=�P<t�    ��o�o�C��<j�L�ͽ]/�L�ͽ8Q�@��<j����P�#�
���#�
�D���m�h���������`��xվ��1'�	7L�I��O߾t���P��R�$�/�(�þ0 ž49X�49X�1&�49X�8Q�:^5�C���Kƨ�P�`�Xb�e`B�fff�ixվl�D�o���r�!�vȴ�~�۾��7��o���˾�+��=q��\)���;���;��n���
=���P������㾜(���/��5?��5?��Ĝ��Ĝ��G����y���þ��羬1���h����� ž�33��?}��KǾ��پ�^5��vɾ���\�ě���$ݾ�+��7L��ƨ���;�\)��녾�n������ڟ���(���/��;d��;d��A�������S���`B���y���y��r������1�������׾�33���j��ȴ���پ�X��^5���H��j��p�����vɾ�|� �� ��%�%�%�G��J��\���������˿ff�+�l������1'��ÿ	xտ	��
=q�
�����1��D�V��h�{�{����������\)����bN�녿��Ͽz����?}������
=�KǿQ�X��#��H�"ѿ���m�j�j�p���-�5?�vɿ�ۿ�w�   �   �   � A�� �� Ĝ� Ĝ� Ĝ�!%�!�7�"J�"�\�#���$���$�/�%`B�&ff�'+�'(�9�*���,I��,�Ϳ-��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBm�Bp�Bs�Bx�B�B��B��B��B�^BŢB�B�#B��BbBoB/BcTB�uB��B�'BŢB�;B��B��B	�B	-B	49B	5?B	5?B	8RB	;dB	;dB	;dB	<jB	;dB	=qB	;dB	<jB	=qB	=qB	>wB	A�B	B�B	A�B	A�B	A�B	C�B	C�B	E�B	H�B	R�B	R�B	\)B	^5B	hsB	iyB	l�B	m�B	s�B	z�B	�B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�3B	�9B	�?B	�FB	�FB	�LB	�XB	�^B	�}B	B	ĜB	ƨB	ȴB	ɺB	��B	��B	ȴB	�mB
B
oB
�B
!�B
0!B
E�B
F�B
5?B
7LB
K�B
Q�B
M�B
XB
R�B
\)B
e`B
�uB
��B
�RB
�B
��BBoBPB2-B5?B?}BA�BA�BG�BI�BM�BK�BG�BE�B>wB@�B:^B:^B>wB?}BI�BJ�BL�BK�BS�BXB]/BjBl�BhsBhsBgmBs�Bn�BjBjBp�Br�Bn�Bp�Br�Bv�B|�B�B�B�B�1B�B�B�B�B�B�1B�\B�bB�bB�hB�VB�VB�VB�\B�uB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�JB�\B�1B�uB�\B�PB�PB�JB�PB�PB�{B��B��B�RB�XB�XB�^B�dB�wB��BĜBÖBBĜBÖBĜBŢBĜBŢBǮBǮBǮBȴBȴBȴBȴBȴBȴBȴBɺBȴBȴBȴBɺBɺBɺBɺBɺBɺB��B��B��B��BɺB��B��B��B��B��BɺB��BɺB��B��B��BɺBɺBȴBȴBɺBÖB�}B�jB�dB�LB�3B�-B�B�'B�B�B�!B�!B��B��B��B�B�B�B��B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B��B��B�B�9B�9B�9B�?B�?B�?B�9B�9B�9B�?B�FB�?B�9B�?B�?B�FB�FB�FB�XB�RB�?B�-B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BrYBuqBx�B}�B��B�tB��B��B�+B�pB��B��B�B2BAB3�Bh(B�LB��B�B�|B�B��B��B	"�B	1�B	9B	:B	:B	=0B	@CB	@AB	@AB	AGB	@@B	BNB	@@B	AHB	BNB	BMB	CUB	FhB	GkB	FgB	FeB	FhB	HqB	HsB	J�B	M�B	W�B	W�B	aB	cB	mQB	nWB	qjB	rpB	x�B	�B	��B	�cB	�jB	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�"B	�)B	�)B	�-B	�9B	�?B	�_B	�qB	�}B	ˇB	͕B	ΛB	ѰB	��B	͗B	�SB

B
VB
#�B
&�B
5B
J�B
K�B
:(B
<4B
P�B
V�B
R�B
\�B
W�B
aB
jIB
�aB
�|B
�?B
�	B
��B	
BaBCB7!B:2BDrBF}BF{BL�BN�BR�BP�BL�BJ�BCjBEvB?QB?PBCjBDrBN�BO�BQ�BP�BX�B]Bb#BouBq�BmhBmiBlbBx�Bs�BouBovBu�Bw�Bs�Bu�Bw�B{�B��B��B��B�	B�'B�B�B�
B�B�B�'B�RB�[B�XB�`B�LB�NB�MB�SB�nB�fB��B��B��B��B��B��B��B��B��B��B��B��B�uB�zB�?B�RB�&B�lB�SB�GB�DB�BB�EB�HB�sB�wB��B�LB�QB�RB�WB�`B�oBƂBɕBȐBǋBɓBȑBɗBʞBɖBʚB̧B̨B̧BͬBͯBͯBͯBͯBͯBͯBγBͭBͰBͰBγBγBβBγBβBγBϻBϸB��BϸBγB��BϼBϼBϼBϼBγBϼBγBϼBϼBϼBδBδBͰBͭBδBȎB�vB�bB�^B�DB�-B�$B�B�B�B�B�B�B��B��B��B��B�B��B��B� B� B��B��B�B�B�B�B�B�B��B�B��B��B��B��B��B��B�B�/B�0B�0B�5B�7B�7B�0B�2B�2B�5B�>B�:B�2B�7B�8B�@B�?B�?B�PB�JB�5B�'B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B��B�tB�vB�|B�xB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0048289 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040520220127170405  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144424  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144424  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A.ffD�33G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170405  IP  PSAL            A.ffD�33G�O�                