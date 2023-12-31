CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  V   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-08-04T21:26:33Z creation; 2023-08-05T07:55:33Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_030d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8    FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8(   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    88   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8H   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8P   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9    	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9,   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    90   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     94   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9T   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9t   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	X  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  D0   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	X  F�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	X  R8   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  [�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  g@   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  p�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  |H   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  �P   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
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
_FillValue                  ,  �0Argo profile    3.1 1.2 19500101000000  20190804212633  20230805075533  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               ,A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��w'�}'1   @��w���@@R���p������DZ1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 7.6 dbar]                                                      A   A��A   A1��A@  AQ��A^ffAq��A���A���A���A���A�  A���A�ffA�  A�33A�ffA�ffA�ffA�33A���A�33A�33A�ffB��B��B
��B33B  B��B33B   B$  B(  B+��B/��B3��B733B;33B?33BC33BG33BK33BO33BS��BW33B[��B_��Bd  Bh  BlffBo33Bs33Bw��B|  B�33B���B���B�  B�  B�ffB���B�  B�33B���B���B�33B���B���B�  B�  B�33B���B���B�  B���B���B���B�  B�33B���B���B�  B�  B�33B���B���B���B�  B�33B�33B�ffB�  B�ffB�  B֙�B���B�33B���B�ffB���B�33B�  B���B���B�ffC��C��C� C� C	� C� C��C��C��C�3CL�CL�CffCffC� C� C!��C#L�C%L�C'ffC)� C+�3C-ffC/� C1��C3L�C5ffC7� C9L�C;ffC=��C?L�CA� CC�3CEffCG33CIffCK��CML�CO� CQ�3CS� CU33CW� CY��C[ffC]33C_ffCa��CcffCe33CgffCi��CkL�Cm�CoffCq��CsL�Cu�CwL�Cy� C{L�C}� C�3C�� C���C��3C���C��fC�� C�ٚC�� C���C�� C�ٚC�� C��fC���C��3C���C�� C��fC���C��3C�ٚC�� C��fC���C�� C�ٚC�� C��3C��fC���C�� C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC��fC�ٚC���C�� C��3C��fC���C���C���C��3C��3C��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C��3C��fC¦fC���C���C�� C�� C�� C�� C�� C�� C�� C�� C���C�ٚCϦfCг3C�� Cҙ�CӦfC�� Cՙ�C֦fC�� C���C�ٚCڳ3Cی�CܦfCݳ3C�� C�ٚC�3CᙚC�3C�� C�fC�� C�ٚC�� C虚C�3C���C�3C왚C� C�3C���C�� C�3C�C��C�� C��fC�ٚC�� C��3C���C�L�C��fD L�D� D��D�3D&fDs3D�fD��D
,�DffD� D��DY�D��D�3D3DS3D�3D��D3D9�D� D��DfD@ D�fD ��D"  D#33D$�fD%��D&��D(@ D)y�D*��D,  D-9�D.s3D/�3D0�3D29�D3� D4�fD5�3D7@ D8��D9� D:�3D<&fD=� D>�3D@�DAL�DB�fDC�fDE  DF@ DG� DH� DJ  DKFfDL�fDM�fDO�DPL�DQy�DR��DT  DU@ DV� DW�fDYfDZL�D[� D\�3D^  D_9�D`s3Da�fDcfDd@ De�fDf��Dg��DiL�Dj�fDk� Dl��Dn33Dol�Dp��Dq��Ds,�Dts3Du�3Dv�3Dx@ Dy��Dz� D{��D}9�D~y�D� D��3D�&fD�� D�VfD�  D���D�6fD�� D�|�D��D���D�VfD�3D��3D�@ D��3D�vfD��D���D�S3D���D��3D�@ D�ٚD�vfD�3D��3D�P D�� D�� D�0 D��3D�y�D�#3D��3D�` D���D���D�9�D�ٚD�� D�#3D���D�VfD�  D���D�6fD��3D�� D��D���D�` D�  D��fD�9�D�� D�y�D�  D��fD�` D���D���D�9�D�ٚD�y�D��D���D�Y�D�fD��3D�9�D��fD�|�D�&fD��fD�ffD�fD��3D�FfD��D�|�D�3D��fD�\�D�3D���D�9�D�� D�|�D��D��fD�VfD���D���D�<�D��fD�y�D�#3D���D�S3D���D��fD�@ D�ٚD�vfD�3D³3D�S3D��fDę�D�9�D�ٚD�|�D�  D��fD�Y�D���Dɣ3D�@ D��3D�y�D�3D̼�D�Y�D��fDΣ3D�C3D��3DЃ3D�#3D��3D�c3D�3DӖfD�9�D�� D�vfD� Dּ�D�ffD�3Dأ3D�C3D��3Dڃ3D�&fD�� D�Y�D�  DݦfD�C3D�� D߀ D�  D��3D�ffD���D�3D�9�D��3D�|�D�fD�� D�\�D��fD�fD�6fD��3D�3D�#3D��3D�VfD���D��D�C3D��D� D��D�fD�S3D�� D� D�<�D�� D�p D� D���D�VfD���D��3D�<�D���D�|�D�  D��3D�\�D���D�� D�@ D��f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A   A��A   A1��A@  AQ��A^ffAq��A���A���A���A���A�  A���A�ffA�  A�33A�ffA�ffA�ffA�33A���A�33A�33A�ffB��B��B
��B33B  B��B33B   B$  B(  B+��B/��B3��B733B;33B?33BC33BG33BK33BO33BS��BW33B[��B_��Bd  Bh  BlffBo33Bs33Bw��B|  B�33B���B���B�  B�  B�ffB���B�  B�33B���B���B�33B���B���B�  B�  B�33B���B���B�  B���B���B���B�  B�33B���B���B�  B�  B�33B���B���B���B�  B�33B�33B�ffB�  B�ffB�  B֙�B���B�33B���B�ffB���B�33B�  B���B���B�ffC��C��C� C� C	� C� C��C��C��C�3CL�CL�CffCffC� C� C!��C#L�C%L�C'ffC)� C+�3C-ffC/� C1��C3L�C5ffC7� C9L�C;ffC=��C?L�CA� CC�3CEffCG33CIffCK��CML�CO� CQ�3CS� CU33CW� CY��C[ffC]33C_ffCa��CcffCe33CgffCi��CkL�Cm�CoffCq��CsL�Cu�CwL�Cy� C{L�C}� C�3C�� C���C��3C���C��fC�� C�ٚC�� C���C�� C�ٚC�� C��fC���C��3C���C�� C��fC���C��3C�ٚC�� C��fC���C�� C�ٚC�� C��3C��fC���C�� C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC��fC�ٚC���C�� C��3C��fC���C���C���C��3C��3C��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C��3C��fC¦fC���C���C�� C�� C�� C�� C�� C�� C�� C�� C���C�ٚCϦfCг3C�� Cҙ�CӦfC�� Cՙ�C֦fC�� C���C�ٚCڳ3Cی�CܦfCݳ3C�� C�ٚC�3CᙚC�3C�� C�fC�� C�ٚC�� C虚C�3C���C�3C왚C� C�3C���C�� C�3C�C��C�� C��fC�ٚC�� C��3C���C�L�C��fD L�D� D��D�3D&fDs3D�fD��D
,�DffD� D��DY�D��D�3D3DS3D�3D��D3D9�D� D��DfD@ D�fD ��D"  D#33D$�fD%��D&��D(@ D)y�D*��D,  D-9�D.s3D/�3D0�3D29�D3� D4�fD5�3D7@ D8��D9� D:�3D<&fD=� D>�3D@�DAL�DB�fDC�fDE  DF@ DG� DH� DJ  DKFfDL�fDM�fDO�DPL�DQy�DR��DT  DU@ DV� DW�fDYfDZL�D[� D\�3D^  D_9�D`s3Da�fDcfDd@ De�fDf��Dg��DiL�Dj�fDk� Dl��Dn33Dol�Dp��Dq��Ds,�Dts3Du�3Dv�3Dx@ Dy��Dz� D{��D}9�D~y�D� D��3D�&fD�� D�VfD�  D���D�6fD�� D�|�D��D���D�VfD�3D��3D�@ D��3D�vfD��D���D�S3D���D��3D�@ D�ٚD�vfD�3D��3D�P D�� D�� D�0 D��3D�y�D�#3D��3D�` D���D���D�9�D�ٚD�� D�#3D���D�VfD�  D���D�6fD��3D�� D��D���D�` D�  D��fD�9�D�� D�y�D�  D��fD�` D���D���D�9�D�ٚD�y�D��D���D�Y�D�fD��3D�9�D��fD�|�D�&fD��fD�ffD�fD��3D�FfD��D�|�D�3D��fD�\�D�3D���D�9�D�� D�|�D��D��fD�VfD���D���D�<�D��fD�y�D�#3D���D�S3D���D��fD�@ D�ٚD�vfD�3D³3D�S3D��fDę�D�9�D�ٚD�|�D�  D��fD�Y�D���Dɣ3D�@ D��3D�y�D�3D̼�D�Y�D��fDΣ3D�C3D��3DЃ3D�#3D��3D�c3D�3DӖfD�9�D�� D�vfD� Dּ�D�ffD�3Dأ3D�C3D��3Dڃ3D�&fD�� D�Y�D�  DݦfD�C3D�� D߀ D�  D��3D�ffD���D�3D�9�D��3D�|�D�fD�� D�\�D��fD�fD�6fD��3D�3D�#3D��3D�VfD���D��D�C3D��D� D��D�fD�S3D�� D� D�<�D�� D�p D� D���D�VfD���D��3D�<�D���D�|�D�  D��3D�\�D���D�� D�@ D��f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�n�@�S�@��y@�`B@�V@��@R��@-?}@~�@7�;@s�
@vv�@L��@�w?��?�?��?�x�?s33?�(�?�^5?�z�?���@Ĝ@(��@.$�@-�T@,�/@)7L@(  @*��@+�@,Z@-��@-�-@-�-@-�-@-@-@-p�@.ff@/;d@/�@0Q�@0r�@0A�@0A�@0��@0Q�@0 �@/�w@/��@/��@.�y@.E�@-@-�@.V@.ȴ@.�+@.{@-p�@-p�@-�T@-p�@-`B@-�T@/�w@/��@/l�@/l�@/\)@.��@.ȴ@-@,��@,�j@,��@,�@,�j@,��@,�/@-�@-�@,I�@+��@+S�@+@*M�@)hs@)%@(��@(Ĝ@*�@.v�@/
=@.�@-�h@,1@+"�@*n�@(��@(b@&V@$I�@#o@!��@|�@v�@�j@@=q@l�@"�@"�@n�@G�@�P@@�R@�@S�@
�\@�9@V@��@��?�o?��?��?�P?���?���?�1'?͑h?�x�?Ƨ�?�?�dZ?�I�?�V?�S�?��?�?}?䛦?���?�!?��?�\)?ؓu?�j?ش9?�E�?�E�?ӶF?�hs?Ͼw?���?�?���?�Q�?���?���?��y?���?�?���?�o?���?��?�+?�?}?���?��`?���?��h?��H?���?��^?�?�J?�bN?��h?�C�?��^?��F?���?{�m?wK�?s33?p �?s33?n��?mO�?l�D?gl�?`A�?d�/?`�?\(�?["�?[dZ?g�?kƨ?t�j?��7?��/?��j?�`B?��`?w��?rn�?qhs?qhs?r-?o�;?`�?Tz�?Q��?N��?PbN?Qhs?R�!?F��?7�P?6?5?}?4z�?4��?6�+?8�u?7�P?2�?/\)?-��?,��?-��?0�`?3t�?5?}?:��?8��?4�j?333?1hs?/�;?*~�?,I�?=/?J~�?O�;?P��?Q&�?O\)?O\)?O�?N��?J~�?@Ĝ?:^5?1�?1��?0 �?.V?,1?'l�?5??��?
=?-?r�?��?$�?�T?`B?��?�
?o?S�? Ĝ>��>�!>�>��>��
>���>�G�>߾w>�V>��H>�!>� �>���>��>���>ܬ>���>���>�ȴ>�->�S�>���>�>���>�I�>�I�>��>��u>�ƨ>��>q��>aG�>J��>:^5>,1>�->��>�>bN>	7L>1'>$�>+=��=>J>J>V>	7L>�>%=�x�=�
==�9X=��w=�+=Y�=<j=�P<���<49X;D����o���
�u��t���9X��`B�o��P�#�
�49X�L�ͽ]/�q����%��+��O߽�hs���㽣�
��{��^5�Ƨ���`��G���h�����1'�I��t���u��R�&�y�)��+�/��2-�6E��8Q�;dZ�>vɾC���E�˾E�˾F��I�^�M��Q녾S�ϾT���W
=�Y��^5?�aG��dZ�hr��k��m�h�p�׾t�j�vȴ�x���|푾����J�������˾�+��7L��C����;���\)��녾��Ͼ�����+��b�������������㾜����-���R���w��A���G���MӾ��
��Z��`B���T���y���xվ���V������ ž�&龲-���F����ȴ��Q쾹�#���m��p���vɾ�|�����J�ě������$ݾ�1'��=q��C���������\)��녾��Ͼ��׍P�ؓu����ܬ�޸R�߾w��A�����������
��ff����þ����1��V��������33���j��KǾ�Q��X��^5��dZ��p����ۿ A�� Ĝ��7�J�Mӿo���������T��T�����	7L�	xտ	��
�������ƨ�1��D�O߿{�V����;��`�-�33��F�9X����?}�����E���P�b��u����X������������X�X��#��#�^5�dZ���/�/�/�p���5?��R�;d�|�|�   � �� �� Ĝ�!%�!%� Ĝ�#S��$���%��%�˿%`B�%�˿%�T�&ff�&�y�'+�'��'(1'�(r��(�9�(�ÿ)7L�)7L�)xտ)��)��*=q�*=q�*=q�*=q�*���+�+C��+��,I��-V�-O߿-�h1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�n�@�S�@��y@�`B@�V@��@R��@-?}@~�@7�;@s�
@vv�@L��@�w?��?�?��?�x�?s33?�(�?�^5?�z�?���@Ĝ@(��@.$�@-�T@,�/@)7L@(  @*��@+�@,Z@-��@-�-@-�-@-�-@-@-@-p�@.ff@/;d@/�@0Q�@0r�@0A�@0A�@0��@0Q�@0 �@/�w@/��@/��@.�y@.E�@-@-�@.V@.ȴ@.�+@.{@-p�@-p�@-�T@-p�@-`B@-�T@/�w@/��@/l�@/l�@/\)@.��@.ȴ@-@,��@,�j@,��@,�@,�j@,��@,�/@-�@-�@,I�@+��@+S�@+@*M�@)hs@)%@(��@(Ĝ@*�@.v�@/
=@.�@-�h@,1@+"�@*n�@(��@(b@&V@$I�@#o@!��@|�@v�@�j@@=q@l�@"�@"�@n�@G�@�P@@�R@�@S�@
�\@�9@V@��@��?�o?��?��?�P?���?���?�1'?͑h?�x�?Ƨ�?�?�dZ?�I�?�V?�S�?��?�?}?䛦?���?�!?��?�\)?ؓu?�j?ش9?�E�?�E�?ӶF?�hs?Ͼw?���?�?���?�Q�?���?���?��y?���?�?���?�o?���?��?�+?�?}?���?��`?���?��h?��H?���?��^?�?�J?�bN?��h?�C�?��^?��F?���?{�m?wK�?s33?p �?s33?n��?mO�?l�D?gl�?`A�?d�/?`�?\(�?["�?[dZ?g�?kƨ?t�j?��7?��/?��j?�`B?��`?w��?rn�?qhs?qhs?r-?o�;?`�?Tz�?Q��?N��?PbN?Qhs?R�!?F��?7�P?6?5?}?4z�?4��?6�+?8�u?7�P?2�?/\)?-��?,��?-��?0�`?3t�?5?}?:��?8��?4�j?333?1hs?/�;?*~�?,I�?=/?J~�?O�;?P��?Q&�?O\)?O\)?O�?N��?J~�?@Ĝ?:^5?1�?1��?0 �?.V?,1?'l�?5??��?
=?-?r�?��?$�?�T?`B?��?�
?o?S�? Ĝ>��>�!>�>��>��
>���>�G�>߾w>�V>��H>�!>� �>���>��>���>ܬ>���>���>�ȴ>�->�S�>���>�>���>�I�>�I�>��>��u>�ƨ>��>q��>aG�>J��>:^5>,1>�->��>�>bN>	7L>1'>$�>+=��=>J>J>V>	7L>�>%=�x�=�
==�9X=��w=�+=Y�=<j=�P<���<49X;D����o���
�u��t���9X��`B�o��P�#�
�49X�L�ͽ]/�q����%��+��O߽�hs���㽣�
��{��^5�Ƨ���`��G���h�����1'�I��t���u��R�&�y�)��+�/��2-�6E��8Q�;dZ�>vɾC���E�˾E�˾F��I�^�M��Q녾S�ϾT���W
=�Y��^5?�aG��dZ�hr��k��m�h�p�׾t�j�vȴ�x���|푾����J�������˾�+��7L��C����;���\)��녾��Ͼ�����+��b�������������㾜����-���R���w��A���G���MӾ��
��Z��`B���T���y���xվ���V������ ž�&龲-���F����ȴ��Q쾹�#���m��p���vɾ�|�����J�ě������$ݾ�1'��=q��C���������\)��녾��Ͼ��׍P�ؓu����ܬ�޸R�߾w��A�����������
��ff����þ����1��V��������33���j��KǾ�Q��X��^5��dZ��p����ۿ A�� Ĝ��7�J�Mӿo���������T��T�����	7L�	xտ	��
�������ƨ�1��D�O߿{�V����;��`�-�33��F�9X����?}�����E���P�b��u����X������������X�X��#��#�^5�dZ���/�/�/�p���5?��R�;d�|�|�   � �� �� Ĝ�!%�!%� Ĝ�#S��$���%��%�˿%`B�%�˿%�T�&ff�&�y�'+�'��'(1'�(r��(�9�(�ÿ)7L�)7L�)xտ)��)��*=q�*=q�*=q�*=q�*���+�+C��+��,I��-V�-O߿-�h1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B�B��BĜB�BBD�B�bB0!Bu�B��B
oB
_;B
ffB
9XB
33B
�B
!�B
#�B
33B
_;B
��B
��B�B=qB�JB��B��B�!B�B�'B�XB�qB�}BÖBĜBŢBȴBǮBǮBȴB��B��B��B��B�B�
B�
B�B�B�B�B�B�B�B�B�
B�B�B�B�)B�B�B�#B�#B�/B�/B�5B�TB�NB�NB�TB�NB�TB�TB�NB�HB�NB�HB�HB�NB�NB�NB�NB�NB�TB�HB�HB�BB�BB�BB�BB�BB�HB�)B�yB�B�B�B�B�B�B�B�B�yB�mB�fB�fB�ZB�ZB�HB�BB�5B�)B�B�B�B�
B��B��B��B��B�
B�B��B��B��BĜBÖB�jB�jB�XB�?B�-B�B��B��B��B��B��B�B�B�^B��BǮBǮBȴBȴBȴBȴBĜBƨBǮBǮBƨBĜBÖBBB�qB�qB�9B�'B�9B�dBƨBǮBÖBĜBĜBB�wB�jB�dB�9B�9B�'B�!B�-B�-B�9B�!B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�'B�FB�?B�qB�qB�}B�}B�dB�dB�dB�dB�dB�^B�XB�?B�9B�3B�FB�?B�?B�LB�'B�!B�B�!B�!B�'B�'B�-B�-B�'B�-B�-B�3B�9B�FB�FB�^B�XB�XB�RB�RB�RB�RB�^B�qB��BƨBƨBƨB��BǮBƨBƨBŢB��BB��B�wB�}B�wB�qB�dB�dB�RB�FB�FB�3B�3B�9B�-B�-B�-B�9B�3B�3B�-B�-B�!B�!B�!B�!B�B�B�!B�?B�3B�?B�9B�LB�FB�FB�9B�-B�'B�!B�!B�B�B�B�B�B�B�3B�-B�'B�'B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B��B�B��BĜB�BBD�B�bB0!Bu�B��B
oB
_;B
ffB
9XB
33B
�B
!�B
#�B
33B
_;B
��B
��B�B=qB�JB��B��B�!B�B�'B�XB�qB�}BÖBĜBŢBȴBǮBǮBȴB��B��B��B��B�B�
B�
B�B�B�B�B�B�B�B�B�
B�B�B�B�)B�B�B�#B�#B�/B�/B�5B�TB�NB�NB�TB�NB�TB�TB�NB�HB�NB�HB�HB�NB�NB�NB�NB�NB�TB�HB�HB�BB�BB�BB�BB�BB�HB�)B�yB�B�B�B�B�B�B�B�B�yB�mB�fB�fB�ZB�ZB�HB�BB�5B�)B�B�B�B�
B��B��B��B��B�
B�B��B��B��BĜBÖB�jB�jB�XB�?B�-B�B��B��B��B��B��B�B�B�^B��BǮBǮBȴBȴBȴBȴBĜBƨBǮBǮBƨBĜBÖBBB�qB�qB�9B�'B�9B�dBƨBǮBÖBĜBĜBB�wB�jB�dB�9B�9B�'B�!B�-B�-B�9B�!B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�'B�FB�?B�qB�qB�}B�}B�dB�dB�dB�dB�dB�^B�XB�?B�9B�3B�FB�?B�?B�LB�'B�!B�B�!B�!B�'B�'B�-B�-B�'B�-B�-B�3B�9B�FB�FB�^B�XB�XB�RB�RB�RB�RB�^B�qB��BƨBƨBƨB��BǮBƨBƨBŢB��BB��B�wB�}B�wB�qB�dB�dB�RB�FB�FB�3B�3B�9B�-B�-B�-B�9B�3B�3B�-B�-B�!B�!B�!B�!B�B�B�!B�?B�3B�?B�9B�LB�FB�FB�9B�-B�'B�!B�!B�B�B�B�B�B�B�3B�-B�'B�'B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755332023080507553320230805075533  IF  ARFMCODA030d                                                                20190804212633                      G�O�G�O�G�O�                IF  ARGQCOQC4.3                                                                 20190804212641  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.3                                                                 20190804212641  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology 20190821161257  IP  PSAL            A   D��fG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A   D��fG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172539  IP  PSAL            A   D��fG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075533  IP  PSAL            A   D��fG�O�                