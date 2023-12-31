CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  K   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:50Z creation; 2021-06-07T15:42:39Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        	,  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	,  FP   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  O|   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	,  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  Z�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  d    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  fl   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  o�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  q�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  {   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  �<   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 L  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	,  �    HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �    HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �,   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �\   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �\   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �\   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �\Argo profile    3.1 1.2 19500101000000  20190620115650  20210607154239  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               !A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @؉v�8�1   @؉v�8�@TN�`���@+����8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A0  A@  AQ��A`  Ap  A���A�  A�  A�33A�33A���A�33A�  A�  A�  A���A�  A�ffA���A���A���B   B��B  B  B��B��B  B��B��B$  B)33B,  B/33B4  B8��B;��B>��BC33BH  BK33BO��BTffBX  B[33B`  Bd��Bh  Bk33BnffBs��Bx��B|  B33B�ffB���B�ffB�33B�  B���B�ffB�33B���B�ffB�33B�  B���B�ffB�33B���B���B�ffB�33B�  B���B���B�33B���B���B�ffB�33B���B���B�ffB�33B���B�B�ffB�33B�  Bʙ�B�33B���B�ffB�33B�  B���B�ffB�33B�  B�B�ffB�ffB�33C� CffCffCL�C	33C�C� C��C�3C�3C��C� CffCL�CL�C33C!�C#  C%ffC'�fC)��C+��C-��C/��C1�3C3� C5L�C7��C9� C;ffC=L�C?L�CA��CC��CE33CG33CIL�CKffCM� CO� CQ��CSL�CUffCW��CYL�C[ffC]��C_L�CaffCc��CeL�Cg� Ci�3Ck� Cm33CoffCq��CsffCu33Cw� Cy��C{��C}ffC33C�� C��3C�ٚC�� C��3C��fC���C�� C��3C��fC�ٚC���C�� C�� C��3C��3C��fC��fC��fC��fC��fC���C���C���C��fC��fC��fC��3C��3C�� C�� C���C�ٚC��fC��3C�� C���C���C��fC�� C���C�ٚC��3C�� C���C��fC�� C�ٚC��3C�� C��3C�� C�� C�� C���C�ٚC��fC��3C�� C���C�ٚC��3C�� C�ٚC��fC��3C���CæfCĦfC�� Cƙ�Cǳ3C���CɦfC�� C���C̳3C͌�CΦfC�� CЙ�Cѳ3C���CӦfC�� C�ٚC�� Cי�Cس3C���CڦfC�� C���Cݳ3Cތ�CߦfC�� CᙚC�fC�� C�fC�3C���C�fC�� C�ٚC�3C��C�fC�� CC�3C�� C�fC�� C�ٚC��3C�� C���C��3C�� C���C�� C��D @ Ds3D��DfD9�Ds3D�3D��D
,�Dl�D�3D  DFfDs3D�fD� D  DffD�3D  DL�Dy�D� D�3D@ Dy�D ��D!�fD#  D$ffD%�fD&��D(33D)s3D*��D,fD-9�D.ffD/�3D1  D2L�D3� D4�3D6  D7L�D8y�D9�fD:�fD<&fD=l�D>��D?�3DA9�DB�fDC��DD��DF&fDGs3DH��DJfDK@ DLy�DM� DO  DP,�DQy�DR�fDS��DU,�DVffDW� DX� DZ9�D[� D\�fD]��D_9�D`��Da�fDc�Dd9�Del�Df� Dg��Di33Djs3Dk�3Dl�3Dn33Doy�Dp�fDq�3Ds@ Dtl�Du� Dv�3Dx@ Dyl�Dz�fD{��D}L�D~�fD� D�� D�  D���D�\�D���D�� D�C3D��3D��fD��D�� D�S3D���D��3D�@ D�� D��3D�#3D��fD�\�D��3D��fD�<�D�ٚD�vfD�  D�� D�\�D���D���D�9�D���D�|�D�  D��3D�Y�D�  D��3D�<�D��3D�|�D�)�D��fD�c3D�3D��3D�C3D��3D��3D�#3D��fD�ffD���D�� D�6fD�� D�|�D�fD��3D�S3D��3D��fD�9�D�� D��3D��D��3D�Y�D�3D�� D�<�D�� D�� D�#3D��fD�Y�D���D�� D�C3D�ٚD�� D�)�D��3D�c3D�  D�� D�<�D��3D�|�D� D���D�VfD���D��fD�C3D�� D�� D�#3D��3D�ffD���D��fD�<�D�ٚD�y�D��D¼�D�` D��3Dę�D�@ D�ٚD�vfD�3DǼ�D�i�D�	�DɦfD�FfD��fDˉ�D�,�D̼�D�S3D���DΜ�D�C3D��D�|�D� DѼ�D�Y�D�  DӦfD�@ D���D�|�D��Dּ�D�c3D�fDأ3D�@ D�� Dڀ D�#3D۶fD�VfD��fDݙ�D�<�D�� D߃3D��D�3D�Y�D�fD�3D�@ D�� D� D�#3D�ɚD�` D��fD�fD�C3D���D�y�D�fD�fD�Y�D���D� D�9�D��3D�y�D�  D��3D�i�D�  D�fD�<�D��3D�|�D�fD�� D�ffD�  D���D�<�D���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A0  A@  AQ��A`  Ap  A���A�  A�  A�33A�33A���A�33A�  A�  A�  A���A�  A�ffA���A���A���B   B��B  B  B��B��B  B��B��B$  B)33B,  B/33B4  B8��B;��B>��BC33BH  BK33BO��BTffBX  B[33B`  Bd��Bh  Bk33BnffBs��Bx��B|  B33B�ffB���B�ffB�33B�  B���B�ffB�33B���B�ffB�33B�  B���B�ffB�33B���B���B�ffB�33B�  B���B���B�33B���B���B�ffB�33B���B���B�ffB�33B���B�B�ffB�33B�  Bʙ�B�33B���B�ffB�33B�  B���B�ffB�33B�  B�B�ffB�ffB�33C� CffCffCL�C	33C�C� C��C�3C�3C��C� CffCL�CL�C33C!�C#  C%ffC'�fC)��C+��C-��C/��C1�3C3� C5L�C7��C9� C;ffC=L�C?L�CA��CC��CE33CG33CIL�CKffCM� CO� CQ��CSL�CUffCW��CYL�C[ffC]��C_L�CaffCc��CeL�Cg� Ci�3Ck� Cm33CoffCq��CsffCu33Cw� Cy��C{��C}ffC33C�� C��3C�ٚC�� C��3C��fC���C�� C��3C��fC�ٚC���C�� C�� C��3C��3C��fC��fC��fC��fC��fC���C���C���C��fC��fC��fC��3C��3C�� C�� C���C�ٚC��fC��3C�� C���C���C��fC�� C���C�ٚC��3C�� C���C��fC�� C�ٚC��3C�� C��3C�� C�� C�� C���C�ٚC��fC��3C�� C���C�ٚC��3C�� C�ٚC��fC��3C���CæfCĦfC�� Cƙ�Cǳ3C���CɦfC�� C���C̳3C͌�CΦfC�� CЙ�Cѳ3C���CӦfC�� C�ٚC�� Cי�Cس3C���CڦfC�� C���Cݳ3Cތ�CߦfC�� CᙚC�fC�� C�fC�3C���C�fC�� C�ٚC�3C��C�fC�� CC�3C�� C�fC�� C�ٚC��3C�� C���C��3C�� C���C�� C��D @ Ds3D��DfD9�Ds3D�3D��D
,�Dl�D�3D  DFfDs3D�fD� D  DffD�3D  DL�Dy�D� D�3D@ Dy�D ��D!�fD#  D$ffD%�fD&��D(33D)s3D*��D,fD-9�D.ffD/�3D1  D2L�D3� D4�3D6  D7L�D8y�D9�fD:�fD<&fD=l�D>��D?�3DA9�DB�fDC��DD��DF&fDGs3DH��DJfDK@ DLy�DM� DO  DP,�DQy�DR�fDS��DU,�DVffDW� DX� DZ9�D[� D\�fD]��D_9�D`��Da�fDc�Dd9�Del�Df� Dg��Di33Djs3Dk�3Dl�3Dn33Doy�Dp�fDq�3Ds@ Dtl�Du� Dv�3Dx@ Dyl�Dz�fD{��D}L�D~�fD� D�� D�  D���D�\�D���D�� D�C3D��3D��fD��D�� D�S3D���D��3D�@ D�� D��3D�#3D��fD�\�D��3D��fD�<�D�ٚD�vfD�  D�� D�\�D���D���D�9�D���D�|�D�  D��3D�Y�D�  D��3D�<�D��3D�|�D�)�D��fD�c3D�3D��3D�C3D��3D��3D�#3D��fD�ffD���D�� D�6fD�� D�|�D�fD��3D�S3D��3D��fD�9�D�� D��3D��D��3D�Y�D�3D�� D�<�D�� D�� D�#3D��fD�Y�D���D�� D�C3D�ٚD�� D�)�D��3D�c3D�  D�� D�<�D��3D�|�D� D���D�VfD���D��fD�C3D�� D�� D�#3D��3D�ffD���D��fD�<�D�ٚD�y�D��D¼�D�` D��3Dę�D�@ D�ٚD�vfD�3DǼ�D�i�D�	�DɦfD�FfD��fDˉ�D�,�D̼�D�S3D���DΜ�D�C3D��D�|�D� DѼ�D�Y�D�  DӦfD�@ D���D�|�D��Dּ�D�c3D�fDأ3D�@ D�� Dڀ D�#3D۶fD�VfD��fDݙ�D�<�D�� D߃3D��D�3D�Y�D�fD�3D�@ D�� D� D�#3D�ɚD�` D��fD�fD�C3D���D�y�D�fD�fD�Y�D���D� D�9�D��3D�y�D�  D��3D�i�D�  D�fD�<�D��3D�|�D�fD�� D�ffD�  D���D�<�D���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��������˿��/��\)���������P��`B���Ͽ~�ۿq녿bMӿQ���D���1hs�z�'�>Ǯ?�Ĝ?�=q?�$�?�bN?���?m��?<�?,�D?(r�?^5?��?J>�ff>��>�&�>�V>�|�>���>���>�X>�  >և+>ݲ->�t�>�33>�p�>�E�?o?�?�?$�?#S�?"J?�?%�?&��?-V?6�+?T9X?]p�?ix�?mV?mV?o��?st�?{��?�  ?�+?��R?�-?�?�v�?��?���?�o?��
?��/?�$�?���?�dZ?�bN?ӕ�?��/?�?}?��/?��?��?���?�J?�t�?�1'?��#@��@Ĝ@�j@V@/@&v�@-?}@,�/@0�u@:~�@<(�@<�j@@�u@AG�@@��@?�P@>E�@;C�@:�@9�@7|�@4I�@3"�@1hs@/��@/+@/��@/;d@-�@,�j@+t�@'�@$I�@$1@$j@%�-@&�R@'�P@(r�@*M�@+o@-`B@,�j@*-@)%@&�y@%�@&$�@%�-@%�@$��@#��@#C�@#33@!��@!X@ ��@ Q�@   @K�@��@ff@�@@@@ȴ@�@�y@ff@�y@v�@$�@��@V@?}@?}@�@�@/@O�@O�@?}@�/@z�@1@��@��@&�@�`@Q�@  @1'@�@��@��@�w@��@5?@�@"�@�\@��@��@  @  @  @A�@Q�@Q�@Q�@Q�@Q�@�@Ĝ@�`@&�@7L@7L@�@��@�9@��@  @b@A�@bN@Q�@A�@r�@r�@r�@ �@�;@�P@;d@�y@��@�+@$�@�T@��@p�@p�@p�@`B@�@��@@@�T@�T@��@�@V@�@(�@��@ƨ@t�@�
@9X@(�@I�@I�@��@ƨ@1@Z@1@1@��@�F@S�@"�@
�\@
�@	�#@	��@	&�@��@1'@�;@�y@V@V@p�@�D@�
@o@M�@=q@�\@&�@ bN?��R?���?�1'?���?��T?���?��`?�"�?���?�I�?Ձ?� �?��m?��
?�p�?���?��/?���?��^?�+?��\?�/?��H?�Q�?�ff?�`B?�Ĝ?�~�?�`B?���?}�-?|(�?u?}?l1?f��?co?^v�?Y�#?P �?Fff??�w?;"�?7�P?1��?,�D?'�?"M�?"�?��?33?��?	7L?�7>�^5>�!>��>�+>�E�>�l�>��T>�;d>�n�>��7>hr�>]/>G�>333>�w>+=�h=���=ě�=��
=�\)=�+=aG�=+<���<o��o�8Q�y�#��%���
��E����������`B���+�	7L�hs��+��R�&�y�0 ž6E��:^5�@��P�`�["Ѿ\(��bMӾo���y�#�|푾�J������9��\)��t�����
=��"Ѿ��R��MӾ�xվ�1������-����KǾ��#��j��푾�|��o�š˾Ƨ��7L��I������bN��n���t��Ձ�ؓu�ڟ���5?��G����
���/���T��l���l�������1��� ž�&��-��33��F��?}��E����پ��#��^5��j���ۿ%�J�Mӿ��S��Z��/����9�	�^�
���ƨ�V��h������\)��׿녿�!��Ͽ9X������E��ȴ��ٿQ�������"ѿ����/���|�   � A�� Ĝ�!���"J�"�\�"��#���#�
�#�
�$��$Z�%��%�˿%�T�&ff�&��&�y�'l��'(r��(�ÿ)7L�)xտ)�^�*=q�+�+��+ƨ�,1�,I��,�Ϳ-O߿-��.V�.���/\)�/���0 ſ0bN�0�`�1&�1녿2�!�3t��3�Ͽ49X�49X�4z�5?}�5?}�5��6�6�6E��6ȴ�7
=�7Kǿ7Kǿ7�ٿ8b�8Q�8Q�8�u�8���8���9��9X�9���9���:^5�:^5�:^5�:�H�;"ѿ;"ѿ;"ѿ;��<(��;�m�<(��<(��<j�<��<��<��<��=/�=p��=p��=p��=�=�>vɿ>�R�>�R�>�ۿ>�ۿ?|�?�w�@A��@��@��@�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ������˿��/��\)���������P��`B���Ͽ~�ۿq녿bMӿQ���D���1hs�z�'�>Ǯ?�Ĝ?�=q?�$�?�bN?���?m��?<�?,�D?(r�?^5?��?J>�ff>��>�&�>�V>�|�>���>���>�X>�  >և+>ݲ->�t�>�33>�p�>�E�?o?�?�?$�?#S�?"J?�?%�?&��?-V?6�+?T9X?]p�?ix�?mV?mV?o��?st�?{��?�  ?�+?��R?�-?�?�v�?��?���?�o?��
?��/?�$�?���?�dZ?�bN?ӕ�?��/?�?}?��/?��?��?���?�J?�t�?�1'?��#@��@Ĝ@�j@V@/@&v�@-?}@,�/@0�u@:~�@<(�@<�j@@�u@AG�@@��@?�P@>E�@;C�@:�@9�@7|�@4I�@3"�@1hs@/��@/+@/��@/;d@-�@,�j@+t�@'�@$I�@$1@$j@%�-@&�R@'�P@(r�@*M�@+o@-`B@,�j@*-@)%@&�y@%�@&$�@%�-@%�@$��@#��@#C�@#33@!��@!X@ ��@ Q�@   @K�@��@ff@�@@@@ȴ@�@�y@ff@�y@v�@$�@��@V@?}@?}@�@�@/@O�@O�@?}@�/@z�@1@��@��@&�@�`@Q�@  @1'@�@��@��@�w@��@5?@�@"�@�\@��@��@  @  @  @A�@Q�@Q�@Q�@Q�@Q�@�@Ĝ@�`@&�@7L@7L@�@��@�9@��@  @b@A�@bN@Q�@A�@r�@r�@r�@ �@�;@�P@;d@�y@��@�+@$�@�T@��@p�@p�@p�@`B@�@��@@@�T@�T@��@�@V@�@(�@��@ƨ@t�@�
@9X@(�@I�@I�@��@ƨ@1@Z@1@1@��@�F@S�@"�@
�\@
�@	�#@	��@	&�@��@1'@�;@�y@V@V@p�@�D@�
@o@M�@=q@�\@&�@ bN?��R?���?�1'?���?��T?���?��`?�"�?���?�I�?Ձ?� �?��m?��
?�p�?���?��/?���?��^?�+?��\?�/?��H?�Q�?�ff?�`B?�Ĝ?�~�?�`B?���?}�-?|(�?u?}?l1?f��?co?^v�?Y�#?P �?Fff??�w?;"�?7�P?1��?,�D?'�?"M�?"�?��?33?��?	7L?�7>�^5>�!>��>�+>�E�>�l�>��T>�;d>�n�>��7>hr�>]/>G�>333>�w>+=�h=���=ě�=��
=�\)=�+=aG�=+<���<o��o�8Q�y�#��%���
��E����������`B���+�	7L�hs��+��R�&�y�0 ž6E��:^5�@��P�`�["Ѿ\(��bMӾo���y�#�|푾�J������9��\)��t�����
=��"Ѿ��R��MӾ�xվ�1������-����KǾ��#��j��푾�|��o�š˾Ƨ��7L��I������bN��n���t��Ձ�ؓu�ڟ���5?��G����
���/���T��l���l�������1��� ž�&��-��33��F��?}��E����پ��#��^5��j���ۿ%�J�Mӿ��S��Z��/����9�	�^�
���ƨ�V��h������\)��׿녿�!��Ͽ9X������E��ȴ��ٿQ�������"ѿ����/���|�   � A�� Ĝ�!���"J�"�\�"��#���#�
�#�
�$��$Z�%��%�˿%�T�&ff�&��&�y�'l��'(r��(�ÿ)7L�)xտ)�^�*=q�+�+��+ƨ�,1�,I��,�Ϳ-O߿-��.V�.���/\)�/���0 ſ0bN�0�`�1&�1녿2�!�3t��3�Ͽ49X�49X�4z�5?}�5?}�5��6�6�6E��6ȴ�7
=�7Kǿ7Kǿ7�ٿ8b�8Q�8Q�8�u�8���8���9��9X�9���9���:^5�:^5�:^5�:�H�;"ѿ;"ѿ;"ѿ;��<(��;�m�<(��<(��<j�<��<��<��<��=/�=p��=p��=p��=�=�>vɿ>�R�>�R�>�ۿ>�ۿ?|�?�w�@A��@��@��@�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBXB]/Bv�B�'B�BBJB�B �B6FBVBq�B�VB��B�LB�BS�B��B��B��B�B	�B	=qB	[#B	�B	�JB	�JB	��B	�LB	�^B	ɺB	�B	�B
  B	��B
B
bB
�B
�B
�B
!�B
$�B
:^B
A�B
?}B
D�B
K�B
P�B
R�B
\)B
bNB
dZB
e`B
hsB
n�B
t�B
�B
�1B
�uB
��B
��B
��B
��B
��B
�!B
�FB
ÖB
�jB
�B
�BB
�B
�B
��B
��B
��B
��BB+BJBbBVBoB{BhB"�B%�B(�B&�B6FB�BZBu�B�+B�JB�bB��B�B�3B�LB��B�B�B�TB�mB�fB�`B�ZB�TB�BB�BB�/B�B�B�B�B�B�B�B�
B��B��B��B��B��B��B��B�
B�B�)B�;B�TB�fB�`B�ZB�ZB�BB�5B�NB�BB�BB�;B�5B�/B�/B�5B�#B�#B�#B�B�B�B�B�B�B�B�B�)B�#B�/B�/B�/B�/B�5B�5B�)B�)B�)B�/B�)B�/B�/B�/B�/B�)B�/B�)B�#B�B�B�B�
B�
B�
B�B�
B�
B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBƨBƨBɺBȴBǮBŢB��B�}B�qB�dB�XB�RB�RB�FB�?B�3B�3B�9B�-B�3B�3B�?B�?B�9B�-B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B[�B`�BzUB��B�B�B�BB$QB9�BY�Bu6B��B�DB��B�BW�B�MB��B�xB�B	B	@�B	^�B	��B	��B	��B	�>B	��B	��B	�FB	�B	�*B
�B
 zB
�B
�B
B
!>B
!>B
%WB
(iB
=�B
EB
C	B
H(B
OSB
TqB
V~B
_�B
e�B
g�B
h�B
k�B
r$B
xHB
��B
��B
�B
�DB
�8B
�>B
�cB
�|B
��B
��B
�"B
��B
ۜB
��B
�*B
�0B
�gB
�sB
�gB zB�B
�B�B�B�B�BB�B&]B)oB,�B*uB9�B"DB]�ByOB��B��B��B�WB��B��B��B�xBِBۜB��B��B��B��B��B��B��B��B�BݩBݩBݩBۜBِBܣBۜBږB؊B؊B�kB�_B�YB�eB�xBږBܣBߵB��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��BޯBޯBޯBݩBݩBۜBܣBۜBܣBܣBܣBߵBޯB�B�B�B�B��B��BߵBߵBߵB�BߵB�B�B�B�BߵB�BߵBޯBܣBۜBۜBږBږBږBۜBږBږBږB؊BׄBׄB�qB�qB�eB�eB�YB�YB�_B�SB�eB�eB�_B�eB�eB�_B�kB�kB�kB�kB�kB�kB�kB�kB�eB�eB�kB�kB�kB�kB�kB�qB�kB�qB�qB�kB�kB�eB�_B�_B�_B�_B�_B�_B�_B�_B�YB�YB�_B�_B�_B�_B�_B�eB�_B�_B�_B�YB�SB�SB�SB�SB�YB�_B�_B�_B�_B�YB�_B�eB�kB�kB�kB�qB�kB�kB�kB�kB�eB�kB�kB�kB�eB�eB�eB�eB�_B�YB�YB�YB�SB�MB�MB�SB�_B�SB�SB�FB�@B�:B�4B�4B�FB�@B�:B�.B�B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B��B��B��B��B��B�|B�|B�|B�|B�uB�oB�iB�iB�cB�WB�QB�JB�JB�JB�QB�JB�DB�>B�>B�>B�8B�8B�>B�>B�>B�>B�>B�DB�DB�DB�>B�>B�>B�8B�8B�8B�2B�2B�2B�8B�8B�8B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�DB�DB�DB�DB�QB�JB�JB�JB�DB�JB�JB�QB�JB�JB�JB�JB�JB�JB�JB�JB�QB�QB�QB�QB�QB�QB�QB�QB�QB�WB�WB�QB�WB�WB�WB�WB�WB�]B�WB�]B�WB�]B�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�iB�iB�iB�oB�oB�iB�iB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542392021060715423920210607154239  IF  ARFMCODA029d                                                                20190620115650                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115720  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.2                                                                 20190620115720  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154239  IP  PSAL            A0  D���G�O�                