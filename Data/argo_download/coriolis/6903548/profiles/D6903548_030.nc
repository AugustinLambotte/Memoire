CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  R   	N_HISTORY          N_CALIB          
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
resolution        =���   axis      Z        	H  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  D    PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	H  Ft   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	H  R   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  [X   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  f�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  p<   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  {�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  �    PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  �t   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	H  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �(   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �,   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �0   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �4   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �X   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20190620115650  20210607154239  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @؇7����1   @؇7����@T\)t�V{@*�1=g�8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      @���@�33A��A��AffA,��A<��AQ��Aa��Aq��A���A�  A�33A�33A�33A���A�  A�ffA���A�  A�ffA�  A�  A���A�  A���B ffB  B  BffB33B  BffBffB   B$  B(ffB,  B/��B3��B8ffB;��B@��BC33BG��BL  BO��BT  BXffB\  B_��BdffBh  Bk��Bq33Bt  Bx��B|  B�  B�  B�  B�  B�33B�  B���B���B�33B�  B���B�33B�33B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B���B�  B�  B�33B���BǙ�B�33B�ffB�ffB���B�ffB���B�  B���B�ffB�33B�  B�  B�  B�  C� CffCffCL�C	��C��C��C� C� C� CffCffCffCffCffC� C!� C#� C%��C'��C)��C+�3C-ffC/ffC1��C3L�C5� C7��C9L�C;� C=�3C?ffCA33CC� CE�3CG� CIL�CK� CM��CO� CQffCS33CU� CW�3CY��C[ffC]L�C_33Ca� Cc��Ce�3Cg� CiffCkL�Cm33Co�Cq  CsffCu�3Cw��Cy� C{� C}ffCL�C���C���C���C��3C��3C��fC��fC�ٚC�ٚC���C���C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C�� C��3C��3C��3C�� C�� C���C���C���C�ٚC�ٚC��fC��fC��fC��3C�� C�� C���C���C���C���C��fC��3C��3C�� C�� C���C���C�ٚC��fC��fC��3C�� C�� C���C���C���C��fC��3C�� C���C�ٚC��fC��fC�� C CÙ�CĦfCų3C���C�ٚCȳ3Cɀ Cʌ�C�� C���C���C���Cϙ�CЦfCѦfCҳ3Cӳ3C�� C���C֦fCצfCس3C�� C�ٚC۳3C܌�CݦfC�� CߦfC�� C�ٚC�� C�fC�� C�ٚC�� C�fC��C�3C�ٚC�� C�3C홚C��C�3C��fC���C�� C�3C��fC���C���C��3C�ٚC���C�� C��3D 33Ds3D��D�fD&fDl�D�3D�3D
33Dl�D��D�fD&fDffD� D��DL�D� D�3D  DS3D�fD��D�3D33Dy�D �fD"�D#33D$` D%�3D'  D(L�D)y�D*�3D+��D-&fD.y�D/�3D0��D2&fD3s3D4� D5�3D7,�D8s3D9��D;  D<@ D=�fD>�fD@�DA9�DBffDC�3DE  DF33DG� DH�3DJ  DK9�DLl�DM� DN��DPY�DQ��DR��DS� DU  DVffDW��DX�3DZ@ D[� D\�fD^�D_L�D`� Da� DcfDd9�Del�Df� Dh�DiL�Dj��Dk��Dl�fDn33Do� Dp�3Dq��Ds@ Dt� Du� DwfDx,�Dyl�Dz�3D|  D}FfD~s3D��D��3D��D��fD�` D�fD��3D�@ D���D�|�D��D���D�\�D���D���D�<�D��3D�|�D��D���D�Y�D���D���D�9�D���D�� D��D��fD�c3D�3D��3D�C3D�ٚD�|�D�&fD�� D�Y�D��fD��3D�33D��3D�s3D�fD���D�\�D�  D���D�@ D��fD�� D��D���D�c3D�  D���D�<�D���D�� D�  D��3D�i�D�  D��fD�9�D�� D�y�D�fD�� D�c3D�  D��3D�@ D�� D�� D�#3D���D�` D�fD�� D�9�D�ٚD�y�D��D���D�Y�D���D���D�@ D��3D��3D�#3D���D�VfD���D��3D�C3D�� D��3D�&fD���D�S3D���D�� D�9�D��fD�� D�  D�� D�` D��fDĜ�D�@ D�ٚD�s3D��D��3D�` D���DɖfD�33D�� D�|�D�)�D�ɚD�ffD�fDΩ�D�I�D���D�|�D� DѶfD�VfD���DӠ D�9�D��fD�p D��D�ɚD�` D��3DؖfD�<�D��fDڀ D��D۹�D�VfD�3Dݣ3D�@ D�� D߀ D�  D��3D�ffD���D�3D�9�D��3D� D��D�3D�` D�	�D�3D�@ D�ٚD�vfD�3D�3D�S3D��3D�3D�33D��fD�vfD��D﹚D�\�D�3D��D�9�D��3D�|�D�&fD��3D�ffD���D��fD�@ D��fD�� D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @���@�33A��A��AffA,��A<��AQ��Aa��Aq��A���A�  A�33A�33A�33A���A�  A�ffA���A�  A�ffA�  A�  A���A�  A���B ffB  B  BffB33B  BffBffB   B$  B(ffB,  B/��B3��B8ffB;��B@��BC33BG��BL  BO��BT  BXffB\  B_��BdffBh  Bk��Bq33Bt  Bx��B|  B�  B�  B�  B�  B�33B�  B���B���B�33B�  B���B�33B�33B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B���B�  B�  B�33B���BǙ�B�33B�ffB�ffB���B�ffB���B�  B���B�ffB�33B�  B�  B�  B�  C� CffCffCL�C	��C��C��C� C� C� CffCffCffCffCffC� C!� C#� C%��C'��C)��C+�3C-ffC/ffC1��C3L�C5� C7��C9L�C;� C=�3C?ffCA33CC� CE�3CG� CIL�CK� CM��CO� CQffCS33CU� CW�3CY��C[ffC]L�C_33Ca� Cc��Ce�3Cg� CiffCkL�Cm33Co�Cq  CsffCu�3Cw��Cy� C{� C}ffCL�C���C���C���C��3C��3C��fC��fC�ٚC�ٚC���C���C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C�� C��3C��3C��3C�� C�� C���C���C���C�ٚC�ٚC��fC��fC��fC��3C�� C�� C���C���C���C���C��fC��3C��3C�� C�� C���C���C�ٚC��fC��fC��3C�� C�� C���C���C���C��fC��3C�� C���C�ٚC��fC��fC�� C CÙ�CĦfCų3C���C�ٚCȳ3Cɀ Cʌ�C�� C���C���C���Cϙ�CЦfCѦfCҳ3Cӳ3C�� C���C֦fCצfCس3C�� C�ٚC۳3C܌�CݦfC�� CߦfC�� C�ٚC�� C�fC�� C�ٚC�� C�fC��C�3C�ٚC�� C�3C홚C��C�3C��fC���C�� C�3C��fC���C���C��3C�ٚC���C�� C��3D 33Ds3D��D�fD&fDl�D�3D�3D
33Dl�D��D�fD&fDffD� D��DL�D� D�3D  DS3D�fD��D�3D33Dy�D �fD"�D#33D$` D%�3D'  D(L�D)y�D*�3D+��D-&fD.y�D/�3D0��D2&fD3s3D4� D5�3D7,�D8s3D9��D;  D<@ D=�fD>�fD@�DA9�DBffDC�3DE  DF33DG� DH�3DJ  DK9�DLl�DM� DN��DPY�DQ��DR��DS� DU  DVffDW��DX�3DZ@ D[� D\�fD^�D_L�D`� Da� DcfDd9�Del�Df� Dh�DiL�Dj��Dk��Dl�fDn33Do� Dp�3Dq��Ds@ Dt� Du� DwfDx,�Dyl�Dz�3D|  D}FfD~s3D��D��3D��D��fD�` D�fD��3D�@ D���D�|�D��D���D�\�D���D���D�<�D��3D�|�D��D���D�Y�D���D���D�9�D���D�� D��D��fD�c3D�3D��3D�C3D�ٚD�|�D�&fD�� D�Y�D��fD��3D�33D��3D�s3D�fD���D�\�D�  D���D�@ D��fD�� D��D���D�c3D�  D���D�<�D���D�� D�  D��3D�i�D�  D��fD�9�D�� D�y�D�fD�� D�c3D�  D��3D�@ D�� D�� D�#3D���D�` D�fD�� D�9�D�ٚD�y�D��D���D�Y�D���D���D�@ D��3D��3D�#3D���D�VfD���D��3D�C3D�� D��3D�&fD���D�S3D���D�� D�9�D��fD�� D�  D�� D�` D��fDĜ�D�@ D�ٚD�s3D��D��3D�` D���DɖfD�33D�� D�|�D�)�D�ɚD�ffD�fDΩ�D�I�D���D�|�D� DѶfD�VfD���DӠ D�9�D��fD�p D��D�ɚD�` D��3DؖfD�<�D��fDڀ D��D۹�D�VfD�3Dݣ3D�@ D�� D߀ D�  D��3D�ffD���D�3D�9�D��3D� D��D�3D�` D�	�D�3D�@ D�ٚD�vfD�3D�3D�S3D��3D�3D�33D��fD�vfD��D﹚D�\�D�3D��D�9�D��3D�|�D�&fD��3D�ffD���D��fD�@ D��fD�� D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����Ĝ�Ͼw���H��Z��33��녿��;���;��󶿌(��$����1��w>vȴ?%?t�?"�\?.�?e`B?z�?v�+?��?�I�?ə�?�
=?���?�A�?�hs?Ͳ-?�j?��#?Ł?��w?�/?�p�?o�;?n�?r-?|�?���?�X?���?��?z��?s��?}�?�n�?�+?��?�^5?���?��?��?��7?�`B?���?�ff?���?��?��?ļj?�
=?��y?�\?�r�?�M�?���?�v�@X@ƨ@v�@+@�@�@r�@$�@�7@ 1'?���?���@��@n�@&�@^5@�@=q@hs@X@�@9X@�j@ff@ff@{@{@��@@ff@
^5@Z@��@��@��@ƨ@`B@��@"��@$1@(��@(A�@'\)@)��@+S�@*�H@*�\@.v�@.�+@-�T@-/@-p�@.�R@,��@*�@)X@)�^@+�
@+��@*�!@'\)@'��@(A�@*n�@)7L@)�@*M�@+C�@,z�@,9X@,Z@-/@.@-�@-��@.$�@-�@+ƨ@*^5@)�7@)7L@(��@'�@&��@&��@&��@&v�@&E�@&E�@%�T@&��@&v�@&�+@&ff@&5?@%��@%�h@%`B@$�@$��@#"�@"�\@!��@!��@!��@!��@!�#@"J@"�@"-@"-@"M�@"=q@"�@!�^@ Q�@p�@�@�
@9X@9X@j@z�@z�@Z@dZ@J@r�@��@Z@�F@��@�
@�m@�@�@9X@Z@z�@�@�@�@�@�@�@1@�!@�!@^5@�^@%@�y@{@�T@��@�T@�T@�T@�T@�@E�@v�@;d@�P@�P@�@ȴ@��@E�@�@ȴ@�R@��@�y@+@+@;d@
=@��@��@�R@ff@�@�-@p�@��@��@�@
�@	�^@�@�@�@��@	&�@��@�@b@�@+@E�@�T@`B@��@�
@^5@x�@ Ĝ?�\)?�O�?���?�7L?�Q�?�K�?��/?���?�hs?�;d?�v�?��?�V?��H?�=q?陚?�?}?�|�?ߝ�?��?���?ա�?��?��;?͑h?ɺ^?�
=?�t�?�bN?��m?�=q?�$�?��?�5??��#?��?���?��;?���?��?���?�I�?��y?�%?w
=?o\)?j��?fff?_|�?WK�?N�?K�?H��?E�?;dZ?9��?9�#?/��?(��?"M�?��?{?��?�\?�7>�p�>�33>�~�>�"�>ȴ9>��F>�&�>�>�Z>��R>���>��>�ƨ>o��>I�^>1&�>"��>t�>   =��=\=��T=�7L=q��=49X=C�<ě�<t�    �49X�ě��t��D���y�#��\)��1��^5�Ƨ��
=����C��\)����%�T�-V�49X�F��N��["Ѿgl��q���|푾��7������+��I���V��񪾖��
=��"Ѿ�;d���
���
��`B��l���xվ�������� ž�&龵?}���پ��#������\�Õ���$ݾǮ�ȴ9��=q��ƨ��O߾�hs��t������և+��b�����"Ѿ޸R������Z���y��xվ���V�����&��-��?}��KǾ�Q��X��^5���m����|� Ĝ��7��\��
��
��/��T��y�+�	7L�	xտ���1�1��Ϳ�����\)����bN�&����녿-��!�33��Ͽz�����E��ȴ�
=�b����X�^5��H���(����/�5?�vɿ�R��R�;d�|��w� A�� Ĝ�!�7�"J�"�\�#o�#S��#�
�$��$���$�/�%��%`B�%�˿%�T�&$ݿ&��&��&�y�'l��'��(1'�(r��(�ÿ(�ÿ)xտ)�^�)��*=q�+�+C��+ƨ�,1�,1�,I��,�Ϳ-V�-O߿-�h�-��.{�.V�.���.��/\)�/\)�/���0 ſ0bN�0�׿0�`�1&�1hs�1녿2n��2�!�2�3�Ͽ3�F�3�Ͽ4���5?}�6E��6�+�6�+�7
=�7
=�7�P�8b�8Q�8Q�8Q�8�u�9��9X�9���9�#�:��:���:���:�H�;dZ�;�m�;�m�<(��<j�<��<푿=/�=p��=�-�=�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ��Ĝ�Ͼw���H��Z��33��녿��;���;��󶿌(��$����1��w>vȴ?%?t�?"�\?.�?e`B?z�?v�+?��?�I�?ə�?�
=?���?�A�?�hs?Ͳ-?�j?��#?Ł?��w?�/?�p�?o�;?n�?r-?|�?���?�X?���?��?z��?s��?}�?�n�?�+?��?�^5?���?��?��?��7?�`B?���?�ff?���?��?��?ļj?�
=?��y?�\?�r�?�M�?���?�v�@X@ƨ@v�@+@�@�@r�@$�@�7@ 1'?���?���@��@n�@&�@^5@�@=q@hs@X@�@9X@�j@ff@ff@{@{@��@@ff@
^5@Z@��@��@��@ƨ@`B@��@"��@$1@(��@(A�@'\)@)��@+S�@*�H@*�\@.v�@.�+@-�T@-/@-p�@.�R@,��@*�@)X@)�^@+�
@+��@*�!@'\)@'��@(A�@*n�@)7L@)�@*M�@+C�@,z�@,9X@,Z@-/@.@-�@-��@.$�@-�@+ƨ@*^5@)�7@)7L@(��@'�@&��@&��@&��@&v�@&E�@&E�@%�T@&��@&v�@&�+@&ff@&5?@%��@%�h@%`B@$�@$��@#"�@"�\@!��@!��@!��@!��@!�#@"J@"�@"-@"-@"M�@"=q@"�@!�^@ Q�@p�@�@�
@9X@9X@j@z�@z�@Z@dZ@J@r�@��@Z@�F@��@�
@�m@�@�@9X@Z@z�@�@�@�@�@�@�@1@�!@�!@^5@�^@%@�y@{@�T@��@�T@�T@�T@�T@�@E�@v�@;d@�P@�P@�@ȴ@��@E�@�@ȴ@�R@��@�y@+@+@;d@
=@��@��@�R@ff@�@�-@p�@��@��@�@
�@	�^@�@�@�@��@	&�@��@�@b@�@+@E�@�T@`B@��@�
@^5@x�@ Ĝ?�\)?�O�?���?�7L?�Q�?�K�?��/?���?�hs?�;d?�v�?��?�V?��H?�=q?陚?�?}?�|�?ߝ�?��?���?ա�?��?��;?͑h?ɺ^?�
=?�t�?�bN?��m?�=q?�$�?��?�5??��#?��?���?��;?���?��?���?�I�?��y?�%?w
=?o\)?j��?fff?_|�?WK�?N�?K�?H��?E�?;dZ?9��?9�#?/��?(��?"M�?��?{?��?�\?�7>�p�>�33>�~�>�"�>ȴ9>��F>�&�>�>�Z>��R>���>��>�ƨ>o��>I�^>1&�>"��>t�>   =��=\=��T=�7L=q��=49X=C�<ě�<t�    �49X�ě��t��D���y�#��\)��1��^5�Ƨ��
=����C��\)����%�T�-V�49X�F��N��["Ѿgl��q���|푾��7������+��I���V��񪾖��
=��"Ѿ�;d���
���
��`B��l���xվ�������� ž�&龵?}���پ��#������\�Õ���$ݾǮ�ȴ9��=q��ƨ��O߾�hs��t������և+��b�����"Ѿ޸R������Z���y��xվ���V�����&��-��?}��KǾ�Q��X��^5���m����|� Ĝ��7��\��
��
��/��T��y�+�	7L�	xտ���1�1��Ϳ�����\)����bN�&����녿-��!�33��Ͽz�����E��ȴ�
=�b����X�^5��H���(����/�5?�vɿ�R��R�;d�|��w� A�� Ĝ�!�7�"J�"�\�#o�#S��#�
�$��$���$�/�%��%`B�%�˿%�T�&$ݿ&��&��&�y�'l��'��(1'�(r��(�ÿ(�ÿ)xտ)�^�)��*=q�+�+C��+ƨ�,1�,1�,I��,�Ϳ-V�-O߿-�h�-��.{�.V�.���.��/\)�/\)�/���0 ſ0bN�0�׿0�`�1&�1hs�1녿2n��2�!�2�3�Ͽ3�F�3�Ͽ4���5?}�6E��6�+�6�+�7
=�7
=�7�P�8b�8Q�8Q�8Q�8�u�9��9X�9���9�#�:��:���:���:�H�;dZ�;�m�;�m�<(��<j�<��<푿=/�=p��=�-�=�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB{B�B)�B>wBA�BG�B\)B}�B�1B�/Bu�B��B�BK�B��B�NB��B+B7LBC�Bm�B�7B�dB�B�sB		7B	2-B	R�B	R�B	��B	��B	�jB	��B	�
B
1B
(�B
9XB
H�B
YB
p�B
� B
�B
�B
�B
�7B
�hB
��B
��B
��B
�?B
B
ǮB
��B
��B
��B
�/B
�TB
�sB
�B
��B
��B
��BVB�B%�B2-B:^B?}BB�BI�BM�BQ�BR�BQ�BS�BR�BW
BN�BO�BN�B33B\)B\)BW
BT�BXBZB[#B`BBdZBgmBhsBn�Bn�Bm�Bo�Bp�Bv�B}�B�%B�hB��B��B��B�B�-B�dB�qBȴBɺBǮB��B��B��B��B�B�5B�#B�)B�)B�BB�5B�
B�B�)B�;B�BB�5B�#B�)B�#B�HB�HB�NB�ZB�fB�sB�yB�B�B�B�B�B�B�B�B�B�B�yB�yB�sB�sB�mB�mB�sB�sB�yB�sB�B�B�B�B�B�B�B�B�yB�sB�sB�fB�`B�`B�fB�`B�fB�fB�fB�mB�mB�mB�mB�fB�`B�NB�HB�;B�;B�BB�BB�BB�BB�BB�;B�BB�5B�B�B�B��B��B��B�B�B�
B�B�B�
B�
B�
B�
B�
B�
B�
B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBȴBɺBȴBŢBŢBŢBĜBŢBĜBÖBBB�wB�}B�}B�qB�dB�jB�dB�dB�^B�jB�dB�^B�^B�XB�^B�XB�RB�FB�FB�?B�9B�-B�-B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BBB-�BBBEBK:B_�B��B��B�ByOB�eB"DBOSB�8B��B�gB
�B:�BG"BqB��B��BۜB��B	�B	5�B	V~B	V~B	�B	�8B	��B	�xB	ږB
�B
,�B
<�B
L@B
\�B
t0B
��B
��B
��B
��B
��B
��B
�&B
�WB
�WB
��B
�B
�:B
�YB
�qB
؊B
�B
��B
��B
�*B
�UB
�mB zB�B"DB)oB5�B=�BC	BFBMFBQ_BUxBV~BUxBW�BV~BZ�BReBSkBReB6�B_�B_�BZ�BX�B[�B]�B^�Bc�Bg�Bj�Bk�Br$Br$BqBs*Bt0BzUB��B��B��B�>B�QB��B��B��B��B��B�@B�FB�:B�_B�xB�qBׄBِB��BޯBߵBߵB��B��BږBݩBߵB��B��B��BޯBߵBޯB��B��B��B��B��B��B�B�B�B�$B�$B�*B�$B�$B�B�B�B�B�B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BܣBۜBِB؊B؊B؊BِBِBږBِBِBږBږBږBږBږBږBږBږBׄBׄB؊B�~B؊B�eB�eB�eB�eB�eB�eB�eB�eB�eB�eB�kB�qB�xB�qB�xB�qB�xB�eB�xB�xB�~B�xB�~BׄBׄBׄBِB؊BׄBׄB؊BׄBׄBׄBׄB�~B�qB�qB�qB�eB�_B�eB�qB�qB�qB�qB�qB�qB�qB�kB�kB�kB�eB�_B�_B�_B�SB�MB�MB�FB�FB�@B�@B�FB�@B�.B�.B�.B�(B�.B�(B�"B�B�B�B�	B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�oB�oB�|B�|B�cB�oB�|B�uB�|B�oB�cB�WB�WB�WB�]B�cB�cB�]B�cB�QB�JB�QB�QB�JB�QB�JB�JB�DB�DB�8B�8B�8B�8B�2B�2B�8B�8B�8B�8B�8B�8B�2B�2B�8B�2B�2B�2B�8B�8B�8B�>B�>B�>B�>B�>B�>B�DB�DB�DB�DB�DB�JB�DB�JB�DB�DB�>B�>B�>B�>B�DB�>B�>B�>B�>B�>B�>B�DB�DB�DB�DB�DB�DB�JB�DB�JB�JB�QB�QB�QB�QB�QB�QB�QB�QB�WB�WB�QB�WB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�]B�cB�iB�cB�cB�iB�cB�iB�iB�cB�cB�iB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B�|B��B��B��B��B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542392021060715423920210607154239  IF  ARFMCODA029d                                                                20190620115650                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115719  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.2                                                                 20190620115719  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154239  IP  PSAL            @���D��3G�O�                