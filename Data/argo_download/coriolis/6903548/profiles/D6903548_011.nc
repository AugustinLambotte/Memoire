CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  S   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:49Z creation; 2021-06-07T15:42:38Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        	L  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  D$   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	L  Fx   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	L  R   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  [d   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  g   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  pP   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  {�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  �<   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  �0   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �(   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �8   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �<   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �L   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �P   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �T   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �X   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �|   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20190620115649  20210607154238  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @�x�l�l1   @�x�$�8�@TjQ�u�@,8V�%1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 5.2 dbar]                                                      @���@�33@�33A��A��AffA.ffAC33AP  Ac33Aq��A�  A�  A�33A���A���A���A�33A�ffA�ffA���A�ffA�  A�33A�  A���A�  B   B  B  B  B  B  B��BffB   B$��B'33B*��B.��B333B7��B<  B@ffBD  BG��BLffBO��BS33BW��B[��B_��Bc��Bh  BlffBo��BtffBx  B{��B��B�  B�ffB���B�ffB�ffB�ffB���B�  B���B���B�33B�33B�33B�33B���B���B���B���B���B���B�  B�  B�33B���B���B���B�  B�33B���B���B�  B�ffB���B���B�  B�33B�33B���B�ffB���B�ffB�  B���B晚B�ffB���B�B�ffB���B���C33C� CffCL�C	� C� CffCL�CL�C��C��C��C��C��C��C33C!L�C#ffC%ffC'� C)��C+��C-L�C/ffC1� C333C5L�C7� C933C;L�C=� C?33CAL�CCffCE� CG�3CIffCK� CM��COffCQ� CS�3CUffCW33CYffC[� C]L�C_� Ca��CcffCe� Cg�3Ci� Ck33Cm� Co�3Cq� CsL�Cu�CwffCy��C{ffC}33C� C�ٚC���C��3C���C���C��3C��fC���C�� C��3C��fC��fC���C���C��3C��fC��fC���C���C�� C��3C��fC��fC���C���C���C�� C��3C��fC�ٚC�ٚC�ٚC�ٚC���C���C���C���C���C���C���C�ٚC��fC��fC��3C��3C�� C�� C���C���C���C��fC�� C���C�ٚC��fC��3C�� C���C��fC��3C�� C�ٚC��fC�� C�� C���C³3C�� C�� C�� C���C���CȦfCɦfCʳ3C�� C���C͙�CΦfCϳ3C�� C���CҦfCӦfCԳ3C�� C���Cי�CئfCٳ3C�� C���Cܙ�CݦfC޳3C�� C���C�fC�3C���C�fC�3C�� C癚C�3C�� C�fC�� C�ٚC��3CC�3C���C�fC��C�fC�� C��fC���C��fC�� C��fC�s3C��3D 9�Dy�D��D  D@ D�fD�3D��D
L�D� D��D��D&fDffD�fD��D,�Dl�D�3D��D@ Ds3D�3D  D33Dl�D �fD!��D#33D$y�D%�fD&��D(,�D)� D*� D+��D-9�D.y�D/� D1  D2FfD3l�D4� D63D7L�D8�fD9��D:�3D<,�D=l�D>��D?�fDA&fDBffDC�fDD�3DF9�DG�fDH�3DI�fDK9�DL��DM� DN�3DPL�DQ�fDR�3DS��DU9�DV��DW��DY�DZFfD[�fD\�fD]�3D_&fD`s3Da� Dc�Dd@ Del�Df� Dg�3Di&fDjy�Dk�3Dl��DnFfDo�fDp�fDrfDsFfDt�fDu��Dv�3Dx9�DyffDz��D|fD}9�D~l�D��D��fD�  D���D�VfD�� D�� D�0 D�� D�s3D�fD���D�\�D�  D��3D�9�D�� D�� D�)�D�ɚD�\�D��3D���D�@ D�ٚD�vfD�#3D��3D�VfD���D�� D�6fD�� D���D�#3D�� D�\�D���D���D�9�D�ٚD�|�D�  D��3D�\�D�  D��fD�<�D��D�|�D� D��3D�VfD���D�� D�9�D��fD��3D�#3D�� D�` D�3D��fD�<�D�� D��fD��D��fD�\�D�	�D��fD�C3D��3D��3D�#3D��3D�c3D�	�D���D�33D�ٚD�y�D��D��3D�\�D���D��3D�<�D��fD���D��D��3D�Y�D�  D��fD�@ D�ٚD��fD�#3D�� D�\�D���D���D�@ D��3D�y�D�  D��fD�\�D��3DĜ�D�FfD�� D�y�D�fDǳ3D�P D�� Dɐ D�0 D�� D�p D�3D̳3D�S3D��fDΙ�D�@ D��fD�y�D� Dѹ�D�c3D���DӜ�D�@ D�� DՀ D�  D��3D�\�D��3D؜�D�<�D�ٚD�y�D�fD�� D�` D�  Dݠ D�<�D�� D߀ D�#3D�fD�\�D�3D��D�9�D��3D�|�D�,�D�ɚD�l�D��D��D�<�D�� D�s3D��D��D�ffD�  D왚D�6fD���DD�)�D��D�P D��3D�D�C3D�� D�y�D��D��fD�VfD��fD��3D�C3D�� D��3D��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @���@�33@�33A��A��AffA.ffAC33AP  Ac33Aq��A�  A�  A�33A���A���A���A�33A�ffA�ffA���A�ffA�  A�33A�  A���A�  B   B  B  B  B  B  B��BffB   B$��B'33B*��B.��B333B7��B<  B@ffBD  BG��BLffBO��BS33BW��B[��B_��Bc��Bh  BlffBo��BtffBx  B{��B��B�  B�ffB���B�ffB�ffB�ffB���B�  B���B���B�33B�33B�33B�33B���B���B���B���B���B���B�  B�  B�33B���B���B���B�  B�33B���B���B�  B�ffB���B���B�  B�33B�33B���B�ffB���B�ffB�  B���B晚B�ffB���B�B�ffB���B���C33C� CffCL�C	� C� CffCL�CL�C��C��C��C��C��C��C33C!L�C#ffC%ffC'� C)��C+��C-L�C/ffC1� C333C5L�C7� C933C;L�C=� C?33CAL�CCffCE� CG�3CIffCK� CM��COffCQ� CS�3CUffCW33CYffC[� C]L�C_� Ca��CcffCe� Cg�3Ci� Ck33Cm� Co�3Cq� CsL�Cu�CwffCy��C{ffC}33C� C�ٚC���C��3C���C���C��3C��fC���C�� C��3C��fC��fC���C���C��3C��fC��fC���C���C�� C��3C��fC��fC���C���C���C�� C��3C��fC�ٚC�ٚC�ٚC�ٚC���C���C���C���C���C���C���C�ٚC��fC��fC��3C��3C�� C�� C���C���C���C��fC�� C���C�ٚC��fC��3C�� C���C��fC��3C�� C�ٚC��fC�� C�� C���C³3C�� C�� C�� C���C���CȦfCɦfCʳ3C�� C���C͙�CΦfCϳ3C�� C���CҦfCӦfCԳ3C�� C���Cי�CئfCٳ3C�� C���Cܙ�CݦfC޳3C�� C���C�fC�3C���C�fC�3C�� C癚C�3C�� C�fC�� C�ٚC��3CC�3C���C�fC��C�fC�� C��fC���C��fC�� C��fC�s3C��3D 9�Dy�D��D  D@ D�fD�3D��D
L�D� D��D��D&fDffD�fD��D,�Dl�D�3D��D@ Ds3D�3D  D33Dl�D �fD!��D#33D$y�D%�fD&��D(,�D)� D*� D+��D-9�D.y�D/� D1  D2FfD3l�D4� D63D7L�D8�fD9��D:�3D<,�D=l�D>��D?�fDA&fDBffDC�fDD�3DF9�DG�fDH�3DI�fDK9�DL��DM� DN�3DPL�DQ�fDR�3DS��DU9�DV��DW��DY�DZFfD[�fD\�fD]�3D_&fD`s3Da� Dc�Dd@ Del�Df� Dg�3Di&fDjy�Dk�3Dl��DnFfDo�fDp�fDrfDsFfDt�fDu��Dv�3Dx9�DyffDz��D|fD}9�D~l�D��D��fD�  D���D�VfD�� D�� D�0 D�� D�s3D�fD���D�\�D�  D��3D�9�D�� D�� D�)�D�ɚD�\�D��3D���D�@ D�ٚD�vfD�#3D��3D�VfD���D�� D�6fD�� D���D�#3D�� D�\�D���D���D�9�D�ٚD�|�D�  D��3D�\�D�  D��fD�<�D��D�|�D� D��3D�VfD���D�� D�9�D��fD��3D�#3D�� D�` D�3D��fD�<�D�� D��fD��D��fD�\�D�	�D��fD�C3D��3D��3D�#3D��3D�c3D�	�D���D�33D�ٚD�y�D��D��3D�\�D���D��3D�<�D��fD���D��D��3D�Y�D�  D��fD�@ D�ٚD��fD�#3D�� D�\�D���D���D�@ D��3D�y�D�  D��fD�\�D��3DĜ�D�FfD�� D�y�D�fDǳ3D�P D�� Dɐ D�0 D�� D�p D�3D̳3D�S3D��fDΙ�D�@ D��fD�y�D� Dѹ�D�c3D���DӜ�D�@ D�� DՀ D�  D��3D�\�D��3D؜�D�<�D�ٚD�y�D�fD�� D�` D�  Dݠ D�<�D�� D߀ D�#3D�fD�\�D�3D��D�9�D��3D�|�D�,�D�ɚD�l�D��D��D�<�D�� D�s3D��D��D�ffD�  D왚D�6fD���DD�)�D��D�P D��3D�D�C3D�� D�y�D��D��fD�VfD��fD��3D�C3D�� D��3D��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�T@ ��@�/@Q�?�l�?�7L?���?�b?t�j?[dZ?,1���B�忄Z��S����m�X�u�X�u�W�ٿR�!�:��� ſ#S��+��;"ѿq녿����������푿�"ѿ�����P��E��y���>�ۿhs��?}� ��!�7�   �5?�#�
��ۿ녾�?}���~����������և+��bN���H��/�N��\)���#��`B=aG�>6E�>�t�>�;d>��>�\)>��m?	�^?��?��?MO�?q�?}�-?���?��?�"�?�j?���?��?�/?�5??�G�?���?�bN?�o?�E�?��?�^5?�l�?և+?���?��T?�?�?�dZ@(�@�@ƨ@p�@�@\)@S�@�-?��?���?��?��R@%@~�@��@�@��@v�@b@	�^@
�@S�@V@V@��@O�@�@�^@33@ r�@#�@#ƨ@%p�@'l�@'�P@#�
@"~�@"�!@"��@!�#@"=q@ b@|�@"-@#��@#��@$9X@$�j@%�h@&{@&$�@';d@&ff@&$�@'�@'�@'�w@'|�@'�;@(bN@(1'@'�w@%�@#ƨ@#�@$I�@$z�@$��@$��@&$�@&$�@&$�@&$�@&{@%�@!%@�+@V@E�@V@ȴ@�@;d@V@�m@dZ@o@�\@�^@�@�`@�@|�@�w@l�@
=@�y@��@;d@|�@��@  @�@E�@�@�/@�
@"�@��@n�@�@�^@��@�^@�#@9X@5?@v�@��@�y@K�@K�@�@  @\)@l�@\)@l�@l�@\)@\)@�P@�P@b@�@Q�@�@l�@��@v�@�y@��@\)@�@�+@5?@V@z�@�/@�D@1@S�@�F@	�#@	��@(�@��@��@1@t�@
�H@
�\@
M�@	X@�@
=@@I�@��?�"�?��u?��?�E�?���?���?�|�?���?��m?���?��?��/?��?޸R?�p�?���?���?�~�?ؓu?�ff?��?��?�?}?��?�7L?׍P?��?�z�?�bN?��?�?��?��?�j?��?�b?��?���?���?�33?�&�?��-?��/?�&�?�/?���?�~�?�7L?��9?���?���?�l�?}�-?{�m?z^5?g�?\�?T�j?P �?Rn�?M�h?E��?@  ?4z�?)�^?"�\?   ?X?z�?�?��?Z>�v�>���>ڟ�>��;>��7>�?}>�>�hs>�$�>��7>|�>u>j~�>T��>>v�>-V>#�
>�w>n�>�=�S�=�j=�t�=q��=8Q�=\)<�j<49X�ě���C��ě���h��w�y�#�����Ƨ�������m�C���+�&�y�,1�2-�5?}�:^5�D���F��E�˾Kƨ�Y��_;d�cS��w�پ�  ���7������$ݾ��^���;��n���������/�������
���y���羬1���D�������F����ȴ�������H���7��1'��O߾�녾���t���b�ڟ��ܬ�޸R�߾w��A����
��`B���������ȴ���H���m��p�����|� A��%����o�Z�����/��˿ff�l��r��	7L�	��
=q����V�{���\)�\)��;��`�hs�녿-���F��Ͽ9X�z����?}��+��P�b��u����������dZ�(����/�푿/�p��p���5?�vɿ�R�;d�|��w� �� Ĝ�!%�!�7�"J�"�\�"�\�"��#S��$Z�$�/�%`B�%`B�%�˿&$ݿ&��'l��'��(1'�(�9�)xտ)��*���+�+�+ƨ�,1�,�D�-V�-�h�.{�/\)�0bN�0bN�0bN�0�׿0�`�1hs�1녿2n��2�3t��3�F�3�Ͽ49X�49X�4z�4z�49X�49X�4z�4�j�4�j�4���5?}�5?}�5?}�5��5��6E��6�+�6ȴ�6ȴ�7
=�7Kǿ7�P�7�ٿ8b�8Q�8�u�9��9X�9���:��:^5�:���;"ѿ;dZ�;�m�<(��<j�<��<푿=/�=/�=/�=p��=�-�=�-�=�-�>5?�>5?�>vɿ>�R�>�R�>�ۿ?;d�?�w�?�w�@  �@  �@A��@A��@A��@A��@�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�T@ ��@�/@Q�?�l�?�7L?���?�b?t�j?[dZ?,1���B�忄Z��S����m�X�u�X�u�W�ٿR�!�:��� ſ#S��+��;"ѿq녿����������푿�"ѿ�����P��E��y���>�ۿhs��?}� ��!�7�   �5?�#�
��ۿ녾�?}���~����������և+��bN���H��/�N��\)���#��`B=aG�>6E�>�t�>�;d>��>�\)>��m?	�^?��?��?MO�?q�?}�-?���?��?�"�?�j?���?��?�/?�5??�G�?���?�bN?�o?�E�?��?�^5?�l�?և+?���?��T?�?�?�dZ@(�@�@ƨ@p�@�@\)@S�@�-?��?���?��?��R@%@~�@��@�@��@v�@b@	�^@
�@S�@V@V@��@O�@�@�^@33@ r�@#�@#ƨ@%p�@'l�@'�P@#�
@"~�@"�!@"��@!�#@"=q@ b@|�@"-@#��@#��@$9X@$�j@%�h@&{@&$�@';d@&ff@&$�@'�@'�@'�w@'|�@'�;@(bN@(1'@'�w@%�@#ƨ@#�@$I�@$z�@$��@$��@&$�@&$�@&$�@&$�@&{@%�@!%@�+@V@E�@V@ȴ@�@;d@V@�m@dZ@o@�\@�^@�@�`@�@|�@�w@l�@
=@�y@��@;d@|�@��@  @�@E�@�@�/@�
@"�@��@n�@�@�^@��@�^@�#@9X@5?@v�@��@�y@K�@K�@�@  @\)@l�@\)@l�@l�@\)@\)@�P@�P@b@�@Q�@�@l�@��@v�@�y@��@\)@�@�+@5?@V@z�@�/@�D@1@S�@�F@	�#@	��@(�@��@��@1@t�@
�H@
�\@
M�@	X@�@
=@@I�@��?�"�?��u?��?�E�?���?���?�|�?���?��m?���?��?��/?��?޸R?�p�?���?���?�~�?ؓu?�ff?��?��?�?}?��?�7L?׍P?��?�z�?�bN?��?�?��?��?�j?��?�b?��?���?���?�33?�&�?��-?��/?�&�?�/?���?�~�?�7L?��9?���?���?�l�?}�-?{�m?z^5?g�?\�?T�j?P �?Rn�?M�h?E��?@  ?4z�?)�^?"�\?   ?X?z�?�?��?Z>�v�>���>ڟ�>��;>��7>�?}>�>�hs>�$�>��7>|�>u>j~�>T��>>v�>-V>#�
>�w>n�>�=�S�=�j=�t�=q��=8Q�=\)<�j<49X�ě���C��ě���h��w�y�#�����Ƨ�������m�C���+�&�y�,1�2-�5?}�:^5�D���F��E�˾Kƨ�Y��_;d�cS��w�پ�  ���7������$ݾ��^���;��n���������/�������
���y���羬1���D�������F����ȴ�������H���7��1'��O߾�녾���t���b�ڟ��ܬ�޸R�߾w��A����
��`B���������ȴ���H���m��p�����|� A��%����o�Z�����/��˿ff�l��r��	7L�	��
=q����V�{���\)�\)��;��`�hs�녿-���F��Ͽ9X�z����?}��+��P�b��u����������dZ�(����/�푿/�p��p���5?�vɿ�R�;d�|��w� �� Ĝ�!%�!�7�"J�"�\�"�\�"��#S��$Z�$�/�%`B�%`B�%�˿&$ݿ&��'l��'��(1'�(�9�)xտ)��*���+�+�+ƨ�,1�,�D�-V�-�h�.{�/\)�0bN�0bN�0bN�0�׿0�`�1hs�1녿2n��2�3t��3�F�3�Ͽ49X�49X�4z�4z�49X�49X�4z�4�j�4�j�4���5?}�5?}�5?}�5��5��6E��6�+�6ȴ�6ȴ�7
=�7Kǿ7�P�7�ٿ8b�8Q�8�u�9��9X�9���:��:^5�:���;"ѿ;dZ�;�m�<(��<j�<��<푿=/�=/�=/�=p��=�-�=�-�=�-�>5?�>5?�>vɿ>�R�>�R�>�ۿ?;d�?�w�?�w�@  �@  �@A��@A��@A��@A��@�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B�B.B.B �B/B�B�B7LB�B>wBr�B�B8RB�;B?}BDBJBhB)�Bn�BƨB+B��B�B`BB�DBɺB�B	B	
=B	�B	�B	1'B	VB	t�B	�PB	�bB	��B	��B	�B	�B	�'B	�LB	�jB	�}B	��B	�)B	�5B	�`B	�B	�B	��B

=B
�B
$�B
/B
2-B
B�B
M�B
W
B
e`B
ffB
n�B
t�B
�B
�JB
�{B
��B
�'B
�^B
��B
�
B
�)B
�B
�
B
�#B
�)B
�/B
�BB
�NB
�mB
��B
��BBPB�B8RB:^B>wBM�BN�BS�Bp�Bl�Bo�Bs�B{�B�1B�7B�7B}�Bz�B}�B~�B�B�%B�JB�PB�VB�oB��B��B��B��B��B�B�B�3B�dB��BƨB��B��B�/B�)B�5B�)B�TB�5B�)B�/B�/B�/B�#B�#B�#B�BB�fB�fB�fB�fB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�B�B�B�B�B�B�B�B�B�B�fB�sB�ZB�ZB�ZB�NB�`B�`B�ZB�yB�BB�;B�;B�5B�/B�)B�5B�#B�/B�B�#B�B�#B�#B�)B�#B�)B�#B�B�B�
B�B��B�B��B��B��B��B��B��B�B�#B�#B�/B�5B�)B�#B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��BǮBƨBB��B��B�wB�wB�wB�wB�wB�qB�dB�RB�^B�RB�FB�FB�XB�XB�XB�XB�LB�RB�dB��B��B��B��B�wB��B�^B�^B�RB�RB�9B�3B�3B�-B�B�B�B��B�B��B��B��B��B��B��B��B�B�!B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�[BB1�B1�B$QB2�BB#JB:�BBBBv<BݩB;�B��BC	B�B�B�B-�Br$B�4B
�B�B,Bc�B��B�FB�B	�B	�B	2B	"DB	4�B	Y�B	xHB	��B	��B	�B	�QB	��B	��B	��B	��B	��B	�	B	؊B	ߵB	��B	��B	�$B	�<B	�gB
�B
B
(iB
2�B
5�B
FB
Q_B
Z�B
h�B
i�B
r$B
xHB
��B
��B
�B
�B
��B
��B
ׄB
ږB
ߵB
ۜB
ږB
ޯB
ߵB
�B
��B
��B
��B
�UB
�gB�B�B&B;�B=�BBBQ_BReBW�Bt0BpBs*BwBBsB��B��B��B��B~mB��B��B��B��B��B��B��B��B�B�B�DB�JB�B��B��B��B��B�B�4B�MB�~B�BߵB��BߵB��B��BߵB�B�B�BޯBޯBޯB��B��B��B��B��B�B�B�B�B�B�B�$B�$B�*B�*B�*B�6B�0B�*B�B�B�B�B�B�B�B�*B�*B�*B�*B�$B�B��B��B��B��B��B��B��B��B��B�B��B��B��B��B�BߵB��BޯB�BݩBޯBݩBޯBޯBߵBޯBߵBޯBܣBݩBږBۜB؊BِBׄB�~B�~B�~B�~BׄBܣBޯBޯB�B��BߵBޯBِBِB�qB�xB�xB�qB�xB�qB�xB�qB�xB�qB�~B؊B�~B�~B�xB�xB�~BׄBׄB�~B�~B�~B�xB�kB�qB�qB�qB�kB�eB�_B�eB؊BِBِB؊BׄB�~B�~B�xB�qB�qB�eB�_B�MB�:B�4B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�oB�oB�oB�iB�uB�|B��B��B��B��B��B��B��B�uB�oB�oB��B��B�|B��B�|B�uB�uB�uB�uB�oB�oB�iB�iB�iB�cB�cB�]B�WB�WB�QB�JB�DB�DB�DB�JB�DB�DB�>B�>B�>B�>B�>B�>B�>B�>B�8B�>B�>B�>B�>B�>B�>B�>B�>B�>B�DB�>B�8B�8B�>B�2B�8B�2B�2B�2B�2B�2B�2B�8B�>B�DB�DB�>B�DB�DB�DB�>B�DB�DB�DB�JB�DB�DB�DB�DB�DB�DB�JB�DB�JB�JB�JB�JB�QB�JB�JB�JB�QB�QB�JB�QB�QB�QB�QB�QB�QB�WB�WB�WB�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�iB�cB�iB�iB�iB�iB�iB�iB�oB�oB�iB�iB�iB�oB�oB�oB�oB�oB�oB�uB�uB�oB�uB�oB�oB�uB�uB�oB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�uB�uB�|B�uB�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B�|B�|B�|B�|B�|B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542382021060715423820210607154238  IF  ARFMCODA029d                                                                20190620115649                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115706  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115706  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154238  IP  PSAL            @���D��G�O�                