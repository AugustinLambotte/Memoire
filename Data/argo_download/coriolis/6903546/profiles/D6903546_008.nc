CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  V   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2018-08-09T22:09:51Z creation; 2023-09-02T14:39:23Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8    DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8    PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8(   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8h   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     9   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9,   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9L   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9l   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9p   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9x   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9|   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
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
_FillValue                 X  D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	X  F`   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	X  R   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  [h   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  g   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  pp   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  |    PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  �x   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  �(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	X  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �4   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �8   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �<   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �@   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �D   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �Argo profile    3.1 1.2 19500101000000  20180809220951  20230902143923  6903546 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU003                  5900A04                         844 @�xv����1   @�xw����@R��;rj����1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.0 dbar]                                                      @ff@Fff@�33@�  @���@�  A33A  A!��A333AA��ANffAa��Aq��A�  A���A���A���A���A�ffA�33A�33A�  A�33A�ffA�  A���A�ffA홚A�33B ffB��B  BffBffBffB��B��B��B$  B'��B+33B0  B4  B7��B<  B@  BDffBG33BK33BO��BT  BX��B[��B^��Bc33Bh  Bk33Bo��BtffBw��Bz��B33B�  B���B�  B�33B���B�ffB���B�  B���B���B�33B���B���B�33B���B���B�  B�ffB���B���B�  B�33B���B�  B�33B���B���B�33B���B���B�  B���B���B�  Bř�B���B���B�ffB�  B�  B���B�33B�33B���BꙚB�ffB���B���B�ffB���CL�C33C� CffC	L�C��C��C� CffCffCffCffCffCffCffCffC!ffC#ffC%� C'� C)��C+L�C-L�C/ffC1� C3��C5L�C7� C9��C;L�C=� C?��CAL�CC� CE�3CGffCI33CKffCM��COffCQ�CSffCU��CWffCY33C[� C]�3C_ffCa33Cc� Ce��CgffCi33Ck� Cm�3Co� CqL�Cs�CuffCw��CyffC{33C}  CffC���C��3C���C�� C�ٚC�� C��3C���C���C��3C�ٚC���C�� C��3C��fC���C���C�� C��3C��fC��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C���C���C���C�ٚC�ٚC��fC��fC��fC��3C��3C��3C�� C�� C���C���C���C���C��fC���C��fC��fC��fC��fC��3C��3C��3C�� C�� C���C���C��fC��3C��3C�� C���C¦fCó3C�� C�ٚCƳ3C�� C�ٚCɳ3Cʙ�C˳3C�� CͦfC�� C�ٚCг3Cь�Cҳ3C���CԦfCՌ�Cֳ3C���Cس3Cٙ�Cڌ�C۳3C���C�� CަfCߌ�C�� C�ٚC���C�3C䙚C��C�3C�ٚC���C�� C�fC뙚C��C�� C��3C��fC�ٚC���C�� C�3C��3C��3C��fC��fC��fC��fC�Y�C�ٚD &fDffD�fD� D&fDffD��D�3D
9�D� D��D��D  D� D� D��D33Dy�D��DfD33DffD�3DfD9�Dl�D �fD"fD#33D$y�D%�fD&�3D(&fD)y�D*�3D,�D-L�D.��D/�3D1�D2S3D3�3D4��D5� D79�D8s3D9� D;�D<FfD=� D>��D@  DAL�DBy�DC�fDD��DF9�DG��DH�fDJ  DK@ DL�fDM�fDN�fDP,�DQs3DR��DTfDU33DVffDW�3DYfDZ9�D[l�D\� D]��D_Y�D`�3Da��DcfDdFfDe�fDf� DhfDiL�Djy�Dk�fDl��Dn33Doy�Dp�fDr  Ds33Dtl�Du��Dv��Dx33Dy� Dz�3D{��D}@ D~� D�3D�y�D�fD���D�Y�D���D��3D�FfD�� D�vfD��D��fD�` D���D���D�6fD��3D�p D� D�� D�P D�� D��3D�33D��fD�|�D�&fD��3D�\�D���D��fD�6fD��fD�y�D��D�� D�c3D�fD�� D�9�D��3D�� D��D���D�\�D�  D��fD�<�D��3D�� D��D��fD�VfD��fD���D�9�D���D�� D��D���D�` D���D��fD�33D��3D�vfD��D���D�\�D�3D���D�6fD�� D�|�D��D���D�Y�D���D���D�@ D��3D�vfD��D�� D�ffD���D��3D�9�D�� D�y�D�3D���D�i�D�fD��fD�C3D��3D��3D�#3D��3D�c3D�3D��fD�I�D���D�|�D� D¹�D�VfD�  DĦfD�C3D�� D�|�D��DǼ�D�\�D�  Dɣ3D�<�D��fDˀ D��D̹�D�Y�D��fDΣ3D�FfD���D�s3D��D��fD�` D���DӜ�D�9�D�ٚD�vfD��D�� D�Y�D���Dؙ�D�9�D�� D�y�D�fD۶fD�VfD�3Dݣ3D�C3D��fD�y�D�  D��fD�\�D��3D��D�C3D�ٚD�s3D��D�ɚD�ffD�fD�3D�C3D��fD鉚D�,�D��D�\�D���D� D�FfD�� D� D�  D�� D�` D�3D�3D�FfD���D�s3D��D�� D�Y�D��fD�� D�<�D�� D�� D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@Fff@�33@�  @���@�  A33A  A!��A333AA��ANffAa��Aq��A�  A���A���A���A���A�ffA�33A�33A�  A�33A�ffA�  A���A�ffA홚A�33B ffB��B  BffBffBffB��B��B��B$  B'��B+33B0  B4  B7��B<  B@  BDffBG33BK33BO��BT  BX��B[��B^��Bc33Bh  Bk33Bo��BtffBw��Bz��B33B�  B���B�  B�33B���B�ffB���B�  B���B���B�33B���B���B�33B���B���B�  B�ffB���B���B�  B�33B���B�  B�33B���B���B�33B���B���B�  B���B���B�  Bř�B���B���B�ffB�  B�  B���B�33B�33B���BꙚB�ffB���B���B�ffB���CL�C33C� CffC	L�C��C��C� CffCffCffCffCffCffCffCffC!ffC#ffC%� C'� C)��C+L�C-L�C/ffC1� C3��C5L�C7� C9��C;L�C=� C?��CAL�CC� CE�3CGffCI33CKffCM��COffCQ�CSffCU��CWffCY33C[� C]�3C_ffCa33Cc� Ce��CgffCi33Ck� Cm�3Co� CqL�Cs�CuffCw��CyffC{33C}  CffC���C��3C���C�� C�ٚC�� C��3C���C���C��3C�ٚC���C�� C��3C��fC���C���C�� C��3C��fC��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C�� C���C���C���C�ٚC�ٚC��fC��fC��fC��3C��3C��3C�� C�� C���C���C���C���C��fC���C��fC��fC��fC��fC��3C��3C��3C�� C�� C���C���C��fC��3C��3C�� C���C¦fCó3C�� C�ٚCƳ3C�� C�ٚCɳ3Cʙ�C˳3C�� CͦfC�� C�ٚCг3Cь�Cҳ3C���CԦfCՌ�Cֳ3C���Cس3Cٙ�Cڌ�C۳3C���C�� CަfCߌ�C�� C�ٚC���C�3C䙚C��C�3C�ٚC���C�� C�fC뙚C��C�� C��3C��fC�ٚC���C�� C�3C��3C��3C��fC��fC��fC��fC�Y�C�ٚD &fDffD�fD� D&fDffD��D�3D
9�D� D��D��D  D� D� D��D33Dy�D��DfD33DffD�3DfD9�Dl�D �fD"fD#33D$y�D%�fD&�3D(&fD)y�D*�3D,�D-L�D.��D/�3D1�D2S3D3�3D4��D5� D79�D8s3D9� D;�D<FfD=� D>��D@  DAL�DBy�DC�fDD��DF9�DG��DH�fDJ  DK@ DL�fDM�fDN�fDP,�DQs3DR��DTfDU33DVffDW�3DYfDZ9�D[l�D\� D]��D_Y�D`�3Da��DcfDdFfDe�fDf� DhfDiL�Djy�Dk�fDl��Dn33Doy�Dp�fDr  Ds33Dtl�Du��Dv��Dx33Dy� Dz�3D{��D}@ D~� D�3D�y�D�fD���D�Y�D���D��3D�FfD�� D�vfD��D��fD�` D���D���D�6fD��3D�p D� D�� D�P D�� D��3D�33D��fD�|�D�&fD��3D�\�D���D��fD�6fD��fD�y�D��D�� D�c3D�fD�� D�9�D��3D�� D��D���D�\�D�  D��fD�<�D��3D�� D��D��fD�VfD��fD���D�9�D���D�� D��D���D�` D���D��fD�33D��3D�vfD��D���D�\�D�3D���D�6fD�� D�|�D��D���D�Y�D���D���D�@ D��3D�vfD��D�� D�ffD���D��3D�9�D�� D�y�D�3D���D�i�D�fD��fD�C3D��3D��3D�#3D��3D�c3D�3D��fD�I�D���D�|�D� D¹�D�VfD�  DĦfD�C3D�� D�|�D��DǼ�D�\�D�  Dɣ3D�<�D��fDˀ D��D̹�D�Y�D��fDΣ3D�FfD���D�s3D��D��fD�` D���DӜ�D�9�D�ٚD�vfD��D�� D�Y�D���Dؙ�D�9�D�� D�y�D�fD۶fD�VfD�3Dݣ3D�C3D��fD�y�D�  D��fD�\�D��3D��D�C3D�ٚD�s3D��D�ɚD�ffD�fD�3D�C3D��fD鉚D�,�D��D�\�D���D� D�FfD�� D� D�  D�� D�` D�3D�3D�FfD���D�s3D��D�� D�Y�D��fD�� D�<�D�� D�� D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�7L@�7L@�&�@ӶF@��@ϥ�@�o@��y@�t�@�X@���@�x�@��@��@��`@��@�^5@��-@��@���@�I�@�@��P@�9X@y��@x  @rn�@h�u@_l�@Z�@Yhs@S�@J��@;@0��@!��@�@j@	�@1@ r�?��?��T?��j?���?��m?�=q?��?��`?���?�l�?У�?�I�?��H?�Q�?�Z?� �?��R?���?���?�?���?���?��F?� �?�Z?��?���?��?͑h?���?��?�33?Ұ!?��;?̬?Ƈ+?��?��?���?���?�Q�?���?�I�?�9X?�v�?���?���?��!?�t�?��j?�z�?�V?�Ĝ?�t�?�J?��?��h?��7?��R?�7L?�1?���?�?�o?�G�?��?�/?�V?�V?��H?�ff?�ff?���?�o?�G�?|�?v�+?t�j?v�+?pbN?g�?f$�?`Ĝ?Z�?Y�?X��?U?}?P�`?N{?Ix�?F��?B��?;dZ?;"�?8Q�?6ȴ?6ȴ?7K�??|�?E�?<(�?7��?=�?C�
?<(�?5?}?*~�?&$�?$Z?|�???}?X?��??}?j?�-?�H?X?b?
=?ȴ?j?�-?j?�H?��?�j?��?C�?C�?
=q?�>�^5>��>�ff>�`B>�V>�&�>�1>߾w>�bN>�ƨ>�z�>޸R>�/>�
=>���>�hs>���>�X>�K�>�ȴ>��>�Ĝ>�Ĝ>�A�>��>��>��->��->��->�/>���>���>��>���>��>��u>���>��>��>�\)>���>���>���>��^>�$�>��>��>���>��>�1'>�7L>�1'>��>���>���>��7>�%>~��>}�>{�m>y�#>vȴ>s�F>q��>n��>hr�>fff>fff>e`B>dZ>M��>F��>F��>F��>E��>G�>@�><j>8Q�>7K�>9X>:^5>;dZ>:^5>:^5>:^5>;dZ><j>>v�><j>:^5><j>@�>F��>H�9>H�9>A�7>@�>>v�>:^5>;dZ>9X>6E�>2->1&�>-V>,1>,1>,1>-V>!��>��>�>�+>�u>��>�u>��>�>�->�u>z�>n�>bN>I�>$�=��#=�F=�=�"�=ȴ9=�{=��-=�7L=T��=�w<��
<D��;��
    ��o�u��j�����t��@��y�#��C����㽧�E��ě�����"ѽ�/��`B����#�J�+�C��O߾I��\)����u�����R�#�
�%�T�'+�0 ž333�49X�6E��;dZ�?|�C���H�9�J���L�;P�`�S�ϾV�Xb�Z��]/�_;d�`A��aG��cS��dZ�fff�fff�gl��ixվk��l�D�m�h�p�׾s�F�u�w�پy�#�z�H�|푾�  ��%��J��o���������$ݾ���+������9��7L���^����������ƨ��I����;�O߾�V�����\)���;��bN��hs��n���񪾓t����Ͼ�zᾔ����������+���+���+��
=���P��b���u�������������"Ѿ�"Ѿ��㾛�㾜(������������/���-��5?���R��;d��;d���w��A���A���A���G��������徢�徣�
��Z��`B���y���r���xվ�~������V��{����� ž�&龳�F���j����ȴ��KǾ�Q쾹�#���m��p���vɾ�|����\�Õ���������$ݾ�1'��7L������ƨ��O߾����bN����z��
=�ؓu�ڟ��ۥ�ܬ��;d��G������S���Z��`B��ff����þ����1��V��h������ ž�׾�-��F��9X��?}��E���KǾ��پ�X���#���H���m��푾�vɾ�|� ��G��J�������
��/�`B��T�$ݿ��l��r��	�^�
=q�
=q�
=q�
���C����1�I���D��D��ͿO߿�h���{�V������\)�\)�����;��;��׿&�hs����n���33��F��Ͽ9X��j����?}��E��ȴ�
=��P��ٿ�ٿb�Q��u�X�����#�^5�"ѿdZ������m�j���p��5?�v�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�7L@�7L@�&�@ӶF@��@ϥ�@�o@��y@�t�@�X@���@�x�@��@��@��`@��@�^5@��-@��@���@�I�@�@��P@�9X@y��@x  @rn�@h�u@_l�@Z�@Yhs@S�@J��@;@0��@!��@�@j@	�@1@ r�?��?��T?��j?���?��m?�=q?��?��`?���?�l�?У�?�I�?��H?�Q�?�Z?� �?��R?���?���?�?���?���?��F?� �?�Z?��?���?��?͑h?���?��?�33?Ұ!?��;?̬?Ƈ+?��?��?���?���?�Q�?���?�I�?�9X?�v�?���?���?��!?�t�?��j?�z�?�V?�Ĝ?�t�?�J?��?��h?��7?��R?�7L?�1?���?�?�o?�G�?��?�/?�V?�V?��H?�ff?�ff?���?�o?�G�?|�?v�+?t�j?v�+?pbN?g�?f$�?`Ĝ?Z�?Y�?X��?U?}?P�`?N{?Ix�?F��?B��?;dZ?;"�?8Q�?6ȴ?6ȴ?7K�??|�?E�?<(�?7��?=�?C�
?<(�?5?}?*~�?&$�?$Z?|�???}?X?��??}?j?�-?�H?X?b?
=?ȴ?j?�-?j?�H?��?�j?��?C�?C�?
=q?�>�^5>��>�ff>�`B>�V>�&�>�1>߾w>�bN>�ƨ>�z�>޸R>�/>�
=>���>�hs>���>�X>�K�>�ȴ>��>�Ĝ>�Ĝ>�A�>��>��>��->��->��->�/>���>���>��>���>��>��u>���>��>��>�\)>���>���>���>��^>�$�>��>��>���>��>�1'>�7L>�1'>��>���>���>��7>�%>~��>}�>{�m>y�#>vȴ>s�F>q��>n��>hr�>fff>fff>e`B>dZ>M��>F��>F��>F��>E��>G�>@�><j>8Q�>7K�>9X>:^5>;dZ>:^5>:^5>:^5>;dZ><j>>v�><j>:^5><j>@�>F��>H�9>H�9>A�7>@�>>v�>:^5>;dZ>9X>6E�>2->1&�>-V>,1>,1>,1>-V>!��>��>�>�+>�u>��>�u>��>�>�->�u>z�>n�>bN>I�>$�=��#=�F=�=�"�=ȴ9=�{=��-=�7L=T��=�w<��
<D��;��
    ��o�u��j�����t��@��y�#��C����㽧�E��ě�����"ѽ�/��`B����#�J�+�C��O߾I��\)����u�����R�#�
�%�T�'+�0 ž333�49X�6E��;dZ�?|�C���H�9�J���L�;P�`�S�ϾV�Xb�Z��]/�_;d�`A��aG��cS��dZ�fff�fff�gl��ixվk��l�D�m�h�p�׾s�F�u�w�پy�#�z�H�|푾�  ��%��J��o���������$ݾ���+������9��7L���^����������ƨ��I����;�O߾�V�����\)���;��bN��hs��n���񪾓t����Ͼ�zᾔ����������+���+���+��
=���P��b���u�������������"Ѿ�"Ѿ��㾛�㾜(������������/���-��5?���R��;d��;d���w��A���A���A���G��������徢�徣�
��Z��`B���y���r���xվ�~������V��{����� ž�&龳�F���j����ȴ��KǾ�Q쾹�#���m��p���vɾ�|����\�Õ���������$ݾ�1'��7L������ƨ��O߾����bN����z��
=�ؓu�ڟ��ۥ�ܬ��;d��G������S���Z��`B��ff����þ����1��V��h������ ž�׾�-��F��9X��?}��E���KǾ��پ�X���#���H���m��푾�vɾ�|� ��G��J�������
��/�`B��T�$ݿ��l��r��	�^�
=q�
=q�
=q�
���C����1�I���D��D��ͿO߿�h���{�V������\)�\)�����;��;��׿&�hs����n���33��F��Ͽ9X��j����?}��E��ȴ�
=��P��ٿ�ٿb�Q��u�X�����#�^5�"ѿdZ������m�j���p��5?�v�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�fB�fB�B	�!B
T�B
gmB
�7B
��B
��B
�TB
�sB
�B
�B
�B
��B
��BB�B>wB33B5?B:^B9XB?}BG�BL�BO�BP�BR�BT�BS�BJ�BS�B[#B\)BhsBk�Bl�BffBl�Bq�Bs�Bu�Bw�B{�B}�B~�B�B�B�B�+B�7B�7B�PB�DB�VB�DB�PB�VB�hB��B��B��B��B�{B��B�JB�B�B�LB�jB�XB��B�qB��B�wB�qB�qB�jB�qB�dB�qB�RB�FB�^B�?B�B�RB�dB�XB�qB�qBȴB�dB�qB�}B�wB�}B�jB�}B��B��B�}B�wB�qB�jB�jB�jB�wB��B��B�dB�jB�jB�qB�jB�jB�jB�}B�jB�jB�^B�XB�XB�RB�RB�^B�RB�XB�LB�?B�?B�?B�9B�3B�9B�3B�3B�LB�?B�XB�LB�RB�LB�dB�XB�RB�RB�?B�?B�FB�9B�'B�FB�9B�3B�RB�XB�RB�XB�XB�XB�XB�dB�jB�^B�^B�^B�^B�XB�RB�RB�RB�RB�?B�9B�-B�9B�?B�?B�9B�^B�-B�-B�3B�?B�?B�?B�9B�9B�FB�'B�'B�!B�B�B�B�B�B�B�!B�!B�!B�!B�!B�!B�'B�'B�'B�'B�'B�'B�'B�!B�'B�!B�'B�'B�'B�!B�!B�!B�'B�-B�-B�-B�-B�'B�-B�'B�'B�-B�'B�-B�-B�-B�-B�'B�'B�'B�'B�'B�-B�!B�-B�!B�!B�!B�B�!B�!B�!B�B�B�!B�!B�!B�B�!B�!B�!B�'B�!B�'B�'B�!B�-B�-B�3B�-B�'B�-B�'B�-B�!B�-B�3B�-B�!B�!B�!B�'B�'B�'B�'B�'B�!B�!B�!B�'B�!B�'B�'B�'B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�fB�fB�B	�!B
T�B
gmB
�7B
��B
��B
�TB
�sB
�B
�B
�B
��B
��BB�B>wB33B5?B:^B9XB?}BG�BL�BO�BP�BR�BT�BS�BJ�BS�B[#B\)BhsBk�Bl�BffBl�Bq�Bs�Bu�Bw�B{�B}�B~�B�B�B�B�+B�7B�7B�PB�DB�VB�DB�PB�VB�hB��B��B��B��B�{B��B�JB�B�B�LB�jB�XB��B�qB��B�wB�qB�qB�jB�qB�dB�qB�RB�FB�^B�?B�B�RB�dB�XB�qB�qBȴB�dB�qB�}B�wB�}B�jB�}B��B��B�}B�wB�qB�jB�jB�jB�wB��B��B�dB�jB�jB�qB�jB�jB�jB�}B�jB�jB�^B�XB�XB�RB�RB�^B�RB�XB�LB�?B�?B�?B�9B�3B�9B�3B�3B�LB�?B�XB�LB�RB�LB�dB�XB�RB�RB�?B�?B�FB�9B�'B�FB�9B�3B�RB�XB�RB�XB�XB�XB�XB�dB�jB�^B�^B�^B�^B�XB�RB�RB�RB�RB�?B�9B�-B�9B�?B�?B�9B�^B�-B�-B�3B�?B�?B�?B�9B�9B�FB�'B�'B�!B�B�B�B�B�B�B�!B�!B�!B�!B�!B�!B�'B�'B�'B�'B�'B�'B�'B�!B�'B�!B�'B�'B�'B�!B�!B�!B�'B�-B�-B�-B�-B�'B�-B�'B�'B�-B�'B�-B�-B�-B�-B�'B�'B�'B�'B�'B�-B�!B�-B�!B�!B�!B�B�!B�!B�!B�B�B�!B�!B�!B�B�!B�!B�!B�'B�!B�'B�'B�!B�-B�-B�3B�-B�'B�-B�'B�-B�!B�-B�3B�-B�!B�!B�!B�'B�'B�'B�'B�'B�!B�!B�!B�'B�!B�'B�'B�'B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202309021439232023090214392320230902143923  IF  ARFMCODA023a                                                                20180809220951                      G�O�G�O�G�O�                IF  ARGQCOQC3.4                                                                 20180809220954  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.4                                                                 20180809220954  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology 20190821161306  IP  PSAL            @ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607145130  IP  PSAL            @ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20220330173425  IP  PSAL            @ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230114155654  IP  PSAL            @ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230902143923  IP  PSAL            @ffD��3G�O�                