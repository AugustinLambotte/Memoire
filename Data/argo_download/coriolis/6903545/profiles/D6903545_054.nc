CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  P   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:29Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        	@  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	@  Fh   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	@  Q�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  [8   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  dx   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  f�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  p   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  rX   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  {�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  �(   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �h   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	@  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �T   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �X   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �\   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �`   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �d   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �(   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �(   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �(   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �(Argo profile    3.1 1.2 19500101000000  20200829092229  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               6A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��v-��.1   @��v-��.@R�Ý���%ƿ�K��8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A0  A;33AL��A`  AnffA|��A�ffA�ffA�  A�  A���A�  A���A�33A���A�33Aՙ�A�33A���A�33A�ffB   B33BffB��B33BffB��B��B   B#��B&��B,  B0��B4ffB7��B;33B>��BC��BH��BLffBP  BS��BW33BZ��B^ffBc��Bi33Bl��BpffBtffBx  B|  B��B���B���B�ffB�ffB�33B�33B�  B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�Bę�Bƙ�Bș�B˙�Bϙ�Bә�Bי�Bۙ�Bߙ�B㙚B癚B뙚B���B���B���B���B�  C  C  C�C�C	33CL�CffCffCffC� C��C�3C��C�fC� C  C!�C#33C%33C'L�C)L�C+L�C-L�C/L�C1ffC3ffC5� C7� C9��C;33C=L�C?ffCAffCC� CE��CG�3CIffCK� CM��COL�CQffCS��CUL�CWffCY��C[L�C]ffC_� Ca33CcffCe� CgL�Ci� Ck�3Cm� CoL�Cq�CsL�Cu��CwffCy33C{  C}ffC��C�� C��fC���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C���C���C��3C��3C��fC��fC�ٚC���C���C���C�� C�� C��3C��fC��fC��fC���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C�� C�� C�s3C�� C�� C�� C�� C�� C���C���C���C���C���C���C���C���C�CÌ�Cę�Cř�Cƙ�CǦfCȦfCɦfCʦfC˳3C̳3CͦfCγ3Cϳ3Cг3Cѳ3Cҳ3Cӳ3CԳ3C�� C�� C�� C���C���Cڙ�CۦfCܳ3Cݳ3C�� C�� C���C�fC�fC�� C���C�ٚC�3C�� C�ٚC�3C��C�fC�� C��fC��C�fC�� C�fC��C�� C�ٚC�� C��3C���C�� C��3C���C��3D ,�D` D�3D�D@ Dy�D�3D��D
,�Dl�D��D��D33Dy�D�fD3D@ D� D�fD��D,�DffD� D�fD,�Dy�D �fD!�3D#&fD$y�D%��D&�3D(FfD)�fD*�fD,fD-,�D.s3D/��D1  D2FfD3s3D4�fD5�3D7FfD8y�D9�3D:�fD<&fD=ffD>�fD?��DA33DBs3DC��DE  DFFfDG��DH�3DI��DK  DLs3DM�fDN�3DP  DQffDR��DT3DUS3DVy�DW�fDX�fDZ,�D[y�D\��D^  D_9�D`y�Da��Dc  Dd&fDes3Df��Dg�fDi9�Dj��Dk� Dl��Dn33Dol�Dp��Dq��Ds&fDtffDu��Dv��Dx,�Dys3Dz��D{�3D}&fD~` D��D��fD�&fD��3D�c3D�3D��fD�<�D��3D�y�D�#3D�� D�\�D���D��fD�6fD�ٚD�y�D��D�� D�VfD���D��3D�<�D��fD�s3D�3D�� D�S3D��3D��fD�9�D���D�y�D�  D�ɚD�c3D�3D��3D�FfD���D�s3D��D�� D�\�D���D���D�6fD�ٚD�|�D��D��3D�Y�D���D��fD�33D�� D�p D� D�� D�P D��3D���D�<�D��3D���D��D�� D�` D���D���D�6fD��3D�s3D�3D���D�\�D�� D��3D�6fD���D��3D��D��fD�\�D���D��fD�@ D�ٚD�vfD�#3D�� D�\�D�  D�� D�@ D�� D��3D�fD�� D�ffD�  D���D�6fD��fD�vfD�fD¼�D�\�D�  Dģ3D�FfD��fD�|�D� DǶfD�Y�D�  DɦfD�C3D�� Dˀ D�#3D̼�D�VfD�  DΙ�D�6fD�� D�|�D��DѼ�D�\�D�  DӖfD�<�D��3D�|�D�fDּ�D�i�D�3Dأ3D�C3D�� D�y�D��D��3D�` D���DݖfD�33D��3D�p D� D�3D�VfD���D�3D�<�D��fD�3D�#3D�� D�c3D���D��D�FfD���D�vfD�  D�ɚD�c3D�  D� D�<�D�ٚD�y�D��D��D�` D�fD�D�0 D���D�|�D��D��fD�VfD��fD��fD�9�D���D�� D��D��3D�ffD��f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A0  A;33AL��A`  AnffA|��A�ffA�ffA�  A�  A���A�  A���A�33A���A�33Aՙ�A�33A���A�33A�ffB   B33BffB��B33BffB��B��B   B#��B&��B,  B0��B4ffB7��B;33B>��BC��BH��BLffBP  BS��BW33BZ��B^ffBc��Bi33Bl��BpffBtffBx  B|  B��B���B���B�ffB�ffB�33B�33B�  B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�Bę�Bƙ�Bș�B˙�Bϙ�Bә�Bי�Bۙ�Bߙ�B㙚B癚B뙚B���B���B���B���B�  C  C  C�C�C	33CL�CffCffCffC� C��C�3C��C�fC� C  C!�C#33C%33C'L�C)L�C+L�C-L�C/L�C1ffC3ffC5� C7� C9��C;33C=L�C?ffCAffCC� CE��CG�3CIffCK� CM��COL�CQffCS��CUL�CWffCY��C[L�C]ffC_� Ca33CcffCe� CgL�Ci� Ck�3Cm� CoL�Cq�CsL�Cu��CwffCy33C{  C}ffC��C�� C��fC���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C���C���C��3C��3C��fC��fC�ٚC���C���C���C�� C�� C��3C��fC��fC��fC���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C�� C�� C�s3C�� C�� C�� C�� C�� C���C���C���C���C���C���C���C���C�CÌ�Cę�Cř�Cƙ�CǦfCȦfCɦfCʦfC˳3C̳3CͦfCγ3Cϳ3Cг3Cѳ3Cҳ3Cӳ3CԳ3C�� C�� C�� C���C���Cڙ�CۦfCܳ3Cݳ3C�� C�� C���C�fC�fC�� C���C�ٚC�3C�� C�ٚC�3C��C�fC�� C��fC��C�fC�� C�fC��C�� C�ٚC�� C��3C���C�� C��3C���C��3D ,�D` D�3D�D@ Dy�D�3D��D
,�Dl�D��D��D33Dy�D�fD3D@ D� D�fD��D,�DffD� D�fD,�Dy�D �fD!�3D#&fD$y�D%��D&�3D(FfD)�fD*�fD,fD-,�D.s3D/��D1  D2FfD3s3D4�fD5�3D7FfD8y�D9�3D:�fD<&fD=ffD>�fD?��DA33DBs3DC��DE  DFFfDG��DH�3DI��DK  DLs3DM�fDN�3DP  DQffDR��DT3DUS3DVy�DW�fDX�fDZ,�D[y�D\��D^  D_9�D`y�Da��Dc  Dd&fDes3Df��Dg�fDi9�Dj��Dk� Dl��Dn33Dol�Dp��Dq��Ds&fDtffDu��Dv��Dx,�Dys3Dz��D{�3D}&fD~` D��D��fD�&fD��3D�c3D�3D��fD�<�D��3D�y�D�#3D�� D�\�D���D��fD�6fD�ٚD�y�D��D�� D�VfD���D��3D�<�D��fD�s3D�3D�� D�S3D��3D��fD�9�D���D�y�D�  D�ɚD�c3D�3D��3D�FfD���D�s3D��D�� D�\�D���D���D�6fD�ٚD�|�D��D��3D�Y�D���D��fD�33D�� D�p D� D�� D�P D��3D���D�<�D��3D���D��D�� D�` D���D���D�6fD��3D�s3D�3D���D�\�D�� D��3D�6fD���D��3D��D��fD�\�D���D��fD�@ D�ٚD�vfD�#3D�� D�\�D�  D�� D�@ D�� D��3D�fD�� D�ffD�  D���D�6fD��fD�vfD�fD¼�D�\�D�  Dģ3D�FfD��fD�|�D� DǶfD�Y�D�  DɦfD�C3D�� Dˀ D�#3D̼�D�VfD�  DΙ�D�6fD�� D�|�D��DѼ�D�\�D�  DӖfD�<�D��3D�|�D�fDּ�D�i�D�3Dأ3D�C3D�� D�y�D��D��3D�` D���DݖfD�33D��3D�p D� D�3D�VfD���D�3D�<�D��fD�3D�#3D�� D�c3D���D��D�FfD���D�vfD�  D�ɚD�c3D�  D� D�<�D�ٚD�y�D��D��D�` D�fD�D�0 D���D�|�D��D��fD�VfD��fD��fD�9�D���D�� D��D��3D�ffD��f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����+��p��n�>�$�>��?$Z?�X?���?��?��;?��?���?�hs?��?�~�?��?��@ 1'@z�@{@�\@�7@�@�@`B@v�@ff@v�@O�@v�@  @��@
�@b@�^@l�@�@��@1@ff@ 1'@|�@|�@';d@&�R@8��@/�P@)��@(�9@)&�@/�P@1x�@3�@9%@:n�@:�H@:�@;�
@;�m@:n�@:��@:�@:~�@:-@9&�@8b@7\)@7�@6ȴ@6ff@4��@49X@3dZ@2�\@2=q@1�#@0��@/�w@/�;@/l�@/
=@/|�@/l�@.v�@.@-��@-/@-/@,�j@+dZ@*^5@'�@&�@%�h@#@;d@��@�H@�^@��@
=@�@��@��@�;@��@�/@�@
�@�P@@�@�j@O�@{@�D@-@ �?��R?���?�ƨ?���?��?�%?���?���?�+?���?��?��`?��H?Լj?�o?�G�?Ͼw?�dZ?Ƨ�?ě�?Ł?�Z?�M�?�I�?�ȴ?�$�?�t�?���?� �?�I�?���?�+?�?���?�t�?�-?�|�?�{?�p�?��D?�"�?���?��+?�b?�r�?�l�?�`B?��?���?���?��P?�ff?���?��?v�+?l��?i��?e`B?f$�?m��?q&�?p��?p��?n�?nV?l��?lI�?bJ?XQ�?T�j?T9X?Q��?P �?N��?J��?H1'?Fff?C�
?>v�?9�#?9X?:�?>��?AG�?A%??;d?=/?:�H?9X?9X?9X?9�#?8��?8��?7
=?4��?3�F?2-?1hs?1��?1�?1hs?0��?0��?/\)?,I�?/\)?.�?-��?*��?(r�?&�y?$�/?!%?;d?(�?�H?��?��?��?�m?"�?�#?�?X?��?�P?�P?�+?�j?�F?hs?��?��?��?
=q?+?�T?S�?   >���>�j>�j>�^5>�Q�>�33>�h>��
>׍P>��`>�C�>ě�>�J>�dZ>�K�>�>��j>���>�x�>�ff>�S�>���>�
=>��>�hs>��>�C�>�C�>�=q>�+>�J>x��>l�D>gl�>V>H�9>=p�>0 �>$�/>V>+>J=�=�=�x�=�"�=��`=\=�9X=��T=�t�=}�=aG�=T��=0 �=+=�P=t�<�h<���<���<49X<o;�o:�o<o;��
;o�49X�u�\)�8Q�T���ixսq���Y��}󶽋C���O߽��w�� Ž�vɽ�������;d��xս���%���	7L�O߾bN�z�������-�"��%�T�$�/�&�y�+�/��7KǾ:^5�=p��>vɾB�\�F��G��H�9�Kƨ�M��P�`�S�ϾT���Xb�["Ѿ^5?�`A��dZ�gl��j~��l�D�p�׾r�!�s�F�vȴ�x���}󶾀���J�����������𾇮��7L���^��C���I����;�V���`��n���t��������+��b��������"Ѿ�(�������-���R���w��Ĝ��Ĝ��G���Z���y���þ�xվ����D���D����� ž�&龳33��9X����KǾ�����^5��dZ���m��푾���  �����J�Õ��š˾�$ݾƧ��1'��=q��I���O߾�����;���`��t���z�և+��b�ٙ���"Ѿܬ�ݲ-��5?�߾w��G���MӾ��
���/���T��ff���y��r������~���D��h�����׾�-��33��9X����KǾ�X���H��j��p���vɾ�|� A��G�����Mӿ����������$ݿff����y�l��r��	7L�	��
~��
���C����1�I��V��������\)���� ſ�׿�`�hs����n��n���!�33�33��F��F��Ͽz�9X�������E���+�
=��ٿb��u������#�^5�����H�dZ����m���/��-��vɿ�R�;d�   �   � �� Ĝ� Ĝ�!G��!G��!���"Mӿ"Mӿ"�\�#o�#S��#���#�
�$Z�$�/�$�/�%��%�˿%�˿&ff�&ff�&ff�&�y�'+�'��(1'�(r��(�ÿ)7L�)�^�)��)��*~��*~��*���+C�1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111���+��p��n�>�$�>��?$Z?�X?���?��?��;G�O�G�O�?�hs?��?�~�?��?��@ 1'@z�@{@�\@�7@�@�@`B@v�@ff@v�@O�@v�@  @��@
�@b@�^@l�@�@��@1@ff@ 1'@|�@|�@';d@&�R@8��@/�P@)��@(�9@)&�@/�P@1x�@3�@9%@:n�@:�H@:�@;�
@;�m@:n�@:��@:�@:~�@:-@9&�@8b@7\)@7�@6ȴ@6ff@4��@49X@3dZ@2�\@2=q@1�#@0��@/�w@/�;@/l�@/
=@/|�@/l�@.v�@.@-��@-/@-/@,�j@+dZ@*^5@'�@&�@%�h@#@;d@��@�H@�^@��@
=@�@��@��@�;@��@�/@�@
�@�P@@�@�j@O�@{@�D@-@ �?��R?���?�ƨ?���?��?�%?���?���?�+?���?��?��`?��H?Լj?�o?�G�?Ͼw?�dZ?Ƨ�?ě�?Ł?�Z?�M�?�I�?�ȴ?�$�?�t�?���?� �?�I�?���?�+?�?���?�t�?�-?�|�?�{?�p�?��D?�"�?���?��+?�b?�r�?�l�?�`B?��?���?���?��P?�ff?���?��?v�+?l��?i��?e`B?f$�?m��?q&�?p��?p��?n�?nV?l��?lI�?bJ?XQ�?T�j?T9X?Q��?P �?N��?J��?H1'?Fff?C�
?>v�?9�#?9X?:�?>��?AG�?A%??;d?=/?:�H?9X?9X?9X?9�#?8��?8��?7
=?4��?3�F?2-?1hs?1��?1�?1hs?0��?0��?/\)?,I�?/\)?.�?-��?*��?(r�?&�y?$�/?!%?;d?(�?�H?��?��?��?�m?"�?�#?�?X?��?�P?�P?�+?�j?�F?hs?��?��?��?
=q?+?�T?S�?   >���>�j>�j>�^5>�Q�>�33>�h>��
>׍P>��`>�C�>ě�>�J>�dZ>�K�>�>��j>���>�x�>�ff>�S�>���>�
=>��>�hs>��>�C�>�C�>�=q>�+>�J>x��>l�D>gl�>V>H�9>=p�>0 �>$�/>V>+>J=�=�=�x�=�"�=��`=\=�9X=��T=�t�=}�=aG�=T��=0 �=+=�P=t�<�h<���<���<49X<o;�o:�o<o;��
;o�49X�u�\)�8Q�T���ixսq���Y��}󶽋C���O߽��w�� Ž�vɽ�������;d��xս���%���	7L�O߾bN�z�������-�"��%�T�$�/�&�y�+�/��7KǾ:^5�=p��>vɾB�\�F��G��H�9�Kƨ�M��P�`�S�ϾT���Xb�["Ѿ^5?�`A��dZ�gl��j~��l�D�p�׾r�!�s�F�vȴ�x���}󶾀���J�����������𾇮��7L���^��C���I����;�V���`��n���t��������+��b��������"Ѿ�(�������-���R���w��Ĝ��Ĝ��G���Z���y���þ�xվ����D���D����� ž�&龳33��9X����KǾ�����^5��dZ���m��푾���  �����J�Õ��š˾�$ݾƧ��1'��=q��I���O߾�����;���`��t���z�և+��b�ٙ���"Ѿܬ�ݲ-��5?�߾w��G���MӾ��
���/���T��ff���y��r������~���D��h�����׾�-��33��9X����KǾ�X���H��j��p���vɾ�|� A��G�����Mӿ����������$ݿff����y�l��r��	7L�	��
~��
���C����1�I��V��������\)���� ſ�׿�`�hs����n��n���!�33�33��F��F��Ͽz�9X�������E���+�
=��ٿb��u������#�^5�����H�dZ����m���/��-��vɿ�R�;d�   �   � �� Ĝ� Ĝ�!G��!G��!���"Mӿ"Mӿ"�\�#o�#S��#���#�
�$Z�$�/�$�/�%��%�˿%�˿&ff�&ff�&ff�&�y�'+�'��(1'�(r��(�ÿ)7L�)�^�)��)��*~��*~��*���+C�1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB~�B�5B	dZB	��B	�XB	�fB
%�B
6FB
e`B
l�B
v�B
T�B
�oB
��B
ŢB
��B
�B
�NB
��B1B
=BJBuB{B�B�B�B&�B(�B,B+B5?B7LB?}BN�B`BB\)B_;BgmBjBu�Bx�B�B�hB��B�'B��B�FB�'B�?B�jB��B�?B�B�BB�TB�ZB�yB�sB�B�sB�B�B�B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�fB�`B�`B�ZB�ZB�NB�HB�B�5B�B�B�
B�
B��B��B��B��B��B��B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBŢBŢBĜBÖBĜB��BÖB��B��B��B�wB�jB�jB�dB�dB�^B�^B�XB�XB�9B�?B�LB�FB�LB�FB�?B�FB�LB�LB�LB�RB�^B�RB�XB�XB�^B�?B�9B�9B�3B�'B�'B�B�B�B�B�'B�-B�3B�3B�?B�?B�?B�FB�9B�9B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�B�3B�3B�3B�9B�3B�-B�3B�3B�9B�9B�9B�9B�?B�9B�9B�?B�FB�FB�FB�LB�LB�LB�FB�RB�?B�LB�LB�LB�LB�FB�?B�?B�FB�9B�9B�?B�?B�FB�?B�FB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�LB�LB�LB�FB�FB�?B�FB�FB�?B�9B�?B�?B�9B�3B�3B�'B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B~�B�5B	dZB	��B	�XB	�fB
%�B
6FB
e`B
l�G�O�G�O�B
�oB
��B
ŢB
��B
�B
�NB
��B1B
=BJBuB{B�B�B�B&�B(�B,B+B5?B7LB?}BN�B`BB\)B_;BgmBjBu�Bx�B�B�hB��B�'B��B�FB�'B�?B�jB��B�?B�B�BB�TB�ZB�yB�sB�B�sB�B�B�B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�fB�`B�`B�ZB�ZB�NB�HB�B�5B�B�B�
B�
B��B��B��B��B��B��B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBŢBŢBĜBÖBĜB��BÖB��B��B��B�wB�jB�jB�dB�dB�^B�^B�XB�XB�9B�?B�LB�FB�LB�FB�?B�FB�LB�LB�LB�RB�^B�RB�XB�XB�^B�?B�9B�9B�3B�'B�'B�B�B�B�B�'B�-B�3B�3B�?B�?B�?B�FB�9B�9B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�B�3B�3B�3B�9B�3B�-B�3B�3B�9B�9B�9B�9B�?B�9B�9B�?B�FB�FB�FB�LB�LB�LB�FB�RB�?B�LB�LB�LB�LB�FB�?B�?B�FB�9B�9B�?B�?B�FB�?B�FB�LB�LB�RB�RB�RB�RB�RB�RB�RB�RB�LB�LB�LB�FB�FB�?B�FB�FB�?B�9B�?B�?B�9B�3B�3B�'B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092229                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092320  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092320  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A0  D��fG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172540  IP  PSAL            A0  D��fG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A0  D��fG�O�                