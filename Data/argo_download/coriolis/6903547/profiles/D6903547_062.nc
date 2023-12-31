CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  W   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:30:27Z creation; 2022-03-27T13:14:32Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        	\  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  D4   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	\  F�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	\  R@   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  [�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  gP   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  p�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  s   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  |`   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  �p   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	\  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �    	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �$   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �T   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �T   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �T   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �TArgo profile    3.1 1.2 19500101000000  20200829093027  20220327131432  6903547 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               >A   IF                                  2C  D   ARVOR                           AI2600-18EU004                  5900A04                         844 @���8�91   @����ƻ`@R���\�@	��ڲ�1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 3.8 dbar]                                                      @�  @�  @�  @�33A   A  AffA1��A<��AL��A\��Ap  A~ffA���A�  A�33A���A���A�  A�  A�  A���A�33A�  AᙚA�  A�ffA���A�33B33BffB  B��B��BffBffB  B#��B)��B-33B0��B4ffB8  B;33B>��BD  BI33BL��BPffBT  BW��B[33B^��Bb��BfffBj��BnffBrffBvffBzffB~ffB�ffB�ffB�ffB���B���B�  B�33B�ffB�ffB���B���B�  B�33B�ffB�ffB���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B�ffB�33B�33B���Bș�B�33B�  B�ffB�  Bۙ�Bߙ�B�ffB�ffB�33B���B�B�33B���B�33CffC33CffC�3C	�3C�3CffC  CffC�3CffC  C�C33C33C33C!33C#�C%�C'�C)  C+  C-  C/ffC1�3C3� C5ffC7L�C933C;�C=�C?  CA  CC  CE  CG  CI  CK  CMffCO�3CQ��CSffCU33CW�CYffC[�3C]� C_ffCa33Cc� Ce�3Cg� CiL�Ck� Cm�3CoffCq�CsL�Cu� Cw33Cy� C{�3C}� CL�C���C��3C�ٚC�� C��3C��fC���C���C��3C��fC�ٚC�ٚC���C�� C��3C��3C��fC���C���C��3C��fC�ٚC�ٚC���C���C�� C�� C��3C��fC���C���C���C��3C��fC��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C�� C��3C��fC��fC��fC��fC��fC��fC��fC��3C��3C��3C�� C�� C�� C���C���C���C�ٚC�ٚC�� Cĳ3CŦfCƦfC���C�ٚCɦfC�� C���C̳3C͙�CΌ�Cϳ3C�ٚC�� CҦfCә�CԌ�C�� C��fC���Cس3Cٌ�CڦfC�� Cܙ�Cݳ3C�� CߦfC�3C���C�3C㙚C�� C��fC�ٚC���C�� C�fCꙚC�3C�ٚC���C�3C��C�fC�� C�C�fC��3C���C��3C���C��3C��fC�Y�C�  D FfDs3D�fD�3DL�D�fD� D	  D
@ D�fD��D�DFfD� D��D�3D33Ds3D�3D�3D33Ds3D�3D�3D,�Dl�D ��D!�3D#9�D$l�D%��D'fD(@ D)� D*� D+��D-9�D.y�D/� D0�3D2,�D3� D4� D6  D7FfD8s3D9�fD:�3D<FfD=� D>��D?�3DA33DBs3DC��DEfDF@ DGy�DH��DI��DK9�DLy�DM��DN��DP@ DQ��DR� DS�3DU,�DVl�DW��DX�3DZ9�D[y�D\�fD]��D_9�D`��Da� Db��Dd9�Dey�Df��DhfDiS3Dj� Dk��Dl�3Dn33Doy�Dp�3Dq��DsFfDt�fDu��Dv�3Dx9�Dy�fDz��D{��D}&fD~y�D��D��3D�  D���D�\�D�  D�� D�C3D��3D��fD�&fD���D�c3D��fD���D�FfD��3D�|�D��D���D�\�D���D��3D�9�D��3D�|�D��D��fD�` D���D���D�9�D��fD�y�D��D�� D�VfD�  D���D�@ D��fD�|�D�&fD��fD�ffD�fD��3D�@ D�� D��3D�&fD��fD�i�D�	�D���D�0 D��3D�vfD�fD��fD�VfD��fD��fD�6fD���D�� D�#3D��fD�ffD�  D���D�<�D�ٚD�y�D�fD��3D�P D���D���D�<�D��fD�|�D�  D��fD�` D���D�� D�9�D��3D�y�D�#3D���D�\�D���D���D�<�D�� D��3D�fD���D�c3D���D�� D�@ D��fD�s3D� D¼�D�i�D�fDģ3D�C3D��Dƃ3D��D�� D�i�D�fDɣ3D�@ D���D�|�D��D��3D�\�D��fDΜ�D�6fD��fD�y�D�  D��fD�` D�  DӠ D�@ D��3DՃ3D�&fD��fD�Y�D�� Dؠ D�C3D���D�vfD�  Dۼ�D�P D��3DݖfD�@ D���D�y�D�&fD��3D�` D�3D�fD�<�D��D�fD�#3D�� D�` D�  D� D�<�D��fD�|�D�  D�� D�` D�  D�3D�<�D��3D�|�D�&fD��fD�i�D�3D��D�@ D��3D�vfD��D���D�Y�D�3D��fD�C3D���D�� D�#3D��fD�i�D�  D�Vf11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�  @�  @�  @�33A   A  AffA1��A<��AL��A\��Ap  A~ffA���A�  A�33A���A���A�  A�  A�  A���A�33A�  AᙚA�  A�ffA���A�33B33BffB  B��B��BffBffB  B#��B)��B-33B0��B4ffB8  B;33B>��BD  BI33BL��BPffBT  BW��B[33B^��Bb��BfffBj��BnffBrffBvffBzffB~ffB�ffB�ffB�ffB���B���B�  B�33B�ffB�ffB���B���B�  B�33B�ffB�ffB���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B���B���B���B�ffB�33B�33B���Bș�B�33B�  B�ffB�  Bۙ�Bߙ�B�ffB�ffB�33B���B�B�33B���B�33CffC33CffC�3C	�3C�3CffC  CffC�3CffC  C�C33C33C33C!33C#�C%�C'�C)  C+  C-  C/ffC1�3C3� C5ffC7L�C933C;�C=�C?  CA  CC  CE  CG  CI  CK  CMffCO�3CQ��CSffCU33CW�CYffC[�3C]� C_ffCa33Cc� Ce�3Cg� CiL�Ck� Cm�3CoffCq�CsL�Cu� Cw33Cy� C{�3C}� CL�C���C��3C�ٚC�� C��3C��fC���C���C��3C��fC�ٚC�ٚC���C�� C��3C��3C��fC���C���C��3C��fC�ٚC�ٚC���C���C�� C�� C��3C��fC���C���C���C��3C��fC��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C�� C��3C��fC��fC��fC��fC��fC��fC��fC��3C��3C��3C�� C�� C�� C���C���C���C�ٚC�ٚC�� Cĳ3CŦfCƦfC���C�ٚCɦfC�� C���C̳3C͙�CΌ�Cϳ3C�ٚC�� CҦfCә�CԌ�C�� C��fC���Cس3Cٌ�CڦfC�� Cܙ�Cݳ3C�� CߦfC�3C���C�3C㙚C�� C��fC�ٚC���C�� C�fCꙚC�3C�ٚC���C�3C��C�fC�� C�C�fC��3C���C��3C���C��3C��fC�Y�C�  D FfDs3D�fD�3DL�D�fD� D	  D
@ D�fD��D�DFfD� D��D�3D33Ds3D�3D�3D33Ds3D�3D�3D,�Dl�D ��D!�3D#9�D$l�D%��D'fD(@ D)� D*� D+��D-9�D.y�D/� D0�3D2,�D3� D4� D6  D7FfD8s3D9�fD:�3D<FfD=� D>��D?�3DA33DBs3DC��DEfDF@ DGy�DH��DI��DK9�DLy�DM��DN��DP@ DQ��DR� DS�3DU,�DVl�DW��DX�3DZ9�D[y�D\�fD]��D_9�D`��Da� Db��Dd9�Dey�Df��DhfDiS3Dj� Dk��Dl�3Dn33Doy�Dp�3Dq��DsFfDt�fDu��Dv�3Dx9�Dy�fDz��D{��D}&fD~y�D��D��3D�  D���D�\�D�  D�� D�C3D��3D��fD�&fD���D�c3D��fD���D�FfD��3D�|�D��D���D�\�D���D��3D�9�D��3D�|�D��D��fD�` D���D���D�9�D��fD�y�D��D�� D�VfD�  D���D�@ D��fD�|�D�&fD��fD�ffD�fD��3D�@ D�� D��3D�&fD��fD�i�D�	�D���D�0 D��3D�vfD�fD��fD�VfD��fD��fD�6fD���D�� D�#3D��fD�ffD�  D���D�<�D�ٚD�y�D�fD��3D�P D���D���D�<�D��fD�|�D�  D��fD�` D���D�� D�9�D��3D�y�D�#3D���D�\�D���D���D�<�D�� D��3D�fD���D�c3D���D�� D�@ D��fD�s3D� D¼�D�i�D�fDģ3D�C3D��Dƃ3D��D�� D�i�D�fDɣ3D�@ D���D�|�D��D��3D�\�D��fDΜ�D�6fD��fD�y�D�  D��fD�` D�  DӠ D�@ D��3DՃ3D�&fD��fD�Y�D�� Dؠ D�C3D���D�vfD�  Dۼ�D�P D��3DݖfD�@ D���D�y�D�&fD��3D�` D�3D�fD�<�D��D�fD�#3D�� D�` D�  D� D�<�D��fD�|�D�  D�� D�` D�  D�3D�<�D��3D�|�D�&fD��fD�i�D�3D��D�@ D��3D�vfD��D���D�Y�D�3D��fD�C3D���D�� D�#3D��fD�i�D�  D�Vf11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@;��@;��@;�F@;��@;�
@;�m@;��@;ƨ@;��@;�@;��@;��@;�@;t�@;t�@;�@<(�@<��@<�D@<�D@<��@<�D@<�@<��@<��@<��@<��@<��@<��@<�@<�@<�j@<��@<�@<�D@<��@<(�@<9X@<�j@<j@<I�@<(�@<9X@<Z@<�@<�j@<z�@<z�@<�D@<�D@<�D@<�D@<�D@<�@<�@<�j@<�@<�@<�@<�@<�@<�D@<�D@<z�@<Z@<9X@<9X@<9X@<I�@<Z@<z�@<�D@<(�@<9X@<�D@<��@<�D@<��@<�D@<��@<�D@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<�@<��@<��@<�@<��@<�D@<�@<�D@<�D@<�D@<��@<�D@<��@<��@<��@<�j@<�@<�@<��@<��@<�@<�@<�@<�@<�/@<��@<�j@<�@<�@<�@<�@<�j@<�j@<�j@<�D@<z�@<�D@<�D@<z�@<��@<�@<��@<�@<�@<��@<�D@<��@<��@<�@<�j@<��@<z�@<�D@<z�@<z�@<z�@<�D@<�D@<�D@<z�@<��@<��@<�@<�j@<�j@<�@<��@<�D@<�D@<j@;��@;��@;C�@;��@;�@:�@:-@3t�@2n�@2J@)��@(r�@$��@"n�@!G�@ b@@O�@�@�D@j@I�@�@�
@S�@�H@�#@�u@�@(�@;d@z�@|�@@�@Z@�@�/@1@dZ@��@�\@-@hs@7L@�@ ��@ ��?��;?��h?�+?��?��y?Ͼw?ȓu?�Ĝ?�Q�?�Z?� �?���?�{?�1?���?�x�?�Q�?�
=?��+?���?��?�-?���?��T?�G�?}p�?{dZ?w�P?rn�?l��?_;d?W�P?B�\?+?(1'?$�/?��?��?/?v�?��?   ?�R?�m?��?K�?��?�+?�?�!?-?�`?V?	�^?r�?�9?�T?��>�j>�Ĝ>ܬ>ٙ�>�n�>ȴ9>ě�>Õ�>��>�`B>��>�\)>���>�O�>�O�>���>���>�7L>�7L>�1'>�%>q��>\(�>V>D��>8Q�>,1>#�
>bN>C�>%=�/=�1=�\)=ix�=D��=#�
<��<�/<�9X<T��;D����`B�#�
�49X�T����C����49X�T���ixսm�h�}�y�#�e`B��O߽��P���㽛�㽰 Ž��ͽ�/��G���;d���`�����G���xս���ٽ��#�o�+�
=q�C��O߾\)�bN�hs�z��+��P�����w�!���$�/�%�T�')��,1�.{�0 ž49X�8Q�:^5�;dZ�=p��C���H�9�O�;�R�T���W
=�\(��^5?�_;d�`A��cS��e`B�fff�ixվk��k��l�D�n���s�F�s�F�s�F�u�u�vȴ�x���{�m��%��o������˾���+��+��+��������1'��1'��1'���^������ƨ��ƨ��ƨ��ƨ��I����;��;��n���t����Ͼ��Ͼ��+�����;d��Ĝ������S���ff����þ�~�������D��{�����������F��?}���پ�Q쾷KǾ�ȴ���#��푾�푾��۾��7�\����š˾�1'��ƨ���`��녾�n���t������Ձ�Ձ�և+��
=��b�ٙ��ڟ���/��5?��A���A���Ĝ�����`B���T��`B���/��l���xվ�~���1��{��!��F��!��9X����ȴ��Q��Q��Q���#���H���m��j���۾�|� A��J����
��������$ݿ+�	xտ
~��
=q�
~��
������D�O߿V�{�{�V����`����-�n���!��!��!�t���F�9X�?}�E��ȴ�
=��������#���dZ����m�/�5?�vɿ�ۿ�ۿ;d�|��w�   � A�� ��!%�!���#o�#���#�
�$��$��$Z�$���$�/�%��%��%��%��%`B�%�T�&ff�&�y�&�y�'+�'l��'(1'�(r��(r��(r��(�ÿ)��*���+�*���*���*���+C��+�+�+C��+C��+��,�D�-��/���/�;�0bN�1���1���1녿2-�2-�333�3t��33311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @;��@;��@;�F@;��@;�
@;�m@;��@;ƨ@;��@;�@;��@;��@;�@;t�@;t�@;�@<(�@<��@<�D@<�D@<��@<�D@<�@<��@<��@<��@<��@<��@<��@<�@<�@<�j@<��@<�@<�D@<��@<(�@<9X@<�j@<j@<I�@<(�@<9X@<Z@<�@<�j@<z�@<z�@<�D@<�D@<�D@<�D@<�D@<�@<�@<�j@<�@<�@<�@<�@<�@<�D@<�D@<z�@<Z@<9X@<9X@<9X@<I�@<Z@<z�@<�D@<(�@<9X@<�D@<��@<�D@<��@<�D@<��@<�D@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<��@<�@<��@<��@<�@<��@<�D@<�@<�D@<�D@<�D@<��@<�D@<��@<��@<��@<�j@<�@<�@<��@<��@<�@<�@<�@<�@<�/@<��@<�j@<�@<�@<�@<�@<�j@<�j@<�j@<�D@<z�@<�D@<�D@<z�@<��@<�@<��@<�@<�@<��@<�D@<��@<��@<�@<�j@<��@<z�@<�D@<z�@<z�@<z�@<�D@<�D@<�D@<z�@<��@<��@<�@<�j@<�j@<�@<��@<�D@<�D@<j@;��@;��@;C�@;��@;�@:�@:-@3t�@2n�@2J@)��@(r�@$��@"n�@!G�@ b@@O�@�@�D@j@I�@�@�
@S�@�H@�#@�u@�@(�@;d@z�@|�@@�@Z@�@�/@1@dZ@��@�\@-@hs@7L@�@ ��@ ��?��;?��h?�+?��?��y?Ͼw?ȓu?�Ĝ?�Q�?�Z?� �?���?�{?�1?���?�x�?�Q�?�
=?��+?���?��?�-?���?��T?�G�?}p�?{dZ?w�P?rn�?l��?_;d?W�P?B�\?+?(1'?$�/?��?��?/?v�?��?   ?�R?�m?��?K�?��?�+?�?�!?-?�`?V?	�^?r�?�9?�T?��>�j>�Ĝ>ܬ>ٙ�>�n�>ȴ9>ě�>Õ�>��>�`B>��>�\)>���>�O�>�O�>���>���>�7L>�7L>�1'>�%>q��>\(�>V>D��>8Q�>,1>#�
>bN>C�>%=�/=�1=�\)=ix�=D��=#�
<��<�/<�9X<T��;D����`B�#�
�49X�T����C����49X�T���ixսm�h�}�y�#�e`B��O߽��P���㽛�㽰 Ž��ͽ�/��G���;d���`�����G���xս���ٽ��#�o�+�
=q�C��O߾\)�bN�hs�z��+��P�����w�!���$�/�%�T�')��,1�.{�0 ž49X�8Q�:^5�;dZ�=p��C���H�9�O�;�R�T���W
=�\(��^5?�_;d�`A��cS��e`B�fff�ixվk��k��l�D�n���s�F�s�F�s�F�u�u�vȴ�x���{�m��%��o������˾���+��+��+��������1'��1'��1'���^������ƨ��ƨ��ƨ��ƨ��I����;��;��n���t����Ͼ��Ͼ��+�����;d��Ĝ������S���ff����þ�~�������D��{�����������F��?}���پ�Q쾷KǾ�ȴ���#��푾�푾��۾��7�\����š˾�1'��ƨ���`��녾�n���t������Ձ�Ձ�և+��
=��b�ٙ��ڟ���/��5?��A���A���Ĝ�����`B���T��`B���/��l���xվ�~���1��{��!��F��!��9X����ȴ��Q��Q��Q���#���H���m��j���۾�|� A��J����
��������$ݿ+�	xտ
~��
=q�
~��
������D�O߿V�{�{�V����`����-�n���!��!��!�t���F�9X�?}�E��ȴ�
=��������#���dZ����m�/�5?�vɿ�ۿ�ۿ;d�|��w�   � A�� ��!%�!���#o�#���#�
�$��$��$Z�$���$�/�%��%��%��%��%`B�%�T�&ff�&�y�&�y�'+�'l��'(1'�(r��(r��(r��(�ÿ)��*���+�*���*���*���+C��+�+�+C��+C��+��,�D�-��/���/�;�0bN�1���1���1녿2-�2-�333�3t��33311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�sB�sB�`B�TB�5B�;B�#B�/B�5B�BB�;B�BB�;B�BB�BB�;B�BB�;B�5B�/B�/B�/B�B�
BƨBŢB��B��B�B�jB�jB�^B�dB�^B�XB�XB�XB�^B�RB�RB�RB�XB�?B�FB�RB�FB�?B�?B�?B�3B�B�B�B�B�B�-B�B�B�B�B�B�B�B��B�B��B�B��B�B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�B�!B�B�!B�!B�!B�B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�'B�!B�!B�!B�'B�!B�!B�!B�'B�'B�!B�'B�'B�'B�'B�'B�'B�'B�'11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�RB�LB�EB�LB�LB�EB�LB�LB�LB�EB�EB�EB�LB�LB�LB�LB�?B�EB�EB�EB�EB�LB�EB�EB�EB�EB�EB�LB�EB�EB�EB�EB�EB�EB�LB�-B�LB�EB�EB�^B�EB�LB�LB�LB�LB�EB�EB�EB�EB�LB�LB�EB�LB�EB�EB�EB�EB�EB�EB�EB�LB�EB�EB�?B�LB�EB�LB�LB�EB�RB�LB�EB�EB�?B�EB�EB�EB�EB�LB�EB�EB�EB�LB�EB�EB�EB�EB�EB�EB�LB�EB�LB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�?B�EB�EB�?B�EB�EB�EB�?B�EB�LB�EB�EB�EB�EB�EB�EB�?B�EB�EB�EB�EB�EB�EB�EB�?B�EB�EB�EB�EB�EB�EB�?B�EB�EB�EB�EB�EB�EB�EB�?B�?B�EB�?B�?B�?B�?B�?B�EB�?B�EB�EB�EB�EB�?B�3B�9B�?B�?B�?B�3B�3B�'B�B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BܾBڲBՓB֙BҁBԍBՓBנB֙BנB֙BנBנB֙BנB֙BՓBԍBԍBԍB�{B�hB�B� B��B��B�rB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�rB�rB�lB�rB�yB��B�yB�lB�fB�lB�yB�lB�fB�NB�lB�NB�fB�ZB�lB�`B�`B�ZB�`B�TB�TB�ZB�ZB�GB�NB�GB�GB�AB�NB�GB�;B�AB�GB�;B�GB�AB�;B�AB�AB�;B�;B�5B�;B�AB�;B�5B�)B�5B�5B�5B�;B�5B�5B�AB�5B�;B�5B�;B�;B�;B�;B�5B�;B�;B�;B�;B�;B�AB�;B�;B�;B�AB�;B�;B�;B�;B�;B�AB�;B�NB�AB�GB�GB�GB�NB�AB�GB�GB�GB�NB�TB�NB�TB�TB�TB�TB�TB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�`B�`B�`B�`B�`B�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�lB�fB�lB�lB�lB�lB�lB�lB�lB�lB�lB�lB�lB�lB�rB�lB�lB�lB�rB�lB�lB�rB�lB�lB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�yB�yB�yB�yB�yB�rB�rB�rB�yB�rB�rB�rB�rB�rB�rB�lB�rB�rB�rB�rB�rB�lB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�lB�rB�lB�rB�rB�rB�lB�rB�rB�rB�lB�rB�rB�rB�rB�lB�rB�rB�rB�rB�rB�rB�rB�lB�lB�lB�rB�rB�lB�rB�lB�rB�rB�lB�rB�lB�lB�rB�rB�lB�lB�lB�rB�rB�rB�lB�rB�rB�lB�rB�lB�lB�lB�lB�lB�lB�lB�rB�lB�rB�rB�rB�rB�lB�rB�lB�lB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�rB�yB�rB�rB�rB�rB�rB�rB�rB�rB�yB�rB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�B�yB�yB�B�yB�B�B�B�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B��B�B�B�B��B��B�B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 0.99978, vertically averaged dS= -0.0084303                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202203271314322022032713143220220327131432  IF  ARFMCODA035h                                                                20200829093027                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829093129  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200829093129  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200918073012  IP  PSAL            @�  D�VfG�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20210128114650  IP  PSAL            @�  D�VfG�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20220327131432  IP  PSAL            @�  D�VfG�O�                