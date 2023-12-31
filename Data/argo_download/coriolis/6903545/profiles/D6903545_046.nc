CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  S   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:28Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
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
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200829092228  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               .A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��vO��O1   @��wW:�@R�����N4ϣ�1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 8.8 dbar]                                                      A��A!��A0  A@  AK33A[33Aq��A���A�  A�ffA�33A�33A�ffA�33A�33A�  A�  Aљ�A�  A�33A�ffA홚A���A���B��BffB
��B��B��BffBffBffB"ffB&  B*  B.  B2  B7��B=��BA��BE33BI33BL��BP��BT��BXffB\ffB`  Bd  Bh  Bk��Bo��Bs��Bw33B{33B33B���B�ffB�ffB�33B�33B�  B���B���B�ffB�33B�  B���B���B�ffB�33B���B���B�ffB�33B�  B���B���B�ffB���B���B�ffB�33B�33B�  B���B���B�ffB�33B�33B�  B���B˙�B�ffB�33B�  B���B���B���B晚BꙚBB�B���B���B���CL�CL�CL�CL�C	L�CffCffC� C��C��C�3C�3C��C��C�fCffC!  C#  C%�C'33C)33C+33C-� C/� C1ffC3ffC5ffC7� C9� C;� C=��C?��CA33CCL�CEL�CGL�CIffCK� CM�3COffCQffCS��CUffCW�CYL�C[� C]L�C_� Ca�3CcffCeL�Cg�CiffCk��Cm� CoL�Cq33Cs�CuffCw��Cy�3C{��C}� CffC��fC���C���C�� C��3C��fC�ٚC�ٚC���C���C�� C�� C��3C��fC��fC���C���C�� C��3C��fC�ٚC�ٚC���C�� C��3C��fC��fC��fC���C���C���C���C���C���C�� C�� C�� C�� C�� C��3C��3C��3C��3C��fC��3C��3C��3C��3C�� C�� C�� C�� C�� C���C���C���C���C��fC��fC��3C��3C�� C�� C���C�ٚC��fC��3C�� CĀ Cŀ Cƀ CǦfCȦfCɦfCʳ3C˳3C�� C�� C�ٚCϳ3C�� C���CҦfC�� C���CզfC�� C���CئfC�� C�ٚC۳3Cܙ�C�� C�ٚC�� C�fC��C�3C���C�� C�fC晚C�� C��fC���C�� C�3C왚C��C�3C��fC�ٚC���C�� C�3C��fC���C���C���C���C�� C�s3C�33D Y�Dy�D� D� D&fDs3D�fD	  D
33D� D��DfD@ D� D��D��D9�Ds3D��D  DFfD��D�3D  D33D� D �3D!��D#FfD$�fD%�fD'fD(33D)� D*��D+��D-33D.l�D/�fD0�fD2&fD3ffD4�fD5�fD7&fD8l�D9��D;fD<L�D=s3D>� D?��DA33DBl�DC� DE�DFL�DG�fDH�fDJfDKFfDL��DM��DN�fDP33DQ� DR�3DS��DU@ DV� DW��DX��DZ9�D[� D\�fD]�3D_FfD`y�Da��Db� Dd9�De�3Df��Dh�DiFfDj� Dk� Dm  DnFfDo�fDp�3Dr�Ds9�Dt` Du��Dw�Dx@ Dyy�Dz��D{�3D},�D~l�D��D�y�D�  D��3D�\�D��fD�� D�<�D�ٚD�y�D�fD��fD�Y�D���D���D�C3D�ٚD�� D�&fD�� D�Y�D��fD��3D�33D��3D�s3D�3D��fD�Y�D�  D��fD�@ D���D�|�D��D���D�` D�3D���D�9�D�� D�|�D��D���D�VfD��fD���D�<�D��3D��fD��D��fD�P D���D���D�@ D�� D�p D� D��3D�S3D�  D���D�9�D��3D�y�D�#3D���D�\�D���D�� D�C3D��3D��fD�  D��fD�Y�D�  D���D�6fD�� D�� D�  D�� D�` D�3D���D�@ D��D��3D��D��fD�S3D���D���D�L�D���D�|�D� D��3D�Y�D��3D���D�FfD��fD���D��D¶fD�\�D�fDĠ D�<�D��3Dƃ3D�  D�� D�` D�3DɖfD�9�D���Dˀ D�#3D̼�D�VfD�� DΜ�D�I�D��DЉ�D�&fD��3D�c3D�  DӠ D�C3D��3DՃ3D�fDֹ�D�` D��fDؓ3D�<�D��fDچfD�&fD��fD�i�D���Dݓ3D�6fD���D߀ D�#3D๚D�VfD�  D♚D�6fD�� D� D��D�� D�c3D�fD癚D�@ D��fD�|�D�fD�� D�ffD�  D왚D�33D�� D��D�)�D��fD�c3D�3D�3D�C3D��fD�D��D���D�c3D�  D���D�<�D���D�� D�  D���D�S3D�3D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A��A!��A0  A@  AK33A[33Aq��A���A�  A�ffA�33A�33A�ffA�33A�33A�  A�  Aљ�A�  A�33A�ffA홚A���A���B��BffB
��B��B��BffBffBffB"ffB&  B*  B.  B2  B7��B=��BA��BE33BI33BL��BP��BT��BXffB\ffB`  Bd  Bh  Bk��Bo��Bs��Bw33B{33B33B���B�ffB�ffB�33B�33B�  B���B���B�ffB�33B�  B���B���B�ffB�33B���B���B�ffB�33B�  B���B���B�ffB���B���B�ffB�33B�33B�  B���B���B�ffB�33B�33B�  B���B˙�B�ffB�33B�  B���B���B���B晚BꙚBB�B���B���B���CL�CL�CL�CL�C	L�CffCffC� C��C��C�3C�3C��C��C�fCffC!  C#  C%�C'33C)33C+33C-� C/� C1ffC3ffC5ffC7� C9� C;� C=��C?��CA33CCL�CEL�CGL�CIffCK� CM�3COffCQffCS��CUffCW�CYL�C[� C]L�C_� Ca�3CcffCeL�Cg�CiffCk��Cm� CoL�Cq33Cs�CuffCw��Cy�3C{��C}� CffC��fC���C���C�� C��3C��fC�ٚC�ٚC���C���C�� C�� C��3C��fC��fC���C���C�� C��3C��fC�ٚC�ٚC���C�� C��3C��fC��fC��fC���C���C���C���C���C���C�� C�� C�� C�� C�� C��3C��3C��3C��3C��fC��3C��3C��3C��3C�� C�� C�� C�� C�� C���C���C���C���C��fC��fC��3C��3C�� C�� C���C�ٚC��fC��3C�� CĀ Cŀ Cƀ CǦfCȦfCɦfCʳ3C˳3C�� C�� C�ٚCϳ3C�� C���CҦfC�� C���CզfC�� C���CئfC�� C�ٚC۳3Cܙ�C�� C�ٚC�� C�fC��C�3C���C�� C�fC晚C�� C��fC���C�� C�3C왚C��C�3C��fC�ٚC���C�� C�3C��fC���C���C���C���C�� C�s3C�33D Y�Dy�D� D� D&fDs3D�fD	  D
33D� D��DfD@ D� D��D��D9�Ds3D��D  DFfD��D�3D  D33D� D �3D!��D#FfD$�fD%�fD'fD(33D)� D*��D+��D-33D.l�D/�fD0�fD2&fD3ffD4�fD5�fD7&fD8l�D9��D;fD<L�D=s3D>� D?��DA33DBl�DC� DE�DFL�DG�fDH�fDJfDKFfDL��DM��DN�fDP33DQ� DR�3DS��DU@ DV� DW��DX��DZ9�D[� D\�fD]�3D_FfD`y�Da��Db� Dd9�De�3Df��Dh�DiFfDj� Dk� Dm  DnFfDo�fDp�3Dr�Ds9�Dt` Du��Dw�Dx@ Dyy�Dz��D{�3D},�D~l�D��D�y�D�  D��3D�\�D��fD�� D�<�D�ٚD�y�D�fD��fD�Y�D���D���D�C3D�ٚD�� D�&fD�� D�Y�D��fD��3D�33D��3D�s3D�3D��fD�Y�D�  D��fD�@ D���D�|�D��D���D�` D�3D���D�9�D�� D�|�D��D���D�VfD��fD���D�<�D��3D��fD��D��fD�P D���D���D�@ D�� D�p D� D��3D�S3D�  D���D�9�D��3D�y�D�#3D���D�\�D���D�� D�C3D��3D��fD�  D��fD�Y�D�  D���D�6fD�� D�� D�  D�� D�` D�3D���D�@ D��D��3D��D��fD�S3D���D���D�L�D���D�|�D� D��3D�Y�D��3D���D�FfD��fD���D��D¶fD�\�D�fDĠ D�<�D��3Dƃ3D�  D�� D�` D�3DɖfD�9�D���Dˀ D�#3D̼�D�VfD�� DΜ�D�I�D��DЉ�D�&fD��3D�c3D�  DӠ D�C3D��3DՃ3D�fDֹ�D�` D��fDؓ3D�<�D��fDچfD�&fD��fD�i�D���Dݓ3D�6fD���D߀ D�#3D๚D�VfD�  D♚D�6fD�� D� D��D�� D�c3D�fD癚D�@ D��fD�|�D�fD�� D�ffD�  D왚D�33D�� D��D�)�D��fD�c3D�3D�3D�C3D��fD�D��D���D�c3D�  D���D�<�D���D�� D�  D���D�S3D�3D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�+@�J@�v�@c@t��@fv�@ko@>E�@E��@?�@=?}@=��@�?ٙ�?��?��;?��?���@+@dZ@�^?�ff?�33?��#?�9X?щ7?�ƨ?��@
�!@�@��@��@p�@9�^@3"�@;ƨ@;S�@@  @?��@:�@5��@2J@1�^@0�9@4��@6�y@6E�@6V@7+@7�@7;d@7�@5V@2�@.{@,�@*�@&�@&{@$�j@$�@#S�@"^5@"�H@"=q@%��@$(�@#�@ Q�@|�@�@��@ff@z�@�\@7L@��@�P@5?@�@�F@@33@�@~�@�#@�^@��@7L@bN@�@  @�u@1'@�@�@|�@��@J?��?���?�I�?�O�?�C�?�"�?���?��?��??�D?��?�"�?�?�1'?���?���?�v�?��H?�r�?��
?щ7?�&�?�Ĝ?Ͳ-?���?�Z?�S�?�S�?��7?�C�?�?}?�A�?���?���?�ƨ?�?��?��F?��?�b?�ff?��!?�%?��?��?��h?�O�?�|�?���?�C�?��u?�K�?��y?���?���?��?�bN?|�?y�#?st�?m��?l��?h��?h1'?f�y?fff?d�?`A�?^��?^5??_|�?c��?f$�?h��?f��?e�T?b�\?^5??Y��?U?Tz�?S��?Q&�?O�;?O��?P �?P �?M�h?I�^?G+?F$�?DZ??�w?5�?4�j?4z�?;dZ?>v�?@�?BJ?A%?@A�?>v�?=�-?;�m?;�m?;�m?;��?;"�?9��?8Q�?7��?7�P?7K�?6ȴ?7
=?6ȴ?6ȴ?6ȴ?6E�?3�F?0��?0��?0�`?1&�?0��?/\)?,�D?+�?-��?.�?.��?,�D?*��?+ƨ?*��?(�9?&�y?&$�?!�7?/?dZ?�j?��?I�?ƨ?
=q?1'?�y?�/?M�>���>��H>�Q�>�>�->�V>�ff>�;d>ڟ�>�hs>ě�>�ȴ>���>�G�>��>���>���>��^>���>��>���>���>��7>vȴ>r�!>n��>ix�>ix�>j~�>w��>��\>�ƨ>�V>�V>�V>��^>��>��>}�>y�#>s�F>j~�>gl�>dZ>bM�>bM�>V>=p�>7K�>B�\>D��>7K�>,1>��>�+>I�>$�=��>%=��m=�h=�"�=�"�=��=�^5=�v�=ě�=��
=�hs=y�#=@�=8Q�=P�`=D��=�P<���<��
;o�D����`B;��
;�o    �D���u���m�h�q���H�9�m�h���
��Q콼j����"ѽ��������������������ͽ����;d��G������o�I��z��P�����R�"��(�þ333�0 ž6E��=p��@��B�\�E�˾I�^�Kƨ�N��W
=�Y��Xb�_;d�cS��cS��fff�gl��hr��l�D�o���r�!�t�j�w�پ{�m�|푾}󶾁%������˾�1'��7L���^��ƨ��O߾���V�����\)��bN��bN��hs��n����Ͼ�����+���P�����"Ѿ����;d��;d��A��������徤Z��`B��ff��l����þ�~��������h�������!���j��ȴ������^5��^5���m��p���vɾ����%���7��o����ȴ9������ƨ������;��hs��녾��Ͼ����Ձ�׍P����ٙ������"Ѿۥ��(���;d�������/��`B��ff��r��������1��V����&��33��?}��ȴ��KǾ��#���#��dZ��dZ��j��vɿ ��J�Mӿ��o����Z����˿��+�l��1'��9��9��ÿ�ÿ	�^�
=q����I���D�O߿�h���V���\)���� ſ&����-�33��!��!��33��F�9X�9X��j�?}��ȴ�ȴ�Kǿb�Q����X��#�^5��H�dZ���(����푿/��5?��R�|�   � Ĝ�!%�!G��!�7�#o�#S��#�
�$Z�$�/�%��%�˿%�T�&$ݿ&��'l��'(�9�(�ÿ)7L�)xտ)�^�)��*=q�*���+�+�+�+C��+��,1�,I��,�D�-V�-V�-O߿-�h�-��.{�.V�.��/\)�/\)�/\)�/�;�/�;1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�+@�J@�v�@c@t��@fv�@ko@>E�@E��@?�G�O�G�O�@�?ٙ�?��?��;?��?���@+@dZ@�^?�ff?�33?��#?�9X?щ7?�ƨ?��@
�!@�@��@��@p�@9�^@3"�@;ƨ@;S�@@  @?��@:�@5��@2J@1�^@0�9@4��@6�y@6E�@6V@7+@7�@7;d@7�@5V@2�@.{@,�@*�@&�@&{@$�j@$�@#S�@"^5@"�H@"=q@%��@$(�@#�@ Q�@|�@�@��@ff@z�@�\@7L@��@�P@5?@�@�F@@33@�@~�@�#@�^@��@7L@bN@�@  @�u@1'@�@�@|�@��@J?��?���?�I�?�O�?�C�?�"�?���?��?��??�D?��?�"�?�?�1'?���?���?�v�?��H?�r�?��
?щ7?�&�?�Ĝ?Ͳ-?���?�Z?�S�?�S�?��7?�C�?�?}?�A�?���?���?�ƨ?�?��?��F?��?�b?�ff?��!?�%?��?��?��h?�O�?�|�?���?�C�?��u?�K�?��y?���?���?��?�bN?|�?y�#?st�?m��?l��?h��?h1'?f�y?fff?d�?`A�?^��?^5??_|�?c��?f$�?h��?f��?e�T?b�\?^5??Y��?U?Tz�?S��?Q&�?O�;?O��?P �?P �?M�h?I�^?G+?F$�?DZ??�w?5�?4�j?4z�?;dZ?>v�?@�?BJ?A%?@A�?>v�?=�-?;�m?;�m?;�m?;��?;"�?9��?8Q�?7��?7�P?7K�?6ȴ?7
=?6ȴ?6ȴ?6ȴ?6E�?3�F?0��?0��?0�`?1&�?0��?/\)?,�D?+�?-��?.�?.��?,�D?*��?+ƨ?*��?(�9?&�y?&$�?!�7?/?dZ?�j?��?I�?ƨ?
=q?1'?�y?�/?M�>���>��H>�Q�>�>�->�V>�ff>�;d>ڟ�>�hs>ě�>�ȴ>���>�G�>��>���>���>��^>���>��>���>���>��7>vȴ>r�!>n��>ix�>ix�>j~�>w��>��\>�ƨ>�V>�V>�V>��^>��>��>}�>y�#>s�F>j~�>gl�>dZ>bM�>bM�>V>=p�>7K�>B�\>D��>7K�>,1>��>�+>I�>$�=��>%=��m=�h=�"�=�"�=��=�^5=�v�=ě�=��
=�hs=y�#=@�=8Q�=P�`=D��=�P<���<��
;o�D����`B;��
;�o    �D���u���m�h�q���H�9�m�h���
��Q콼j����"ѽ��������������������ͽ����;d��G������o�I��z��P�����R�"��(�þ333�0 ž6E��=p��@��B�\�E�˾I�^�Kƨ�N��W
=�Y��Xb�_;d�cS��cS��fff�gl��hr��l�D�o���r�!�t�j�w�پ{�m�|푾}󶾁%������˾�1'��7L���^��ƨ��O߾���V�����\)��bN��bN��hs��n����Ͼ�����+���P�����"Ѿ����;d��;d��A��������徤Z��`B��ff��l����þ�~��������h�������!���j��ȴ������^5��^5���m��p���vɾ����%���7��o����ȴ9������ƨ������;��hs��녾��Ͼ����Ձ�׍P����ٙ������"Ѿۥ��(���;d�������/��`B��ff��r��������1��V����&��33��?}��ȴ��KǾ��#���#��dZ��dZ��j��vɿ ��J�Mӿ��o����Z����˿��+�l��1'��9��9��ÿ�ÿ	�^�
=q����I���D�O߿�h���V���\)���� ſ&����-�33��!��!��33��F�9X�9X��j�?}��ȴ�ȴ�Kǿb�Q����X��#�^5��H�dZ���(����푿/��5?��R�|�   � Ĝ�!%�!G��!�7�#o�#S��#�
�$Z�$�/�%��%�˿%�T�&$ݿ&��'l��'(�9�(�ÿ)7L�)xտ)�^�)��*=q�*���+�+�+�+C��+��,1�,I��,�D�-V�-V�-O߿-�h�-��.{�.V�.��/\)�/\)�/\)�/�;�/�;1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��BC�BE�B�7B	B	�B	u�B	q�B	��B	�B	ǮB	�{B
�B	�B
\B
5?B
dZB
ZB
�B
�jB
�B
ŢB
ĜB
�-B
��B
�B
�;B
��B49BB�B;dBN�Bo�B�B�!BŢBÖB�B�5B�/B�B��B��B��B�B�fB�`B�fB�fB�sB�yB�B�yB�`B�HB�ZB�sB�)B�B�)B�5B�B�)B�;B�;B�HB�ZB�ZB�5B�;B�;B�5B�ZB�NB�/B�B�B�B�B�
B�B�B�
B�B�B�
B�B�B�B�B�B�#B�/B�)B�/B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBǮBŢBŢBȴBƨBŢBĜBĜBÖBB�}B�wB�qB�qB�jB�jB�dB�dB�^B�qB�9B�?B�3B�9B�?B�9B�?B�9B�RB�XB�RB�LB�LB�LB�FB�FB�FB�?B�9B�3B�3B�3B�-B�3B�-B�3B�-B�3B�3B�3B�9B�3B�FB�LB�XB�^B�XB�^B�RB�RB�LB�LB�LB�RB�FB�LB�RB�RB�XB�RB�LB�LB�LB�FB�9B�3B�?B�LB�^B�^B�jB�dB�jB�jB�jB�jB�jB�qB�qB�qB�qB�qB�qB�wB�wB�wB�wB�}B�wB�}B�}B�wB�wB�wB�}B�}B�}B�}B�wB�}B�}B��B��B��B��B��B��B�}B��B�}B�wB�jB�dB�^B�^B�RB�RB�RB�RB�LB�LB�LB�LB�FB�FB�?B�?B�9B�9B�-B�-B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B��BC�BE�B�7B	B	�B	u�B	q�B	��B	�G�O�G�O�B
�B	�B
\B
5?B
dZB
ZB
�B
�jB
�B
ŢB
ĜB
�-B
��B
�B
�;B
��B49BB�B;dBN�Bo�B�B�!BŢBÖB�B�5B�/B�B��B��B��B�B�fB�`B�fB�fB�sB�yB�B�yB�`B�HB�ZB�sB�)B�B�)B�5B�B�)B�;B�;B�HB�ZB�ZB�5B�;B�;B�5B�ZB�NB�/B�B�B�B�B�
B�B�B�
B�B�B�
B�B�B�B�B�B�#B�/B�)B�/B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBǮBŢBŢBȴBƨBŢBĜBĜBÖBB�}B�wB�qB�qB�jB�jB�dB�dB�^B�qB�9B�?B�3B�9B�?B�9B�?B�9B�RB�XB�RB�LB�LB�LB�FB�FB�FB�?B�9B�3B�3B�3B�-B�3B�-B�3B�-B�3B�3B�3B�9B�3B�FB�LB�XB�^B�XB�^B�RB�RB�LB�LB�LB�RB�FB�LB�RB�RB�XB�RB�LB�LB�LB�FB�9B�3B�?B�LB�^B�^B�jB�dB�jB�jB�jB�jB�jB�qB�qB�qB�qB�qB�qB�wB�wB�wB�wB�}B�wB�}B�}B�wB�wB�wB�}B�}B�}B�}B�wB�}B�}B��B��B��B��B��B��B�}B��B�}B�wB�jB�dB�^B�^B�RB�RB�RB�RB�LB�LB�LB�LB�FB�FB�?B�?B�9B�9B�-B�-B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
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
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092228                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092312  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200829092312  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A��D�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172539  IP  PSAL            A��D�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A��D�� G�O�                