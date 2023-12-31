CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
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
resolution        =���   axis      Z        �  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  BD   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D<   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  L   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g|   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  it   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  qH   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  y   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  {   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �    HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �`   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �p   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �t   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170405  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               =A   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @��I��J1   @��I��J@S`c<�@l���Pj�8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A0  A>ffAQ��A`  Ap  A�  A�33A�  A���A�  A���A�ffA�33A���A�33A�33A���A�ffA�33A�  A�33A�ffB  B��B��BffB  B  B  B��B#��B'33B+��B/33B333B733B;��B?33BC��BG��BL  BO��BTffBV��B[33B_33Bc33Bg��Bl  Bp��Bs��Bv��B{��B�  B���B�  B�ffB�  B���B�ffB���B�ffB�33B���B�ffB�  B���B�ffB�  B���B���B�ffB���B���B�ffB�33B�  B���B���B�ffB�33B���B���B�ffB�33B�  B���BÙ�B�ffB�  B�ffB�  Bҙ�B�33B���B�  B♚B�  B�ffB�  B�ffB�  B�ffB�  C33C� C�3C� C	L�C�CffC��CffCL�C�CffC��C��C� CL�C!33C#�C%ffC'�3C)�3C+��C-� C/� C1L�C3L�C533C733C933C;33C=33C?�CA�CCL�CE��CG� CIL�CK��CM� COffCQL�CSL�CUL�CWL�CYL�C[L�C]� C_� Ca�3CcffCe� Cg��CiffCk�CmL�Co� Cq33Cs  CuL�Cw� CyffC{33C}� C�3C���C��3C��fC���C��3C�ٚC���C��3C���C���C��3C���C��3C���C���C��3C���C��3C���C�� C�ٚC���C��3C��fC���C��3C��fC�ٚC���C�� C��3C��3C��fC���C���C�� C�� C��3C��3C��fC��fC��fC��fC��fC�ٚC��fC��fC��3C��3C�� C�� C�� C�� C���C���C���C���C���C��fC��fC��3C��3C�� C���C�ٚC��fC��fC�� CĀ CŌ�Cƙ�CǙ�CȦfCɳ3C�� C�ٚC�ٚC��fCγ3Cπ CЌ�Cљ�CҦfCӦfCԳ3C�� C���C���C�ٚC��fC��3C�� C܀ Cݙ�Cޙ�C߳3C�� C�� C�ٚC��fC�� C� C��C�fC�3C�� C���C�ٚC��3C��3C�� C��C�C�fC�3C�� C�� C�ٚC��fC��fC��3C�� C�L�C�ٚD @ D��D� D� D33Ds3D�3D��D
33DffD� D��D33Ds3D��D��D@ D��D��DfDS3D� D�3D�fD  Dy�D �3D"3D#FfD$y�D%� D&�3D(&fD)ffD*��D+�3D-@ D.y�D/�3D0��D2@ D3�fD4�3D5��D7FfD8y�D9�3D:�3D<33D=s3D>��D@fDAL�DB�3DC� DD� DF&fDGl�DH�3DI��DK9�DLy�DM��DO  DPL�DQy�DR�fDT  DU9�DVs3DW��DX�fDZ,�D[s3D\� D^�D_S3D`��Da�fDcfDd@ De� Df�fDhfDiL�Djy�Dk�fDl��Dn33Dos3Dp�3Dq��DsFfDty�Du��Dv��DxL�Dy� Dz�3D{�fD}9�D~��D� D�y�D�3D���D�ffD�3D�� D�<�D��3D�y�D�&fD��3D�c3D�fD��fD�FfD�� D�|�D��D���D�Y�D�  D��fD�FfD��fD�|�D�fD��3D�ffD���D��fD�6fD��fD�|�D�  D��3D�\�D��3D���D�L�D��D���D�)�D�ɚD�ffD���D���D�C3D�� D��3D��D��fD�` D���D��fD�@ D�� D�|�D�  D�� D�c3D�  D�� D�6fD�ٚD�|�D�#3D��fD�` D���D���D�33D���D��fD�#3D��3D�` D�  D��3D�FfD��D���D�  D��3D�VfD��3D���D�<�D�� D�vfD��D�� D�Y�D��fD�� D�<�D��fD�� D�  D��fD�Y�D��fD���D�@ D�ٚD��fD��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A0  A>ffAQ��A`  Ap  A�  A�33A�  A���A�  A���A�ffA�33A���A�33A�33A���A�ffA�33A�  A�33A�ffB  B��B��BffB  B  B  B��B#��B'33B+��B/33B333B733B;��B?33BC��BG��BL  BO��BTffBV��B[33B_33Bc33Bg��Bl  Bp��Bs��Bv��B{��B�  B���B�  B�ffB�  B���B�ffB���B�ffB�33B���B�ffB�  B���B�ffB�  B���B���B�ffB���B���B�ffB�33B�  B���B���B�ffB�33B���B���B�ffB�33B�  B���BÙ�B�ffB�  B�ffB�  Bҙ�B�33B���B�  B♚B�  B�ffB�  B�ffB�  B�ffB�  C33C� C�3C� C	L�C�CffC��CffCL�C�CffC��C��C� CL�C!33C#�C%ffC'�3C)�3C+��C-� C/� C1L�C3L�C533C733C933C;33C=33C?�CA�CCL�CE��CG� CIL�CK��CM� COffCQL�CSL�CUL�CWL�CYL�C[L�C]� C_� Ca�3CcffCe� Cg��CiffCk�CmL�Co� Cq33Cs  CuL�Cw� CyffC{33C}� C�3C���C��3C��fC���C��3C�ٚC���C��3C���C���C��3C���C��3C���C���C��3C���C��3C���C�� C�ٚC���C��3C��fC���C��3C��fC�ٚC���C�� C��3C��3C��fC���C���C�� C�� C��3C��3C��fC��fC��fC��fC��fC�ٚC��fC��fC��3C��3C�� C�� C�� C�� C���C���C���C���C���C��fC��fC��3C��3C�� C���C�ٚC��fC��fC�� CĀ CŌ�Cƙ�CǙ�CȦfCɳ3C�� C�ٚC�ٚC��fCγ3Cπ CЌ�Cљ�CҦfCӦfCԳ3C�� C���C���C�ٚC��fC��3C�� C܀ Cݙ�Cޙ�C߳3C�� C�� C�ٚC��fC�� C� C��C�fC�3C�� C���C�ٚC��3C��3C�� C��C�C�fC�3C�� C�� C�ٚC��fC��fC��3C�� C�L�C�ٚD @ D��D� D� D33Ds3D�3D��D
33DffD� D��D33Ds3D��D��D@ D��D��DfDS3D� D�3D�fD  Dy�D �3D"3D#FfD$y�D%� D&�3D(&fD)ffD*��D+�3D-@ D.y�D/�3D0��D2@ D3�fD4�3D5��D7FfD8y�D9�3D:�3D<33D=s3D>��D@fDAL�DB�3DC� DD� DF&fDGl�DH�3DI��DK9�DLy�DM��DO  DPL�DQy�DR�fDT  DU9�DVs3DW��DX�fDZ,�D[s3D\� D^�D_S3D`��Da�fDcfDd@ De� Df�fDhfDiL�Djy�Dk�fDl��Dn33Dos3Dp�3Dq��DsFfDty�Du��Dv��DxL�Dy� Dz�3D{�fD}9�D~��D� D�y�D�3D���D�ffD�3D�� D�<�D��3D�y�D�&fD��3D�c3D�fD��fD�FfD�� D�|�D��D���D�Y�D�  D��fD�FfD��fD�|�D�fD��3D�ffD���D��fD�6fD��fD�|�D�  D��3D�\�D��3D���D�L�D��D���D�)�D�ɚD�ffD���D���D�C3D�� D��3D��D��fD�` D���D��fD�@ D�� D�|�D�  D�� D�c3D�  D�� D�6fD�ٚD�|�D�#3D��fD�` D���D���D�33D���D��fD�#3D��3D�` D�  D��3D�FfD��D���D�  D��3D�VfD��3D���D�<�D�� D�vfD��D�� D�Y�D��fD�� D�<�D��fD�� D�  D��fD�Y�D��fD���D�@ D�ٚD��fD��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����Ĝ��5?��񪿰���~�������������^5��+�������/���׿|(��vE��s�F��u��7���þ���vɾ�KǾ�=q��+��(�þt�����\�}�Y��t�;��
=�{=��#>��=��#=���<u��vɾ	7L�r�!��hs��z��z��n���dZ�7Kǽ�P=�l�=>!��>fff>�bN>�5?>�p�>�1>��
>�M�>�=q>�V>��>��>��>��D?�?Vȴ?V�+?� �?�o?�z�?�J?h��?Hr�?\j?o�;?o�;?n{?�33?�$�?���?��9?��-?��u?���?��m?�;d?��?���?��?��?�V?�S�?��?��?�&�?�b?��F?��u?�Q�?�&�?���?�?���?�x�?�/?�~�?�=q?���?���?���?�\)?��?���?���?Ƈ+?�\)?��P?�V@	�@
-@�j@{@Z@
^5@\)@�@�9@�@X@�@J@��@r�@v�@/@�@��@��@��@��@j@�h@�R@�@��@o@(�@��@�-@5?@�@�P@  @  @  @ �@�9@�@X@X@��@��@�@n�@n�@=q@x�@hs@x�@hs@x�@x�@��@�u@r�@A�@b@�w@�@��@@t�@X@��@1'@1'@Ĝ@��@r�@@�@
J@	�@�@�+@{@�u@l�@I�@G�@ Ĝ@ 1'?��?�dZ?��H?�X?��u?��u?��u?��?�
=?�7L?���?�"�?���@ ��@��@�H@��@C�@t�@�@��@��@��@��@��@@�7@ �u?�{?��-?�I�?���?��9?�1'?�1'?�Q�?�x�?�7L?��^?��#?��#?��H?���?�b?�b?��P?��y?�ȴ?��y?�ff?��?��?�F?���?��?�O�?��?�1?�=q?�r�?׍P?��?�Z?���?�  ?��m?ě�?���?��R?�|�?�G�?�Z?�z�?�t�?�o?���?�M�?� �?���?�C�?�7L?�J?�j?��#?�&�?�?�$�?��?�`B?|j?q��?j~�?d�?vȴ?�ff?��?v?cS�?bM�?fff?lI�?MO�?D�?>5??2-??ff?G�?�?$�/?��>�D>�%>\>š�>�(�?��>��#>�>�x�>��>���>�!>�{>�Z>��>�=q>��7>��H>��h>�1>��h>�p�>�V>�l�>�{>�Z>�"�>�hs>�hs>�hs>�o>��/>�o>E��>8Q�>,1>+>2->��> Ĝ>,1>&�y>C�=�F=�`B=��=�F=�h=�=�S�=�v�=��=�%=D��=P�`=49X<��<�C�;o<�o<�`B<���<�j<��
<49X;D���ě���C���`B�o���,1�,1�49X�49X�P�`�y�#��%��7L���T��E���vɽ�-������`�������$ݾC��hs����-�"��A�7�Kƨ�V�_;d�dZ�k��p�׾t�j�x���z�H�{�m��  ��  ���7���\�������𾆧𾆧�$ݾ��������bN��zᾘ�u���㾟�w��G��������徤Z��ff��l���1��V���h���h���h��{�����������׾��!��KǾ��پ�Q쾺^5���H��j���۾\�Õ������$ݾ�+��+��1'��1'��=q���;�\)���`��녾�t����ϾՁ�׍P�ۥ�ݲ-�޸R��;d��;d��A���A���Ĝ���
��`B��ff���T���y��~���� ž����&��׾� ž� ž�� ž� ž�� �111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ��Ĝ��5?��񪿰���~�������������^5��+�������/���׿|(��vE��s�F��u��7���þ���vɾ�KǾ�=q��+��(�þt�����\�}�Y��t�;��
=�{=��#>��=��#=���<u��vɾ	7L�r�!��hs��z��z��n���dZ�7Kǽ�P=�l�=>!��>fff>�bN>�5?>�p�>�1>��
>�M�>�=q>�V>��>��>��>��D?�?Vȴ?V�+?� �?�o?�z�?�J?h��?Hr�?\j?o�;?o�;?n{?�33?�$�?���?��9?��-?��u?���?��m?�;d?��?���?��?��?�V?�S�?��?��?�&�?�b?��F?��u?�Q�?�&�?���?�?���?�x�?�/?�~�?�=q?���?���?���?�\)?��?���?���?Ƈ+?�\)?��P?�V@	�@
-@�j@{@Z@
^5@\)@�@�9@�@X@�@J@��@r�@v�@/@�@��@��@��@��@j@�h@�R@�@��@o@(�@��@�-@5?@�@�P@  @  @  @ �@�9@�@X@X@��@��@�@n�@n�@=q@x�@hs@x�@hs@x�@x�@��@�u@r�@A�@b@�w@�@��@@t�@X@��@1'@1'@Ĝ@��@r�@@�@
J@	�@�@�+@{@�u@l�@I�@G�@ Ĝ@ 1'?��?�dZ?��H?�X?��u?��u?��u?��?�
=?�7L?���?�"�?���@ ��@��@�H@��@C�@t�@�@��@��@��@��@��@@�7@ �u?�{?��-?�I�?���?��9?�1'?�1'?�Q�?�x�?�7L?��^?��#?��#?��H?���?�b?�b?��P?��y?�ȴ?��y?�ff?��?��?�F?���?��?�O�?��?�1?�=q?�r�?׍P?��?�Z?���?�  ?��m?ě�?���?��R?�|�?�G�?�Z?�z�?�t�?�o?���?�M�?� �?���?�C�?�7L?�J?�j?��#?�&�?�?�$�?��?�`B?|j?q��?j~�?d�?vȴ?�ff?��?v?cS�?bM�?fff?lI�?MO�?D�?>5??2-??ff?G�?�?$�/?��>�D>�%>\>š�>�(�?��>��#>�>�x�>��>���>�!>�{>�Z>��>�=q>��7>��H>��h>�1>��h>�p�>�V>�l�>�{>�Z>�"�>�hs>�hs>�hs>�o>��/>�o>E��>8Q�>,1>+>2->��> Ĝ>,1>&�y>C�=�F=�`B=��=�F=�h=�=�S�=�v�=��=�%=D��=P�`=49X<��<�C�;o<�o<�`B<���<�j<��
<49X;D���ě���C���`B�o���,1�,1�49X�49X�P�`�y�#��%��7L���T��E���vɽ�-������`�������$ݾC��hs����-�"��A�7�Kƨ�V�_;d�dZ�k��p�׾t�j�x���z�H�{�m��  ��  ���7���\�������𾆧𾆧�$ݾ��������bN��zᾘ�u���㾟�w��G��������徤Z��ff��l���1��V���h���h���h��{�����������׾��!��KǾ��پ�Q쾺^5���H��j���۾\�Õ������$ݾ�+��+��1'��1'��=q���;�\)���`��녾�t����ϾՁ�׍P�ۥ�ݲ-�޸R��;d��;d��A���A���Ĝ���
��`B��ff���T���y��~���� ž����&��׾� ž� ž�� ž� ž�� �111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�B	�B	 �B	$�B	'�B	.B	�B	1'B	0!B	<jB	A�B	@�B	J�B	H�B	\)B	_;B	�bB	��B	�-B	�}B	ǮB	��B	�)B	�)B	�yB
  B
B
\B
{B
�B
"�B
&�B
:^B
;dB
C�B
G�B
P�B
J�B
5?B
E�B
G�B
?}B
C�B
E�B
G�B
L�B
XB
`BB
p�B
}�B
�B
�+B
�7B
�uB
��B
��B
��B
��B
��B
�B
�B
�FB
�!B
��B
�RB
�B
�yB
�B
��B
�B
��B
�B
�`B
�B
��B
��BB%BhB�B{B!�B�BhBbBoB{B�B�B�B�B'�B2-B5?B33B33B/B6FB@�B49B1'B-B0!B1'B1'B7LB5?B33B5?B2-B33B2-B2-B>wBffBZBy�Bt�B�%B�+B��B��B��B��B��B�B�B�B�B�B�-B�9B�9B�3B�FB�9B�3B�3B�FB�LB�LB�RB�qB�}B��BƨB��B��B��B��B��B��B�
B��B�
B�B�B�B�B�B�#B�#B�#B�#B�)B�/B�/B�/B�/B�)B�)B�/B�/B�)B�#B�B�)B�)B�#B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBǮBŢB�}B�wB�}B�qB�dB�jB�dB�dB�jB�^B�FB�jB�wB��B��BÖBȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺB��BɺBȴBȴBǮBǮBɺB��BɺBɺB��BɺB��B��B��B��B��B��B��B��BɺBɺBB�wB�dB�dB�XB�LB�RB�XB�FB�?B�9B�?B�-B�!B�'B��B��B��B��B��B�B�B�B�B�B�B�B��B��B��B��B��B��B�{B�bB�1B�7B�7B�%B�B�B�B�bB��B��B��B�VB�hB�uB��B�VB�DB�=B�+B}�B{�B|�B~�B�B�DBz�Bv�Bw�By�B~�B�1B�+B�%B�%B�1B�7B�JB�JB�JB�JB�DB�7B�7B�1B�+B�=B�VB�{B��B��B��B��B��B��B��B��B�uB�JB�=B�=B�7B�DB�VB�JB�JB�\B�bB�bB�VB�\B�bB�hB�hB�uB�uB�oB�oB�oB�bB�bB�hB�bB�hB�bB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B	#�B	�B	%�B	)�B	,�B	3B	�B	6$B	5B	AeB	F�B	E}B	O�B	M�B	a#B	d6B	�bB	��B	�/B	�}B	̮B	��B	�*B	�*B	�{B
B
B
_B
~B
"�B
'�B
+�B
?dB
@kB
H�B
L�B
U�B
O�B
:DB
J�B
L�B
D�B
H�B
J�B
L�B
Q�B
]B
eJB
u�B
��B
�&B
�4B
�?B
�B
��B
��B
��B
��B
��B
�B
�B
�QB
�+B
ƕB
�^B
�"B
�B
�B
��B
��B
��B
��B
�oB
��B
��B�BB2ByB�B�B&�B�BvBoB~B�B�B �B$�B �B-B7>B:OB8DB8EB4,B;UBE�B9JB66B2B5/B67B68B<]B:OB8EB:QB7=B8DB7<B7>BC�BkyB_0B~�By�B�:B�@B��B��B��B��B��B�B�"B�!B�)B�4B�DB�OB�PB�HB�[B�SB�KB�IB�\B�cB�dB�iBBēBśB��B��B��B��B�
B�B�B� B�B�"B�*B�(B�/B�-B�5B�:B�:B�=B�=B�AB�GB�GB�GB�GB�AB�AB�GB�IB�AB�:B�7B�?B�CB�;B�1B�)B�5B�B��B��B��B�B�B�B�	B��B��B��B��B��B�B�B��B��BʹBĒBÍBĕBB�{B��B�zB�}B�B�uB�_B�BÎBŘBơBȭB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BǨBÎB�{B�|B�kB�bB�lB�pB�[B�VB�LB�WB�DB�:B�AB�B�B�B�B�B�&B�*B�(B�(B�(B�(B�"B�B�B�B��B��B��B��B�sB�DB�JB�JB�6B�,B�&B�+B�vB��B��B��B�jB�}B��B��B�lB�VB�NB�=B�B��B��B�B�3B�XB�B{�B|�B~�B�B�CB�?B�9B�8B�CB�JB�]B�]B�^B�ZB�UB�LB�IB�CB�>B�PB�hB��B��B��B��B��B��B��B��B��B��B�]B�QB�LB�KB�VB�iB�\B�^B�nB�vB�sB�kB�mB�uB�{B�zB��B��B��B��B��B�tB�sB�|B�tB�zB�qB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0049507 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040520220127170405  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144421  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144421  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A0  D��fG�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170405  IP  PSAL            A0  D��fG�O�                