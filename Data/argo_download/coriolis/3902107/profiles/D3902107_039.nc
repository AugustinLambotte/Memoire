CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T15:12:16Z creation; 2020-11-23T11:33:26Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z          :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   L�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       V�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ^�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       `�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   h�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       j�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       s   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   {   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       }$   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   �8   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       �@   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �    HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �$   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �(   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �,   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �0   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �T   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828151216  20220127170814  3902107 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               'A   IF                                  2C  D   ARVOR                           AI2600-18EU007                  5900A04                         844 @�'7��I�1   @�'7��I�@S�O� D�@?��LS8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A333AA��AL��A^ffAt��A���A���A�  A�ffA�33A�  A���A���A�33A���A�33A���A�33A陚A�ffA�  A���B��B��B  B��B��B  B33B��B$ffB(  B+��B0  B3��B733B;��B@ffBC��BG33BK��BPffBS��BW33B\  B_��Bc��BhffBlffBp  Bt  Bw��B{33B33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�33B���B���B���B���B���B���B�  BÙ�Bř�B���B���B���B�33B���B�ffB���B�  B���BꙚB�33B�  B���B���B���CffCffCffCffC	� C� C��CL�CL�C� C33CL�CffC33CffC��C!ffC#� C%�3C'� C)L�C+� C-�3C/� C1L�C3�C5ffC7��C9ffC;33C=� C?�3CA� CCL�CE�CGffCI�3CK� CML�CO�CQffCS�3CU��CWffCYL�C[33C]�C_ffCa��Cc�3Ce��Cg� CiffCkL�Cm33Co33Cq�Cs�Cu�Cw�Cy�C{33C}33C33C���C���C���C���C���C���C���C���C���C���C�� C��3C��3C��fC�ٚC�ٚC���C���C�� C��3C��3C��fC��fC���C���C�� C��fC��fC��fC��3C�� C�� C�� C���C��fC��fC��3C���C��fC��3C���C��3C�� C��fC�� C���C��3C�� C�ٚC��3C�� C���C��3C���C��fC�� C���C��3C���C��3C���C��fC�� C��fC�� C�ٚC�� CÙ�C�� C�ٚC�� CǦfC�� C�ٚC�� C˦fČ�Cͳ3C�ٚC�� Cг3Cљ�C�� C�ٚC���Cճ3C֦fCי�C�� C��fC���C�� Cܳ3Cݙ�Cތ�C�� C��fC�ٚC���C�� C�3C�fC晚C��C� C�3C��fC�ٚC�ٚC���C�� C�3C�3C�fC�C�C��C��3C�� C���C���C���C�@ C���D ,�Dl�D�3D��DFfD� D��D�3D
33Ds3D��D�3D33Ds3D�3D�3D,�Ds3D��DfDL�Ds3D��D  D@ D� D ��D!��D#L�D$��D%��D&�fD(,�D)y�D*�fD,  D-9�D.y�D/��D1  D2,�D3s3D4� D5��D733D8� D9��D:��D<L�D=�fD>�3D?��DA@ DB�fDC� DE�DF9�DGl�DH��DJ�DKL�DLy�DM�fDN��DP33DQl�DR� DS��DU33DVl�DW�3DX��DZ,�D[� D\��D^fD_9�D`l�Da�fDb��DdFfDey�Df��Dh  Di33DjffDk��Dm�DnL�Do�3Dp��Dq�fDs,�Dty�Du�fDv��Dx,�Dyy�Dz�fD{��D},�D~y�D�fD��fD�  D��fD�\�D�  D��fD�<�D��fD�� D��D��fD�` D�  D���D�<�D���D�|�D�  D�� D�S3D���D�� D�<�D�ٚD�vfD�fD���D�\�D�  D��3D�<�D��3D�y�D�  D���D�VfD���D���D�9�D�ٚD�y�D��D�� D�ffD���D�� D�<�D���D�|�D�#3D���D�S3D���D���D�9�D�ٚD�y�D��D���D�c3D���D��fD�@ D��3D�� D�  D�� D�` D�3D���D�33D���D��fD�#3D��3D�c3D�  D��3D�6fD�ٚD�� D��D��3D�S3D�� D�� D�0 D��3D�s3D�3D��3D�S3D��3D��fD�9�D���D�� D�&fD���D�P D��fD���D�FfD�� D�y�D�&fD�� D�\�D��fDĐ D�<�D��fD�vfD�#3D�� D�` D�3Dɣ3D�C3D��3DˆfD��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A333AA��AL��A^ffAt��A���A���A�  A�ffA�33A�  A���A���A�33A���A�33A���A�33A陚A�ffA�  A���B��B��B  B��B��B  B33B��B$ffB(  B+��B0  B3��B733B;��B@ffBC��BG33BK��BPffBS��BW33B\  B_��Bc��BhffBlffBp  Bt  Bw��B{33B33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�33B���B���B���B���B���B���B�  BÙ�Bř�B���B���B���B�33B���B�ffB���B�  B���BꙚB�33B�  B���B���B���CffCffCffCffC	� C� C��CL�CL�C� C33CL�CffC33CffC��C!ffC#� C%�3C'� C)L�C+� C-�3C/� C1L�C3�C5ffC7��C9ffC;33C=� C?�3CA� CCL�CE�CGffCI�3CK� CML�CO�CQffCS�3CU��CWffCYL�C[33C]�C_ffCa��Cc�3Ce��Cg� CiffCkL�Cm33Co33Cq�Cs�Cu�Cw�Cy�C{33C}33C33C���C���C���C���C���C���C���C���C���C���C�� C��3C��3C��fC�ٚC�ٚC���C���C�� C��3C��3C��fC��fC���C���C�� C��fC��fC��fC��3C�� C�� C�� C���C��fC��fC��3C���C��fC��3C���C��3C�� C��fC�� C���C��3C�� C�ٚC��3C�� C���C��3C���C��fC�� C���C��3C���C��3C���C��fC�� C��fC�� C�ٚC�� CÙ�C�� C�ٚC�� CǦfC�� C�ٚC�� C˦fČ�Cͳ3C�ٚC�� Cг3Cљ�C�� C�ٚC���Cճ3C֦fCי�C�� C��fC���C�� Cܳ3Cݙ�Cތ�C�� C��fC�ٚC���C�� C�3C�fC晚C��C� C�3C��fC�ٚC�ٚC���C�� C�3C�3C�fC�C�C��C��3C�� C���C���C���C�@ C���D ,�Dl�D�3D��DFfD� D��D�3D
33Ds3D��D�3D33Ds3D�3D�3D,�Ds3D��DfDL�Ds3D��D  D@ D� D ��D!��D#L�D$��D%��D&�fD(,�D)y�D*�fD,  D-9�D.y�D/��D1  D2,�D3s3D4� D5��D733D8� D9��D:��D<L�D=�fD>�3D?��DA@ DB�fDC� DE�DF9�DGl�DH��DJ�DKL�DLy�DM�fDN��DP33DQl�DR� DS��DU33DVl�DW�3DX��DZ,�D[� D\��D^fD_9�D`l�Da�fDb��DdFfDey�Df��Dh  Di33DjffDk��Dm�DnL�Do�3Dp��Dq�fDs,�Dty�Du�fDv��Dx,�Dyy�Dz�fD{��D},�D~y�D�fD��fD�  D��fD�\�D�  D��fD�<�D��fD�� D��D��fD�` D�  D���D�<�D���D�|�D�  D�� D�S3D���D�� D�<�D�ٚD�vfD�fD���D�\�D�  D��3D�<�D��3D�y�D�  D���D�VfD���D���D�9�D�ٚD�y�D��D�� D�ffD���D�� D�<�D���D�|�D�#3D���D�S3D���D���D�9�D�ٚD�y�D��D���D�c3D���D��fD�@ D��3D�� D�  D�� D�` D�3D���D�33D���D��fD�#3D��3D�c3D�  D��3D�6fD�ٚD�� D��D��3D�S3D�� D�� D�0 D��3D�s3D�3D��3D�S3D��3D��fD�9�D���D�� D�&fD���D�P D��fD���D�FfD�� D�y�D�&fD�� D�\�D��fDĐ D�<�D��fD�vfD�#3D�� D�` D�3Dɣ3D�C3D��3DˆfD��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��j@�1'@��
@�|�@�S�@���@��\@��T@�G�@�Ĝ@� �@�+@�{@��P@�t�@���@��9@��@�+@���@��\@�{@f��@a�#@T�/@S�F@X1'@^v�@g�@W��@S�@Q��@J~�@G+@J�\@J��@G�@HbN@H�9@H��@G�P@A�@?l�@@1'@<(�@<��@=/@;��@<�D@<1@:��@:�@;o@;33@:�@:-@9x�@8��@5�-@333@1hs@0bN@0  @0Q�@1��@4�D@333@3@.v�@/�P@/l�@.�y@.�R@.E�@-�T@-?}@+�F@*�H@,�@.�y@-�T@,(�@,�@3@:-@A��@G|�@E�@D9X@Cƨ@C�@Bn�@@bN@?K�@>E�@;��@8A�@6��@6ff@6�R@6��@4�/@0��@-�@1��@1�7@2^5@1�#@1�@2��@2�\@1X@/��@.�@-��@,Z@+ƨ@+"�@*n�@*J@+@*�@(Ĝ@(�@(Ĝ@&5?@&�+@&�+@&ȴ@&�@'�@&��@&�+@&V@&��@'+@'�@'��@'�@'�P@'K�@&�@&E�@%�T@%�-@&�+@%�-@%�@%��@ ��@33@��@|�@;d@"�@n�@J@�;@ff@O�@j@��@+@j@ A�@ Ĝ@;d@$�@K�@ƨ@M�@�@x�@�#@�@hs@��@  @�@��@?}@�@
-@�/@X@ A�?���?�7L@ 1'@��@�?�Ĝ?��m?�?�?�{?�"�?���?ٙ�?�~�?��m?׮?׍P?ա�?�o?ӶF?�l�?�X?�I�?��?�?陚?���?�"�?�~�?��#?�u?�Q�?�P?�E�?�Z?�M�?�J?�J?�hs?�A�?޸R?�p�?�I�?���?�r�?��T?ԛ�?�9X?�t�?�n�?�\)?�5??Ͳ-?�p�?̬?�j?�1?���?ȴ9?�ȴ?�Z?�S�?�n�?�-?��7?�bN?�O�?��D?�=q?��P?�E�?��
?�n�?���?��?��h?�I�?���?���?�r�?��u?���?���?�^5?��#?��?��^?��?��+?��?��?�A�?�  ?�  ?�\)?��R?�5??�5??�{?�I�?�~�?�l�?�ff?���?���?�Q�?rn�?`�?Y��?["�?]/?bM�?VE�?X��?U?}?Kƨ?:^5?,��?$�/?"��?X? �?1?�7>�>��7>��>�M�>�=q>l�D>H�9><j>!��>   =�
==���=\=�t�=q��=]/=L��=8Q�=<j=H�9<�h<�C�;�`B;�`B<�j=�w=C�<T����o��/��h��1�t����t��0 Ž�%��O߽����9X������xվo�
=q�O߾�P��-�"��-V�?|�O�;�Z��_;d�ixվs�F�~�۾��7�����+���9���;�hs������P�����Ĝ���y�������F���#���۾���Ƨ���;�և+�����Ĝ��Z���y��xվ�V��&��F��X���ۿG��Mӿo���ff�	��C����I���h������녿t��?}�E���P�b��u���dZ�j�vɿ   � Ĝ�!G��!���"�\�"��"��"��#S��$��%��&ff�'l��'(�9�)�^�*~��+C��+ƨ�,1�,I��,�D�.{�/���0 ſ0�׿1녿2n��2�!�2�!�3�Ͽ49X�4�j�5?}�5�6�+�7
=�7�P�8Q�9��9X�9���9�#�:��:^5�:���:�H�;"ѿ;dZ�<(��=/�=p��=�-�=�>vɿ>�R�>�ۿ?|�?�w�@  �@A��@��@��@Ĝ�A%�A�7�A�7�A�7�A���B�\�Co�C���C���C���D��D��DZ�DZ�D��D��D��D��C�
�D��D��D��DZ�DZ�D�/�D�/�D�/�E��E��D�/�D�/�D�/�D�/�D���D��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��j@�1'@��
@�|�@�S�@���@��\@��T@�G�@�Ĝ@� �@�+@�{@��P@�t�@���@��9@��@�+@���@��\@�{@f��@a�#@T�/@S�F@X1'@^v�@g�@W��@S�@Q��@J~�@G+@J�\@J��@G�@HbN@H�9@H��@G�P@A�@?l�@@1'@<(�@<��@=/@;��@<�D@<1@:��@:�@;o@;33@:�@:-@9x�@8��@5�-@333@1hs@0bN@0  @0Q�@1��@4�D@333@3@.v�@/�P@/l�@.�y@.�R@.E�@-�T@-?}@+�F@*�H@,�@.�y@-�T@,(�@,�@3@:-@A��@G|�@E�@D9X@Cƨ@C�@Bn�@@bN@?K�@>E�@;��@8A�@6��@6ff@6�R@6��@4�/@0��@-�@1��@1�7@2^5@1�#@1�@2��@2�\@1X@/��@.�@-��@,Z@+ƨ@+"�@*n�@*J@+@*�@(Ĝ@(�@(Ĝ@&5?@&�+@&�+@&ȴ@&�@'�@&��@&�+@&V@&��@'+@'�@'��@'�@'�P@'K�@&�@&E�@%�T@%�-@&�+@%�-@%�@%��@ ��@33@��@|�@;d@"�@n�@J@�;@ff@O�@j@��@+@j@ A�@ Ĝ@;d@$�@K�@ƨ@M�@�@x�@�#@�@hs@��@  @�@��@?}@�@
-@�/@X@ A�?���?�7L@ 1'@��@�?�Ĝ?��m?�?�?�{?�"�?���?ٙ�?�~�?��m?׮?׍P?ա�?�o?ӶF?�l�?�X?�I�?��?�?陚?���?�"�?�~�?��#?�u?�Q�?�P?�E�?�Z?�M�?�J?�J?�hs?�A�?޸R?�p�?�I�?���?�r�?��T?ԛ�?�9X?�t�?�n�?�\)?�5??Ͳ-?�p�?̬?�j?�1?���?ȴ9?�ȴ?�Z?�S�?�n�?�-?��7?�bN?�O�?��D?�=q?��P?�E�?��
?�n�?���?��?��h?�I�?���?���?�r�?��u?���?���?�^5?��#?��?��^?��?��+?��?��?�A�?�  ?�  ?�\)?��R?�5??�5??�{?�I�?�~�?�l�?�ff?���?���?�Q�?rn�?`�?Y��?["�?]/?bM�?VE�?X��?U?}?Kƨ?:^5?,��?$�/?"��?X? �?1?�7>�>��7>��>�M�>�=q>l�D>H�9><j>!��>   =�
==���=\=�t�=q��=]/=L��=8Q�=<j=H�9<�h<�C�;�`B;�`B<�j=�w=C�<T����o��/��h��1�t����t��0 Ž�%��O߽����9X������xվo�
=q�O߾�P��-�"��-V�?|�O�;�Z��_;d�ixվs�F�~�۾��7�����+���9���;�hs������P�����Ĝ���y�������F���#���۾���Ƨ���;�և+�����Ĝ��Z���y��xվ�V��&��F��X���ۿG��Mӿo���ff�	��C����I���h������녿t��?}�E���P�b��u���dZ�j�vɿ   � Ĝ�!G��!���"�\�"��"��"��#S��$��%��&ff�'l��'(�9�)�^�*~��+C��+ƨ�,1�,I��,�D�.{�/���0 ſ0�׿1녿2n��2�!�2�!�3�Ͽ49X�4�j�5?}�5�6�+�7
=�7�P�8Q�9��9X�9���9�#�:��:^5�:���:�H�;"ѿ;dZ�<(��=/�=p��=�-�=�>vɿ>�R�>�ۿ?|�?�w�@  �@A��@��@��@Ĝ�A%�A�7�A�7�A�7�A���B�\�Co�C���C���C���D��D��DZ�DZ�D��D��D��D��C�
�D��D��D��DZ�DZ�D�/�D�/�D�/�E��E��D�/�D�/�D�/�D�/�D���D��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B�uB�bB�VB�oB�Bt�Br�BD�BI�B:^B5?B?}B^5BVBVBN�B[#BF�BJ�BQ�BXBYB_;BaHBdZBffBgmBjBp�Bt�Bx�Bv�Bx�B|�B~�B�B�B�B�B�JB�PB�VB�VB�JB�+B�+B�B�%B�7B�=B�hB�DB�hB�\B��B��B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B�}B��B�B�B�
B��B�B�B�B�B�B��B��B��B��B�B��B��B��B��B��B��B��B�
B�
B�B�B�B�B��B�B��B��B��B��B��B�B�B��B�B��B��B��B��B��B�
B�B�B�B�#B�#B�#B�;B�5B�BB�;B�;B�;B�5B�;B�5B�NB�HB�BB�/B��B��B��B��BɺB��B�B��B��B��B��B��B�B��B�/B�`B�ZB�TB�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�jB�qB��B��B��B�wB�RB�9B�-B�'B�B��B��B��B��B��B��B��B��B�B�B�B�B�FB�qB��BÖBĜBĜBƨBĜBŢBĜBŢBĜBƨBÖBĜBŢBĜBBBÖB��B��B��B�}B�}B��B��B��B�wB�}B�}B�wB�wB�}B�wB�wB�qB�}B��B�qB�qB�wB�wB�wB�qB�jB�dB�dB�dB�XB�RB�FB�LB�XB�?B�FB�LB�FB�LB�RB�^B�^B�dB��B��B�}B�wB�qB�qB�qB�wB�}B�}B�}B��B��B�}B��B�}B��B��B��B�jB�RB�3B��B�B��B�'B�9B�3B�9B�LB�?B�-B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B�uB�bB�VB�oB�Bt�Br�BD�BI�B:^B5?B?}B^5BVBVBN�B[#BF�BJ�BQ�BXBYB_;BaHBdZBffBgmBjBp�Bt�Bx�Bv�Bx�B|�B~�B�B�B�B�B�JB�PB�VB�VB�JB�+B�+B�B�%B�7B�=B�hB�DB�hB�\B��B��B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B�}B��B�B�B�
B��B�B�B�B�B�B��B��B��B��B�B��B��B��B��B��B��B��B�
B�
B�B�B�B�B��B�B��B��B��B��B��B�B�B��B�B��B��B��B��B��B�
B�B�B�B�#B�#B�#B�;B�5B�BB�;B�;B�;B�5B�;B�5B�NB�HB�BB�/B��B��B��B��BɺB��B�B��B��B��B��B��B�B��B�/B�`B�ZB�TB�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�jB�qB��B��B��B�wB�RB�9B�-B�'B�B��B��B��B��B��B��B��B��B�B�B�B�B�FB�qB��BÖBĜBĜBƨBĜBŢBĜBŢBĜBƨBÖBĜBŢBĜBBBÖB��B��B��B�}B�}B��B��B��B�wB�}B�}B�wB�wB�}B�wB�wB�qB�}B��B�qB�qB�wB�wB�wB�qB�jB�dB�dB�dB�XB�RB�FB�LB�XB�?B�FB�LB�FB�LB�RB�^B�^B�dB��B��B�}B�wB�qB�qB�qB�wB�}B�}B�}B��B��B�}B��B�}B��B��B��B�jB�RB�3B��B�B��B�'B�9B�3B�9B�LB�?B�-B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               No adjustment was necessary. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                                                      202011231133262022012717081420220127170814  IF  ARFMCODA035h                                                                20200828151216                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828151253  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828151253  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201123113326  IP  PSAL            A333D��G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170814  IP  PSAL            A333D��G�O�                