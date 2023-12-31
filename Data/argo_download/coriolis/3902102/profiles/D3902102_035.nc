CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  M   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:22Z creation; 2020-11-17T12:19:02Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        	4  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  C�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	4  E�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  O(   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	4  Qx   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  Z�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  c�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  f0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  od   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  q�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  z�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  �l   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 P  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	4  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
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
_FillValue                  ,  �TArgo profile    3.1 1.2 19500101000000  20200828144322  20220127170400  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               #A   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @��֌So�1   @��֌So�@Sa�5����(�±��8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 10.0 dbar]                                                     A&ffA1��AA��AP  A^ffAnffA|��A�  A�ffA�ffA�  A�33A�ffA�  A�33A�ffA�  A�  A�33A�33A�ffA�  A�33B��BffB  B  B��B��B33B33B#33B'33B+33B/33B333B733B;��B?33BC��BG��BL  BP  BTffBW33B[33B_��Bc��Bh  Bl  Bp  Bt  BxffB{33B33B���B���B�  B�  B�33B���B���B�  B���B�  B�  B�33B���B���B���B���B���B�  B�33B���B���B���B���B�  B�  B�33B���B���B�  B�  B�33B���B���B�  B�  B�ffB�33B���B�ffB���B�33B�  B♚B�33B�  B���B�B���B�33B�33C��C��C33CL�C	L�CL�CffCffC� C33CL�CffC� C�3CffC� C!�3C#ffC%�C'L�C)� C+L�C-� C/�3C1� C333C5� C7��C9� C;L�C=�C?ffCA��CCffCE33CG�CIffCK��CM� COL�CQ33CS�CUffCW��CY�3C[��C]ffC_L�Ca33Cc�CeffCg��Ci��Ck��Cm� CoffCqL�Cs33Cu�Cw  CyffC{�fC}�fC��C��fC�ٚC���C���C�� C�� C�� C��3C��3C��3C��3C��3C��3C�� C�� C�� C���C�ٚC��fC��3C�� C���C��fC��3C���C�ٚC��fC�� C�� C���C��3C���C��fC��3C���C���C��3C���C��3C�� C��fC�� C�ٚC�� C���C��fC�� C�ٚC��fC�� C���C���C��3C�� C���C���C�ٚC��fC��3C�� C�� C���C���C��3C�� C���C�ٚC��fCĳ3Cŀ Cƀ CǙ�CȦfC�� C���C�ٚC��3C�� C΀ CϦfCг3Cѳ3C�� C�� C���Cՙ�C֦fC׳3Cس3C�� C���C�ٚCܦfCݳ3C�� Cߙ�C�fC�3C�3C�� C䙚C�fC�3C�3C���C陚C�fC�3C�� C�ٚC�fC�3C�� C�fC�fC�� C���C��fC�� C���C��fC�� C�� C��fD ,�D� D�3D��D33Ds3D�3D�3D
33Ds3D��D  DFfDs3D� D�3DFfD� D��D�3D,�Dl�D��D�3D@ D�fD ��D!��D#9�D$� D%��D'  D(@ D)�fD*��D,  D-33D.�fD/� D0��D29�D3s3D4��D5�3D79�D8� D9��D;  D<33D=` D>��D@�DAFfDBy�DC��DD��DF@ DG� DH�fDJ�DKS3DL� DM� DN��DP33DQ� DR��DS��DU&fDV� DW��DX��DZ&fD[y�D\��D]��D_&fD`s3Da� Db�3Dd,�De� Df� Dg��Di@ Djl�Dk� Dl�3Dn@ Dos3Dp��Dq�fDs9�Dt��Du��Dv� Dx&fDy�fDz�fD{��D},�D~ffD� D�p D� D��3D�VfD���D�� D�9�D��fD�� D��D��fD�VfD��fD��fD�9�D���D�s3D��D�� D�Y�D��fD��3D�33D��3D�s3D�3D��fD�Y�D���D���D�9�D�� D�y�D�3D���D�ffD�fD��fD�I�D���D�s3D�fD���D�c3D���D���D�C3D�� D�|�D�  D���D�` D�  D��fD�<�D��fD�� D�&fD��3D�` D���D���D�6fD��fD�s3D�fD��fD�\�D�  D��3D�I�D���D�s3D�#3D��3D�\�D���D���D�6fD��fD�y�D��D�� D�ffD���D��fD�<�D��3D��fD�  D��fD�Y�D�  D��fD�@ D���D�y�D�fD��3D�ffD���D���D�C3D�ٚD�� D�&fD�� D�VfD��3DĜ�D�I�D��DƆfD�&fD��3D�c3D���Də�D�<�D��fD�|�D�&fD��fD�ffD�	�DΩ�D�<�D��3D�vfD�fDѼ�D�c3D���Dә�D�C3D�� D�|�D��D�� D�c3D���Dؠ D�FfD���D�s3D��D�ɚD�c3D�3Dݠ D�@ D�� D�|�D�fD๚D�\�D��fD��D�I�D���D�s3D�fD幚D�c3D���D癚D�9�D�ٚD�y�D�  D��fD�` D���D�3D�33D�� D�l�D� D� D�P D��3D�fD�<�D��D�fD�)�D�� D�VfD���D��3D�<�D�ٚD�L�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A&ffA1��AA��AP  A^ffAnffA|��A�  A�ffA�ffA�  A�33A�ffA�  A�33A�ffA�  A�  A�33A�33A�ffA�  A�33B��BffB  B  B��B��B33B33B#33B'33B+33B/33B333B733B;��B?33BC��BG��BL  BP  BTffBW33B[33B_��Bc��Bh  Bl  Bp  Bt  BxffB{33B33B���B���B�  B�  B�33B���B���B�  B���B�  B�  B�33B���B���B���B���B���B�  B�33B���B���B���B���B�  B�  B�33B���B���B�  B�  B�33B���B���B�  B�  B�ffB�33B���B�ffB���B�33B�  B♚B�33B�  B���B�B���B�33B�33C��C��C33CL�C	L�CL�CffCffC� C33CL�CffC� C�3CffC� C!�3C#ffC%�C'L�C)� C+L�C-� C/�3C1� C333C5� C7��C9� C;L�C=�C?ffCA��CCffCE33CG�CIffCK��CM� COL�CQ33CS�CUffCW��CY�3C[��C]ffC_L�Ca33Cc�CeffCg��Ci��Ck��Cm� CoffCqL�Cs33Cu�Cw  CyffC{�fC}�fC��C��fC�ٚC���C���C�� C�� C�� C��3C��3C��3C��3C��3C��3C�� C�� C�� C���C�ٚC��fC��3C�� C���C��fC��3C���C�ٚC��fC�� C�� C���C��3C���C��fC��3C���C���C��3C���C��3C�� C��fC�� C�ٚC�� C���C��fC�� C�ٚC��fC�� C���C���C��3C�� C���C���C�ٚC��fC��3C�� C�� C���C���C��3C�� C���C�ٚC��fCĳ3Cŀ Cƀ CǙ�CȦfC�� C���C�ٚC��3C�� C΀ CϦfCг3Cѳ3C�� C�� C���Cՙ�C֦fC׳3Cس3C�� C���C�ٚCܦfCݳ3C�� Cߙ�C�fC�3C�3C�� C䙚C�fC�3C�3C���C陚C�fC�3C�� C�ٚC�fC�3C�� C�fC�fC�� C���C��fC�� C���C��fC�� C�� C��fD ,�D� D�3D��D33Ds3D�3D�3D
33Ds3D��D  DFfDs3D� D�3DFfD� D��D�3D,�Dl�D��D�3D@ D�fD ��D!��D#9�D$� D%��D'  D(@ D)�fD*��D,  D-33D.�fD/� D0��D29�D3s3D4��D5�3D79�D8� D9��D;  D<33D=` D>��D@�DAFfDBy�DC��DD��DF@ DG� DH�fDJ�DKS3DL� DM� DN��DP33DQ� DR��DS��DU&fDV� DW��DX��DZ&fD[y�D\��D]��D_&fD`s3Da� Db�3Dd,�De� Df� Dg��Di@ Djl�Dk� Dl�3Dn@ Dos3Dp��Dq�fDs9�Dt��Du��Dv� Dx&fDy�fDz�fD{��D},�D~ffD� D�p D� D��3D�VfD���D�� D�9�D��fD�� D��D��fD�VfD��fD��fD�9�D���D�s3D��D�� D�Y�D��fD��3D�33D��3D�s3D�3D��fD�Y�D���D���D�9�D�� D�y�D�3D���D�ffD�fD��fD�I�D���D�s3D�fD���D�c3D���D���D�C3D�� D�|�D�  D���D�` D�  D��fD�<�D��fD�� D�&fD��3D�` D���D���D�6fD��fD�s3D�fD��fD�\�D�  D��3D�I�D���D�s3D�#3D��3D�\�D���D���D�6fD��fD�y�D��D�� D�ffD���D��fD�<�D��3D��fD�  D��fD�Y�D�  D��fD�@ D���D�y�D�fD��3D�ffD���D���D�C3D�ٚD�� D�&fD�� D�VfD��3DĜ�D�I�D��DƆfD�&fD��3D�c3D���Də�D�<�D��fD�|�D�&fD��fD�ffD�	�DΩ�D�<�D��3D�vfD�fDѼ�D�c3D���Dә�D�C3D�� D�|�D��D�� D�c3D���Dؠ D�FfD���D�s3D��D�ɚD�c3D�3Dݠ D�@ D�� D�|�D�fD๚D�\�D��fD��D�I�D���D�s3D�fD幚D�c3D���D癚D�9�D�ٚD�y�D�  D��fD�` D���D�3D�33D�� D�l�D� D� D�P D��3D�fD�<�D��D�fD�)�D�� D�VfD���D��3D�<�D�ٚD�L�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?��-?��P@ �u@@�`@|�@
=@C�@�y@( �@2�@=V@A%@B�\@I�#@Kƨ@M��@Nv�@N$�@N�R@M@M��@MO�@K�m@L�@M�@N@M�T@M`B@M/@MV@L�@JJ@I�@J^5@J�@K"�@I��@I�7@I��@H��@H�@H �@HA�@Hb@G�;@Hr�@L�/@M��@N{@N�@Ol�@Ol�@L�j@K��@J��@IX@H �@H  @IX@J��@K@K�@K��@J��@I��@Hr�@Hb@G+@F@E�T@Ep�@D�/@C��@CS�@B�!@BJ@@b@<�j@:�\@9��@8�9@8�9@6V@5�@49X@3��@1�7@/�;@/��@/�@1��@3S�@4��@5��@6V@6�y@7;d@7\)@7+@6��@4��@4�@3�F@1��@/\)@/
=@0A�@/�w@/|�@-O�@*J@)�7@(�`@(�@'�@'�P@'+@&v�@$Z@#"�@"��@"�@"��@&��@(Ĝ@)G�@)�7@)��@)��@'+@#t�@!�@!X@   @�@I�@�!@�\@��@��@�@�D@��@�D@ƨ@M�@��@A�@�u@&�@��@b@�P@K�@�y@�+@V@E�@5?@`B@��@X@1'@�;@��@�@E�@`B@��@�D@��@1@��@
��@	��@��@bN@  @;d@�@S�@7L?�v�?���?��m?���?��H?���?��u?��/?�t�?�t�?�S�?�33?�J?���?�/?�I�?���?�"�?��#?�j?�"�?�^?�?���?�n�?�7?�&�?� �?�v�?�1'?�?���?Ѓ?�v�?�j?ʟ�?�x�?���?ȴ9?�V?�;d?�A�?��`?���?���?�  ?�\)?ʟ�?ɺ^?�1'?�ȴ?�ȴ?�%?�Ĝ?���?�G�?���?�5??�p�?�V?��?�1?��^?��?�l�?���?���?�33?��?��?��?���?�?�Z?���?��!?�J?���?�  ?���?�5??��h?�(�?���?��^?�+?�?�Z?���?�M�?���?��7?�Ĝ?�|�?�V?�p�?��D?��?�=q?�x�?�x�?���?���?��T?��j?�z�?���?���?���?|�?z^5?q�?e��?cS�?]/?St�?N��?Gl�?@�?9�?0��?,1?)�^?-O�?.V?(r�?#o? �? A�?�R?5??^5?V?�? Ĝ>�33>��y>��y>ܬ>�I�>��F>��y>�;d>�"�>�
=>��`>�7L>{�m>cS�>Xb>D��>:^5>,1>O�=��#=�`B=��`=�j=�hs=�o=]/=49X='�=#�
<��<�1<D��    ���
�D�����
��h�+��P�@��]/�y�#��+��t����㽬1�\�����G�������\)��P����w�$�/�(�þ/��49X�6E��9X�>vɾB�\�G��L�;P�`�R�V�["Ѿ]/�`A��cS��gl��k��n���q���t�j�u�x���|푾����J��o������˾���1'���9��=q��ƨ��O߾�V��\)���`��񪾔�������
=��b������"Ѿ�(����-���R��;d��A�������MӾ��
��Z���T���xվ�~���1���h������ ž������!���j��ȴ������X���#��dZ��p���  �\����Ƨ�Ǯ�ɺ^��ƨ���;����bN��n���z�Ձ�׍P��b�����"Ѿ�/��/�޸R�߾w��A���MӾ�Z���y����þ�~���D��h�����׾����!��33���j��ȴ��Q��X��dZ��j��p����۾�vɾ�|�%������J�Mӿ��o��������/������/����T�	��	7L�r��1'�r��	xտ
~��I��bN� ſ&�hs�-�33������t��n��33�n��33��+�
=�Q�������������������������/��-�5?��R�����(�����m�/��R��w��w�|�|�"J�"�\�"��"��#o�#���$��#�
�$�/�$���%��%�T�&ff�&�y�&�y�'��(r��(�9�(�ÿ(�ÿ)xտ)xտ)�^�*=q�*=q�*���+C��+��+ƨ�,1�,�D�,��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ?��-?��P@ �u@@�`@|�@
=@C�@�y@( �@2�@=V@A%@B�\@I�#@Kƨ@M��@Nv�@N$�@N�R@M@M��@MO�@K�m@L�@M�@N@M�T@M`B@M/@MV@L�@JJ@I�@J^5@J�@K"�@I��@I�7@I��@H��@H�@H �@HA�@Hb@G�;@Hr�@L�/@M��@N{@N�@Ol�@Ol�@L�j@K��@J��@IX@H �@H  @IX@J��@K@K�@K��@J��@I��@Hr�@Hb@G+@F@E�T@Ep�@D�/@C��@CS�@B�!@BJ@@b@<�j@:�\@9��@8�9@8�9@6V@5�@49X@3��@1�7@/�;@/��@/�@1��@3S�@4��@5��@6V@6�y@7;d@7\)@7+@6��@4��@4�@3�F@1��@/\)@/
=@0A�@/�w@/|�@-O�@*J@)�7@(�`@(�@'�@'�P@'+@&v�@$Z@#"�@"��@"�@"��@&��@(Ĝ@)G�@)�7@)��@)��@'+@#t�@!�@!X@   @�@I�@�!@�\@��@��@�@�D@��@�D@ƨ@M�@��@A�@�u@&�@��@b@�P@K�@�y@�+@V@E�@5?@`B@��@X@1'@�;@��@�@E�@`B@��@�D@��@1@��@
��@	��@��@bN@  @;d@�@S�@7L?�v�?���?��m?���?��H?���?��u?��/?�t�?�t�?�S�?�33?�J?���?�/?�I�?���?�"�?��#?�j?�"�?�^?�?���?�n�?�7?�&�?� �?�v�?�1'?�?���?Ѓ?�v�?�j?ʟ�?�x�?���?ȴ9?�V?�;d?�A�?��`?���?���?�  ?�\)?ʟ�?ɺ^?�1'?�ȴ?�ȴ?�%?�Ĝ?���?�G�?���?�5??�p�?�V?��?�1?��^?��?�l�?���?���?�33?��?��?��?���?�?�Z?���?��!?�J?���?�  ?���?�5??��h?�(�?���?��^?�+?�?�Z?���?�M�?���?��7?�Ĝ?�|�?�V?�p�?��D?��?�=q?�x�?�x�?���?���?��T?��j?�z�?���?���?���?|�?z^5?q�?e��?cS�?]/?St�?N��?Gl�?@�?9�?0��?,1?)�^?-O�?.V?(r�?#o? �? A�?�R?5??^5?V?�? Ĝ>�33>��y>��y>ܬ>�I�>��F>��y>�;d>�"�>�
=>��`>�7L>{�m>cS�>Xb>D��>:^5>,1>O�=��#=�`B=��`=�j=�hs=�o=]/=49X='�=#�
<��<�1<D��    ���
�D�����
��h�+��P�@��]/�y�#��+��t����㽬1�\�����G�������\)��P����w�$�/�(�þ/��49X�6E��9X�>vɾB�\�G��L�;P�`�R�V�["Ѿ]/�`A��cS��gl��k��n���q���t�j�u�x���|푾����J��o������˾���1'���9��=q��ƨ��O߾�V��\)���`��񪾔�������
=��b������"Ѿ�(����-���R��;d��A�������MӾ��
��Z���T���xվ�~���1���h������ ž������!���j��ȴ������X���#��dZ��p���  �\����Ƨ�Ǯ�ɺ^��ƨ���;����bN��n���z�Ձ�׍P��b�����"Ѿ�/��/�޸R�߾w��A���MӾ�Z���y����þ�~���D��h�����׾����!��33���j��ȴ��Q��X��dZ��j��p����۾�vɾ�|�%������J�Mӿ��o��������/������/����T�	��	7L�r��1'�r��	xտ
~��I��bN� ſ&�hs�-�33������t��n��33�n��33��+�
=�Q�������������������������/��-�5?��R�����(�����m�/��R��w��w�|�|�"J�"�\�"��"��#o�#���$��#�
�$�/�$���%��%�T�&ff�&�y�&�y�'��(r��(�9�(�ÿ(�ÿ)xտ)xտ)�^�*=q�*=q�*���+C��+��+ƨ�,1�,�D�,��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�jB	�fB

=B
.B
D�B
jB
��B
�XB
��B1BL�B� B��B�-B�HB�B��B��B��B��B��B��BB��B��BBBBBBBB��B��B��BBBB��B  B��B��B��B��B��B  B��BBJBPB\BDBJB+BB%BBBB+B
=BDBJBPBPB
=B+B1B+B%BBBB��B��B��B  B��B��B�B�B�B�B��B�B�B�sB�fB�sB�mB�mB�B�B�B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�fB�mB�`B�`B�ZB�ZB�ZB�BB�;B�BB�BB�NB�B�B�B�B�B�B�B�sB�fB�TB�`B�;B�BB�5B�5B�HB�HB�TB�NB�TB�NB�HB�NB�BB�NB�NB�sB�`B�`B�fB�`B�`B�`B�ZB�`B�ZB�TB�TB�HB�;B�BB�5B�BB�BB�;B�;B�BB�HB�TB�NB�NB�NB�HB�;B�;B�5B�/B�B�B�B��B�B��B�B�B�
B��B��B�B��B�
B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��BǮBǮBƨBȴBƨBƨBƨB��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBȴBȴBȴBǮBǮBǮBƨBƨBĜB��B��B�wB�qB�jB�dB�XB�RB�FB�FB�FB�RB�XB�XB�RB�RB�RB�RB�XB�LB�FB�3B�-B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B	�gB	�dB
=B
4B
J�B
p�B
��B
�`B
��B=BR�B�B�
B�AB�`B��B��BBBBBB
0BBB	-B
2B
0B%B%B	+BBBBB	-B
2B
2BBBBBBBBBBB`BhBvB^BbBAB6B<B	+B
0B	-BABTB^BdBhBjBVBABHBCB<B
2B
2BBBBBBB�B��B��B��B��B��B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B�B�B�zB�B�xB�vB�pB�oB�pB�YB�RB�XB�XB�eB�B�B��B��B��B��B��B�B�|B�hB�sB�SB�YB�IB�KB�\B�^B�kB�dB�lB�dB�^B�dB�YB�dB�dB�B�wB�wB�zB�tB�wB�tB�oB�wB�oB�lB�jB�_B�PB�WB�IB�XB�XB�SB�NB�WB�]B�iB�fB�cB�fB�\B�QB�QB�JB�HB�0B�+B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B� B� B�B��B�B�B�B�B��B��B�B�B��B��B��B��B��B��B��B̺B��B̾B̻B̻B��B��B��B��B�
B�B�B�B�B�B�B�B��B��B��B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B̻B̻BʯBǚBƕBćBÄB�{B�tB�gB�cB�ZB�XB�UB�dB�iB�gB�eB�eB�eB�eB�hB�]B�VB�DB�=B�-B�#B�-B�$B�B�B� B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0002 (+/- 2e-05) , vertically averaged dS =0.0059273 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219022022012717040020220127170400  IF  ARFMCODA035h                                                                20200828144322                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144403  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144403  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121902  IP  PSAL            A&ffD�L�G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170400  IP  PSAL            A&ffD�L�G�O�                