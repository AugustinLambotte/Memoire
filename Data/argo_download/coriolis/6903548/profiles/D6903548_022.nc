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
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20190620115650  20210607154239  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @؁6ffff1   @؁7�݀@Tf��O0x@,�+c�C�1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 9.1 dbar]                                                      AffAffA.ffA>ffANffAa��AnffA|��A�33A�ffA�ffA�33A�33A�ffA�  A�33A�33A���A�  A�ffA���A���A���B   BffB  BffB��B  B��B33B��B$ffB(  B,  B0  B4  B8  B;33B?33BC��BG��BK��BO��BT  BX  B\ffB_33Bc33Bg��Bk��Bp  Bt  BxffB{33B33B���B���B���B���B�  B�  B�33B���B���B���B���B���B���B�  B�ffB���B���B���B�  B�  B�33B���B���B�  B���B���B�  B���B���B�  B���B���B�  B�ffB���B�33B�33B�33B�33B�33Bڙ�Bޙ�B���B�  B�33BB���B�33B���B�33CL�C��CffC33C	� C�3C� CffCL�C�C  CffC�3C� CffC33C!�C#ffC%��C'ffC)33C+�C-ffC/�3C1� C3L�C5�C7ffC9��C;ffC=33C?� CA��CCffCE33CG� CI�3CK� CML�CO�CQL�CS� CUL�CW� CY�3C[� C]L�C_�CaffCc�3Ce� CgL�Ci�Ck� Cm��Co�3Cq��Cs� CuffCwL�CyL�C{33C}�C  C��3C��fC��fC��fC��fC��fC�ٚC�ٚC�ٚC�ٚC�ٚC��fC��fC��3C��3C�� C�� C���C���C���C���C��fC��3C��3C�� C���C�ٚC��fC��fC��3C�� C�� C���C���C��3C�� C���C�ٚC��fC�� C�� C���C���C��fC��3C�� C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��3C�� C���C��fC�� C���C��fC��3C�� C�Có3Cĳ3C�� C�� C�� C���CɦfCʦfC˳3C�� C���CΙ�CϦfCг3C�� C�ٚCӳ3C�� C�ٚCֳ3Cי�CئfC�� Cڙ�C۳3C�� CݦfC޳3C���C�fC�� C���C�fC�� C���C�3C��C�fC�� C�fC�� C�ٚC�� CC�� C�ٚC�� C�fC��C��fC�� C��fC���C��3C���C�ffC�  D FfDl�D�3D��D&fDs3D�fD��D
9�Dy�D� D�D9�D` D�3D��DFfDs3D� D�3DFfDy�D��D  DS3D�fD �3D"fD#@ D$s3D%��D&�fD(33D)y�D*��D,fD-9�D.s3D/�3D0��D29�D3��D4� D6  D7@ D8l�D9�3D;�D<@ D=l�D>��D@�DA9�DBffDC�3DEfDF33DG� DH��DJ  DK33DLffDM��DO�DP@ DQy�DR��DS��DU,�DVl�DW�fDY  DZ33D[l�D\��D]��D_,�D`y�Da�fDc  Dd33De�fDf� Dg�3Di33Djs3Dk� Dl��Dn9�Do�fDp�3Dr  DsL�Dty�Du��Dv� Dx9�Dy��Dz�fD|fD}FfD~� D�fD��3D�&fD�ɚD�\�D�� D��3D�33D��fD�|�D�  D��fD�` D�3D��fD�C3D���D�y�D�3D��3D�VfD���D�� D�FfD��3D�� D�  D�� D�c3D��fD���D�C3D�ٚD�s3D��D��fD�c3D�  D�� D�<�D���D�� D�  D��3D�c3D�3D��fD�<�D��3D�y�D�  D��fD�c3D�3D��3D�C3D��fD�y�D�3D���D�` D���D��fD�@ D�� D�|�D��D���D�` D�  D��3D�9�D���D��3D��D��fD�S3D�� D���D�L�D��D���D�)�D���D�S3D���D�� D�C3D��fD�|�D�3D��3D�` D���D�� D�FfD�� D�|�D��D���D�\�D�  D��3D�FfD���D�s3D�fD¹�D�` D�3Dę�D�33D���D�|�D��Dǹ�D�Y�D���Dɓ3D�9�D���D�vfD�  D��fD�c3D�  DΜ�D�6fD��3D�s3D�3Dѳ3D�VfD���DӜ�D�C3D��fDՌ�D�  Dֳ3D�Y�D���Dؙ�D�<�D��3Dڀ D��D۳3D�\�D�fDݦfD�<�D��3D�vfD��D๚D�VfD�  D� D�@ D��fD�|�D��D�3D�S3D��3D�3D�33D��fD�y�D�  D��fD�\�D��3D�fD�6fD��fD�y�D��D�ɚD�c3D���D�3D�<�D���D�s3D�fD���D�\�D���D���D�<�D�� D�y�D��D�� D�` D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  AffAffA.ffA>ffANffAa��AnffA|��A�33A�ffA�ffA�33A�33A�ffA�  A�33A�33A���A�  A�ffA���A���A���B   BffB  BffB��B  B��B33B��B$ffB(  B,  B0  B4  B8  B;33B?33BC��BG��BK��BO��BT  BX  B\ffB_33Bc33Bg��Bk��Bp  Bt  BxffB{33B33B���B���B���B���B�  B�  B�33B���B���B���B���B���B���B�  B�ffB���B���B���B�  B�  B�33B���B���B�  B���B���B�  B���B���B�  B���B���B�  B�ffB���B�33B�33B�33B�33B�33Bڙ�Bޙ�B���B�  B�33BB���B�33B���B�33CL�C��CffC33C	� C�3C� CffCL�C�C  CffC�3C� CffC33C!�C#ffC%��C'ffC)33C+�C-ffC/�3C1� C3L�C5�C7ffC9��C;ffC=33C?� CA��CCffCE33CG� CI�3CK� CML�CO�CQL�CS� CUL�CW� CY�3C[� C]L�C_�CaffCc�3Ce� CgL�Ci�Ck� Cm��Co�3Cq��Cs� CuffCwL�CyL�C{33C}�C  C��3C��fC��fC��fC��fC��fC�ٚC�ٚC�ٚC�ٚC�ٚC��fC��fC��3C��3C�� C�� C���C���C���C���C��fC��3C��3C�� C���C�ٚC��fC��fC��3C�� C�� C���C���C��3C�� C���C�ٚC��fC�� C�� C���C���C��fC��3C�� C���C�ٚC��fC��3C�� C���C���C��fC��3C�� C�ٚC��3C�� C���C��fC�� C���C��fC��3C�� C�Có3Cĳ3C�� C�� C�� C���CɦfCʦfC˳3C�� C���CΙ�CϦfCг3C�� C�ٚCӳ3C�� C�ٚCֳ3Cי�CئfC�� Cڙ�C۳3C�� CݦfC޳3C���C�fC�� C���C�fC�� C���C�3C��C�fC�� C�fC�� C�ٚC�� CC�� C�ٚC�� C�fC��C��fC�� C��fC���C��3C���C�ffC�  D FfDl�D�3D��D&fDs3D�fD��D
9�Dy�D� D�D9�D` D�3D��DFfDs3D� D�3DFfDy�D��D  DS3D�fD �3D"fD#@ D$s3D%��D&�fD(33D)y�D*��D,fD-9�D.s3D/�3D0��D29�D3��D4� D6  D7@ D8l�D9�3D;�D<@ D=l�D>��D@�DA9�DBffDC�3DEfDF33DG� DH��DJ  DK33DLffDM��DO�DP@ DQy�DR��DS��DU,�DVl�DW�fDY  DZ33D[l�D\��D]��D_,�D`y�Da�fDc  Dd33De�fDf� Dg�3Di33Djs3Dk� Dl��Dn9�Do�fDp�3Dr  DsL�Dty�Du��Dv� Dx9�Dy��Dz�fD|fD}FfD~� D�fD��3D�&fD�ɚD�\�D�� D��3D�33D��fD�|�D�  D��fD�` D�3D��fD�C3D���D�y�D�3D��3D�VfD���D�� D�FfD��3D�� D�  D�� D�c3D��fD���D�C3D�ٚD�s3D��D��fD�c3D�  D�� D�<�D���D�� D�  D��3D�c3D�3D��fD�<�D��3D�y�D�  D��fD�c3D�3D��3D�C3D��fD�y�D�3D���D�` D���D��fD�@ D�� D�|�D��D���D�` D�  D��3D�9�D���D��3D��D��fD�S3D�� D���D�L�D��D���D�)�D���D�S3D���D�� D�C3D��fD�|�D�3D��3D�` D���D�� D�FfD�� D�|�D��D���D�\�D�  D��3D�FfD���D�s3D�fD¹�D�` D�3Dę�D�33D���D�|�D��Dǹ�D�Y�D���Dɓ3D�9�D���D�vfD�  D��fD�c3D�  DΜ�D�6fD��3D�s3D�3Dѳ3D�VfD���DӜ�D�C3D��fDՌ�D�  Dֳ3D�Y�D���Dؙ�D�<�D��3Dڀ D��D۳3D�\�D�fDݦfD�<�D��3D�vfD��D๚D�VfD�  D� D�@ D��fD�|�D��D�3D�S3D��3D�3D�33D��fD�y�D�  D��fD�\�D��3D�fD�6fD��fD�y�D��D�ɚD�c3D���D�3D�<�D���D�s3D�fD���D�\�D���D���D�<�D�� D�y�D��D�� D�` D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?�J@�@Ĝ@+��@1�^@:�@D�@Hr�@I�#@J-@K��@Pb@N{@G��@G\)@Hb@H��@G;d@E�@G��@G�@Fff@A%@-�h@V@�?��?׮?�j@�/@�@��@%�@+�
@0b@<1@q�@r��@sdZ@s�F@t�@}@}p�@zn�@t1@r^5@o+@j�@` �@^��@]@^��@_�@^�+@\(�@Yhs@Y%@PA�@O;d@V5?@Q&�@Q7L@J�!@I�@D��@E�T@Gl�@H�`@K��@OK�@R�@X �@\�@]?}@[33@_�@g\)@f��@f5?@e�@e`B@e��@e��@e�@e�-@d�/@co@` �@_�@[��@X��@W;d@^@b�\@a��@`��@`A�@_;d@]p�@[�m@Z-@Y&�@Z=q@Y��@Y��@X�9@T�@Q��@P�@N@G��@F5?@B��@>{@>��@Ax�@@�u@?�@>$�@<I�@9�^@9��@8�`@8Q�@8bN@8�9@9&�@7�@6��@2�!@-`B@+�F@,(�@,1@%O�@!7L@ 1'@
=@�@ A�@!G�@ A�@�;@�P@\)@�w@ r�@�@�+@�@9X@1@ƨ@@�F@�m@�@t�@C�@�@5?@�+@�R@v�@{@�@��@��@�@?}@�@@��@�T@��@�T@�T@�@��@�-@p�@�/@"�@-@�#@��@x�@x�@x�@hs@X@7L@�@��@�`@Ĝ@�R@��@7L@�@��@�@G�@7L@G�@Ĝ@1'@�@ �@A�@Ĝ@%@x�@�@b@ȴ@ff@?}@V@z�@z�@��@�j@�/@�@��@�+@��@��@��@�y@5?@
=@
=@�@�R@ff@E�@5?@5?@5?@5?@E�@E�@V@V@v�@ff@5?@�h@O�@��@(�@9X@I�@Z@z�@�D@�j@��@V@�@/@V@��@��@�@�j@�@9X@1@�
@��@t�@C�@
�!@
�@	��@	hs@	�@bN@�@�y@ff@�@�j@1@"�@M�@7L@ �u?��?�{?���?�r�?�ȴ?���?�M�?�p�?�x�?�A�?�E�?�5??�~�?Ƨ�?Õ�?��?��?��
?��h?�Q�?�t�?��w?�V?�7L?�ff?��?��F?�Ĝ?�Q�?��?r�!?j~�?b��?Z��?R�?L�D?I�^?G�?A%?=/?9�?/\)?$�/?;d?b?�F?��?r�>�|�>�?}>�/>�7L>�^5>�~�>��w>�7L>���>���>p��>W
=>O�;>2->!��>�+>	7L>   =�j=���=��=q��=49X=C�<�<�/<���<���<D��;o��o���P�`��C����㽟�w���T��vɽ������`��`B���پ$ݾ\)��w�-V�5?}�C���L�;V�`A��dZ�e`B�gl��n���w�پ�����������1'���^��C���ƨ��O߾�bN��t���
=������R��G����T��V�������j��Q쾻dZ����  �\�Õ�����Ƨ��7L�����bN��t�����b�����"Ѿ�/��5?�߾w�������
��ff��l����þ����V�����&��F����KǾ�Q���#���H��p���|� �� Ĝ����o������˿����r��	xտ
~������D��Ϳ�h�{�������;�bN�hs�-�33��Ͽz��j�������E��Kǿ�P��ٿQ�������#�^5�"ѿ��(����푿p��5?��w� A��!%�!G��!G��!�7�"�\�#S��#���#�
�$Z�%��%�˿&$ݿ&�y�'+�'+�'l��'l��'(r��(�ÿ)�^�*~��*���+�+C��+��,1�,�D�,�Ϳ-V�-V�-V�-O߿-��-��.{�.{�.{�.V�.V�.��/�;�0 ſ0 ſ0�`�1hs�1녿2-�2�!�2�333�333�3�Ͽ4z�4�j�4�j�4�j�5?}�5��6�6E��6E��6�+�6ȴ�7�P�7�ٿ7�ٿ7�ٿ8b�8Q�8Q�8�u�8�u�8�u�8�u�8���9��9���:^5�:���:�H�:�H�:�H�:�H�;"ѿ;dZ�;dZ�;��;�m�<(��<��<푿=/�=p��=p��=p�111111111111111111111111111111111144111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ?�J@�@Ĝ@+��@1�^@:�@D�@Hr�@I�#@J-@K��@Pb@N{@G��@G\)@Hb@H��@G;d@E�@G��@G�@Fff@A%@-�h@V@�?��?׮?�j@�/@�@��@%�@+�
G�O�G�O�@q�@r��@sdZ@s�F@t�@}@}p�@zn�@t1@r^5@o+@j�@` �@^��@]@^��@_�@^�+@\(�@Yhs@Y%G�O�G�O�@V5?@Q&�@Q7L@J�!@I�@D��@E�T@Gl�@H�`@K��@OK�@R�@X �@\�@]?}@[33@_�@g\)@f��@f5?@e�@e`B@e��@e��@e�@e�-@d�/@co@` �@_�@[��@X��@W;d@^@b�\@a��@`��@`A�@_;d@]p�@[�m@Z-@Y&�@Z=q@Y��@Y��@X�9@T�@Q��@P�@N@G��@F5?@B��@>{@>��@Ax�@@�u@?�@>$�@<I�@9�^@9��@8�`@8Q�@8bN@8�9@9&�@7�@6��@2�!@-`B@+�F@,(�@,1@%O�@!7L@ 1'@
=@�@ A�@!G�@ A�@�;@�P@\)@�w@ r�@�@�+@�@9X@1@ƨ@@�F@�m@�@t�@C�@�@5?@�+@�R@v�@{@�@��@��@�@?}@�@@��@�T@��@�T@�T@�@��@�-@p�@�/@"�@-@�#@��@x�@x�@x�@hs@X@7L@�@��@�`@Ĝ@�R@��@7L@�@��@�@G�@7L@G�@Ĝ@1'@�@ �@A�@Ĝ@%@x�@�@b@ȴ@ff@?}@V@z�@z�@��@�j@�/@�@��@�+@��@��@��@�y@5?@
=@
=@�@�R@ff@E�@5?@5?@5?@5?@E�@E�@V@V@v�@ff@5?@�h@O�@��@(�@9X@I�@Z@z�@�D@�j@��@V@�@/@V@��@��@�@�j@�@9X@1@�
@��@t�@C�@
�!@
�@	��@	hs@	�@bN@�@�y@ff@�@�j@1@"�@M�@7L@ �u?��?�{?���?�r�?�ȴ?���?�M�?�p�?�x�?�A�?�E�?�5??�~�?Ƨ�?Õ�?��?��?��
?��h?�Q�?�t�?��w?�V?�7L?�ff?��?��F?�Ĝ?�Q�?��?r�!?j~�?b��?Z��?R�?L�D?I�^?G�?A%?=/?9�?/\)?$�/?;d?b?�F?��?r�>�|�>�?}>�/>�7L>�^5>�~�>��w>�7L>���>���>p��>W
=>O�;>2->!��>�+>	7L>   =�j=���=��=q��=49X=C�<�<�/<���<���<D��;o��o���P�`��C����㽟�w���T��vɽ������`��`B���پ$ݾ\)��w�-V�5?}�C���L�;V�`A��dZ�e`B�gl��n���w�پ�����������1'���^��C���ƨ��O߾�bN��t���
=������R��G����T��V�������j��Q쾻dZ����  �\�Õ�����Ƨ��7L�����bN��t�����b�����"Ѿ�/��5?�߾w�������
��ff��l����þ����V�����&��F����KǾ�Q���#���H��p���|� �� Ĝ����o������˿����r��	xտ
~������D��Ϳ�h�{�������;�bN�hs�-�33��Ͽz��j�������E��Kǿ�P��ٿQ�������#�^5�"ѿ��(����푿p��5?��w� A��!%�!G��!G��!�7�"�\�#S��#���#�
�$Z�%��%�˿&$ݿ&�y�'+�'+�'l��'l��'(r��(�ÿ)�^�*~��*���+�+C��+��,1�,�D�,�Ϳ-V�-V�-V�-O߿-��-��.{�.{�.{�.V�.V�.��/�;�0 ſ0 ſ0�`�1hs�1녿2-�2�!�2�333�333�3�Ͽ4z�4�j�4�j�4�j�5?}�5��6�6E��6E��6�+�6ȴ�7�P�7�ٿ7�ٿ7�ٿ8b�8Q�8Q�8�u�8�u�8�u�8�u�8���9��9���:^5�:���:�H�:�H�:�H�:�H�;"ѿ;dZ�;dZ�;��;�m�<(��<��<푿=/�=p��=p��=p�111111111111111111111111111111111144111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�By�B��BjB��BB�B+BVBhB�B#�B33B'�B&�B2-B+B2-B9XBA�BC�BC�BcTB	\B	��B	�/B
B
)�B
`BB
]/B
�B
�qB
�fB
�B
��B
ŢBVB��B��B�B�B��B��B�dB�-B�'B��B�B�-B��B��B��B��B��B��B��B��B��Bm�B��B��B��B��B��B��B��B��B��B�B�9BĜB��B��B�mB�`B�yBB��B1B��B%BB1B%B%BBB  B��B��B��BBJBuB{BuBuBoBhBVBJBJBJBhBbB\B+BBB��B��B��B�B�B�B�B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�ZB�NB�HB�;B�HB��B��B��B��B��B��B��B��B��B��B��B�
B�B��B��B��B��B��B��B��B�B�
B�B�B�B�)B�5B�/B�5B�/B�)B�/B�/B�)B�)B�/B�5B�/B�5B�;B�;B�5B�5B�;B�;B�5B�5B�/B�B�#B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBȴBȴBȴBȴBȴBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BȴBȴBȴB��B��B��B��BɺBĜBB�qB�jB�jB�dB�jB�jB�XB�FB�9B�9B�?B�?B�?B�9B�?B�3B�-B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111144111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BB}gB�mBnB�B�B�BB
�B�B�B2B'cB6�B+|B*uB5�B.�B5�B<�BEBG"BG"Bf�B	�B	�qB	�B
�B
-�B
c�B
`�B
��B
��B
��B
�$G�O�G�O�BY�B�QB�|B��B��B�B�B��B��B��B��B��B��B�QB�>B�]B�]B�iB�WB�B�uG�O�G�O�B�]B�cB�oB�2B�DB�,B�QB�WB�JB��B��B�(B�_B�xB��B��B�B�B�B�B�B	�B�B�B	�B	�B�B�B�B zB�[B�B�B�BBBBB�B�B�B�B�B�B�B�B�B
�B�B�B�sB�[B�UB�$B�B�B�$B�UB�OB�HB�<B�$B�$B�*B�$B�$B�*B�0B�6B�0B�0B��B��B��B��B��B�kB�eB�kB�kB�~B؊B؊BׄBׄB؊B؊BږBۜB�~BׄB�~B�~B�~B�~BׄBِBږBِBِBݩBߵB��B�B��B�BߵB�B�BߵBߵB�B��B�B��B��B��B��B��B��B��B��B��B�BݩBޯBݩBݩBܣBܣBݩBܣBܣBܣBۜBۜBۜBِB�xB�eB�eB�eB�eB�kB�qB�eB�YB�eB�_B�eB�eB�xB�kB�kB�eB�YB�SB�FB�FB�@B�@B�@B�@B�@B�@B�FB�FB�SB�SB�YB�YB�_B�YB�SB�eB�eB�eB�_B�eB�eB�eB�eB�eB�eB�eB�eB�kB�kB�kB�kB�_B�_B�_B�YB�YB�_B�_B�eB�_B�_B�eB�kB�kB�kB�kB�qB�qB�qB�qB�qB�qB�xB�qB�qB�qB�qB�qB�qB�kB�kB�kB�eB�_B�_B�_B�YB�YB�SB�SB�MB�MB�MB�SB�@B�@B�@B�MB�SB�SB�MB�FB�(B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�uB�iB�oB�oB�oB�iB�iB�iB�iB�iB�cB�]B�cB�cB�iB�iB�iB�cB�iB�WB�QB�JB�>B�,B�2B�8B�8B�2B�,B�2B�,B�2B�2B�2B�,B�,B�&B�&B�&B�&B�&B�&B�,B�2B�2B�2B�,B�,B�2B�,B�2B�2B�2B�2B�2B�2B�8B�>B�>B�>B�>B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�QB�QB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�QB�JB�QB�QB�QB�QB�WB�WB�QB�QB�QB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�uB�uB�|B�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B�|B��B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111144111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542392021060715423920210607154239  IF  ARFMCODA029d                                                                20190620115650                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115713  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115713  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154239  IP  PSAL            AffD��3G�O�                