CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS   ^   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       06-May-2016 13:52:41Zcreation      
references        (http://www.argodatamgt.org/Documentation   comment       bThis netCDF file is generated using BODC's argoReader and netCDF writer software (argo@bodc.ac.uk)     user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    <@   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    <P   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    <T   REFERENCE_DATE_TIME                	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    <X   DATE_CREATION                  	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    <h   DATE_UPDATE                	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    <x   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    <�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  <�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  <�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  =   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        =@   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    =D   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    =H   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     =L   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    =l   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    =p   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     =t   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     =�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     =�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    =�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   axis      T      
_FillValue        A.�~       
resolution        X           =�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    =�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        X           =�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   	valid_min         �V�        	valid_max         @V�        axis      Y      
_FillValue        @�i�            =�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    	valid_min         �f�        	valid_max         @f�        axis      X      
_FillValue        @�i�            =�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    =�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    >    VERTICAL_SAMPLING_SCHEME                   	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    >   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ?   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ?   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        axis      Z      
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     x  ?   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  @�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  B   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  `  C�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  `  C�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  `  D@   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     x  D�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  F   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  G�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  `  I   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  `  Ih   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  `  I�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     x  J(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  K�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     x  M   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  0  N�   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    N�   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    Q�   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    T�   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  ,  W�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    W�   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    W�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    X   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    X   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                    X,   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  8  Y,   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    Yd   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  @  Yt   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        Y�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        Y�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        Y�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  @  Y�Argo profile    3.1 1.2 19500101000000  20170202062757  20170202062757  6900603 Argo UK                                                         Jon Turton                                                      PRES            TEMP            PSAL               �A   BO  30511                           2C  D   PROVOR                          OIN 07UK-S3-01                  V4.65                           841 @���8��1   @���8��@O"���m�>��E���1   ARGOS   Primary sampling: averaged                                                                                                                                                                                                                                      ����A   A   A   ?�  @�  Ap  A�  B  B8  B`  B�  B�  B�  B�  B�  B�  B�  C  C  C  C&  C0  C:  CD  CU  Cn  C�  C�� C�� C�� C�  C�  C΀ Cۀ C�  C�  D @ D  D  D@ D@ D� D%� D,  D2@ D8� D>� DE  DK@ DQ� DX  D^  Dd@ Dk  Dp� Dw  D}@ D�� D�� D�  D�  D�@ D�� D�� D�� D�� D�� D�  D�  D�` D�` D�� D�� D�  D�� D�  D�@ D�@ D�` Dƀ D�� D�� D�� D�  D�@ D�@ D܀ Dߠ D�� D�� D�  D�  D�  D�@ D�` D�� D�` B�7B�7B�%BdZBZBP�BP�BN�BN�BN�BM�BJ�BF�BD�B@�B<jB7LB0!B+B&�B#�B"�B�B��B  BJBJBB��B�B�B�sB�sB�B�B�B��B�B�ZB�/B�B�B�HB�5BɺBɺB�BɺB�-B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�qB��B�ZB  BB��B��B�ZB�B��B��BɺBƨBĜBÖBÖBBÖBB��B��B�}B�wB�jB�dB�jB�jB�qB�qB�}B��@߅@ߍP@�K�@���@�`B@�&�@� �@��y@�ȴ@ư!@Ƨ�@���@���@�9X@�;d@��@���@�K�@��T@�Ĝ@���@���@���@�hs@��/@�n�@���@�33@��@���@�dZ@�Ĝ@��@�  @�I�@��
@��m@�p�@��D@�O�@�~�@�bN@��@�n�@��@�;d@��@�Ĝ@���@~��@zM�@v�y@u��@so@p��@o��@ol�@n�y@n��@n��@n��@o
=@p�u@t�@tz�@|�/@��w@��F@�5?@�ff@{dZ@t�D@q��@o��@mV@j��@i��@g�@f�R@e?}@dZ@b�H@a�^@_�@]�T@\9X@Z=q@X�u@W�w@V�R@U`B@Tz�@S�
@S��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ?�  @�  Ap  A�  B  B8  B`  B�  B�  B�  B�  B�  B�  B�  C  C  C  C&  C0  C:  CD  CU  Cn  C�  C�� C�� C�� C�  C�  C΀ Cۀ C�  C�  D @ D  D  D@ D@ D� D%� D,  D2@ D8� D>� DE  DK@ DQ� DX  D^  Dd@ Dk  Dp� Dw  D}@ D�� D�� D�  D�  D�@ D�� D�� D�� D�� D�� D�  D�  D�` D�` D�� D�� D�  D�� D�  D�@ D�@ D�` Dƀ D�� D�� D�� D�  D�@ D�@ D܀ Dߠ D�� D�� D�  D�  D�  D�@ D�` D�� D�` B�7B�7B�%BdZBZBP�BP�BN�BN�BN�BM�BJ�BF�BD�B@�B<jB7LB0!B+B&�B#�B"�B�B��B  BJBJBB��B�B�B�sB�sB�B�B�B��B�B�ZB�/B�B�B�HB�5BɺBɺB�BɺB�-B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�qB��B�ZB  BB��B��B�ZB�B��B��BɺBƨBĜBÖBÖBBÖBB��B��B�}B�wB�jB�dB�jB�jB�qB�qB�}B��@߅@ߍP@�K�@���@�`B@�&�@� �@��y@�ȴ@ư!@Ƨ�@���@���@�9X@�;d@��@���@�K�@��T@�Ĝ@���@���@���@�hs@��/@�n�@���@�33@��@���@�dZ@�Ĝ@��@�  @�I�@��
@��m@�p�@��D@�O�@�~�@�bN@��@�n�@��@�;d@��@�Ĝ@���@~��@zM�@v�y@u��@so@p��@o��@ol�@n�y@n��@n��@n��@o
=@p�u@t�@tz�@|�/@��w@��F@�5?@�ff@{dZ@t�D@q��@o��@mV@j��@i��@g�@f�R@e?}@dZ@b�H@a�^@_�@]�T@\9X@Z=q@X�u@W�w@V�R@U`B@Tz�@S�
@S��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            PSAL            TEMP            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                       TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                             ds=0                                                                                                                                                                                                                                                           none                                                                                                                                                                                                                                                            Calibration error is manufacturer specified accuracy. Float that auto-corrects pressure. The quoted error is 2.4 dBar.                                                                                                                                          OW(2008) analysis. Mapping scales 8/4 and 4/2. OW chosen theta levels. Reference data CTD_for_DMQC_2010V2. Small drifts < 0.01 in salinity magnitude suggested. Conclusion, no correction.                                                                      The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   201108011022272011080110250720110801102507  BO  BO  BO  BO  ARGQARGQARSQARSQRTQCSCUTnullOW  2.0 2.0 null1.2                                                                                                                                                                                                                                                                 20110505191202201105161054552011080110222720110801102507QCP$QCP$IP  IP                                                                  ?�  ?�  G�O�G�O�D�` D�` G�O�G�O�F�^ F�^ F�^ F�^ 1688849860362238131072                                          