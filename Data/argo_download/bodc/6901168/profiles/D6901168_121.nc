CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS   9   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       10-Mar-2023 03:22:00Zcreation      
references        (http://www.argodatamgt.org/Documentation   comment       bThis netCDF file is generated using BODC's argoReader and netCDF writer software (argo@bodc.ac.uk)     user_manual_version       3.4    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    <H   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    <X   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    <\   REFERENCE_DATE_TIME                	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    <`   DATE_CREATION                  	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    <p   DATE_UPDATE                	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    <�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    <�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  <�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  <�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  =   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        =H   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    =L   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    =P   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     =T   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    =t   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    =x   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     =|   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     =�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     =�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    =�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   axis      T      
_FillValue        A.�~       
resolution        >�E�vQ�        =�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    =�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�E�vQ�        =�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   	valid_min         �V�        	valid_max         @V�        axis      Y      
_FillValue        @�i�            =�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    	valid_min         �f�        	valid_max         @f�        axis      X      
_FillValue        @�i�            =�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    >   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    >   VERTICAL_SAMPLING_SCHEME                   	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    >   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ?   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ?   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        axis      Z      
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������      �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  @   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  @�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  <  A�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  <  B   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  <  BD   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������      �  B�   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  Cd   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  DH   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  <  E,   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  <  Eh   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  <  E�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������      �  E�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  F�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���      �  G�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  �  H�   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                 	   I   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                 	   R   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                 	   [   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  �  d   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    d�   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    d�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    d�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    d�   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  d�   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  T  f|   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    f�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  f�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        gH   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        g`   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        gx   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  `  g�Argo profile    3.1 1.2 19500101000000  20230310032206  20230310032206  6901168 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               yA   BO  97091                           2C  D   APEX                            6610                            121012                          846 @�0��_1�1   @�0��_1�@N��G�{�C7�͞��1   ARGOS   Primary sampling: discrete                                                                                                                                                                                                                                      ����A   F   A   @���A!��Aq��A���A�33A�ffB��B��BG33Bo��B�33B�33B�33B���B���B���C  C  C33C   C)�3C4�C>  CH33C[�Co�fC�L�C�  C�L�C��3C�33C�  C��C�  D	ffD��D"��D.�3D;�3DHfDT� D`��DmffDy�3D�6fD��fD��fD��D�6fD�|�D�� D�fD��3D���Dԃ3D���D�ffB5?B33B49B49B49B33B-B!�B{�BǮBÖB�B��B�B�7B��B��B�jB�FB��B��B�}B�XB��BB�3B�3B�qB�LB�^B�RB�?B��B�B�B��B��B�DB��B��B��B�hB�hB�B�VB�PB�hB�=B�JB�DB��B��B�3B�dB�PB�XB�@��H@��@���@�@�@�
=@�o@�
=@׮@���@�x�@ϥ�@�/@�O�@�?}@��F@�M�@�O�@���@��/@�n�@��m@���@��m@��@�@��;@�ff@�33@�$�@�dZ@���@�ȴ@�n�@��@�@{��@v��@q�@o|�@m�@lj@jM�@i&�@fE�@fȴ@d�/@cdZ@c@b�!@ax�@e��@jJ@h�@^ff@U�-@O�111111111111111111111111111111111111111111111111111111111   333333333333333333333333333333333333333333333333333333333   111111111111111111111111111111111111111111111111111111111   @�  A$��At��A�fgA���A�  B��B��BH  BpfgB���B���B���B�33B�33B�33C33C33CffC 33C)�fC4L�C>33CHffC[L�Cp�C�fgC��C�fgC��C�L�C��C�34C��D	s3D��D"��D/  D;� DH3DT��D`��Dms3Dz  D�<�D���D���D�3D�<�D��3D��fD��D���D�3Dԉ�D�3D�l�B?B<�B=�B=�B=�B<�B6�B+�B��B�pB�XB��B��B��B��BڧB١B�,B�B��B�KB�?B�B�EB�QB��B��B�3B�B� B�B�B��B��B��B��B��B�B�CB�OB�CB�*B�*B��B�B�B�+B� B�B�B��B��B��B�'B�B�B��@��H@��@���@�@�@�
=@�o@�
=@׮@���@�x�@ϥ�@�/@�O�@�?}@��F@�M�@�O�@���@��/@�n�@��m@���@��m@��@�@��;@�ff@�33@�$�@�dZ@���@�ȴ@�n�@��@�@{��@v��@q�@o|�@m�@lj@jM�@i&�@fE�@fȴ@d�/@cdZ@c@b�!@ax�@e��@jJ@h�@^ff@U�-@O�111111111111111111111111111111111111111111111111111111111   333333333333333333333333333333333333333333333333333333333   111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                   N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.2                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             e_time assumes 0.09dBar/s ascent rate, alpha=0.0267, tau=18.6s                                                                                                                                                                                                  N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=-0.0096321                                                                                                                                                                                                                                                  Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                    N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 2.5/0.8 LAT 1.5/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 1/5. MAP_P_DELTA 100                                                                                                                                                   202303090956542023030911043220230309095654202303090956542023030911043220230309110132202303090956542023030911043220230309110432  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARSQARSQARSQRTQCRTSPSCUTnullnullOW  2.0 1.0 2.0 nullnull0.1                                                                                                                                                                                                                                                                                                                                                                                                 201710261306222017102613062320190314120254202303090956542023030911013220230309110432QCP$CV  QCP$IP  IP  IP                                                                                                  @���G�O�@���G�O�G�O�G�O�D�ffG�O�D�ffG�O�G�O�G�O�G���G���G���G���G���G���6389758                         131072                                                          