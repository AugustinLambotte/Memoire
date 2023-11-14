CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS   9   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       10-Mar-2023 03:16:07Zcreation      
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
_FillValue                     d�   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                     d�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                     d�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                     d�   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                    e   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  p  g   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                     g�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  �  g�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�         h,   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�         hL   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�         hl   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  �  h�Argo profile    3.1 1.2 19500101000000  20230310031611  20230310031611  6901168 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               lA   BO  92837                           2C  D   APEX                            6610                            121012                          846 @�[6;�1   @�[6;�@N�q�i��A(H��1   ARGOS   Primary sampling: discrete                                                                                                                                                                                                                                      ����A   F   B   @���A$��Ap  A�33Ař�A�  B��B"  BHffBn��B�ffB���B���B���B�  B�33C� C�fC�fC 33C*� C4  C=�fCH�C\ffCp  C�  C�  C��C�  C���C�ٚC�ٚC��fD	s3D�fD"�fD/fD;�fDH�DT� D`�fDm�fDz�D�33D�\�D�� D�3D�C3D���D���D�	�D�� D�	�D�y�D�fD�y�BA�B)�B
�HB
%B
8RB	�sB
�%B	�{B
��B
��Be`BR�Bu�B�Bx�BiyB�{B�bBv�B|�By�B�{B��B��B�DB�hB��B�VB��B��B�\B��B��B��B��B�{B�bB�JB�7B�DB�DB�%Bx�B�B�+B�1B�=B�Bv�B|�B�hB��B�B��B�?B�jB�L@��#@�^@�h@��@�@��;@�l�@ɡ�@��-@�t�@��w@���@���@���@��@��/@�@��^@�%@�@�(�@�%@��h@�@�Ĝ@�ff@��-@��`@��@�A�@{S�@vff@sC�@p  @m�T@l�/@k33@iX@hb@f��@fE�@d�j@c"�@`r�@aG�@a�#@a�^@aG�@[dZ@Z�H@`r�@c��@ko@iX@fȴ@a��@]V111111111111111111111111111111111111111111111111111111111   344433443333333333333333333333333333333333333333333333333   111411111111111111111111111111111111111111111111111111111   @�fgA)��At��A���A�  A�ffB  B#33BI��Bp  B�  B�fgB�34B�fgBݙ�B���C��C33C33C � C*��C4L�C>33CHfgC\�3CpL�C�&fC�&fC�33C�&fC��3C�  C�  C��D	�fD��D"��D/�D;��DH,�DT�3D`��Dm��Dz,�D�<�D�fgD�ɚD��D�L�D��4D��gD�4D���D�4Dԃ4D� D�4BJ�G�O�G�O�G�O�B
A�B	�G�O�G�O�B
�BBn�B\%B~�B�QB�Br�B��B��B�B� B�B��B�B��B�vB��B��B��B��B��B��B��B��B��B��B��B��B�}B�jB�wB�wB�XB�	B�@B�_B�eB�qB�:B�B�"B��B��B�PB�B�uBŠB��@��#@�^@�hG�O�@�@��;@�l�@ɡ�@��-@�t�@��w@���@���@���@��@��/@�@��^@�%@�@�(�@�%@��h@�@�Ĝ@�ff@��-@��`@��@�A�@{S�@vff@sC�@p  @m�T@l�/@k33@iX@hb@f��@fE�@d�j@c"�@`r�@aG�@a�#@a�^@aG�@[dZ@Z�H@`r�@c��@ko@iX@fȴ@a��@]V111111111111111111111111111111111111111111111111111111111   344433443333333333333333333333333333333333333333333333333   111411111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
G�O�G�O�G�O�<#�
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
;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                   N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.3                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             e_time assumes 0.09dBar/s ascent rate, alpha=0.0267, tau=18.6s                                                                                                                                                                                                  N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=-0.0091391                                                                                                                                                                                                                                                  Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                    N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 2.5/0.8 LAT 1.5/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 1/5. MAP_P_DELTA 100                                                                                                                                                   202303090956542023030911043220230309095654202303090956542023030911043220230309110132202303090956542023030911043220230309110432  BO  BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQARSQRTSPSCUTSCUTSCUTSCUTnullnullOW  1.0 2.0 2.0 2.0 2.0 nullnull0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 2017061813132820190314120254201903141202542019031412025420190314120254202303090956542023030911013220230309110432CV  QCF$QCF$QCF$QCF$IP  IP  IP                  TEMP            TEMP            TEMP            TEMP                                                            G�O�E� E� Ap  A$��G�O�G�O�G�O�G�O�E�H E�H Ap  A$��G�O�G�O�G�O�G�R�G�R�G�R�G�R�G�R�G�R�G�R�G�R�                131072          131072          131072          131072                                                          