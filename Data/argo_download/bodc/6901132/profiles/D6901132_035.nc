CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS   F   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       07-Aug-2021 07:37:36Zcreation      
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
resolution        ?�������       ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       @8   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       AP   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  Bh   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  B�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������       C@   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       DX   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       Ep   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  H  F�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  H  F�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                  H  G   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������       G`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       Hx   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���       I�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  J�   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    K   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    Q   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    W   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ]   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ]\   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ]x   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ]�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ]�   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ]�   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  _�   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    _�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  `   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        `|   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        `�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        `�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  `�Argo profile    3.1 1.2 19500101000000  20210807073739  20210807073739  6901132 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               #A   BO  124072                          2C  D   APEX                            6236                            120210                          846 @֣.���1   @֣.���@R+S����)Ƨ1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            A   A   A   @���AG�Ad  A�A�G�A�(�B��Bp�B*�RBF=qB[33Bl�B�  B�=qB�(�B��)B��RB�
=B��3B���B�8RB�C ^�C	��C�C�)C(��C2nC<�)CF޸CPs3CZ��Cd��Cn��CxQ�C��C�C���C�+�C�)C��C��C�
C��{C�fC�Y�C�aHC�C�C��C�B�C��C�1�C���D	�D�{D"{D.��DG��D`e�Dy��D�G
D��=D�>fD��=D�2=D��=D�33D�fD�5�D��B�?Bn�BخB OBUgBabB��B�]B�B	4B	l�B	��B
?B
�uB
��B
�sB
�B"BI�Bp�B�5B�|B��B��B�oB��BHBxB}"B��B��B�2B�wB�B�WB�QB��B�WB�tB�2B��B��B�TB��B�B��B�hB��B��B��B�2B�,B��B��B��B�TB��B��B�'B�pB��B�]B�B�B��B�\B��B�hB��B�N?hK^?���?*�L������&��ο��������󶿚#:�ra|�:C��Z��>��>��??�W?)��?D�f?��b?���?��?�?ІY?��~?�_?��*?���?S�
?Jں?��?W>�?�x?� �?���?�)_?G�?0��?/˒?A��?>p;?4Ɇ?%��?c�?�?
)�?�A>�H>�u�>��z>�E�>�,�>�*0>cS�=��=/O���.���¾He���:���<����Ͽ�Y�S����|�'1�S�9�~�?�DɆ1111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111  @�fgA�AjffA���A�z�A�\)B34B
>B,Q�BG�B\��Bn�B���B�
>B���B���B��B��
B�� B�B�B�Q�C �C

=CW
C�C)8RC2�zC=�CGECPٙCZ�RCe)Cn�3Cx�RC�:�C�8RC�.C�^�C�O\C�:�C�<)C�J=C��C�9�C���C��{C�w
C�H�C�u�C�NC�eC�.D	34D�D"D.��DGˆD`]Dy�HD�S�D��
D�K3D��
D�?
D��
D�@ D��3D�B�D�(�B��BwB��BBY{BelB�B��B��B	5�B	nB	�fB
?�B
��B
��B
�|B
�B@BJ�BsCB�WB��B��B�aB��B�B��B|�B�oB��B��B�B�B�CB�2B�B��B�LB��B��B�{B�kB��B�B�tB�LB��B�9B�[B��B��B��B�UB�6B��B��B�JB�B�|B��B�MB��B�\B�\B�;B��B�GB��B�2B��?hK^?���?*�L������&��ο��������󶿚#:�ra|�:C��Z��>��>��??�W?)��?D�f?��b?���?��?�?ІY?��~?�_?��*?���?S�
?Jں?��?W>�?�x?� �?���?�)_?G�?0��?/˒?A��?>p;?4Ɇ?%��?c�?�?
)�?�A>�H>�u�>��z>�E�>�,�>�*0>cS�=��=/O���.���¾He���:���<����Ͽ�Y�S����|�'1�S�9�~�?�DɆ1111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�
<#�
<��<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.4                                                                                                                                                                                                                                                        none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=-0.0041881                                                                                                                                                                                                                                                  Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 2/0.5.MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2021v01 and ARGO2020v03 ref. data.                                                                                                                    202108051535122021080614372520210805153512202108051535122021080614372520210806143725BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190828115032201908281150322019082811503620190828115041202006241321372021080515351220210806143725  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@���@���@���G�O�G�O�G�O�G�O�D��D��D��G�O�G�O�G�T G�T G�T G�T G�T G�T G�T                                 6389758         6389758         131072                                          