CDF   
   
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   V   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       [2019-06-20T12:01:21Z creation; 2020-01-21T11:41:39Z last update (coriolis COCQ (V3.2) tool)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8    DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    80   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8@   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8P   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8X   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    94   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    98   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     9<   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9\   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9|   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
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
resolution        =���   axis      Z        X  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  <8   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        X  <�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  =�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     X  >@   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  ?�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  @�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  AH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  B�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  B�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  DP   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  E�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  F    PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  GX   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  G�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    Rd   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    Rh   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    Rl   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    Rp   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  Rt   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    R�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    R�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    R�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         R�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         R�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        R�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    R�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  I   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    I8   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    L8   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    O8   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  R8Argo profile    3.1 1.2 19500101000000  20190620120121  20211006133609  6902729 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               ;A   IF                                  2C  D   ARVOR                           AI2600-16FR312                  5900A04                         844 @�Yf����1   @�Yf����@T[?�W�z@/��p��|8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 25 dbar average from 500 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                       A!��A+33AA��AP  Ac33Aq��A~ffA�  A�  A�  A�33A�33A���A�  A�33A�33A�ffA�ffA�ffA�ffA�ffA�33B   B33B��BffB33B��B  BffB   B#��B'��B+��B/��B3��B8  B<  B@��BD  BG33B[33B�33B���B�  B�33B�ffB���B���CL�C� C� C%ffC/�3C9ffCC33CM� CWL�Ca�Ck� Cu�3C��C�� C��3C���C���C��fC���C��3C���C���C��3C���C�� C�ٚC�� Cʙ�Cϳ3C�ٚC�� C�� C�� C�fC��fC��fC�L�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A!��A+33AA��AP  Ac33Aq��A~ffA�  A�  A�  A�33A�33A���A�  A�33A�33A�ffA�ffA�ffA�ffA�ffA�33B   B33B��BffB33B��B  BffB   B#��B'��B+��B/��B3��B8  B<  B@��BD  BG33B[33B�33B���B�  B�33B�ffB���B���CL�C� C� C%ffC/�3C9ffCC33CM� CWL�Ca�Ck� Cu�3C��C�� C��3C���C���C��fC���C��3C���C���C��3C���C�� C�ٚC�� Cʙ�Cϳ3C�ٚC�� C�� C�� C�fC��fC��fC�L�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����/��p��Ͳ-��X��$ݿ����T�����T��z῾�ۿ�j��1'��(���Z���˿�+��?}��ff���h��;d��J��F�j~�������MӾ�Q�333��`B<�j���T���>L��>�~�?	7L?Fff?S��?M�h?WK�?�x�?��7?��y?st�?��+?�v�?Ł?��?���@��@�;@�@��@�@
^5@t�@C�@	G�@	G�@Ĝ@�@�@V@��@�@1'@o@	G�@�-@V@o@��@��@ A�?��m?�Q�?���?� �?�b?��?�hs?�o?��?�G�?�+?�~�?��911111111111111111111144141411114414411444111111111111111111111111111111111111111111111  ��/��p��Ͳ-��X��$ݿ����T�����T��z῾�ۿ�j��1'��(���Z���˿�+��?}��ff���h��;dG�O�G�O��j~�G�O���M�G�O��333��`B<�j���TG�O�G�O�>�~�G�O�G�O�?S��?M�hG�O�G�O�G�O�?��y?st�?��+?�v�?Ł?��?���@��@�;@�@��@�@
^5@t�@C�@	G�@	G�@Ĝ@�@�@V@��@�@1'@o@	G�@�-@V@o@��@��@ A�?��m?�Q�?���?� �?�b?��?�hs?�o?��?�G�?�+?�~�?��911111111111111111111144141411114414411444111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;oG�O�;oG�O�;o;o;o;oG�O�G�O�;oG�O�G�O�;o;oG�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	L�B	M�B	O�B	P�B	S�B	VB	VB	YB	ZB	aHB	aHB	�B	��B	�9B	��B	ÖB	ŢB	��B	��B	�NB	�B
!�B	�B
#�B
H�B
%�B
cTB
I�B
P�B
cTB
w�B
n�B
JB
��B
�B
�B
�3B
�wB
�NB
�B
�B
�)B
�#B
��B>wBO�BjB�VB�dBȴB��BƨB�wB�^B�wB�qB�XB�^B�^B��BBǮB��B��BǮB�}B�qBƨBĜB��B�RB�-B�'B�B��B�B��B��B��B��B��B��B�B�^BŢBȴ11111111111111111111144141411114414411444111111111111111111111111111111111111111111111  B	L�B	M�B	O�B	P�B	S�B	VB	VB	YB	ZB	aHB	aHB	�B	��B	�9B	��B	ÖB	ŢB	��B	��B	�NB	�G�O�G�O�B
#�G�O�B
%�G�O�B
I�B
P�B
cTB
w�G�O�G�O�B
��G�O�G�O�B
�3B
�wG�O�G�O�G�O�B
�)B
�#B
��B>wBO�BjB�VB�dBȴB��BƨB�wB�^B�wB�qB�XB�^B�^B��BBǮB��B��BǮB�}B�qBƨBĜB��B�RB�-B�'B�B��B�B��B��B��B��B��B��B�B�^BŢBȴ11111111111111111111144141411114414411444111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
G�O�<#�
G�O�<#�
<#�
<#�
<#�
G�O�G�O�<#�
G�O�G�O�<#�
<#�
G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                                                    202004090850042021100613360920211006133609  IF  ARFMCODA029d                                                                20190620120121                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620120151  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.2                                                                 20190620120151  QCF$                G�O�G�O�G�O�0000000000004000IF      SCOO0.48                                                                20200121112537  CF  TEMP            B��B��?�                  IF      SCOO0.48                                                                20200121112537  CF  PSAL            B33B33?�                  IF      SCOO0.48                                                                20200121112537  CF  PSAL            B��B��?�                  IF      SCOO0.48                                                                20200121112537  CF  TEMP            B33B33?�                  IF      COCQ3.2                                                                 20200121114139                      G�O�G�O�G�O�                IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200409085004  IP  PSAL            A!��C�L�G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211006133609  IP  PSAL            A!��C�L�G�O�                