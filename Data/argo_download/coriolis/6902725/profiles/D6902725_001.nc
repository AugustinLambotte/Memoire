CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   l   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       [2017-11-16T08:31:47Z creation; 2020-01-21T15:41:25Z last update (coriolis COCQ (V3.2) tool)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8(   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  80   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8p   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     9   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     94   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9T   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9t   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9x   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >��	4E�   
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  <h   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  <�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  >�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  >�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  @�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  BP   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  B�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  Dl   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  D�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  F�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  H8   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  H�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  l  JT   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  J�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    U�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    U�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    U�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    U�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  U�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    V   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    V,   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    V0   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         V@   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         VD   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        VH   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    VL   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  Lp   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    L�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    O�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    R�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  U�Argo profile    3.1 1.2 19500101000000  20171116083147  20210930141421  6902725 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-16FR308                  5900A04                         844 @�&��6�1   @�&�ax9�@R�~u	O���p�X1   GPS     A   A   B   Primary sampling: averaged [10 sec sampling, 25 dbar average from 1000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 6.7 dbar]                                                      @�33@���A  AffA.ffAA��AS33A`  Aq��A���A���A�  A���A�33A�33A���A���A���A�  A���A�  Aݙ�A�ffA�33A�33A�ffB��B��BffB  B  B��B33B ffB$ffB'33B+33B/��B3��B7��B;33B?��BC��BH  B\ffB�33B���B�33B���B�33B�  B�ffCL�CL�C�C%ffC/� C9� CC33CM33CW33Ca  Ck  Cu  CffC���C�� C�� C��fC��fC�� C�ٚC���C�ٚC��3C��fC��fC��3C���C���Cϳ3CԳ3Cٳ3C޳3C���C�� C���C�ٚC�� D @ D�fD��D  D,�DffD%��D+�3D233D8` D>��DE&fDK,�DQ` DW��D]�3Dd9�Dj� Dp��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�33@���A  AffA.ffAA��AS33A`  Aq��A���A���A�  A���A�33A�33A���A���A���A�  A���A�  Aݙ�A�ffA�33A�33A�ffB��B��BffB  B  B��B33B ffB$ffB'33B+33B/��B3��B7��B;33B?��BC��BH  B\ffB�33B���B�33B���B�33B�  B�ffCL�CL�C�C%ffC/� C9� CC33CM33CW33Ca  Ck  Cu  CffC���C�� C�� C��fC��fC�� C�ٚC���C�ٚC��3C��fC��fC��3C���C���Cϳ3CԳ3Cٳ3C޳3C���C�� C���C�ٚC�� D @ D�fD��D  D,�DffD%��D+�3D233D8` D>��DE&fDK,�DQ` DW��D]�3Dd9�Dj� Dp��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@o�;@vff@�A�@��@���@��D@�ȴ@�V@vv�@D1@��@�?���?�5?? A�>�u>�?��=�P�}󶾠A���&�r�!��|��?}���>S��?XQ�?^v�?d�/?r�?|(�?��?���?�Z?��y?��?�Q�?�
=?ɺ^?�O�?θR?�=q?�n�?�7@6�R@m�-@\j@_�w@]/@X�@SC�@J��@D(�@<��@5�@)X@$�@{@@�
?�{?�;d?޸R?Ͳ-?�"�?��?�{?�Ĝ?���?xQ�?Fff?.�?=�?=p�?"��>�^5>��?&$�?R�?N�?E`B?49X?!%?{?C�>�&�>�"�>��>�Q�>��u>���>���>~��>["�=�j=�7L<������
�H�9��9X��`B�V�!���@��T���bMӾ���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@o�;@vff@�A�@��@���@��D@�ȴ@�V@vv�@D1@��@�?���?�5?? A�>�u>�?��=�P�}󶾠A���&�r�!��|��?}���>S��?XQ�?^v�?d�/?r�?|(�?��?���?�Z?��y?��?�Q�?�
=?ɺ^?�O�?θR?�=q?�n�?�7@6�R@m�-@\j@_�w@]/@X�@SC�@J��@D(�@<��@5�@)X@$�@{@@�
?�{?�;d?޸R?Ͳ-?�"�?��?�{?�Ĝ?���?xQ�?Fff?.�?=�?=p�?"��>�^5>��?&$�?R�?N�?E`B?49X?!%?{?C�>�&�>�"�>��>�Q�>��u>���>���>~��>["�=�j=�7L<������
�H�9��9X��`B�V�!���@��T���bMӾ���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oA��wB 1'B �'B)�B�FB�B;dB\)BuB�XB�B�B&�BL�BA�B�B�B�;B��BhB+B?}BiyBu�B�JB��B�B	�B	W
B	dZB	u�B	}�B	�DB	��B	�B	�B	�B	ŢB
"�B
-B
1'B
C�B
S�B
cTB
�bB�B�B�B33B9XB;dB8RB1'B+B#�B�BPB	7BB��B�mB�mB�ZB�)B�
B��B��BƨB��B��B�^B�B��B�B�B��B��B��B�3B�}B��B��B�dB�LB�3B�FB�'B�!B�B�B�B�-B�9B�9B�3B�B�B�B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111A��wB 1'B �'B)�B�FB�B;dB\)BuB�XB�B�B&�BL�BA�B�B�B�;B��BhB+B?}BiyBu�B�JB��B�B	�B	W
B	dZB	u�B	}�B	�DB	��G�O�G�O�G�O�G�O�B
"�B
-B
1'B
C�B
S�B
cTB
�bB�B�B�B33B9XB;dB8RB1'B+B#�B�BPB	7BB��B�mB�mB�ZB�)B�
B��B��BƨB��B��B�^B�B��B�B�B��B��B��B�3B�}B��B��B�dB�LB�3B�FB�'B�!B�B�B�B�-B�9B�9B�3B�B�B�B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                                                    202004090832202021093014142120210930141421  IF  ARFMCODA016c                                                                20171116083147                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171116083157  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.1                                                                 20171116083157  QCF$                G�O�G�O�G�O�0000000000004000IF      SCOO0.48                                                                20200121152216  CF  PSAL            B$ffB'33?�                  IF      SCOO0.48                                                                20200121152216  CF  TEMP            B+33B/��@�                  IF      COCQ3.2                                                                 20200121154125                      G�O�G�O�G�O�                IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200409083220  IP  PSAL            @�33Dp��G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20210930141421  IP  PSAL            @�33Dp��G�O�                