CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       [2019-06-20T12:01:25Z creation; 2020-01-21T11:41:42Z last update (coriolis COCQ (V3.2) tool)    
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
resolution        =���   axis      Z        p  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  =P   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  =�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  @\   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  @�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  Ch   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  E�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  Ft   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  I�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  K�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  N`   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  N�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Ql   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  R   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ]�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ]�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ]�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ]�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ]�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ^$   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ^4   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ^8   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ^H   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ^L   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ^P   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ^T   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  Tx   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    T�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    W�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Z�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ]�Argo profile    3.1 1.2 19500101000000  20190620120125  20211006133609  6902729 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-16FR312                  5900A04                         844 @؏f���1   @؏g`T�@@T~T���@2�r2L�f1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 25 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.2 dbar]                                                      @��@Fff@�33@�  @�  @�33@���AffAffA.ffA@  AS33A`  Ak33A|��A�ffA�  A���A�  A���A�  A���A�33A�33A�  A���A�33A�  A�  A���B ffB  B  B��B33B  B��B33B ffB$ffB'33B+33B/��B3��B8  B:��B?33BC��BF��B\ffB�33B�ffB�ffB�  B�33B�33B�ffC33CL�C��C%��C/� C9L�CCffCM� CW33CaffCkffCuL�CL�C��3C��fC��3C���C��fC��3C���C�� C��fC��fC��3C���C�� C�ٚC�� Cϳ3C��fC��fC��fC�� C�fC�� C�3C���D ,�DY�D�fD�fD33D��D%�fD+�3D233D8s3D>�3DEfDK` DQy�DW��D^fDd9�Dj�fDp��Dw  D}@ D���D���D�	�D�#3D�@ D�vfD�y�D�� D���D��fD��3D�33D�<�D�VfD�s3D��3D�� D��3D���D�#3D�C3D�\�D�y�Dɜ�D̳3D�ٚD��3D��D�<�D�S3D�vfD� D��D���D� D��D�0 D�VfD�y�D�� D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@Fff@�33@�  @�  @�33@���AffAffA.ffA@  AS33A`  Ak33A|��A�ffA�  A���A�  A���A�  A���A�33A�33A�  A���A�33A�  A�  A���B ffB  B  B��B33B  B��B33B ffB$ffB'33B+33B/��B3��B8  B:��B?33BC��BF��B\ffB�33B�ffB�ffB�  B�33B�33B�ffC33CL�C��C%��C/� C9L�CCffCM� CW33CaffCkffCuL�CL�C��3C��fC��3C���C��fC��3C���C�� C��fC��fC��3C���C�� C�ٚC�� Cϳ3C��fC��fC��fC�� C�fC�� C�3C���D ,�DY�D�fD�fD33D��D%�fD+�3D233D8s3D>�3DEfDK` DQy�DW��D^fDd9�Dj�fDp��Dw  D}@ D���D���D�	�D�#3D�@ D�vfD�y�D�� D���D��fD��3D�33D�<�D�VfD�s3D��3D�� D��3D���D�#3D�C3D�\�D�y�Dɜ�D̳3D�ٚD��3D��D�<�D�S3D�vfD� D��D���D� D��D�0 D�VfD�y�D�� D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����9������b���P��1'��1'���T��Mӿ�녿U?}�D�/�7Kǿ6E��#�
�vɿb�bN�`B�%�$ݿ A����j���+�J����v�;��
=�->�b?��?��?F��?�\)?�l�?��R?�"�?���?�9?�ff?�O�?��-@ Q�@b@
�@p�@|�@�7@��@�-@ �9@8r�@e?}@v�@o;d@m�@m@j�@h�`@jM�@pQ�@t�@p�@mV@iG�@e`B@_+@W�@R��@OK�@I�@BM�@@A�@<�j@:�!@9�@6��@5�-@3�
@,�@'�;@%p�@!�7@�@J@E�@�P@E�@��@�P@r�@
=@p�@��@�@��@9X@	%?��w?�+?��T?��?��?�M�?o�?E�?'+?>�/>�l�>l�D>49X=��=�%<t���P��E��\)�0 žO�;�u���^������Ĝ�����پ��7��ƨ�����r���-��^5�G�����	xտV� ſ����ٿ��(��vɿ A��"Mӿ#o�$�/�&��(1'�)��+��-O߿/��0�`�333�49X�5��6�+111111111111141111111144111111414411111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111���9������b���P��1'��1'���T��Mӿ�녿U?}�D�/�7Kǿ6E�G�O��vɿb�bN�`B�%�$ݿ A����jG�O�G�O���v�;��
=�->�b?��?��G�O�?�\)G�O�G�O�?�"�?���?�9?�ff?�O�?��-@ Q�@b@
�@p�@|�@�7G�O�G�O�@ �9@8r�@e?}@v�@o;d@m�@m@j�@h�`@jM�@pQ�@t�@p�@mV@iG�@e`B@_+@W�@R��@OK�@I�@BM�@@A�@<�j@:�!@9�@6��@5�-@3�
@,�@'�;@%p�@!�7@�@J@E�@�P@E�@��@�P@r�@
=@p�@��@�@��@9X@	%?��w?�+?��T?��?��?�M�?o�?E�?'+?>�/>�l�>l�D>49X=��=�%<t���P��E��\)�0 žO�;�u���^������Ĝ�����پ��7��ƨ�����r���-��^5�G�����	xտV� ſ����ٿ��(��vɿ A��"Mӿ#o�$�/�&��(1'�)��+��-O߿/��0�`�333�49X�5��6�+111111111111141111111144111111414411111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;oG�O�;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�XB�RB�RB�^B�RB�LB�qB��B�BB,B=qBP�BR�BO�Bo�Bw�Bp�B�B�hB�bB�uB�'B�B�?B�sB�B�B	1B	9XB	]/B	L�B	��B	�RB	��B	�B
  B
PB
%�B
D�B
A�B
5?B
K�B
aHB
_;B
l�B
{�B
k�B
��B
��B
��B�B��B�#B�HB�fB�B�B��B	7B�B�B�B�B�B\B%BB��B��B�B�B�B�B�B�B�B�B�B�ZB�`B�BB�5B�B��B�B�B��B��B��B��B��B��B��B��B��B��BɺB��B�wB�XB�LB�3B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111141111111144111111414411111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B�XB�RB�RB�^B�RB�LB�qB��B�BB,B=qBP�BR�G�O�Bo�Bw�Bp�B�B�hB�bB�uB�'G�O�G�O�B�sB�B�B	1B	9XB	]/G�O�B	��G�O�G�O�B	�B
  B
PB
%�B
D�B
A�B
5?B
K�B
aHB
_;B
l�B
{�G�O�G�O�B
��B
��B�B��B�#B�HB�fB�B�B��B	7B�B�B�B�B�B\B%BB��B��B�B�B�B�B�B�B�B�B�B�ZB�`B�BB�5B�B��B�B�B��B��B��B��B��B��B��B��B��B��BɺB��B�wB�XB�LB�3B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111141111111144111111414411111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�<#�
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
G�O�<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                                                    202004090850442021100613360920211006133609  IF  ARFMCODA029d                                                                20190620120125                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620120219  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620120219  QCF$                G�O�G�O�G�O�0000000000004000IF      SCOO0.48                                                                20200121112917  CF  TEMP            B ffB ff?�                  IF      SCOO0.48                                                                20200121112917  CF  TEMP            Ak33Ak33?�                  IF      SCOO0.48                                                                20200121112917  CF  PSAL            B?33BC��?�                  IF      SCOO0.48                                                                20200121112917  CF  PSAL            B ffB ff?�                  IF      SCOO0.48                                                                20200121112917  CF  TEMP            B?33BC��?�                  IF      SCOO0.48                                                                20200121112917  CF  PSAL            Ak33Ak33?�                  IF      COCQ3.2                                                                 20200121114142                      G�O�G�O�G�O�                IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200409085044  IP  PSAL            @��D��3G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211006133609  IP  PSAL            @��D��3G�O�                