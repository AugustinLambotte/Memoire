CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   v   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2019-06-20T11:53:05Z creation; 2019-06-20T11:54:14Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8,   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8<   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8L   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8T   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    90   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    94   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     98   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9X   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9x   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  <�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  =,   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  ?   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  ?|   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  AT   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  C,   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  C�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  E|   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  E�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  G�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  I�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  J   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  x  K�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Ll   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    W�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    W�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    W�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    W�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  W�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    W�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    X    HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    X   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         X   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         X   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        X   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    X    	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ND   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Nt   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Qt   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Tt   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  WtArgo profile    3.1 1.2 19500101000000  20190620115305  20200409090956  6902730 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               pA   IF                                  2C  D   ARVOR                           AI2600-16FR313                  5900A04                         844 @�}��>��1   @�}觡�@Q33333�3+�Z��1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 25 dbar average from 1300 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.0 dbar]                                                      @ff@Fff@�  @�33@�  @�  A��AffA   A333A@  ANffAc33AnffA|��A�ffA�  A�33A���A�ffA�ffA�33A�33A�  A�  A�  A�  A���A�  A�  B33BffB  B33B33B33B��B��B ffB#33B'��B,  B/��B4ffB733B;33B?��BC��BHffB^  B�  B�  B�33B���Bҙ�B�ffB�  C� C33C33C%�C/��C9� CC33CM�CW33Ca33Ck�CuffCffC��3C��3C���C���C���C��3C��fC���C��3C�� C���C���C���CŦfCʙ�C���C�ٚC�� C޳3C���C�3C��fC�3C��3D ,�Dl�D��D�DffDl�D%��D+�fD2&fD8` D>� DD��DK,�DQl�DW�fD^  Dd@ Dj�3Dp�3Dw  D}FfD���D�ɚD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@Fff@�  @�33@�  @�  A��AffA   A333A@  ANffAc33AnffA|��A�ffA�  A�33A���A�ffA�ffA�33A�33A�  A�  A�  A�  A���A�  A�  B33BffB  B33B33B33B��B��B ffB#33B'��B,  B/��B4ffB733B;33B?��BC��BHffB^  B�  B�  B�33B���Bҙ�B�ffB�  C� C33C33C%�C/��C9� CC33CM�CW33Ca33Ck�CuffCffC��3C��3C���C���C���C��3C��fC���C��3C�� C���C���C���CŦfCʙ�C���C�ٚC�� C޳3C���C�3C��fC�3C��3D ,�Dl�D��D�DffDl�D%��D+�fD2&fD8` D>� DD��DK,�DQl�DW�fD^  Dd@ Dj�3Dp�3Dw  D}FfD���D�ɚD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�(�@���@�o@��H@�t�@˥�@˾w@�  @�A�@̃@͡�@���@�M�@Η�@���@ΰ!@ΰ!@Χ�@��@���@��
@���@��@�1@q%@T��@@A�@:M�@@t�@Ĝ?���?�ff?��?޸R?���?���?˅?ȓu?Ƨ�?��7?��?�J?��H?�+?��h?�(�?�z�?��w?�1?���?�?}?z�?o�?b�\?L�D?@Ĝ?0bN?(�9?'+?%��?&ff?'+?(1'?+�?+�?,��?.{?.V?&�y?�R?�?��?V?	7L?��>��m>�r�>��>���>�E�>���>�Z>�(�>��`>�$�>u>_;d>Xb>E��>8Q�>&�y=�S�=�7L=#�
<o�#�
�#�
�����bN�,1�=p��L�;S�Ͼe`B�w�پ����C�����n�������MӾ�l������ȴ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@�(�@���@�o@��H@�t�@˥�@˾w@�  @�A�@̃@͡�@���@�M�@Η�@���@ΰ!@ΰ!@Χ�@��@���@��
@���@��@�1@q%@T��@@A�@:M�@@t�@Ĝ?���?�ff?��?޸R?���?���?˅?ȓu?Ƨ�?��7?��?�J?��H?�+?��h?�(�?�z�?��w?�1?���?�?}?z�?o�?b�\?L�D?@Ĝ?0bN?(�9?'+?%��?&ff?'+?(1'?+�?+�?,��?.{?.V?&�y?�R?�?��?V?	7L?��>��m>�r�>��>���>�E�>���>�Z>�(�>��`>�$�>u>_;d>Xb>E��>8Q�>&�y=�S�=�7L=#�
<o�#�
�#�
�����bN�,1�=p��L�;S�Ͼe`B�w�پ����C�����n�������MӾ�l������ȴ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBBBDB49B��B� B��B��B��B�B�LBÖB�;B�B��B  B	7B
=BDBhB�B(�Bz�B�{B�'B'�B��B�-B�B	�B	�XB	�B
�B
l�B
v�B
z�B
�+B
�VB
��B
��B
��B
��B
�'B
�'B
��B
�#B
�TB
��B
��BuB%�B8RBG�BT�BhsBm�Bq�Br�Br�Br�Bv�Bz�B{�B{�B}�B� B�B�+B�VB�hB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BBBDB49B��B� B��B��B��B�B�LBÖB�;B�B��B  B	7B
=BDBhB�B(�Bz�B�{B�'B'�B��B�-B�B	�B	�XB	�B
�B
l�B
v�B
z�B
�+B
�VB
��B
��B
��B
��B
�'B
�'B
��B
�#B
�TB
��B
��BuB%�B8RBG�BT�BhsBm�Bq�Br�Br�Br�Bv�Bz�B{�B{�B}�B� B�B�+B�VB�hB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 1.1,  -CTD2018V01 & ARGO2018V01 -                                                                                                                                    202004090909562020040909095620200409090956  IF  ARFMCODA029d                                                                20190620115305                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115414  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115414  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200409090956  IP  PSAL            @ffD���G�O�                