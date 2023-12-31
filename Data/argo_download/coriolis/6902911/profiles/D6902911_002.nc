CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2020-08-29T01:28:45Z creation; 2022-06-28T15:34:54Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_050f      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        h  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  =D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        h  =�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  @H   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     h  @�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  CL   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  E�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  FP   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  IT   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  K�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  N$   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  N�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Q(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     h  Q�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ]�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ]�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ]�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ]�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ]�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ]�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ]�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ]�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ]�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ^    HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ^   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ^   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  T,   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    T\   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    W\   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Z\   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ]\Argo profile    3.1 1.2 19500101000000  20200829012845  20230130142503  6902911 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18FR002                  5900A04                         844 @�x&So�1   @�x'��@Q�����'�)��n1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 25 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.1 dbar]                                                      @��@@  @�33@�  @�  @�33A   A  A   A1��A>ffAT��Ac33AnffA|��A�33A���A���A�  A�  A���A���A���A�33A�ffA�33A�33A�33A���A�  B ffBffB  B33B��BffBffB  B ffB#33B'��B,  B/33B4  B8��B<ffB?��BC33BF��B\ffB���B�33B�  B�ffBҙ�B�ffB���C�3CL�C33C%33C/ffC9��CC33CM33CW�CaffCkL�CuL�C� C�� C���C���C���C��fC���C�ٚC�ٚC��fC��3C��3C���C��fC�� C��3C��fC��fC�ٚC޳3C� C��C� C�3C��3D L�Ds3D� D��D  Ds3D%� D,�D2,�D8l�D>��DE  DKL�DQs3DW�fD]�3Dd,�Dj� Dp��Dw  D}L�D���D��3D��fD�fD�33D�i�D�vfD���D�� D���D���D�  D�9�D�P D�|�D�� D���D�ٚD���D��D�C3D�Y�D�p Dɩ�D̶fD���D�3D�  D�<�D�c3D߀ D��D�� D���D���D��D�C3D�ffD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@@  @�33@�  @�  @�33A   A  A   A1��A>ffAT��Ac33AnffA|��A�33A���A���A�  A�  A���A���A���A�33A�ffA�33A�33A�33A���A�  B ffBffB  B33B��BffBffB  B ffB#33B'��B,  B/33B4  B8��B<ffB?��BC33BF��B\ffB���B�33B�  B�ffBҙ�B�ffB���C�3CL�C33C%33C/ffC9��CC33CM33CW�CaffCkL�CuL�C� C�� C���C���C���C��fC���C�ٚC�ٚC��fC��3C��3C���C��fC�� C��3C��fC��fC�ٚC޳3C� C��C� C�3C��3D L�Ds3D� D��D  Ds3D%� D,�D2,�D8l�D>��DE  DKL�DQs3DW�fD]�3Dd,�Dj� Dp��Dw  D}L�D���D��3D��fD�fD�33D�i�D�vfD���D�� D���D���D�  D�9�D�P D�|�D�� D���D�ٚD���D��D�C3D�Y�D�p Dɩ�D̶fD���D�3D�  D�<�D�c3D߀ D��D�� D���D���D��D�C3D�ffD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@Ӿw@ԛ�@��@�Ĝ@�I�@��@ϥ�@��
@�5?@�{@�C�@�9X@�{@�p�@�A�@�z�@�(�@���@�9X@�1'@�O�@�?}@��@�/@�7L@�33@^{@A�#@�?���?�T>j~�>Z�>�  >�?}?�?'+?1��?DZ?f��?[dZ?f��?�A�?��?��/?��7?��?��?���?�(�@M�@�@�@�@	��@O�@��?�o?��#?�A�?�?�ȴ?���?� �?��h?�x�?�ff?�ff?���?�5??["�?2�?Q�?��?��>�ȴ>�`B>�/>��>��`>�J>��>��h>��D>� �>��#>�=q>���>��/>�K�>��?G�?�\? �>�->��>���>�;d>��;>^5?>)��=�`B=���=m�h<��t��C��q����j���O߾)��I�^�bMӾw�پ��9��zᾞ5?��ff���׾�Q쾾�۾Ƨ��V�׍P������xվ�&����� A�������	�^��D���-��j��+�X���vɿ!%�#S��$�/�'l��)xտ,I��.���0 ſ1녿3�Ͽ5?}�6�+�7�P1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @Ӿw@ԛ�@��@�Ĝ@�I�@��@ϥ�@��
@�5?@�{@�C�@�9X@�{@�p�@�A�@�z�@�(�@���@�9X@�1'@�O�@�?}@��@�/@�7L@�33@^{@A�#@�?���?�T>j~�>Z�>�  >�?}?�?'+?1��?DZ?f��?[dZ?f��?�A�?��?��/?��7?��?��?���?�(�@M�@�@�@�@	��@O�@��?�o?��#?�A�?�?�ȴ?���?� �?��h?�x�?�ff?�ff?���?�5??["�?2�?Q�?��?��>�ȴ>�`B>�/>��>��`>�J>��>��h>��D>� �>��#>�=q>���>��/>�K�>��?G�?�\? �>�->��>���>�;d>��;>^5?>)��=�`B=���=m�h<��t��C��q����j���O߾)��I�^�bMӾw�پ��9��zᾞ5?��ff���׾�Q쾾�۾Ƨ��V�׍P������xվ�&����� A�������	�^��D���-��j��+�X���vɿ!%�#S��$�/�'l��)xտ,I��.���0 ſ1녿3�Ͽ5?}�6�+�7�P1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB+BZBŢB�B%B(�B8RBF�B[#By�B�qB�TB�B��BbB(�BF�B[#BgmB�B��B�B�wB��B	M�B	�9B	q�B	]/B	x�B	H�B	6FB	G�B	k�B	x�B	�B	�3B	�/B	�B
  B
&�B
,B
5?B
W
B
t�B
�bB
��B
��B
��B
ŢB
��BM�B^5Bv�B}�B|�B�+BhsBgmBjBl�BiyBdZBcTBcTBcTBcTBdZBe`BaHBQ�B?}B2-B/B.B.B+B+B+B+B,B)�B'�B(�B,B7LB7LB49B;dBM�BS�BR�BW
B[#B]/B_;BaHBcTB`BBl�B`BBiyBbNB`BBaHBcTBhsBk�Br�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B+BZBŢB�B%B(�B8RBF�B[#By�B�qB�TB�B��BbB(�BF�B[#BgmB�B��B�B�wB��B	M�B	�9B	q�B	]/B	x�B	H�B	6FB	G�B	k�B	x�B	�B	�3B	�/B	�B
  B
&�B
,B
5?B
W
B
t�B
�bB
��B
��B
��B
ŢB
��BM�B^5Bv�B}�B|�B�+BhsBgmBjBl�BiyBdZBcTBcTBcTBcTBdZBe`BaHBQ�B?}B2-B/B.B.B+B+B+B+B,B)�B'�B(�B,B7LB7LB49B;dBM�BS�BR�BW
B[#B]/B_;BaHBcTB`BBl�B`BBiyBbNB`BBaHBcTBhsBk�Br�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                                                                    202301301425032023013014250320230130142503  IF  ARFMCODA035h                                                                20200829012845                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829012916                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829012916                      G�O�G�O�G�O�                IF  ARFMCODA050c                                                                20220608152657                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220608152752                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220608152752                      G�O�G�O�G�O�                IF  ARFMCODA050d                                                                20220618152744                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220618152925                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220618152925                      G�O�G�O�G�O�                IF  ARFMCODA050f                                                                20220628153333                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220628153454  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC5.8                                                                 20220628153454  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230130142503  IP  PSAL            @��D���G�O�                