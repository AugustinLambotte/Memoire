CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   r   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2020-08-29T01:34:27Z creation; 2022-06-26T15:27:47Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  <�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  =   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  >�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  ?T   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  A   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  B�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  CX   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  E    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  E�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  G\   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  I$   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  I�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  K`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  K�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    V�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    V�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    W    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    W   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  W   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    WH   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    WX   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    W\   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Wl   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Wp   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Wt   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Wx   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  M�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    M�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    P�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    S�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  V�Argo profile    3.1 1.2 19500101000000  20200829013427  20230130143206  6902912 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-18FR003                  5900A04                         844 @�{h8�91   @�{h8�9@Q�v���2�$��8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 25 dbar average from 1500 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                      A�  A�33A���A���A���A���A�  A�  A�  A���A���A���A�  A�  A�A�ffA�33B33B  B  BffB��B33B  B   B$  B(ffB+33B/33B4  B733B;��B@ffBC��BF��B\ffB���B�  B�  B���B�33B�33B���CffC� C��C%�3C/ffC933CC�CM33CW��Ca� Ck�CuffC�fC��3C���C��3C��fC��3C��fC�ٚC��fC��fC�� C��3C���C��3C�� C���Cϙ�CԳ3C���C���C��fC�� C��fC�3C�� D S3D� D� D��D@ Ds3D%� D+�3D2,�D8` D>��DD��DKS3DQ�fDW�fD^fDdL�Dj� Dp� Dv��D}` D���D��3D���D�fD�<�D�` D�� D���D�� D�� D���D�3D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�  A�33A���A���A���A���A�  A�  A�  A���A���A���A�  A�  A�A�ffA�33B33B  B  BffB��B33B  B   B$  B(ffB+33B/33B4  B733B;��B@ffBC��BF��B\ffB���B�  B�  B���B�33B�33B���CffC� C��C%�3C/ffC933CC�CM33CW��Ca� Ck�CuffC�fC��3C���C��3C��fC��3C��fC�ٚC��fC��fC�� C��3C���C��3C�� C���Cϙ�CԳ3C���C���C��fC�� C��fC�3C�� D S3D� D� D��D@ Ds3D%� D+�3D2,�D8` D>��DD��DKS3DQ�fDW�fD^fDdL�Dj� Dp� Dv��D}` D���D��3D���D�fD�<�D�` D�� D���D�� D�� D���D�3D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�ƨ@���@���@F�+@9�7@U�@_\)@Nv�@/�@��@r�@r�?�Ĝ?��\?��?�1?�-?��?öF?�=q?�V?߾w?���?���?�bN?��?�Q�?�5??�=q?���?��;?��?Լj?�V?�\@E�@33@1'@\)@@�;?�~�?�j?�1'?���?�?��?��?���?���?�=q?��+?�r�?��?�5??�1'?���?�5??�dZ?���?���?�?��?�ff?}�-?��F?;d?y��?u?}?hr�?Y�?L��??�w?(�9?�?X?�F?I�?$�>>�>�Q�>��->�7L>gl�>=p�>�=�S�=��-=Y�<T�����
�P�`���-��l��n��,1�B�\�Y��n�����7���9��zᾠA������F���۾�=q�׍P��S���׿�7��/111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111  @���@�ƨ@���@���@F�+@9�7@U�@_\)@Nv�@/�@��@r�@r�?�Ĝ?��\?��?�1?�-?��?öF?�=q?�V?߾w?���?���?�bN?��?�Q�?�5?G�O�G�O�?��;?��?Լj?�V?�\@E�@33@1'@\)@@�;?�~�?�j?�1'?���?�?��?��?���?���?�=q?��+?�r�?��?�5??�1'?���?�5??�dZ?���?���?�?��?�ff?}�-?��F?;d?y��?u?}?hr�?Y�?L��??�w?(�9?�?X?�F?I�?$�>>�>�Q�>��->�7L>gl�>=p�>�=�S�=��-=Y�<T�����
�P�`���-��l��n��,1�B�\�Y��n�����7���9��zᾠA������F���۾�=q�׍P��S���׿�7��/111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB ĜB%B�B�B%B(�B��BaHB+BBL�B�B��B	�B	�B	�B	/B	jB	o�B	�}B	�;B
B
�B
33B
9XB
O�B
T�B
jB
n�B
�B
ZB
�1B
��B
�!B
�?B
�5BN�B��B�'B�qBĜBB�XB�9B�FB�B��B��B�B��B��B��B��B�B��B��B��B��B��B��B��B��B�B�B�B�B�'B�'B�'B�'B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111  B ĜB%B�B�B%B(�B��BaHB+BBL�B�B��B	�B	�B	�B	/B	jB	o�B	�}B	�;B
B
�B
33B
9XB
O�B
T�B
jB
n�G�O�G�O�B
�1B
��B
�!B
�?B
�5BN�B��B�'B�qBĜBB�XB�9B�FB�B��B��B�B��B��B��B��B�B��B��B��B��B��B��B��B��B�B�B�B�B�'B�'B�'B�'B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                                                                    202301301432062023013014320620230130143206  IF  ARFMCODA035h                                                                20200829013427                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829013504                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829013504                      G�O�G�O�G�O�                IF  ARFMCODA049c                                                                20220527153145                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220527153245                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220527153245                      G�O�G�O�G�O�                IF  ARFMCODA050a                                                                20220606152902                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220606153000                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220606153000                      G�O�G�O�G�O�                IF  ARFMCODA050d                                                                20220616152710                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220616152816                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220616152816                      G�O�G�O�G�O�                IF  ARFMCODA050f                                                                20220626152648                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220626152747  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC5.8                                                                 20220626152747  QCF$                G�O�G�O�G�O�0000000000004000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230130143206  IP  PSAL            A�  D��fG�O�                