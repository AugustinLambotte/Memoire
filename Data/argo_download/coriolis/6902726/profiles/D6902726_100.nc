CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   v   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2020-05-11T14:45:34Z creation; 2020-05-11T14:47:09Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_034d      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
_FillValue                  ,  WtArgo profile    3.1 1.2 19500101000000  20200511144534  20211004092454  6902726 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               dA   IF                                  2C  D   ARVOR                           AI2600-16FR309                  5900A04                         844 @�s�I��J1   @�s����@Q1�����3��g��1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 25 dbar average from 1300 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 9.0 dbar]                                                      AffA!��A0  AA��AP  A`  Aq��A�  A���A�ffA�ffA���A�  A�33A�  A�  A�  A�  Aٙ�A�33A�  A�  A�33A���B��BffBffB  B  B  B  B   B$  B(  B,  B/��B3��B7��B;��B?��BD  BG��B^  B���B�33B���B���B�33B�  B�  C��C��C� C%33C/33C9L�CCffCM33CWffCa�3CkL�Cu  CffC�� C��fC��fC��3C��fC���C���C��fC��fC�� C�ٚC���C�� CŦfCʙ�Cϙ�CԳ3Cٙ�CަfC�fC虚C��fC�ٚC�ٚD 9�Dl�D��D��Dl�Dy�D%��D+�3D2&fD8y�D>�fDD��DK  DQ` DW�fD]��Dd33Djl�Dp�fDw3D}s3D���D��D�fD��D�6fD�VfD�vfD��fD�� D�&f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  AffA!��A0  AA��AP  A`  Aq��A�  A���A�ffA�ffA���A�  A�33A�  A�  A�  A�  Aٙ�A�33A�  A�  A�33A���B��BffBffB  B  B  B  B   B$  B(  B,  B/��B3��B7��B;��B?��BD  BG��B^  B���B�33B���B���B�33B�  B�  C��C��C� C%33C/33C9L�CCffCM33CWffCa�3CkL�Cu  CffC�� C��fC��fC��3C��fC���C���C��fC��fC�� C�ٚC���C�� CŦfCʙ�Cϙ�CԳ3Cٙ�CަfC�fC虚C��fC�ٚC�ٚD 9�Dl�D��D��Dl�Dy�D%��D+�3D2&fD8y�D>�fDD��DK  DQ` DW�fD]��Dd33Djl�Dp�fDw3D}s3D���D��D�fD��D�6fD�VfD�vfD��fD�� D�&f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�ƨ@��@��@�bN@��^@jM�@ZM�@O�w@J�H@D9X@,�@!�^@�h@��@M�@
�H@l�?��?��`?��;?�hs?��?��`?�~�?���?��y?�
=?�\)?�ȴ?�?}?��;?�\)?��?�\)?��w?�V?�1?�l�?�p�?�?��?�O�?���?]/?E�?<j?<�?4z�?K?CS�?>5??@  ?7��?:�H?6ȴ?2-?,�D?(��?,1?0�`?7��?0��?*=q?&ff?$��?!�7?/?��?�D?�9?�>�dZ>��>�l�>ڟ�>�hs>š�>�Q�>���>�>�"�>�n�>��9>}�>n��>^5?>Q�>7K�>�u=�=���=49X<���t��t�����ȴ9�   ��w�B�\�Z��r�!��o�������w���9X��vɾ�+��n���"Ѿ�Z�� ž�X��\�	x�1111111111111111111111111111111111111111111111111111111111111444111111111111111111111111111111111111111111111111111111  @�ƨ@��@��@�bN@��^@jM�@ZM�@O�w@J�H@D9X@,�@!�^@�h@��@M�@
�H@l�?��?��`?��;?�hs?��?��`?�~�?���?��y?�
=?�\)?�ȴ?�?}?��;?�\)?��?�\)?��w?�V?�1?�l�?�p�?�?��?�O�?���?]/?E�?<j?<�?4z�?K?CS�?>5??@  ?7��?:�H?6ȴ?2-?,�D?(��?,1?0�`?7��G�O�G�O�G�O�?$��?!�7?/?��?�D?�9?�>�dZ>��>�l�>ڟ�>�hs>š�>�Q�>���>�>�"�>�n�>��9>}�>n��>^5?>Q�>7K�>�u=�=���=49X<���t��t�����ȴ9�   ��w�B�\�Z��r�!��o�������w���9X��vɾ�+��n���"Ѿ�Z�� ž�X��\�	x�1111111111111111111111111111111111111111111111111111111111111444111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�Bk�B�BD�BZB��Bu�B��B��BN�B��BL�B�B��B	\B	=qB	S�B	�=B	�!B	�^B	ɺB
  B
1B
�B
L�B
n�B
o�B
�B
�B
��B
��B
��B
��B
��B
��B
��B
�B
�FB
�wB
�wB
B
��B
�;B
��B%B\B�B)�B=qBC�BE�BK�BM�BR�BS�BW
BXBYB]/BaHBffBgmB�B�B��BJB-BG�B1'BD�BQ�B<jBI�BE�BE�B6FB9XB49B7LB8RB5?B0!B%�B-Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bv�Bw�Bx�By�By�Bz�Bz�B{�B{�B{�B{�B|�B|�B|�B|�B}�B}�B}�B}�B}�B}�B}�B}�B~�B~�B~�B~�1111111111111111111111111111111111111111111111111111111111111444444444444444444444441111111111111111111111111111111111  B=2B�%BnBfvB|?BNB�oB��B�Bq�B�9Bp7B��B�PB	2�B	aB	w�B	��B	��B	�%B	�B
#�B
,B
B�B
p�B
��B
��B
�B
�B
ÿB
��B
��B
��B
ýB
��B
��B
�B
�KB
�}B
�{B
�B
��BHB�B*;B3sBB�BNBa�Bg�Bi�Bo�Bq�BwBxB{3B|6B}AB�[B�tB��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�"B�!B�!B�!1111111111111111111111111111111111111111111111111111111111111444444444444444444444441111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0009 (+/- 7e-06) , vertically averaged dS =0.035011 (+/- 0.01)                                                                                                                                                                                             No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              Salinity drift or offset detected - OWC fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                          202005120935392021100409245420211004092454  IF  ARFMCODA034d                                                                20200511144534                      G�O�G�O�G�O�                IF  ARGQCOQC4.5                                                                 20200511144709  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.5                                                                 20200511144709  QCF$                G�O�G�O�G�O�0000000000004000IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200512093539  IP  PSAL            AffD�&fG�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211004092454  IP  PSAL            AffD�&fG�O�                