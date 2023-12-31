CDF       
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2017-08-20T00:05:40Z creation; 2018-06-23T14:07:31Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8X   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9    JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9(   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >��	4E�   
_FillValue        A.�~            9,   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           94   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9<   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9D   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9H   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9P   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9T   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9X   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9\   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :\   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        H  :`   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  <�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        H  =<   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  ?�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     H  @   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  B`   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  D�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  E<   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  G�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  H   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  J`   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  L�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  M<   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  O�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     H  P   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    [�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    [�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    [�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    [�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  [�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    \   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    \   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    \    HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         \0   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         \4   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        \8   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    \<   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  R`   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    R�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    U�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    X�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  [�  [�             ,  [�Argo profile    3.1 1.2 19500101000000  20170820000540  20180623140731  3901874 MOCCA-EU                                                        Romain Cancouet                                                 PRES            TEMP            PSAL               "A   IF                                  2C  D   ARVOR                           AR2600-16FR037                  5900A00                         844 @�(+<M`1   @�(.��@S��:�����E�l2   IRIDIUM A   A   A   Primary sampling: averaged [10 sec sampling, 4 dbar average from 2000 dbar to 1400 dbar; 10 sec sampling, 2 dbar average from 1400 dbar to 250 dbar; 10 sec sampling, 1 dbar average from 250 dbar to 2.5 dbar]                                                    @陚@���AffA!��A0  A@  AQ��A\��AnffA���A�33A�  A���A�33A�  A���A�33A�  Aə�A�  A�ffAݙ�A�  A�ffA���B   B��B��B33B��BffB  B��B!33B%33B(��B,��B0��B4��B8��B<��B@��BD��BH��BLffBPffBTffBXffB\��B`��BdffBh��Bl��Bp��Bu33By33B}33B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�33B�33B�33B�33B�ffB�ffB�ffB���B���B���B���B�  B�  B�33B�33B�33B�33B�ffB�ffB���B�B���B���B�  B�  B�  B�33B�33B�33B�ffB�ffBٙ�Bۙ�B���B���B���B�  B�  B�  B�33B�ffB�ffB�ffB�B���B���B���B�  B�  B�  C ��C��C�3C�3C��C�fC�fC  C	�C
�C�C33CL�CffCffC� C� C��C��C��C  C33C  C��C�3C�fC�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @陚@���AffA!��A0  A@  AQ��A\��AnffA���A�33A�  A���A�33A�  A���A�33A�  Aə�A�  A�ffAݙ�A�  A�ffA���B   B��B��B33B��BffB  B��B!33B%33B(��B,��B0��B4��B8��B<��B@��BD��BH��BLffBPffBTffBXffB\��B`��BdffBh��Bl��Bp��Bu33By33B}33B���B���B���B���B���B���B���B���B�  B�  B�  B�  B�33B�33B�33B�33B�ffB�ffB�ffB���B���B���B���B�  B�  B�33B�33B�33B�33B�ffB�ffB���B�B���B���B�  B�  B�  B�33B�33B�33B�ffB�ffBٙ�Bۙ�B���B���B���B�  B�  B�  B�33B�ffB�ffB�ffB�B���B���B���B�  B�  B�  C ��C��C�3C�3C��C�fC�fC  C	�C
�C�C33CL�CffCffC� C� C��C��C��C  C33C  C��C�3C�fC�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�Ĝ@�M�@��\@��@�33@���@�&�@���@�V@�bN@�7L@ǥ�@�ȴ@�?}@�o@�;d@�E�@��@�|�@�`B@�@�
=@�C�@�@� �@��@�5?@�hs@�dZ@�"�@�I�@�%@�1'@�G�@~�@~{@|�D@}?}@}@}`B@~{@}?}@z�H@{��@{�F@{C�@{@z��@xĜ@w�;@xA�@w�;@w��@w
=@v�@vȴ@v�+@v5?@v@u��@tI�@qG�@p��@q�^@pr�@pA�@j�!@f��@c�
@co@b�\@ax�@`A�@]@]?}@\�@\j@\j@\9X@\Z@\I�@\j@\Z@\z�@\j@\j@\j@\�D@\��@]O�@]p�@]V@\�@\��@\�@[��@X�u@XbN@XA�@YX@ZM�@[@[C�@Z��@Z-@Y�@Y�7@X�9@W�;@V��@V�R@V�+@VE�@W+@V��@V5?@V@U�T@U�T@U�@T�@T1@R�@RM�@Q7L@P�u@Pr�@O�P@O|�@O+@Nff@M�@M�T@M�@M`B@KdZ@J�@J^5@Ko@Lz�@Lj@L�@L�@L1@L�@LZ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�Ĝ@�M�@��\@��@�33@���@�&�@���@�V@�bN@�7L@ǥ�@�ȴ@�?}@�o@�;d@�E�@��@�|�@�`B@�@�
=@�C�@�@� �@��@�5?@�hs@�dZ@�"�@�I�@�%@�1'@�G�@~�@~{@|�D@}?}@}@}`B@~{@}?}@z�H@{��@{�F@{C�@{@z��@xĜ@w�;@xA�@w�;@w��@w
=@v�@vȴ@v�+@v5?@v@u��@tI�@qG�@p��@q�^@pr�@pA�@j�!@f��@c�
@co@b�\@ax�@`A�@]@]?}@\�@\j@\j@\9X@\Z@\I�@\j@\Z@\z�@\j@\j@\j@\�D@\��@]O�@]p�@]V@\�@\��@\�@[��@X�u@XbN@XA�@YX@ZM�@[@[C�@Z��@Z-@Y�@Y�7@X�9@W�;@V��@V�R@V�+@VE�@W+@V��@V5?@V@U�T@U�T@U�@T�@T1@R�@RM�@Q7L@P�u@Pr�@O�P@O|�@O+@Nff@M�@M�T@M�@M`B@KdZ@J�@J^5@Ko@Lz�@Lj@L�@L�@L1@L�@LZ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB\B�{B`BB	}�B
�hB
�XBB�B��B�;B��B%B�B�B!�B�B�B/B�B��B��B�BȴB�B��B�ZB�B�ZB�B%B+BPBuB�B �B,B9XB<jB@�B@�BG�BH�BF�BH�BH�BI�BH�BH�BN�BI�BM�BM�BO�BO�BO�BO�BP�BP�BO�BN�BO�BJ�BJ�BH�BH�BI�BG�BG�B:^B6FB7LB6FB7LB49B49B5?B5?B8RB9XB;dB;dB<jB<jB=qB<jB>wB=qB>wB>wBB�BB�BC�BD�BE�BB�BB�BC�B?}B?}B@�BC�BF�BF�BG�BG�BG�BI�BD�BF�BC�BD�BB�BC�BD�BD�BD�BC�BD�BD�BE�BD�BB�BB�B>wB>wB=qB=qB>wB<jB;dB:^B<jB9XB9XB7LB7LB8RB49B5?B5?B<jB9XB;dB:^B;dB;dB=q11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B\B�{B`BB	}�B
�hB
�XBB�B��B�;B��B%B�B�B!�B�B�B/B�B��B��B�BȴB�B��B�ZB�B�ZB�B%B+BPBuB�B �B,B9XB<jB@�B@�BG�BH�BF�BH�BH�BI�BH�BH�BN�BI�BM�BM�BO�BO�BO�BO�BP�BP�BO�BN�BO�BJ�BJ�BH�BH�BI�BG�BG�B:^B6FB7LB6FB7LB49B49B5?B5?B8RB9XB;dB;dB<jB<jB=qB<jB>wB=qB>wB>wBB�BB�BC�BD�BE�BB�BB�BC�B?}B?}B@�BC�BF�BF�BG�BG�BG�BI�BD�BF�BC�BD�BB�BC�BD�BD�BD�BC�BD�BD�BE�BD�BB�BB�B>wB>wB=qB=qB>wB<jB;dB:^B<jB9XB9XB7LB7LB8RB49B5?B5?B<jB9XB;dB:^B;dB;dB=q11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              201806231407312018062314073120180623140731  IF  ARFMCODA012d                                                                20170820000540                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170820000550                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170820000550                      G�O�G�O�G�O�                IF  ARFMCODA012d                                                                20170831155245                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170831155254                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170831155254                      G�O�G�O�G�O�                IF  ARFMCODA013a                                                                20170908134305                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170908134314                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170908134314                      G�O�G�O�G�O�                IF  ARFMCODA013d                                                                20170918134741                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170918134752                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170918134752                      G�O�G�O�G�O�                IF  ARFMCODA013f                                                                20170928133732                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170928133741                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170928133741                      G�O�G�O�G�O�                IF  ARFMCODA014g                                                                20171008124030                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171008124041                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171008124041                      G�O�G�O�G�O�                IF  ARFMCODA014g                                                                20171008133704                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171008133713                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171008133713                      G�O�G�O�G�O�                IF  ARFMCODA014g                                                                20171018133815                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171018133825                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171018133825                      G�O�G�O�G�O�                IF  ARFMCODA014g                                                                20171028133854                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171028133901                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171028133901                      G�O�G�O�G�O�                IF  ARFMCODA015a                                                                20171107133915                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20171107133928  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.1                                                                 20171107133928  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.0 CTD2017V1                                                       20171114123404  IP  PSAL            @陚C�G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20180623140731  IP  PSAL            @陚C�G�O�                