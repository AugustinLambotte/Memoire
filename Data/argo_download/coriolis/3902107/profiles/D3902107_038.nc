CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T15:12:16Z creation; 2020-11-23T11:33:26Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8(   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8h   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9,   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         90   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    98   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9<   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9D   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9L   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9T   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9X   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9`   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9d   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9h   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9l   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :l   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        <  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        <  C|   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  J�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     <  L�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  S�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  [    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  \�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  d   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  e�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  m   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  tT   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  v$   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  }`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     <  0   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �(   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �,   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �<   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �@   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �D   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �H   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �l   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828151216  20220127170814  3902107 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               &A   IF                                  2C  D   ARVOR                           AI2600-18EU007                  5900A04                         844 @�$�q�1   @�$�q�@SϿ���@��( �8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A   A0  A>ffAQ��Aa��AnffA~ffA�  A���A���A�33A���A���A�ffA�  A�ffA���A�33Aݙ�A�  A�A�  A�ffB  BffB��B��B  B  BffB33B#33B'��B+��B/��B3��B8  B;��B@  BD  BHffBK33BO33BS��BW��B\ffB_33Bc��Bh  Bk33Bp  BtffBw��Bz��B33B�  B���B�  B���B�33B���B���B�  B���B�ffB�  B���B���B�33B���B�ffB�33B�  B���B���B�ffB���B���B�ffB�33B���B���B�ffB�33B���B���B�33B�33B���Bř�B�ffB���B�33B���B�ffB�  B�33B���B�ffB�  B�ffB�  B���B�33B���C��C� CL�C33C	� C��C��C� CL�C33C  CffC�3C��C� CffC!L�C#33C%33C'�C)�C+  C-ffC/��C1�3C3� C5� C7ffC9ffC;ffC=L�C?L�CA33CC��CE��CG��CI��CK��CM��COL�CQL�CSL�CU� CW��CYL�C[L�C]� C_� Ca�3CcL�CeffCg��Ci33CkffCm� CoL�Cq� Cs�3CuffCw33CyffC{��C}� CL�C���C��3C���C�� C��fC���C��3C���C�� C��fC���C�� C��3C�ٚC���C��3C���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C���C�� C��3C��fC��fC��fC��fC��fC��fC��3C��3C��3C��3C�� C�� C�� C�� C���C���C�ٚC��fC��3C��3C�� C���C��fC��fC��3C�� C�ٚC��3C�� C���C��fC³3C�� Cę�Cų3C�� CǦfC�� C���Cʳ3Cˌ�C̦fC�� CΦfCό�Cг3C���CҦfCӌ�CԦfC�� Cֳ3Cי�C�� C��fC���C�� CܦfC݌�C޳3C�ٚC�� C�3C�fC��C� C�3C��fC�ٚC���C�� C�3C�fC왚C홚C��C� C�3C��fC�ٚC���C�� C�� C��3C���C�� C��3C�L�C��3D L�D��D��DfD@ Dy�D�3D��D
&fD` D�fD��D9�D�fD��D��DFfD�fD� D  D9�D� D�fD�3D9�DffD �3D"  D#33D$ffD%�3D'fD(@ D)y�D*��D+�3D-,�D.� D/�fD0��D2,�D3ffD4��D6�D7S3D8y�D9�fD:�fD<,�D=l�D>�3D@  DA33DBl�DC� DD��DF9�DGy�DH��DJ  DK,�DLy�DM� DN�3DP&fDQs3DR�fDS��DU33DVy�DW��DX�3DZ,�D[ffD\�fD]�fD_9�D`�3Da��Db�fDd,�De� Df�3Dg��Di,�Djl�Dk�3Dl��Dn@ Do�fDp��Dq�3Ds@ Dt��Du��Dv��Dx@ Dys3Dz��D{�fD}  D~y�D�3D�|�D�fD���D�Y�D��3D�� D�<�D�ٚD��3D�&fD��3D�ffD���D���D�@ D��fD�|�D�#3D���D�Y�D��fD��3D�<�D��fD�|�D�3D��fD�` D���D��fD�@ D�� D�� D�  D��3D�Y�D���D��3D�FfD�� D�y�D�fD�� D�c3D�3D��fD�FfD��D�|�D�3D���D�` D���D��3D�@ D���D�y�D�fD��fD�Y�D���D���D�@ D��3D�y�D�  D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A   A0  A>ffAQ��Aa��AnffA~ffA�  A���A���A�33A���A���A�ffA�  A�ffA���A�33Aݙ�A�  A�A�  A�ffB  BffB��B��B  B  BffB33B#33B'��B+��B/��B3��B8  B;��B@  BD  BHffBK33BO33BS��BW��B\ffB_33Bc��Bh  Bk33Bp  BtffBw��Bz��B33B�  B���B�  B���B�33B���B���B�  B���B�ffB�  B���B���B�33B���B�ffB�33B�  B���B���B�ffB���B���B�ffB�33B���B���B�ffB�33B���B���B�33B�33B���Bř�B�ffB���B�33B���B�ffB�  B�33B���B�ffB�  B�ffB�  B���B�33B���C��C� CL�C33C	� C��C��C� CL�C33C  CffC�3C��C� CffC!L�C#33C%33C'�C)�C+  C-ffC/��C1�3C3� C5� C7ffC9ffC;ffC=L�C?L�CA33CC��CE��CG��CI��CK��CM��COL�CQL�CSL�CU� CW��CYL�C[L�C]� C_� Ca�3CcL�CeffCg��Ci33CkffCm� CoL�Cq� Cs�3CuffCw33CyffC{��C}� CL�C���C��3C���C�� C��fC���C��3C���C�� C��fC���C�� C��3C�ٚC���C��3C���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��fC�ٚC���C���C�� C��3C��fC��fC��fC��fC��fC��fC��3C��3C��3C��3C�� C�� C�� C�� C���C���C�ٚC��fC��3C��3C�� C���C��fC��fC��3C�� C�ٚC��3C�� C���C��fC³3C�� Cę�Cų3C�� CǦfC�� C���Cʳ3Cˌ�C̦fC�� CΦfCό�Cг3C���CҦfCӌ�CԦfC�� Cֳ3Cי�C�� C��fC���C�� CܦfC݌�C޳3C�ٚC�� C�3C�fC��C� C�3C��fC�ٚC���C�� C�3C�fC왚C홚C��C� C�3C��fC�ٚC���C�� C�� C��3C���C�� C��3C�L�C��3D L�D��D��DfD@ Dy�D�3D��D
&fD` D�fD��D9�D�fD��D��DFfD�fD� D  D9�D� D�fD�3D9�DffD �3D"  D#33D$ffD%�3D'fD(@ D)y�D*��D+�3D-,�D.� D/�fD0��D2,�D3ffD4��D6�D7S3D8y�D9�fD:�fD<,�D=l�D>�3D@  DA33DBl�DC� DD��DF9�DGy�DH��DJ  DK,�DLy�DM� DN�3DP&fDQs3DR�fDS��DU33DVy�DW��DX�3DZ,�D[ffD\�fD]�fD_9�D`�3Da��Db�fDd,�De� Df�3Dg��Di,�Djl�Dk�3Dl��Dn@ Do�fDp��Dq�3Ds@ Dt��Du��Dv��Dx@ Dys3Dz��D{�fD}  D~y�D�3D�|�D�fD���D�Y�D��3D�� D�<�D�ٚD��3D�&fD��3D�ffD���D���D�@ D��fD�|�D�#3D���D�Y�D��fD��3D�<�D��fD�|�D�3D��fD�` D���D��fD�@ D�� D�� D�  D��3D�Y�D���D��3D�FfD�� D�y�D�fD�� D�c3D�3D��fD�FfD��D�|�D�3D���D�` D���D��3D�@ D���D�y�D�fD��fD�Y�D���D���D�@ D��3D�y�D�  D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���k����j�ȴ��j?$�?g+?h1'?ix�?���?�~�?�;d?��?�$�?��?�j?�V?͑h?Ͼw?ѩ�?щ7?�?ڟ�?�!@�m@�H@~�@z�@�@bN@��@@|�@#�F@%`B@(��@"�\@ Q�@%��@,�/@0��@8��@8��@?+@?�P@@��@Cƨ@D�j@A��@>��@;�
@8�9@=��@>�y@=�h@?�w@>E�@;"�@8�u@;"�@7�@8bN@8bN@9��@;��@<��@>@=@<�/@=`B@;�
@:-@;@<I�@;��@:-@8b@5O�@2�!@333@3"�@2��@2�\@2^5@3"�@4�@3t�@2=q@.$�@*=q@(�9@(�9@*�\@/
=@0 �@1%@/l�@/�@/\)@.��@)��@+dZ@,�/@-��@-O�@,��@+�@)�@)�7@'�@�;@S�@��@bN@�@��@?}@p�@�@�@p�@z�@@`B@��@��@�+@b@�j@
=@�j@v�@ 1'@!%@!G�@!G�@!�^@ȴ@/@�+@?}@Z@/@�h@�@ �`@!&�@!�#@"=q@!��@"�H@!��@��@\)@K�@��@�\@V@@p�@C�@hs@G�@�
@ff@7L@bN@�@?}@�@�\@��@&�@�9@�@��@�7@G�@�u@Ĝ@x�@��@Q�@K�@�@O�@/@�@�j@C�@
�@	��@Q�@l�@ff@��@p�@/@�@�F@C�@�!@~�@^5@J@�@7L?���?�p�?��m?�"�?��^?�l�?�?���?�t�?�!?�7?� �?���?�O�?�D?��?�^5?�x�?�r�?���?��?��?�Z?㕁?�Ĝ?�  ?�5??ܬ?�j?�dZ?�~�?�X?և+?�?}?���?�S�?�J?� �?Η�?Ͳ-?̬?˅?��?ɺ^?���?Ƈ+?š�?��/?öF?�G�?���?�V?���?��?�
=?�ff?��
?�%?���?��?�{?�V?���?�^5?��u?��P?�E�?���?�9X?�-?�&�?��?��h?���?�=q?��?��/?�Ĝ?�G�?�bN?�V?�?�Q�?��y?��?�S�?;d?r�!?n��?i��?dZ?[��?Q��?J��?C�
?7��?&$�?�?��?p�?33?V?�T>��>��>�1>�x�>��>��H>��+>fff>F��>?|�>:^5>��>�=��=� �=�C�=aG�=<j=}�=�7L=�t�=��=t���1��w�o�ě����
�t��y�#��t���{���������S����#�$ݾz�!���(�þ+�2-�;dZ�B�\�G��S�ϾY��dZ�n���vȴ��  ��o����7L��O߾�녾�񪾖��"Ѿ�A����y��1�������!��Q쾽�\�ě��Ƨ��1'�������;�����޸R���
����������ٿJ������9�ƨ����!�����j��R��w�!�7�#���$Z�$Z�%��'l��(�ÿ)7L�+�-V�-��-��.{�/��0bN�0�׿1&�3t��4���5?}�6�6ȴ�7Kǿ7�ٿ7�P�7�P�7�P�8b�9��9�#�:�H�9�#�:�H�:���<(��<푿>vɿ@  �AG��BJ�B�\�B�\�Co�C�
�D��D��DZ�D��C���C���DZ�E��E�˿E�˿E�˿E`B�E`B�E`B�E`B�E��E��E�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 �k����j�ȴ��j?$�?g+?h1'?ix�?���?�~�?�;d?��?�$�?��?�j?�V?͑h?Ͼw?ѩ�?щ7?�?ڟ�?�!@�m@�H@~�@z�@�@bN@��@@|�@#�F@%`B@(��@"�\@ Q�@%��@,�/@0��@8��@8��@?+@?�P@@��@Cƨ@D�j@A��@>��@;�
@8�9@=��@>�y@=�h@?�w@>E�@;"�@8�u@;"�@7�@8bN@8bN@9��@;��@<��@>@=@<�/@=`B@;�
@:-@;@<I�@;��@:-@8b@5O�@2�!@333@3"�@2��@2�\@2^5@3"�@4�@3t�@2=q@.$�@*=q@(�9@(�9@*�\@/
=@0 �@1%@/l�@/�@/\)@.��@)��@+dZ@,�/@-��@-O�@,��@+�@)�@)�7@'�@�;@S�@��@bN@�@��@?}@p�@�@�@p�@z�@@`B@��@��@�+@b@�j@
=@�j@v�@ 1'@!%@!G�@!G�@!�^@ȴ@/@�+@?}@Z@/@�h@�@ �`@!&�@!�#@"=q@!��@"�H@!��@��@\)@K�@��@�\@V@@p�@C�@hs@G�@�
@ff@7L@bN@�@?}@�@�\@��@&�@�9@�@��@�7@G�@�u@Ĝ@x�@��@Q�@K�@�@O�@/@�@�j@C�@
�@	��@Q�@l�@ff@��@p�@/@�@�F@C�@�!@~�@^5@J@�@7L?���?�p�?��m?�"�?��^?�l�?�?���?�t�?�!?�7?� �?���?�O�?�D?��?�^5?�x�?�r�?���?��?��?�Z?㕁?�Ĝ?�  ?�5??ܬ?�j?�dZ?�~�?�X?և+?�?}?���?�S�?�J?� �?Η�?Ͳ-?̬?˅?��?ɺ^?���?Ƈ+?š�?��/?öF?�G�?���?�V?���?��?�
=?�ff?��
?�%?���?��?�{?�V?���?�^5?��u?��P?�E�?���?�9X?�-?�&�?��?��h?���?�=q?��?��/?�Ĝ?�G�?�bN?�V?�?�Q�?��y?��?�S�?;d?r�!?n��?i��?dZ?[��?Q��?J��?C�
?7��?&$�?�?��?p�?33?V?�T>��>��>�1>�x�>��>��H>��+>fff>F��>?|�>:^5>��>�=��=� �=�C�=aG�=<j=}�=�7L=�t�=��=t���1��w�o�ě����
�t��y�#��t���{���������S����#�$ݾz�!���(�þ+�2-�;dZ�B�\�G��S�ϾY��dZ�n���vȴ��  ��o����7L��O߾�녾�񪾖��"Ѿ�A����y��1�������!��Q쾽�\�ě��Ƨ��1'�������;�����޸R���
����������ٿJ������9�ƨ����!�����j��R��w�!�7�#���$Z�$Z�%��'l��(�ÿ)7L�+�-V�-��-��.{�/��0bN�0�׿1&�3t��4���5?}�6�6ȴ�7Kǿ7�ٿ7�P�7�P�7�P�8b�9��9�#�:�H�9�#�:�H�:���<(��<푿>vɿ@  �AG��BJ�B�\�B�\�Co�C�
�D��D��DZ�D��C���C���DZ�E��E�˿E�˿E�˿E`B�E`B�E`B�E`B�E��E��E�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB!�B\)B�5Bq�B��B��B��B�B33B��B	 �B	XB	�!B	��B	��B	��B	�B	�5B	�TB	�B	��B
{B
A�B
iyB
x�B
y�B
� B
�\B
��B
B
�B
�BB�B �B�B�B2-BE�BJ�B]/BffBl�By�B{�B� B�DB�DB�DB�B�7B��B��B�uB��B��B��B��B��B��B��B��B��B��B�B�!B�?B�?B�RB�^B�FB�XB�jB�qB�?B�RB�wB�9B�LB�RB�LB�RB�RB�RB�^B�RB�LB�RB�?B�'B�B�3B�dB�wB��B��B��BĜB��B��B��BĜBƨBɺBƨBĜBBBB�9B�B��B�B��B�B�B�B�!B�B�!B�3B�B�3B�FB�XB�XB�qBÖBȴBŢBŢB��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�/B�BB�5B�5B�)B�#B�B�B��B��B��BɺB��BǮBƨB�?B��BƨB��B��B�
B��B��B��B��B��B��B��B��B��B��B��B�B�
B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBȴBɺBȴBɺBȴBɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBȴBǮBǮBǮBȴBȴBȴBɺBɺBɺBɺBǮBȴBȴBɺBǮBǮBȴBƨBƨBĜBŢBƨBŢBÖBĜBŢBĜBĜBĜBĜBÖBĜBŢBĜBÖBĜBĜBĜBĜBŢBBBƨBB��B�dBĜB��B�}B��B��B�}B�qB�wB�jB�qB�dB�XB�^B�LB�FB�3B�!B�B�!B�3B�-B�'B�B�B�B�!B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B!�B\)B�5Bq�B��B��B��B�B33B��B	 �B	XB	�!B	��B	��B	��B	�B	�5B	�TB	�B	��B
{B
A�B
iyB
x�B
y�B
� B
�\B
��B
B
�B
�BB�B �B�B�B2-BE�BJ�B]/BffBl�By�B{�B� B�DB�DB�DB�B�7B��B��B�uB��B��B��B��B��B��B��B��B��B��B�B�!B�?B�?B�RB�^B�FB�XB�jB�qB�?B�RB�wB�9B�LB�RB�LB�RB�RB�RB�^B�RB�LB�RB�?B�'B�B�3B�dB�wB��B��B��BĜB��B��B��BĜBƨBɺBƨBĜBBBB�9B�B��B�B��B�B�B�B�!B�B�!B�3B�B�3B�FB�XB�XB�qBÖBȴBŢBŢB��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�/B�BB�5B�5B�)B�#B�B�B��B��B��BɺB��BǮBƨB�?B��BƨB��B��B�
B��B��B��B��B��B��B��B��B��B��B��B�B�
B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBȴBɺBȴBɺBȴBɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBȴBǮBǮBǮBȴBȴBȴBɺBɺBɺBɺBǮBȴBȴBɺBǮBǮBȴBƨBƨBĜBŢBƨBŢBÖBĜBŢBĜBĜBĜBĜBÖBĜBŢBĜBÖBĜBĜBĜBĜBŢBBBƨBB��B�dBĜB��B�}B��B��B�}B�qB�wB�jB�qB�dB�XB�^B�LB�FB�3B�!B�B�!B�3B�-B�'B�B�B�B�!B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               No adjustment was necessary. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                                                      202011231133262022012717081420220127170814  IF  ARFMCODA035h                                                                20200828151216                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828151252  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828151252  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201123113326  IP  PSAL            A   D���G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170814  IP  PSAL            A   D���G�O�                