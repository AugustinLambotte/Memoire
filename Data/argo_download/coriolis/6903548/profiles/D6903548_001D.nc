CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  U   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:48Z creation; 2021-06-07T15:42:38Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8    FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8(   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    88   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8H   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8P   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9    	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9,   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    90   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     94   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9T   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9t   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
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
resolution        =���   axis      Z        T  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  @,   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        T  A�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  F�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     T  H0   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  M�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  R�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  T0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  Y�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  Z�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  `0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  e�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  f�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 X  l0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  m�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    |8   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    |<   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    |@   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    |D   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  |H   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    |�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    |�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    |�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         |�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         |�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        |�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    |�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  r�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    s   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    v   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    y   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  |Argo profile    3.1 1.2 19500101000000  20190620115648  20210607154238  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @�p��-��1   @�p�B�\ @TN�=��B@.W�﫭Q1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 100 dbar; 10 sec sampling, 2 dbar average from 100 dbar to 500 dbar; 10 sec sampling, 5 dbar average from 500 dbar to 1000 dbar]                                                       A0  AH  Ad��A�33A�ffA�33A�33A�33A�  A���A�33B��B	33B��B��B!33B(ffB133B:ffBBffBJffBR  BZ  Ba��Bi33Br  Bz  B���B�  B�  B�33B�  B���B���B���B���B���B���B���B���B���B���B�  B���Bę�B�ffB���B�33B�ffBؙ�B���B�ffB䙚B���B�ffB�B���B�ffB���C L�C� C33C33CL�C
ffCL�CL�C� CL�C� CL�CffCL�CffCL�C ffC"33C$33C&ffC(ffC*L�C,� C.L�C0ffC2ffC4L�C6L�C8L�C:33C<ffC>ffC@ffCBffCDffCF� CHffCJL�CLffCNL�CPL�CRL�CTffCVffCXL�CZ33C\ffC^33C`� Cb��Cd��Cf� ChffCjL�ClL�Cn33CpffCr��Ct� CvL�CxffCz��C|� C~L�C�@ C�L�C�33C�@ C�@ C�@ C�L�C�&fC�33C�33C�33C�33C�33C�@ C�L�C�33C�&fC�33C�33C�L�C�33C�&fC�33C�@ C�@ C�@ C�@ C�33C�33C�&fC�@ C�@ C��C�@ C�&fC�&fC��C��C�@ C�&fC�33C�@ C�&fC��C�@ C�L�C�33C�33C�33C��C�&fC�33C�@ C�&fC�33C��C�@ C�L�C�&fC��C�33C�L�C�33C�&fC�33C�33C��C�&fC��C�@ C��C��C�L�C�&fC�&fC�33C�33C�&fC�L�C��C��C�L�C�L�C�@ C�33C�&fC��C�33C��C�&fC�&fC�&fC��C�33C�@ C�33C�&fC��C�&fC�33C�&fC��C�33C�&fC�33C�&fC�&fC�33C�33C�&fC�&fC��C�&fC�&fC�&fC��C�33C�33C��C�&fC�&fC�&fC�ٚC�� C�ٚD&fDy�D�3D  D33Dy�D�fD	��D@ Ds3D�3D��D33Dy�D�3D�3D,�Dl�D��D�3D9�Dy�D�3D�3D9�D y�D!��D"��D$33D%l�D&��D'�3D)33D*� D+��D,��D.9�D/s3D0��D2  D3FfD4y�D5�3D6��D89�D9l�D:��D;�3D=,�D>y�D?��D@�3DB33DCl�DD��DE��DG33DHy�DI��DJ��DL9�DMs3DN� DO��DQ9�DRy�DS�3DT�3DV@ DWy�DX��DY�3D[33D\l�D]��D^�3D`33Day�Db��Dc�3De@ Dfy�Dg��Dh�3Dj33Dks3Dl��Dm�3Do@ Dp� Dq�3Dr�3Dt�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A0  AH  Ad��A�33A�ffA�33A�33A�33A�  A���A�33B��B	33B��B��B!33B(ffB133B:ffBBffBJffBR  BZ  Ba��Bi33Br  Bz  B���B�  B�  B�33B�  B���B���B���B���B���B���B���B���B���B���B�  B���Bę�B�ffB���B�33B�ffBؙ�B���B�ffB䙚B���B�ffB�B���B�ffB���C L�C� C33C33CL�C
ffCL�CL�C� CL�C� CL�CffCL�CffCL�C ffC"33C$33C&ffC(ffC*L�C,� C.L�C0ffC2ffC4L�C6L�C8L�C:33C<ffC>ffC@ffCBffCDffCF� CHffCJL�CLffCNL�CPL�CRL�CTffCVffCXL�CZ33C\ffC^33C`� Cb��Cd��Cf� ChffCjL�ClL�Cn33CpffCr��Ct� CvL�CxffCz��C|� C~L�C�@ C�L�C�33C�@ C�@ C�@ C�L�C�&fC�33C�33C�33C�33C�33C�@ C�L�C�33C�&fC�33C�33C�L�C�33C�&fC�33C�@ C�@ C�@ C�@ C�33C�33C�&fC�@ C�@ C��C�@ C�&fC�&fC��C��C�@ C�&fC�33C�@ C�&fC��C�@ C�L�C�33C�33C�33C��C�&fC�33C�@ C�&fC�33C��C�@ C�L�C�&fC��C�33C�L�C�33C�&fC�33C�33C��C�&fC��C�@ C��C��C�L�C�&fC�&fC�33C�33C�&fC�L�C��C��C�L�C�L�C�@ C�33C�&fC��C�33C��C�&fC�&fC�&fC��C�33C�@ C�33C�&fC��C�&fC�33C�&fC��C�33C�&fC�33C�&fC�&fC�33C�33C�&fC�&fC��C�&fC�&fC�&fC��C�33C�33C��C�&fC�&fC�&fC�ٚC�� C�ٚD&fDy�D�3D  D33Dy�D�fD	��D@ Ds3D�3D��D33Dy�D�3D�3D,�Dl�D��D�3D9�Dy�D�3D�3D9�D y�D!��D"��D$33D%l�D&��D'�3D)33D*� D+��D,��D.9�D/s3D0��D2  D3FfD4y�D5�3D6��D89�D9l�D:��D;�3D=,�D>y�D?��D@�3DB33DCl�DD��DE��DG33DHy�DI��DJ��DL9�DMs3DN� DO��DQ9�DRy�DS�3DT�3DV@ DWy�DX��DY�3D[33D\l�D]��D^�3D`33Day�Db��Dc�3De@ Dfy�Dg��Dh�3Dj33Dks3Dl��Dm�3Do@ Dp� Dq�3Dr�3Dt�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?��9?�"�?���?���?��?�bN?���?���@�@
~�@�j@��?���?�  ?�33?У�?�1'@/@&�@"J@+��@4�@5�@%��@�\@33@��@�@�u@�9@1@�m@O�@@�@��@�@  @)hs@.�+@-�@ 1'@ �`@!��@!�^@"��@( �@)�@)x�@*J@+��@1��@49X@4z�@4(�@2n�@0��@0��@0�`@0�`@0�`@/��@.�R@-�-@,��@,1@)G�@';d@'\)@&�y@)�^@.5?@/;d@.��@.E�@-�@,�@+�F@*�@*J@)��@)��@)hs@)��@*^5@*^5@*^5@*^5@)��@)�7@'��@&�+@&ff@&V@&$�@&@%�@%�@$(�@"��@"~�@"�!@#t�@$�D@%p�@%�@&�@';d@'�;@(1'@( �@(b@(b@( �@( �@(b@(Q�@'�@&{@%O�@!�@�@I�@(�@1@�@O�@�
@�y@1@~�@�@��@��@+@ �`@!hs@�R@�\@r�@l�@;d@�j@�\@J@��@��@\)@ �@A�@ �@��@hs@�@�H@"�@o@o@��@��@�D@ff@��@v�@p�@��@�@A�@��@��@��@�@z�@��@��@�7@�^@�^@�@�#@��@�7@�7@x�@�@�@�w@|�@|�@|�@��@�w@b@�;@��@\)@E�@�-@?}@�@�-@V@ȴ@�R@�@;d@�@�@b@1'@A�@b@b@ �@+@�-@�@9X@I�@dZ@
n�@	��@
J@
�!@r�@b@�@�m@��@ �u?��?�/?��?�"�?���?���?�"�?�(�?�^5?��9?�ȴ?�?�?}?��j?���?�?�Z?�;d?��?�?�Q�?�7?ش9?���?̋D?�t�?�?���?�\)?�^5?��P?��?�33?� �?���?���?�Q�?��?��
?�Ĝ?��h?��?�33?�bN?xb?q&�?co?V?MO�?G�?>5??6�+?2n�?.��?*=q?%�T?|�?X??&�?�D?	��?�
>�>�1>��
>��>Ǯ>�>��/>�M�>�;d>���>��9>��7>r�!>_;d>R�>6E�>"��>
=q=��=�
==�-=��-=L��<��ͼ49X��/��w�]/�aG��ixսm�h��O߽�1�Ƨ���������h�+�z������$�/�/��7KǾ?|�Kƨ�W
=�aG�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ?��9?�"�?���?���?��?�bN?���?���@�@
~�@�j@��?���?�  ?�33?У�?�1'@/@&�@"J@+��@4�@5�@%��@�\@33@��@�@�u@�9@1@�m@O�@@�@��@�@  @)hs@.�+@-�@ 1'@ �`@!��@!�^@"��@( �@)�@)x�@*J@+��@1��@49X@4z�@4(�@2n�@0��@0��@0�`@0�`@0�`@/��@.�R@-�-@,��@,1@)G�@';d@'\)@&�y@)�^@.5?@/;d@.��@.E�@-�@,�@+�F@*�@*J@)��@)��@)hs@)��@*^5@*^5@*^5@*^5@)��@)�7@'��@&�+@&ff@&V@&$�@&@%�@%�@$(�@"��@"~�@"�!@#t�@$�D@%p�@%�@&�@';d@'�;@(1'@( �@(b@(b@( �@( �@(b@(Q�@'�@&{@%O�@!�@�@I�@(�@1@�@O�@�
@�y@1@~�@�@��@��@+@ �`@!hs@�R@�\@r�@l�@;d@�j@�\@J@��@��@\)@ �@A�@ �@��@hs@�@�H@"�@o@o@��@��@�D@ff@��@v�@p�@��@�@A�@��@��@��@�@z�@��@��@�7@�^@�^@�@�#@��@�7@�7@x�@�@�@�w@|�@|�@|�@��@�w@b@�;@��@\)@E�@�-@?}@�@�-@V@ȴ@�R@�@;d@�@�@b@1'@A�@b@b@ �@+@�-@�@9X@I�@dZ@
n�@	��@
J@
�!@r�@b@�@�m@��@ �u?��?�/?��?�"�?���?���?�"�?�(�?�^5?��9?�ȴ?�?�?}?��j?���?�?�Z?�;d?��?�?�Q�?�7?ش9?���?̋D?�t�?�?���?�\)?�^5?��P?��?�33?� �?���?���?�Q�?��?��
?�Ĝ?��h?��?�33?�bN?xb?q&�?co?V?MO�?G�?>5??6�+?2n�?.��?*=q?%�T?|�?X??&�?�D?	��?�
>�>�1>��
>��>Ǯ>�>��/>�M�>�;d>���>��9>��7>r�!>_;d>R�>6E�>"��>
=q=��=�
==�-=��-=L��<��ͼ49X��/��w�]/�aG��ixսm�h��O߽�1�Ƨ���������h�+�z������$�/�/��7KǾ?|�Kƨ�W
=�aG�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBM�B]/BQ�B��B�)B��B	Q�B	x�B	��B	��B	��B
_;B
�hB
�B
u�B
z�B
��B
�/B
�yBhB49BN�Be`BffBcTBM�BL�BK�BK�BQ�B\)Bk�Bx�B�PB�B|�By�B�B�^B��B��B�wB�wB��BB��B��B��B�B�B�fB�sB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�TB�ZB�`B�mB�B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�`B�mB�mB�sB�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�mB�5B�BB�BB�BB�HB�NB�NB��B��B��B��B�)B�HB�fB�mB�B�fB�/B�/B�#B�B�
B��B��B��B��B�B�)B�/B�#B�/B�5B�BB�BB�TB�TB�TB�TB�BB�ZB�mB�B�B�sB�mB�`B�;B�;B�5B�/B�#B�B�B��B��B��B��B�B��B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�
B�
B�
B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺB��BȴBÖB��B��B�}B�qB�^B�RB�LB�LB�LB�RB�LB�FB�FB�9B�?B�LB�9B�3B�-B�-B�!B�!B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   BQ_B`�BUxB�,BߵB	 zB	UxB	|aB	��B	�~B	�gB
b�B
��B
��B
yOB
~mB
�]B
�B
�B�B7�BReBh�Bi�Bf�BQ_BPYBOSBOSBUxB_�BoB|aB��B��B�zB}gB��B��B�SB�YB�B�B�B�B�B�qBׄBݩBۜB��B��B�*B�6B�<B�BB�0B�*B�0B�6B�6B�0B�B�B�B�0B�B��B��B��B��B�B�BB�OB�OB�BB�<B�6B�6B�6B�0B�0B�0B�*B�6B�6B�6B�6B�<B�*B�6B�$B�B�B�B�B�B�B�B��B��B��B��B�B�B�B�B�0B�<B�<B�BB�<B�<B�BB�BB�BB�BB�HB�<B�*B��B��B��B��B��B��B��B��BׄB�xB�qB؊BߵB��B��B��B�B��B�B�BޯBݩBږB؊B�xB�qB�xBۜBߵB�BޯB�B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B�BޯBܣBۜB؊B؊B؊B؊BِB؊B؊BِB؊BِBׄB�~B�~B�xB�xB�xB�~B�~BׄBׄBׄB�~B�xB�qB�eB�eB�kB�eB�xB�~B�xBׄBׄBِBِBږBږBږBږBۜBܣBׄB�~B�kB�~B�~B�xB�qB�kBׄB�eB�qBׄB�eB�SB�MB�SB�FB�FB�SB�MB�MB�SB�eB�eB�_B�YB�YB�SB�SB�MB�MB�MB�FB�FB�FB�MB�@B�"B�B�B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�|B�|B�uB�|B�uB�uB�oB�uB�oB�oB�oB�oB�iB�iB�cB�cB�oB�]B�]B�]B�cB�iB�]B�]B�WB�QB�WB�JB�JB�DB�DB�DB�JB�DB�JB�DB�>B�DB�DB�DB�>B�>B�DB�DB�DB�DB�DB�DB�DB�>B�DB�DB�DB�DB�DB�DB�DB�DB�>11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542382021060715423820210607154238  IF  ARFMCODA029d                                                                20190620115648                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115701  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115701  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154238  IP  PSAL            A0  Dt�G�O�                