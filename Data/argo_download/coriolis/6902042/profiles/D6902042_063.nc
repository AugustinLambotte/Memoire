CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   U   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       S2015-03-25T14:53:58Z creation; 2022-08-24T14:00:31Z last update (BSH ARSQ software)    comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.     comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8\   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9    WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9    JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T      
resolution        >�EȠ�Q)        9$   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9,   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�EȠ�Q)        90   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           98   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9@   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9H   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9L   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9T   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :T   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    :X   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    :\   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :`   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        T  :d   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  ;�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        T  <   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  =d   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     T  =�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  ?   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  @d   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  @�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  B   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  Bh   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  C�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  E   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  Eh   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  F�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  G   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  Hh   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    H�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    K�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    N�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  Q�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    Q�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    Q�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    Q�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    Q�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  Q�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    R   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    R$   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    R(   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         R8   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         R<   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        R@   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    RDArgo profile    3.1 1.2 19500101000000  20150325145358  20220824140031  6902042 E-AIMS                                                          Waldemar Walczowski                                             PRES            TEMP            PSAL               ?A   IF  39212438                        2C  D   NEMO                            298                             n/a                             860 @�B��z1   @�B��z@Si���a@5loTK1   GPS     Primary sampling: discrete []                                                                                                                                                                                                                                      B   B   B               @,��>L��>L��=���@�ffA��A[33A���A�  B33BFffBk33B�  B���B���B�  Bڙ�B���C  C33CffC33C)�C3L�C=33CGL�CQ�C[�3Ce� Cn��Cx� C��3C���C��3C�s3C�� C�� C�s3C��C���C�Y�C�@ C�Y�C�Y�C�&fC�3C��C�33DfD	33D,�D�fD&fD"L�D(�3D.��D;33DG� DS�fD`� DmY�Dyl�D�3D�` D�|�D���D���D�C3D��fD�ٚD���D�l�D���D���D�3D�P Dډ�D��3D�  D�3D�3D�34444444111111111111111111111111111111111111111111111111111111111111111111111111111144   G�O�G�O�G�O�G�O�G�O�G�O�G�O�@�ffA��A[33A���A�  B33BFffBk33B�  B���B���B�  Bڙ�B���C  C33CffC33C)�C3L�C=33CGL�CQ�C[�3Ce� Cn��Cx� C��3C���C��3C�s3C�� C�� C�s3C��C���C�Y�C�@ C�Y�C�Y�C�&fC�3C��C�33DfD	33D,�D�fD&fD"L�D(�3D.��D;33DG� DS�fD`� DmY�Dyl�D�3D�` D�|�D���D���D�C3D��fD�ٚD���D�l�D���D���D�3D�P Dډ�D��3D�  D�3G�O�G�O�4444444111111111111111111111111111111111111111111111111111111111111111111111111111144   G�O�G�O�G�O�G�O�G�O�G�O�G�O�@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�BŢB��B[#B@�BE�Bu�B\BK�BK�BK�BK�BK�BK�BK�BK�BL�BK�BK�BK�BL�BK�BG�BH�BH�BH�BI�BG�BI�BJ�BG�BK�BM�BM�BM�BM�BL�BI�BC�BC�BE�BB�B<jB8RB8RB5?B2-B0!B$�BDB��B��B�B��B��B��B��B�B�B�B��B�jB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��4444444111111111111111111111111111111111111111111111111111111111111111111111111111111   G�O�G�O�G�O�G�O�G�O�G�O�G�O�BK�BK�BK�BK�BK�BK�BK�BK�BL�BK�BK�BK�BL�BK�BG�BH�BH�BH�BI�BG�BI�BJ�BG�BK�BM�BM�BM�BM�BL�BI�BC�BC�BE�BB�B<kB8PB8OB5?B2,B0!B$�BCB��B��B�B��B��B��B��B�B�B�B��B�jB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�4444444111111111111111111111111111111111111111111111111111111111111111111111111111144   G�O�G�O�G�O�G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�?�7?BM�?ǍP@9�@8bN@6��@.ȴ@:^5@:J@:-@9��@:=q@9�@:�\@:�\@:��@:�H@:�!@:�!@;��@;�@8Ĝ@9&�@9x�@9�#@9��@97L@9hs@;S�@9%@;33@<j@;�F@<��@<�/@<��@4I�@4(�@4I�@4I�@4(�@-�@+�@'K�@&��@#t�@ Ĝ@�?�?}?�V?�M�?�\)?�9X?���?�`B?���?{�m?y�?k�>�Q�>`A��%�y�#�������;�ؓu��l���Z�n������R�&��/\)�.V�-O߿4�j�8b�E��C�
�CS��A���@Ĝ�@��@�4444444111111111111111111111111111111111111111111111111111111111111111111111111111111   G�O�G�O�G�O�G�O�G�O�G�O�G�O�@:^5@:J@:-@9��@:=q@9�@:�\@:�\@:��@:�H@:�!@:�!@;��@;�@8Ĝ@9&�@9x�@9�#@9��@97L@9hs@;S�@9%@;33@<j@;�F@<��@<�/@<��@4I�@4(�@4I�@4I�@4(�@-�@+�@'K�@&��@#t�@ Ĝ@�?�?}?�V?�M�?�\)?�9X?���?�`B?���?{�m?y�?k�>�Q�>`A��%�y�#�������;�ؓu��l���Z�n������R�&��/\)�.V�-O߿4�j�8b�E��C�
�CS��A���@Ĝ�@��@�4444444111111111111111111111111111111111111111111111111111111111111111111111111111111   G�O�G�O�G�O�G�O�G�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Surface pressure = 0 dbar                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202208241400312022082414003120220824140031  IF  ARGQCOAR1.0                                                                 20150325152854  QCF$                G�O�G�O�G�O�004100          IF  ARGQCOAR1.0                                                                 20150325152854  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20150325152854  QCC$                G�O�G�O�G�O�004100          IF  CORTCOOA6.2 RTQCGL01                                                        20150326110148  QCF$TEMP            D/  Dz  G�O�6               IF  CORTCOOA6.2 RTQCGL01                                                        20150326111611  QCP$PSAL            G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20161226172805  QCP$PSAL            G�O�G�O�G�O�                IF  CODMCOOA6.2 DMQCGL01                                                        20161226172920  QCF$TEMP            D/  Dz  G�O�6               IF      CORA4.3                                                                 20170804062926  QCP$                G�O�G�O�G�O�                IF  ARGQMIMAR3.2                                                                20190128142326  QCP$                G�O�G�O�G�O�                IF      SCOO0.49                                                                20190925162808  CF  TEMP                    ?�                  IF      SCOO0.49                                                                20190925162808  CF  TEMP            >L��>L��?�                  IF      SCOO0.49                                                                20190925162808  CF  TEMP            @�ff@�ff@�                  IF      SCOO0.49                                                                20190925162808  CF  PSAL                    ?�                  IF      SCOO0.49                                                                20190925162808  CF  PSAL            >L��>L��?�                  IF      SCOO0.49                                                                20190925162808  CF  PSAL            @�ff@�ff@�                  IF      SCOO0.49                                                                20190925162808  CF  PRES                    ?�                  IF      SCOO0.49                                                                20190925162808  CF  PRES            @,��@,��?�                  IF      COFC3.2                                                                 20220822141606                      G�O�G�O�G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20220824140031  IP  PSAL                D�3G�O�                