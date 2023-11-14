CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2020-05-11T14:52:45Z creation; 2020-05-11T14:54:23Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z        L  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  =(   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        L  =�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  @   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     L  @�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  B�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  E4   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  E�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  H�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  J�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M@   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  M�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  P    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     L  P�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    \\   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    \`   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    \d   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    \h   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  \l   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    \�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    \�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    \�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         \�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         \�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        \�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    \�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  S    SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    S0   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    V0   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Y0   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  \0Argo profile    3.1 1.2 19500101000000  20200511145245  20211005090932  6902728 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               VA   IF                                  2C  D   ARVOR                           AI2600-16FR311                  5900A04                         844 @�B''�}(1   @�B''�}(@R1麮t4�,��-�8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 25 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                      A&ffA1��A>ffAP  A`  Ap  A�  A�  A�ffA�ffA�  A���A���A�33A�  A�33A�33A�  A���A�  A�33A�ffB   B  B��B33B��B  BffB33B��B#��B(  B,ffB/��B333B8ffB=33B@  BE33BH  B^  B�  B�  B�33B�ffB�ffB�33B�33C��C� CL�C%L�C/L�C9ffCC��CM�3CW�3Ca��CkffCuL�C��C��3C���C���C���C���C��fC���C��3C��3C�� C���C��3C��3Cř�CʦfCϳ3Cԙ�C٦fC���C�� C�� C�ٚC���C�� D @ Dl�D�3D�3D&fDy�D%� D+��D2�D8` D>� DE�DK,�DQY�DW��D^3Dd33Djy�Dp�3Dv��D}9�D���D���D���D��D�9�D�` D�� D���D�� D��fD���D��D�C3D�\�D���D��fD��3D��fD�� D��D�C3D�s3D�vfDə�D̹�D���D���D�fD�@ D�y�D�s3D��D��D�� D���D�,�D�6fD�P D��3D�Y�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A&ffA1��A>ffAP  A`  Ap  A�  A�  A�ffA�ffA�  A���A���A�33A�  A�33A�33A�  A���A�  A�33A�ffB   B  B��B33B��B  BffB33B��B#��B(  B,ffB/��B333B8ffB=33B@  BE33BH  B^  B�  B�  B�33B�ffB�ffB�33B�33C��C� CL�C%L�C/L�C9ffCC��CM�3CW�3Ca��CkffCuL�C��C��3C���C���C���C���C��fC���C��3C��3C�� C���C��3C��3Cř�CʦfCϳ3Cԙ�C٦fC���C�� C�� C�ٚC���C�� D @ Dl�D�3D�3D&fDy�D%� D+��D2�D8` D>� DE�DK,�DQY�DW��D^3Dd33Djy�Dp�3Dv��D}9�D���D���D���D��D�9�D�` D�� D���D�� D��fD���D��D�C3D�\�D���D��fD��3D��fD�� D��D�C3D�s3D�vfDə�D̹�D���D���D�fD�@ D�y�D�s3D��D��D�� D���D�,�D�6fD�P D��3D�Y�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����J��J��J��-��Mӿ�-��-��Mӿ�-��-��-��n���Mӿ�Mӿ�������vɿ�V��I���X��ff���T�ӕ���-��&�ѩ�������o�� ſ��#��l�����`B�����w���!���m��~���1'���#��9X��t��O��?\j@�@6v�@H�9@EO�@=�@17L@"�H@@�@��?�dZ?�t�?��?�Ĝ?��H?���?�b?�K�?� �?vȴ?l1?Z�?I��?<�?1�?*=q?"�\?/?�P?&�?�h?�>�!>�l�>ۥ�>�1'>�j>�K�>��>�S�>��>�n�>t�j>9X>o=�1=]/<t��t��8Q콇+��F�hs�%�T�;dZ�I�^�Y��k��~�۾��9��V��b�������羳�F��^5�Õ������������1���j��p��Mӿ�/�C����hs��F��P����/��w�"��%��')xտ+��-O߿.��2-�49X�6E��8b�9��:^5�;��=�-�?;d�@A��@Ĝ�A%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��J��J��J��-��Mӿ�-��-��Mӿ�-��-��-��n���Mӿ�Mӿ�������vɿ�V��I���X��ff���T�ӕ���-��&�ѩ�������o�� ſ��#��l�����`B�����w���!���m��~���1'���#��9X��t��O��?\j@�@6v�@H�9@EO�@=�@17L@"�H@@�@��?�dZ?�t�?��?�Ĝ?��H?���?�b?�K�?� �?vȴ?l1?Z�?I��?<�?1�?*=q?"�\?/?�P?&�?�h?�>�!>�l�>ۥ�>�1'>�j>�K�>��>�S�>��>�n�>t�j>9X>o=�1=]/<t��t��8Q콇+��F�hs�%�T�;dZ�I�^�Y��k��~�۾��9��V��b�������羳�F��^5�Õ������������1���j��p��Mӿ�/�C����hs��F��P����/��w�"��%��')xտ+��-O߿.��2-�49X�6E��8b�9��:^5�;��=�-�?;d�@A��@Ĝ�A%111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�JB	�JB	�JB	�JB	�JB	�DB	�JB	�JB	�JB	�JB	�JB	�PB	�PB	�PB	�bB	�oB	�{B	�{B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�9B	�-B	�B	�3B	�dB	ȴB	��B	��B	��B	�)B	�HB	�yB
+B
�wBk�B�HBoB�BuB	7B��B�B�B�sB��B��B��BŢB�wB�3B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B	�JB	�JB	�JB	�JB	�JB	�DB	�JB	�JB	�JB	�JB	�JB	�PB	�PB	�PB	�bB	�oB	�{B	�{B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�9B	�-B	�B	�3B	�dB	ȴB	��B	��B	��B	�)B	�HB	�yB
+B
�wBk�B�HBoB�BuB	7B��B�B�B�sB��B��B��BŢB�wB�3B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <4` <4` <4` <4� <4@ <4` <4@ <4  <4` <4` <4` <4� <4@ <4@ <4@ <4` <4  <4� <4  <4  <4` <4� <4� <4� <4� <4` <4� <4@ <4� <4� <4� <4� <4� <4� <4� <4� <4� <4� <4� <5  <4� <5  <5  <6� <7� <8� <8� <8� <8� <8� <8� <8� <8� <8` <8@ <8` <8  <8@ <8  <7� <7� <8  <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7� <7` <7� <7� <7� <7` <7� <7� <7� <7� <7� <7� <7� <7` <7� <7@ <7` <7� <7` <7� <7` <7� <7` <7@ <7@ <7@ <7� <7� <7� <7` <7` <7` <7@ <7` <7� <7� <7� <7� <7` <7@ <7@ <7  <7� <7` <7@ <7� <7` <7` <7� <7` <7` <7@ <7@ <7@ <7� <7@ <7@ <7@ <7` <7@ <7` <7` PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                                                    202005120939182021100509093220211005090932  IF  ARFMCODA034d                                                                20200511145245                      G�O�G�O�G�O�                IF  ARGQCOQC4.5                                                                 20200511145423  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.5                                                                 20200511145423  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200512093918  IP  PSAL            A&ffD�Y�G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211005090932  IP  PSAL            A&ffD�Y�G�O�                