CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2019-09-14T15:28:32Z creation; 2019-09-14T15:28:50Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_030f      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        `  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  =<   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        `  =�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  @4   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     `  @�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  C,   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  E�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  F$   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  H�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  I   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  K|   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  Nt   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  P�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  Ql   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ](   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ],   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ]0   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ]4   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ]8   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ]x   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ]�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ]�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ]�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ]�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ]�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ]�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  S�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    S�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    V�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Y�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  \�Argo profile    3.1 1.2 19500101000000  20190914152832  20211006133609  6902729 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               �A   IF                                  2C  D   ARVOR                           AI2600-16FR312                  5900A04                         844 @��i>��?1   @��i>��?@To�4�[.@3���4�8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 25 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                      @�33@�  @�  @ٙ�@���A��A!��A1��AA��AQ��Aa��AnffA���A���A�  A�  A�  A�33A�33A�  A���A���A�33A���A���A�  A�33A�  A�33B33B33B33B��B��B33B��B��B$ffB)33B,��B/��B4ffB8  B;��B@  BC��BG33B\��B�33B�  B���B�33B�ffB�ffB�  CffC  C  C%�C/33C9��CC��CMffCWffCa� CkffCuL�CffC�� C���C�� C��3C��3C�� C�� C���C��3C��fC��3C�� C�� C�� Cʙ�Cϙ�CԦfCٳ3C���C�� C�fC�� C�� C�� D @ D� D�fD��DFfDs3D%�3D,  D2,�D8` D>�fDE�DK,�DQl�DW� D]��DdS3Djy�Dp�3Dw  D}33D��3D���D�3D��D�6fD�i�D�|�D���D��fD��fD��fD�&fD�<�D�\�D�y�D�� D��3D�� D��fD�  D�6fD�c3D�|�Də�D̼�D�� D�fD��D�9�D�Y�D�vfD�3D�fD�ٚD�  D��D�6fD�|�D�  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�33@�  @�  @ٙ�@���A��A!��A1��AA��AQ��Aa��AnffA���A���A�  A�  A�  A�33A�33A�  A���A���A�33A���A���A�  A�33A�  A�33B33B33B33B��B��B33B��B��B$ffB)33B,��B/��B4ffB8  B;��B@  BC��BG33B\��B�33B�  B���B�33B�ffB�ffB�  CffC  C  C%�C/33C9��CC��CMffCWffCa� CkffCuL�CffC�� C���C�� C��3C��3C�� C�� C���C��3C��fC��3C�� C�� C�� Cʙ�Cϙ�CԦfCٳ3C���C�� C�fC�� C�� C�� D @ D� D�fD��DFfDs3D%�3D,  D2,�D8` D>�fDE�DK,�DQl�DW� D]��DdS3Djy�Dp�3Dw  D}33D��3D���D�3D��D�6fD�i�D�|�D���D��fD��fD��fD�&fD�<�D�\�D�y�D�� D��3D�� D��fD�  D�6fD�c3D�|�Də�D̼�D�� D�fD��D�9�D�Y�D�vfD�3D�fD�ٚD�  D��D�6fD�|�D�  11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����bN��A���A��� ſ�  �ϝ��Ͼw�θR�˥㿫"Ѿ�l��-V>�C�?[��?X�u?;��>��=�%��  �X���F������푿��Ϳ��7��vɿ��+��|��p����#�ٺ^��"ѿ�
=�ٺ^����ߝ���G���-��n���Mӿ����
=�щ7��
=�щ7��~���hs���#�kƨ��w?I�^?�(�?�z�@1@2M�@E@B~�@C"�@;33@,�@0 �@1G�@.$�@0�9@01'@.@+��@+C�@( �@$�@"J@{@~�@�w@�-@��@��@t�@	��@�y@9X@ b?�/?��9?��?�j?�K�?��?���?�r�?щ7?��?�
=?���?�Q�?�G�?�  ?_�w?I��?,��?r�>�~�>���>n��>%�T=��=aG�<T����h��+�ȴ9����F��e`B�~�۾�7L������þ�9X��Q�Ǯ���`��(���ff��1��-���H�G���/�1'�����녿�Ͽ�ٿ"ѿ�ۿ Ĝ�"��$���&ff�&��'��)xտ+C��,�D�.{�/\)�0bN�2-�3t��3��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111��bN��A���A��� ſ�  �ϝ��Ͼw�θR�˥㿫"Ѿ�l��-V>�C�?[��?X�u?;��>��=�%��  �X���F������푿��Ϳ��7��vɿ��+��|��p����#�ٺ^��"ѿ�
=�ٺ^����ߝ���G���-��n���Mӿ����
=�щ7��
=�щ7��~���hs���#�kƨ��w?I�^?�(�?�z�@1@2M�@E@B~�@C"�@;33@,�@0 �@1G�@.$�@0�9@01'@.@+��@+C�@( �@$�@"J@{@~�@�w@�-@��@��@t�@	��@�y@9X@ b?�/?��9?��?�j?�K�?��?���?�r�?щ7?��?�
=?���?�Q�?�G�?�  ?_�w?I��?,��?r�>�~�>���>n��>%�T=��=aG�<T����h��+�ȴ9����F��e`B�~�۾�7L������þ�9X��Q�Ǯ���`��(���ff��1��-���H�G���/�1'�����녿�Ͽ�ٿ"ѿ�ۿ Ĝ�"��$���&ff�&��'��)xտ+C��,�D�.{�/\)�0bN�2-�3t��3��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oA��uA���A���A���A���A��A��A��RA��TA��\B5?B�B�B
=BA�B{�B�B#�B(�B@�BD�BO�B_;Bp�Bx�B�B�B�1B�B�7B�bB��B��B��B��B��B�B�!B�-B�9B�FB�qB�wBBǮB��B��B�B	F�B	�?B
1'B
��B
��B)�B��BŢBBǮBĜB�B�wBĜBŢB��B��B��B��B�B�B��B��B��B��B��B��BǮBÖBĜBĜBĜBĜBĜBÖBÖB�wB�}B�}B�}B�}B��B�qB�wB�}B�qB�dB�9B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A��uA���A���A���A���A��A��A��RA��TA��\B5?B�B�B
=BA�B{�B�B#�B(�B@�BD�BO�B_;Bp�Bx�B�B�B�1B�B�7B�bB��B��B��B��B��B�B�!B�-B�9B�FB�qB�wBBǮB��B��B�B	F�B	�?B
1'B
��B
��B)�B��BŢBBǮBĜB�B�wBĜBŢB��B��B��B��B�B�B��B��B��B��B��B��BǮBÖBĜBĜBĜBĜBĜBÖBÖB�wB�}B�}B�}B�}B��B�qB�wB�}B�qB�dB�9B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                                                    202004090851152021100613360920211006133609  IF  ARFMCODA030f                                                                20190914152832                      G�O�G�O�G�O�                IF  ARGQCOQC4.3                                                                 20190914152850  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC4.3                                                                 20190914152850  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200409085115  IP  PSAL            @�33D�  G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211006133609  IP  PSAL            @�33D�  G�O�                