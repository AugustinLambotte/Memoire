CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:29Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z           :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z           C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  J�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���        L�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        S�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Z�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        \�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  c�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        e�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        l�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        u�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  |�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        ~�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �T   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �d   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �h   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �x   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �|   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200829092229  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               <A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��w��O�1   @��w��O�@Q���ᛸ�1�&U{h8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A.ffA@  AQ��Ac33AnffA���A�33A���A�  A�  A���A�  A�33A�ffA�  Aљ�A�ffA�33A�  A�ffA�33B   B��B33B33B  B33B  BffB   B#33B'��B+33B0ffB2��B6��B;��B?��BDffBI33BL  BP  BR��BW33B[33B_��Bc��Bg33BlffBp  Bs��Bw��B{��B33B�  B���B���B�ffB�  B���B���B���B���B���B���B���B�  B�ffB���B���B�  B�  B�33B���B�  B���B���B�33B�ffB���B���B�ffB�  B�ffB�  B���B�ffB���B�ffB�  Bʙ�BΙ�B���B�ffB���B�  B���B�ffB���B�ffB�33B���B���B���C33C33C33C�C	�C33C33C33CL�CL�CL�CffCL�CL�CffCffC!� C#��C%�3C'�3C)�3C+��C-��C/�fC1� C3ffC5� C7  C9�C;L�C=33C?33CA33CC33CE�CG�CI  CKffCM�fCO�fCQ�fCS��CU��CW��CY�3C[��C]� C_L�Ca33Cc�CeffCg��Ci�3Ck��Cm� CoffCqL�Cs33Cu33Cw�Cy�C{  C}ffC�fC��3C��3C��fC��fC��fC��fC��fC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C�� C���C�� C��3C���C���C�� C��3C��3C��3C��fC��fC���C���C���C���C���C���C�� C�� C�� C���C���C���C���C���C���C���C���C���C���C��fC��fC��fC��fC��fC��fC��fC��3C��3C���C��fC��3C�� C���C��fC�� C���C�� C���Có3Cę�Cų3C�ٚC�� CȦfCɌ�CʦfC�� C̳3C͙�C�� C�ٚC���Cѳ3CҦfCӌ�CԀ Cճ3C�ٚC���C�� Cٳ3Cڳ3CۦfCܙ�C݌�Cހ C߳3C�3CᙚC�� C�3C�3C�3C�3C�3C�3C�3C�� C���C왚C��fC�3C�� C�ٚC�3C�� C�ٚC��3C���C��fC�� C��fC�� C�s3C���D 9�Ds3D�fD��DL�D� D��D��D
@ Ds3D�fD�3DFfDy�D��D  D33Dl�D�fDfD,�Dl�D�3D  D@ Dl�D ��D"fD#@ D$y�D%�3D&�3D(,�D)l�D*� D+��D-,�D.ffD/�fD0�fD2&fD3l�D4� D5��D733D8s3D9��D;  D<,�D=s3D>�fD?��DA33DBs3DC��DD��DF33DGs3DH��DJ  DK@ DL�fDM�fDO�DPS3DQy�DR�fDS�3DU@ DV�3DWٚDYfDZ@ D[y�D\�3D]�fD_9�D`�3Da��Db� Dd&fDey�Df��Dg��Di,�Djs3Dk� Dm�Dn9�Dos3Dp�3Dq��Ds@ DtffDu�3Dv��Dx,�Dyy�Dz��D|  D}9�D~y�D�3D�y�D�fD��fD�Y�D���D���D�C3D��D�� D��D���D�Y�D��fD���D�9�D�ٚD�y�D��D�� D�c3D�fD���D�33D�ٚD�� D��D��fD�` D���D���D�<�D���D�� D�fD���D�` D���D��3D�<�D��fD��3D�#3D��3D�c3D�fD��3D�<�D�� D��3D��D��fD�S3D�� D�� D�33D��3D�vfD��D�� D�ffD�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A.ffA@  AQ��Ac33AnffA���A�33A���A�  A�  A���A�  A�33A�ffA�  Aљ�A�ffA�33A�  A�ffA�33B   B��B33B33B  B33B  BffB   B#33B'��B+33B0ffB2��B6��B;��B?��BDffBI33BL  BP  BR��BW33B[33B_��Bc��Bg33BlffBp  Bs��Bw��B{��B33B�  B���B���B�ffB�  B���B���B���B���B���B���B���B�  B�ffB���B���B�  B�  B�33B���B�  B���B���B�33B�ffB���B���B�ffB�  B�ffB�  B���B�ffB���B�ffB�  Bʙ�BΙ�B���B�ffB���B�  B���B�ffB���B�ffB�33B���B���B���C33C33C33C�C	�C33C33C33CL�CL�CL�CffCL�CL�CffCffC!� C#��C%�3C'�3C)�3C+��C-��C/�fC1� C3ffC5� C7  C9�C;L�C=33C?33CA33CC33CE�CG�CI  CKffCM�fCO�fCQ�fCS��CU��CW��CY�3C[��C]� C_L�Ca33Cc�CeffCg��Ci�3Ck��Cm� CoffCqL�Cs33Cu33Cw�Cy�C{  C}ffC�fC��3C��3C��fC��fC��fC��fC��fC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C�� C���C�� C��3C���C���C�� C��3C��3C��3C��fC��fC���C���C���C���C���C���C�� C�� C�� C���C���C���C���C���C���C���C���C���C���C��fC��fC��fC��fC��fC��fC��fC��3C��3C���C��fC��3C�� C���C��fC�� C���C�� C���Có3Cę�Cų3C�ٚC�� CȦfCɌ�CʦfC�� C̳3C͙�C�� C�ٚC���Cѳ3CҦfCӌ�CԀ Cճ3C�ٚC���C�� Cٳ3Cڳ3CۦfCܙ�C݌�Cހ C߳3C�3CᙚC�� C�3C�3C�3C�3C�3C�3C�3C�� C���C왚C��fC�3C�� C�ٚC�3C�� C�ٚC��3C���C��fC�� C��fC�� C�s3C���D 9�Ds3D�fD��DL�D� D��D��D
@ Ds3D�fD�3DFfDy�D��D  D33Dl�D�fDfD,�Dl�D�3D  D@ Dl�D ��D"fD#@ D$y�D%�3D&�3D(,�D)l�D*� D+��D-,�D.ffD/�fD0�fD2&fD3l�D4� D5��D733D8s3D9��D;  D<,�D=s3D>�fD?��DA33DBs3DC��DD��DF33DGs3DH��DJ  DK@ DL�fDM�fDO�DPS3DQy�DR�fDS�3DU@ DV�3DWٚDYfDZ@ D[y�D\�3D]�fD_9�D`�3Da��Db� Dd&fDey�Df��Dg��Di,�Djs3Dk� Dm�Dn9�Dos3Dp�3Dq��Ds@ DtffDu�3Dv��Dx,�Dyy�Dz��D|  D}9�D~y�D�3D�y�D�fD��fD�Y�D���D���D�C3D��D�� D��D���D�Y�D��fD���D�9�D�ٚD�y�D��D�� D�c3D�fD���D�33D�ٚD�� D��D��fD�` D���D���D�<�D���D�� D�fD���D�` D���D��3D�<�D��fD��3D�#3D��3D�c3D�fD��3D�<�D�� D��3D��D��fD�S3D�� D�� D�33D��3D�vfD��D�� D�ffD�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��� �� �� �� �� A�� A�� A��;d� A��   � A�� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ��w� ��� Ĝ�!G��"Mӿ$��$���$�/�$���#���#o�#S��#�
�#���#�
�$Z�%��%��%`B�$��$�/�%��%`B�%`B�%`B�#S��$��$Z�%��%`B�$�/�&$ݿ&ff�&ff�&ff�&ff���Kǿb��+����ȴ��P�33�hs����hs��1'�J�$ݿ����;���پ�bN����>49X>�9X?ff?5?}?�33?���?�?}?ӶF?��?���?��@ �u?�{?��?�dZ?�$�?��?�|�?�"�?�7L?���?�?�G�?߾w?��;?�;d?���?�/?�j?��m?�"�?�=q?�b?�`B?�M�?ϝ�?��?�{?̋D?�"�?ə�?��y?�ff?���?��?�\)?���?�p�?���?��j?��\?��?�A�?�;d?�V?�/?�C�?��y?��
?�%?� �?���?���?�
=?���?���?��?���?�Ĝ?��w?��R?��h?��?�X?�ȴ?�ff?��F?��\?�A�?�w?}/?v?kC�?gl�?h1'?g�?h1'?hr�?g�?f$�?c��?_;d?[�m?WK�?U�?Q��?L�D?L1?KC�?G+?Fff?E�T?D�/?DZ?C��?Co?@Ĝ?@A�??�w?=�-?;�m?;��?8��?3t�?0 �?.�?.V?.V?-�h?-V?,1?)�^?(1'?'l�?'l�?&�y?"J? �?��?�?/?dZ?�H?�#?X?�u?��?
=?ȴ?9X?�?�;?��?C�?
=q?��?�y?��? A�>��H>�X>��j>�x�>�Z>ڟ�>և+>�n�>��;>�I�>�ƨ>�1'>��>��>�^5>�->���>�l�>��/>�Ĝ>���>�n�>��`>�\)>�I�>��>��>|�>q��>o��>m�h>aG�>T��>S��>R�>O�;>Kƨ>J��>A�7>>v�>@�>@�>@�>>v�><j><j>>v�>;dZ>7K�>333>0 �>+>(��>'�>&�y>$�/>!��>��>z�>V>C�>�>J=���=��=�=��=�
==��`=Ƨ�=�Q�=�^5=�-=���=���=�O�=aG�=L��=,1=C�<�<�/<�9X<�C��o�e`B��`B���t��,1�49X�D����%��o��7L��hs���P���w��1��-��9X��Q�Ƨ�ȴ9���ͽ���
=�����xս�����C��O߾�u����R�$�/�+�-V�/��0 ž49X�7KǾ>vɾ@��D���G��I�^�J���J���L�;O�;�W
=�Xb�Z��["Ѿ]/�e`B�hr��l�D�n���q���r�!�t�j�|푾�J��o��o���������9��C���I�����hs��t���zᾖ�+���P�������㾜���/���R���w���w��G��������徣�
��Z��Z��`B���T��r���~����D���h�����&龲�!���F��9X��9X��ȴ���#��^5���H���m��  �\�Õ��Ƨ�ȴ9��C���\)���`��t��Ձ��
=�׍P�ؓu�������ۥ�ݲ-��;d��A���Ĝ������MӾ��T��������D��V��{���� ž�!��9X��Q���H��푿G����o�o���111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111� �� �� �� �� A�� A�� A��;d� A��   � A�� �� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ� Ĝ��w� ��� Ĝ�!G��"Mӿ$��$���$�/�$���#���#o�#S��#�
�#���#�
�$Z�%��%��%`B�$��$�/�%��%`B�%`B�%`B�#S��$��$Z�%��%`B�$�/�&$ݿ&ff�&ff�&ff�&ff���Kǿb��+����ȴ��P�33�hs����hs��1'�J�$�G�O�G�O����پ�bN����>49X>�9X?ff?5?}?�33?���?�?}?ӶF?��?���?��@ �u?�{?��?�dZ?�$�?��?�|�?�"�?�7L?���?�?�G�?߾w?��;?�;d?���?�/?�j?��m?�"�?�=q?�b?�`B?�M�?ϝ�?��?�{?̋D?�"�?ə�?��y?�ff?���?��?�\)?���?�p�?���?��j?��\?��?�A�?�;d?�V?�/?�C�?��y?��
?�%?� �?���?���?�
=?���?���?��?���?�Ĝ?��w?��R?��h?��?�X?�ȴ?�ff?��F?��\?�A�?�w?}/?v?kC�?gl�?h1'?g�?h1'?hr�?g�?f$�?c��?_;d?[�m?WK�?U�?Q��?L�D?L1?KC�?G+?Fff?E�T?D�/?DZ?C��?Co?@Ĝ?@A�??�w?=�-?;�m?;��?8��?3t�?0 �?.�?.V?.V?-�h?-V?,1?)�^?(1'?'l�?'l�?&�y?"J? �?��?�?/?dZ?�H?�#?X?�u?��?
=?ȴ?9X?�?�;?��?C�?
=q?��?�y?��? A�>��H>�X>��j>�x�>�Z>ڟ�>և+>�n�>��;>�I�>�ƨ>�1'>��>��>�^5>�->���>�l�>��/>�Ĝ>���>�n�>��`>�\)>�I�>��>��>|�>q��>o��>m�h>aG�>T��>S��>R�>O�;>Kƨ>J��>A�7>>v�>@�>@�>@�>>v�><j><j>>v�>;dZ>7K�>333>0 �>+>(��>'�>&�y>$�/>!��>��>z�>V>C�>�>J=���=��=�=��=�
==��`=Ƨ�=�Q�=�^5=�-=���=���=�O�=aG�=L��=,1=C�<�<�/<�9X<�C��o�e`B��`B���t��,1�49X�D����%��o��7L��hs���P���w��1��-��9X��Q�Ƨ�ȴ9���ͽ���
=�����xս�����C��O߾�u����R�$�/�+�-V�/��0 ž49X�7KǾ>vɾ@��D���G��I�^�J���J���L�;O�;�W
=�Xb�Z��["Ѿ]/�e`B�hr��l�D�n���q���r�!�t�j�|푾�J��o��o���������9��C���I�����hs��t���zᾖ�+���P�������㾜���/���R���w���w��G��������徣�
��Z��Z��`B���T��r���~����D���h�����&龲�!���F��9X��9X��ȴ���#��^5���H���m��  �\�Õ��Ƨ�ȴ9��C���\)���`��t��Ձ��
=�׍P�ؓu�������ۥ�ݲ-��;d��A���Ĝ������MӾ��T��������D��V��{���� ž�!��9X��Q���H��푿G����o�o���111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
2-B
2-B
2-B
1'B
1'B
2-B
2-B
1'B
2-B
49B
2-B
1'B
1'B
2-B
2-B
2-B
1'B
1'B
2-B
2-B
33B
33B
33B
+B
.B
/B
/B
.B
/B
/B
.B
/B
/B
/B
/B
.B
.B
.B
-B
.B
.B
0!B
.B
-B
.B
/B
-B
/B
-B
.B
/B
.B
.B
-B
.B
0!B
-B
9XB
8RB
5?B
9XB
33B
<jB
5?B
9XB
49B
A�B
:^B
D�B
B�B
S�B
D�B
uB
jB
y�B
��B
ÖB
B
�TB
��BB�BbNB��B��B��B�;B��BƨBƨBƨBǮBƨBĜBƨBBBB��BĜB��B��BBBB��B��BBBBÖBBBB��B��B��B��B�}B�qB�}B�}B�wB�wB�qB�dB�dB�dB�jB�^B�^B�^B�XB�RB�^B�FB�LB�FB�LB�LB�FB�FB�FB�?B�FB�?B�?B�9B�9B�3B�?B�9B�-B�9B�-B�9B�3B�-B�'B�!B�B�'B�-B�3B�-B�3B�-B�3B�'B�9B�3B�3B�3B�3B�3B�-B�3B�3B�3B�-B�-B�-B�9B�3B�-B�3B�9B�-B�-B�'B�-B�-B�3B�3B�-B�-B�3B�-B�-B�-B�3B�-B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�9B�9B�9B�9B�9B�?B�9B�9B�9B�9B�9B�?B�-B�3B�3B�-B�-B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B
2-B
2-B
2-B
1'B
1'B
2-B
2-B
1'B
2-B
49B
2-B
1'B
1'B
2-B
2-B
2-B
1'B
1'B
2-B
2-B
33B
33B
33B
+B
.B
/B
/B
.B
/B
/B
.B
/B
/B
/B
/B
.B
.B
.B
-B
.B
.B
0!B
.B
-B
.B
/B
-B
/B
-B
.B
/B
.B
.B
-B
.B
0!B
-B
9XB
8RB
5?B
9XB
33B
<jB
5?B
9XB
49B
A�B
:^B
D�B
B�B
S�G�O�G�O�B
jB
y�B
��B
ÖB
B
�TB
��BB�BbNB��B��B��B�;B��BƨBƨBƨBǮBƨBĜBƨBBBB��BĜB��B��BBBB��B��BBBBÖBBBB��B��B��B��B�}B�qB�}B�}B�wB�wB�qB�dB�dB�dB�jB�^B�^B�^B�XB�RB�^B�FB�LB�FB�LB�LB�FB�FB�FB�?B�FB�?B�?B�9B�9B�3B�?B�9B�-B�9B�-B�9B�3B�-B�'B�!B�B�'B�-B�3B�-B�3B�-B�3B�'B�9B�3B�3B�3B�3B�3B�-B�3B�3B�3B�-B�-B�-B�9B�3B�-B�3B�9B�-B�-B�'B�-B�-B�3B�3B�-B�-B�3B�-B�-B�-B�3B�-B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�9B�9B�9B�9B�9B�?B�9B�9B�9B�9B�9B�?B�-B�3B�3B�-B�-B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092229                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092325  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092325  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A.ffD�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172540  IP  PSAL            A.ffD�� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A.ffD�� G�O�                