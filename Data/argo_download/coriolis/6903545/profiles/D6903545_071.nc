CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:30Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  Cd   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  J8   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  K�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  R�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Y�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  [P   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  b$   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  c�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  j�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  q�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  s<   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  z   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  {�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �H   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �X   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �\   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �l   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �p   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �t   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �x   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200829092230  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               GA   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @�򈈈�1   @�򈈈�@P��US�8ɸ����8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A6ffAH  AY��Ak33A{33A�ffA�ffA�ffA�ffA���A���A�  A���A���A�33Aՙ�A�ffA�  A�ffA�  B ffB��B33B
��B��B��B��B33B��B$  B(ffB-33B1��B6ffB;��B@ffBE��BJ��BP  BU33BZ��B`ffBe��Bk��Bq33Bv��B|��B�ffB�ffB�ffB�ffB�ffB�ffB�ffB���B���B���B�  B�33B�ffB���B�  B�33B�ffB���B�  B�33B���B�  B�ffB���B�33B���B���B�33Bߙ�B�33B晚B�  B홚B���B�33B���B�  CL�C�C�fC�3C	ffC  C�fC��CffC33C�C�fC�3CffC33C   C!��C#��C%ffC'33C)�C*�fC,�3C/ffC233C4  C5��C7��C9� C;L�C=�C>�fC@�3CC� CF33CH  CI��CK��CMffCO33CP�fCR�3CUffCX�CY�fC[��C]ffC_33Ca  Cb��Ce� Ch33Cj  Ck�fCm�3Co� CqffCs33Cu�Cv�fCx�3C{ffC~�C��C���C��3C���C�� C�ffC��3C��C�  C�ٚC��3C���C�s3C�� C��C��fC�� C���C�s3C�� C��C��fC�� C���C�s3C��3C��3C���C��fC�s3C��3C��3C���C��fC�� C�� C��C�ٚC��3C���C�ffC��fC��3C���C��fC�� C�� C��C��3C���C��fC�� C�ffC��fC��3C���C��3C���C�s3C��3C��C��3C���C��3C���C�s3C�� C��C�  C��fC���CƳ3CǙ�CȀ C�ffC�� C��C�  C��fC���Cϳ3CЙ�Cь�C�s3C�Y�C�� C�&fC��C�  C��fC���C�� CۦfCܙ�C݌�C�s3C�ffC�L�C�3C�&fC��C��C��3C��fC�ٚC���C�� C�3C�fC왚C��C� C�s3C�s3C�ffC�Y�C�L�C��3C�&fC��C��C��C��C�� C��fD 33Dy�D��D  D@ D�fD��D�3D
  Dl�D� D�D3DffD�3DfDS3DffD��DfD3D` D��DfDY�Dl�D � D"�D#,�D$� D%ٚD&��D(@ D)�3D*�fD+��D-L�D.` D/��D13D2&fD3y�D4�3D5�fD7@ D8� D9�3D:��D<33D=��D>�3D?�3DA33DB��DC�3DD�3DF@ DG��DH�3DI�3DK3DL� DM��DO�DP,�DQS3DR� DT,�DUS3DVy�DW��DX� DZ33D[� D\�fD]��D_�D`y�Da��Dc3Dd9�De` Df�fDg��Di` Dj�fDk��Dl�3Dn@ Do��Dp�3Dq��Ds�Dt�fDu��Dw3Dx33DyY�Dz� D{�3D}Y�D~� D�fD�ffD��D�� D�c3D��fD���D�C3D��fD���D�  D��3D�I�D���D���D�S3D��D�� D�fD�� D�FfD�� D���D�VfD�� D��3D��D�� D�C3D�  D���D�S3D��fD�|�D�fD���D�C3D�  D���D�S3D���D��3D��D�� D�I�D�� D���D�\�D��fD�� D�&fD��3D�\�D��fD���D�&fD��3D�\�D��D���D�p D��D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A6ffAH  AY��Ak33A{33A�ffA�ffA�ffA�ffA���A���A�  A���A���A�33Aՙ�A�ffA�  A�ffA�  B ffB��B33B
��B��B��B��B33B��B$  B(ffB-33B1��B6ffB;��B@ffBE��BJ��BP  BU33BZ��B`ffBe��Bk��Bq33Bv��B|��B�ffB�ffB�ffB�ffB�ffB�ffB�ffB���B���B���B�  B�33B�ffB���B�  B�33B�ffB���B�  B�33B���B�  B�ffB���B�33B���B���B�33Bߙ�B�33B晚B�  B홚B���B�33B���B�  CL�C�C�fC�3C	ffC  C�fC��CffC33C�C�fC�3CffC33C   C!��C#��C%ffC'33C)�C*�fC,�3C/ffC233C4  C5��C7��C9� C;L�C=�C>�fC@�3CC� CF33CH  CI��CK��CMffCO33CP�fCR�3CUffCX�CY�fC[��C]ffC_33Ca  Cb��Ce� Ch33Cj  Ck�fCm�3Co� CqffCs33Cu�Cv�fCx�3C{ffC~�C��C���C��3C���C�� C�ffC��3C��C�  C�ٚC��3C���C�s3C�� C��C��fC�� C���C�s3C�� C��C��fC�� C���C�s3C��3C��3C���C��fC�s3C��3C��3C���C��fC�� C�� C��C�ٚC��3C���C�ffC��fC��3C���C��fC�� C�� C��C��3C���C��fC�� C�ffC��fC��3C���C��3C���C�s3C��3C��C��3C���C��3C���C�s3C�� C��C�  C��fC���CƳ3CǙ�CȀ C�ffC�� C��C�  C��fC���Cϳ3CЙ�Cь�C�s3C�Y�C�� C�&fC��C�  C��fC���C�� CۦfCܙ�C݌�C�s3C�ffC�L�C�3C�&fC��C��C��3C��fC�ٚC���C�� C�3C�fC왚C��C� C�s3C�s3C�ffC�Y�C�L�C��3C�&fC��C��C��C��C�� C��fD 33Dy�D��D  D@ D�fD��D�3D
  Dl�D� D�D3DffD�3DfDS3DffD��DfD3D` D��DfDY�Dl�D � D"�D#,�D$� D%ٚD&��D(@ D)�3D*�fD+��D-L�D.` D/��D13D2&fD3y�D4�3D5�fD7@ D8� D9�3D:��D<33D=��D>�3D?�3DA33DB��DC�3DD�3DF@ DG��DH�3DI�3DK3DL� DM��DO�DP,�DQS3DR� DT,�DUS3DVy�DW��DX� DZ33D[� D\�fD]��D_�D`y�Da��Dc3Dd9�De` Df�fDg��Di` Dj�fDk��Dl�3Dn@ Do��Dp�3Dq��Ds�Dt�fDu��Dw3Dx33DyY�Dz� D{�3D}Y�D~� D�fD�ffD��D�� D�c3D��fD���D�C3D��fD���D�  D��3D�I�D���D���D�S3D��D�� D�fD�� D�FfD�� D���D�VfD�� D��3D��D�� D�C3D�  D���D�S3D��fD�|�D�fD���D�C3D�  D���D�S3D���D��3D��D�� D�I�D�� D���D�\�D��fD�� D�&fD��3D�\�D��fD���D�&fD��3D�\�D��D���D�p D��D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����翪������=q��^5����������"ѿ��������㿬���vɿ�녿�����&鿱hs��$ݿ�����p���%�����H���#��^5��^5�����ə��Ƨ��$ݿ��T����ȴ��Kǿ�l���l���Kǿ�
=�Ƨ�Ƈ+�Ł���`��  ������V��/����^5��xտ������  ���w���D��C���r������t���vɿ�{���-��j��dZ��dZ���m������H������񪿏|�ۿ�"ѿ��^���}�|��z^5�k��^�ۿZ��@Ĝ�6�+�5?}�3t��0 ſ-���w�^5��`����n������1��S�������˾t�j�["Ѿ@��<j�2-���J���������+�aG���9X<t�=8Q�=��=�+=�%=Y�=8Q�=�j>J>�>5?}>:^5>Kƨ>^5?>l�D>|�>�%>��\>�7L>��`>��>�(�>�Z>��h>���>���>�Q�>\>ɺ^>�ƨ>��>�\)>��>ܬ>ݲ->���>��y>�V>�h>��>��>��P>��>�z�>�bN>t�=ix�=�{=���=�
==�9X=���=�O�=��=m�h=Y�=�O�>�>�>$�=��=��T=ix�<��
<���<u�t��#�
�ě�=#�
=H�9=\)=<j=T����t���O߽�-����Ƨ���`����"ѽ���h�����9X����<j��`B<�1=T��>�>6E�>H�9>�"�>��
>��+>��>��>�ff>���>��>�|�>�?}>��/>���?�`?�? �>��H>�p�>��#>�>�;d>��T>���>��>�n�>�%>�v�>�t�>��T>�h>�(�>�Q�>���>�K�>��m?M�?��>�|�>�dZ>�X>�->>�l�>ܬ>�z�>��>��>�^5>�Q�>�9X>��j>��F>��h>�ff>�Z>�Z>��
>��
>���>���>���>���>�Ĝ>�G�>���>���>�M�>��R>��>�\)>��>���>���>��T>���>�>�>�ff>�ff>�`B>��
>��
>�G�>���>�M�>���>�|�>Ǯ>�7L>���>�5?>�"�>��>��>��D>�V>�1>���>��F>�9X>�p�>���>�C�>�ƨ>�$�>��>�v�>�j>�E�>� �>�~�>��y>�Z>�Ĝ>�5?>���>�b>���>���>��>�=q>��>~��>z�H>{�m>ix�>`A�>aG�>`A�>W
=>T��>P�`>I�^><j>1&�>�w>�>I�>   =�"�=ě�=��
=��=u=m�h=D��=��=t�<�/<�1<�o<49X;��
    �o�o�#�
���ͼ�`B���0 ŽD���P�`�u��+��\)���w���T�� Ž�j������xս��o�z���(�þ+�2-�7KǾ?|�I�^�T���\(��hr��l�D�l�D�u�~�۾�%�����7L��I���\)��񪾖�+����������r����羮{������&龵���#��푾����J�Õ������$ݾ������;�bN�������5?��;d�������/��~��� ž�9X��E����پ�KǾ��پ�����dZ��j� A��ff�r��r��r���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ���翪������=q��^5����������"ѿ��������㿬���vɿ�녿�����&鿱hs��$ݿ�����p���%�����H���#��^5��^5�����ə��Ƨ��$ݿ��T����ȴ��Kǿ�l���l���Kǿ�
=�Ƨ�Ƈ+�Ł���`��  ������V��/����^5��xտ������  ���w���D��C���r������t���vɿ�{���-��j��dZ��dZ���m������H������񪿏|�ۿ�"ѿ��^���}�|��z^5�k��^�ۿZ��@Ĝ�6�+�5?}�3t��0 ſ-���w�^5��`����n������1��S�������˾t�j�["Ѿ@��<j�2-���J���������+�aG���9X<t�=8Q�=��=�+=�%=Y�=8Q�=�j>J>�>5?}>:^5>Kƨ>^5?>l�D>|�>�%>��\>�7L>��`>��>�(�>�Z>��h>���>���>�Q�>\>ɺ^>�ƨ>��>�\)>��>ܬ>ݲ->���>��y>�V>�h>��>��>��P>��>�z�>�bN>t�=ix�=�{=���=�
==�9X=���=�O�=��=m�h=Y�=�O�>�>�>$�=��=��T=ix�<��
<���<u�t��#�
�ě�=#�
=H�9=\)=<j=T����t���O߽�-����Ƨ���`����"ѽ���h�����9X����<j��`B<�1=T��>�>6E�>H�9>�"�>��
>��+>��>��>�ff>���>��>�|�>�?}>��/>���?�`?�? �>��H>�p�>��#>�>�;d>��T>���>��>�n�>�%>�v�>�t�>��T>�h>�(�>�Q�>���>�K�>��m?M�?��>�|�>�dZ>�X>�->>�l�>ܬ>�z�>��>��>�^5>�Q�>�9X>��j>��F>��h>�ff>�Z>�Z>��
>��
>���>���>���>���>�Ĝ>�G�>���>���>�M�>��R>��>�\)>��>���>���>��T>���>�>�>�ff>�ff>�`B>��
>��
>�G�>���>�M�>���>�|�>Ǯ>�7L>���>�5?>�"�>��>��>��D>�V>�1>���>��F>�9X>�p�>���>�C�>�ƨ>�$�>��>�v�>�j>�E�>� �>�~�>��y>�Z>�Ĝ>�5?>���>�b>���>���>��>�=q>��>~��>z�H>{�m>ix�>`A�>aG�>`A�>W
=>T��>P�`>I�^><j>1&�>�w>�>I�>   =�"�=ě�=��
=��=u=m�h=D��=��=t�<�/<�1<�o<49X;��
    �o�o�#�
���ͼ�`B���0 ŽD���P�`�u��+��\)���w���T�� Ž�j������xս��o�z���(�þ+�2-�7KǾ?|�I�^�T���\(��hr��l�D�l�D�u�~�۾�%�����7L��I���\)��񪾖�+����������r����羮{������&龵���#��푾����J�Õ������$ݾ������;�bN�������5?��;d�������/��~��� ž�9X��E����پ�KǾ��پ�����dZ��j� A��ff�r��r��r���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB �B"�B!�B!�B"�B$�B"�B&�B&�B%�B/B;dBM�B��B��B�BgmB(�B@�Bo�B�+B�wBDB�B49BH�BZB}�B��B�FB��BŢB��B�;B�B��B��B��B	  B	B	
=B	�B	�B	!�B	-B	7LB	=qB	?}B	@�B	@�B	<jB	K�B	L�B	Q�B	R�B	W
B	[#B	^5B	e`B	bNB	ffB	ffB	hsB	jB	p�B	r�B	w�B	z�B	�B	�7B	�=B	�\B	��B	��B	�B	�!B	�-B	�qB	��B	��B	��B	�5B	�BB	�ZB	�HB	�sB	�B	��B
B
JB
#�B
'�B
(�B
6FB
;dB
:^B
I�B
Q�B
ZB
`BB
e`B
m�B
v�B
z�B
y�B
� B
�%B
�PB
�{B
��B
��B
��B
��B
��B
��B
��B
�FB
�dB
��B
B
�dB
ǮB
��B
��B
��B
�
B
�B
�#B
�BB
�ZB
�TB
�B
�B
�B
�B
��B
��B
��B
��B
��B
��B  B  BBBB
��B
��BB
�B
�B
�B
��B
��B
�sB
��B
�B
�TB
��B
��B
��B
��B
��B
��B
��B
��BBBBBBBBBB%BB1BDBJBPB1BPB	7B1B+B+B1B+B1B1B	7BJBPBoBoBuB�B"�B)�B+B�B5?B2-B.B9XB7LBG�BB�BH�BN�BQ�BQ�BO�BW
BT�BS�BR�BP�BS�BZBVBVBVBXBVBZBW
BN�BZBcTBm�BaHBiyBm�Bl�Bm�Bm�Bo�Bm�Bl�Bn�Bm�Bo�Bm�Bl�Bl�Bl�Bl�Bk�Bm�Bl�Bp�Bn�Bn�Bm�Bm�Bl�Bm�Bn�Bn�Bn�Bm�Bn�Bo�Bp�Bq�Bq�Br�Br�Bt�Br�Bx�Bx�Bx�B|�B|�B~�B~�B� B~�B� B�B�B�B�B�+B�=B�\B�VB�oB��B��B�{B�oB�hB�bB�hB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B�B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B �B"�B!�B!�B"�B$�B"�B&�B&�B%�B/B;dBM�B��B��B�BgmB(�B@�Bo�B�+B�wBDB�B49BH�BZB}�B��B�FB��BŢB��B�;B�B��B��B��B	  B	B	
=B	�B	�B	!�B	-B	7LB	=qB	?}B	@�B	@�B	<jB	K�B	L�B	Q�B	R�B	W
B	[#B	^5B	e`B	bNB	ffB	ffB	hsB	jB	p�B	r�B	w�B	z�B	�B	�7B	�=B	�\B	��B	��B	�B	�!B	�-B	�qB	��B	��B	��B	�5B	�BB	�ZB	�HB	�sB	�B	��B
B
JB
#�B
'�B
(�B
6FB
;dB
:^B
I�B
Q�B
ZB
`BB
e`B
m�B
v�B
z�B
y�B
� B
�%B
�PB
�{B
��B
��B
��B
��B
��B
��B
��B
�FB
�dB
��B
B
�dB
ǮB
��B
��B
��B
�
B
�B
�#B
�BB
�ZB
�TB
�B
�B
�B
�B
��B
��B
��B
��B
��B
��B  B  BBBB
��B
��BB
�B
�B
�B
��B
��B
�sB
��B
�B
�TB
��B
��B
��B
��B
��B
��B
��B
��BBBBBBBBBB%BB1BDBJBPB1BPB	7B1B+B+B1B+B1B1B	7BJBPBoBoBuB�B"�B)�B+B�B5?B2-B.B9XB7LBG�BB�BH�BN�BQ�BQ�BO�BW
BT�BS�BR�BP�BS�BZBVBVBVBXBVBZBW
BN�BZBcTBm�BaHBiyBm�Bl�Bm�Bm�Bo�Bm�Bl�Bn�Bm�Bo�Bm�Bl�Bl�Bl�Bl�Bk�Bm�Bl�Bp�Bn�Bn�Bm�Bm�Bl�Bm�Bn�Bn�Bn�Bm�Bn�Bo�Bp�Bq�Bq�Br�Br�Bt�Br�Bx�Bx�Bx�B|�B|�B~�B~�B� B~�B� B�B�B�B�B�+B�=B�\B�VB�oB��B��B�{B�oB�hB�bB�hB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B�B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092230                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092337  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092337  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134655  IP  PSAL            A6ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172541  IP  PSAL            A6ffD��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A6ffD��3G�O�                