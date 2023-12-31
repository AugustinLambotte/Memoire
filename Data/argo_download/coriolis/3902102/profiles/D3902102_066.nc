CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:24Z creation; 2020-11-17T12:19:03Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        �  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  B   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  Mp   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  \�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  f    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  h   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  o�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  w8   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  y    PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �8   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �h   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �h   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �h   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �hArgo profile    3.1 1.2 19500101000000  20200828144324  20220127170406  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               BA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @�V��>�1   @�V��>�@R�<��
�#�M=#S8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A333AA��ANffA^ffAl��A�  A�  A�33A�33A�33A�33A�33A�  A���A�ffA�33A���A�ffA�33A�  A���A�ffB33B33B��B��BffB  B  B��B#��B'��B,ffB0  B3��B7��B;��B@  BDffBH  BK��BO33BS33BW33B[33B_��Bc��Bh  Bk33Bo33Bs��Bx  B|��B��B�ffB���B�ffB�  B���B���B���B�ffB���B�ffB�ffB�ffB�ffB�ffB���B���B���B���B�  B�  B�33B�33B�ffB���B���B���B���B���B���B�  B�  B�  B�33B�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�33B�  B���B�ffB�33B�  B���B�ffC  CffC��CffC	L�C�C� C��C�3C��C� CffCL�C33C�C  C!ffC#��C%�3C'��C)� C+ffC-L�C/33C1�C3  C5ffC7�fC9��C;�3C=�3C?�3CA�3CC�3CE�3CG��CI�fCKffCM  CO  CQL�CSL�CUffCWffCYffC[� C]� C_��CaL�CcL�CeffCg� Ci� Ck33CmL�CoL�CqffCs� Cu� Cw��Cy33C{L�C}ffC� C���C��fC��3C�� C���C��3C�� C���C��fC��3C�� C���C��fC��3C�� C���C��fC�� C��fC�� C�ٚC�� C���C�� C��3C�ٚC�� C��fC���C���C��3C�ٚC���C�� C��3C��fC��fC���C���C���C��3C��3C��fC��fC�ٚC���C���C���C�� C�� C���C���C���C���C���C�ٚC�ٚC�ٚC�ٚC�ٚC��fC��fC��fC�ٚC�ٚC�ٚC���C�� Cĳ3Cų3CƦfCǦfCȦfCɦfCʦfC˦fC̳3Cͳ3Cγ3C�� C�� C�ٚCҦfCӳ3C���Cՙ�Cֳ3C�� Cؙ�C٦fC�� Cۙ�Cܳ3C�� C޳3Cߌ�C�3C�� C�3C㙚C�3C�ٚC�� C�fC虚C� C�3C���C�� C��3C�fCC���C�3C��3C��fC�ٚC���C�� C�� C��3C��fC�L�C�� D 9�D�3D�3D3DS3Dy�D� D��D
9�Dy�D� D3DL�Dy�D�fD�3DFfD� D��D��D9�D�fD�3D�fD9�D��D � D!�3D#&fD$` D%��D'3D(S3D)��D*�fD,fD-@ D.y�D/�3D0��D2  D3y�D4�3D6�D7FfD8�fD9� D;  D<@ D=�fD>�fD@�DAS3DBy�DC� DD�fDF33DGy�DH�fDJ3DK9�DLffDM��DN�fDP&fDQffDR�fDS��DU33DV� DW�fDY�DZ9�D[ffD\�3D]��D_@ D`�fDa��Db�3DdFfDey�Df��Dh  DiS3Dj� Dk�3Dl� Dn33Do�fDp��Dq�3Ds,�DtffDu�fDv� Dx  Dyy�Dz�3D|�D}L�D~��D�3D�|�D� D��3D�S3D��fD���D�FfD��fD��fD�&fD�ɚD�` D���D���D�<�D���D�y�D��D��fD�P D���D���D�I�D��D��fD��D�� D�S3D��3D���D�C3D��3D��fD�  D���D�c3D�  D���D�@ D�� D�� D�  D��3D�\�D��3D���D�FfD��fD��fD�&fD��fD�i�D���D��fD�9�D���D�� D�&fD�ɚD�\�D��3D���D�@ D��3D���D�  D��fD�\�D�3D��fD�<�D�� D�s3D��D�� D�\�D���D��fD�6fD��fD�y�D�#3D�� D�\�D���D���D�<�D�ٚD��3D��D���D�Vf111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A333AA��ANffA^ffAl��A�  A�  A�33A�33A�33A�33A�33A�  A���A�ffA�33A���A�ffA�33A�  A���A�ffB33B33B��B��BffB  B  B��B#��B'��B,ffB0  B3��B7��B;��B@  BDffBH  BK��BO33BS33BW33B[33B_��Bc��Bh  Bk33Bo33Bs��Bx  B|��B��B�ffB���B�ffB�  B���B���B���B�ffB���B�ffB�ffB�ffB�ffB�ffB���B���B���B���B�  B�  B�33B�33B�ffB���B���B���B���B���B���B�  B�  B�  B�33B�33B�33B�33B�ffB�ffB�ffB�ffB�ffB�33B�  B���B�ffB�33B�  B���B�ffC  CffC��CffC	L�C�C� C��C�3C��C� CffCL�C33C�C  C!ffC#��C%�3C'��C)� C+ffC-L�C/33C1�C3  C5ffC7�fC9��C;�3C=�3C?�3CA�3CC�3CE�3CG��CI�fCKffCM  CO  CQL�CSL�CUffCWffCYffC[� C]� C_��CaL�CcL�CeffCg� Ci� Ck33CmL�CoL�CqffCs� Cu� Cw��Cy33C{L�C}ffC� C���C��fC��3C�� C���C��3C�� C���C��fC��3C�� C���C��fC��3C�� C���C��fC�� C��fC�� C�ٚC�� C���C�� C��3C�ٚC�� C��fC���C���C��3C�ٚC���C�� C��3C��fC��fC���C���C���C��3C��3C��fC��fC�ٚC���C���C���C�� C�� C���C���C���C���C���C�ٚC�ٚC�ٚC�ٚC�ٚC��fC��fC��fC�ٚC�ٚC�ٚC���C�� Cĳ3Cų3CƦfCǦfCȦfCɦfCʦfC˦fC̳3Cͳ3Cγ3C�� C�� C�ٚCҦfCӳ3C���Cՙ�Cֳ3C�� Cؙ�C٦fC�� Cۙ�Cܳ3C�� C޳3Cߌ�C�3C�� C�3C㙚C�3C�ٚC�� C�fC虚C� C�3C���C�� C��3C�fCC���C�3C��3C��fC�ٚC���C�� C�� C��3C��fC�L�C�� D 9�D�3D�3D3DS3Dy�D� D��D
9�Dy�D� D3DL�Dy�D�fD�3DFfD� D��D��D9�D�fD�3D�fD9�D��D � D!�3D#&fD$` D%��D'3D(S3D)��D*�fD,fD-@ D.y�D/�3D0��D2  D3y�D4�3D6�D7FfD8�fD9� D;  D<@ D=�fD>�fD@�DAS3DBy�DC� DD�fDF33DGy�DH�fDJ3DK9�DLffDM��DN�fDP&fDQffDR�fDS��DU33DV� DW�fDY�DZ9�D[ffD\�3D]��D_@ D`�fDa��Db�3DdFfDey�Df��Dh  DiS3Dj� Dk�3Dl� Dn33Do�fDp��Dq�3Ds,�DtffDu�fDv� Dx  Dyy�Dz�3D|�D}L�D~��D�3D�|�D� D��3D�S3D��fD���D�FfD��fD��fD�&fD�ɚD�` D���D���D�<�D���D�y�D��D��fD�P D���D���D�I�D��D��fD��D�� D�S3D��3D���D�C3D��3D��fD�  D���D�c3D�  D���D�@ D�� D�� D�  D��3D�\�D��3D���D�FfD��fD��fD�&fD��fD�i�D���D��fD�9�D���D�� D�&fD�ɚD�\�D��3D���D�@ D��3D���D�  D��fD�\�D�3D��fD�<�D�� D�s3D��D�� D�\�D���D��fD�6fD��fD�y�D�#3D�� D�\�D���D���D�<�D�ٚD��3D��D���D�Vf111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����������������?}���������������?}�����z�䛦�䛦��Z��9X��9X������Ͽ�F��t�����!��n���녿�hs��Ĝ��bN���ۿ�1�ա˿ѩ����������E�������w��^5���翸�9�����ɺ^��?}�g��9���"��V��C����w�پ:^5��w�J�t��B�\��񪾔zᾗ
=�N��n��%�P�`=�E�>Xb>t�j>�
=>�^5?!%?-�h?0 �?P �?o��?��?���?�l�?�Q�?��?�$�?ȓu?�(�?�{?�?ڟ�?���?�/?Լj?��;?Ł?�ff?��?v�+?�bN?��H?�Z?���?���?�hs?�n�?�9X?���?�?}?�?}?���?��/?���?�o?���?��;?�;d?���?�p�?�V?��?�(�?�1?�1?��?�bN?�33?�?�
=?�+?�+?�K�?���?�Q�?���?��+?�`B?�bN?�~�?���?�{?�~�?�1'?�dZ?�I�?�v�?��;?��R?��?��`?�G�?�&�?�bN?���?��R?��?�V?��h?�/?��?��?���?�b?�K�?��!?}�-?|j?x�u?~5??�S�?���?~v�?s�F?t9X?r�!?r�!?qhs?p�`?r�!?s�F?st�?s��?v�+?w�P?{�m?}�-?��?��`?�%?��?�  ?�w?�M�?�33?���?� �?v?st�?q�?rn�?q��?n{?e�T?bJ?a%?\(�?Z�?Z^5?Z�?Z�H?["�?[dZ?["�?YX?VE�?Tz�?S��?St�?S33?S33?R�?R�!?R�?R�?S�F?S�F?O\)?Ix�?Gl�?DZ?F�y?A��?>5???�w?A�7?BJ?Co?C��?DZ?E�T?H1'?G�?G+?F��?F��?G+?Hr�?G�?Hr�?I��?J=q?H�9?DZ?A��?>��?>�R??�w?@�?@Ĝ?A�7?AG�?A�7?BJ?B�\?CS�?D�?DZ?E�?E�T?Gl�?G+?E`B?E�?E�T?F�y?G�?I��?J��?J��?Ix�?Gl�?E�T?E��?E`B?E�?A��?AG�?@  ?>v�?;�m?;dZ?;dZ?;��?<(�?;�m?<j?9�?{?M�>�?}>�F>� �>���>� �>���>���?1?��?��?	x�?%?G�?��?r�?`B>�^5>�~�>�V>�1>��>��j>�h>�`B>�n�>�p�>���>�\)>�;d>�l�>�>��>��
>�5?>�(�>��>�M�>�Z>�(�>�t�>���>���>���>��>��>��>��>���>�o>l�D>^5?>Z�>O�;>D��>6E�>/�>0 �>/�>'�>�>	7L=���=���=�\)=L��;�o�ě����
        ��C��'+�+�o�e`B�ě�����w�,1��w��7L�u��+���T���`������G���������ٽ��#���m�   �	7L�O߾hs����w�%�T�+�0 ž49X�7KǾ>vɾ?|�I�^�L�;Q녾Z��_;d�dZ�ixվp�׾u�x���z�H�{�m�~�۾�%���\���������𾇮�����1'���9������C���C���ƨ�����`��n���t����Ͼ����
=�����������㾝/��;d��A���G����徤Z���/���/���/���/��`B��ff��r���1���h���h��{�������������׾��!���F���F��9X��E����پ��پ��#��dZ��j��|��  ���7����Ǯ�ȴ9��=q��C���I����;�����`���`��녾��ϾՁ�և+�և+�׍P�ؓu�ٙ��ڟ�111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ���������������?}���������������?}�����z�䛦�䛦��Z��9X��9X������Ͽ�F��t�����!��n���녿�hs��Ĝ��bN���ۿ�1�ա˿ѩ����������E�������w��^5���翸�9�����ɺ^��?}G�O�G�O��"��V��C����w�پ:^5��w�J�t��B�\��񪾔zᾗ
=�N��n��%�P�`=�E�>Xb>t�j>�
=>�^5?!%?-�h?0 �?P �?o��?��?���?�l�?�Q�?��?�$�?ȓu?�(�?�{?�?ڟ�?���?�/?Լj?��;?Ł?�ff?��?v�+?�bN?��H?�Z?���?���?�hs?�n�?�9X?���?�?}?�?}?���?��/?���?�o?���?��;?�;d?���?�p�?�V?��?�(�?�1?�1?��?�bN?�33?�?�
=?�+?�+?�K�?���?�Q�?���?��+?�`B?�bN?�~�?���?�{?�~�?�1'?�dZ?�I�?�v�?��;?��R?��?��`?�G�?�&�?�bN?���?��R?��?�V?��h?�/?��?��?���?�b?�K�?��!?}�-?|j?x�u?~5??�S�?���?~v�?s�F?t9X?r�!?r�!?qhs?p�`?r�!?s�F?st�?s��?v�+?w�P?{�m?}�-?��?��`?�%?��?�  ?�w?�M�?�33?���?� �?v?st�?q�?rn�?q��?n{?e�T?bJ?a%?\(�?Z�?Z^5?Z�?Z�H?["�?[dZ?["�?YX?VE�?Tz�?S��?St�?S33?S33?R�?R�!?R�?R�?S�F?S�F?O\)?Ix�?Gl�?DZ?F�y?A��?>5???�w?A�7?BJ?Co?C��?DZ?E�T?H1'?G�?G+?F��?F��?G+?Hr�?G�?Hr�?I��?J=q?H�9?DZ?A��?>��?>�R??�w?@�?@Ĝ?A�7?AG�?A�7?BJ?B�\?CS�?D�?DZ?E�?E�T?Gl�?G+?E`B?E�?E�T?F�y?G�?I��?J��?J��?Ix�?Gl�?E�T?E��?E`B?E�?A��?AG�?@  ?>v�?;�m?;dZ?;dZ?;��?<(�?;�m?<j?9�?{?M�>�?}>�F>� �>���>� �>���>���?1?��?��?	x�?%?G�?��?r�?`B>�^5>�~�>�V>�1>��>��j>�h>�`B>�n�>�p�>���>�\)>�;d>�l�>�>��>��
>�5?>�(�>��>�M�>�Z>�(�>�t�>���>���>���>��>��>��>��>���>�o>l�D>^5?>Z�>O�;>D��>6E�>/�>0 �>/�>'�>�>	7L=���=���=�\)=L��;�o�ě����
        ��C��'+�+�o�e`B�ě�����w�,1��w��7L�u��+���T���`������G���������ٽ��#���m�   �	7L�O߾hs����w�%�T�+�0 ž49X�7KǾ>vɾ?|�I�^�L�;Q녾Z��_;d�dZ�ixվp�׾u�x���z�H�{�m�~�۾�%���\���������𾇮�����1'���9������C���C���ƨ�����`��n���t����Ͼ����
=�����������㾝/��;d��A���G����徤Z���/���/���/���/��`B��ff��r���1���h���h��{�������������׾��!���F���F��9X��E����پ��پ��#��dZ��j��|��  ���7����Ǯ�ȴ9��=q��C���I����;�����`���`��녾��ϾՁ�և+�և+�׍P�ؓu�ٙ��ڟ�111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	-B	-B	.B	.B	.B	.B	-B	.B	.B	.B	-B	.B	.B	/B	/B	/B	0!B	/B	/B	/B	0!B	0!B	2-B	33B	5?B	9XB	=qB	@�B	I�B	XB	m�B	u�B	�B	�VB	��B	�B	�FB	ŢB	ŢB	�mB	�B	��B
  B
7LB	��B
5?B
33B
E�B
E�B
G�B
P�B
]/B
YB
[#B
t�B
_;B
e`B
o�B
hsB
hsB
u�B
t�B
�B
�7B
��B
�B
��B
�3B
��B
�B
�
B
�B
�B\B�B,B/B1'B5?B5?B:^B1'BC�BB�B>wB<jB:^B;dB�B �B�B&�B=qBF�BK�BM�BM�BO�BR�BT�BXBXBYBYBZB\)B^5B_;BbNB`BB`BBaHBbNBcTBcTBdZBgmBk�Bq�Bt�Bw�Bx�Bx�B{�B|�B}�B~�B~�B|�B{�By�Bx�Bu�Bs�Bv�Bs�Bs�B}�B~�B|�B}�B� B�B~�B� B�B�B�B�B�B�B�B�B�B�B�B�B~�B�B�B}�B�1B�JB�JB�1B�1B�1B�7B�JB�JB�JB�\B�bB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�!B�'B�-B�3B�9B�9B�3B�9B�9B�9B�3B�3B�9B�9B�3B�3B�3B�3B�3B�9B�3B�3B�-B�XB�'B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�'B�B�B�B�B�B�B�!B�!B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B	1�B	1�B	2�B	2�B	2�B	2�B	1�B	2�B	2�B	2�B	1�B	2�B	2�B	3�B	3�B	3�B	4�B	3�B	3�B	3�B	4�B	4�B	6�B	7�B	:	B	>!B	B7B	EOB	N�B	\�B	r]B	z�B	��B	�#B	��B	��B	�B	�qB	�qB	�?B	�cB	��B
�G�O�G�O�B
:B
8B
JyB
JxB
L�B
U�B
bB
]�B
_�B
y�B
dB
j6B
ttB
mIB
mKB
z�B
y�B
��B
�B
��B
��B
��B
�B
ӳB
��B
��B
��B
��B7BaB0�B3�B6B: B:"B??B6BHwBGpBCYBALB?@B@DBsB%�B �B+�BBQBK�BP�BR�BR�BT�BW�BY�B\�B\�B]�B]�B_BaBcBdBg0Be!Be%Bf*Bg2Bh8Bh8Bi;BlPBpeBv�By�B|�B}�B}�B��B��B��B��B��B��B��B~�B}�Bz�Bx�B{�Bx�Bx�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�,B�/B�B�B�B�B�-B�-B�.B�?B�GB�XB�`B�`B�~B�zB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�<B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B� B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�oB��B��B��B�}B��B��B��B��B��B�}B��B��B��B�uB�}B�}B�B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0047612 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040620220127170406  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144426  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144426  QCF$                G�O�G�O�G�O�0000000000004000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A333D�VfG�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170406  IP  PSAL            A333D�VfG�O�                