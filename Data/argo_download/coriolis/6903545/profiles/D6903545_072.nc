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
resolution        =���   axis      Z          :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  J�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       L�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       S�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Z�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       \�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  c�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       e\   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       lp   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       uL   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  |`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ~(   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �<   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �l   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �l   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �l   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �lArgo profile    3.1 1.2 19500101000000  20200829092230  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               HA   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @�w�[�1   @�x|e�@@P�̛����9O�G�1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 3.8 dbar]                                                      @�33@�33@�  @���A   A��A!��A0  A;33AL��A^ffAnffA~ffA���A�  A���A�33A�33A���A�33A�33A�  A�33A���A���A�33A�33A���A�33B33B  B  B  B��B��B33B33B#33B'��B,  B0ffB4  B7��B<  BA33BD  BH  BL  BP  BTffBXffB\  B^��Bb��Bh  Bk��Bq33Bs��Bw33Bz��B~��B���B���B�  B�33B�  B���B�  B���B���B�33B�33B�  B�33B���B���B�  B�  B�ffB���B���B�33B���B���B�  B�  B�33B�  B�33B���B���B�33B���B�  B�33B���B�ffB�ffB�33B�  B�  B�33Bޙ�B♚B晚BꙚB���B�  B�ffB���B�  C�3C� CL�C�C	ffC�3C��C��C��C��C��C� CffCL�CL�C33C!� C#��C%��C'��C)L�C+L�C-ffC/ffC1� C3� C5�3C7ffC9� C;��C=33C?L�CA� CC33CEffCG�3CI� CKL�CM�COffCQ��CSL�CU33CW33CYL�C[� C]ffC_L�Ca��Cc�3CeffCg� Ci��CkL�CmL�Co� CqL�Cs�CuffCw��Cy� C{ffC}L�CL�C���C���C���C���C���C���C���C��fC��fC���C���C���C��fC��fC��fC��fC���C���C���C���C���C���C���C���C���C���C���C���C��fC���C��fC��fC��fC���C���C���C���C���C���C���C���C�� C�� C�� C�� C�� C�� C���C���C���C���C�� C���C���C���C���C���C��fC��fC��fC��3C��3C�� C���C���C��3C�� CÙ�Cĳ3C���CƦfCǳ3C���CɦfCʳ3C���C̳3C͌�Cγ3C���CЦfC�� C�ٚCӳ3Cԙ�CզfCֳ3C׌�CئfCٳ3Cڌ�CۦfC�� CݦfC޳3C���C�3C��C�fC�� C�fC��C�3C���C�3C陚C�� C�ٚC�� C��3CC�� C�ٚC�� C�3C�C��C�� C��3C��fC�ٚC���C�ffC�ٚD &fD` D��D3DS3D��D��D	3D
Y�D� D�fD�3D,�DffD� D��DL�D��D�fD�D9�DffD��D��DFfD�fD � D"  D#9�D$� D%� D'�D(@ D)l�D*�fD+�fD-&fD.l�D/�fD0�3D29�D3l�D4�fD5�fD733D8s3D9�fD:��D<33D=�fD>�fD@  DA9�DBs3DC��DEfDF,�DGs3DH��DI��DK@ DL�3DM��DOfDPFfDQ�fDR��DT�DU9�DV` DW�3DY  DZ33D[ffD\��D^fD_33D`` Da��Db�3Dd@ De��Df��Dg� Di9�Dj�fDk��Dl�fDn9�Dos3Dp�fDr  DsL�Dt� Du�3Dw�Dx9�Dys3Dz�3D{�fD}&fD~y�D�3D���D�)�D�ɚD�i�D���D��3D�6fD�ٚD�|�D�#3D��fD�\�D��3D�� D�<�D���D�� D�  D��3D�VfD���D��3D�9�D��3D�p D��D�ɚD�c3D�3D���D�C3D�� D�|�D�3D���D�i�D���D��3D�9�D�� D��3D��D��3D�Y�D�  D���D�<�D�ٚD�vfD��D�s3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�33@�33@�  @���A   A��A!��A0  A;33AL��A^ffAnffA~ffA���A�  A���A�33A�33A���A�33A�33A�  A�33A���A���A�33A�33A���A�33B33B  B  B  B��B��B33B33B#33B'��B,  B0ffB4  B7��B<  BA33BD  BH  BL  BP  BTffBXffB\  B^��Bb��Bh  Bk��Bq33Bs��Bw33Bz��B~��B���B���B�  B�33B�  B���B�  B���B���B�33B�33B�  B�33B���B���B�  B�  B�ffB���B���B�33B���B���B�  B�  B�33B�  B�33B���B���B�33B���B�  B�33B���B�ffB�ffB�33B�  B�  B�33Bޙ�B♚B晚BꙚB���B�  B�ffB���B�  C�3C� CL�C�C	ffC�3C��C��C��C��C��C� CffCL�CL�C33C!� C#��C%��C'��C)L�C+L�C-ffC/ffC1� C3� C5�3C7ffC9� C;��C=33C?L�CA� CC33CEffCG�3CI� CKL�CM�COffCQ��CSL�CU33CW33CYL�C[� C]ffC_L�Ca��Cc�3CeffCg� Ci��CkL�CmL�Co� CqL�Cs�CuffCw��Cy� C{ffC}L�CL�C���C���C���C���C���C���C���C��fC��fC���C���C���C��fC��fC��fC��fC���C���C���C���C���C���C���C���C���C���C���C���C��fC���C��fC��fC��fC���C���C���C���C���C���C���C���C�� C�� C�� C�� C�� C�� C���C���C���C���C�� C���C���C���C���C���C��fC��fC��fC��3C��3C�� C���C���C��3C�� CÙ�Cĳ3C���CƦfCǳ3C���CɦfCʳ3C���C̳3C͌�Cγ3C���CЦfC�� C�ٚCӳ3Cԙ�CզfCֳ3C׌�CئfCٳ3Cڌ�CۦfC�� CݦfC޳3C���C�3C��C�fC�� C�fC��C�3C���C�3C陚C�� C�ٚC�� C��3CC�� C�ٚC�� C�3C�C��C�� C��3C��fC�ٚC���C�ffC�ٚD &fD` D��D3DS3D��D��D	3D
Y�D� D�fD�3D,�DffD� D��DL�D��D�fD�D9�DffD��D��DFfD�fD � D"  D#9�D$� D%� D'�D(@ D)l�D*�fD+�fD-&fD.l�D/�fD0�3D29�D3l�D4�fD5�fD733D8s3D9�fD:��D<33D=�fD>�fD@  DA9�DBs3DC��DEfDF,�DGs3DH��DI��DK@ DL�3DM��DOfDPFfDQ�fDR��DT�DU9�DV` DW�3DY  DZ33D[ffD\��D^fD_33D`` Da��Db�3Dd@ De��Df��Dg� Di9�Dj�fDk��Dl�fDn9�Dos3Dp�fDr  DsL�Dt� Du�3Dw�Dx9�Dys3Dz�3D{�fD}&fD~y�D�3D���D�)�D�ɚD�i�D���D��3D�6fD�ٚD�|�D�#3D��fD�\�D��3D�� D�<�D���D�� D�  D��3D�VfD���D��3D�9�D��3D�p D��D�ɚD�c3D�3D���D�C3D�� D�|�D�3D���D�i�D���D��3D�9�D�� D��3D��D��3D�Y�D�  D���D�<�D�ٚD�vfD��D�s3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���q녿lI��h�9�e`B�b��]�Xb�W�ٿX�u�Z��[��[�m�]/�^�R�_�w�c���c���d��e`B�g+�h1'�pbN��V��o���
���˿�������T��������?}���;��G�������Mӿ�񪿤�/���^5�����S����T���#��푿�;d��  �����%��Ĝ��A���;d���ۿ��R��;d���`��t��öF������z��`B�ǍP��
=��ȴ�Ƨ��E�����E��Ƈ+�Ł���T��z�öF��t���o������Ϳ����������������Ͽ��F��t���S���t���t���t�����������S���o��񪿲�!��Mӿ�Mӿ��!��Mӿ�S���  ��A��� ſ�vɿ����=q���;���忑&鿐Ĝ��\)��ƨ�������~�ۿ}p��|(��y��s�Ͽg_�w�Z^5�W�ٿP�`�@  �.V���������{������H���ؓu��hs������H��&龟;d�S��:�o>_;d>��>��>���>�J>���>��>�7L>�C�>�ƨ>t�j>��7>��\>���>�\)>��>���>�\)>�j?�7?Z?`B?�?�/? A�>�X>�?}>�ȴ?Z?��?ff?��>�|�>��H>�dZ>�`B>ܬ>�ƨ>�5??��?�#?"�?�?p�?v�?!��?"��?#��?#��?#S�?"��?"�\?5??��?��?��?"�?�?�#?��?�?��?�u?��?�+???��?n�?��?&�?�`?��?bN?bN? �?�?��?V?V?{?{?{?��?O�?�D?�9?`B?`B?��?S�?�\?M�?�7?   >���>�p�>��>��>���? A�? Ĝ? Ĝ? �>�|�>�|�>���>�v�>��>�v�?   >�v�>��>�v�>�p�>��m>�^5>��H>�^5>��m>��m>��H>��H>��H>��>�j>��H>�K�>�F>� �>�{>�D>�ff>�z�>�>�t�>�n�>�bN>��;>��;>��>�\)>��`>��>�t�>�
=>�;d>ݲ->�/>ܬ>ܬ>�(�>ۥ�>��>׍P>׍P>׍P>߾w>���>���>�J>��!>��h>��h>�~�>��/>�;d>��->���>��P>�O�>��^>��^>�7L>�+>�o>�J>���>��!>�C�>��>��>�b>�z�>�bN>ɺ^>�C�>�I�>Ձ>�M�>�;d>ؓu>�>��`>���>�C�>�7L>�o>�  >�v�>��H>�9X>�9X>�9X>��F>��h>���>��T>�M�>��>�>�hs>���>�+>w��>o��>fff>`A�>Y�>N�>D��>?|�>5?}>,1>(��>"��>�u>O�>+=�=�G�=ȴ9=�Q�=�-=�1=��=m�h=@�=49X<�`B<��
<49X;ě���o��`B�T���ě��\)�'H�9�}󶽕�������^5���ͽ�;d���m�O߾��$�/�(�þ'+�7KǾI�^�R�^5?�fff�m�h�z�H��%���9��I���\)��녾�������-��G����
��Z���y��1����� ž�&龶ȴ���۾�o�����+�ɺ^���`��
=�ݲ-������þ���{��33���پ�Q��Q��Q��Q��Q�������#��X��^5��dZ���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   �q녿lI��h�9�e`B�b��]�Xb�W�ٿX�u�Z��[��[�m�]/�^�R�_�w�c���c���d��e`B�g+�h1'�pbN��V��o���
���˿�������T��������?}���;��G�������Mӿ�񪿤�/���^5�����S����T���#��푿�;d��  �����%��Ĝ��A���;d���ۿ��R��;d���`��t��öF������z��`B�ǍP��
=��ȴ�Ƨ��E�����E��Ƈ+�Ł���T��z�öF��t���o������Ϳ����������������Ͽ��F��t���S���t���t���t�����������S���o��񪿲�!��Mӿ�Mӿ��!��Mӿ�S���  ��A��� ſ�vɿ����=q���;���忑&鿐Ĝ��\)��ƨ�������~�ۿ}p��|(��y��s�Ͽg_�w�Z^5�W�ٿP�`�@  �.V���������{������H���ؓu��hs������H��&龟;d�S��:�o>_;d>��>��>���>�J>���>��>�7L>�C�>�ƨ>t�j>��7>��\>���>�\)>��>���>�\)>�j?�7?Z?`B?�?�/? A�>�X>�?}>�ȴ?Z?��?ff?��>�|�>��H>�dZ>�`B>ܬ>�ƨ>�5??��?�#?"�?�?p�?v�?!��?"��?#��?#��?#S�?"��?"�\?5??��?��?��?"�?�?�#?��?�?��?�u?��?�+???��?n�?��?&�?�`?��?bN?bN? �?�?��?V?V?{?{?{?��?O�?�D?�9?`B?`B?��?S�?�\?M�?�7?   >���>�p�>��>��>���? A�? Ĝ? Ĝ? �>�|�>�|�>���>�v�>��>�v�?   >�v�>��>�v�>�p�>��m>�^5>��H>�^5>��m>��m>��H>��H>��H>��>�j>��H>�K�>�F>� �>�{>�D>�ff>�z�>�>�t�>�n�>�bN>��;>��;>��>�\)>��`>��>�t�>�
=>�;d>ݲ->�/>ܬ>ܬ>�(�>ۥ�>��>׍P>׍P>׍P>߾w>���>���>�J>��!>��h>��h>�~�>��/>�;d>��->���>��P>�O�>��^>��^>�7L>�+>�o>�J>���>��!>�C�>��>��>�b>�z�>�bN>ɺ^>�C�>�I�>Ձ>�M�>�;d>ؓu>�>��`>���>�C�>�7L>�o>�  >�v�>��H>�9X>�9X>�9X>��F>��h>���>��T>�M�>��>�>�hs>���>�+>w��>o��>fff>`A�>Y�>N�>D��>?|�>5?}>,1>(��>"��>�u>O�>+=�=�G�=ȴ9=�Q�=�-=�1=��=m�h=@�=49X<�`B<��
<49X;ě���o��`B�T���ě��\)�'H�9�}󶽕�������^5���ͽ�;d���m�O߾��$�/�(�þ'+�7KǾI�^�R�^5?�fff�m�h�z�H��%���9��I���\)��녾�������-��G����
��Z���y��1����� ž�&龶ȴ���۾�o�����+�ɺ^���`��
=�ݲ-������þ���{��33���پ�Q��Q��Q��Q��Q�������#��X��^5��dZ���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBv�B�uB��B�RBŢB�BM�B^5BdZBr�B�B�=B�{B��B��B�3B�?B�^B�}BƨBȴB�B�B)�B)�B.BA�BffBq�Bz�B�B�uB�-BBǮB��B�fBB�B)�Bs�B�1B��B��B�
B�B��B��B+BuB�B"�B%�B0!B49B>wBH�BJ�BXB�B��B�/B�ZB�`B�sB�B�B��B  B�B8RBC�BI�BW
B_;Br�B�B�+B�VB�uB��B�B�3B�LBBǮBȴBȴB��B��B��B��B�
B�#B�#B�5B�;B�TB�B	B	VB	�B	�B	 �B	'�B	6FB	B�B	J�B	T�B	YB	hsB	o�B	x�B	�B	�7B	�oB	��B	��B	�B	�LB	�qB	ĜB	��B	�B	�`B	�B	��B
B
%B
DB
%B
hB
�B
�B
"�B
+B
/B
:^B
D�B
VB
l�B
�B
�PB
�oB
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
�B
�?B
�qB
ƨB
ǮB
��B
�B
�B
�#B
�B
�B
�B
�#B
�#B
�/B
�BB
�NB
�fB
�sB
�mB
�mB
�ZB
�`B
�HB
�yB
�B
��B
��BBBBB1B
=BDBJBPBVB\BhBuB�B�B�B�B�B�B�B�B�B�B�B�B!�B#�B'�B(�B)�B)�B+B+B,B,B.B.B.B/B/B1'B2-B2-B2-B33B7LB7LB7LB8RB:^B;dB;dB;dB<jB=qB=qB=qB>wB?}BA�BA�BA�BA�BC�BC�BD�BE�BF�BF�BH�BJ�BK�BK�BK�BL�BL�BL�BM�BN�BN�BO�BP�BQ�BR�BR�BS�BT�BVBVBVBVBT�BVBW
BW
BYBYBYBYBXBYBYBXBZB[#B]/B^5B^5B^5B_;B_;B_;BaHBaHBaHBbNBe`Be`Be`BdZBcTBe`BffBffBgmBgmBffBffBgmBgmBiyBhsBiyBk�Bl�Bn�Bp�Bz�B�B�B�+B�+B�+B�+B�7B�=B�PB�bB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   Bv�B�uB��B�RBŢB�BM�B^5BdZBr�B�B�=B�{B��B��B�3B�?B�^B�}BƨBȴB�B�B)�B)�B.BA�BffBq�Bz�B�B�uB�-BBǮB��B�fBB�B)�Bs�B�1B��B��B�
B�B��B��B+BuB�B"�B%�B0!B49B>wBH�BJ�BXB�B��B�/B�ZB�`B�sB�B�B��B  B�B8RBC�BI�BW
B_;Br�B�B�+B�VB�uB��B�B�3B�LBBǮBȴBȴB��B��B��B��B�
B�#B�#B�5B�;B�TB�B	B	VB	�B	�B	 �B	'�B	6FB	B�B	J�B	T�B	YB	hsB	o�B	x�B	�B	�7B	�oB	��B	��B	�B	�LB	�qB	ĜB	��B	�B	�`B	�B	��B
B
%B
DB
%B
hB
�B
�B
"�B
+B
/B
:^B
D�B
VB
l�B
�B
�PB
�oB
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
�B
�?B
�qB
ƨB
ǮB
��B
�B
�B
�#B
�B
�B
�B
�#B
�#B
�/B
�BB
�NB
�fB
�sB
�mB
�mB
�ZB
�`B
�HB
�yB
�B
��B
��BBBBB1B
=BDBJBPBVB\BhBuB�B�B�B�B�B�B�B�B�B�B�B�B!�B#�B'�B(�B)�B)�B+B+B,B,B.B.B.B/B/B1'B2-B2-B2-B33B7LB7LB7LB8RB:^B;dB;dB;dB<jB=qB=qB=qB>wB?}BA�BA�BA�BA�BC�BC�BD�BE�BF�BF�BH�BJ�BK�BK�BK�BL�BL�BL�BM�BN�BN�BO�BP�BQ�BR�BR�BS�BT�BVBVBVBVBT�BVBW
BW
BYBYBYBYBXBYBYBXBZB[#B]/B^5B^5B^5B_;B_;B_;BaHBaHBaHBbNBe`Be`Be`BdZBcTBe`BffBffBgmBgmBffBffBgmBgmBiyBhsBiyBk�Bl�Bn�Bp�Bz�B�B�B�+B�+B�+B�+B�7B�=B�PB�bB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092230                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092337  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200829092337  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134655  IP  PSAL            @�33D�s3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172541  IP  PSAL            @�33D�s3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            @�33D�s3G�O�                