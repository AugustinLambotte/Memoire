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
resolution        =���   axis      Z        �  C4   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  I�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  K�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  R@   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Z�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  aL   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  pX   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  r   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  x�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  z`   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �l   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �p   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �t   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �x   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �|   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �@   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �@   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �@   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �@Argo profile    3.1 1.2 19500101000000  20200829092230  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               EA   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��ffff1   @��ffff@Q�S�
T�7��e}U\8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A8  ANffAfffA|��A���A���A���A���A���Ař�Aљ�A���A���A���B ��B��B��B��B��B33B%33B+��B2  B8ffB>��BE33BK��BR  BXffB^��Be33Bk��Br  Bx��B33B�  B�33B���B���B�33B���B���B�33B���B�  B�ffB���B�33B���B�  B�ffB���B�33B���B�  B�ffB�  B�  B�33Bי�B�  B�ffB�  B�  B�  B�ffB�  B�ffB���B���C�C��C� C33C	  C� C33C�fC��CffC�C�fC� C33C�fC�3C!ffC#33C$�fC'� C*33C,  C-�3C/� C133C3  C5� C833C9�fC;��C=L�C?�C@��CCffCF  CG��CI� CK33CM  CO� CR33CT  CU�3CWffCY33CZ�fC]� C`�Ca�fCc��CeL�Cg�Ch��CkL�Cm�fCo��CqL�Cs  Cu� Cx�Cy��C{ffC}�C��C��C��fC��3C���C�ffC��fC��fC��3C�� C�� C�  C���C��fC�s3C��fC��fC��3C�� C�� C��3C�� C���C�� C��3C�� C���C���C�  C���C��fC�s3C��fC��fC��3C���C�ffC��fC��fC��3C���C�ffC��fC��fC���C��fC�� C�� C��C��3C�ٚC��3C���C�s3C�� C��C��3C�ٚC��3C���C�� C�ffC��3C��C��3C���C��3C�C�s3C�� C��C��3C�ٚC�� Cə�Cʀ C�ffC̳3C��C��3C���Cг3Cљ�C�s3C�Y�CԳ3C��C��3C�ٚC�� C٦fCڌ�Cۀ C�ffC�� C��C�  C��fC���C�� C�fC��C� C�ffC�L�C�3C�&fC��C��C��3C��fC�ٚC���C�3C�fC�C��C�s3C�Y�C�� C�&fC��C��3C�Y�C��fD 33D�3D��DٚD9�D�3D� D� D
S3D� D�3D� D�D��D� D�3D,�D` D�3DfD9�Ds3D�fD� D�D��D �fD!��D#,�D$ffD%� D&ٚD(L�D)�fD*� D+��D-33D.s3D/�3D0�3D29�D3y�D4��D5��D79�D8y�D9��D:��D<  D=��D>�3D@�DAFfDB� DC��DD�3DF,�DGffDH� DI� DK  DL` DM� DN�fDP&fDQs3DR�3DT  DUL�DVS3DW�fDX�3DZFfD[��D\��D^  D_S3D`ffDa��DcfDd3De` Df��Dg��Di33Dj��Dk��Dl��Dn@ DoS3Dp��DrfDs  Dt� Du� Dv�3Dx�Dys3DzٚD{�3D}�D~l�D�3D�vfD�fD��fD�i�D���D��fD�6fD��fD�vfD�#3D�ٚD�l�D���D�� D�#3D�ٚD�� D�#3D��fD�L�D�� D���D�L�D�� D�s3D�fD�� D�vfD�	�D���D�0 D�� D�vfD�,�D�� D�P D�� D��fD�L�D���D�l�D�  D���D�p D�3D��fD�,�D�� D�|�D�9�D�� D�ffD���D��fD�,�D��fD�\�D��D��fD�p D�fD�� D�6fD�ɚD�|�D�6fD�ٚ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A8  ANffAfffA|��A���A���A���A���A���Ař�Aљ�A���A���A���B ��B��B��B��B��B33B%33B+��B2  B8ffB>��BE33BK��BR  BXffB^��Be33Bk��Br  Bx��B33B�  B�33B���B���B�33B���B���B�33B���B�  B�ffB���B�33B���B�  B�ffB���B�33B���B�  B�ffB�  B�  B�33Bי�B�  B�ffB�  B�  B�  B�ffB�  B�ffB���B���C�C��C� C33C	  C� C33C�fC��CffC�C�fC� C33C�fC�3C!ffC#33C$�fC'� C*33C,  C-�3C/� C133C3  C5� C833C9�fC;��C=L�C?�C@��CCffCF  CG��CI� CK33CM  CO� CR33CT  CU�3CWffCY33CZ�fC]� C`�Ca�fCc��CeL�Cg�Ch��CkL�Cm�fCo��CqL�Cs  Cu� Cx�Cy��C{ffC}�C��C��C��fC��3C���C�ffC��fC��fC��3C�� C�� C�  C���C��fC�s3C��fC��fC��3C�� C�� C��3C�� C���C�� C��3C�� C���C���C�  C���C��fC�s3C��fC��fC��3C���C�ffC��fC��fC��3C���C�ffC��fC��fC���C��fC�� C�� C��C��3C�ٚC��3C���C�s3C�� C��C��3C�ٚC��3C���C�� C�ffC��3C��C��3C���C��3C�C�s3C�� C��C��3C�ٚC�� Cə�Cʀ C�ffC̳3C��C��3C���Cг3Cљ�C�s3C�Y�CԳ3C��C��3C�ٚC�� C٦fCڌ�Cۀ C�ffC�� C��C�  C��fC���C�� C�fC��C� C�ffC�L�C�3C�&fC��C��C��3C��fC�ٚC���C�3C�fC�C��C�s3C�Y�C�� C�&fC��C��3C�Y�C��fD 33D�3D��DٚD9�D�3D� D� D
S3D� D�3D� D�D��D� D�3D,�D` D�3DfD9�Ds3D�fD� D�D��D �fD!��D#,�D$ffD%� D&ٚD(L�D)�fD*� D+��D-33D.s3D/�3D0�3D29�D3y�D4��D5��D79�D8y�D9��D:��D<  D=��D>�3D@�DAFfDB� DC��DD�3DF,�DGffDH� DI� DK  DL` DM� DN�fDP&fDQs3DR�3DT  DUL�DVS3DW�fDX�3DZFfD[��D\��D^  D_S3D`ffDa��DcfDd3De` Df��Dg��Di33Dj��Dk��Dl��Dn@ DoS3Dp��DrfDs  Dt� Du� Dv�3Dx�Dys3DzٚD{�3D}�D~l�D�3D�vfD�fD��fD�i�D���D��fD�6fD��fD�vfD�#3D�ٚD�l�D���D�� D�#3D�ٚD�� D�#3D��fD�L�D�� D���D�L�D�� D�s3D�fD�� D�vfD�	�D���D�0 D�� D�vfD�,�D�� D�P D�� D��fD�L�D���D�l�D�  D���D�p D�3D��fD�,�D�� D�|�D�9�D�� D�ffD���D��fD�,�D��fD�\�D��D��fD�p D�fD�� D�6fD�ɚD�|�D�6fD�ٚ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���a%�Z�H�Vȴ�Q녿I7L�;��7
=�1���/��(�9�%`B�p�����G��ڟ���V���!������>/�>���>�K�?�D?E�?��?�?,1?A%?T��?e�?wK�?{��?���?�33?���?��+?���?��`?��!?�z�?��?�7L?��H?��?��?<�?&��?hs?#o?>5??J=q?BJ?D��?=p�?W
=?[dZ?e�?yX?f�y?g+?f$�?d�?co?W��?&�y?5??
=?p�?"M�?.V?/\)?,I�?,�D?/��?5?}?AG�?W�P?j~�?s�F?y��?y�?yX?z�H?{"�?}p�?~v�?���?�Ĝ?�%?��7?��\?�o?��F?��j?��j?�`B?��+?�+?��?���?�Q�?���?�^5?��#?���?�?��/?��F?�o?��F?�o?��7?� �?|�?{dZ?x�u?t9X?s�F?pbN?bJ?X��?X��?YX?X�u?VE�?S��?O�?KC�?9��?;d?\)?
~�?��>��>�?}>���>��>���>��>T��>C��>H�9>0 �>J>Y�>�  >�\)>�S�>���>ؓu>���>�&�>�>�V>�O�>�\)>�bN>��>�t�>��>�`B>�/>��?	��?
=q?�7?�
?	��?�h?1???	��?	7L?�?�h?��?��?�?�`?t�?t�?��?��?9X?z�??}??}?E�?E�?
=?�P?X?��?��? Ĝ?$Z?$��?%`B?'+?,I�?�w?�H??K�?#S�?)�^?/�??�w?:��?'�?�D>��`>���>�%>\>��>�7L>�r�>���>�M�>��>��->��w>�Ĝ>�Ĝ>�Ĝ>���>��/>��y>�r�>���>���>���>���>�~�>��>�V>�->�j>�v�>�  >�  >��>�dZ>��#>�9X>��j>�?}>�K�>�^5>�v�>�  >�$�>�=q>�ƨ>���>���>�\)>�$�>���>���>���>�J>�J>�o>���>�$�>��9>�C�>��`>�n�>���>���>��>��+>��u>�"�>�/>��R>�A�>�5?>���>���>�->�E�>��#>�v�>�%>�%>��7>\>�  >�%>�J>\>�=q>�hs>��`>�V>ɺ^>Ǯ>�C�>�ƨ>�C�>���>Ƨ�>�J>��m>��H>�dZ>�dZ>��H>�>�>�S�>�;d>��>�z�>��>�ƨ>���>s�F>cS�>P�`>=p�>.{>!��>�>�>z�>n�>bN>O�=���=�;d=��=���=��T=���=��=��=��-=��w=��w=���=��
=��T=��-=]/=\)<�C�;ě�:�o�ě���t����+�#�
�<j�H�9�y�#��\)���㽛�㽛�㽝�-�� Ž�����`��S���h�����#�   �$ݾ
=q�hs��u�$�/�)��0 ž333�7KǾ;dZ�B�\�B�\�O�;�V�W
=�^5?�`A��fff�o���y�#�}󶾀  ��  ��o������V���`��녾�
=��"Ѿ�Ĝ���T���羫���{������9X��X��dZ�Õ��Ƨ��7L��C���V��t������S���D��&���j���H��p���^511111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111�a%�Z�H�Vȴ�Q녿I7L�;��7
=�1���/��(�9�%`B�p�����G��ڟ���V���!������>/�>���>�K�?�D?E�?��?�?,1?A%?T��?e�?wK�?{��?���?�33?���?��+?���?��`?��!?�z�?��?�7L?��H?��?��?<�?&��?hs?#o?>5??J=q?BJ?D��?=p�?W
=?[dZ?e�?yX?f�y?g+?f$�?d�?co?W��?&�y?5??
=?p�?"M�?.V?/\)?,I�?,�D?/��?5?}?AG�?W�P?j~�?s�F?y��?y�?yX?z�H?{"�?}p�?~v�?���?�Ĝ?�%?��7?��\?�o?��F?��j?��j?�`B?��+?�+?��?���?�Q�?���?�^5?��#?���?�?��/?��F?�o?��F?�o?��7?� �?|�?{dZ?x�u?t9X?s�F?pbN?bJ?X��?X��?YX?X�u?VE�?S��?O�?KC�?9��?;d?\)?
~�?��>��>�?}>���>��>���>��>T��>C��>H�9>0 �>J>Y�>�  >�\)>�S�>���>ؓu>���>�&�>�>�V>�O�>�\)>�bN>��>�t�>��>�`B>�/>��?	��?
=q?�7?�
?	��?�h?1???	��?	7L?�?�h?��?��?�?�`?t�?t�?��?��?9X?z�??}??}?E�?E�?
=?�P?X?��?��? Ĝ?$Z?$��?%`B?'+?,I�?�w?�H??K�?#S�?)�^?/�??�w?:��?'�?�D>��`>���>�%>\>��>�7L>�r�>���>�M�>��>��->��w>�Ĝ>�Ĝ>�Ĝ>���>��/>��y>�r�>���>���>���>���>�~�>��>�V>�->�j>�v�>�  >�  >��>�dZ>��#>�9X>��j>�?}>�K�>�^5>�v�>�  >�$�>�=q>�ƨ>���>���>�\)>�$�>���>���>���>�J>�J>�o>���>�$�>��9>�C�>��`>�n�>���>���>��>��+>��u>�"�>�/>��R>�A�>�5?>���>���>�->�E�>��#>�v�>�%>�%>��7>\>�  >�%>�J>\>�=q>�hs>��`>�V>ɺ^>Ǯ>�C�>�ƨ>�C�>���>Ƨ�>�J>��m>��H>�dZ>�dZ>��H>�>�>�S�>�;d>��>�z�>��>�ƨ>���>s�F>cS�>P�`>=p�>.{>!��>�>�>z�>n�>bN>O�=���=�;d=��=���=��T=���=��=��=��-=��w=��w=���=��
=��T=��-=]/=\)<�C�;ě�:�o�ě���t����+�#�
�<j�H�9�y�#��\)���㽛�㽛�㽝�-�� Ž�����`��S���h�����#�   �$ݾ
=q�hs��u�$�/�)��0 ž333�7KǾ;dZ�B�\�B�\�O�;�V�W
=�^5?�`A��fff�o���y�#�}󶾀  ��  ��o������V���`��녾�
=��"Ѿ�Ĝ���T���羫���{������9X��X��dZ�Õ��Ƨ��7L��C���V��t������S���D��&���j���H��p���^511111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�B	!�B	/B	C�B	_;B	o�B	z�B	�+B	�=B	�\B	��B	��B	��B	�B	�B	�qB	ƨB	��B	��B
�B
$�B
H�B
VB
YB
_;B
bNB
_;B
r�B
s�B
�B
�JB
�VB
�VB
�{B
��B
��B
��B
��B
��B
��B
�uB
��B
��B
��B
��B
hsB
�B
�%B
u�B
�B
�B
�VB
�PB
�bB
�{B
�{B
��B
��B
��B
��B
�{B
�uB
��B
�DB
�=B
�B
�DB
�B
�DB
�\B
�bB
�bB
�hB
�{B
��B
��B
��B
�{B
�B
�B
�B
�B
�!B
�!B
�3B
�'B
�3B
�-B
�-B
�3B
�9B
�?B
�?B
�?B
�FB
�9B
�FB
�RB
�XB
�XB
�dB
�dB
�qB
�^B
�RB
�^B
�XB
�LB
�qB
�^B
�^B
�jB
�dB
�qB
�jB
�dB
�^B
�dB
�RB
�XB
�^B
�^B
�dB
�dB
�^B
�dB
�^B
�?B
�XB
�-B
�9B
�3B
�9B
�XB
�^B
�FB
�9B
�FB
�XB
�LB
�3B
�3B
�qB
�dB
ŢB
�-B
��B
��B
�)B
�;B
�B
�B
�B
�/B
�;B
�BB
�NB
�TB
�fB
�mB
�B
��B
��BB+B1B\B{B�B{B�BuB�B�B�B�B�B�B�B�B!�B!�B"�B#�B&�B$�B(�B+B,B-B.B.B-B/B0!B33B6FB7LB8RB6FB:^B;dB;dB;dB>wB5?BF�BF�BH�BO�BJ�BB�B6FB8RB7LB9XB7LB8RB8RB5?B6FB5?B6FB8RB9XB9XB9XB;dB;dB;dB<jB<jB=qB=qB>wB=qB>wB>wB@�BG�BI�BJ�BJ�BK�BK�BL�BM�BN�BN�BP�BP�BT�BT�BXBXBYBZBVBS�BT�BT�BT�BS�BVBW
BYBZBZB\)B[#B^5BaHBcTBbNBbNBbNBcTBe`Be`BgmBiyBk�Bo�Bq�Bs�Bu�Bw�Bz�B|�B|�B}�B~�B�B�B�B�%B�DB�VB�\B�hB�oB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B��B�B��B��B��B��B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B	�B	!�B	/B	C�B	_;B	o�B	z�B	�+B	�=B	�\B	��B	��B	��B	�B	�B	�qB	ƨB	��B	��B
�B
$�B
H�B
VB
YB
_;B
bNB
_;B
r�B
s�B
�B
�JB
�VB
�VB
�{B
��B
��B
��B
��B
��B
��B
�uB
��B
��B
��B
��B
hsB
�B
�%B
u�B
�B
�B
�VB
�PB
�bB
�{B
�{B
��B
��B
��B
��B
�{B
�uB
��B
�DB
�=B
�B
�DB
�B
�DB
�\B
�bB
�bB
�hB
�{B
��B
��B
��B
�{B
�B
�B
�B
�B
�!B
�!B
�3B
�'B
�3B
�-B
�-B
�3B
�9B
�?B
�?B
�?B
�FB
�9B
�FB
�RB
�XB
�XB
�dB
�dB
�qB
�^B
�RB
�^B
�XB
�LB
�qB
�^B
�^B
�jB
�dB
�qB
�jB
�dB
�^B
�dB
�RB
�XB
�^B
�^B
�dB
�dB
�^B
�dB
�^B
�?B
�XB
�-B
�9B
�3B
�9B
�XB
�^B
�FB
�9B
�FB
�XB
�LB
�3B
�3B
�qB
�dB
ŢB
�-B
��B
��B
�)B
�;B
�B
�B
�B
�/B
�;B
�BB
�NB
�TB
�fB
�mB
�B
��B
��BB+B1B\B{B�B{B�BuB�B�B�B�B�B�B�B�B!�B!�B"�B#�B&�B$�B(�B+B,B-B.B.B-B/B0!B33B6FB7LB8RB6FB:^B;dB;dB;dB>wB5?BF�BF�BH�BO�BJ�BB�B6FB8RB7LB9XB7LB8RB8RB5?B6FB5?B6FB8RB9XB9XB9XB;dB;dB;dB<jB<jB=qB=qB>wB=qB>wB>wB@�BG�BI�BJ�BJ�BK�BK�BL�BM�BN�BN�BP�BP�BT�BT�BXBXBYBZBVBS�BT�BT�BT�BS�BVBW
BYBZBZB\)B[#B^5BaHBcTBbNBbNBbNBcTBe`Be`BgmBiyBk�Bo�Bq�Bs�Bu�Bw�Bz�B|�B|�B}�B~�B�B�B�B�%B�DB�VB�\B�hB�oB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B��B�B��B��B��B��B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092230                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092334  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092334  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134655  IP  PSAL            A8  D�ٚG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172541  IP  PSAL            A8  D�ٚG�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A8  D�ٚG�O�                