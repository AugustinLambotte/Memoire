CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2016-12-19T13:55:46Z creation; 2023-09-07T16:10:45Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8    DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8    PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8(   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8h   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     9   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9,   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9L   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9l   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9p   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9x   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >��	4E�   
_FillValue        A.�~            9|   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Bh   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  DX   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  L   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N    TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]p   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _`   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  p�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  xx   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  zh   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �$   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �(   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �,   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �0   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �4   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �t   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20161219135546  20230907161045  3901875 MOCCA-EU                                                        Romain Cancouet                                                 PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR                           AR2600-16FR038                  5900A00                         844 @��a��O�1   @��\�O� @R?-�      1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 250 dbar; 10 sec sampling, 2 dbar average from 250 dbar to 1400 dbar; 10 sec sampling, 4 dbar average from 1400 dbar to 1000 dbar]                                                     A+33AH  Ac33A�ffA�  A���A���A�  A�ffAᙚA�B33B
  B��B��B!33B)33B1��B:ffBA��BI33BR  BZffBa��Bh��Bq33By��B�  B�  B���B���B���B�  B�  B���B���B���B���B���B���B�33B���B���B���Bę�Bș�B�  B���B���Bؙ�Bܙ�B�  B�  B���B���B���B�  B�  B�33C L�CL�CffCffC� C
� C� CffC� C��C� CL�C� C�3C� CL�C 33C"ffC$� C&L�C(ffC*� C,L�C.L�C0� C2L�C4ffC6��C8ffC:� C<��C>ffC@33CBffCD� CFL�CH33CJL�CL� CNL�CP� CR�3CT� CVffCX33CZffC\��C^ffC`L�Cb33CdffCf��Ch� CjffClL�Cn� Cp�3Cr� CtffCv33CxffCz� C|ffC~L�C��C�33C�L�C�@ C�33C��C�33C�@ C�&fC�33C�@ C�@ C�L�C�&fC�33C�33C�&fC�@ C�33C��C�33C�33C��C�&fC��C�&fC��C�33C�&fC�33C�L�C�33C�@ C��C�@ C�33C�33C�&fC�L�C�L�C�L�C�&fC�&fC�&fC��C�@ C�&fC��C�33C��C�@ C�&fC��C�&fC�33C�&fC�&fC��C�33C�@ C�@ C�@ C�@ C�33C�&fC�33C�&fC�33C�&fC�&fC�33C��C�&fC�&fC�33C�&fC�&fC�&fC�&fC��C�@ C�&fC�@ C��C�&fC��C�&fC�@ C�@ C�@ C�@ C�33C�&fC��C��C�33C��C�&fC�&fC�&fC��C��C�&fC�33C�@ C�33C�33C��C�33C�@ C�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC��C��C�&fC�33C�@ C�L�C�33C�&fC�&fD 3D �3DfD��D  D�fD  D� D�D��D�D��DfD��D3D�3D�D�3D	�D	� D
3D
��D3D��D�D� D  D�3D3D�3D3D�3D3D��D  D� D�D�3D  D�3D�D�3D�D� D3D��D�D��D3D� D3D��D�D��D  D��D3D�3D3D��D�D�3D�D��D   D � D!  D!�3D"�D"��D#�D#�3D$�D$�3D%�D%��D&�D&�3D'3D'��D(�D(�3D)�D)��D*�D*��D+fD+��D,3D,� D-�D-��D.3D.��D/3D/��D03D0�3D1�D1�fD2  D2�fD3&fD3��D43D4��D5&fD5�fD6  D6��D73D7� D83D8� D93D9��D:�D:�3D;3D;��D<�D<��D=�D=��D>�D>��D?  D?� D@3D@��DA3DA� DB3DB�3DC3DC� DD�DD��DE�DE��DF3DF�3DG�DG��DH�DH�3DI�DI��DJ�DJ��DK3DK�3DL�DL��DM3DM��DN3DN�3DO3DO�3DP�DP��DQ�DQ�3DR  DR�3DS3DS��DT�DT��DU3DU��DV�DV�3DW3DW�3DXfDX��DY3DY�3DZ�DZ��D[3D[��D\3D\��D]�D]��D^3D^�3D_�D_��D`�D`�3Da3Da� Db3Db�3Dc�Dc�3Dd3Dd��De�De��Df3Df�3DgfDg��Dh3Dh� Di  Di��Dj3Dj�3Dk�Dk��Dl�Dl��Dm3Dm��Dn3Dn�3Do3Do�3Dp�Dp�3Dq3Dq��Dr  Dr�3Ds�Ds�3Dt3Dt��Du3Du��Dv�Dv��Dw�Dw��Dx3Dx��Dx� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A+33AH  Ac33A�ffA�  A���A���A�  A�ffAᙚA�B33B
  B��B��B!33B)33B1��B:ffBA��BI33BR  BZffBa��Bh��Bq33By��B�  B�  B���B���B���B�  B�  B���B���B���B���B���B���B�33B���B���B���Bę�Bș�B�  B���B���Bؙ�Bܙ�B�  B�  B���B���B���B�  B�  B�33C L�CL�CffCffC� C
� C� CffC� C��C� CL�C� C�3C� CL�C 33C"ffC$� C&L�C(ffC*� C,L�C.L�C0� C2L�C4ffC6��C8ffC:� C<��C>ffC@33CBffCD� CFL�CH33CJL�CL� CNL�CP� CR�3CT� CVffCX33CZffC\��C^ffC`L�Cb33CdffCf��Ch� CjffClL�Cn� Cp�3Cr� CtffCv33CxffCz� C|ffC~L�C��C�33C�L�C�@ C�33C��C�33C�@ C�&fC�33C�@ C�@ C�L�C�&fC�33C�33C�&fC�@ C�33C��C�33C�33C��C�&fC��C�&fC��C�33C�&fC�33C�L�C�33C�@ C��C�@ C�33C�33C�&fC�L�C�L�C�L�C�&fC�&fC�&fC��C�@ C�&fC��C�33C��C�@ C�&fC��C�&fC�33C�&fC�&fC��C�33C�@ C�@ C�@ C�@ C�33C�&fC�33C�&fC�33C�&fC�&fC�33C��C�&fC�&fC�33C�&fC�&fC�&fC�&fC��C�@ C�&fC�@ C��C�&fC��C�&fC�@ C�@ C�@ C�@ C�33C�&fC��C��C�33C��C�&fC�&fC�&fC��C��C�&fC�33C�@ C�33C�33C��C�33C�@ C�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC��C��C�&fC�33C�@ C�L�C�33C�&fC�&fD 3D �3DfD��D  D�fD  D� D�D��D�D��DfD��D3D�3D�D�3D	�D	� D
3D
��D3D��D�D� D  D�3D3D�3D3D�3D3D��D  D� D�D�3D  D�3D�D�3D�D� D3D��D�D��D3D� D3D��D�D��D  D��D3D�3D3D��D�D�3D�D��D   D � D!  D!�3D"�D"��D#�D#�3D$�D$�3D%�D%��D&�D&�3D'3D'��D(�D(�3D)�D)��D*�D*��D+fD+��D,3D,� D-�D-��D.3D.��D/3D/��D03D0�3D1�D1�fD2  D2�fD3&fD3��D43D4��D5&fD5�fD6  D6��D73D7� D83D8� D93D9��D:�D:�3D;3D;��D<�D<��D=�D=��D>�D>��D?  D?� D@3D@��DA3DA� DB3DB�3DC3DC� DD�DD��DE�DE��DF3DF�3DG�DG��DH�DH�3DI�DI��DJ�DJ��DK3DK�3DL�DL��DM3DM��DN3DN�3DO3DO�3DP�DP��DQ�DQ�3DR  DR�3DS3DS��DT�DT��DU3DU��DV�DV�3DW3DW�3DXfDX��DY3DY�3DZ�DZ��D[3D[��D\3D\��D]�D]��D^3D^�3D_�D_��D`�D`�3Da3Da� Db3Db�3Dc�Dc�3Dd3Dd��De�De��Df3Df�3DgfDg��Dh3Dh� Di  Di��Dj3Dj�3Dk�Dk��Dl�Dl��Dm3Dm��Dn3Dn�3Do3Do�3Dp�Dp�3Dq3Dq��Dr  Dr�3Ds�Ds�3Dt3Dt��Du3Du��Dv�Dv��Dw�Dw��Dx3Dx��Dx� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�Z@�Z@�Q�@� �@��@�b@Ͼw@ΰ!@ȓu@���@���@��R@1�#@@Q�?�(�?�?}?���?���?ɺ^?ě�?���?� �?��P?�5??��h?�I�?���?��?�1?��?���?~5??vE�?mO�?e��?]�-?T9X?O�;?Co?8��?2�!?1hs?+C�?&��?   ?j?dZ?��?��?
~�?��?�?�?�?+?ff?��?S�?   >��>��>�X>�?}>��>�h>�h>���>���>ݲ->�(�>��>ۥ�>�5?>ۥ�>�b>�n�>�O�>�ƨ>�7L>�ƨ>�O�>��>ܬ>�A�>߾w>ۥ�>���>��
>ۥ�>��>ɺ^>�E�>�Ĝ>���>�"�>�r�>��>���>��`>���>���>��9>�1'>���>���>�
=>�
=>�>�t�>��>��>�J>}�>u>z�H>�J>��>�o>��>w��>fff>]/>\(�>["�>\(�>Xb>Q�>M��>@�>5?}>+>�R>z�>\)>I�>�=��=��m=��#=�l�=�
==��=��`=ě�=�Q�=�Q�=�E�=�-=���=��
=��-=��=�7L=�o=}�=m�h=]/=L��=49X=�w=C�<�h<�1<�t�<u<T��<D��<D��<49X<o;ě�;ě�;��
;D��;D��:�o    �o���
��`B�e`B��C����㼬1���
��1��9X��9X��j��9X�ě���9X���
���㼛�㼣�
���
���
��1��j�ě����ͼ�/��/��/��`B���o�\)��P��w�#�
�0 Ž0 Ž8Q�D���D���H�9�L�ͽT���]/�e`B�aG��m�h�y�#�}󶽅���O߽�hs������P���㽝�-���w���
���T���T��1�� Ž�9X��E���vɽ���\�ě��Ƨ������`����
=��;d��`B��xս�h��������#���m���m�   �%�o�o�J�J�����$ݾ$ݾ+�+�1'�1'�C��O߾\)�V�\)�bN�t���+�����R�!���#�
�$�/�$�/�%�T�(�þ)��+�,1�,1�.{�2-�49X�6E��6E��7KǾ8Q�9X�9X�:^5�:^5�;dZ�=p��=p��>vɾ@��B�\�C���D���D���E�˾F��G��I�^�L�;M��O�;�Q녾R�T���V�W
=�Xb�Xb�Y��Z��\(��]/�]/�^5?�_;d�`A��bMӾbMӾe`B�fff�gl��hr��hr��hr��ixվixվk��k��n���p�׾q���r�!�s�F�u�vȴ�x���z�H�|푾}�}�~�۾�  ��������%��%���\���\��������������˾�$ݾ���+�����1'��7L���^���^��=q������C���C���C���C���ƨ��ƨ��I����;�O߾�O߾�O߾���V������;���`��hs��n����Ͼ�������������
=���P���u��������������"Ѿ��㾜(������/���-���-���R���R��;d���w��A���Ĝ������MӾ��徣S����
��Z��`B���T���T���y���r����þ�xվ��義�義�羪~��������1��{��������������׾�&龱����-���!��33���F���j��?}����E���E���KǾ�KǾ��پ�Q쾸�����#���#���#��^5���H��dZ���m��j��푾�󶾾vɾ��۾�  ��������%���7�\��o�Õ��������š˾�$ݾ�+�Ǯ�Ǯ��7L11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�Z@�Z@�Q�@� �@��@�b@Ͼw@ΰ!@ȓu@���@���@��R@1�#@@Q�?�(�?�?}?���?���?ɺ^?ě�?���?� �?��P?�5??��h?�I�?���?��?�1?��?���?~5??vE�?mO�?e��?]�-?T9X?O�;?Co?8��?2�!?1hs?+C�?&��?   ?j?dZ?��?��?
~�?��?�?�?�?+?ff?��?S�?   >��>��>�X>�?}>��>�h>�h>���>���>ݲ->�(�>��>ۥ�>�5?>ۥ�>�b>�n�>�O�>�ƨ>�7L>�ƨ>�O�>��>ܬ>�A�>߾w>ۥ�>���>��
>ۥ�>��>ɺ^>�E�>�Ĝ>���>�"�>�r�>��>���>��`>���>���>��9>�1'>���>���>�
=>�
=>�>�t�>��>��>�J>}�>u>z�H>�J>��>�o>��>w��>fff>]/>\(�>["�>\(�>Xb>Q�>M��>@�>5?}>+>�R>z�>\)>I�>�=��=��m=��#=�l�=�
==��=��`=ě�=�Q�=�Q�=�E�=�-=���=��
=��-=��=�7L=�o=}�=m�h=]/=L��=49X=�w=C�<�h<�1<�t�<u<T��<D��<D��<49X<o;ě�;ě�;��
;D��;D��:�o    �o���
��`B�e`B��C����㼬1���
��1��9X��9X��j��9X�ě���9X���
���㼛�㼣�
���
���
��1��j�ě����ͼ�/��/��/��`B���o�\)��P��w�#�
�0 Ž0 Ž8Q�D���D���H�9�L�ͽT���]/�e`B�aG��m�h�y�#�}󶽅���O߽�hs������P���㽝�-���w���
���T���T��1�� Ž�9X��E���vɽ���\�ě��Ƨ������`����
=��;d��`B��xս�h��������#���m���m�   �%�o�o�J�J�����$ݾ$ݾ+�+�1'�1'�C��O߾\)�V�\)�bN�t���+�����R�!���#�
�$�/�$�/�%�T�(�þ)��+�,1�,1�.{�2-�49X�6E��6E��7KǾ8Q�9X�9X�:^5�:^5�;dZ�=p��=p��>vɾ@��B�\�C���D���D���E�˾F��G��I�^�L�;M��O�;�Q녾R�T���V�W
=�Xb�Xb�Y��Z��\(��]/�]/�^5?�_;d�`A��bMӾbMӾe`B�fff�gl��hr��hr��hr��ixվixվk��k��n���p�׾q���r�!�s�F�u�vȴ�x���z�H�|푾}�}�~�۾�  ��������%��%���\���\��������������˾�$ݾ���+�����1'��7L���^���^��=q������C���C���C���C���ƨ��ƨ��I����;�O߾�O߾�O߾���V������;���`��hs��n����Ͼ�������������
=���P���u��������������"Ѿ��㾜(������/���-���-���R���R��;d���w��A���Ĝ������MӾ��徣S����
��Z��`B���T���T���y���r����þ�xվ��義�義�羪~��������1��{��������������׾�&龱����-���!��33���F���j��?}����E���E���KǾ�KǾ��پ�Q쾸�����#���#���#��^5���H��dZ���m��j��푾�󶾾vɾ��۾�  ��������%���7�\��o�Õ��������š˾�$ݾ�+�Ǯ�Ǯ��7L11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
~�B
� B
� B
� B
�B
�B
�B
�PB
��B�B?}B;dBK�B@�BS�BffBw�BffBy�B{�B}�B�B�+B|�B� B�%B�PB�VB�VB�PB�PB�=B�VB�PB�oB�hB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B
~�B
� B
� B
� B
�B
�B
�B
�PB
��B�B?}B;dBK�B@�BS�BffBw�BffBy�B{�B}�B�B�+B|�B� B�%B�PB�VB�VB�PB�PB�=B�VB�PB�oB�hB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202309071610452023090716104520230907161045  IF  ARFMCODA008c                                                                20161219135546                      G�O�G�O�G�O�                IF  ARGQCOQC2.9                                                                 20161219135553  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC2.9                                                                 20161219135553  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.0 CTD2017V1                                                       20171114124206  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20180625141209  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20190314142625  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology 20191009151724  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200227203607  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20201009145903  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210531124731  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20220223174922  IP  PSAL            A+33Dx� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230907161045  IP  PSAL            A+33Dx� G�O�                