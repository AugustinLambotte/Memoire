CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2016-12-19T13:55:08Z creation; 2023-09-07T12:38:35Z last update (BSH ARSQ software)    
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
_FillValue                 �  B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  Dt   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  LD   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N8   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  V   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  q`   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  y0   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  {$   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �$   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �d   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �t   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �x   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20161219135508  20230907123835  3901872 MOCCA-EU                                                        Romain Cancouet                                                 PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR                           AR2600-16FR035                  5900A00                         844 @��y'�}(1   @��u}'Ҁ@Q�5?|����x���1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 250 dbar; 10 sec sampling, 2 dbar average from 250 dbar to 1400 dbar; 10 sec sampling, 4 dbar average from 1400 dbar to 1000 dbar]                                                     A0  AI��Ad��A�33A�33A���A���A�ffAљ�A���A���B  B
��B��B  B"  B*ffB2��B:ffBB  BI��BQ33BZ  Ba��Bh��Bq33Bz  B�33B�  B���B�  B���B�ffB���B�33B�  B���B���B���B���B���B���B���B���B�  B�33B�ffBЙ�B�  B�33Bܙ�B���B�  B�ffB�  B�ffB���B�33B���C �CffC��CffC�C
L�C� CL�C�CL�C� CffC33CffC��CffC �C"ffC$��C&L�C(� C*��C,ffC.� C0��C2L�C4� C6�3C8ffC:� C<��C>L�C@ffCBffCD� CF��CH33CJL�CLL�CNL�CPffCR� CT��CV33CXL�CZL�C\ffC^ffC`ffCb� CdffCfffChffCjffClL�CnL�CpL�Cr��Ct��Cv� CxffCzL�C|��C~ffC�&fC�@ C�33C��C�@ C�Y�C�@ C�&fC�@ C�&fC��C�33C�L�C�33C�&fC�@ C�&fC��C�&fC�33C�L�C�&fC��C�&fC�33C�@ C�L�C�Y�C�33C��C��C�&fC�33C�@ C�@ C�L�C�Y�C�Y�C�Y�C�33C��C��C�&fC�&fC�&fC�33C�33C�33C�33C�33C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�33C�33C�&fC��C��C�33C�Y�C�Y�C�L�C�L�C�L�C�@ C�33C�33C�&fC��C�33C�Y�C�@ C�33C��C�33C�33C�&fC�@ C�L�C�33C�&fC�@ C�L�C�@ C�33C�@ C�L�C�33C�@ C�L�C�33C�33C�@ C�@ C�L�C�&fC�&fC�33C�&fC�&fC�&fC�L�C�@ C�33C�33C�&fC�@ C�33C��C�&fC�@ C�L�C�@ C�&fC�33C�L�C�L�C�@ C�&fD 3D ��D�D� D  D� D  D��D3D�3D  D�3DfD��D3D��D  D� D	fD	��D
3D
�fD  D�3D3D�3D3D��D&fD� D�D��D  D�3D�D��D3D��D�D��D3D�3D3D��DfD��D�D��D3D��D3D��D3D��D3D��D  D��D3D�3D�D�3D  D�3D �D � D!3D!�3D"3D"��D#�D#� D$3D$��D%�D%�3D&�D&�3D'  D'��D(�D(��D)3D)��D*�D*��D+�D+��D,�D,��D-�D-�3D.�D.��D/�D/�3D0�D0��D13D1� D23D2��D3�D3�3D4�D4�3D53D5�3D6�D6��D7�D7��D8�D8��D9  D9�3D:3D:��D;3D;�3D<3D<��D=  D=��D>3D>� D?�D?�3D@3D@��DA�DA�3DB�DB��DC3DC�3DD�DD�3DE�DE��DF3DF��DG  DG��DH�DH��DI�DI��DJ  DJ� DK3DK�3DL�DL��DM  DM��DN�DN�3DO�DO��DP�DP��DQ�DQ��DR3DR� DS3DS��DT3DT��DU3DU��DV�DV��DW�DW��DX3DX�3DY�DY�fDZ3DZ�3D[3D[�3D\  D\��D]�D]�3D^  D^�3D_  D_��D`3D`��Da�Da�3Db�Db��Dc  Dc��Dd3Dd��De  De�3Df�Df��Dg  Dg� Dh�Dh�3Di�Di��Dj�Dj��Dk�Dk�3Dl3Dl��Dm3Dm��Dn�Dn�3Do�Do��Dp  Dp��Dq�Dq��Dr�Dr�3Ds�Ds� Dt�Dt��Du�Du��Dv�Dv��Dw3Dw��Dx�Dx�3Dy�Dy�3Dz3Dz�3D{�D{��D{� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A0  AI��Ad��A�33A�33A���A���A�ffAљ�A���A���B  B
��B��B  B"  B*ffB2��B:ffBB  BI��BQ33BZ  Ba��Bh��Bq33Bz  B�33B�  B���B�  B���B�ffB���B�33B�  B���B���B���B���B���B���B���B���B�  B�33B�ffBЙ�B�  B�33Bܙ�B���B�  B�ffB�  B�ffB���B�33B���C �CffC��CffC�C
L�C� CL�C�CL�C� CffC33CffC��CffC �C"ffC$��C&L�C(� C*��C,ffC.� C0��C2L�C4� C6�3C8ffC:� C<��C>L�C@ffCBffCD� CF��CH33CJL�CLL�CNL�CPffCR� CT��CV33CXL�CZL�C\ffC^ffC`ffCb� CdffCfffChffCjffClL�CnL�CpL�Cr��Ct��Cv� CxffCzL�C|��C~ffC�&fC�@ C�33C��C�@ C�Y�C�@ C�&fC�@ C�&fC��C�33C�L�C�33C�&fC�@ C�&fC��C�&fC�33C�L�C�&fC��C�&fC�33C�@ C�L�C�Y�C�33C��C��C�&fC�33C�@ C�@ C�L�C�Y�C�Y�C�Y�C�33C��C��C�&fC�&fC�&fC�33C�33C�33C�33C�33C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�33C�33C�&fC��C��C�33C�Y�C�Y�C�L�C�L�C�L�C�@ C�33C�33C�&fC��C�33C�Y�C�@ C�33C��C�33C�33C�&fC�@ C�L�C�33C�&fC�@ C�L�C�@ C�33C�@ C�L�C�33C�@ C�L�C�33C�33C�@ C�@ C�L�C�&fC�&fC�33C�&fC�&fC�&fC�L�C�@ C�33C�33C�&fC�@ C�33C��C�&fC�@ C�L�C�@ C�&fC�33C�L�C�L�C�@ C�&fD 3D ��D�D� D  D� D  D��D3D�3D  D�3DfD��D3D��D  D� D	fD	��D
3D
�fD  D�3D3D�3D3D��D&fD� D�D��D  D�3D�D��D3D��D�D��D3D�3D3D��DfD��D�D��D3D��D3D��D3D��D3D��D  D��D3D�3D�D�3D  D�3D �D � D!3D!�3D"3D"��D#�D#� D$3D$��D%�D%�3D&�D&�3D'  D'��D(�D(��D)3D)��D*�D*��D+�D+��D,�D,��D-�D-�3D.�D.��D/�D/�3D0�D0��D13D1� D23D2��D3�D3�3D4�D4�3D53D5�3D6�D6��D7�D7��D8�D8��D9  D9�3D:3D:��D;3D;�3D<3D<��D=  D=��D>3D>� D?�D?�3D@3D@��DA�DA�3DB�DB��DC3DC�3DD�DD�3DE�DE��DF3DF��DG  DG��DH�DH��DI�DI��DJ  DJ� DK3DK�3DL�DL��DM  DM��DN�DN�3DO�DO��DP�DP��DQ�DQ��DR3DR� DS3DS��DT3DT��DU3DU��DV�DV��DW�DW��DX3DX�3DY�DY�fDZ3DZ�3D[3D[�3D\  D\��D]�D]�3D^  D^�3D_  D_��D`3D`��Da�Da�3Db�Db��Dc  Dc��Dd3Dd��De  De�3Df�Df��Dg  Dg� Dh�Dh�3Di�Di��Dj�Dj��Dk�Dk�3Dl3Dl��Dm3Dm��Dn�Dn�3Do�Do��Dp  Dp��Dq�Dq��Dr�Dr�3Ds�Ds� Dt�Dt��Du�Du��Dv�Dv��Dw3Dw��Dx�Dx�3Dy�Dy�3Dz3Dz�3D{�D{��D{� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�S�@��@���@�~�@�=q@���@֟�@ϕ�@��@Ĭ@�5?@�ȴ@��^@Q�@4�@�@1'@
J?�n�?�  ?\?���?�+?�+?��^?��\?wK�?k�?]�-?Ix�?A%?9�#?1��?)x�?"J?��?-??��?J>��m>�?}>>�`B>׍P>�ƨ>\>�%>��>�9X>��/>�/>��u>�\)>�=q>���>�J>vȴ>p��>q��>w��>�J>���>���>~��>u>s�F>m�h>q��>q��>m�h>ix�>_;d>Y�>Y�>N�>E��>B�\>?|�>H�9>Q�>M��>H�9>D��>>v�>6E�>333>/�>-V>(��>�->�u>V>$�>�>�>�>o>�>J>   =��#=��=�l�=�/=�
==���=�^5=��=�\)=�+=u=}�=�t�=��=�hs=�o=ix�=,1=�P=,1=T��=Y�=H�9=<j=,1=#�
=t�=o<�`B<���<�t�<�o<�t�<�9X<���<u<D��<D��<t�;o:�o:�o;D��    �o�D���D����o���
��o��o�D�����
��`B��`B�ě����
�o�o�D����o��o��o��o��o��o��o��o��o��o��o��o��o��o�D����o�D���D���D���D���D���D���D����o    ;D��;D��;��
;ě�<o<t�<t�<t�<t�<#�
<49X<49X<#�
<#�
<#�
<49X<#�
<#�
<t�;�`B;�`B<o<t�<t�<#�
<t�<o<o<t�<t�<t�;�`B;�`B;�`B;�`B<o<o<t�<t�<o;��
;��
;��
;��
;��
;��
;��
;�o;D��;D��;D��;o:�o            �o��o��o�D����o�ě���`B�o��`B�ě���`B�o�49X�D���u��o��o��o�u��C���t���1��j�ě�������`B��h�o�+�t����#�
�'0 Ž49X�8Q�@��H�9�T���Y��e`B�m�h�y�#��%�����7L��O߽�\)������P���㽡�����
���{��9X��j�\�Ƨ�ȴ9�������
=��/��G���l�������#�   �o�����	7L�I��V�hs�n���+��u���������-�"��%�T�(�þ+�-V�0 ž1&�333�49X�5?}�7KǾ7KǾ9X�;dZ�;dZ�<j�>vɾ?|�@��C���E�˾F��H�9�Kƨ�M��O�;�Q녾T���W
=�Xb�Z��]/�_;d�_;d�`A��cS��cS��e`B�gl��gl��ixվixվk��l�D�n���q���s�F�t�j�u�w�پw�پy�#�{�m�}󶾀���%���7���7���7���\���������$ݾ�+��7L���^������ƨ��O߾�V��\)���`��n���񪾓�Ͼ�zᾔzᾔzᾕ���
=��b�������������������"Ѿ�(�������-���-���R��;d���w��Ĝ��G����徣�
��`B��ff��l���r����þ�xվ��羪~���������1��1���D���h��������� ž� ž� ž�&龱�����!���F��?}��E���ȴ���پ�X���#��^5���H���m��푾�vɾ�|������7��o�Õ��������š˾�+�Ǯ�ȴ9��7L������C����;�V������;��bN��hs��녾�n�����t����Ͼ�z��z��������
=�׍P�׍P��b�ٙ�����ڟ���"Ѿ�"Ѿۥ�ܬ�ݲ-11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@���@�S�@��@���@�~�@�=q@���@֟�@ϕ�@��@Ĭ@�5?@�ȴ@��^@Q�@4�@�@1'@
J?�n�?�  ?\?���?�+?�+?��^?��\?wK�?k�?]�-?Ix�?A%?9�#?1��?)x�?"J?��?-??��?J>��m>�?}>>�`B>׍P>�ƨ>\>�%>��>�9X>��/>�/>��u>�\)>�=q>���>�J>vȴ>p��>q��>w��>�J>���>���>~��>u>s�F>m�h>q��>q��>m�h>ix�>_;d>Y�>Y�>N�>E��>B�\>?|�>H�9>Q�>M��>H�9>D��>>v�>6E�>333>/�>-V>(��>�->�u>V>$�>�>�>�>o>�>J>   =��#=��=�l�=�/=�
==���=�^5=��=�\)=�+=u=}�=�t�=��=�hs=�o=ix�=,1=�P=,1=T��=Y�=H�9=<j=,1=#�
=t�=o<�`B<���<�t�<�o<�t�<�9X<���<u<D��<D��<t�;o:�o:�o;D��    �o�D���D����o���
��o��o�D�����
��`B��`B�ě����
�o�o�D����o��o��o��o��o��o��o��o��o��o��o��o��o��o�D����o�D���D���D���D���D���D���D����o    ;D��;D��;��
;ě�<o<t�<t�<t�<t�<#�
<49X<49X<#�
<#�
<#�
<49X<#�
<#�
<t�;�`B;�`B<o<t�<t�<#�
<t�<o<o<t�<t�<t�;�`B;�`B;�`B;�`B<o<o<t�<t�<o;��
;��
;��
;��
;��
;��
;��
;�o;D��;D��;D��;o:�o            �o��o��o�D����o�ě���`B�o��`B�ě���`B�o�49X�D���u��o��o��o�u��C���t���1��j�ě�������`B��h�o�+�t����#�
�'0 Ž49X�8Q�@��H�9�T���Y��e`B�m�h�y�#��%�����7L��O߽�\)������P���㽡�����
���{��9X��j�\�Ƨ�ȴ9�������
=��/��G���l�������#�   �o�����	7L�I��V�hs�n���+��u���������-�"��%�T�(�þ+�-V�0 ž1&�333�49X�5?}�7KǾ7KǾ9X�;dZ�;dZ�<j�>vɾ?|�@��C���E�˾F��H�9�Kƨ�M��O�;�Q녾T���W
=�Xb�Z��]/�_;d�_;d�`A��cS��cS��e`B�gl��gl��ixվixվk��l�D�n���q���s�F�t�j�u�w�پw�پy�#�{�m�}󶾀���%���7���7���7���\���������$ݾ�+��7L���^������ƨ��O߾�V��\)���`��n���񪾓�Ͼ�zᾔzᾔzᾕ���
=��b�������������������"Ѿ�(�������-���-���R��;d���w��Ĝ��G����徣�
��`B��ff��l���r����þ�xվ��羪~���������1��1���D���h��������� ž� ž� ž�&龱�����!���F��?}��E���ȴ���پ�X���#��^5���H���m��푾�vɾ�|������7��o�Õ��������š˾�+�Ǯ�ȴ9��7L������C����;�V������;��bN��hs��녾�n�����t����Ͼ�z��z��������
=�׍P�׍P��b�ٙ�����ڟ���"Ѿ�"Ѿۥ�ܬ�ݲ-11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	VB	��B	�qB	�`B
+B
'�B
<jB
e`B
u�B
��B
�FB
ɺB
�fB
�#B,BbB{B'�B#�B<jB;dBA�BA�B?}BS�B[#B`BBbNBe`BiyBp�Bo�Br�Bt�Bv�Bx�Bz�B|�B}�B~�B�B�B�B�B�%B�=B�=B�7B�DB�JB�JB�PB�JB�JB�PB�PB�VB�\B�\B�bB�uB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B	VB	��B	�qB	�`B
+B
'�B
<jB
e`B
u�B
��B
�FB
ɺB
�fB
�#B,BbB{B'�B#�B<jB;dBA�BA�B?}BS�B[#B`BBbNBe`BiyBp�Bo�Br�Bt�Bv�Bx�Bz�B|�B}�B~�B�B�B�B�B�%B�=B�=B�7B�DB�JB�JB�PB�JB�JB�PB�PB�VB�\B�\B�bB�uB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202309071238352023090712383520230907123835  IF  ARFMCODA008c                                                                20161219135508                      G�O�G�O�G�O�                IF  ARGQCOQC2.9                                                                 20161219135516  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC2.9                                                                 20161219135516  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.0 CTD2017V1                                                       20171114115739  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20180621164100  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20190314105017  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology 20190920141310  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20201026173348  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210613120958  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20220225161333  IP  PSAL            A0  D{� G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230907123835  IP  PSAL            A0  D{� G�O�                