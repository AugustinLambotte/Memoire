CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   \   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2022-06-24T11:05:03Z creation; 2022-08-25T12:59:51Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_050f      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    9�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    9�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :    REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :$   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :4   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  :<   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  :|   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  :�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        :�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;    DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;(   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;,   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;0   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     ;P   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     ;p   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    ;�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           ;�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    ;�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            ;�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           ;�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           ;�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    ;�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    ;�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ;�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    ;�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     �  <�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  ?�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  @   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  A�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  CP   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  C�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  E   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  F�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  F�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  HX   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  H�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  J$   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  K�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  K�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  M`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  M�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    [�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    [�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    [�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    [�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  [�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    [�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    \   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    \   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         \   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         \   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        \    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    \$   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  O,   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Ol   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Sl   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Wl   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  [lArgo profile    3.1 1.2 19500101000000  20220624110503  20220825125951  6900884 AWI                                                             Gerd ROHARDT                                                    MTIME           PRES            TEMP            PSAL               �A   IF                                  2C  D   NEMO                            230                             23-May-2012                     860 @׋�A��1   @׋�A��@S~M�0��@"�T�y�1   GPS         B   A   B   Primary sampling: discrete []                                                                                                                                                                                                                                      �\��E   �^W:�   �bd�   �d2�   �f/�f�  �hEȤ   �j�j��  �k�S   �m���  �oI���  �q3�   �qr(5@  �r�a�   �s?V@  �t�Շ�  �u�I�@  �v�k�  �v����  �x�w   �y\�%�  �zO��   �z��:@  �{���   �|���  �}��@  �~��~�  ��n^   ���r�  ����   �����  ��ax:�  ���&O   ��ޠ   ��&N!`  ��z@  ����.�  ����l�  ��y\�  ����   ���m:p  ��[f��  ���&N�  ��[��  ���s�   ��N��P  ��ƻ[  �����  ��@�tp  ��;*@  �����  ��3�a�  ���u1   ��>� �  ����,�  ��\(��  ����b�  ����X  ��#Eg�  ���u1   ��w`�  ���?�  ���b��  ��\�  ��+<�  ���b:�  ���Sp  ��Vx��  ����@  �����0  ���Q�  �����  ��fff�  ��<M^�  ��
���  ����8  ���z�  ����>�  ��fO`  ��%��t  ���>2�  ���W�  ��Ѻ��  ���H8  ��ٱ�p  ���|��  ���$h�  ���j�  ����  ���>��  ����  ���B  ���X�2  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000=���?L��@`  @�33@�33A��A>ffAH  Ak33A�  A�ffA�  A���A�ffA�ffA�  A���B��B33BffB��B!��B*��B2ffB;33BG33BX��BnffB~ffB�  B���B�ffB�  B�ffB�  B�33B�ffB�ffC��C��C�3C  C)33C3ffC=33CG�3CQ� C[  Ce��Co��Cx��C�ffC�� C�ffC�Y�C�� C�ٚC��fC�Y�C��fC�� C��3C�33C�L�C³3CǦfC�Y�C���C�ٚC���D�D	s3D��DٚD&fD"` D(�3D.��D;FfDG� DL� DL� DL��DL��DL��DL��DL��DL�3DL�3DL�3DL�3DL�311111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444>L��?ffg@fff@�ff@�ffA34A@  AI��Al��A���A�33A���A�fgA�33A�33A���A���B  B��B��B33B"  B+33B2��B;��BG��BY33Bn��B~��B�33B�  B���B�33B���B�33B�ffBי�B홙C�4C�4C��C�C)L�C3� C=L�CG��CQ��C[�Ce�gCo�gCx�4C�s3C���C�s3C�fgC���C��gC��3C�fgC��3C���C�  C�@ C�Y�C�� Cǳ3C�fgC�ٚC��gC��gD3D	y�D� D� D,�D"ffD(��D.�3D;L�DG�fDL�fG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@�dZ@��
@��;@��@�1@��@��m@��m@��@�  @�(�@�b@�1'@�1'@�(�@�9X@�1'@� �@� �@� �@� �@�(�@�(�@�1'@�Q�@�9X@�Z@�A�@�bN@�I�@�9X@�1'@�9X@�9X@� �@���@���@���@�j@�M�@�ȴ@�hs@�J@��!@�@��F@��@�x�@�&�@�j@�=q@��!@�v�@�=q@�z�@��u@��/@�|�@��P@�(�@�+@���@�7L@�t�@�n�@�Q�@��\@��m@��@�A�@~5?@u/@tz�@l�D@dZ@SdZ@Ahs@'�?�v�?��m?�C�?�?�?���?�C�?�?�"�?�?�?��H?�~�?ۅ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�dZ@��
@��;@��@�1@��@��m@��m@��@�  @�(�@�b@�1'@�1'@�(�@�9X@�1'@� �@� �@� �@� �@�(�@�(�@�1'@�Q�@�9X@�Z@�A�@�bN@�I�@�9X@�1'@�9X@�9X@� �@���@���@���@�j@�M�@�ȴ@�hs@�J@��!@�@��F@��@�x�@�&�@�j@�=q@��!@�v�@�=q@�z�@��u@��/@�|�@��P@�(�@�+@���@�7L@�t�@�n�@�Q�@��\@��m@��@�A�@~5?@u/@tz�@l�D@dZ@SdZ@Ahs@'�?�v�?��m?�C�?�?�?���?�C�?�?�"�?�?�?��H?�~�?ۅ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBM�B?}B@�B?}B?}B?}B?}B?}B>wB?}B@�B@�B@�B@�B@�B@�B?}B?}B?}B@�B?}B?}B?}B@�B@�B@�B@�B@�B@�B@�B?}B?}BA�BC�B@�BJ�BE�B,B��B�5B0!B@�BF�BK�BM�BR�BQ�BS�BYB[#BYB^5B^5BYB5?B(�BA�BXBZBbNBcTBcTBbNB_;B`BB[#B\)BK�BI�B<jB=qB5?BC�B>wB7LB1'B(�BoB�B�TB�sB�sB�mB�sB�mB�mB�mB�sB�mB�sB�sB�y11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111BM�B?pB@wB?pB?pB?sB?nB?pB>hB?sB@wB@tB@tB@tB@wB@wB?nB?nB?nB@tB?nB?nB?nB@wB@tB@tB@wB@tB@tB@tB?nB?nBA~BC�B@wBJ�BE�B+�B��B�&B0B@vBF�BK�BM�BR�BQ�BS�BYB[BYB^)B^&BYB56B(�BA~BXBZBbABcEBcFBbAB_/B`2B[B\BK�BI�B<_B=dB51BC�B>iB7=B1B(�BaB�B�DB�hG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�MTIME           PRES            TEMP            PSAL            Not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Not applicable                                                                                                                                                                                                                                                  Surface pressure = -0.1 dbar                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Not applicable                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              20220825125951202208251259512022082512595120220825125951IF  ARFMCODA050f                                                                20220624110503                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220624110711  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC5.8                                                                 20220624110711  QCF$                G�O�G�O�G�O�0000000000000100GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20220825125951  IP  PSAL            =���DL�3G�O�                