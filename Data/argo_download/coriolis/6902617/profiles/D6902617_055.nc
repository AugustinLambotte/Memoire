CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:21Z creation; 2018-10-26T19:00:47Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    6�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7T   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    7�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    7�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     7�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     88   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           8\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            8h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  @`   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  B   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  H�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  J�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  QP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Y�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  `�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b@   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  o�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  q|   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  x@   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  y�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
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
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181026190021  20181119104030  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               7A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @����t�1   @����t�@S�_�ʺ@ ����9T8   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A>ffAP  A`  Ap  A~ffA�33A�  A���A�  A���A���A�  A�  A�  A�  A�  A�  A�33A�33A�  B   B��B  BffB  B  B  B  B ffB$  B(  B+��B/��B4  B8ffB<  B@  BD  BG��BL  BPffBTffBX  B[��B_��Bc��Bh  Bl  Bp  Bs��Bw��B{��B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  C  C� C  C	� C  C��C  C� C  C� C�C� C   C"� C%  C'� C*  C,ffC.�fC1ffC4  C6� C9  C;� C>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW�CY��C\�C^� Ca  Cc� Cf  Ch� Ck  CmffCp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�s3C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C��3C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�33Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�L�C��C�� C��3C�@ C� C�� C��C�L�C���C���C��C���C��D � DfD@ D� D� D  D@ D	� D
�fD  D@ D� D� D  D@ D� D� D  D@ D� D� DfD@ D� D� D   D!@ D"� D#� D$��D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@�fDA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS9�DT� DU� DW  DX@ DYy�DZ��D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh�fDi� Dj��Dl9�Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�  D��3D�C3D��3D�� D�  D�� D�c3D�  D�� D�@ D���D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D��3D�#3D��3D�c3D�3D��fD�i�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @ff@@  @�  @�  @�  @�  A   A  A   A0  A>ffAP  A`  Ap  A~ffA�33A�  A���A�  A���A���A�  A�  A�  A�  A�  A�  A�33A�33A�  B   B��B  BffB  B  B  B  B ffB$  B(  B+��B/��B4  B8ffB<  B@  BD  BG��BL  BPffBTffBX  B[��B_��Bc��Bh  Bl  Bp  Bs��Bw��B{��B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  C  C� C  C	� C  C��C  C� C  C� C�C� C   C"� C%  C'� C*  C,ffC.�fC1ffC4  C6� C9  C;� C>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW�CY��C\�C^� Ca  Cc� Cf  Ch� Ck  CmffCp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�s3C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C��3C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�33Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�L�C��C�� C��3C�@ C� C�� C��C�L�C���C���C��C���C��D � DfD@ D� D� D  D@ D	� D
�fD  D@ D� D� D  D@ D� D� D  D@ D� D� DfD@ D� D� D   D!@ D"� D#� D$��D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@�fDA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS9�DT� DU� DW  DX@ DYy�DZ��D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh�fDi� Dj��Dl9�Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�  D��3D�C3D��3D�� D�  D�� D�c3D�  D�� D�@ D���D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D��3D�#3D��3D�c3D�3D��fD�i�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@hr�@h�@h�u@h�@h�u@h�@hbN@hbN@h�@h�@h�u@h�u@h��@h��@h�9@hĜ@h��@h��@hĜ@hĜ@hĜ@h��@h��@h��@i%@i%@i%@i%@i%@i%@h��@h��@h�u@h�`@i�@i�@i�@i�@i%@h��@h�9@h�9@h��@h�u@h�u@h�u@h�u@h�u@h�u@hbN@hQ�@hA�@hA�@hQ�@hA�@hQ�@hbN@h�@hr�@hr�@hbN@hr�@h��@h��@h�@hbN@h��@i%@i7L@i7L@i7L@i&�@i�@h�9@h�9@h��@h�9@h��@h��@h��@h�9@h��@h�@hĜ@h�`@h��@hĜ@h��@hĜ@h�`@h��@hbN@hA�@hA�@hQ�@h1'@h1'@hA�@hA�@hA�@hb@g�@h  @h  @g�@g�@g�;@g�;@g�@g�w@g�@g�w@g�w@g�@g�@g��@g�@g�w@h �@g�w@g
=@g�@f$�@f5?@f�@g
=@g+@d(�@b-@a�^@ahs@b�@bJ@c�
@b�@c�m@c��@c"�@bJ@b^5@a&�@a&�@_;d@]p�@Y�#@]/@Y��@R=q@K�@K33@J��@G�P@E��@Ep�@E�T@E`B@D��@E�@C��@AX@@�u@<z�@1��@,�j@&��@$��@!�^@�j@�P@��@~�@��@v�@��@�R@��@��@?}?��?���?��?�D?��?���@�P@dZ@O�@E�@I�@��@V@5?@�@�D@?}@�@5?@E�@{@��@�@�m@�`@b@V@��@��@  @�
@dZ?�1?�X?���?���?�?}?��;?�u?�&�?ؓu?�%?�z�?�$�?���?�=q?���?�X?��P?�/?�ȴ?��?���?�E�?��?�x�?�C�?���?���?�33?���?�-?��?��7?���?��?�ff?� �?j~�?e�?_;d?_;d?^�R?_|�?^��?\�?X��?L��?8b? Ĝ??�
>�j>� �>ڟ�>���>Ǯ>Õ�>���>��T>���>�/>�V>|�>n��>k�>gl�>hr�>aG�>V>R�>8Q�>&�y>I�>�=�=��=��=�\)=�+=�%=q��=\)<D����o�u��t��������+�#�
�0 Ž@���o��O߽�����-���
��9X�Ƨ��;d����o�C���P�$�/�333�F��R�Xb�\(��`A��bMӾbMӾe`B�fff�m�h�t�j���9��C���I����;����`��񪾒񪾓t������u��������������������(����-��Ĝ��Z���/���y����������j���پ�^5��󶾾�۾�  ��J�����+�Ǯ��1'��=q��O߾��;��z�ڟ��ݲ-�����`B���þ���V��{������H���ۿ%�J��/�$ݿ+���1'�1'��9�
~����D��׿�!����E���+��u�Q��u�����#�^5����dZ��H�^5�dZ��m��m�"ѿj����-��-�/�푿/�푿��/�p���-��-����-��R�!G��!G��!�7�"Mӿ"�\�"J�!���"J1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @hr�@h�@h�u@h�@h�u@h�@hbN@hbN@h�@h�@h�u@h�u@h��@h��@h�9@hĜ@h��@h��@hĜ@hĜ@hĜ@h��@h��@h��@i%@i%@i%@i%@i%@i%@h��@h��@h�u@h�`@i�@i�@i�@i�@i%@h��@h�9@h�9@h��@h�u@h�u@h�u@h�u@h�u@h�u@hbN@hQ�@hA�@hA�@hQ�@hA�@hQ�@hbN@h�@hr�@hr�@hbN@hr�@h��@h��@h�@hbN@h��@i%@i7L@i7L@i7L@i&�@i�@h�9@h�9@h��@h�9@h��@h��@h��@h�9@h��@h�@hĜ@h�`@h��@hĜ@h��@hĜ@h�`@h��@hbN@hA�@hA�@hQ�@h1'@h1'@hA�@hA�@hA�@hb@g�@h  @h  @g�@g�@g�;@g�;@g�@g�w@g�@g�w@g�w@g�@g�@g��@g�@g�w@h �@g�w@g
=@g�@f$�@f5?@f�@g
=@g+@d(�@b-@a�^@ahs@b�@bJ@c�
@b�@c�m@c��@c"�@bJ@b^5@a&�@a&�@_;d@]p�@Y�#@]/@Y��@R=q@K�@K33@J��@G�P@E��@Ep�@E�T@E`B@D��@E�@C��@AX@@�u@<z�@1��@,�j@&��@$��@!�^@�j@�P@��@~�@��@v�@��@�R@��@��@?}?��?���?��?�D?��?���@�P@dZ@O�@E�@I�@��@V@5?@�@�D@?}@�@5?@E�@{@��@�@�m@�`@b@V@��@��@  @�
@dZ?�1?�X?���?���?�?}?��;?�u?�&�?ؓu?�%?�z�?�$�?���?�=q?���?�X?��P?�/?�ȴ?��?���?�E�?��?�x�?�C�?���?���?�33?���?�-?��?��7?���?��?�ff?� �?j~�?e�?_;d?_;d?^�R?_|�?^��?\�?X��?L��?8b? Ĝ??�
>�j>� �>ڟ�>���>Ǯ>Õ�>���>��T>���>�/>�V>|�>n��>k�>gl�>hr�>aG�>V>R�>8Q�>&�y>I�>�=�=��=��=�\)=�+=�%=q��=\)<D����o�u��t��������+�#�
�0 Ž@���o��O߽�����-���
��9X�Ƨ��;d����o�C���P�$�/�333�F��R�Xb�\(��`A��bMӾbMӾe`B�fff�m�h�t�j���9��C���I����;����`��񪾒񪾓t������u��������������������(����-��Ĝ��Z���/���y����������j���پ�^5��󶾾�۾�  ��J�����+�Ǯ��1'��=q��O߾��;��z�ڟ��ݲ-�����`B���þ���V��{������H���ۿ%�J��/�$ݿ+���1'�1'��9�
~����D��׿�!����E���+��u�Q��u�����#�^5����dZ��H�^5�dZ��m��m�"ѿj����-��-�/�푿/�푿��/�p���-��-����-��R�!G��!G��!�7�"Mӿ"�\�"J�!���"J1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB`BB`BB`BB`BB`BBaHBaHBaHBaHB`BBaHB`BBaHB`BBaHB`BBaHB`BB`BB`BBaHB`BBaHB`BB`BBaHBaHB`BBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B^5B_;B^5B_;B^5B]/B\)BZBYBYBYB[#B[#B[#BZBYBYBYBW
BW
BT�BS�BQ�BN�BL�BG�B>wB<jB<jB:^B7LB5?B5?B5?B33B33B2-B0!B-B)�B#�B�B�BoBVBDB	7B%BB%B
=BPBPBJB	7BB��B��B�B�B��B��BDB%�B+B.B/B/B0!B49B49B2-B1'B49B5?B49B49B49B5?B6FB49B2-B1'B/B-B(�B#�B�B�B�B�B�B�B�B�BoBPB+BB��B�B�B�yB�sB�fB�`B�BB�5B�5B�;B�HB�NB�ZB�`B�ZB�TB�TB�TB�TB�NB�NB�BB�5B�#B�B��B��B��B��B��B��B��B��B��B��BĜB��B�wB�dB�^B�XB�RB�LB�FB�FB�?B�9B�3B�-B�'B�'B�'B�'B�'B�'B�!B�'B�!B�!B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   BaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB`BB`BB`BB`BB`BBaHBaHBaHBaHB`BBaHB`BBaHB`BBaHB`BBaHB`BB`BB`BBaHB`BBaHB`BB`BBaHBaHB`BBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B^5B_;B^5B_;B^5B]/B\)BZBYBYBYB[#B[#B[#BZBYBYBYBW
BW
BT�BS�BQ�BN�BL�BG�B>wB<jB<jB:^B7LB5?B5?B5?B33B33B2-B0!B-B)�B#�B�B�BoBVBDB	7B%BB%B
=BPBPBJB	7BB��B��B�B�B��B��BDB%�B+B.B/B/B0!B49B49B2-B1'B49B5?B49B49B49B5?B6FB49B2-B1'B/B-B(�B#�B�B�B�B�B�B�B�B�BoBPB+BB��B�B�B�yB�sB�fB�`B�BB�5B�5B�;B�HB�NB�ZB�`B�ZB�TB�TB�TB�TB�NB�NB�BB�5B�#B�B��B��B��B��B��B��B��B��B��B��BĜB��B�wB�dB�^B�XB�RB�LB�FB�FB�?B�9B�3B�-B�'B�'B�'B�'B�'B�'B�!B�'B�!B�!B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040302018111910403020181119104030  IF  ARFMCODA024c                                                                20181026190021                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190047  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181026190047  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104030  IP  PSAL            @ffD�i�G�O�                