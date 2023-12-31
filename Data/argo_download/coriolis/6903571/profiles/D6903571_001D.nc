CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       �2022-11-16T16:26:47Z creation; 2022-11-16T16:27:59Z last update (coriolis COQC software); 2023-02-15T11:51:39ZZ DMQC performed on CORE variables; 2023-03-14T14:46:07ZZ DMQC performed on CORE variables   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.41   Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_053c      comment_dmqc_operator         iPRIMARY | http://orcid.org/0000-0003-2516-6106 | Jan Even Øie Nilsen, Institute of Marine Research (IMR)         C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    ;L   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    ;\   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    ;`   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    ;d   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    ;t   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    ;�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    ;�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  ;�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  ;�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  <   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        <\   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    <`   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    <d   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     <h   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    <�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    <�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     <�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     <�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     <�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    <�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         <�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    <�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            =    LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           =   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           =   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    =   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    =   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    =$   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    =(   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    =,   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    =0   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    =4   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        >4   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     �  >8   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  L�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        X  N�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        X  W�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _H   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     X  a    TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  hx   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  o�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  q�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  y    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  z�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  �0   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  �`   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     X  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �`   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �d   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �h   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �l   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �p   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �(   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �(   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �(   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  �(Argo profile    3.1 1.2 19500101000000  20221116162647  20230314144607  6903571 NorArgo2                                                        Kjell Arne Mork                                                 MTIME           PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR_D                         AD2700-19NO001                  5608A13                         838 @�P`�1   @�C��R�@RZ��Ͽ�HI�g��1   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 10 dbar; 10 sec sampling, 2 dbar average from 10 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 1000 dbar]                                                       ?HEȢ   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?cW���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?pN��@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?vx��   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?|(���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�)V��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��\�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�Sʠ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���8�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��#��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���
@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�j�|  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�2��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���F�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?� �.P  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�W�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�r(  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�"�@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�(3�P  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�9D�`  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�Y�S�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�y�H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��&N   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�؎�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�"�@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�d�	�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?����H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�[fǀ  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�"�9H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�|e�8  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�v�I4  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��\(�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�
���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�ZC�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��?�0  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���kX  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�e�8$  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���t  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�\,  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��J��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�Vٱ�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?� <�v  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��/hL  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�%�	|  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��`�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�d��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�W:�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�q�r  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��6�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�o��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?��}'�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�<�u�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?Ǯu�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�I��J  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�wwwx  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���8  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�^З�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?̬�X  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?���j2  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?Ϫ���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ?�v���  09999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990999999099999909999990  AVffAs33A���A�  A�33A���A�  A�33A�33A�33B  B
  BffB  B$  B*  B2��B;33BD  BL��BT  BZffBd  Bl��Bt  Bz��B���B���B���B�  B���B�33B���B�  B�ffB�  B���B�  B���B�  B���B�  B�33B�ffB���B�  B�ffB�ffBڙ�B���B���B晚BꙚBB�ffB�ffB�ffB�ffCL�CffCL�C33C	�C�C�C  C  C  C�fC��C��C� C�fCL�C!33C#33C%�C'  C(�fC*�3C-  C/L�C133C3  C4�fC6�3C9  C;33C=  C>��CA  CC33CD�fCG  CI33CJ�fCL��CN�3CP�fCR�fCU  CW�CY33CZ��C\�fC_  Ca  Cc�Ce�Cg�Ci�Ck33Cl��Cn�fCp�fCr�fCu  Cw  Cy�C{  C}�C~�3C�ffC�ffC�ffC�ffC�ffC�Y�C�Y�C�� C���C�� C�� C�� C�� C�� C���C�� C�� C�� C�s3C�s3C���C�� C�ffC���C�� C�s3C���C�s3C�ffC�� C���C�� C�s3C���C�� C�ffC�� C�s3C�Y�C�s3C���C�s3C�Y�C�Y�C�Y�C�s3C�� C�� C���C���C�s3C�L�C�Y�C�ffC�Y�C�Y�C�Y�C�Y�C�ffC�ffC�ffC�ffC�ffC�Y�C�s3C���C�CÀ C�ffCŀ Cƀ Cǀ CȀ Cɀ C�s3C�ffC̀ C͌�C΀ C�ffC�ffC�ffC�ffC�Y�C�ffC�ffC�s3C׀ C؀ C�ffC�Y�C�ffC�Y�C�ffC�s3C�ffC�ffC�ffC� C�Y�C� C� C��C��C�Y�C�Y�C�s3C�Y�C�ffC�Y�C�s3C�Y�C�s3C�s3C�s3C�ffC�ffC�s3C�ffC�ffC�� C�Y�C�Y�C�s3C�Y�C�s3C�ffC�s3D 33D �3D33D��D9�D�3D@ D��D33D� D,�D��D9�D��DL�D��DFfD� D	9�D	��D
@ D
�fD9�D�3D9�D� D@ D��D33D�3D33D�3D@ D��D9�D��D,�D��D@ D�fD@ D��D9�D� D,�D�3D9�D��D33D��D,�D�3D,�D� D33D�3D9�D� D@ D��DFfD�fD@ D� D 9�D � D!9�D!�3D"FfD"��D#,�D#�3D$9�D$�fD%9�D%�3D&,�D&�3D'33D'��D(9�D(� D)33D)�3D*9�D*�3D+9�D+�3D,&fD,�3D-9�D-��D.9�D.�3D/9�D/��D033D0�3D1,�D1��D233D2�3D333D3�3D433D4�3D5,�D5��D6,�D6�fD79�D7�fD8@ D8��D933D9��D:33D:� D;33D;� D<@ D<�3D=33D=�3D>33D>� D?9�D?� D@9�D@��DA9�DA� DB@ DB� DCFfDC�fDDFfDD��DE,�DE��DF9�DF��DG,�DG� DH@ DH�fDI,�DI��DJ,�DJ��DK9�DK� DL@ DL� DM33DM��DN9�DN��DO9�DO�3DP,�DP��DQ,�DQ��DR@ DR� DS@ DS��DT9�DT��DU9�DU��DV33DV��DW9�DW�fDX&fDX�3DY33DY��DZ@ DZ�3D[9�D[��D\&fD\��D]9�D]�3D^9�D^�3D_9�D_� D`33D`� Da,�Da�3Db33Db�3Dc,�Dc��Dd@ Dd� De,�De�3Df33Df��Dg9�Dg�3Dh33Dh��Di33Di��Dj9�Dj��Dk9�Dk�3Dl33Dl��Dm,�Dm��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  AVffAs33A���A�  A�33A���A�  A�33A�33A�33B  B
  BffB  B$  B*  B2��B;33BD  BL��BT  BZffBd  Bl��Bt  Bz��B���B���B���B�  B���B�33B���B�  B�ffB�  B���B�  B���B�  B���B�  B�33B�ffB���B�  B�ffB�ffBڙ�B���B���B晚BꙚBB�ffB�ffB�ffB�ffCL�CffCL�C33C	�C�C�C  C  C  C�fC��C��C� C�fCL�C!33C#33C%�C'  C(�fC*�3C-  C/L�C133C3  C4�fC6�3C9  C;33C=  C>��CA  CC33CD�fCG  CI33CJ�fCL��CN�3CP�fCR�fCU  CW�CY33CZ��C\�fC_  Ca  Cc�Ce�Cg�Ci�Ck33Cl��Cn�fCp�fCr�fCu  Cw  Cy�C{  C}�C~�3C�ffC�ffC�ffC�ffC�ffC�Y�C�Y�C�� C���C�� C�� C�� C�� C�� C���C�� C�� C�� C�s3C�s3C���C�� C�ffC���C�� C�s3C���C�s3C�ffC�� C���C�� C�s3C���C�� C�ffC�� C�s3C�Y�C�s3C���C�s3C�Y�C�Y�C�Y�C�s3C�� C�� C���C���C�s3C�L�C�Y�C�ffC�Y�C�Y�C�Y�C�Y�C�ffC�ffC�ffC�ffC�ffC�Y�C�s3C���C�CÀ C�ffCŀ Cƀ Cǀ CȀ Cɀ C�s3C�ffC̀ C͌�C΀ C�ffC�ffC�ffC�ffC�Y�C�ffC�ffC�s3C׀ C؀ C�ffC�Y�C�ffC�Y�C�ffC�s3C�ffC�ffC�ffC� C�Y�C� C� C��C��C�Y�C�Y�C�s3C�Y�C�ffC�Y�C�s3C�Y�C�s3C�s3C�s3C�ffC�ffC�s3C�ffC�ffC�� C�Y�C�Y�C�s3C�Y�C�s3C�ffC�s3D 33D �3D33D��D9�D�3D@ D��D33D� D,�D��D9�D��DL�D��DFfD� D	9�D	��D
@ D
�fD9�D�3D9�D� D@ D��D33D�3D33D�3D@ D��D9�D��D,�D��D@ D�fD@ D��D9�D� D,�D�3D9�D��D33D��D,�D�3D,�D� D33D�3D9�D� D@ D��DFfD�fD@ D� D 9�D � D!9�D!�3D"FfD"��D#,�D#�3D$9�D$�fD%9�D%�3D&,�D&�3D'33D'��D(9�D(� D)33D)�3D*9�D*�3D+9�D+�3D,&fD,�3D-9�D-��D.9�D.�3D/9�D/��D033D0�3D1,�D1��D233D2�3D333D3�3D433D4�3D5,�D5��D6,�D6�fD79�D7�fD8@ D8��D933D9��D:33D:� D;33D;� D<@ D<�3D=33D=�3D>33D>� D?9�D?� D@9�D@��DA9�DA� DB@ DB� DCFfDC�fDDFfDD��DE,�DE��DF9�DF��DG,�DG� DH@ DH�fDI,�DI��DJ,�DJ��DK9�DK� DL@ DL� DM33DM��DN9�DN��DO9�DO�3DP,�DP��DQ,�DQ��DR@ DR� DS@ DS��DT9�DT��DU9�DU��DV33DV��DW9�DW�fDX&fDX�3DY33DY��DZ@ DZ�3D[9�D[��D\&fD\��D]9�D]�3D^9�D^�3D_9�D_� D`33D`� Da,�Da�3Db33Db�3Dc,�Dc��Dd@ Dd� De,�De�3Df33Df��Dg9�Dg�3Dh33Dh��Di33Di��Dj9�Dj��Dk9�Dk�3Dl33Dl��Dm,�Dm��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��>�dZ>�dZ>�^5>�^5>���>���>���>���>�Q�>���>���>�K�>�K�>���>�K�>�ȴ>�>�9X>�33>��!>��!>��!>�->�->�->�->��!>��!>��!>��F>�33>��!>�33>��F>�>�ȴ>�>�E�>�?}>�E�>�Q�>��#>�^5>��H>�dZ>�dZ>�dZ>�dZ>��H>��#>�X>�X>�^5>�^5>��H>��H>��#>�K�>�E�>�?}>�>��j>�33>�->���>�&�>�->��!>���>� �>���>���>��>���>���>��D>��>�ff>���>�(�>���>�I�>���>y�#>_;d>=p�>��=��=��=� �=q��=#�
=C�<�9X    �T����1��`B��h�C��'0 ŽaG��y�#��7L����������
=��xս��F�J�	7L�\)����u������ Ĝ�%�T�'+�.{�/��0 ž0 ž2-�2-�2-�2-�2-�2-�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�2-�333�49X�5?}�6E��7KǾ9X�9X�:^5�;dZ�>vɾ?|�A�7�B�\�C���B�\�D���D���E�˾F��G��H�9�I�^�L�;M��O�;�O�;�Q녾R�S�ϾT���W
=�]/�aG��e`B�gl��j~��m�h�o���s�F�t�j�u�vȴ�w�پw�پw�پx���x���y�#�y�#�z�H�z�H�z�H�{�m�z�H�{�m�{�m�{�m�{�m�|푾|푾|푾|푾|푾|푾{�m�y�#�x���vȴ�u�t�j�t�j�t�j�s�F�s�F�s�F�s�F�s�F�r�!�r�!�s�F�t�j�r�!�vȴ�y�#�z�H�{�m�{�m�|푾{�m�}�}�~�۾�  �����������%��J��o�������˾���+��+��+�����+�����+�����C�������ƨ��ƨ��C���ƨ��I���O߾�O߾��;��;��;��;��;�O߾�O߾�����V��n���zᾕ���b���u��b���u��b��b��b��b���+��������zᾔzᾑ녾�t���������������Ͼ��Ͼ��Ͼ���������+���P�����"Ѿ��㾛�㾜(���(���(���(���(���(���(����㾛�㾜(����㾛�㾛�㾛�㾛"Ѿ�"Ѿ�����-��5?���R���R���R��;d���R���R��5?���R��;d���w��A���A����R��5?��5?���-�����(����㾛"Ѿ�������������������u���u���u���u��b��b���P���P��
=��
=��
=��
=��
=����������������������������������������������������������������zᾔzᾔzᾔzᾔzᾔ����������+��
=��
=��
=��
=��
=��
=���P���P���P��b��b��b���P���P���P��b���u�������������㾝/���-���-��5?��5?��5?���-��5?��5?���R���R���R��5?���R��5?��5?��5?��5?��5?��5?��5?���-���-�����/��/�����/������������㾜���(���(���(������������������������/�����(���(���(����㾛�㾛�㾛�㾛"Ѿ�"Ѿ�"�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  >�dZ>�dZ>�^5>�^5>���>���>���>���>�Q�>���>���>�K�>�K�>���>�K�>�ȴ>�>�9X>�33>��!>��!>��!>�->�->�->�->��!>��!>��!>��F>�33>��!>�33>��F>�>�ȴ>�>�E�>�?}>�E�>�Q�>��#>�^5>��H>�dZ>�dZ>�dZ>�dZ>��H>��#>�X>�X>�^5>�^5>��H>��H>��#>�K�>�E�>�?}>�>��j>�33>�->���>�&�>�->��!>���>� �>���>���>��>���>���>��D>��>�ff>���>�(�>���>�I�>���>y�#>_;d>=p�>��=��=��=� �=q��=#�
=C�<�9X    �T����1��`B��h�C��'0 ŽaG��y�#��7L����������
=��xս��F�J�	7L�\)����u������ Ĝ�%�T�'+�.{�/��0 ž0 ž2-�2-�2-�2-�2-�2-�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�1&�2-�333�49X�5?}�6E��7KǾ9X�9X�:^5�;dZ�>vɾ?|�A�7�B�\�C���B�\�D���D���E�˾F��G��H�9�I�^�L�;M��O�;�O�;�Q녾R�S�ϾT���W
=�]/�aG��e`B�gl��j~��m�h�o���s�F�t�j�u�vȴ�w�پw�پw�پx���x���y�#�y�#�z�H�z�H�z�H�{�m�z�H�{�m�{�m�{�m�{�m�|푾|푾|푾|푾|푾|푾{�m�y�#�x���vȴ�u�t�j�t�j�t�j�s�F�s�F�s�F�s�F�s�F�r�!�r�!�s�F�t�j�r�!�vȴ�y�#�z�H�{�m�{�m�|푾{�m�}�}�~�۾�  �����������%��J��o�������˾���+��+��+�����+�����+�����C�������ƨ��ƨ��C���ƨ��I���O߾�O߾��;��;��;��;��;�O߾�O߾�����V��n���zᾕ���b���u��b���u��b��b��b��b���+��������zᾔzᾑ녾�t���������������Ͼ��Ͼ��Ͼ���������+���P�����"Ѿ��㾛�㾜(���(���(���(���(���(���(����㾛�㾜(����㾛�㾛�㾛�㾛"Ѿ�"Ѿ�����-��5?���R���R���R��;d���R���R��5?���R��;d���w��A���A����R��5?��5?���-�����(����㾛"Ѿ�������������������u���u���u���u��b��b���P���P��
=��
=��
=��
=��
=����������������������������������������������������������������zᾔzᾔzᾔzᾔzᾔ����������+��
=��
=��
=��
=��
=��
=���P���P���P��b��b��b���P���P���P��b���u�������������㾝/���-���-��5?��5?��5?���-��5?��5?���R���R���R��5?���R��5?��5?��5?��5?��5?��5?��5?���-���-�����/��/�����/������������㾜���(���(���(������������������������/�����(���(���(����㾛�㾛�㾛�㾛"Ѿ�"Ѿ�"�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�!B�B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�'B�'B�'B�'B�'B�'B�'B�!B�'B�'B�'B�'B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�B�B�!B�!B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�!B�B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�'B�'B�'B�'B�'B�'B�'B�!B�'B�'B�'B�'B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�!B�B�B�!B�!B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES.                                                                                                                                                                                                                                           TEMP_ADJUSTED = TEMP.                                                                                                                                                                                                                                           PSAL_ADJUSTED = PSAL.                                                                                                                                                                                                                                           not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  The quoted error is manufacturer specified accuracy.                                                                                                                                                                                                            The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity offset or drift detected. The quoted error is max[0.01, statistical uncertainty] in PSS-78.                                                                                                                                             20230314102334202303141023342023031410233420230314102334IF  ARFMCODA053c                                                                20221116162647                      G�O�G�O�G�O�                IF  ARGQCOQC5.9                                                                 20221116162759  QCP$                G�O�G�O�G�O�000000000288F37EIF  ARGQCOQC5.9                                                                 20221116162759  QCF$                G�O�G�O�G�O�0000000000000000