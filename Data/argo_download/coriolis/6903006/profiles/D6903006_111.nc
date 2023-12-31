CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  G   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:23Z creation; 2023-01-04T00:23:00Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_054b      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    :�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    :�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  ;   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  ;D   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  ;�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        ;�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     <   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     <8   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    <X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         <\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    <d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            <h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           <p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           <x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    <�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    <�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    <�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    <�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        =�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     
8  =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  G�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          I    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  N<   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          O�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  T�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       U�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       [   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  `    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ah   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  f�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       g�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       l�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  r   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       sL   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 H  xh   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       y�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �D   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �H   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �L   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �P   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �T   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  ~�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                       SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  �Argo profile    3.1 1.2 19500101000000  20230104002023  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               oA   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�wwww1   @�wwww@S8D�2��q���8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 2000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �V�l   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �w:��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���%�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���8�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���H(  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����>�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����h  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��'qf  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��~�/  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���%��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��]�RL  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��� �,  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����t  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��u��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��2q�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �ę���  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���>��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��Pg(�  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���Q)V  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��9D�Z  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��[�j1  999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990 A)��A>ffAS33A`  Al��A�  A���A�33A�  A���A�  A���A���A�ffA�ffA���A�33A�33A噚A�ffA���A�33BffB��B��B  BffBffBffB33B#��B(ffB,  B.��B333B8  B;��B?33BDffBHffBT  Bg��B{��B���B���B�ffB���B�  B���B�33B�  B�  B���B뙚B�  C �C�fC	��C��C��C  C�3C#  C(L�C-33C233C7L�C<  C@�3CF  CK  CO�fCT�fCY�fC_  Cc��Ch��Cn�Cr�fCw�3C|��C��3C���C�  C�ffC�  C���C�  C�ffC���C�s3C��C���C�  C�s3C�  C�s3C��fC�ffC��fC�ffC��3C�s3C�ٚC�s3C��C�� C��3C�s3C�  CɌ�C��C�s3C���C�s3C��fC�Y�C��fC�s3C��C�s3C�  C��C��fC�L�C��fC��C�  C�ffC��fC�ffC�ٚD ,�Dl�D�fD� D  D` D� D� D
  DffD� D��DL�D� D�3D� D33D�fD�3D��DFfD� D� D�DL�D��D ��D!��D#@ D$y�D%��D&��D(33D)y�D*��D,  D-,�D.y�D/ٚD1�D2@ D3` D4��D5��D7L�D8� D9�3D;�D<FfD=y�D>�fD@3DA@ DBffDC�3DEfDF9�DGl�DH�fDI�fDK&fDLffDM�fDN��DP33DQl�DR��DS��DU,�DVy�DW��DYfDZL�D[� D\� D^  D_@ D`� Da� Dc  Dd@ Dey�Df�3Dg�3Di9�Dj�fDk��Dl��Dn&fDo` Dp��Dr3DsL�Dt�fDu�fDwfDxL�Dy�3D{` D}� D�0 D�s3D��fD�� D�0 D�l�D��3D�� D�33D�i�D�� D���D�<�D�p D��fD��D�  D�c3D��3D���D�0 D�i�D���D���D�0 D�l�D��fD��fD�&fD�l�D��3D��D�  D�l�D�� D��fD�#3D�ffD���D��3D�0 D�p D���D��D�&fD�i�D���D��3D�  D�c3D�� D�� D�)�D�ffDæfD���D�33D�l�Dȩ�D��fD�)�D�i�Dͩ�D���D�&fD�c3DҦfD�� D�0 D�c3D׬�D���D�<�D�l�Dܣ3D��3D�33D�l�D�fD��fD�,�D�s3D橚D��fD�&fD�l�D�fD��3D�,�D�ffD�fD��D�33D�y�D�I�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A)��A>ffAS33A`  Al��A�  A���A�33A�  A���A�  A���A���A�ffA�ffA���A�33A�33A噚A�ffA���A�33BffB��B��B  BffBffBffB33B#��B(ffB,  B.��B333B8  B;��B?33BDffBHffBT  Bg��B{��B���B���B�ffB���B�  B���B�33B�  B�  B���B뙚B�  C �C�fC	��C��C��C  C�3C#  C(L�C-33C233C7L�C<  C@�3CF  CK  CO�fCT�fCY�fC_  Cc��Ch��Cn�Cr�fCw�3C|��C��3C���C�  C�ffC�  C���C�  C�ffC���C�s3C��C���C�  C�s3C�  C�s3C��fC�ffC��fC�ffC��3C�s3C�ٚC�s3C��C�� C��3C�s3C�  CɌ�C��C�s3C���C�s3C��fC�Y�C��fC�s3C��C�s3C�  C��C��fC�L�C��fC��C�  C�ffC��fC�ffC�ٚD ,�Dl�D�fD� D  D` D� D� D
  DffD� D��DL�D� D�3D� D33D�fD�3D��DFfD� D� D�DL�D��D ��D!��D#@ D$y�D%��D&��D(33D)y�D*��D,  D-,�D.y�D/ٚD1�D2@ D3` D4��D5��D7L�D8� D9�3D;�D<FfD=y�D>�fD@3DA@ DBffDC�3DEfDF9�DGl�DH�fDI�fDK&fDLffDM�fDN��DP33DQl�DR��DS��DU,�DVy�DW��DYfDZL�D[� D\� D^  D_@ D`� Da� Dc  Dd@ Dey�Df�3Dg�3Di9�Dj�fDk��Dl��Dn&fDo` Dp��Dr3DsL�Dt�fDu�fDwfDxL�Dy�3D{` D}� D�0 D�s3D��fD�� D�0 D�l�D��3D�� D�33D�i�D�� D���D�<�D�p D��fD��D�  D�c3D��3D���D�0 D�i�D���D���D�0 D�l�D��fD��fD�&fD�l�D��3D��D�  D�l�D�� D��fD�#3D�ffD���D��3D�0 D�p D���D��D�&fD�i�D���D��3D�  D�c3D�� D�� D�)�D�ffDæfD���D�33D�l�Dȩ�D��fD�)�D�i�Dͩ�D���D�&fD�c3DҦfD�� D�0 D�c3D׬�D���D�<�D�l�Dܣ3D��3D�33D�l�D�fD��fD�,�D�s3D橚D��fD�&fD�l�D�fD��3D�,�D�ffD�fD��D�33D�y�D�I�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����n���J��-�Ұ!��S���J��녿ѩ���J��-��녿щ7�ѩ���hs�����
��Z�ӶF���
���Ͽӕ��Ұ!��t��ӕ���녿�녿�Z���m��p���1'��ȴ���ٿ�
=�$�/�T��>���>�"�>š�?�?=�?\�?|�?�^5?���?�j?Ł?�/?��D@�#@�h@�+@
��@t�@��@@	hs@ȴ@$�@�T@�h@/@�m@��@x�@�!@9X@z�@S�@~�@ �9?��w?�5??�I�?�ƨ?���?�$�?�z�?�7?�I�?�^5?��?�E�?�&�?���?��#?�^5?ٙ�?���?�Ĝ?��/?߾w?�V?ݑh?ӕ�?�j?�
=?ǍP?��?�?�Z?�%?�1?�bN?���?�&�?ϝ�?�X?���?°!?��?�"�?���?�K�?��?���?��\?{�m?s33?rn�?_|�?N{?O\)?J��?:�H?+?&$�?$��? �?5??"�?X?b?�?��?�?�?��?��?�
>�Q�>��>>>���>��>޸R>���>š�>���>�G�>��>�/>�;d>�
=>���>�$�>|�>z�H>e`B>\(�>Q�>D��>?|�><j>49X>��>\)=���==��=�G�=��`=���=���=�;d=�`B=��=�;d=�O�=T��<�/��o���
�ě���o��C��ě���`B�\)�8Q�P�`�ixս�7L���w��9X��j�ě���������G���l���xս�����ٽ��ٽ��ٽ�����%�t���P����� Ĝ� Ĝ��w��R�!���$�/�#�
�'0 ž5?}�@��I�^�I�^�R�W
=�Y��aG��q���vȴ�w�پ|푾��7������˾�1'��\)���Ͼ�������;d��MӾ�`B�������33��j��%�Ƨ��1'�ɺ^���;�V��V��bN�և+�ڟ��޸R��G��������y��r���������&��׾�F�����پ��H��j��p�� A���\����T�ff��9�
~��C����V���&����#��-�   �"Mӿ$���$���$Z�&�y�'+�'l��)xտ+C��,�D�-��/��0 ſ0bN�0�׿0�׿0�`�0�`�0�׿1&�0�`�0bN�/��.��.��.��/��0bN�2-�49X�5��6�6E��6E��6�5111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��n���J��-�Ұ!��S���J��녿ѩ���J��-��녿щ7�ѩ���hs�����
��Z�ӶF���
���Ͽӕ��Ұ!��t��ӕ���녿�녿�Z���m��p���1'��ȴ���ٿ�
=�$�/�T��>���>�"�>š�?�?=�?\�?|�?�^5?���?�j?Ł?�/?��D@�#@�h@�+@
��@t�@��@@	hs@ȴ@$�@�T@�h@/@�m@��@x�@�!@9X@z�@S�@~�@ �9?��w?�5??�I�?�ƨ?���?�$�?�z�?�7?�I�?�^5?��?�E�?�&�?���?��#?�^5?ٙ�?���?�Ĝ?��/?߾w?�V?ݑh?ӕ�?�j?�
=?ǍP?��?�?�Z?�%?�1?�bN?���?�&�?ϝ�?�X?���?°!?��?�"�?���?�K�?��?���?��\?{�m?s33?rn�?_|�?N{?O\)?J��?:�H?+?&$�?$��? �?5??"�?X?b?�?��?�?�?��?��?�
>�Q�>��>>>���>��>޸R>���>š�>���>�G�>��>�/>�;d>�
=>���>�$�>|�>z�H>e`B>\(�>Q�>D��>?|�><j>49X>��>\)=���==��=�G�=��`=���=���=�;d=�`B=��=�;d=�O�=T��<�/��o���
�ě���o��C��ě���`B�\)�8Q�P�`�ixս�7L���w��9X��j�ě���������G���l���xս�����ٽ��ٽ��ٽ�����%�t���P����� Ĝ� Ĝ��w��R�!���$�/�#�
�'0 ž5?}�@��I�^�I�^�R�W
=�Y��aG��q���vȴ�w�پ|푾��7������˾�1'��\)���Ͼ�������;d��MӾ�`B�������33��j��%�Ƨ��1'�ɺ^���;�V��V��bN�և+�ڟ��޸R��G��������y��r���������&��׾�F�����پ��H��j��p�� A���\����T�ff��9�
~��C����V���&����#��-�   �"Mӿ$���$���$Z�&�y�'+�'l��)xտ+C��,�D�-��/��0 ſ0bN�0�׿0�׿0�`�0�`�0�׿1&�0�`�0bN�/��.��.��.��/��0bN�2-�49X�5��6�6E��6E��6�5111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�PB	�DB	�JB	�DB	�PB	�=B	�DB	�DB	�=B	�DB	�DB	�DB	�DB	�=B	�PB	�7B	�7B	�1B	�=B	�JB	�1B	�DB	�DB	�=B	�%B	�DB	��B	��B	ȴB	��B	��B	�B	�B
$�B
aHB
_;B
q�B
gmB
��B
�'B
�XB
��B
�BVB%�B8RBZBx�B}�B�DB�VB��B��B��B��B��B�B�B�B�B�B�B�B�B�9B�^B�jB�qB�wB�}B��B�}B�wB�wB�wB�jB�jB�jB�XB�XB�LB�FB�?B�B�B�B�B�-B�^B�qB�dB�XB�^B�3B�!B�B�B�RB�jB�RB�?B��B�B�!B�dB�}B�XB�RB�RB�LB�LB�LB�9B�B��B��B�{B�uB�oB�\B�7B�DB�DB�1B�B�%B�%B�%B�%B�+B�+B�%B�%B�+B�1B�1B�1B�+B�1B�7B�=B�=B�=B�DB�\B�bB�VB�JB�JB�7B�=B�=B�DB�DB�JB�JB�PB�VB�PB�PB�VB�VB�\B�bB�bB�\B�\B�\B�bB�VB�\B�VB�\B�bB�oB�{B��B��B�{B�oB�hB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�uB�oB�uB�{B�{B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B	�PB	�DB	�JB	�DB	�PB	�=B	�DB	�DB	�=B	�DB	�DB	�DB	�DB	�=B	�PB	�7B	�7B	�1B	�=B	�JB	�1B	�DB	�DB	�=B	�%B	�DB	��B	��B	ȴB	��B	��B	�B	�B
$�B
aHB
_;B
q�B
gmB
��B
�'B
�XB
��B
�BVB%�B8RBZBx�B}�B�DB�VB��B��B��B��B��B�B�B�B�B�B�B�B�B�9B�^B�jB�qB�wB�}B��B�}B�wB�wB�wB�jB�jB�jB�XB�XB�LB�FB�?B�B�B�B�B�-B�^B�qB�dB�XB�^B�3B�!B�B�B�RB�jB�RB�?B��B�B�!B�dB�}B�XB�RB�RB�LB�LB�LB�9B�B��B��B�{B�uB�oB�\B�7B�DB�DB�1B�B�%B�%B�%B�%B�+B�+B�%B�%B�+B�1B�1B�1B�+B�1B�7B�=B�=B�=B�DB�\B�bB�VB�JB�JB�7B�=B�=B�DB�DB�JB�JB�PB�VB�PB�PB�VB�VB�\B�bB�bB�\B�\B�\B�bB�VB�\B�VB�\B�bB�oB�{B��B��B�{B�oB�hB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�uB�oB�uB�{B�{B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002023                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002300  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002300  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A)��D�I�G�O�                