CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   L   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2012-08-16T14:42:00Z creation; 2019-06-03T17:52:29Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7(   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    78   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7<   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7@   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7P   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7`   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7p   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7x   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8(   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8,   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    80   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     84   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8T   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8X   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8\   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8|   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T      
resolution        >�EȠ�Q)        8�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�EȠ�Q)        8�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        0  :    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  ;0   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        0  ;|   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  <�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     0  <�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  >(   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  ?X   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  ?�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  @�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  A    TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  BP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  C�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  C�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  L  D�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     0  EH   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  Fx   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    F�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    I�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    L�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  O�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    O�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    O�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    O�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    O�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  O�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    P$   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    P4   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    P8   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         PH   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         PL   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        PP   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    PTArgo profile    3.1 1.2 19500101000000  20120816144200  20190603175229  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            TEMP            PSAL               �A   IF  26871679                        2C  D   APEX                            2573                            n/a                             846 @�VQk��~1   @�VR5~@R�n��P�&�&�x��1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   B   A   @�ffAffAk33A�  A�33A���B
  BffB1��BD��BXffBl��B���B�33B���B�  B���B�ffB�  B�  B���B�ffC�C33CL�CffC)� C3  C=ffCGL�CQ33C[33CeffCo�Cy� C�� C��3C���C���C��fC��fC�� C���C��3C��3C���C���C��fC�C�� C̦fCљ�C֙�CۦfC�� C�fC��C��C���C���C�� D�fD@ D��D	9�D�fD,�D� DS3D�3DL�D� D,�DٚD"S3D$ٚ2222222222222222222222222222222222222222222222222222222222222222222222222222����@���A��AfffA�ffA�  A�34B
  B33B0ffBD  BXffBm34B�  B�ffB���B���B�33B���B���BЙ�B�33B�  C�C33CL�C$ffC-�fC8L�CB33CL�CV�C`L�Cj  CtffC~ffC�&fC�@ C�  C��C��C�33C�  C�&fC�&fC��C��C��C�  C�33C��C��C��C��C�33C��C�  C�  C��C��C�33D � D��Ds3D�3D
� D�fDy�D�Dl�DfDy�D�fD�3D!�D#�32222222222222222222222222222222222222222222222222222222222222222222222222222@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�p�B��BB	+B	��B
bB
�'Br�B��B��B�TB��BVB\B{B�B �B�B�BoBJB	7BB��B��B�B�yB�ZB�/B�B��BB�LB�LB�LB�9B�-B�?B�-B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B��B��B��1411111111111111111111111111111111111111111111111111111111111111111111111111A�p�G�O�BB	+B	��B
bB
�'Br�B��B��B�TB��BVB\B{B�B �B�B�BoBJB	7BB��B��B�B�yB�ZB�/B�B��BB�LB�LB�LB�9B�-B�?B�-B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B��B��B��1411111111111111111111111111111111111111111111111111111111111111111111111111<��
G�O�<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
<��
=D��@�r�@��@���@���@Ol�@J�\@hbN@V��@O��@G��@E�@G;d@A�7@:~�@5��@0��@-�T@';d@!%@Q�@|�@
�H@?�C�?�"�?ݲ-?�ȴ?��;?�J?���?�M�?nV?j=q?a%?O�;?F�y?H��?:^5?'+?ƨ>��>�>���>�x�>��D>��/>�Ĝ>m�h>_;d>bM�>F��>t�j>49X>.{>6E�>o=�9X=��w=�v�=�E�=u<D��<�`B=��w=�%=���=��T=�t�=T��=#�
<��
<49X�e`B�o�H�91111111111111111111111111111111111111111111111111111111111111111111111111111=D��@�r�@��@���@���@Ol�@J�\@hbN@V��@O��@G��@E�@G;d@A�7@:~�@5��@0��@-�T@';d@!%@Q�@|�@
�H@?�C�?�"�?ݲ-?�ȴ?��;?�J?���?�M�?nV?j=q?a%?O�;?F�y?H��?:^5?'+?ƨ>��>�>���>�x�>��D>��/>�Ĝ>m�h>_;d>bM�>F��>t�j>49X>.{>6E�>o=�9X=��w=�v�=�E�=u<D��<�`B=��w=�%=���=��T=�t�=T��=#�
<��
<49X�e`B�o�H�91111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           201709071514272017090715142720170907151426  IF  CORTCOOA5.2 RTQCGL01                                                        20120817025552  QCP$TEMP            G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20120817030251  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO1.4                                                                 20120817155939  QC                  G�O�G�O�G�O�                IF  ARGQCOAR1.0                                                                 20120816162325  QCP$                G�O�G�O�G�O�DEBFC           IF  ARGQCOAR1.0                                                                 20120816162325  QCF$                G�O�G�O�G�O�08200           IF      SCOO1.4                                                                 20121226104120  QC                  G�O�G�O�G�O�                IF  ARSQOW  1.0 CTD2016V1                                                       20170907151427  IP  PSAL            @�ffD$ٚG�O�                IF      COFC3.2                                                                 20190603175229                      G�O�G�O�G�O�                