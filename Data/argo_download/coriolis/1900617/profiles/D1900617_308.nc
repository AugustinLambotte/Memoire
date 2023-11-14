CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   Q   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      history       X2014-10-15T18:37:46Z creation; 2019-06-03T17:53:07Z last update (coriolis COFC software)   comment       bThe profile number used to assign the CONFIG_MISSION_NUMBER has not been check against ANDRO data.        @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        D  :    PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ;D   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        D  ;�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  <�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     D  =0   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  >t   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  ?�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  @   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  AP   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  A�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  B�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  D,   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  D�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  T  E�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     D  F   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  G\   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    G�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    J�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    M�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  P�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    P�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    P�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    P�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    P�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  P�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    Q   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    Q   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    Q   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Q,   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Q0   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Q4   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Q8Argo profile    3.1 1.2 19500101000000  20141015183746  20190603175307  1900617 CONGAS                                                          Alain SERPETTE                                                  PRES            PSAL            TEMP              4A   IF  34441454                        2C  D   APEX                            2573                            n/a                             846 @����t1   @����t@Qz�G��%��x���1   ARGOS   Primary sampling: discrete []                                                                                                                                                                                                                                      A   A   A   @�  A��AfffA�  A���A�33B
ffB��B2ffBFffBY33Bm��B�  B�33B�  B�33B�33B���B�  B�33Bڙ�BC� C��CL�CffC)��C3ffC=��CG� CQ  C[ffCe33CoffCy�C�� C�s3C�� C��3C���C�� C��3C��fC�� C�� C���C�� C���C³3Cǳ3C�� Cь�Cֳ3C۳3C�� C��C�3CC��C��fC�� D� DY�D� D	9�D�fDY�D�fD@ D� D@ D��DY�D� D"33D$�3D'S3D)��D,9�D.�3D/��222222222222222222222222222222222222222222222222222222222222222222222222222222222   �L��@���A��Ac34A�33A���A�ffB��B33B133BD  BXffBl��B34B�ffB���B���B�33B�ffB���B�  B�  B�ffCL�C  C�C$L�C.�C8L�CB33CK�3CV�C_�fCj�Cs��C}�4C���C�ٚC��C��fC�ٚC��C�  C��C��C��fC��C��fC��C��C��C��fC��C��C��C��fC��C��3C��fC�  C�ٚD ��DfD��D�fD
s3DfDs3D��Dl�D��Dy�DfDl�D � D#� D&  D(Y�D*�fD-� D.ff222222222222222222222222222222222222222222222222222222222222222222222222222222222   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B
�XB
�XB
�^B
�wB
�qB
�jB
�^B
�RB
�B
�7B
�)B
�NB
�mB
�B
�B
��B
��B%B1BPB�B%�B0!B>wBJ�BbNBS�BZBcTBq�Bx�B|�B~�B�B�B�%B�+B�1B�1B�+B�1B�1B�+B�1B�1B�1B�1B�1B�1B�7B�=B�JB�DB�=B�DB�DB�=B�DB�DB�DB�DB�DB�JB�PB�PB�PB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�bB�\B�\B�bB�\111111111111111111111111111111111111111111111111111111111111111111111111111111111   B
�XB
�XB
�^B
�wB
�qB
�jB
�^B
�RB
�B
�7B
�)B
�NB
�mB
�B
�B
��B
��B%B1BPB�B%�B0!B>wBJ�BbNBS�BZBcTBq�Bx�B|�B~�B�B�B�%B�+B�1B�1B�+B�1B�1B�+B�1B�1B�1B�1B�1B�1B�7B�=B�JB�DB�=B�DB�DB�=B�DB�DB�DB�DB�DB�JB�PB�PB�PB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�bB�\B�\B�bB�\111111111111111111111111111111111111111111111111111111111111111111111111111111111   <��
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
<��
<��
<��
<��
<��
<��
<��
@�X@�7L@١�@�5?@�$�@�@ّh@֏\@��m@�o@`b?��j?��P?kƨ?U�?>�R?#S�?t�?��>ȴ9>���>��9>"��>C��>C��>�ƨ>$�/<u<��
==��`=�%=49X<�1<u<o;D���e`B��w�ixսy�#���P���`��`B���#���\)��+���"��)����)��B�\�L�;R�cS��k��r�!�s�F�vȴ�{�m�~�۾�o���9������\)��t�������5?��MӾ��1�������!��X���\��=q��V���;111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�X@�7L@١�@�5?@�$�@�@ّh@֏\@��m@�o@`b?��j?��P?kƨ?U�?>�R?#S�?t�?��>ȴ9>���>��9>"��>C��>C��>�ƨ>$�/<u<��
==��`=�%=49X<�1<u<o;D���e`B��w�ixսy�#���P���`��`B���#���\)��+���"��)����)��B�\�L�;R�cS��k��r�!�s�F�vȴ�{�m�~�۾�o���9������\)��t�������5?��MӾ��1�������!��X���\��=q��V���;111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            PSAL            TEMP            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted for offset by using surface pressure, following the DM pressure adjustment procedure described in the Argo quality control manual; No significant pressure drift detected - Calibration error is manufacturer specified accuracy              No correction - Method OW : Weighted least squares - Error = maximum [ statistical uncertainty, 0.01]                                                                                                                                                           No significant temperature drift detected - Calibration error is manufacturer specified accuracy                                                                                                                                                                201709071532282017090715322820170907153229  IF  ARGQCOAR1.0                                                                 20141015072639  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20141015072639  QCF$                G�O�G�O�G�O�008000          IF  ARGQCOAR1.0                                                                 20141015072639  QCC$                G�O�G�O�G�O�008000          IF  CORTCOOA5.2 RTQCGL01                                                        20141015181747  QCP$PSAL            G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20141015181421  QCP$TEMP            G�O�G�O�G�O�                IF  ARSQOW  1.0 CTD2016V1                                                       20170907153229  IP  PSAL            @�  D/��G�O�                IF      COFC3.2                                                                 20190603175307                      G�O�G�O�G�O�                