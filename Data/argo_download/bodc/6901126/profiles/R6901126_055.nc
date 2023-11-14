CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   D   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
_FillValue                    0�   FORMAT_VERSION                 comment       File format version    
_FillValue                    0�   HANDBOOK_VERSION               comment       Data handbook version      
_FillValue                    0�   REFERENCE_DATE_TIME                 comment       !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    1    DATE_CREATION                   comment       Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    1   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    1    PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    10   PROJECT_NAME                  comment       Name of the project    
_FillValue                  @  18   PI_NAME                   comment       "Name of the principal investigator     
_FillValue                  @  1x   STATION_PARAMETERS           	            conventions       Argo reference table 3     	long_name         ,List of available parameters for the station   
_FillValue                  0  1�   CYCLE_NUMBER               	long_name         Float cycle number     
_FillValue         ��   conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle        1�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    1�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    1�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     1�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    2   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    2   INST_REFERENCE                    	long_name         Instrument type    conventions       Brand, type, serial number     
_FillValue                  @  2   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    2\   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    conventions       8Relative julian days with decimal part (as parts of day)   units         "days since 1950-01-01 00:00:00 UTC     
_FillValue        A.�~            2`   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    2h   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            2l   LATITUDE               	long_name         &Latitude of the station, best estimate     
_FillValue        @�i�       units         degree_north   	valid_min         �V�        	valid_max         @V�             2t   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    
_FillValue        @�i�       units         degree_east    	valid_min         �f�        	valid_max         @f�             2|   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    2�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    2�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    2�   PRES         
      	   	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  3�   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       3�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  4�   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       5<   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    6L   PSAL         
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       6P   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  7`   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       7�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  8�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       8�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       :   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  ;   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ;`   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  D  <p   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       <�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  =�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    =�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    @�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    C�   CALIBRATION_DATE            	             
_FillValue                  ,  F�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    G    HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    G$   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    G(   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    G,   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  G0   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    Gp   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    G�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    G�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   units         decibar    
_FillValue        G�O�        G�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    units         decibar    
_FillValue        G�O�        G�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        G�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    G�ARGO profile    2.2 1.2 19500101000000  20131206041750  20131206041750  6901126 Argo UK                                                         Jon Turton                                                      PRES            PSAL            TEMP               7A   BO  58293                           2B+ A   APEX-SBE 6238                                                   846 @��|/7�@1   @��|/7�@@R�A�7Kǿ���l�C�1   ARGOS   A   @�G�@���A$(�Au�A�{A�ffA�{B
B�
B3��BH��B^  Bq  B�ǮB�.B�\B�G�B��qB�Q�B�8RB�{B۔{B�8RC �qC�CJ=C  C)��C3T{C=^�CF��CQ&fC\�Ce�
Co^�CyY�C�� C���C���C��{C�q�C���C���C��\C��C�� C�z�C�s3C�  Cѫ�Cۃ�C��{C�fC�\)D	J�D�\D"N�D.�
DG�)D`�Dy��D�R�D��qD�\{D��\D�_
D��
D�Uq11111111111111111111111111111111111111111111111111111111111111111111@��H@�\A(��Az�RA�z�A���A�z�B��B 
=B5  BI�
B_33Br33B�aHB�ǮB���B��HB�W
B��B���BǮB�.B���CJ=CffC�
CL�C)�{C3�HC=��CF�RCQs3C\Y�Cf#�Co��Cy�fC�fC���C��RC���C��RC��\C��3C���C��{C��fC��HC���C�&fC���C۪=C���C���C���D	^D�D"a�D.�=DG�\D`�RDy�D�\)D��
D�fD� �D�h�D� �D�_
11111111111111111111111111111111111111111111111111111111111111111111G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�F   Bq�Bp�BpoBq[Bq'Bq�Br�Br�Br|BrBr�BshBr|BsBr�Br�Bs�Br�BtTBs�BvB�B��B��B��B��B��B�?B��B�B�tB�]B�B��B�TB�%B��B�nB�nB�TB��B��B��B��B�%B�ZB��B�B�MB��B�UB�/B�/B�qB��B�0B��B�
B�fB�fB��B�2B��B�@B�TB�B�vB�44444444444444444444444444444444444444444444444444444444444444444444Bq�Bp�BpoBq[Bq'Bq�Br�Br�Br|BrBr�BshBr|BsBr�Br�Bs�Br�BtTBs�BvB�B��B��B��B��B��B�?B��B�B�tB�]B�B��B�TB�%B��B�nB�nB�TB��B��B��B��B�%B�ZB��B�B�MB��B�UB�/B�/B�qB��B�0B��B�
B�fB�fB��B�2B��B�@B�TB�B�vB�11111111111111111111111111111111111111111111111111111111111111111111G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A   ��(��녾Ұ!��,=��̾�*0��觾���
����A�����g8��?��F��%���`B���������ʾ�Z�\�v��V�M���>��,����%���	l�����7Y�	�'�hs�e���͟���X�ȴ9��|������˾;�&�,"h�'��-�h�:xl�F?�c�f�w�4��h
�� i��	����!��<6��M������y���� [��h
����4�K^�BĿ�����6z11111111111111111111111111111111111111111111111111111111111111111111��(��녾Ұ!��,=��̾�*0��觾���
����A�����g8��?��F��%���`B���������ʾ�Z�\�v��V�M���>��,����%���	l�����7Y�	�'�hs�e���͟���X�ȴ9��|������˾;�&�,"h�'��-�h�:xl�F?�c�f�w�4��h
�� i��	����!��<6��M������y���� [��h
����4�K^�BĿ�����6z11111111111111111111111111111111111111111111111111111111111111111111G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            PSAL                            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE from current cycle.                                                                                                                                                                                     PSAL_ADJSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                                                                                                                                                                                                                                                                                    dP = -0.3                                                                                                                                                                                                                                                       e_time assumes 0.09dBar/s ascent rate, alpha=0.0141, tau=6.68s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                                                                                                                                                                                                                                                                                    2013120604151820131206041518                BO  ARGQPREX2.0                                                                 20131206041518  CV                  G�O�G�O�G�O�                BO  ARGQRTSP1.0                                                                 20131206041518  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20131206041518  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20131206041518  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20131206041518  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20131206041518  CR                  G�O�G�O�G�O�                BO  ARGQRTQC2.0                                                                 20131206041519  QCP$                G�O�G�O�G�O�6389758         BO  ARGQRTQC2.0                                                                 20131206041524  QCP$                G�O�G�O�G�O�6389758         