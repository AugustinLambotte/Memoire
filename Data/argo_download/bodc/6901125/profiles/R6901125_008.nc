CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   F   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
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
resolution        =���       2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  3�   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       3�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  5   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       5T   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    6l   PSAL         
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       6p   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  7�   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       7�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  8�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       90   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :H   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       :L   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  ;d   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ;�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  <�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       =   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  >$   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    >T   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    AT   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    DT   CALIBRATION_DATE            	             
_FillValue                  ,  GT   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    G�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    G�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    G�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    G�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  G�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    G�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    G�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    G�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   units         decibar    
_FillValue        G�O�        G�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    units         decibar    
_FillValue        G�O�        G�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        G�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    H ARGO profile    2.2 1.2 19500101000000  20170206005630  20170206005630  6901125 Argo UK                                                         Jon Turton                                                      PRES            PSAL            TEMP               A   BO  38177                           2B+ A   APEX-SBE 6237                                                   846 @�W���a�1   @�W���a�@R�=p��
�"vȴ9X1   ARGOS   A   @��
A�\A,��Av�HA���AŮA��B	�RB\)B1�HBE�\B[��Bn=qB�� B��B�p�B��HB��=B��)B�  B�
=B�W
B���C)C
�C��C\)C)T{C3B�C=8RCF��CP�fC[Y�Ce(�Cn�3Cx�C��C���C��
C��HC���C���C��C���C���C���C���C�b�C�g�Cѵ�C�  C��CC��\D	?\D�RD"b�D.θDG�HD`� Dy�\D�W\D���D�W
D���D�s�D��D�W\D���D�r�1111111111111111111111111111111111111111111111111111111111111111111111  @���A
=A1G�A{\)A�
=A��A�(�B
�
B z�B3  BF�B\�RBo\)B�\B���B�  B�p�B��B�k�B��\Bƙ�B��fB�aHCc�C5�C�=C��C)�)C3�=C=� CF�RCQ.C[�HCep�Co:�Cy�C�{C��qC���C�C��fC��fC���C��C��\C��C���C��fCǋ�C�ٚC�#�C���C�\C��3D	QHD�=D"t{D.�DG�3D`��Dy�HD�`RD��D�` D���D�|{D���D�`RD���D�{�1111111111111111111111111111111111111111111111111111111111111111111111  G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A   Bi�B B
R�BB|�B�B�B�B$�B'mB)_B,�B-)B)�B,"B*KB%FB#:B#�B$B!�B�B�B�BPBB�B�XB�B�aB�RB�YB�[B��B�B��B��B��B� B��B�oB��B�hB�-B��B��B��B�zB�TB�B��B�B�TB�GB��B�2B��B��B�NB�!B��B��B�CB�7B��B�B�KB��B�B��1111111111111111111111111111111111111111111111111111111111111111111111  Bi�B B
R�BB|�B�B�B�B$�B'mB)_B,�B-)B)�B,"B*KB%FB#:B#�B$B!�B�B�B�BPBB�B�XB�B�aB�RB�YB�[B��B�B��B��B��B� B��B�oB��B�hB�-B��B��B��B�zB�TB�B��B�B�TB�GB��B�2B��B��B�NB�!B��B��B�CB�7B��B�B�KB��B�B��1111111111111111111111111111111111111111111111111111111111111111111111  G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B   @�V@�\)@��@�x�@`�?@L�@MrG@Jxl@F�2@B��@>��@=}�@:��@8��@5��@0��@-4@*��@)T�@'�]@'�@![W@�n@F�@{J@s�?��?̨�?�<6?���?��]?�H�?z�1?Z��?7E9?($?}�?�6??��>�>�e>�8�>�i�>��;>�G�>���>��>C,�=��j=��=��I=F?;�t����e�澧E9��4׾�����1'���?���s��m�S��������C��/���64441144111111111111111111111111111111144111144111111111111111111111111  G�O�G�O�G�O�@�x�@`�?G�O�G�O�@Jxl@F�2@B��@>��@=}�@:��@8��@5��@0��@-4@*��@)T�@'�]@'�@![W@�n@F�@{J@s�?��?̨�?�<6?���?��]?�H�?z�1?Z��?7E9?($?}�?�6G�O�G�O�>�>�e>�8�>�i�G�O�G�O�>���>��>C,�=��j=��=��I=F?;�t����e�澧E9��4׾�����1'���?���s��m�S��������C��/���64441144111111111111111111111111111111144111144111111111111111111111111  G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            PSAL                            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE from current cycle.                                                                                                                                                                                     PSAL_ADJSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                                                                                                                                                                                                                                                                                    dP = -0.28                                                                                                                                                                                                                                                      e_time assumes 0.09dBar/s ascent rate, alpha=0.0141, tau=6.68s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                                                                                                                                                                                                                                                                                    2012082319124320120823191243                BO  ARGQPREX2.0                                                                 20120823191243  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20120823191243  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20120823191243  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20120823191243  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20120823191243  CR                  G�O�G�O�G�O�                BO  ARGQRTSP1.0                                                                 20120823191243  CV                  G�O�G�O�G�O�                BO  ARGQRTQC2.0                                                                 20120823191246  QCP$                G�O�G�O�G�O�6389758         BO  ARGQRTQC2.0                                                                 20120823191252  QCP$                G�O�G�O�G�O�6389758         BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            C���C���?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            C���C���?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            C��HC��H?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            C��
C��
?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            A��A��?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            AŮAŮ?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            A,��A,��?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            A�\A�\?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$SIGT            @��
@��
?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCP$                G�O�G�O�G�O�131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            @��
@��
?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            A�\A�\?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            A,��A,��?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            AŮAŮ?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            A��A��?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            C��
C��
?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            C��HC��H?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            C���C���?�  131072          BO  ARGQSCUT2.0                                                                 20120829113209  QCF$TEMP            C���C���?�  131072          