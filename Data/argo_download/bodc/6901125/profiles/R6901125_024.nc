CDF   (   
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   G   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
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
resolution        =���       2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  3�   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       3�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  5   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   units         decibar    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       5\   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    6x   PSAL         
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       6|   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  7�   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       7�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  8�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     units         psu    
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       9D   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :`   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       :d   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  ;�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ;�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  <�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   units         degree_Celsius     
_FillValue        G�O�   comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       =,   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  >H   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    >x   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Ax   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Dx   CALIBRATION_DATE            	             
_FillValue                  ,  Gx   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    G�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    G�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    G�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    G�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  G�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    G�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    H   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    H   HISTORY_START_PRES                    	long_name          Start pressure action applied on   units         decibar    
_FillValue        G�O�        H   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    units         decibar    
_FillValue        G�O�        H   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        H    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    H$ARGO profile    2.2 1.2 19500101000000  20170206004210  20170206004210  6901125 Argo UK                                                         Jon Turton                                                      PRES            PSAL            TEMP               A   BO  46038                           2B+ A   APEX-SBE 6237                                                   846 @��@��1   @��@��@R�n��P��K]�c�1   ARGOS   A   @��@�p�A�RAi��A�=qA�p�A��B�
B#{B5=qBH�B\�BpG�B�{B�33B���B�B�B�ǮB���B�\)B�8RB�  B�qC��C=qC��C��C)�RC3�
C=!HCG��CQ��C[ffCe��CoQ�Cx�C�eC���C���C��\C��C��C��\C��C���C�c�C���C�k�C�G�CѥC���C��C�˅C���D	L)D��D"mqD.�RDGǮD`�Dy��D�a�D��{D�d)D���D�[�D��fD�YHD��D�o�D��11111111111111111111111111111111111111111111111111111111111111111111111 @���@��HA!p�AlQ�A���A���A�G�B�B#B5�BI33B]��Bp��B�k�B��=B�  B���B��B�(�B��3BƏ\B�W
B�{C�{Ch�C�
C�HC)��C4�C=L�CG�
CQ�RC[��Ce�RCo}qCy)C�z�C��C�˅C��C��qC���C��C���C��{C�y�C��\C��HC�]qCѺ�C��C�fC��HC��HD	W
D��D"xRD.�3DGҏDa  Dy�{D�g
D���D�i�D��RD�`�D��D�^�D��D�uD���11111111111111111111111111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B   B
��B
�NB
ԕB
�QB
�yB
רB
֡B
��B
��B
�B
�B)yBOBS�Bc�B}�B�PB�HB�B�KB�B��B��B��B�B��B��B��B��B�:B��B�B��B�B�sB�8B�$B��B��B�@B�`B�LB�B�B�sB��B��B��B��B�eB��B�
B�fB��B�TB��B��B��B��B�B�IB��B�B��B�KB�B�1B�kB�WB��B�511111114111114444144444141111111111111111111111111111111111111111111111 B
��B
�NB
ԕB
�QB
�yB
רB
֡B
��B
��B
�B
�B)yBOBS�Bc�B}�B�PB�HB�B�KB�B��B��B��B�B��B��B��B��B�:B��B�B��B�B�sB�8B�$B��B��B�@B�`B�LB�B�B�sB��B��B��B��B�eB��B�
B�fB��B�TB��B��B��B��B�B�IB��B�B��B�KB�B�1B�kB�WB��B�511111111111111111111111111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B   ��~��8����p���o�Ξ�����ԕ��Y����.�޾����f_پpoi�O�����>w�k;o>��U>H�=k�=���;��4>���>�=q>�u�=e`B��p��L0��+Խۋ��%+Ծs�=[��=��=��=��=���:�IR��j��!���3���Ƽ�1��
ں���� �I� �	��o ��X������6�F�ZkQ�u?}������� ������������F�a|�	*0���mƿ*��~�������替!hs�&�11111114111114444114444111111111111111111111111111111111111111111111111 ��~��8����p���o�Ξ�����ԕG�O����.�޾����f_پpoiG�O�G�O�G�O�G�O�>��U>H�G�O�G�O�G�O�G�O�>�=q>�u�=e`B��p��L0��+Խۋ��%+Ծs�=[��=��=��=��=���:�IR��j��!���3���Ƽ�1��
ں���� �I� �	��o ��X������6�F�ZkQ�u?}������� ������������F�a|�	*0���mƿ*��~�������替!hs�&�11111114111114444114444111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            PSAL                            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE from current cycle.                                                                                                                                                                                     PSAL_ADJSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                                                                                                                                                                                                                                                                                    dP = -0.17                                                                                                                                                                                                                                                      e_time assumes 0.09dBar/s ascent rate, alpha=0.0141, tau=6.68s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                                                                                                                                                                                                                                                                                    2013031414482620130314144826                BO  ARGQRTSP1.0                                                                 20130314144826  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20130314144826  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20130314144826  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20130314144826  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20130314144826  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20130314144826  CR                  G�O�G�O�G�O�                BO  ARGQRTQC2.0                                                                 20130314144837  QCP$                G�O�G�O�G�O�6389758         BO  ARGQRTQC2.0                                                                 20130314144858  QCP$                G�O�G�O�G�O�6389758         BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�{B�{?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�
B�
?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�qB�q?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�33B�33?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B���B���?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�B�B�B�?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B���B���?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�\)B�\)?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�8RB�8R?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�  B�  ?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            C=qC=q?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$PSAL            B�qB�q?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCP$                G�O�G�O�G�O�131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�qB�q?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�
B�
?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�8RB�8R?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�\)B�\)?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B���B���?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�B�B�B�?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B���B���?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�33B�33?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�{B�{?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�
B�
?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�{B�{?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�  B�  ?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�33B�33?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B���B���?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�B�B�B�?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�\)B�\)?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$TEMP            B�8RB�8R?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            C=qC=q?�  131072          BO  ARGQSCUT2.0                                                                 20130320104352  QCF$SIGT            B�  B�  ?�  131072          