CDF   	   
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
_FillValue                    H$ARGO profile    2.2 1.2 19500101000000  20170206003228  20170206003228  6901125 Argo UK                                                         Jon Turton                                                      PRES            PSAL            TEMP               �A   BO  75120                           2B+ A   APEX-SBE 6237                                                   846 @�U���� 1   @�U���� @SKI�^@#t�j~�1   ARGOS   A   @��A�HA:ffA�Q�A�
=A�Q�A��B�HB#B8z�BN��Bb
=Bs��B�33B���B�Q�B�8RB�33B�W
B��3B�=qB޽qB�33C  C�{C�C!
=C+�C4��C>�RCIL�CS��C\�
Cg��CqL�C{�C�ٚC��C�~�C��)C�~�C���C��HC��HC�}qC���C���C���CȨ�C��C��{C氤C�C���D	��D@�D"��D/NDH1HDap�Dz^D���D�%qD���D��D���D�)�DԮD� RD��
D�R11111111111111111111111111111111111111111111111111111111111111111111111 @��
A
=A>�\A�ffA��A�ffA�33B�B$��B9�BP  Bc{Bt�
B��RB�W
B��
B��qB��RB��)B�8RB�B�B�B�RCB�C�
C�C!L�C+ǮC4�C?:�CI�\CS�\C\ٚCg�RCq�\C{ǮC���C��\C�� C��qC�� C���C���C���C���C���C��
C���C��=C��C���C���C���C��)D	�{DQHD"�fD/^�DHA�Da�HDzn�D��)D�-�D��=D�D��=D�2=DԶfD�(�D��\D��11111111111111111111111111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A   B^�B^�B^jB^5B]~B^OB]�B_�B_;Ba�B^B_!B\)B\CB]�B]�B\xBY�BV9BT�BU2BS�BP�BP}BQ�BR�BR�BQ�BN�BL�BLJBK�BL~BN�BN�BNBMjBP.BOBBOBBNVBOvBM�BO\BN�BL�BM6BP.BN�BK�BE�B>�B2�B+B�B�+B��B�AB��B�BB�B��B�KB��B��B�EB��B�B��B�=B�W11111111111111111111111111111111111111111111111111111111111111111111111 B^�B^�B^jB^5B]~B^OB]�B_�B_;Ba�B^B_!B\)B\CB]�B]�B\xBY�BV9BT�BU2BS�BP�BP}BQ�BR�BR�BQ�BN�BL�BLJBK�BL~BN�BN�BNBMjBP.BOBBOBBNVBOvBM�BO\BN�BL�BM6BP.BN�BK�BE�B>�B2�B+B�B�+B��B�AB��B�BB�B��B�KB��B��B�EB��B�B��B�=B�W11111111111111111111111111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A   @�y>@�h�@�q@�z@���@�]d@�i�@��@���@�N�@�3�@�v`@�2�@�YK@���@�hs@���@�
@}��@y�j@x�@t�|@p��@l��@i�@j�!@h��@go@c�4@`��@^�@\��@[��@Z��@Y�C@X�@W'�@W�Q@W�@V��@U��@Uo @T�O@TɆ@T_@R�@R�F@R��@Q�N@O{J@J��@C]�@:+k@2�'@YK?��6?�8�?h1'>��8���辂J�ǻ0��Y���6�*6�4���>i��D�Gs�J�11111111111111111111111111111111111111111111111111111111111111111111111 @�y>@�h�@�q@�z@���@�]d@�i�@��@���@�N�@�3�@�v`@�2�@�YK@���@�hs@���@�
@}��@y�j@x�@t�|@p��@l��@i�@j�!@h��@go@c�4@`��@^�@\��@[��@Z��@Y�C@X�@W'�@W�Q@W�@V��@U��@Uo @T�O@TɆ@T_@R�@R�F@R��@Q�N@O{J@J��@C]�@:+k@2�'@YK?��6?�8�?h1'>��8���辂J�ǻ0��Y���6�*6�4���>i��D�Gs�J�11111111111111111111111111111111111111111111111111111111111111111111111 G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            PSAL                            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE from current cycle.                                                                                                                                                                                     PSAL_ADJSTED = celltm_sbe41(PSAL,TEMP,PRES,e_time,alpha,tau)                                                                                                                                                                                                                                                                                                                                                                                                                                                                    dP = -0.26                                                                                                                                                                                                                                                      e_time assumes 0.09dBar/s ascent rate, alpha=0.0141, tau=6.68s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Salinity corrected for Cell Thermal Mass (CTM), Johnson et al.(2007), JAOT & effects of pressure adjustments                                                                                                                                                                                                                                                                                                                                                                                                                    2015060504111820150605041118                BO  ARGQPREX2.0                                                                 20150605041118  CR                  G�O�G�O�G�O�                BO  ARGQRTSP1.0                                                                 20150605041118  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20150605041118  CV                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20150605041118  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20150605041118  CR                  G�O�G�O�G�O�                BO  ARGQPREX2.0                                                                 20150605041118  CR                  G�O�G�O�G�O�                BO  ARGQRTQC2.0                                                                 20150605041119  QCP$                G�O�G�O�G�O�6389758         BO  ARGQRTQC2.0                                                                 20150605041122  QCP$                G�O�G�O�G�O�6389758         BO  ARGQRTQC2.0                                                                 20150605041126  QCP$                G�O�G�O�G�O�6389758         