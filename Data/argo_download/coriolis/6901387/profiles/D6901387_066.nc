CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   U   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
_FillValue                    0�   FORMAT_VERSION                 comment       File format version    
_FillValue                    0�   HANDBOOK_VERSION               comment       Data handbook version      
_FillValue                    0�   REFERENCE_DATE_TIME                 comment       !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    1    DATE_CREATION                   comment       Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    1   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    1    PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    10   PROJECT_NAME                  comment       Name of the project    
_FillValue                  @  18   PI_NAME                   comment       "Name of the principal investigator     
_FillValue                  @  1x   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  1�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        1�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    1�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    1�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     1�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    2   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    2   INST_REFERENCE                    	long_name         Instrument type    conventions       Brand, type, serial number     
_FillValue                  @  2   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    2\   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            2`   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    2h   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            2l   LATITUDE               	long_name         &Latitude of the station, best estimate     units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�             2t   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�             2|   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    2�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    2�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    2�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    2�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    2�   PRES         
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     T  2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  3�   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     T  4H   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  5�   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     T  5�   PSAL         
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  7H   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  8�   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  8�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  :H   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  :�   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  ;�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  =H   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  =�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  X  >�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     T  ?L   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  @�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    @�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    C�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    F�   CALIBRATION_DATE            	             
_FillValue                  ,  I�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    I�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    J    HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    J   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    J   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  J   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    JL   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    J\   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    J`   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Jp   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Jt   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Jx   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    J|Argo profile    2.2 1.2 19500101000000  20120316154510  20121128094853  6901387 ARGO POLAND                                                     Waldemar WALCZOWSKI                                             PRES            PSAL            TEMP               BA   IF  25685536                        2C  D   NEMO Profiling Float                                            860 @�����G1   @�����G@S�\��N<@T�w[�01   GPS     B   B   B   =���        >L��>L��>���>L��@�ffA  A^ffA���A噚B��BC��Bk33B���B�33B���B�  B�  B���C �fC� C33C��C)  C3ffC=�CG  CQ�C[L�CeL�Co  Cx�C��3C�� C���C�L�C���C�L�C���C�s3C�L�C���C���C���C�� C��C�� C�� C��3D��D	9�D��D�fD3D"` D(y�D.�3D;  DG��DT@ D`�3DmFfDy�fD�fD�Y�D���D�ٚD�fD�c3D��3D��D�#3D�l�D��3D��3D�,�D�` Dڣ3D�� D�&fD�l�D�D�� 1441414111111111111111111111111111111111111111111111111111111111111111111111111111111   =���G�O�G�O�>L��G�O�>���G�O�@�ffA  A^ffA���A噚B��BC��Bk33B���B�33B���B�  B�  B���C �fC� C33C��C)  C3ffC=�CG  CQ�C[L�CeL�Co  Cx�C��3C�� C���C�L�C���C�L�C���C�s3C�L�C���C���C���C�� C��C�� C�� C��3D��D	9�D��D�fD3D"` D(y�D.�3D;  DG��DT@ D`�3DmFfDy�fD�fD�Y�D���D�ٚD�fD�c3D��3D��D�#3D�l�D��3D��3D�,�D�` Dڣ3D�� D�&fD�l�D�D�� 1441414111111111111111111111111111111111111111111111111111111111111111111111111111111   @��G�O�G�O�@��G�O�@��G�O�@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B��B
�PB�B@ޟ�A-@^@��BBÖB��B��B��BĜBǮB�BBoB{BhBJB1BDBDBJBuB\BPB1BB	7BJBJBDBDBJBuBoB�B�BuB�B�B �B$�B)�B$�B1'B)�B{BB
=B�wB��B�B�?BǮB��BÖB�B��B�9B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1144444411111111111111111111111111111111111111111111111111111111111111111111111111111   G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�BÖB��B��B��BĜBǮB�BBoB{BhBJB1BDBDBJBuB\BPB1BB	7BJBJBDBDBJBuBoB�B�BuB�B�B �B$�B)�B$�B1'B)�B{BB
=B�wB��B�B�?BǮB��BÖB�B��B�9B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��4444444411111111111111111111111111111111111111111111111111111111111111111111111111111   G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
?�ȴ?o\)?W��?���?� �?�1?��?��^?�?��@ ��?�dZ?���?���@	�^@)��@)��@(�@#S�@!��@!G�@ A�@�@"-@ 1'@V@dZ@�m@v�@
=@ȴ@�@@��@�R@�@�@�u@�@��@�P@"�@�@�@�#@E�@%@t�?��j?�l�?�"�?Gl�?bJ?^��?~��?vE�?Tz�?1>��>�b=�vɽ�hs�8Q쾏�;�����Ǯ�ۥ��h���ٿC����n��
=�^5�#o�!���,�Ϳ3t��4z�=/�=�-�?�w�Co�D��1144441111111111111111111111111111111111111111111111111111111111111111111111111111111   ?�ȴG�O�G�O�G�O�G�O�G�O�G�O�?��^?�?��@ ��?�dZ?���?���@	�^@)��@)��@(�@#S�@!��@!G�@ A�@�@"-@ 1'@V@dZ@�m@v�@
=@ȴ@�@@��@�R@�@�@�u@�@��@�P@"�@�@�@�#@E�@%@t�?��j?�l�?�"�?Gl�?bJ?^��?~��?vE�?Tz�?1>��>�b=�vɽ�hs�8Q쾏�;�����Ǯ�ۥ��h���ٿC����n��
=�^5�#o�!���,�Ϳ3t��4z�=/�=�-�?�w�Co�D��1444444111111111111111111111111111111111111111111111111111111111111111111111111111111   ;oG�O�G�O�G�O�G�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            PSAL            TEMP            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          201211280948532012112809485320121128094853  IF  ARGQCOAR1.0                                                                 20120316145416  QCP$                G�O�G�O�G�O�DEB6C           IF  ARGQCOAR1.0                                                                 20120316145416  QCF$                G�O�G�O�G�O�04B20           IF      SCOO1.4                                                                 20120316161009  QC                  G�O�G�O�G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2012V01 + ARGO climatology 20121128094853  IP  PSAL            =���D�� G�O�                