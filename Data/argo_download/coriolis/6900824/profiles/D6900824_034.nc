CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   F   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
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
_FillValue                    2�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    2�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    2�   PRES         
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  3�   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       3�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  5   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       5\   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       6t   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  7�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       7�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  8�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       94   PSAL         
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       :L   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  ;d   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ;�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  <�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
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
_FillValue                    G�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         G�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         G�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        G�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    H Argo profile     2.21.2 19500101000000  20101016034022  20120927141725  6900824 BSH                                                             Birgit KLEIN                                                    TEMP            PRES            PSAL               "A   IF  14632602                        2C  D   NEMO Profiling Float                                            860 @ծl��1   @ծl��@N��
=p��@#333331   ARGOS   A   A   A   @���A��Aq��A�33B��BD  Bj��B�ffB�ffB�ffB�ffB�33B�  C �fC33C�C  C)  C2��C<��CF��CK��C[33Ce�Cn�3Cx�fC���C�s3C�� C�ffC���C���C�� C��3C���C�� C�s3C�ffC�ffC�� C�Y�D	9�D��D"@ D.� D:��DG� DT33D`�fDmFfDy��D�#3D�` D���D��3D�fD�c3D��3D��3D�  D�c3D���D�fD�Y�Dڜ�D���D��D�P D�3D��f1111111111111111111111111111111111111111111111111111111111111111111111  @�  A��Ad��A���B��B@��Bg��B���B���B���B���Bؙ�B�ffC �C
ffCL�C33C(33C2  C<  CF  CJ��CZffCdL�Cm�fCx�C�&gC��C��C�  C�&gC�34C��C�L�C�&gC��C��C�  C�  C��C��4D	gDy�D"�D.��D:�gDG��DT  D`�3Dm3Dy�gD�	�D�FfD�� D�əD���D�I�D���D���D�fD�I�D��3D���D�@ Dڃ3D��3D�  D�6fD�D���1111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A$�A$�A$�A$��A�w@���@���@�E�@��\@� �@��-@�G�@��9@�j@��
@���@��@�D@�?}@��@��^@��@�E�@���@�`B@�I�@�n�@�C�@�7@�{@�=q@㝲@���@޸R@�(�@�+@�`B@��m@��y@Ͼw@�V@�x�@��m@���@�t�@��\@�K�@��y@��@���@��P@���@� �@|��@y&�@u�T@qx�@v�y@w�P@w�@w
=@r~�@q&�@p��@m�@mO�@h  @dj@a��@_
=1111111111111111111111111111111111111111111111111111111111111111111111  A$�A$�A$�A$��A�w@���@���@�E�@��\@� �@��-@�G�@��9@�j@��
@���@��@�D@�?}@��@��^@��@�E�@���@�`B@�I�@�n�@�C�@�7@�{@�=q@㝲@���@޸R@�(�@�+@�`B@��m@��y@Ͼw@�V@�x�@��m@���@�t�@��\@�K�@��y@��@���@��P@���@� �@|��@y&�@u�T@qx�@v�y@w�P@w�@w
=@r~�@q&�@p��@m�@mO�@h  @dj@a��@_
=1111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
��B
��B
��B
�BdZB0!B>wBJ�BbNBw�B�Br�Bs�Bv�Bv�By�B{�B�B�+B�=B�VB�hB�uB�hB�bB�DB�Bu�Bm�B`BBP�BF�B<jB49B)�B&�B �B�B�BJB1B  B�B��B�qB�B�!B��B��B�bB�B�B}�Bx�Bt�Br�Bm�B�B�7B�hB��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111  B3B3B3BB��BN�B]%BiqB�B��B��B�iB�pB��B��B��B��B��B��B��B�B�(B�6B�(B�"B�B��B�{B�GB~�Bo�BeSB[BR�BH�BE�B?fB9?B9?B*�B&�B�B.B�sB��BˉBΜB�aB� B��B��B�xB�^B�>B�#B�B��B�tB��B��B�B�$B�*B�CB�0B�bB�5B�B�B�1111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
TEMP            PRES            PSAL            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         PSAL_ADJUSTED=PSAL (re-calculated by using PRES_ADJUSTED) + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                           none                                                                                                                                                                                                                                                            Surface pressure = 0.8 dbar                                                                                                                                                                                                                                     OW : r=0.99991 , vertically averaged dS =0.029857                                                                                                                                                                                                               No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Pressure adjusted by using pressure offset at the sea surface. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                     Significant salinity drift present (salinity adjusted for pressure offset) - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                              201107211650232011072116502320110721165023  IF  ARGQCOAR1.0                                                                 20101016033138  QCF$                G�O�G�O�G�O�02000           IF  ARGQCOAR1.0                                                                 20101016033138  QCP$                G�O�G�O�G�O�DEB6C           IF  CORTCOOA5.2 RTQCGL01                                                        20101103075256  QCP$PSAL            G�O�G�O�G�O�                IF      SCOO1.3                                                                 20101102151121  QC                  G�O�G�O�G�O�                IF  CORTCOOA5.2 RTQCGL01                                                        20101103074325  QCP$TEMP            G�O�G�O�G�O�                GE  ARSQOW  1.0 ARGO climatology, Version: ARGO_for_DMQC_2011V02, 07/11         20110721165023  IP  PSAL            @���D��fG�O�                