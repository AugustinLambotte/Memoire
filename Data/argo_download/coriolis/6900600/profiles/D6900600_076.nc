CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   F   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       2012-08-29T12:48:03Z creation      
references        (http://www.argodatamgt.org/Documentation   comment           user_manual_version       2.4    Conventions       Argo-2.4 CF-1.6    featureType       trajectoryProfile         >   	DATA_TYPE                  	long_name         	Data type      
_FillValue                    2�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    2�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    2�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    2�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    3   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    3   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    3,   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  34   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  3t   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  3�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        3�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    3�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    3�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     3�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    4   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    4   INST_REFERENCE                    	long_name         Instrument type    conventions       Brand, type, serial number     
_FillValue                  @  4   FIRMWARE_VERSION                  	long_name         Instrument version     conventions           
_FillValue                    4X   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    4h   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T           4l   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    4t   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            4x   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           4�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           4�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    4�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    4�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    4�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    4�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    4�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    4�   PRES         
      	   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       5�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  6�   PRES_ADJUSTED            
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       7   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  8    PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       8h   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       9�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  :�   TEMP_ADJUSTED            
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       :�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  ;�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       <@   PSAL         
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_practical_salinity   
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       =X   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  >p   PSAL_ADJUSTED            
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       >�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  H  ?�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       @   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  A0   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    A`   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    D`   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    G`   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    
_FillValue                  ,  J`   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    J�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    J�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    J�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    J�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  J�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    J�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    J�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    J�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         K    HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         K   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        K   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    KArgo profile    2.3 1.2 19500101000000  20120827044334  20130709105252  6900600 BSH                                                             Birgit KLEIN                                                    PRES            PSAL            TEMP               LA   IF  26950128                        2C  D   NEMO Profiling Float                                                            860 @�Xp/4�1   @�Xp/4�@P�5?|��*����o1   ARGOS   A   A   A   Primary sampling: discrete                                                                                                                                                                                                                                      @�ff@���A�33A�  B  BD��Bn  B�ffB���B�  B���B���B�ffC��C��C� CffC&��C3�3C;�CGffCQffC[� Ce� CjL�Cy�3C���C���C��fC���C�� C��3C��fC�L�C�ٚC��3C��3C���C�� C���C�� C̳3C���C�� CۦfC�fC䙚C���C�ٚC���C��fD  D	S3D� D�3D  D"Y�D.�3D;Y�DG��DS��D`ٚDmY�Dy��D���D�p D���D��D�)�D�i�1111111111111111111111111111111111111111111111111111111111111111111111  @�ff@���A�33A�  B  BD��Bn  B�ffB���B�  B���B���B�ffC��C��C� CffC&��C3�3C;�CGffCQffC[� Ce� CjL�Cy�3C���C���C��fC���C�� C��3C��fC�L�C�ٚC��3C��3C���C�� C���C�� C̳3C���C�� CۦfC�fC䙚C���C�ٚC���C��fD  D	S3D� D�3D  D"Y�D.�3D;Y�DG��DS��D`ٚDmY�Dy��D���D�p D���D��D�)�D�i�1111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��AbNAn�@�E�@��!@���@?W��>��>O�;=@�<u�,1�������.{�t���/���
=o=��
=ě�=P�`=���=�%='�=���=�h=���<�/<o�49X�o�#�
�0 ŽT����7L���E����ͽ�F�%�C��z��-�'6E��=p��7KǾB�\�H�9�["Ѿw�پ�1'���;������b��ff��-��1'��MӾ��H�S�� ſt���+��m�|�'l��,�Ϳ1�1111111111111111111111111111111111111111111111111111111111111111111111  AbNAn�@�E�@��!@���@?W��>��>O�;=@�<u�,1�������.{�t���/���
=o=��
=ě�=P�`=���=�%='�=���=�h=���<�/<o�49X�o�#�
�0 ŽT����7L���E����ͽ�F�%�C��z��-�'6E��=p��7KǾB�\�H�9�["Ѿw�پ�1'���;������b��ff��-��1'��MӾ��H�S�� ſt���+��m�|�'l��,�Ϳ1�1111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
�B
�B
��B
�B
�BbB0!B8RB>wBB�BF�BJ�BJ�BO�BO�BT�B^5BbNBp�Bv�B|�B~�B�B�B�B�\B�uB�{B�bB�bB�hB�bB�bB�hB�oB�uB�oB�uB�{B�{B�uB�{B�uB�{B�uB�{B�{B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111  B
�B
�B
��B
�B
�BcB0"B8SB>xBB�BF�BJ�BJ�BO�BO�BT�B^6BbOBp�Bv�B|�B~�B�B�B�B�\B�uB�{B�bB�bB�hB�bB�bB�hB�oB�uB�oB�uB�{B�{B�uB�{B�uB�{B�uB�{B�{B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            PSAL            TEMP            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            Surface pressure = 0 dbar                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressure adjusted by using pressure offset at the sea surface. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                     No significant salinity drift detected (salinity adjusted for pressure offset). OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                      No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          201307091052522013070910525220130709105252  IF  ARGQCOAR1.0                                                                 20120827033905  QCF$                G�O�G�O�G�O�00000           IF  ARGQCOAR1.0                                                                 20120827033905  QCP$                G�O�G�O�G�O�DEB6C           IF  ARGQSCOO1.4                                                                 20120829135531  CF  POS$            G�O�G�O�@@                  IF      SCOO1.4                                                                 20120829135532  QC                  G�O�G�O�G�O�                GE  ARSQOW  1.0 ARGO CTD reference database, Version: CTD_for_DMQC_2013V01, 3/1320130709105252  IP  PSAL            @�ffD�i�G�O�                