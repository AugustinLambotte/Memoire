CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   N_CALIB       	N_HISTORY             	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       2016-12-28T16:37:42Z creation      
references        (http://www.argodatamgt.org/Documentation   comment           user_manual_version       3.03   Conventions       Argo-3.0 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      
_FillValue                    4�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    4�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    4�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    4�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    4�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    4�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    5   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  5   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  5L   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  5�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        5�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    5�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    5�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     5�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    5�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    5�   PLATFORM_TYPE                     	long_name         Type of float      
_FillValue                     5�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                    6   FIRMWARE_VERSION                  	long_name         Instrument version     
_FillValue                    6    WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    60   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       axis      T           64   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    6<   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            6@   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           6H   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           6P   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    6X   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    6\   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    6d   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    6h   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    6l   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    6p   CONFIG_MISSION_NUMBER                  	long_name         'Float's mission number for each profile    conventions       @0..N, 0 : launch mission (if exists), 1 : first complete mission   
_FillValue         ��        7p   PRES         
      
   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  7t   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  @@   PRES_ADJUSTED            
      
   	long_name         SEA PRESSURE   standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  Bt   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  K@   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  Mt   PSAL         
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  V@   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  _   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  a@   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  j   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  l@   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  u   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  }�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  ��   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    
_FillValue                  ,  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �4   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �8   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �<   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �@   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �D   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��Argo profile    3.0 1.2 19500101000000  20161228163742  20170111113733  6902617 BSH                                                             Birgit KLEIN                                                    PRES            PSAL            TEMP               -A   IF  51395660                        2C  D   NOVA                            SN187                           865 @��1   @��@R�!�D(N@�K�Z�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @   @9��@y��@�  @�  @�  A   A��A   A0  A>ffANffA^ffAp  A~ffA�  A���A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�33A�33A�33A�33B��B��B  B  B  B��B  B ffB$  B'��B,  B0ffB4  B8  B<ffB@ffBDffBHffBLffBP  BT  BXffB\  B`  Bd  Bh  Bl  Bp  Bt  BxffB|  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�33B�33B�33B�33B�  B�  C  C� C  C	� C  C� C  C� C  CffC�fC� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9�C;��C>  C@� CC  CEffCG�fCJ� CM  CO��CR�CT��CW�CY��C\�C^� Ca  Cc��Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw��Cz�C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C��3C��3C�33C�s3C�� C��C�L�C�� C��3C�  C�@ C�� C��3C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�33C�s3Cг3C�  C�L�CԀ C�� C��C�@ Cـ C�� C�  C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C��3C�@ C��C�� C�  C�@ C� C�� C��3C�33C�s3C�� C�  C�� C��3D � D  D@ D� D��D  D@ D	� D
�fDfD@ D� D� D  D@ D� D� D  D@ D�fD� D  D@ D� D� D fD!@ D"y�D#��D%  D&@ D'y�D(� D*fD+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?FfD@�fDA�fDC  DD@ DEy�DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DRfDSFfDT� DU� DW  DX@ DY� DZ��D[��D]@ D^�fD_� Da  Db@ Dc� Dd� Df  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� Dp  Dq9�Dr� Ds�fDu  Dv@ Dw� Dx��Dy��D{@ D|�fD}� D  D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D��D���D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�c3D�  D�� D�C3D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D��3D�` D�  DҠ D�@ D�� D�|�D�  D�� D�\�D�  Dנ D�@ D��3Dـ D�  D�� D�` D���Dܠ D�@ D�� Dރ3D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�3D� D�C3D�� D�|�D�  D�� D�` D�  D� D�@ D�� D� D�  D��3D�` D�  D� D�@ D�� D�|�D��D��D�` D�  D�� D�@ D��3D�� D��D�� D�` D�  D�� D�C3D��fD�vf11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @   @9��@y��@�  @�  @�  A   A��A   A0  A>ffANffA^ffAp  A~ffA�  A���A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�33A�33A�33A�33B��B��B  B  B  B��B  B ffB$  B'��B,  B0ffB4  B8  B<ffB@ffBDffBHffBLffBP  BT  BXffB\  B`  Bd  Bh  Bl  Bp  Bt  BxffB|  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�33B�33B�33B�33B�  B�  C  C� C  C	� C  C� C  C� C  CffC�fC� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9�C;��C>  C@� CC  CEffCG�fCJ� CM  CO��CR�CT��CW�CY��C\�C^� Ca  Cc��Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw��Cz�C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C��3C��3C�33C�s3C�� C��C�L�C�� C��3C�  C�@ C�� C��3C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�33C�s3Cг3C�  C�L�CԀ C�� C��C�@ Cـ C�� C�  C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C��3C�@ C��C�� C�  C�@ C� C�� C��3C�33C�s3C�� C�  C�� C��3D � D  D@ D� D��D  D@ D	� D
�fDfD@ D� D� D  D@ D� D� D  D@ D�fD� D  D@ D� D� D fD!@ D"y�D#��D%  D&@ D'y�D(� D*fD+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D8��D:@ D;� D<� D>  D?FfD@�fDA�fDC  DD@ DEy�DF��DH  DI@ DJ� DK� DM  DN@ DO� DP� DRfDSFfDT� DU� DW  DX@ DY� DZ��D[��D]@ D^�fD_� Da  Db@ Dc� Dd� Df  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� Dp  Dq9�Dr� Ds�fDu  Dv@ Dw� Dx��Dy��D{@ D|�fD}� D  D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D��D���D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�c3D�  D�� D�C3D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D��3D�` D�  DҠ D�@ D�� D�|�D�  D�� D�\�D�  Dנ D�@ D��3Dـ D�  D�� D�` D���Dܠ D�@ D�� Dރ3D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�3D� D�C3D�� D�|�D�  D�� D�` D�  D� D�@ D�� D� D�  D��3D�` D�  D� D�@ D�� D�|�D��D��D�` D�  D�� D�@ D��3D�� D��D�� D�` D�  D�� D�C3D��fD�vf11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��B`BB`BB`BB`BB`BB`BB`BB_;B`BB`BB`BB_;B`BB_;B`BB_;B_;B_;B_;B_;B_;B_;B_;B`BB`BB`BB_;B_;B`BB`BB`BB_;B`BB_;B_;B_;B_;B`BB_;B`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B_;B_;B`BB`BB`BB`BB`BB`BB`BB`BB`BBbNBdZBk�Bo�Bo�Bq�Bs�Bt�Bt�Bt�Bt�Bu�Bt�Bv�Bu�Bt�Bt�Bt�Bt�Bt�Bs�Bs�Bs�Br�Br�Bq�Bp�Bo�Bn�Bl�Bl�BjBiyBgmBe`BbNB`BB\)BZBXBXBXBW
BVBT�BR�BO�BL�BF�BE�BA�B>wB:^B7LB33B0!B-B+B(�B'�B&�B%�B#�B"�B�B�B�B"�B#�B&�B%�B%�B$�B#�B!�B�B�B�B�B�B{BuBhBbBVBPBDB	7B%BB  B��B��B��B��B�B�B�B�B�mB�TB�5B�#B��B��B��B��BɺBƨBÖB��B�}B�wB�qB�jB�^B�XB�LB�?B�9B�-B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B`BB`BB`BB`BB`BB`BB`BB_;B`BB`BB`BB_;B`BB_;B`BB_;B_;B_;B_;B_;B_;B_;B_;B`BB`BB`BB_;B_;B`BB`BB`BB_;B`BB_;B_;B_;B_;B`BB_;B`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B_;B_;B`BB`BB`BB`BB`BB`BB`BB`BB`BBbNBdZBk�Bo�Bo�Bq�Bs�Bt�Bt�Bt�Bt�Bu�Bt�Bv�Bu�Bt�Bt�Bt�Bt�Bt�Bs�Bs�Bs�Br�Br�Bq�Bp�Bo�Bn�Bl�Bl�BjBiyBgmBe`BbNB`BB\)BZBXBXBXBW
BVBT�BR�BO�BL�BF�BE�BA�B>wB:^B7LB33B0!B-B+B(�B'�B&�B%�B#�B"�B�B�B�B"�B#�B&�B%�B%�B$�B#�B!�B�B�B�B�B�B{BuBhBbBVBPBDB	7B%BB  B��B��B��B��B�B�B�B�B�mB�TB�5B�#B��B��B��B��BɺBƨBÖB��B�}B�wB�qB�jB�^B�XB�LB�?B�9B�-B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
@��@�  @�b@�  @�  @�  @�  @�b@� �@���@��@���@�1@�1@�  @��@���@��@���@��m@���@�1@�b@��@�  @���@�  @�  @�1@�b@�b@�b@���@��@��@���@�  @�1@�  @�1@�  @�b@�(�@�(�@�1'@� �@�1@�b@�1@�1@�b@�9X@�9X@�A�@�A�@�A�@�1'@�9X@�9X@�A�@�A�@�1'@� �@�(�@�9X@�A�@�A�@�A�@�A�@�A�@�A�@�A�@�A�@�Q�@�Q�@�Q�@�Q�@�Q�@�Q�@�Z@�Z@�Z@�Z@�bN@�bN@�Z@�I�@�Q�@�I�@�Z@�Q�@�Q�@�Z@�bN@�Z@�j@�Z@�Z@�bN@�bN@�Z@�j@�j@�j@�j@�j@�z�@��u@��u@���@���@���@���@��@���@���@���@���@��@���@��@��@��@��9@��j@�Ĝ@�Ĝ@��j@��j@�Ĝ@�Ĝ@���@���@���@���@��j@��j@�Ĝ@�Ĝ@���@���@���@���@���@���@���@���@���@���@���@���@��/@���@��/@��/@��/@���@�Ĝ@��j@�Z@�bN@�1@�\)@�33@�
=@���@�~�@�-@���@�9X@���@��H@�$�@��@��T@�`B@��@�z�@�Z@�(�@��P@�"�@�5?@�@���@�(�@�@�M�@��@�p�@��@��@��@�~�@~$�@y�@y�@y��@y7L@w;d@v��@u@t��@p�9@j��@f�@d�j@c��@_l�@]�@[��@W�w@N��@L1@K@I&�@HA�@G��@G+@CC�@?
=@<�@;��@<j@<(�@9�^@9hs@9&�@7��@6{@6$�@5?}@0��@)�@&ȴ@$1@#t�@"M�@��@��@�\@��@�`@�P@�/@��@ƨ@�@`B@�m?���?�ȴ?�z�?�5??�7L?��?�?���?�J?�v�?�?}?��T?�V?��H?���?y�?f�y?Y�#?R�?L�D?H��?@  ?5?,I�?��?�h?�/>�->ݲ->Õ�>�r�>�(�>�t�>w��>Kƨ>?|�>8Q�>(��>bN>J>%=��=�F=�x�=��=ȴ9=���=���=�t�=}�=aG�=��<���<t���o�D����1���ͼ�/�+�]/�m�h�ixս�C���1��E��ě����`����/��G���h�1'�hs��P��P����'.{�/��2-�333�:^5�=p��;dZ�<j�:^5�:^5�<j�=p��F��L�;Q녾]/�fff�j~��x���}󶾁%���7���\��J���\���˾�����9��ƨ����bN��hs��b������/���w���/���y���þ�V��&龴�j����KǾ��پ�Q쾹�#��dZ���۾�%��o����ě��\��J�\��$ݾ��;��t��Ձ�ؓu��/�������
���/���y��r�������1��1��h���&��������-��-��������-���j��ȴ���پ�Q��Q���پ�Q��Q��X��X���#���#��X��X������E���E���KǾ�X��^5��푾�p�����vɿ�\��
����ff�
���1�I����&�9X�Q�b����^5�p��p��p��/��w�!G��"��"��#o�"�\�"J�#o�"��"Mӿ"Mӿ"J�"Mӿ"��"�\�"Mӿ#o�#S��#���#S��"Mӿ!�7�"J�"�\�#S��#o�#���$Z�$Z�$���$�/�$�/�$���$���$���%`B�'l��'��'��(1'�(1'�(r��(�9�)�^�*���+��,I��,�D�,�Ϳ-O߿-��.V�.���.��/\)�/�;�0bN�0�׿1&�2�6E��6�+�6�+�6ȴ�7�P�7�ٿ8Q�8�u�9���:��:^5�:^5�:��:^5�:��9�#�9�#�:^5�;"ѿ;dZ�;��;��;��;��;��;��;��;��;�m�<(��<��=/�=/�=�-�=p��=�-�=�-�=�-�=�>5?�>5?�>5?�>5?�>5?�>5?�>5?�>5?�>vɿ>vɿ>�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@�  @�b@�  @�  @�  @�  @�b@� �@���@��@���@�1@�1@�  @��@���@��@���@��m@���@�1@�b@��@�  @���@�  @�  @�1@�b@�b@�b@���@��@��@���@�  @�1@�  @�1@�  @�b@�(�@�(�@�1'@� �@�1@�b@�1@�1@�b@�9X@�9X@�A�@�A�@�A�@�1'@�9X@�9X@�A�@�A�@�1'@� �@�(�@�9X@�A�@�A�@�A�@�A�@�A�@�A�@�A�@�A�@�Q�@�Q�@�Q�@�Q�@�Q�@�Q�@�Z@�Z@�Z@�Z@�bN@�bN@�Z@�I�@�Q�@�I�@�Z@�Q�@�Q�@�Z@�bN@�Z@�j@�Z@�Z@�bN@�bN@�Z@�j@�j@�j@�j@�j@�z�@��u@��u@���@���@���@���@��@���@���@���@���@��@���@��@��@��@��9@��j@�Ĝ@�Ĝ@��j@��j@�Ĝ@�Ĝ@���@���@���@���@��j@��j@�Ĝ@�Ĝ@���@���@���@���@���@���@���@���@���@���@���@���@��/@���@��/@��/@��/@���@�Ĝ@��j@�Z@�bN@�1@�\)@�33@�
=@���@�~�@�-@���@�9X@���@��H@�$�@��@��T@�`B@��@�z�@�Z@�(�@��P@�"�@�5?@�@���@�(�@�@�M�@��@�p�@��@��@��@�~�@~$�@y�@y�@y��@y7L@w;d@v��@u@t��@p�9@j��@f�@d�j@c��@_l�@]�@[��@W�w@N��@L1@K@I&�@HA�@G��@G+@CC�@?
=@<�@;��@<j@<(�@9�^@9hs@9&�@7��@6{@6$�@5?}@0��@)�@&ȴ@$1@#t�@"M�@��@��@�\@��@�`@�P@�/@��@ƨ@�@`B@�m?���?�ȴ?�z�?�5??�7L?��?�?���?�J?�v�?�?}?��T?�V?��H?���?y�?f�y?Y�#?R�?L�D?H��?@  ?5?,I�?��?�h?�/>�->ݲ->Õ�>�r�>�(�>�t�>w��>Kƨ>?|�>8Q�>(��>bN>J>%=��=�F=�x�=��=ȴ9=���=���=�t�=}�=aG�=��<���<t���o�D����1���ͼ�/�+�]/�m�h�ixս�C���1��E��ě����`����/��G���h�1'�hs��P��P����'.{�/��2-�333�:^5�=p��;dZ�<j�:^5�:^5�<j�=p��F��L�;Q녾]/�fff�j~��x���}󶾁%���7���\��J���\���˾�����9��ƨ����bN��hs��b������/���w���/���y���þ�V��&龴�j����KǾ��پ�Q쾹�#��dZ���۾�%��o����ě��\��J�\��$ݾ��;��t��Ձ�ؓu��/�������
���/���y��r�������1��1��h���&��������-��-��������-���j��ȴ���پ�Q��Q���پ�Q��Q��X��X���#���#��X��X������E���E���KǾ�X��^5��푾�p�����vɿ�\��
����ff�
���1�I����&�9X�Q�b����^5�p��p��p��/��w�!G��"��"��#o�"�\�"J�#o�"��"Mӿ"Mӿ"J�"Mӿ"��"�\�"Mӿ#o�#S��#���#S��"Mӿ!�7�"J�"�\�#S��#o�#���$Z�$Z�$���$�/�$�/�$���$���$���%`B�'l��'��'��(1'�(1'�(r��(�9�)�^�*���+��,I��,�D�,�Ϳ-O߿-��.V�.���.��/\)�/�;�0bN�0�׿1&�2�6E��6�+�6�+�6ȴ�7�P�7�ٿ8Q�8�u�9���:��:^5�:^5�:��:^5�:��9�#�9�#�:^5�;"ѿ;dZ�;��;��;��;��;��;��;��;��;�m�<(��<��=/�=/�=�-�=p��=�-�=�-�=�-�=�>5?�>5?�>5?�>5?�>5?�>5?�>5?�>5?�>vɿ>vɿ>�R11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            PSAL            TEMP            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          201701111137342017011111373420170111113734  IF  ARGQCOAR1.0                                                                 20161228150020  QCP$                G�O�G�O�G�O�09EBFC          IF  ARGQCOAR1.0                                                                 20161228150020  QCC$                G�O�G�O�G�O�000000          IF  ARGQCOAR1.0                                                                 20161228150020  QCF$                G�O�G�O�G�O�000000          GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2016V01 + ARGO climatology 20170111113735  IP  PSAL            @   D�vfG�O�                