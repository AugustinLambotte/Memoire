CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:18Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    6�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7T   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    7�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    7�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     7�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     88   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           8\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            8h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    8�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  Bd   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  M`   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  X\   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a$   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  cX   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  nT   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �8   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �<   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �@   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �D   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �H   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090450  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               $A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�S�g���1   @�S�(�x@O���*�@��O�2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @   @9��@�  @�  @�  @�  A   A  A!��A1��A@  AP  A^ffAp  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8ffB<  B@  BD  BH  BL  BO��BT  BX  B\  B`  Bc��Bh  Bl  Bp  BtffBx  B{��B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33C  C� C  C	ffC  C� C  C� C  C��C  C� C   C"ffC%  C'��C*  C,� C/  C1� C4  C6� C8�fC;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT��CW  CY� C\  C^� Ca  CcffCe�fCh� Ck  Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C�  C�@ C�� C�� C�  C�@ Cŀ C�� C��3C�@ Cʌ�C���C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�33Cـ C�� C��3C�@ Cހ C�� C�  C�@ C� C�3C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C�� C�  D �fDfDFfD�fD�fD  D@ D	� D
� D  D@ D� D� D  D9�D� D�fD  D@ D�fD� DfD@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2�fD4fD5@ D6y�D7� D9  D:@ D;� D<�fD>  D?@ D@� DA� DC  DD@ DE� DF��DH  DI@ DJ� DK� DM  DNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DY�fDZ� D\  D]FfD^� D_� Da  Db@ Dc�fDd� Df  Dg@ Dh� Di� Dk  Dl@ Dm�fDn� Dp  Dq@ Dr� Ds�fDufDv@ Dw� Dx� Dz  D{@ D|� D}� D  D��D���D�` D�  D��3D�@ D���D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�  D��3D�@ D���D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D���D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�c3D�3Dã3D�@ D�� Dŀ D�#3D�� D�` D�  DȠ D�@ D���D�|�D�  D�� D�` D���D͠ D�@ D�� Dπ D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  Dռ�D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D�3D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D�� D�` D�3D��3D�@ D�� D�� D��D�� D�` D�3D�� D�I�D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @   @9��@�  @�  @�  @�  A   A  A!��A1��A@  AP  A^ffAp  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8ffB<  B@  BD  BH  BL  BO��BT  BX  B\  B`  Bc��Bh  Bl  Bp  BtffBx  B{��B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33C  C� C  C	ffC  C� C  C� C  C��C  C� C   C"ffC%  C'��C*  C,� C/  C1� C4  C6� C8�fC;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT��CW  CY� C\  C^� Ca  CcffCe�fCh� Ck  Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C�  C�@ C�� C�� C�  C�@ Cŀ C�� C��3C�@ Cʌ�C���C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�33Cـ C�� C��3C�@ Cހ C�� C�  C�@ C� C�3C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C�� C�  D �fDfDFfD�fD�fD  D@ D	� D
� D  D@ D� D� D  D9�D� D�fD  D@ D�fD� DfD@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2�fD4fD5@ D6y�D7� D9  D:@ D;� D<�fD>  D?@ D@� DA� DC  DD@ DE� DF��DH  DI@ DJ� DK� DM  DNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DY�fDZ� D\  D]FfD^� D_� Da  Db@ Dc�fDd� Df  Dg@ Dh� Di� Dk  Dl@ Dm�fDn� Dp  Dq@ Dr� Ds�fDufDv@ Dw� Dx� Dz  D{@ D|� D}� D  D��D���D�` D�  D��3D�@ D���D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�  D��3D�@ D���D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D���D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�c3D�3Dã3D�@ D�� Dŀ D�#3D�� D�` D�  DȠ D�@ D���D�|�D�  D�� D�` D���D͠ D�@ D�� Dπ D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  Dռ�D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D�3D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D�� D�` D�3D��3D�@ D�� D�� D��D�� D�` D�3D�� D�I�D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�%@�%@�%@�V@�V@�V@��@��@�&�@��@�V@��@��@��@��@��@�V@���@��@���@���@���@�z�@��u@�Z@��@�"�@�@��@��@�=q@��@��T@��T@��T@��#@���@���@��@�x�@��@�x�@�x�@�p�@�hs@�hs@�?}@���@���@��j@��j@��@��D@�Q�@� �@��
@���@��P@��P@�l�@�\)@�\)@�\)@�\)@�K�@��@��H@��R@��!@�E�@�{@��-@�`B@�&�@���@���@��u@��F@�K�@�33@��H@���@�n�@�^5@�M�@�5?@��@�J@�@�$�@�$�@�J@��@��@���@���@��7@��@���@���@�o@���@��@��/@��`@��/@�j@��j@�V@���@���@��@��#@�p�@��@��D@�1@��
@���@��@���@��h@��@�?}@��@���@���@��u@�bN@��
@��F@���@���@��P@�dZ@�C�@�
=@�ȴ@�n�@�$�@��T@��h@��`@�bN@�9X@��@�  @��m@��w@���@�|�@�|�@�dZ@�\)@�C�@�o@���@��H@�ȴ@��!@�~�@�v�@�V@�-@�$�@��@��@��^@���@�x�@�x�@�`B@�?}@��@�V@�&�@��`@��/@�Ĝ@��@���@��u@��u@��D@�Z@�A�@�1'@� �@��@��
@��@�|�@�K�@�+@�"�@�o@�@��H@���@��+@�=q@���@���@�p�@�/@���@��j@�I�@��@���@���@�l�@�dZ@�\)@�33@���@���@�~�@�^5@�-@��@���@��h@�?}@�G�@�7L@��@�V@�V@�%@���@�%@�V@��@���@��/@��j@��j@�j@���@�l�@���@��@���@�ff@�n�@�E�@�@���@��7@�p�@�hs@�O�@�?}@�/@��@�V@�%@��@�(�@�K�@���@���@�=q@��^@�O�@���@�A�@���@��
@���@��P@��@��@�K�@�33@�"�@�o@�t�@��F@���@��w@��@�S�@���@���@���@���@���@��\@�~�@�ff@�=q@���@�`B@�V@��@�b@��;@��;@�ƨ@�|�@�t�@���@�
=@��H@�{@���@���@��@��@��@��@�%@���@�V@���@��@�r�@�(�@���@���@��@�;d@�@�ȴ@��!@�V@��+@�v�@�~�@�-@�@�x�@�G�@��@���@���@��@���@���@��@���@�z�@�j@�I�@�;@�w@�;@l�@+@~��@~��@~ff@~@}�-@}p�@}`B@|��@|�@{��@{S�@z�H@z��@z�\@z~�@zn�@z-@y�#@y��@yhs@x��@xQ�@xA�@xA�@x1'@xb@w�;@w�w@w�P@w�;@xA�@xr�@xĜ@x��@yhs@y�@y��@x��@xQ�@w�@v�y@vv�@vV@v{@u�h@u�h@u?}@t��@tj@sƨ@s��@s�
@t1@tj@s��@s��@r�@r��@r��@r�@rM�@r=q@rn�@r�!@r�H@s�F@u�@xbN@y��@x�9@x��@{o@|(�@|�j@}V@}V@|�@|�/@|�j@|��@|��@}/@}O�@|��@|�D@|Z@{33@y��@y��@z�@zn�@zM�@z-@y%@w�w@wK�@vȴ@w��@yhs@y�@x�`@xQ�@x  @w�@v��@vȴ@v�@vff@u�@u�h@uO�@u`B@u`B@u?}@u�@t��@t�/@t9X@s�m@s�F@sdZ@s@r��@r��@r��@rn�@r-@qx�@q&�@qG�@qX@qG�@q7L@p�`@p�u@pA�@pb@p  @o��@o��@o\)@o+@o;d@o;d@o\)@o|�@n��@n��@n�+@nV@nE�@n@mp�@m`B@m?}@l�/@lZ@l1@k�@kC�@j��@j^5@i��@i��@i�@i��@i�^@ihs@i%@hĜ@hbN@h1'@h  @g��@g|�@f��@f��@fff@e@e�h@e�@e`B@e?}@e�@eV@d��@dI�@c��@c�
@c�
@c�m@c��@cdZ@c"�@b�H@b�\@b~�@b^5@a�^@a��@a�7@aX@a7L@`��@`��@`bN@` �@_�P@^��@^ȴ@^�R1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�%@�%@�%@�V@�V@�V@��@��@�&�@��@�V@��@��@��@��@��@�V@���@��@���@���@���@�z�@��u@�Z@��@�"�@�@��@��@�=q@��@��T@��T@��T@��#@���@���@��@�x�@��@�x�@�x�@�p�@�hs@�hs@�?}@���@���@��j@��j@��@��D@�Q�@� �@��
@���@��P@��P@�l�@�\)@�\)@�\)@�\)@�K�@��@��H@��R@��!@�E�@�{@��-@�`B@�&�@���@���@��u@��F@�K�@�33@��H@���@�n�@�^5@�M�@�5?@��@�J@�@�$�@�$�@�J@��@��@���@���@��7@��@���@���@�o@���@��@��/@��`@��/@�j@��j@�V@���@���@��@��#@�p�@��@��D@�1@��
@���@��@���@��h@��@�?}@��@���@���@��u@�bN@��
@��F@���@���@��P@�dZ@�C�@�
=@�ȴ@�n�@�$�@��T@��h@��`@�bN@�9X@��@�  @��m@��w@���@�|�@�|�@�dZ@�\)@�C�@�o@���@��H@�ȴ@��!@�~�@�v�@�V@�-@�$�@��@��@��^@���@�x�@�x�@�`B@�?}@��@�V@�&�@��`@��/@�Ĝ@��@���@��u@��u@��D@�Z@�A�@�1'@� �@��@��
@��@�|�@�K�@�+@�"�@�o@�@��H@���@��+@�=q@���@���@�p�@�/@���@��j@�I�@��@���@���@�l�@�dZ@�\)@�33@���@���@�~�@�^5@�-@��@���@��h@�?}@�G�@�7L@��@�V@�V@�%@���@�%@�V@��@���@��/@��j@��j@�j@���@�l�@���@��@���@�ff@�n�@�E�@�@���@��7@�p�@�hs@�O�@�?}@�/@��@�V@�%@��@�(�@�K�@���@���@�=q@��^@�O�@���@�A�@���@��
@���@��P@��@��@�K�@�33@�"�@�o@�t�@��F@���@��w@��@�S�@���@���@���@���@���@��\@�~�@�ff@�=q@���@�`B@�V@��@�b@��;@��;@�ƨ@�|�@�t�@���@�
=@��H@�{@���@���@��@��@��@��@�%@���@�V@���@��@�r�@�(�@���@���@��@�;d@�@�ȴ@��!@�V@��+@�v�@�~�@�-@�@�x�@�G�@��@���@���@��@���@���@��@���@�z�@�j@�I�@�;@�w@�;@l�@+@~��@~��@~ff@~@}�-@}p�@}`B@|��@|�@{��@{S�@z�H@z��@z�\@z~�@zn�@z-@y�#@y��@yhs@x��@xQ�@xA�@xA�@x1'@xb@w�;@w�w@w�P@w�;@xA�@xr�@xĜ@x��@yhs@y�@y��@x��@xQ�@w�@v�y@vv�@vV@v{@u�h@u�h@u?}@t��@tj@sƨ@s��@s�
@t1@tj@s��@s��@r�@r��@r��@r�@rM�@r=q@rn�@r�!@r�H@s�F@u�@xbN@y��@x�9@x��@{o@|(�@|�j@}V@}V@|�@|�/@|�j@|��@|��@}/@}O�@|��@|�D@|Z@{33@y��@y��@z�@zn�@zM�@z-@y%@w�w@wK�@vȴ@w��@yhs@y�@x�`@xQ�@x  @w�@v��@vȴ@v�@vff@u�@u�h@uO�@u`B@u`B@u?}@u�@t��@t�/@t9X@s�m@s�F@sdZ@s@r��@r��@r��@rn�@r-@qx�@q&�@qG�@qX@qG�@q7L@p�`@p�u@pA�@pb@p  @o��@o��@o\)@o+@o;d@o;d@o\)@o|�@n��@n��@n�+@nV@nE�@n@mp�@m`B@m?}@l�/@lZ@l1@k�@kC�@j��@j^5@i��@i��@i�@i��@i�^@ihs@i%@hĜ@hbN@h1'@h  @g��@g|�@f��@f��@fff@e@e�h@e�@e`B@e?}@e�@eV@d��@dI�@c��@c�
@c�
@c�m@c��@cdZ@c"�@b�H@b�\@b~�@b^5@a�^@a��@a�7@aX@a7L@`��@`��@`bN@` �@_�P@^��@^ȴ@^�R1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B!�B!�B"�B"�B"�B"�B#�B$�B%�B&�B&�B&�B&�B&�B&�B%�B%�B$�B%�B%�B'�B'�B'�B'�B(�B)�B)�B)�B)�B)�B)�B)�B(�B'�B&�B&�B%�B$�B#�B"�B"�B"�B!�B!�B!�B �B �B �B!�B!�B!�B!�B!�B!�B!�B!�B"�B#�B#�B#�B$�B$�B#�B'�B(�B(�B+B,B-B.B+B&�B%�B$�B#�B"�B!�B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBoBoBhB\BVBVBPBPBPBJBJBJBDBDBDBDBDB
=B
=B
=B	7B	7B	7B1B1B1B1B+B+B+B%B%B%BBBB%BBBBBBBBBBBBBBBBBBB  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�mB�fB�ZB�ZB�ZB�TB�TB�NB�NB�HB�BB�BB�BB�;B�;B�;B�;B�5B�5B�5B�#B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBŢBÖBBBB��B��B�}B��B�}B�qB�jB�dB�dB�dB�dB�^B�^B�^B�^B�^B�XB�RB�FB�9B�9B�9B�3B�3B�-B�'B�'B�-B�-B�-B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�9B�LB�FB�LB�wB��BĜBŢBƨBƨBƨBƨBǮBǮBɺB��BɺBɺBȴBŢBÖBĜBŢBǮBǮBǮBŢBB��BBƨB��B��B��BɺBɺBɺBɺBȴBȴBǮBǮBƨBǮBǮBǮBȴBȴBȴBǮBǮBƨBƨBǮBǮBǮBǮBǮBǮBǮBǮBǮBȴBȴBȴBȴBȴBȴBɺBɺBɺBɺBɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBǮBǮBǮBǮBȴBȴBȴBǮBǮBǮBǮBǮBǮBǮBǮBƨBƨBŢBŢBĜBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBƨBƨBƨBƨBƨBƨBƨBƨBŢBŢBŢBŢBŢBŢBƨ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B B B B B!B"B"B"B"B"B"B"B!B B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B B!B!B#B$B%$B&*B#B B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B
�B
�B	�BuBoBoBiBiBiBcBcBcB]B]B]B]B]BVBVBVBPBPBPB JB JB JB JB�DB�DB�DB�>B�>B�>B�9B�9B�9B�>B�9B�9B�9B�9B�9B�9B�9B�9B�3B�3B�3B�,B�,B�,B�&B� B� B� B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B߉BނB�vB�vB�vB�pB�pB�jB�jB�dB�^B�^B�^B�XB�XB�XB�XB�RB�RB�RB�@B�-B�'B�!B�B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�~B�~B�~B�xB�rB�fB�YB�YB�YB�SB�SB�MB�GB�GB�MB�MB�MB�GB�AB�;B�5B�5B�/B�/B�)B�)B�)B�/B�/B�)B�)B�)B�#B�B�B�B�B�B�B�B�B�B�B�B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�
B�
B�B�B�B�B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�)B�ZB�lB�fB�lB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.007704                                                                                                                                                                                                              No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904502018110709045020181107090450  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172718  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172718  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090450  IP  PSAL            @   D���G�O�                