CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:20Z creation; 2018-10-26T19:00:44Z last update (coriolis COQC software)   
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
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181026190020  20181119104027  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               1A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @��3��q1   @�재9@S�v�:@�/�Մ�2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@9��@y��@�  @�  @�  A   A  A!��A0  A@  AP  A`  Ap  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  B   B  B  B  B��B  B��B��B��B$  B(  B,  B0ffB4ffB8ffB<ffB@ffBDffBHffBL  BO��BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|ffB�  B�  B�  B�  B�  B�33B�  B�  B���B���B���B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�33B�33B�33B�33B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  C  C� C  C	� C�C� C  C� C  C� C  C� C �C"��C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>�C@� CC  CE� CH  CJ� CM  CO��CR�CT� CV�fCYffC[�fC^ffCa  Cc� Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C���C���C��C�@ C�� C�� C��C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�33C�� C�� C��3C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C�  C�@ Cʀ C˳3C�  C�@ C�s3C�� C�  C�@ CԀ Cճ3C�  C�33Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C��3C�33C�s3C�3C��3C�@ C��C�� C�  C�@ C� C�� C��C�@ C�s3C�� C��C�s3C��3D ��D��D@ D� D� D  D@ D	� D
� D  D@ D� D� DfD@ Dy�D� D  D9�D� D� D��D@ D� D�fD fD!@ D"� D#� D%  D&@ D'� D(� D*  D+9�D,� D-� D/  D0@ D1� D2�fD4  D59�D6y�D7� D9  D:@ D;� D<��D>  D?@ D@� DA� DCfDD@ DE� DF� DHfDIFfDJ� DK� DM  DN@ DO�fDP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� D`��Db@ Dc� Dd� DffDgFfDh� Di� Dk  Dl@ Dm� Dn�fDp  Dq@ Dr� Ds� DufDvFfDw� Dx��Dz  D{FfD|� D}� D  D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�#3D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D��D���D�\�D�  D��3D�@ D�� D��3D�  D�� D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���Dà D�@ D���Dŀ D�  D�� D�` D���DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  Dм�D�\�D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�#3D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�  D��3D�C3D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@9��@y��@�  @�  @�  A   A  A!��A0  A@  AP  A`  Ap  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  B   B  B  B  B��B  B��B��B��B$  B(  B,  B0ffB4ffB8ffB<ffB@ffBDffBHffBL  BO��BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|ffB�  B�  B�  B�  B�  B�33B�  B�  B���B���B���B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�33B�33B�33B�33B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  C  C� C  C	� C�C� C  C� C  C� C  C� C �C"��C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>�C@� CC  CE� CH  CJ� CM  CO��CR�CT� CV�fCYffC[�fC^ffCa  Cc� Cf�Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C���C���C��C�@ C�� C�� C��C�@ C�� C���C�  C�@ C�s3C�� C�  C�@ C�� C�� C�  C�33C�� C�� C��3C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C���C�  C�@ Cʀ C˳3C�  C�@ C�s3C�� C�  C�@ CԀ Cճ3C�  C�33Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C��3C�33C�s3C�3C��3C�@ C��C�� C�  C�@ C� C�� C��C�@ C�s3C�� C��C�s3C��3D ��D��D@ D� D� D  D@ D	� D
� D  D@ D� D� DfD@ Dy�D� D  D9�D� D� D��D@ D� D�fD fD!@ D"� D#� D%  D&@ D'� D(� D*  D+9�D,� D-� D/  D0@ D1� D2�fD4  D59�D6y�D7� D9  D:@ D;� D<��D>  D?@ D@� DA� DCfDD@ DE� DF� DHfDIFfDJ� DK� DM  DN@ DO�fDP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� D`��Db@ Dc� Dd� DffDgFfDh� Di� Dk  Dl@ Dm� Dn�fDp  Dq@ Dr� Ds� DufDvFfDw� Dx��Dz  D{FfD|� D}� D  D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�#3D��3D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D��D���D�\�D�  D��3D�@ D�� D��3D�  D�� D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���Dà D�@ D���Dŀ D�  D�� D�` D���DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  Dм�D�\�D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�#3D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�  D��3D�C3D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@y��@y�7@y��@zn�@z~�@z~�@z~�@z~�@zn�@z~�@z~�@zM�@zM�@zM�@zn�@z�@z~�@z�\@{o@{"�@z�@{@{@{@{@{@{@{o@{o@{"�@{33@{33@{33@{33@{33@{C�@{C�@{C�@{33@{33@{33@{33@{33@{33@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{S�@{S�@{S�@{S�@{S�@{S�@{S�@{dZ@{dZ@{t�@{��@{��@{��@{��@{��@{��@{�F@{ƨ@{��@{�@{t�@{��@{��@{t�@{t�@{t�@{t�@{�@{t�@{�@{�@{�@{��@{��@{��@{�F@{�F@{��@{�F@{��@{��@{��@{��@{��@{��@{��@{�F@{ƨ@{ƨ@{ƨ@{�m@|1@|1@|�@|�@{�m@{�m@|1@{��@|�@|(�@{�m@{�
@{�
@{ƨ@{�m@{�
@|1@|1@|I�@|I�@|9X@|9X@|9X@|I�@|(�@{��@{ƨ@|�@|1@{�m@|1@{�m@|9X@{��@yx�@y7L@yx�@y�#@y�@y�@y�#@y�^@y�#@y�^@z=q@z~�@z~�@y��@z^5@z~�@s�@hr�@T�D@Gl�@B=q@D��@Ct�@F��@I�@F��@Ct�@?��@>�y@=p�@=?}@<z�@6�@5p�@4��@0�9@/;d@-�@(Ĝ@'l�@$�/@!��@��@V@1@��@
=@��@�h@�?��H?�M�?ѩ�?ɺ^?�%?��?���?�G�?�7L?�K�?�E�?��?��?���?���?��/?�z�?��
?�S�?���?�G�?�{?�I�?��^?��#?���?�?|(�?p �?mV?k?d�/?\j?S�F?J=q??�w?6?3��?333?/\)?*��?'�?|�?K�?K�?�`?\)?��?r�?�? A�>��m>�1>ؓu>׍P>�"�>�V>��j>���>�{>�V>���>�r�>��y>�S�>�G�>��>��>�V>��>�$�>���>cS�>P�`>%�T>+=�G�=��=��=�j=��=�O�=<j=�P=o<���<e`B;o��o���
��`B�o�o�+�#�
�@��ixս�\)���
��{��E���vɽě���
=��;d��F���C��n���R�-V�/��1&�333�5?}�B�\�I�^�R�Xb�["Ѿ^5?�_;d�aG��cS��dZ�dZ�fff�hr��l�D�m�h�r�!�u�y�#�|푾}�}�~�۾~�۾����%��J�������˾�������������������+������I���O߾�\)��bN���`��hs���`���`��hs��n������b���u��������������"Ѿ��㾜(����R��;d���w��G����/���y��r����D��9X��Q쾻dZ��푾�vɾ���Ǯ��C���O߾�bN���ؓu����ڟ��ݲ-��5?��5?�޸R�߾w��A�������S����/��l���~���V���׾�!���j��E���Q��j��p���vɾ�|� A�� A�� Ĝ��7�J��\������l���ÿ	�^�
=q�
~��ƨ��D�O߿�h�����{����\)��������;�bN�&�n���33��F��Ͽz��j������+�Q�������H�"ѿdZ�"ѿ"ѿ"ѿ"ѿdZ��m�푿�-�5?�vɿvɿ;d� A�� Ĝ�!�7�!���"Mӿ#S��#�
�$��$Z�$���$�/�$�/�%��%`B�%`B�%�˿&��'+�'��'(1'�(�9�(�ÿ)7L�)xտ)��)�^�*���+�+C��,I��,�D�-V�-��-��.{�.V�.���/���0 ſ0�׿0�`�0�`�0�`�1hs�1���1녿2n��2n��2�333�3�Ͽ4z�4���5?}�5�5�6�6�6�+�6ȴ�7Kǿ7�P�8�u�8���8���8���9��9���:��:��:��9�#�9�#�9�#�9�#�9�#�9�#�9���9���9���9�#�:��9�#�:��:��:��:^5�:^5�:���:���:^5�:^5�:�H�:�H�:���:�H�:���:���:���:���:���:���:���:���:��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @y��@y�7@y��@zn�@z~�@z~�@z~�@z~�@zn�@z~�@z~�@zM�@zM�@zM�@zn�@z�@z~�@z�\@{o@{"�@z�@{@{@{@{@{@{@{o@{o@{"�@{33@{33@{33@{33@{33@{C�@{C�@{C�@{33@{33@{33@{33@{33@{33@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{C�@{S�@{S�@{S�@{S�@{S�@{S�@{S�@{dZ@{dZ@{t�@{��@{��@{��@{��@{��@{��@{�F@{ƨ@{��@{�@{t�@{��@{��@{t�@{t�@{t�@{t�@{�@{t�@{�@{�@{�@{��@{��@{��@{�F@{�F@{��@{�F@{��@{��@{��@{��@{��@{��@{��@{�F@{ƨ@{ƨ@{ƨ@{�m@|1@|1@|�@|�@{�m@{�m@|1@{��@|�@|(�@{�m@{�
@{�
@{ƨ@{�m@{�
@|1@|1@|I�@|I�@|9X@|9X@|9X@|I�@|(�@{��@{ƨ@|�@|1@{�m@|1@{�m@|9X@{��@yx�@y7L@yx�@y�#@y�@y�@y�#@y�^@y�#@y�^@z=q@z~�@z~�@y��@z^5@z~�@s�@hr�@T�D@Gl�@B=q@D��@Ct�@F��@I�@F��@Ct�@?��@>�y@=p�@=?}@<z�@6�@5p�@4��@0�9@/;d@-�@(Ĝ@'l�@$�/@!��@��@V@1@��@
=@��@�h@�?��H?�M�?ѩ�?ɺ^?�%?��?���?�G�?�7L?�K�?�E�?��?��?���?���?��/?�z�?��
?�S�?���?�G�?�{?�I�?��^?��#?���?�?|(�?p �?mV?k?d�/?\j?S�F?J=q??�w?6?3��?333?/\)?*��?'�?|�?K�?K�?�`?\)?��?r�?�? A�>��m>�1>ؓu>׍P>�"�>�V>��j>���>�{>�V>���>�r�>��y>�S�>�G�>��>��>�V>��>�$�>���>cS�>P�`>%�T>+=�G�=��=��=�j=��=�O�=<j=�P=o<���<e`B;o��o���
��`B�o�o�+�#�
�@��ixս�\)���
��{��E���vɽě���
=��;d��F���C��n���R�-V�/��1&�333�5?}�B�\�I�^�R�Xb�["Ѿ^5?�_;d�aG��cS��dZ�dZ�fff�hr��l�D�m�h�r�!�u�y�#�|푾}�}�~�۾~�۾����%��J�������˾�������������������+������I���O߾�\)��bN���`��hs���`���`��hs��n������b���u��������������"Ѿ��㾜(����R��;d���w��G����/���y��r����D��9X��Q쾻dZ��푾�vɾ���Ǯ��C���O߾�bN���ؓu����ڟ��ݲ-��5?��5?�޸R�߾w��A�������S����/��l���~���V���׾�!���j��E���Q��j��p���vɾ�|� A�� A�� Ĝ��7�J��\������l���ÿ	�^�
=q�
~��ƨ��D�O߿�h�����{����\)��������;�bN�&�n���33��F��Ͽz��j������+�Q�������H�"ѿdZ�"ѿ"ѿ"ѿ"ѿdZ��m�푿�-�5?�vɿvɿ;d� A�� Ĝ�!�7�!���"Mӿ#S��#�
�$��$Z�$���$�/�$�/�%��%`B�%`B�%�˿&��'+�'��'(1'�(�9�(�ÿ)7L�)xտ)��)�^�*���+�+C��,I��,�D�-V�-��-��.{�.V�.���/���0 ſ0�׿0�`�0�`�0�`�1hs�1���1녿2n��2n��2�333�3�Ͽ4z�4���5?}�5�5�6�6�6�+�6ȴ�7Kǿ7�P�8�u�8���8���8���9��9���:��:��:��9�#�9�#�9�#�9�#�9�#�9�#�9���9���9���9�#�:��9�#�:��:��:��:^5�:^5�:���:���:^5�:^5�:�H�:�H�:���:�H�:���:���:���:���:���:���:���:���:��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBffBgmBffBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBffBgmBgmBgmBgmBgmBgmBgmBgmBffBffBgmBffBgmBgmBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBgmBgmBffBffBffBffBffBffBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBe`BffBffBffBffBe`Be`Be`BbNBbNBbNBcTBcTBcTBbNBbNBbNBbNBbNBcTBbNBaHB_;BZBO�B<jB(�B�B�B�B�B�B �B�B�B�B�B�B�B{BoBhB\BVBJB
=B1B%BBBB��B��B��B��B�B�B�sB�ZB�;B�)B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBǮBŢBĜBÖBB��B��B�}B�wB�jB�jB�dB�dB�^B�^B�XB�RB�RB�LB�LB�LB�FB�FB�?B�9B�9B�3B�3B�-B�-B�'B�!B�'B�'B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  BffBgmBffBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBffBgmBgmBgmBgmBgmBgmBgmBgmBffBffBgmBffBgmBgmBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBgmBgmBffBffBffBffBffBffBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBgmBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBffBe`BffBffBffBffBe`Be`Be`BbNBbNBbNBcTBcTBcTBbNBbNBbNBbNBbNBcTBbNBaHB_;BZBO�B<jB(�B�B�B�B�B�B �B�B�B�B�B�B�B{BoBhB\BVBJB
=B1B%BBBB��B��B��B��B�B�B�sB�ZB�;B�)B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBǮBŢBĜBÖBB��B��B�}B�wB�jB�jB�dB�dB�^B�^B�XB�RB�RB�LB�LB�LB�FB�FB�?B�9B�9B�3B�3B�-B�-B�'B�!B�'B�'B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040272018111910402720181119104027  IF  ARFMCODA024c                                                                20181026190020                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190044  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181026190044  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104027  IP  PSAL            @ffD��3G�O�                