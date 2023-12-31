CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  .   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:10Z last update (coriolis COQC software)   
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
_FillValue                 0  BT   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  M<   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  Ol   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  X$   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  `�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  c   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  k�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  m�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  v�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  d   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 0  �L   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �|   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �4   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �d   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �d   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �d   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �d             ,  �dArgo profile    3.1 1.2 19500101000000  20181105172617  20181107090442  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�0�bM��1   @�0��8T�@Oo7���A\�W	�2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A�  A�33A�33A�  A�  A�  A�  A���A�  A�  A�  A�  B   B  B  BffBffBffB  B  B   B$ffB(ffB,ffB0ffB4ffB8  B<  B@ffBD  BHffBLffBP  BT  BX  B\ffB`  Bd  Bh  Bl  Bp  Bs��Bx  B|  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�33B�33B�33B�  B�  B�33B�  B�  B�  B�  B�  B�33B�33B�33B�  B�  B�  B���B���B�  B�  B�  B�33B�33B�  B�  B�  B�33B�33B�  B�  C  C� C  C	� C  C� C  C� C  C� C�fC� C   C"ffC%  C'� C*  C,� C.�fC1� C4�C6� C9  C;ffC>  C@� CC  CE��CH�CJ� CM  CO� CR�CT��CW  CY� C\  C^� Ca  Cc��Cf�Ch��Ck  CmffCp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�33C�� C�� C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�33C�� C�� C�  C�L�C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ Cŀ C�� C��C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C���C��C�@ C� C�� C�  C�@ C� C�� C��C�@ C�s3C�� C��C�@ C�� C�� C�  C�� C�  D � D  DFfD� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D�fD  D9�D� D� D  D@ D� D� D fD!@ D"� D#�fD%  D&@ D'� D(� D*  D+@ D,y�D-� D/fD0FfD1�fD2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>fD?@ D@y�DA� DC  DD@ DE� DF�fDH  DI@ DJ� DK� DM  DN9�DO� DP� DR  DS@ DT� DU� DW  DX9�DYy�DZ� D\  D]9�D^� D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�  D�� D�c3D�  D���D�<�D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D��D���D�` D�3D��3D�C3D�� D�|�D�  D�� D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D��3D�  D�� D�c3D�3D��3D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�#3D��3D�` D���D�� D�@ D�� D��3D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� DŃ3D�#3D�� D�c3D�3Dȣ3D�C3D�� Dʀ D�  D��3D�` D�  D͠ D�@ D�� Dσ3D�  D��3D�` D�  DҠ D�C3D��3DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�c3D�  Dܠ D�C3D�� Dހ D�#3D��3D�` D���D� D�@ D���D� D�  D�� D�` D�3D�3D�@ D��3D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D��D�` D�3D��3D�<�D��3D��3D�#3D��fD�S3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A�  A�33A�33A�  A�  A�  A�  A���A�  A�  A�  A�  B   B  B  BffBffBffB  B  B   B$ffB(ffB,ffB0ffB4ffB8  B<  B@ffBD  BHffBLffBP  BT  BX  B\ffB`  Bd  Bh  Bl  Bp  Bs��Bx  B|  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�33B�33B�33B�  B�  B�33B�  B�  B�  B�  B�  B�33B�33B�33B�  B�  B�  B���B���B�  B�  B�  B�33B�33B�  B�  B�  B�33B�33B�  B�  C  C� C  C	� C  C� C  C� C  C� C�fC� C   C"ffC%  C'� C*  C,� C.�fC1� C4�C6� C9  C;ffC>  C@� CC  CE��CH�CJ� CM  CO� CR�CT��CW  CY� C\  C^� Ca  Cc��Cf�Ch��Ck  CmffCp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�33C�� C�� C��3C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�33C�� C�� C�  C�L�C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ Cŀ C�� C��C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C���C��C�@ C� C�� C�  C�@ C� C�� C��C�@ C�s3C�� C��C�@ C�� C�� C�  C�� C�  D � D  DFfD� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D�fD  D9�D� D� D  D@ D� D� D fD!@ D"� D#�fD%  D&@ D'� D(� D*  D+@ D,y�D-� D/fD0FfD1�fD2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>fD?@ D@y�DA� DC  DD@ DE� DF�fDH  DI@ DJ� DK� DM  DN9�DO� DP� DR  DS@ DT� DU� DW  DX9�DYy�DZ� D\  D]9�D^� D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�  D�� D�c3D�  D���D�<�D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D��D���D�` D�3D��3D�C3D�� D�|�D�  D�� D�` D�  D��3D�C3D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�c3D�  D�� D�@ D�� D��3D�  D�� D�c3D�3D��3D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�C3D�� D�� D�#3D��3D�` D���D�� D�@ D�� D��3D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� DŃ3D�#3D�� D�c3D�3Dȣ3D�C3D�� Dʀ D�  D��3D�` D�  D͠ D�@ D�� Dσ3D�  D��3D�` D�  DҠ D�C3D��3DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�c3D�  Dܠ D�C3D�� Dހ D�#3D��3D�` D���D� D�@ D���D� D�  D�� D�` D�3D�3D�@ D��3D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D��D�` D�3D��3D�<�D��3D��3D�#3D��fD�S3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�G�@�G�@�?}@�O�@�O�@�X@���@��@�hs@���@��#@��T@��#@��@��#@��#@��@��T@��@��T@���@��^@�@��#@��T@��@��@��T@��#@��#@��T@��T@��@��@��@��@���@���@�@�@���@��@��#@��#@��#@��T@��@��#@��#@��@��@��@��@��@��@��@��@��T@��T@��T@��@��@��@��T@��T@��T@���@���@��T@��@��@��@��@��@��@��T@��@��@��T@��@��@���@��T@��T@��#@��#@���@���@���@���@�@�@��^@��-@���@��-@���@���@���@���@���@��h@��h@��h@��h@��7@��7@��@�x�@�O�@�G�@�O�@�/@�7L@�/@���@���@�z�@�z�@�j@�bN@�bN@�bN@�bN@�bN@�bN@�Z@�Z@�j@�bN@�Z@�Z@�A�@�A�@�9X@�(�@�(�@�  @��m@��
@��w@�|�@�l�@�t�@�l�@�dZ@�dZ@�\)@�S�@�K�@�;d@�"�@�@��@��y@���@�~�@�V@��@�{@���@���@���@���@���@��^@���@�7L@�V@���@��9@���@��D@��F@�K�@�@��@���@��@��#@��#@�{@�V@�M�@�|�@�K�@��H@��R@���@�ff@�J@���@�@���@���@���@��9@��D@��j@�r�@�A�@�t�@��R@�J@���@��^@�x�@��@���@��@��D@�9X@��@��F@��@�"�@��y@�~�@�M�@�5?@�@��@��@��#@��^@��-@���@���@��@��@��@�x�@�X@�7L@�Ĝ@���@��u@���@���@���@���@���@���@��u@��D@�I�@��@���@�l�@�K�@�C�@�C�@�C�@�C�@�;d@�;d@�;d@�
=@��\@�n�@�ff@�ff@�=q@��@�@��@��7@�hs@��@���@�r�@�Q�@�1'@�w@\)@~��@~ȴ@~V@}�h@|�@|�j@|��@|Z@{ƨ@{dZ@{33@{��@|��@|I�@{t�@{@x�`@xbN@w�@wl�@w�w@w\)@xb@w�@w|�@w�@vȴ@vff@v5?@v$�@v@u�-@u/@t��@t��@t��@tI�@sƨ@r��@qG�@p�`@p�`@p�@p  @o�@nȴ@oK�@o�w@p1'@pr�@pr�@pr�@p�u@pA�@o|�@n�y@n��@n�y@n��@n��@oK�@n��@o��@p  @o�w@o�;@o�;@pbN@p�9@p�@p1'@pA�@p�@p��@pĜ@q�@q��@r=q@q�@p��@qx�@q��@r�H@r~�@s�m@u/@t�j@t�@uV@t��@u/@uV@t�@t�D@tz�@tj@tj@t�D@tI�@t�@uV@t�/@tj@t(�@s��@r�!@r�@r-@q��@q��@q��@qx�@qhs@qX@q7L@q&�@q%@p�`@pĜ@p�9@p�9@p��@pQ�@p �@pb@o�;@o�w@o��@o|�@ol�@ol�@o|�@o|�@oK�@o;d@o�@o
=@n�@nȴ@n��@n@m�T@m��@n@nE�@mp�@l��@l��@l�/@l��@lI�@k�
@kC�@j�!@j�\@j=q@j-@j-@j-@i��@h�`@hr�@hQ�@h  @g�P@g|�@gl�@gl�@gl�@g\)@fȴ@fv�@fV@f@e�h@e`B@eO�@e?}@d�@d(�@c��@ct�@cdZ@b��@bM�@b-@b=q@a��@ax�@`��@`Ĝ@`�9@`��@`��@`��@a&�@ax�@ax�@a%@a�@`Ĝ@_�@` �@_;d@^E�@^{@]��@]��@]p�@]O�@]p�@]��@^$�@^��@^v�@^$�@]O�@\(�@[�@["�@Z�H@Z��@[@Z�!@ZM�@Z�\@Z^5@Yx�@Y�@X��@X�@X�u@X  @W��@W�;@W�@W��@W��@W�@V��@VE�@U@U�@T�j@Tj@T9X@T�@S��@S��@St�@SS�@S@R�H@R�!@R�!@R��@Rn�@Q��@Q�^@Qx�@Q&�@P��@P�9@PQ�@P �@O�w@Ol�@O
=@O
=@N�@N�+@NE�@M�@M�-@M�@M`B@MO�@L��@L��@M`B@MO�@M�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@�G�@�G�@�?}@�O�@�O�@�X@���@��@�hs@���@��#@��T@��#@��@��#@��#@��@��T@��@��T@���@��^@�@��#@��T@��@��@��T@��#@��#@��T@��T@��@��@��@��@���@���@�@�@���@��@��#@��#@��#@��T@��@��#@��#@��@��@��@��@��@��@��@��@��T@��T@��T@��@��@��@��T@��T@��T@���@���@��T@��@��@��@��@��@��@��T@��@��@��T@��@��@���@��T@��T@��#@��#@���@���@���@���@�@�@��^@��-@���@��-@���@���@���@���@���@��h@��h@��h@��h@��7@��7@��@�x�@�O�@�G�@�O�@�/@�7L@�/@���@���@�z�@�z�@�j@�bN@�bN@�bN@�bN@�bN@�bN@�Z@�Z@�j@�bN@�Z@�Z@�A�@�A�@�9X@�(�@�(�@�  @��m@��
@��w@�|�@�l�@�t�@�l�@�dZ@�dZ@�\)@�S�@�K�@�;d@�"�@�@��@��y@���@�~�@�V@��@�{@���@���@���@���@���@��^@���@�7L@�V@���@��9@���@��D@��F@�K�@�@��@���@��@��#@��#@�{@�V@�M�@�|�@�K�@��H@��R@���@�ff@�J@���@�@���@���@���@��9@��D@��j@�r�@�A�@�t�@��R@�J@���@��^@�x�@��@���@��@��D@�9X@��@��F@��@�"�@��y@�~�@�M�@�5?@�@��@��@��#@��^@��-@���@���@��@��@��@�x�@�X@�7L@�Ĝ@���@��u@���@���@���@���@���@���@��u@��D@�I�@��@���@�l�@�K�@�C�@�C�@�C�@�C�@�;d@�;d@�;d@�
=@��\@�n�@�ff@�ff@�=q@��@�@��@��7@�hs@��@���@�r�@�Q�@�1'@�w@\)@~��@~ȴ@~V@}�h@|�@|�j@|��@|Z@{ƨ@{dZ@{33@{��@|��@|I�@{t�@{@x�`@xbN@w�@wl�@w�w@w\)@xb@w�@w|�@w�@vȴ@vff@v5?@v$�@v@u�-@u/@t��@t��@t��@tI�@sƨ@r��@qG�@p�`@p�`@p�@p  @o�@nȴ@oK�@o�w@p1'@pr�@pr�@pr�@p�u@pA�@o|�@n�y@n��@n�y@n��@n��@oK�@n��@o��@p  @o�w@o�;@o�;@pbN@p�9@p�@p1'@pA�@p�@p��@pĜ@q�@q��@r=q@q�@p��@qx�@q��@r�H@r~�@s�m@u/@t�j@t�@uV@t��@u/@uV@t�@t�D@tz�@tj@tj@t�D@tI�@t�@uV@t�/@tj@t(�@s��@r�!@r�@r-@q��@q��@q��@qx�@qhs@qX@q7L@q&�@q%@p�`@pĜ@p�9@p�9@p��@pQ�@p �@pb@o�;@o�w@o��@o|�@ol�@ol�@o|�@o|�@oK�@o;d@o�@o
=@n�@nȴ@n��@n@m�T@m��@n@nE�@mp�@l��@l��@l�/@l��@lI�@k�
@kC�@j�!@j�\@j=q@j-@j-@j-@i��@h�`@hr�@hQ�@h  @g�P@g|�@gl�@gl�@gl�@g\)@fȴ@fv�@fV@f@e�h@e`B@eO�@e?}@d�@d(�@c��@ct�@cdZ@b��@bM�@b-@b=q@a��@ax�@`��@`Ĝ@`�9@`��@`��@`��@a&�@ax�@ax�@a%@a�@`Ĝ@_�@` �@_;d@^E�@^{@]��@]��@]p�@]O�@]p�@]��@^$�@^��@^v�@^$�@]O�@\(�@[�@["�@Z�H@Z��@[@Z�!@ZM�@Z�\@Z^5@Yx�@Y�@X��@X�@X�u@X  @W��@W�;@W�@W��@W��@W�@V��@VE�@U@U�@T�j@Tj@T9X@T�@S��@S��@St�@SS�@S@R�H@R�!@R�!@R��@Rn�@Q��@Q�^@Qx�@Q&�@P��@P�9@PQ�@P �@O�w@Ol�@O
=@O
=@N�@N�+@NE�@M�@M�-@M�@M`B@MO�@L��@L��@M`B@MO�@M�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B�}B�}B��B�}B�}B��B�}B�}B��B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�wB�}B�}B�}B�}B�}B�}B�}B�wB�}B�wB�wB�}B�}B�}B�}B�wB�}B�wB�}B�}B�}B�wB�wB�}B�wB�wB�}B�wB�}B�wB�}B�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�qB�wB�wB�wB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�^B�^B�XB�XB�XB�XB�XB�RB�RB�RB�RB�RB�RB�LB�LB�LB�LB�FB�FB�?B�9B�9B�3B�3B�-B�'B�-B�3B�9B�?B�^B�^B�^B�^B�^B�XB�XB�XB�^B�XB�XB�FB�FB�FB�LB�FB�9B�3B�!B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�uB�oB�oB�uB��B��B�{B�uB�\B�VB�PB�VB�VB�VB�bB�bB�\B�\B�VB�VB�VB�PB�PB�PB�JB�JB�JB�DB�DB�7B�+B�B�B}�B�B�B�B�B�B�B�%B�+B�+B�1B�1B�1B�%B�%B�%B�+B�+B�+B�7B�7B�=B�DB�JB�JB�PB�VB�bB�bB�\B�hB�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�'B�3B�3B�9B�9B�3B�3B�-B�-B�-B�-B�-B�3B�3B�3B�9B�9B�9B�?B�?B�?B�?B�?B�FB�FB�FB�LB�LB�RB�RB�RB�RB�^B�dB�jB�qB�qB�qB�qB�qB�wB�wB�}B�}B��BB��B�}B��B��B��B��B��B��B�}B�}B�}B��B��B��B��B�wB�qB�qB�qB�wB�wB�wB�wB�wB�wB�wB�}B�}B�}B�}B�}B�}B�wB�wB�qB�qB�qB�qB�qB�qB�qB�qB�qB�jB�qB�qB�qB�qB�qB�}B��B��BB��BB��B��BB��B�wB�}B�}B�}B�}B�}B��B��BÖBƨBƨBƨBĜB��B��B��B�}B��B��B��B��B��B��B�}B�wB�wB�wB�wB�wB�qB�wB�}B�}B�}B�wB�wB�wB�qB�jB�jB�dB�dB�jB�jB�jB�qB�qB�qB�qB�wB�wB�}B�}B�wB�wB�}B�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�~B�~B�xB�xB�qB�kB�kB�eB�eB�_B�YB�_B�eB�kB�qB��B��B��B��B��B��B��B��B��B��B��B�xB�xB�xB�~B�xB�kB�eB�TB�GB�AB�AB�;B�;B�5B�5B�5B�/B�/B�)B�)B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�
B�
B�
B�
B�
B�
B�
B�
B�B�B�B�B�B�B�
B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�~B�xB�xB�kB�`B�TB�TBz)B�TB�NBGB~AB�NB�TB�ZB�`B�`B�fB�fB�fB�ZB�ZB�ZB�`B�`B�`B�kB�kB�qB�xB�~B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�#B�)B�/B�;B�;B�AB�;B�;B�AB�GB�NB�NB�ZB�fB�fB�lB�lB�fB�fB�`B�`B�`B�`B�`B�fB�fB�fB�lB�lB�lB�rB�rB�rB�rB�rB�yB�yB�yB�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0037123                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904432018110709044320181107090443  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172710  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172710  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090443  IP  PSAL            @��D�S3G�O�                