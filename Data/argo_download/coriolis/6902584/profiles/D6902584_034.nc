CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  1   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:17Z last update (coriolis COQC software)   
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
_FillValue                 4  B`   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  MX   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  XP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  cH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  n@   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �$   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �d   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �t   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �x   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090449  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               "A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�N�OĄ`1   @�N���?@O� �ܾ��@��R�P1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B��B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BXffB\  B_��Bd  Bh  Bl  Bp  BtffBx  B|  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B���B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  C�fCffC�fC	� C  C� C  C� C�C� C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH�CJ� CM  CO� CR  CT� CW�CY� C\  C^� Ca  CcffCe�fCh� Ck  Cm� Cp  Cr� Cu  Cw� Cz�C|� C~�fC�� C�  C�@ C�� C��3C�  C�@ C���C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C��3C�  C�L�C���C���C��C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C�  C�L�Cό�C���C�  C�@ CԀ C�� C�  C�L�Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C��D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D�fDfDFfD� D� D  D@ D�fD� D��D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(�fD*  D+9�D,� D-�fD/fD0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>  D?@ D@� DA� DC  DD@ DE�fDF� DH  DI@ DJ� DK� DMfDNFfDO�fDP� DQ��DS9�DT� DU� DW  DX@ DY� DZ�fD\  D]@ D^� D_� Da  DbFfDc�fDd�fDf  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� Dp  DqFfDr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�C3D��3D��3D�  D�� D�` D���D���D�@ D��3D��3D�  D���D�\�D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D��3D�@ D���D�� D�  D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  Dà D�C3D�� Dŀ D�  DƼ�D�\�D���DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�  D�� D�` D���DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D��3D�` D�  D��D�<�D���D� D�  D�� D�` D���D���D�@ D�� D�� D�  D��3D�` D�3D��fD�<�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B��B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BXffB\  B_��Bd  Bh  Bl  Bp  BtffBx  B|  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B���B�  B�  B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  C�fCffC�fC	� C  C� C  C� C�C� C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH�CJ� CM  CO� CR  CT� CW�CY� C\  C^� Ca  CcffCe�fCh� Ck  Cm� Cp  Cr� Cu  Cw� Cz�C|� C~�fC�� C�  C�@ C�� C��3C�  C�@ C���C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C��3C�  C�L�C���C���C��C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C���C�  C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C�  C�L�Cό�C���C�  C�@ CԀ C�� C�  C�L�Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C��D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D�fDfDFfD� D� D  D@ D�fD� D��D@ D� D� D   D!@ D"� D#� D%  D&@ D'� D(�fD*  D+9�D,� D-�fD/fD0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<�fD>  D?@ D@� DA� DC  DD@ DE�fDF� DH  DI@ DJ� DK� DMfDNFfDO�fDP� DQ��DS9�DT� DU� DW  DX@ DY� DZ�fD\  D]@ D^� D_� Da  DbFfDc�fDd�fDf  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� Dp  DqFfDr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�C3D��3D��3D�  D�� D�` D���D���D�@ D��3D��3D�  D���D�\�D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D��3D�@ D���D�� D�  D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  Dà D�C3D�� Dŀ D�  DƼ�D�\�D���DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�  D�� D�` D���DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D��3D�` D�  D��D�<�D���D� D�  D�� D�` D���D���D�@ D�� D�� D�  D��3D�` D�3D��fD�<�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�
=@�+@�33@�|�@�dZ@�dZ@�33@��@�bN@�I�@�1'@���@���@���@�&�@�G�@�X@�O�@�7L@�&�@��@��@��@�%@�%@�%@���@���@��u@��D@�r�@�bN@�I�@�1'@��@�1@�  @�  @���@��@��m@��m@��m@��m@��
@�ƨ@�ƨ@��w@��F@��@��F@��F@��F@���@���@���@��P@��P@��P@��P@��@��@�|�@�|�@�t�@�|�@�|�@�t�@�t�@�l�@�t�@�t�@�l�@�l�@�|�@�|�@��@�|�@�t�@�t�@�l�@�S�@�"�@��@�ȴ@���@�M�@��@���@��@�&�@�%@�Ĝ@��@�9X@�(�@�b@���@�1@��\@�X@�?}@�7L@�&�@�G�@�/@�/@���@���@��@��@�1'@��@�t�@��@��!@�ff@���@�@�@�X@��/@�Ĝ@���@���@��9@�9X@�1@�o@��@���@���@�~�@�^5@�{@���@���@��7@�p�@�x�@�p�@�p�@�p�@�`B@�X@�X@�X@�`B@�`B@�X@�`B@�`B@�X@��9@� �@��
@��@���@��@��@���@�  @�1'@�A�@���@��7@�&�@���@���@��9@��u@�z�@�j@�Q�@�1'@�(�@���@��@��@���@� �@�A�@�A�@��@��@��@�l�@�l�@���@��!@��\@�v�@�ff@�5?@�J@�J@�{@���@��#@���@�p�@�X@�G�@���@�Ĝ@���@�z�@��m@���@�ƨ@�ƨ@��F@��@��P@�dZ@��@�@��y@��@�ȴ@���@�~�@�^5@�{@��h@�X@�7L@��@�%@���@��j@�Z@��@��@��
@���@�|�@�\)@�;d@�
=@���@�^5@���@��h@�X@�/@�7L@��@�j@���@��;@���@��@��@�\)@�K�@��@��!@�~�@�-@��#@��^@��h@�X@��@��u@��m@��F@�|�@�33@��@��@�v�@�5?@�E�@�=q@�J@��^@��@�hs@��@���@���@���@��h@�`B@��@���@��`@�Ĝ@���@�Q�@�b@��w@���@� �@�I�@�I�@�b@�  @��m@���@�S�@�
=@��R@�E�@�p�@�/@�hs@�%@���@�S�@�ȴ@�5?@��7@�G�@���@�Ĝ@��@�9X@�1@�1@�  @�@\)@~v�@}�@}p�@|��@}?}@~E�@~�@~�+@
=@~V@}�T@}?}@|��@{�F@{�F@|(�@|��@{ƨ@{�@|�@|��@|��@|I�@{dZ@z-@y7L@x�9@x�9@x��@yx�@yx�@yhs@x�@xb@x1'@x��@y&�@y�7@{33@{"�@z~�@z�\@z�!@z=q@zM�@z�\@z�@{C�@{C�@{��@|1@|I�@|j@}�@}�-@}@}��@}�T@}�@|�@{S�@y�#@yG�@y7L@x��@x�u@xr�@xbN@xr�@x�9@y%@x��@xA�@x  @w�@w��@w��@w��@w�;@w�;@w�@w�@w�@w�@w�@w�;@xb@x �@xb@w�@w�@xA�@xQ�@xr�@x�u@x�u@x��@y%@yx�@yX@x�u@w�@v��@v��@v��@vȴ@v��@v��@v�y@vȴ@vv�@vff@v@u�h@u?}@t��@t��@t�/@t�j@tz�@t9X@t(�@t9X@t�@s�m@sƨ@s�F@sC�@s@r�H@r�\@r��@so@so@s@s"�@s"�@r�@r�!@r-@q�#@q��@q��@q��@q��@qhs@p��@p��@pbN@o|�@o�@n�@n��@n�+@nff@n5?@n5?@n{@m�@m?}@m/@l��@l��@lj@lI�@k��@k�@k��@kS�@j^5@i��@i�7@iX@iG�@h��@h1'@h �@g�P@g�@f�R@fE�@e�T@e�-@e?}@d��@d�j@dj@dZ@d(�@c�F@c��@ct�@cdZ@c33@b��@b~�@b=q@bJ@a�^@ahs@a%@`Ĝ@`bN@` �@_�@_��@_;d@^�R@^��@^V@^{@]�@]�T@]p�@\�/@\�D@\9X@[��@[�@[S�@[33@Z��@Z~�@Z�@Y��@YX@X�`@X�9@X�u@Xr�@XQ�@Xb@W�@W�P@WK�@W+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�
=@�+@�33@�|�@�dZ@�dZ@�33@��@�bN@�I�@�1'@���@���@���@�&�@�G�@�X@�O�@�7L@�&�@��@��@��@�%@�%@�%@���@���@��u@��D@�r�@�bN@�I�@�1'@��@�1@�  @�  @���@��@��m@��m@��m@��m@��
@�ƨ@�ƨ@��w@��F@��@��F@��F@��F@���@���@���@��P@��P@��P@��P@��@��@�|�@�|�@�t�@�|�@�|�@�t�@�t�@�l�@�t�@�t�@�l�@�l�@�|�@�|�@��@�|�@�t�@�t�@�l�@�S�@�"�@��@�ȴ@���@�M�@��@���@��@�&�@�%@�Ĝ@��@�9X@�(�@�b@���@�1@��\@�X@�?}@�7L@�&�@�G�@�/@�/@���@���@��@��@�1'@��@�t�@��@��!@�ff@���@�@�@�X@��/@�Ĝ@���@���@��9@�9X@�1@�o@��@���@���@�~�@�^5@�{@���@���@��7@�p�@�x�@�p�@�p�@�p�@�`B@�X@�X@�X@�`B@�`B@�X@�`B@�`B@�X@��9@� �@��
@��@���@��@��@���@�  @�1'@�A�@���@��7@�&�@���@���@��9@��u@�z�@�j@�Q�@�1'@�(�@���@��@��@���@� �@�A�@�A�@��@��@��@�l�@�l�@���@��!@��\@�v�@�ff@�5?@�J@�J@�{@���@��#@���@�p�@�X@�G�@���@�Ĝ@���@�z�@��m@���@�ƨ@�ƨ@��F@��@��P@�dZ@��@�@��y@��@�ȴ@���@�~�@�^5@�{@��h@�X@�7L@��@�%@���@��j@�Z@��@��@��
@���@�|�@�\)@�;d@�
=@���@�^5@���@��h@�X@�/@�7L@��@�j@���@��;@���@��@��@�\)@�K�@��@��!@�~�@�-@��#@��^@��h@�X@��@��u@��m@��F@�|�@�33@��@��@�v�@�5?@�E�@�=q@�J@��^@��@�hs@��@���@���@���@��h@�`B@��@���@��`@�Ĝ@���@�Q�@�b@��w@���@� �@�I�@�I�@�b@�  @��m@���@�S�@�
=@��R@�E�@�p�@�/@�hs@�%@���@�S�@�ȴ@�5?@��7@�G�@���@�Ĝ@��@�9X@�1@�1@�  @�@\)@~v�@}�@}p�@|��@}?}@~E�@~�@~�+@
=@~V@}�T@}?}@|��@{�F@{�F@|(�@|��@{ƨ@{�@|�@|��@|��@|I�@{dZ@z-@y7L@x�9@x�9@x��@yx�@yx�@yhs@x�@xb@x1'@x��@y&�@y�7@{33@{"�@z~�@z�\@z�!@z=q@zM�@z�\@z�@{C�@{C�@{��@|1@|I�@|j@}�@}�-@}@}��@}�T@}�@|�@{S�@y�#@yG�@y7L@x��@x�u@xr�@xbN@xr�@x�9@y%@x��@xA�@x  @w�@w��@w��@w��@w�;@w�;@w�@w�@w�@w�@w�@w�;@xb@x �@xb@w�@w�@xA�@xQ�@xr�@x�u@x�u@x��@y%@yx�@yX@x�u@w�@v��@v��@v��@vȴ@v��@v��@v�y@vȴ@vv�@vff@v@u�h@u?}@t��@t��@t�/@t�j@tz�@t9X@t(�@t9X@t�@s�m@sƨ@s�F@sC�@s@r�H@r�\@r��@so@so@s@s"�@s"�@r�@r�!@r-@q�#@q��@q��@q��@q��@qhs@p��@p��@pbN@o|�@o�@n�@n��@n�+@nff@n5?@n5?@n{@m�@m?}@m/@l��@l��@lj@lI�@k��@k�@k��@kS�@j^5@i��@i�7@iX@iG�@h��@h1'@h �@g�P@g�@f�R@fE�@e�T@e�-@e?}@d��@d�j@dj@dZ@d(�@c�F@c��@ct�@cdZ@c33@b��@b~�@b=q@bJ@a�^@ahs@a%@`Ĝ@`bN@` �@_�@_��@_;d@^�R@^��@^V@^{@]�@]�T@]p�@\�/@\�D@\9X@[��@[�@[S�@[33@Z��@Z~�@Z�@Y��@YX@X�`@X�9@X�u@Xr�@XQ�@Xb@W�@W�P@WK�@W+111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB,B+B,B+B)�B(�B'�B'�B&�B%�B&�B(�B+B,B,B-B.B/B/B/B/B/B/B/B.B.B/B/B/B/B/B.B/B/B/B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B.B-B.B.B-B-B-B-B,B+B+B+B)�B)�B(�B'�B'�B'�B&�B%�B$�B$�B$�B"�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuB{BuBhBhBbBbBbBVBVBDB
=B	7B	7B1B1B+B%BBBBBBBBBBBB%B%BB%BBBBBB  B  B  B  BBBBB1B
=B	7B1B1B+B+B%B%B%B%BBBBBB%B%B%BBBBBBB  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�sB�sB�mB�mB�fB�fB�fB�`B�ZB�TB�HB�;B�5B�5B�/B�5B�/B�#B�B�B�B�B�B�B�B�B�
B�
B�B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺB��BɺBɺBȴBǮBǮBȴBɺBɺBɺBȴBȴBǮBǮBǮBƨBƨBĜBÖBBÖBŢBƨBƨBŢBŢBĜBÖBB��B�}B�qB�^B�XB�dB�^B�RB�9B�-B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B��B��B��B��B��B�B��B��B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�'B�-B�'B�-B�3B�3B�9B�?B�FB�LB�RB�XB�^B�dB�jB�}B��B��B��BBÖB��B�wB�dB�dB�dB�dB�^B�dB�jB�qB�wB��B��B�}B�wB�wB�}B�}B��B��B��B��B��B��B��B��BÖBĜBĜBĜBŢBƨBǮBȴBȴBȴBɺB��B��B��B��B��B��BɺBȴBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺB��B��BɺBȴBȴBǮBǮBƨBƨBƨBƨBƨBƨBƨBƨBŢBŢBŢBƨBƨBƨBƨBƨBƨBŢBŢBŢBŢBŢBŢBĜBŢBŢBŢBŢBĜBĜBĜBĜBĜBĜBÖBBBBBBBBBB��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B$�B#�B$�B#�B"�B!�B �B �B�B�B�B!�B#�B$�B$�B%�B&�B'�B'�B'�B'�B'�B'�B'�B&�B&�B'�B'�B'�B'�B'�B&�B'�B'�B'�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B%�B&�B&�B%�B%�B%�B%�B$�B#�B#�B#�B"�B"�B!�B �B �B �B�B�B�B�B�B~BeBMBSBYB_BeBkBkBeBeBeB_BYBSBGBAB5B5B)B#B)B#B
B
B	B	B	BBB�B�B�B�B �B �B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B �B�B�B �B �B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�yB�lB�lB�lB�lB�lB�lB�fB�fB�`B�`B�ZB�ZB�ZB�TB�TB�OB�IB�BB�<B�<B�6B�6B�6B�0B�$B�$B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BϼBϼBζBͰB̫B̫B˥BʟBȒBǌBƆBŀBŀB�zB�zB�mB�mB�tB�mB�mB�gB�bB�bB�gB�mB�mB�mB�gB�gB�bB�bB�bB�\B�\B�PB�JB�CB�JB�VB�\B�\B�VB�VB�PB�JB�CB�7B�1B�%B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�2B�>B�>B�>B�CB�JB�>B�,B�B�B�B�B�B�B�B�&B�,B�8B�8B�2B�,B�,B�2B�2B�8B�8B�>B�>B�>B�>B�>B�>B�KB�QB�QB�QB�VB�\B�bB�hB�hB�hB�nB�uB�{BƇBƇBŁB�uB�nB�hB�nB�uB�{B�{BŁB�{BŁB�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�uB�uB�{B�{BŁBƇBǍBǍBȓBȓBəBəBȓBȓBȓBȓBəBəBəBȓBȓBǍBƇBƇBƇBƇBƇBƇBƇBƇBƇBƇBŁBŁBŁBŁBŁBŁBŁB�|BŁB�|B�vB�oB�oB�vB�vB�oB�iB�iB�cB�cB�]B�]B�]B�]B�]B�]B�]B�]B�WB�WB�WB�]B�]B�]B�]B�]B�]B�WB�WB�WB�WB�WB�WB�QB�WB�WB�WB�WB�QB�QB�QB�QB�QB�QB�KB�EB�EB�EB�EB�EB�EB�EB�EB�EB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0071346                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904492018110709044920181107090449  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172717  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172717  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090449  IP  PSAL            @ffD�<�G�O�                