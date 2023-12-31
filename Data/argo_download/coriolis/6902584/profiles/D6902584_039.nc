CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:19Z last update (coriolis COQC software)   
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
_FillValue                 4  Bh   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  Mh   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Xh   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a4   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ch   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  nh   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w4   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �    PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �4   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �4   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �\   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �`   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �d   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �h   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �l   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �    SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �0   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �0   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �0   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090452  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               'A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�[Q33331   @�[Q��;�@P
����@�l�ô1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A��A!��A0  A@  AP  A`  Ap  A�  A�33A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  B   B  B  BffB  B  BffB  B   B$ffB(  B,  B0  B4  B8  B<  B?��BC��BH  BL  BP  BT  BXffB\ffB`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B���B�  B�33B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C/  C1��C4�C6��C9  C;� C>  C@� CC  CE��CH�CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ CŌ�C�� C�  C�L�Cʀ C˳3C�  C�@ Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�33Cހ C���C�  C�L�C� C�� C�  C�@ C� C�3C�  C�@ C�s3C�3C�  C�@ C� C�� C�  C�@ C�� C�� C�  C���C��D �fD  D9�D� D� D  D@ D	� D
� DfD@ D� D� D  DFfD� D� D  D9�D� D�fD  D@ Dy�D��D��D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9fD:FfD;�fD<�fD>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\fD]@ D^� D_� Da  Db@ Dc� Dd� DffDgFfDh�fDi� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� DzfD{FfD|�fD}�fD  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D���D�|�D��D�� D�` D�  DÜ�D�@ D�� Dŀ D�  D�� D�` D�  DȜ�D�@ D�� D�|�D�  D��3D�` D�  D͠ D�C3D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�3Dף3D�C3D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D��D�� D�` D�  D��D�@ D�� D�|�D�  D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�3D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�C3D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   A��A!��A0  A@  AP  A`  Ap  A�  A�33A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�33A�  A�  A�  B   B  B  BffB  B  BffB  B   B$ffB(  B,  B0  B4  B8  B<  B?��BC��BH  BL  BP  BT  BXffB\ffB`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�  B���B�  B�33B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C/  C1��C4�C6��C9  C;� C>  C@� CC  CE��CH�CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr��Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ CŌ�C�� C�  C�L�Cʀ C˳3C�  C�@ Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�33Cހ C���C�  C�L�C� C�� C�  C�@ C� C�3C�  C�@ C�s3C�3C�  C�@ C� C�� C�  C�@ C�� C�� C�  C���C��D �fD  D9�D� D� D  D@ D	� D
� DfD@ D� D� D  DFfD� D� D  D9�D� D�fD  D@ Dy�D��D��D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9fD:FfD;�fD<�fD>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\fD]@ D^� D_� Da  Db@ Dc� Dd� DffDgFfDh�fDi� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� DzfD{FfD|�fD}�fD  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D���D�|�D��D�� D�` D�  DÜ�D�@ D�� Dŀ D�  D�� D�` D�  DȜ�D�@ D�� D�|�D�  D��3D�` D�  D͠ D�C3D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�3Dף3D�C3D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D��D�� D�` D�  D��D�@ D�� D�|�D�  D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�3D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�C3D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@��#@��#@��#@��#@��T@��#@���@���@��#@���@��#@���@�@��^@�@�@���@���@�Q�@�"�@��@��@�9@�(�@��@�A�@�ƨ@�n�@���@�E�@���@�\)@��@�^5@�Q�@�/@�|�@�\)@���@�hs@�O�@�&�@ش9@��
@��y@�~�@թ�@�Ĝ@��y@�O�@Л�@�~�@�?}@Ͳ-@�=q@Ώ\@ΰ!@��@�ȴ@�-@�@�$�@�ȴ@�ȴ@Ώ\@��@�G�@���@�%@̣�@̬@̛�@̋D@�bN@�(�@�ƨ@˝�@�+@�~�@�ff@ʇ+@�n�@�-@Ɂ@ɉ7@Ɂ@�p�@�O�@�O�@�7L@��@���@���@ȼj@ȋD@�Z@�(�@�b@���@��@�dZ@��@���@�@ģ�@��m@�K�@�V@��@���@�ff@�ff@�(�@�;d@���@�@�
=@�E�@�X@��@�z�@�1'@��
@�t�@�"�@��@�ȴ@��R@�v�@�ff@�n�@�v�@�n�@�M�@�~�@�V@�=q@�5?@���@���@���@��@���@��u@�r�@�Q�@� �@�9X@�A�@�|�@��@�Z@�Ĝ@�Q�@��;@�J@�5?@�^5@��@�n�@��H@���@�X@�(�@�(�@�A�@��@���@�t�@��@�C�@�l�@�S�@��+@��j@�ƨ@�
=@���@�ȴ@�ȴ@���@�v�@�V@��-@��@�Ĝ@��j@��9@�V@�X@���@��@�Z@�1'@�bN@�A�@�b@��m@��w@��P@�dZ@�33@�o@��@�o@��@��\@�n�@�E�@��@�{@��#@��h@�hs@�?}@��@��@�Q�@�9X@��m@��w@�33@��y@���@��@���@�@��7@�X@��@��D@�A�@��@��;@�ƨ@��@���@���@���@�v�@�-@���@���@�`B@��9@�  @�ƨ@��@���@��@�t�@�\)@�S�@�C�@�;d@�+@��H@���@�{@���@�{@���@���@���@�9X@�l�@�o@��@��@��+@�-@��7@�&�@��j@�j@�1'@���@��w@��@�@���@�^5@�5?@��T@��h@��@���@��@�9X@���@�ƨ@��F@�dZ@�ȴ@�ff@�J@�@��@�?}@�&�@��`@��@� �@��w@�dZ@�C�@�33@��@��+@��T@��@��#@���@�G�@��@���@��`@�Ĝ@���@��D@�9X@���@�|�@�C�@��H@��R@���@�M�@�{@���@���@���@�p�@�&�@���@��@��`@���@�Ĝ@��@�(�@��
@��@�l�@�33@��@���@���@�n�@�V@�=q@�@���@���@�@���@�@���@���@�X@�?}@�&�@�&�@�V@���@���@���@��@��j@�bN@�w@�@
=@~�+@~��@~�y@;d@�P@|�@�P@|�@�@~�@~�R@~��@~�+@~�+@~5?@}�T@}�T@}�@~@}�T@}��@}�h@}p�@}V@}?}@}V@|�/@|Z@{�
@{��@z�!@zM�@y��@y%@w��@w�P@w��@w�w@w�@wl�@wK�@w�@v�y@v��@v�@w�@wK�@w\)@wl�@w\)@w\)@w;d@w+@w+@w+@w+@v��@v��@v��@vȴ@vȴ@v�R@v�R@v�R@v��@v��@v��@vv�@vff@vE�@vE�@vV@vff@vff@vff@vff@vE�@v{@v@v@u�@u�T@u�-@v{@vv�@v��@vȴ@v�@v��@vV@v5?@v$�@v{@u�@u��@u�h@up�@u�h@u�@up�@u/@t��@t�@tz�@s��@s��@s��@sS�@s@r~�@rJ@qx�@q&�@pĜ@pr�@p1'@pb@o��@o��@oK�@o
=@nV@n@m��@m�h@m/@l�j@l(�@k�F@k@j�H@j�!@jn�@j^5@i�#@i��@i%@h��@hbN@g�w@g\)@f{@e/@d�@dz�@dj@c��@c@b-@a�@a�^@ahs@`�`@`Q�@_�@_+@^ff@\I�@Z�@ZJ@Yx�@Yx�@Yx�@Yhs@Y��@Y��@YX@X��@XQ�@W��@Vff@U�T@U��@U?}@T�/@TZ@T�@S��@S�@S@R��@R^5@R�@QX@Q%@P�`@O�w@N@L�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @���@��#@��#@��#@��#@��T@��#@���@���@��#@���@��#@���@�@��^@�@�@���@���@�Q�@�"�@��@��@�9@�(�@��@�A�@�ƨ@�n�@���@�E�@���@�\)@��@�^5@�Q�@�/@�|�@�\)@���@�hs@�O�@�&�@ش9@��
@��y@�~�@թ�@�Ĝ@��y@�O�@Л�@�~�@�?}@Ͳ-@�=q@Ώ\@ΰ!@��@�ȴ@�-@�@�$�@�ȴ@�ȴ@Ώ\@��@�G�@���@�%@̣�@̬@̛�@̋D@�bN@�(�@�ƨ@˝�@�+@�~�@�ff@ʇ+@�n�@�-@Ɂ@ɉ7@Ɂ@�p�@�O�@�O�@�7L@��@���@���@ȼj@ȋD@�Z@�(�@�b@���@��@�dZ@��@���@�@ģ�@��m@�K�@�V@��@���@�ff@�ff@�(�@�;d@���@�@�
=@�E�@�X@��@�z�@�1'@��
@�t�@�"�@��@�ȴ@��R@�v�@�ff@�n�@�v�@�n�@�M�@�~�@�V@�=q@�5?@���@���@���@��@���@��u@�r�@�Q�@� �@�9X@�A�@�|�@��@�Z@�Ĝ@�Q�@��;@�J@�5?@�^5@��@�n�@��H@���@�X@�(�@�(�@�A�@��@���@�t�@��@�C�@�l�@�S�@��+@��j@�ƨ@�
=@���@�ȴ@�ȴ@���@�v�@�V@��-@��@�Ĝ@��j@��9@�V@�X@���@��@�Z@�1'@�bN@�A�@�b@��m@��w@��P@�dZ@�33@�o@��@�o@��@��\@�n�@�E�@��@�{@��#@��h@�hs@�?}@��@��@�Q�@�9X@��m@��w@�33@��y@���@��@���@�@��7@�X@��@��D@�A�@��@��;@�ƨ@��@���@���@���@�v�@�-@���@���@�`B@��9@�  @�ƨ@��@���@��@�t�@�\)@�S�@�C�@�;d@�+@��H@���@�{@���@�{@���@���@���@�9X@�l�@�o@��@��@��+@�-@��7@�&�@��j@�j@�1'@���@��w@��@�@���@�^5@�5?@��T@��h@��@���@��@�9X@���@�ƨ@��F@�dZ@�ȴ@�ff@�J@�@��@�?}@�&�@��`@��@� �@��w@�dZ@�C�@�33@��@��+@��T@��@��#@���@�G�@��@���@��`@�Ĝ@���@��D@�9X@���@�|�@�C�@��H@��R@���@�M�@�{@���@���@���@�p�@�&�@���@��@��`@���@�Ĝ@��@�(�@��
@��@�l�@�33@��@���@���@�n�@�V@�=q@�@���@���@�@���@�@���@���@�X@�?}@�&�@�&�@�V@���@���@���@��@��j@�bN@�w@�@
=@~�+@~��@~�y@;d@�P@|�@�P@|�@�@~�@~�R@~��@~�+@~�+@~5?@}�T@}�T@}�@~@}�T@}��@}�h@}p�@}V@}?}@}V@|�/@|Z@{�
@{��@z�!@zM�@y��@y%@w��@w�P@w��@w�w@w�@wl�@wK�@w�@v�y@v��@v�@w�@wK�@w\)@wl�@w\)@w\)@w;d@w+@w+@w+@w+@v��@v��@v��@vȴ@vȴ@v�R@v�R@v�R@v��@v��@v��@vv�@vff@vE�@vE�@vV@vff@vff@vff@vff@vE�@v{@v@v@u�@u�T@u�-@v{@vv�@v��@vȴ@v�@v��@vV@v5?@v$�@v{@u�@u��@u�h@up�@u�h@u�@up�@u/@t��@t�@tz�@s��@s��@s��@sS�@s@r~�@rJ@qx�@q&�@pĜ@pr�@p1'@pb@o��@o��@oK�@o
=@nV@n@m��@m�h@m/@l�j@l(�@k�F@k@j�H@j�!@jn�@j^5@i�#@i��@i%@h��@hbN@g�w@g\)@f{@e/@d�@dz�@dj@c��@c@b-@a�@a�^@ahs@`�`@`Q�@_�@_+@^ff@\I�@Z�@ZJ@Yx�@Yx�@Yx�@Yhs@Y��@Y��@YX@X��@XQ�@W��@Vff@U�T@U��@U?}@T�/@TZ@T�@S��@S�@S@R��@R^5@R�@QX@Q%@P�`@O�w@N@L�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB1'B1'B1'B0!B1'B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B/B/B.B,B.B-B+B)�B-B,B+B(�B'�B&�B&�B1'B/B2-B.B2-B8RB=qB<jB;dB?}BG�BF�BF�BG�BI�BH�BI�BJ�BL�BO�BK�BH�BH�BE�BJ�BO�BS�BT�BVBXBYBXBYB]/B^5B_;B]/B\)B[#B\)B[#B\)B]/B]/B\)B\)B[#BYBYBXBYBZB[#BZBYB]/B]/B^5B^5B^5B]/B^5B]/B]/B]/B]/B]/B\)B]/B]/B^5B]/BYBM�BR�BYBXBT�BQ�BI�BF�BC�BB�B=qB:^B9XB9XB8RB6FB5?B5?B5?B6FB7LB7LB7LB7LB6FB6FB6FB6FB7LB8RB8RB9XB:^B9XB9XB9XB9XB8RB8RB8RB8RB7LB7LB7LB8RB9XB;dB9XB=qB?}BB�BA�B?}B:^B:^B;dB;dB=qB>wB;dB8RB49B5?B6FB5?B6FB5?B33B6FB6FB5?B2-B,B(�B%�B$�B$�B%�B%�B$�B$�B"�B �B�B �B �B"�B$�B#�B"�B"�B"�B#�B#�B"�B"�B!�B!�B �B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBuBoBoBbB\BVBPBJBJBDB
=B	7B1B+B+B%B%BBBBB  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�sB�sB�mB�fB�`B�ZB�TB�HB�HB�BB�;B�5B�)B�)B�#B�B�B�B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBƨBǮBǮBǮBƨBŢBŢBŢBĜBĜBĜBB��B��B�}B�wB�qB�qB�jB�jB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�XB�RB�RB�LB�LB�FB�FB�FB�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�FB�FB�FB�FB�FB�FB�FB�LB�LB�RB�RB�RB�LB�?B�FB�?B�9B�FB�LB�RB�XB�^B�^B�^B�^B�^B�dB�dB�dB�jB�jB�jB�jB�jB�qB�wB�wB�wB�wB�qB�wB�qB�jB�jB�dB�^B�XB�XB�RB�?B�3B�-B�-B�9B�9B�3B�9B�9B�9B�9B�FB�LB�RB�^B�dB�jB�qB�qB�qB�wB�wB�wB�}B�}B��B��B��B��B��B��B��BBBÖBĜBĜBŢBŢBƨBƨBǮBǮBǮBȴBȴBȴBȴBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBƨBƨBƨBƨBŢBŢBŢBĜBÖB��B�}B�wB�}B��BBÖBĜBĜBĜBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBBBBB��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B(XB(XB(XB'RB(XB'RB'RB'RB'RB'RB'RB'RB'RB'RB'RB&LB&LB%EB#9B%EB$?B"4B!.B$@B#:B"4B (B"BBB(YB&MB)_B%FB)_B/�B4�B3�B2�B6�B>�B=�B=�B>�B@�B?�B@�BA�BC�BGBB�B?�B?�B<�BA�BGBK)BL/BM5BOABPHBOABPHBT`BUfBVlBT`BSZBRTBSZBRTBSZBT`BT`BSZBSZBRTBPHBPHBOABPHBQNBRTBQNBPHBT`BT`BUfBUfBUfBT`BUfBT`BT`BT`BT`BT`BSZBT`BT`BUfBT`BPHBEBJ$BPIBOBBL0BIB@�B=�B:�B9�B4�B1�B0�B0�B/�B-zB,sB,sB,tB-zB.�B.�B.�B.�B-{B-{B-{B-{B.�B/�B/�B0�B1�B0�B0�B0�B0�B/�B/�B/�B/�B.�B.�B.�B/�B0�B2�B0�B4�B6�B9�B8�B6�B1�B1�B2�B2�B4�B5�B2�B/�B+nB,tB-{B,tB-{B,tB*hB-{B-{B,tB)bB#>B ,BBBBBBBBB�B�B�B�BBBBBBBBBBBBB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B
�B
�B	�B	�B�B�B�B�B�B�B}BvB pB�jB�dB�dB�^B�^B�XB�RB�LB�FB�:B�4B�.B�.B�.B�.B�(B�(B�(B�(B�(B�!B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B�B߯B߯BީBݢBܜBۖBڐB؅B؅B�B�xB�rB�fB�fB�`B�ZB�TB�TB�MB�GB�BB�<B�6B�0B�0B�#B�#B�B�B�B�B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�tB�nB�nB�zB�zB�tB�zB�zB�zB�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�$B�$B�+B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�,B�,B�,B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0085661                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904522018110709045220181107090452  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172719  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172719  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090452  IP  PSAL            @ffD��3G�O�                