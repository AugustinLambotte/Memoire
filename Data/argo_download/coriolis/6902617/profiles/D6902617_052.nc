CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-10-26T19:00:20Z creation; 2018-10-26T19:00:46Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181026190020  20181119104028  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               4A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @��@��1   @��J�Y@S����@�<�z� 1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @   @@  @�  @�  @�  @���@���AffA   A0  A@  AP  A`  Ap  A���A�  A�  A�  A�  A�  A�  A�33A�33A�  A�  A�  A���A�  A�  A�  B   B  B  B  B  B  B  BffB   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bw��B|  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	��C  C� C  C� C�C� C  C� C   C"��C%  C'� C*  C,� C/  C1� C4  C6� C9�C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  CcffCf  Ch��Ck  Cm� Cp  CrffCu  Cw� Cz  C|� C  C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C��3C��3C�@ C���C�� C�  C�33C�s3C��3C��3C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��C�@ Cπ C�� C�  C�L�CԌ�C�� C��3C�@ Cـ C���C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D � D  DFfD� D� D  D@ D	�fD
� D  D@ D� D� D  D@ D� D�fD  D@ D� D� D  D9�D� D�fD   D!@ D"� D#� D%  D&9�D'� D(�fD*  D+@ D,�fD-� D/  D0@ D1�fD2� D4  D5@ D6� D7� D9  D:9�D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc�fDd� Df  DgFfDh� Di�fDkfDl@ Dm� Dn� Dp  Dq@ Dry�Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�#3D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D���D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�#3D��3D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D��D�� D�` D�  Dܠ D�@ D�� Dހ D�  D��3D�c3D�  D��D�@ D�� D� D�  D�� D�c3D�3D� D�C3D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D��D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D��3D�� D��D�� D�` D�3D�� D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @   @@  @�  @�  @�  @���@���AffA   A0  A@  AP  A`  Ap  A���A�  A�  A�  A�  A�  A�  A�33A�33A�  A�  A�  A���A�  A�  A�  B   B  B  B  B  B  B  BffB   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bw��B|  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	��C  C� C  C� C�C� C  C� C   C"��C%  C'� C*  C,� C/  C1� C4  C6� C9�C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  CcffCf  Ch��Ck  Cm� Cp  CrffCu  Cw� Cz  C|� C  C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C��3C��3C�@ C���C�� C�  C�33C�s3C��3C��3C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��C�@ Cπ C�� C�  C�L�CԌ�C�� C��3C�@ Cـ C���C�  C�@ Cހ C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D � D  DFfD� D� D  D@ D	�fD
� D  D@ D� D� D  D@ D� D�fD  D@ D� D� D  D9�D� D�fD   D!@ D"� D#� D%  D&9�D'� D(�fD*  D+@ D,�fD-� D/  D0@ D1�fD2� D4  D5@ D6� D7� D9  D:9�D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc�fDd� Df  DgFfDh� Di�fDkfDl@ Dm� Dn� Dp  Dq@ Dry�Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�#3D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D���D�` D�3D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�#3D��3D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dף3D�@ D�� Dـ D��D�� D�` D�  Dܠ D�@ D�� Dހ D�  D��3D�c3D�  D��D�@ D�� D� D�  D�� D�c3D�3D� D�C3D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D��D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D��3D�� D��D�� D�` D�3D�� D�@ D��3D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�`B@�p�@�x�@��@��7@�p�@�`B@�O�@�%@��j@��/@�%@��@���@�  @��@��;@��
@���@��F@��@���@���@���@��@�dZ@�\)@�@���@��+@�v�@�E�@�J@�@��h@���@���@��h@���@�p�@�hs@��/@���@��u@�j@� �@���@�ƨ@�dZ@��@��@���@�ȴ@���@���@�n�@�^5@�^5@�^5@�^5@��#@�x�@�X@��`@���@��D@�z�@�1@��;@��F@��P@�K�@�K�@��@�o@��@��@��@�K�@�33@���@��!@�M�@�M�@�-@�$�@�$�@�=q@�M�@�ff@�^5@�ff@�V@�M�@�M�@�@���@��@�p�@��/@��/@�9X@���@�l�@�;d@���@�=q@��@��-@��7@��@���@���@��@��m@�C�@���@���@�~�@�=q@��@��h@�G�@�&�@�V@���@�r�@�  @�;d@��!@�$�@��@�@���@��@�?}@���@���@�z�@�1@}�@}�h@|��@|9X@{dZ@z��@yx�@yG�@y�@y�@x�9@x  @w��@w+@vV@u��@t�D@r��@p��@m�-@kƨ@j�@h��@g\)@g
=@f�@f�@ep�@a%@]�h@\�/@\j@W�P@SC�@PQ�@N5?@J�@IG�@G�w@FV@D��@@�9@<�@;�
@9��@81'@5�T@5V@3o@2�\@0�9@-p�@*�!@'��@%?}@!�#@;d@�@��@�@�m@b@�@
M�@ �@�D?��H?�!?�D?�K�?��?���?�&�?�=q?�33?���?�b?���?�;d?�x�?�K�?���?�{?���?���?�`B?�o?�C�?�n�?|j?y�?ix�?_|�?]/?T9X?MO�?H1'?C�
?@�?:�H?333?0bN?*=q?%�?�w?�?�?�y>���>�X>��j>��>�1>���>�>���>�z�>�n�>�I�>Õ�>�^5>���>�?}>�>�"�>���>�$�>|�>p��>k�>N�>9X>$�/>�P>	7L=�=��`=ě�=�{=�7L=}�=u=m�h=0 �<�/<�t�;��
�ě���j��h�t��P�`�]/�}󶽓t����T��9X��Q����ě�����;d��`B��xս��   ���1'�I��bN�t���+���#�
�-V�1&�6E��<j�D���I�^�J���Kƨ�R�W
=�["Ѿ]/�bMӾe`B�m�h�r�!�w�پz�H�~�۾�  ���7������$ݾ��9��7L������O߾��;���Ͼ��+�������㾜(����-��5?���w��G���S����
���
��Z��Z���T��r����義�羪~������V��{������&龲�!���j��E���X�������=q��V���`��t�����
=�������ۥ��/��Ĝ��MӾ�S���`B���y��r������1������!��9X��KǾ�����^5��j��p��������ۿ Ĝ� Ĝ�G��S��������`B��T����y���r���ÿ
=q�
����C����1��ͿO߿{����������;� ſ�׿����t���Ͽ9X�z��������+�
=�
=��ٿQ���������^5����"ѿdZ����m�(����p��5?��R�|�   � A�� ��!%�!�7�!�7�!���!���"��#S��#�
�$��$Z�$���%��%�˿%�T�&��'��(1'�(r��(�ÿ)�^�)��)��*=q�*=q�*~��*���+C��+ƨ�,1�,�Ϳ,�Ϳ-��.{�.{�.���.��/��.��.��/��/��/\)�0 ſ0�׿0�`�1hs�1녿1녿2-�2n��2-�2�!�2�333�333�3t��3�Ͽ49X�4z�4�j�4���5?}�5?}�5��5�6E��6�+�6ȴ�7Kǿ7�P�7�ٿ8b�8b�8Q�8Q�8�u�8���8���9��9��9X�9X�9X�9X�9X�9X�9���9�#�:��9�#�9�#�9�#�:��:^5�:^5�:�H�;"ѿ;"ѿ;"ѿ;"�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�`B@�p�@�x�@��@��7@�p�@�`B@�O�@�%@��j@��/@�%@��@���@�  @��@��;@��
@���@��F@��@���@���@���@��@�dZ@�\)@�@���@��+@�v�@�E�@�J@�@��h@���@���@��h@���@�p�@�hs@��/@���@��u@�j@� �@���@�ƨ@�dZ@��@��@���@�ȴ@���@���@�n�@�^5@�^5@�^5@�^5@��#@�x�@�X@��`@���@��D@�z�@�1@��;@��F@��P@�K�@�K�@��@�o@��@��@��@�K�@�33@���@��!@�M�@�M�@�-@�$�@�$�@�=q@�M�@�ff@�^5@�ff@�V@�M�@�M�@�@���@��@�p�@��/@��/@�9X@���@�l�@�;d@���@�=q@��@��-@��7@��@���@���@��@��m@�C�@���@���@�~�@�=q@��@��h@�G�@�&�@�V@���@�r�@�  @�;d@��!@�$�@��@�@���@��@�?}@���@���@�z�@�1@}�@}�h@|��@|9X@{dZ@z��@yx�@yG�@y�@y�@x�9@x  @w��@w+@vV@u��@t�D@r��@p��@m�-@kƨ@j�@h��@g\)@g
=@f�@f�@ep�@a%@]�h@\�/@\j@W�P@SC�@PQ�@N5?@J�@IG�@G�w@FV@D��@@�9@<�@;�
@9��@81'@5�T@5V@3o@2�\@0�9@-p�@*�!@'��@%?}@!�#@;d@�@��@�@�m@b@�@
M�@ �@�D?��H?�!?�D?�K�?��?���?�&�?�=q?�33?���?�b?���?�;d?�x�?�K�?���?�{?���?���?�`B?�o?�C�?�n�?|j?y�?ix�?_|�?]/?T9X?MO�?H1'?C�
?@�?:�H?333?0bN?*=q?%�?�w?�?�?�y>���>�X>��j>��>�1>���>�>���>�z�>�n�>�I�>Õ�>�^5>���>�?}>�>�"�>���>�$�>|�>p��>k�>N�>9X>$�/>�P>	7L=�=��`=ě�=�{=�7L=}�=u=m�h=0 �<�/<�t�;��
�ě���j��h�t��P�`�]/�}󶽓t����T��9X��Q����ě�����;d��`B��xս��   ���1'�I��bN�t���+���#�
�-V�1&�6E��<j�D���I�^�J���Kƨ�R�W
=�["Ѿ]/�bMӾe`B�m�h�r�!�w�پz�H�~�۾�  ���7������$ݾ��9��7L������O߾��;���Ͼ��+�������㾜(����-��5?���w��G���S����
���
��Z��Z���T��r����義�羪~������V��{������&龲�!���j��E���X�������=q��V���`��t�����
=�������ۥ��/��Ĝ��MӾ�S���`B���y��r������1������!��9X��KǾ�����^5��j��p��������ۿ Ĝ� Ĝ�G��S��������`B��T����y���r���ÿ
=q�
����C����1��ͿO߿{����������;� ſ�׿����t���Ͽ9X�z��������+�
=�
=��ٿQ���������^5����"ѿdZ����m�(����p��5?��R�|�   � A�� ��!%�!�7�!�7�!���!���"��#S��#�
�$��$Z�$���%��%�˿%�T�&��'��(1'�(r��(�ÿ)�^�)��)��*=q�*=q�*~��*���+C��+ƨ�,1�,�Ϳ,�Ϳ-��.{�.{�.���.��/��.��.��/��/��/\)�0 ſ0�׿0�`�1hs�1녿1녿2-�2n��2-�2�!�2�333�333�3t��3�Ͽ49X�4z�4�j�4���5?}�5?}�5��5�6E��6�+�6ȴ�7Kǿ7�P�7�ٿ8b�8b�8Q�8Q�8�u�8���8���9��9��9X�9X�9X�9X�9X�9X�9���9�#�:��9�#�9�#�9�#�:��:^5�:^5�:�H�;"ѿ;"ѿ;"ѿ;"�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBw�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bw�Bu�Bu�Bu�Bt�Bt�Bu�Bu�Bu�Bu�Bu�Bu�Bt�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bs�Bt�Bt�Bt�Bs�Bs�Br�Bt�Br�Bt�Bs�Bs�Bs�Bs�Bq�Br�Bs�Br�Br�Br�Br�Br�Br�Br�Bq�Br�Bq�Bp�Bq�Bq�Bo�Bp�Bo�Bo�Bo�Bp�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bo�Bm�Bm�Bl�Bk�Bk�BjBjBjBjBjBk�Bk�Bk�Bk�Bk�Bk�BjBjBiyBhsBhsBiyBiyBhsBhsBhsBgmBgmBgmBgmBgmBdZBffBe`BdZBdZBdZBcTBcTBdZBcTBcTBbNB`BB_;B^5B]/B]/B\)B\)BZB[#B\)B]/B]/B\)B[#BZBZBXBVBS�BQ�BN�BL�BK�BI�BG�BG�BF�BE�BB�B>wB<jB:^B8RB5?B2-B/B-B+B(�B'�B%�B#�B�B�B�B�B�B�B�B�B�B�B�BuBbBVBPBDB	7B+BBB��B��B��B��B��B�B�B�B�B�mB�`B�NB�BB�/B�)B�B�B�B�
B�B��B��B��B��B��B��B��BƨBĜB��B��B��B�}B�wB�qB�jB�jB�dB�dB�^B�XB�XB�RB�FB�?B�?B�9B�9B�3B�3B�3B�-B�-B�-B�-B�-B�'B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 Bw�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bw�Bu�Bu�Bu�Bt�Bt�Bu�Bu�Bu�Bu�Bu�Bu�Bt�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bs�Bt�Bt�Bt�Bs�Bs�Br�Bt�Br�Bt�Bs�Bs�Bs�Bs�Bq�Br�Bs�Br�Br�Br�Br�Br�Br�Br�Bq�Br�Bq�Bp�Bq�Bq�Bo�Bp�Bo�Bo�Bo�Bp�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bo�Bm�Bm�Bl�Bk�Bk�BjBjBjBjBjBk�Bk�Bk�Bk�Bk�Bk�BjBjBiyBhsBhsBiyBiyBhsBhsBhsBgmBgmBgmBgmBgmBdZBffBe`BdZBdZBdZBcTBcTBdZBcTBcTBbNB`BB_;B^5B]/B]/B\)B\)BZB[#B\)B]/B]/B\)B[#BZBZBXBVBS�BQ�BN�BL�BK�BI�BG�BG�BF�BE�BB�B>wB<jB:^B8RB5?B2-B/B-B+B(�B'�B%�B#�B�B�B�B�B�B�B�B�B�B�B�BuBbBVBPBDB	7B+BBB��B��B��B��B��B�B�B�B�B�mB�`B�NB�BB�/B�)B�B�B�B�
B�B��B��B��B��B��B��B��BƨBĜB��B��B��B�}B�wB�qB�jB�jB�dB�dB�^B�XB�XB�RB�FB�?B�?B�9B�9B�3B�3B�3B�-B�-B�-B�-B�-B�'B�'B�'B�'B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040282018111910402920181119104028  IF  ARFMCODA024c                                                                20181026190020                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190046  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181026190046  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104029  IP  PSAL            @   D��3G�O�                