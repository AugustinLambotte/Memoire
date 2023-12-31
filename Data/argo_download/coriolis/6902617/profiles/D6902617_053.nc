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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181026190020  20181119104029  6902617 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               5A   IF                                  2C  D   NOVA                            SN187                           n/a                             865 @���8p�k1   @����7�;@S�c.��@cp:$� 1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A���A���A���A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BX  B\  B_��Bd  Bh  Bl  BpffBt  Bx  B|  B�  B�33B�33B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B���B�  B�33B�33B�  B�  B�  B�  B�  C  C� C  C	��C  C� C  C� C  C� C  C� C �C"� C%  C'� C*  C,� C/  C1��C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CV�fCY� C\  C^ffC`�fCcffCe�fChffCk  Cm��Cp  Cr� Cu  Cw� Cz  C|� C  C��3C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C��C�L�C�� C��3C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��3C�@ Cό�C�� C�  C�@ CԀ Cճ3C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D ��D��D@ D� D� DfD@ D	� D
�fDfDFfD� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!FfD"�fD#�fD%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7��D9  D:FfD;� D<��D=��D?9�D@� DA��DC  DD@ DE� DF� DH  DI@ DJ�fDK� DM  DN@ DO� DP� DQ��DS@ DT� DU� DWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�` D���D���D�@ D�� D�� D�  D���D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�\�D�  DȠ D�@ D�� Dʀ D��D�� D�c3D�  D͠ D�@ D��3Dπ D�  Dм�D�` D���DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�c3D�  D� D�C3D��3D� D�  D�� D�c3D�  D��D�@ D��3D� D�  D�� D�` D�  D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D��3D�C3D��3D�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A���A���A���A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BX  B\  B_��Bd  Bh  Bl  BpffBt  Bx  B|  B�  B�33B�33B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B���B�  B�33B�33B�  B�  B�  B�  B�  C  C� C  C	��C  C� C  C� C  C� C  C� C �C"� C%  C'� C*  C,� C/  C1��C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CV�fCY� C\  C^ffC`�fCcffCe�fChffCk  Cm��Cp  Cr� Cu  Cw� Cz  C|� C  C��3C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C��C�L�C�� C��3C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��3C�@ Cό�C�� C�  C�@ CԀ Cճ3C�  C�@ Cـ C�� C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C���C�  C�� C�  D ��D��D@ D� D� DfD@ D	� D
�fDfDFfD� D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!FfD"�fD#�fD%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7��D9  D:FfD;� D<��D=��D?9�D@� DA��DC  DD@ DE� DF� DH  DI@ DJ�fDK� DM  DN@ DO� DP� DQ��DS@ DT� DU� DWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�` D���D���D�@ D�� D�� D�  D���D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D�� D�\�D���D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�\�D�  DȠ D�@ D�� Dʀ D��D�� D�c3D�  D͠ D�@ D��3Dπ D�  Dм�D�` D���DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�c3D�  D� D�C3D��3D� D�  D�� D�c3D�  D��D�@ D��3D� D�  D�� D�` D�  D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D��3D�C3D��3D�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@d1@d1@c��@d�@d(�@c��@c�@c�F@co@`b@^V@\j@ZM�@Y�7@XA�@W�w@W�@V$�@U@U�@Vff@VV@U�T@U�@U�h@U�h@U��@U��@U�h@Up�@Up�@Up�@U/@U`B@U�@U�@U�@U�@U`B@U`B@U`B@U`B@U��@U�-@Up�@T�@S��@S�F@So@Q�@Q�7@Q7L@P�`@PQ�@PQ�@P1'@P  @O\)@M��@M?}@L�@L��@J^5@G�@E��@DZ@C�F@B�@B��@?�;@=�@=O�@<j@;�
@;33@:�@:��@9��@7�@5�T@5��@6�+@7
=@7;d@7�@7
=@6�y@7�@8  @9G�@:~�@:n�@:�\@:-@7K�@3�@3o@0��@0r�@.�y@-�@+dZ@*^5@)�^@(bN@&�y@&{@%��@$�@$�@$�@$j@$�@"J@�w@��@?}@n�@E�@O�@�D@-@�w@�@t�@
-@	X@��@j@t�@n�@�@�7@ ��?�;d?���?��-?�(�?�"�?���?�=q?�^5?�=q?���?�Q�?��u?���?��9?�b?�?}?�!?�J?�Z?�z�?��
?�M�?�M�?�J?�-?�7?�G�?�Ĝ?��;?�A�?�A�?�p�?���?���?�?�(�?��H?�^?�9?�1'?�b?���?�?�P?�K�?�ȴ?��?�E�?�$�?��/?�9X?��`?�v�?�V?���?�V?���?��?�J?�;d?�?��?��?��^?�+?��j?�7L?�?p�`?e��?^�R?Z�?O��?L1?G+?@Ĝ?=p�?<(�?9��?0 �?%�T?�!?��>�^5>�K�>�9X>>��y>�Z>���>�(�>�t�>�O�>ɺ^>Ǯ>Õ�>�dZ>�->���>���>��T>���>�I�>��^>���>��\>~��>x��>w��>s�F>l�D>_;d>^5?>aG�>cS�>Xb>F��>9X>6E�>2->.{>(��>!��>�w>�->z�>+==�x�=��=�j=��-=�O�=u=8Q�<��<�1��o��o���
�t��ě����D���y�#��hs���-��-����Ƨ����"ѽ��m�bN��P��-�&�y�,1�6E��;dZ�?|�B�\�E�˾M��P�`�S�ϾY��]/�bMӾe`B�fff�gl��j~��k��m�h�r�!�t�j�vȴ�z�H�}�~�۾~�۾�  ���7��o������𾈴9���^��C����;����bN��hs��n���zᾕ����P��b����������-��5?��;d��A���Ĝ��MӾ�Z���
���
��`B���þ�~���V���h���D���h���h����������&龶E����#��푾�vɾ��۾��۾�|��o��+��C���녾Ձ�և+���և+��b��(��޸R��Ĝ��G���MӾ��/��xվ������V���!��9X���پ�^5��j���ۿ Ĝ�J���S�����/����˿$ݿ$ݿ$ݿ$ݿff��y���r���9�	��
���1�O߿����bN����녿n���33��F��Ͽ9X�������E���+�
=�Kǿ�ٿb�Q��u����X��#�^5�"ѿ��(����/�p���-�vɿ|� Ĝ�"J�#S��#�
�$��%��%�T�&ff�'+�'(�ÿ)7L�)xտ*~��+C��+��+��,1�,I��,�D�,�Ϳ-V�-O߿-O߿-O߿.V�/\)�/�;�/�;�0 ſ0bN�0�`�1���1녿2-�2�!�2�3t��3�Ͽ49X�4�j�4�j�4�j�5?}�5�6�6�+�6ȴ�7Kǿ7�P�7�P�7�ٿ8b�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�9��9X�9���9���9���9���9�#�9�#�9�#�9�#�:��:���;"ѿ;dZ�;"ѿ;"ѿ;��;��;��;��;dZ�;dZ�;��;dZ�;dZ�;"ѿ;"ѿ;"ѿ;"ѿ;"ѿ:���:^5�:���:�H�:�H�;"ѿ:�H�;"ѿ;"ѿ;dZ�;��;�m�<(��<j�<푿=�-�=�-�=p��=p��=p��=p��=�-�=�-�=�-�=�-�=p�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @d1@d1@c��@d�@d(�@c��@c�@c�F@co@`b@^V@\j@ZM�@Y�7@XA�@W�w@W�@V$�@U@U�@Vff@VV@U�T@U�@U�h@U�h@U��@U��@U�h@Up�@Up�@Up�@U/@U`B@U�@U�@U�@U�@U`B@U`B@U`B@U`B@U��@U�-@Up�@T�@S��@S�F@So@Q�@Q�7@Q7L@P�`@PQ�@PQ�@P1'@P  @O\)@M��@M?}@L�@L��@J^5@G�@E��@DZ@C�F@B�@B��@?�;@=�@=O�@<j@;�
@;33@:�@:��@9��@7�@5�T@5��@6�+@7
=@7;d@7�@7
=@6�y@7�@8  @9G�@:~�@:n�@:�\@:-@7K�@3�@3o@0��@0r�@.�y@-�@+dZ@*^5@)�^@(bN@&�y@&{@%��@$�@$�@$�@$j@$�@"J@�w@��@?}@n�@E�@O�@�D@-@�w@�@t�@
-@	X@��@j@t�@n�@�@�7@ ��?�;d?���?��-?�(�?�"�?���?�=q?�^5?�=q?���?�Q�?��u?���?��9?�b?�?}?�!?�J?�Z?�z�?��
?�M�?�M�?�J?�-?�7?�G�?�Ĝ?��;?�A�?�A�?�p�?���?���?�?�(�?��H?�^?�9?�1'?�b?���?�?�P?�K�?�ȴ?��?�E�?�$�?��/?�9X?��`?�v�?�V?���?�V?���?��?�J?�;d?�?��?��?��^?�+?��j?�7L?�?p�`?e��?^�R?Z�?O��?L1?G+?@Ĝ?=p�?<(�?9��?0 �?%�T?�!?��>�^5>�K�>�9X>>��y>�Z>���>�(�>�t�>�O�>ɺ^>Ǯ>Õ�>�dZ>�->���>���>��T>���>�I�>��^>���>��\>~��>x��>w��>s�F>l�D>_;d>^5?>aG�>cS�>Xb>F��>9X>6E�>2->.{>(��>!��>�w>�->z�>+==�x�=��=�j=��-=�O�=u=8Q�<��<�1��o��o���
�t��ě����D���y�#��hs���-��-����Ƨ����"ѽ��m�bN��P��-�&�y�,1�6E��;dZ�?|�B�\�E�˾M��P�`�S�ϾY��]/�bMӾe`B�fff�gl��j~��k��m�h�r�!�t�j�vȴ�z�H�}�~�۾~�۾�  ���7��o������𾈴9���^��C����;����bN��hs��n���zᾕ����P��b����������-��5?��;d��A���Ĝ��MӾ�Z���
���
��`B���þ�~���V���h���D���h���h����������&龶E����#��푾�vɾ��۾��۾�|��o��+��C���녾Ձ�և+���և+��b��(��޸R��Ĝ��G���MӾ��/��xվ������V���!��9X���پ�^5��j���ۿ Ĝ�J���S�����/����˿$ݿ$ݿ$ݿ$ݿff��y���r���9�	��
���1�O߿����bN����녿n���33��F��Ͽ9X�������E���+�
=�Kǿ�ٿb�Q��u����X��#�^5�"ѿ��(����/�p���-�vɿ|� Ĝ�"J�#S��#�
�$��%��%�T�&ff�'+�'(�ÿ)7L�)xտ*~��+C��+��+��,1�,I��,�D�,�Ϳ-V�-O߿-O߿-O߿.V�/\)�/�;�/�;�0 ſ0bN�0�`�1���1녿2-�2�!�2�3t��3�Ͽ49X�4�j�4�j�4�j�5?}�5�6�6�+�6ȴ�7Kǿ7�P�7�P�7�ٿ8b�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�8�u�9��9X�9���9���9���9���9�#�9�#�9�#�9�#�:��:���;"ѿ;dZ�;"ѿ;"ѿ;��;��;��;��;dZ�;dZ�;��;dZ�;dZ�;"ѿ;"ѿ;"ѿ;"ѿ;"ѿ:���:^5�:���:�H�:�H�;"ѿ:�H�;"ѿ;"ѿ;dZ�;��;�m�<(��<j�<푿=�-�=�-�=p��=p��=p��=p��=�-�=�-�=�-�=�-�=p�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB\)B\)B[#B[#B[#BZBZBXBXBVBT�BVBT�BQ�BR�BP�BQ�BQ�BP�BP�BP�BQ�BP�BP�BP�BP�BP�BP�BP�BP�BP�BP�BP�BO�BP�BP�BP�BP�BP�BP�BO�BO�BO�BO�BN�BN�BL�BK�BJ�BJ�BI�BH�BG�BG�BF�BE�BD�BD�BB�B@�B?}B=qB;dB8RB6FB33B2-B1'B1'B-B+B+B)�B(�B'�B&�B&�B$�B%�B$�B#�B$�B%�B&�B'�B'�B(�B)�B,B-B0!B/B-B+B,B%�B$�B$�B"�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BuBhB\BPBJB
=B1B+BBBBBBB  B  B  B  B  B  B  B  B  B  B  B  B  B  B  B  BBB  B  B��B  BBBBBBBBBBBBB  BBB  B  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�mB�TB�/B�)B�B�B�B��B��B��B��B��BǮBŢBÖBÖB��B��B��B�wB�wB�qB�jB�dB�XB�LB�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�3B�-B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�'B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�'B�!B�!B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B\)B\)B[#B[#B[#BZBZBXBXBVBT�BVBT�BQ�BR�BP�BQ�BQ�BP�BP�BP�BQ�BP�BP�BP�BP�BP�BP�BP�BP�BP�BP�BP�BO�BP�BP�BP�BP�BP�BP�BO�BO�BO�BO�BN�BN�BL�BK�BJ�BJ�BI�BH�BG�BG�BF�BE�BD�BD�BB�B@�B?}B=qB;dB8RB6FB33B2-B1'B1'B-B+B+B)�B(�B'�B&�B&�B$�B%�B$�B#�B$�B%�B&�B'�B'�B(�B)�B,B-B0!B/B-B+B,B%�B$�B$�B"�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BuBhB\BPBJB
=B1B+BBBBBBB  B  B  B  B  B  B  B  B  B  B  B  B  B  B  B  BBB  B  B��B  BBBBBBBBBBBBB  BBB  B  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�mB�TB�/B�)B�B�B�B��B��B��B��B��BǮBŢBÖBÖB��B��B��B�wB�wB�qB�jB�dB�XB�LB�FB�FB�FB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�3B�-B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�'B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�'B�!B�!B�!B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811191040292018111910402920181119104029  IF  ARFMCODA024c                                                                20181026190020                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181026190046  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181026190046  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181119104029  IP  PSAL            @ffD�33G�O�                