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
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090451  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               %A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�VQ ��>1   @�VQ��6�@O��"���?����>2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BW��B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�33B�  B���B�  B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�33B�33B�  B���B�  B�  B�  B�33B�  B�  C  C� C  C	� C  C��C�C� C  C� C  C��C �C"� C%  C'��C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\�C^� C`�fCcffCf  Ch��Ck�Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�L�C���C�� C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�s3C�� C��C�L�C�� C��3C��3C�@ C���C���C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��C�L�Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�� C��3C�@ C�� C�� C�  C���C�  D � D  D@ D� D� D  DFfD	� D
� D  D@ D� D� D  D@ D� D� D  D9�D� D� D  D@ D� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-�fD/fD0FfD1� D2� D4  D5FfD6�fD7�fD9fD:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK��DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^�fD_� Da  DbFfDc�fDd� De��Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|�fD}� D~��D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D���D�|�D��D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�@ D���D�� D�#3D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�c3D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���Dà D�C3D�� Dŀ D�  D�� D�\�D�  DȠ D�@ D�� Dʀ D�  D˼�D�` D�  D͠ D�C3D��3Dσ3D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܣ3D�@ D���D�|�D�  D�� D�c3D�3D� D�@ D��3D� D�  D�� D�` D�  D��D�@ D�� D� D�  D��3D�c3D�3D�3D�C3D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D�� D�` D�  D�� D�<�D�� D��3D�#3D��3D�` D�  D�� D�I�D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BW��B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�33B�  B���B�  B�  B�  B�  B�33B�  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�33B�33B�  B���B�  B�  B�  B�33B�  B�  C  C� C  C	� C  C��C�C� C  C� C  C��C �C"� C%  C'��C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\�C^� C`�fCcffCf  Ch��Ck�Cm� Cp  Cr� Cu  Cw��Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�L�C���C�� C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C�s3C�� C��C�L�C�� C��3C��3C�@ C���C���C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C��C�L�Cπ Cг3C�  C�@ CԀ C�� C�  C�@ Cـ C���C�  C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�� C��3C�@ C�� C�� C�  C���C�  D � D  D@ D� D� D  DFfD	� D
� D  D@ D� D� D  D@ D� D� D  D9�D� D� D  D@ D� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-�fD/fD0FfD1� D2� D4  D5FfD6�fD7�fD9fD:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK��DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^�fD_� Da  DbFfDc�fDd� De��Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|�fD}� D~��D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D���D�|�D��D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�3D�� D�@ D���D�� D�#3D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�c3D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���Dà D�C3D�� Dŀ D�  D�� D�\�D�  DȠ D�@ D�� Dʀ D�  D˼�D�` D�  D͠ D�C3D��3Dσ3D�  Dм�D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܣ3D�@ D���D�|�D�  D�� D�c3D�3D� D�@ D��3D� D�  D�� D�` D�  D��D�@ D�� D� D�  D��3D�c3D�3D�3D�C3D�� D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D�� D�` D�  D�� D�<�D�� D��3D�#3D��3D�` D�  D�� D�I�D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��u@���@�z�@�j@�I�@�Q�@�r�@�z�@�1'@�Z@�1'@�S�@���@�V@���@���@�7L@�V@�V@�V@��@���@���@�z�@�Q�@��;@�;d@���@�M�@�{@��@��-@���@���@��7@�x�@�p�@�p�@�p�@�p�@�hs@�hs@�p�@�`B@�X@�O�@�G�@�?}@�/@�&�@�/@�/@��@��@���@���@��@��@��`@��`@��/@��/@���@��/@���@���@���@��j@��9@��9@���@��u@��u@��D@��@�r�@�r�@�r�@�r�@�bN@�I�@�9X@�1'@�1'@�(�@�(�@� �@�b@�  @���@��;@��
@��
@�ƨ@�ƨ@��w@��F@��F@��@��@��@���@�|�@��@�|�@�|�@�\)@�K�@�;d@�o@�o@�
=@��!@�M�@���@��D@�?}@��@���@���@�O�@���@�(�@��+@���@��@�ȴ@���@��F@�S�@���@��#@���@���@���@�  @���@��@�E�@��#@���@�7L@��j@�I�@�1@��;@��w@���@���@�33@���@�~�@�M�@�$�@�J@��@���@��@��/@��j@��u@�Z@� �@�  @��w@�\)@���@���@���@�^5@�V@�5?@�J@��^@��@�O�@�V@��@���@��D@�A�@�1@�ƨ@�l�@�S�@�K�@�C�@�;d@�+@��@��\@�M�@�{@�J@�@��7@�7L@�/@��@��@�%@���@��u@�Z@�1'@� �@���@��;@��;@���@�|�@�S�@�33@�+@�"�@�"�@�o@�@���@��@��y@��H@��@��@�ȴ@���@���@���@��!@��+@��+@��+@�v�@�n�@�5?@�=q@�5?@�J@�@�@�J@��@��@��T@���@���@��7@��@�G�@�V@�V@�Ĝ@���@�Z@�b@��@��
@��
@��@��@�t�@�C�@�
=@���@�=q@�$�@�@��#@���@�p�@�?}@�7L@��@��@���@��u@��@��@��@��@�r�@�9X@��@���@�dZ@�S�@�K�@�33@��@�
=@��y@���@�~�@�V@��^@�x�@��`@��/@��9@�A�@�9X@��@��;@��@�S�@��@��R@��R@��\@�=q@�J@��^@��@��-@���@��T@��#@��^@���@��7@�X@�?}@�V@���@�j@�(�@�@|�@l�@+@~�y@~��@~E�@}�T@}��@}��@}O�@|��@|��@|�D@|�@{��@{�
@{��@{�@{dZ@{dZ@{33@z�@z��@z=q@z�@y�@y��@y7L@x��@x��@x�@xbN@x �@x  @w|�@w�@vȴ@vȴ@v��@v5?@v$�@vV@vv�@vV@vV@vff@v�R@v�R@vȴ@v��@vE�@v$�@v5?@v5?@v5?@vE�@u�@u��@vV@v��@vv�@u�T@u�h@u�@uO�@uO�@uO�@u?}@u`B@uO�@u�@up�@uO�@uV@uV@tZ@s��@t�@uO�@u�@u��@u��@up�@u?}@v��@vȴ@vV@vff@v��@vE�@vV@v$�@v{@v�+@x  @xbN@y�@y�@y&�@y��@zM�@z�H@z��@z�@yX@y��@zJ@y�^@y&�@x1'@w�@w\)@w+@w�@wK�@x �@y7L@x�`@x��@xĜ@x�@xbN@xQ�@xb@w�w@w|�@wl�@vȴ@v�+@vV@vv�@vff@vV@vE�@v@u�T@u�h@uO�@t��@t��@t�/@tj@t1@s�m@s��@s��@sS�@so@r�H@r��@r=q@q��@q�^@qhs@q�@p��@pĜ@p�@p1'@p1'@p1'@p1'@o��@o��@oK�@o
=@n��@nV@n@m��@m�-@mO�@m�@l��@l�@lj@l�@k��@kC�@k@j��@j�!@jn�@j-@j�@i��@i��@i�7@iG�@h��@h�@hb@g�w@g\)@g
=@f��@f�+@f�+@f�+@f�+@fv�@f5?@e�@e�@d��@d��@c�F@c�@cC�@b�H@b-@a��@a��@a��@a��@ax�@ahs@a%@`��@`�9@`1'@`b@_�@_�P@_l�@^��@^�y@^�@^��@^��@^ȴ@^��@^ff@]�-@]�-1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��u@���@�z�@�j@�I�@�Q�@�r�@�z�@�1'@�Z@�1'@�S�@���@�V@���@���@�7L@�V@�V@�V@��@���@���@�z�@�Q�@��;@�;d@���@�M�@�{@��@��-@���@���@��7@�x�@�p�@�p�@�p�@�p�@�hs@�hs@�p�@�`B@�X@�O�@�G�@�?}@�/@�&�@�/@�/@��@��@���@���@��@��@��`@��`@��/@��/@���@��/@���@���@���@��j@��9@��9@���@��u@��u@��D@��@�r�@�r�@�r�@�r�@�bN@�I�@�9X@�1'@�1'@�(�@�(�@� �@�b@�  @���@��;@��
@��
@�ƨ@�ƨ@��w@��F@��F@��@��@��@���@�|�@��@�|�@�|�@�\)@�K�@�;d@�o@�o@�
=@��!@�M�@���@��D@�?}@��@���@���@�O�@���@�(�@��+@���@��@�ȴ@���@��F@�S�@���@��#@���@���@���@�  @���@��@�E�@��#@���@�7L@��j@�I�@�1@��;@��w@���@���@�33@���@�~�@�M�@�$�@�J@��@���@��@��/@��j@��u@�Z@� �@�  @��w@�\)@���@���@���@�^5@�V@�5?@�J@��^@��@�O�@�V@��@���@��D@�A�@�1@�ƨ@�l�@�S�@�K�@�C�@�;d@�+@��@��\@�M�@�{@�J@�@��7@�7L@�/@��@��@�%@���@��u@�Z@�1'@� �@���@��;@��;@���@�|�@�S�@�33@�+@�"�@�"�@�o@�@���@��@��y@��H@��@��@�ȴ@���@���@���@��!@��+@��+@��+@�v�@�n�@�5?@�=q@�5?@�J@�@�@�J@��@��@��T@���@���@��7@��@�G�@�V@�V@�Ĝ@���@�Z@�b@��@��
@��
@��@��@�t�@�C�@�
=@���@�=q@�$�@�@��#@���@�p�@�?}@�7L@��@��@���@��u@��@��@��@��@�r�@�9X@��@���@�dZ@�S�@�K�@�33@��@�
=@��y@���@�~�@�V@��^@�x�@��`@��/@��9@�A�@�9X@��@��;@��@�S�@��@��R@��R@��\@�=q@�J@��^@��@��-@���@��T@��#@��^@���@��7@�X@�?}@�V@���@�j@�(�@�@|�@l�@+@~�y@~��@~E�@}�T@}��@}��@}O�@|��@|��@|�D@|�@{��@{�
@{��@{�@{dZ@{dZ@{33@z�@z��@z=q@z�@y�@y��@y7L@x��@x��@x�@xbN@x �@x  @w|�@w�@vȴ@vȴ@v��@v5?@v$�@vV@vv�@vV@vV@vff@v�R@v�R@vȴ@v��@vE�@v$�@v5?@v5?@v5?@vE�@u�@u��@vV@v��@vv�@u�T@u�h@u�@uO�@uO�@uO�@u?}@u`B@uO�@u�@up�@uO�@uV@uV@tZ@s��@t�@uO�@u�@u��@u��@up�@u?}@v��@vȴ@vV@vff@v��@vE�@vV@v$�@v{@v�+@x  @xbN@y�@y�@y&�@y��@zM�@z�H@z��@z�@yX@y��@zJ@y�^@y&�@x1'@w�@w\)@w+@w�@wK�@x �@y7L@x�`@x��@xĜ@x�@xbN@xQ�@xb@w�w@w|�@wl�@vȴ@v�+@vV@vv�@vff@vV@vE�@v@u�T@u�h@uO�@t��@t��@t�/@tj@t1@s�m@s��@s��@sS�@so@r�H@r��@r=q@q��@q�^@qhs@q�@p��@pĜ@p�@p1'@p1'@p1'@p1'@o��@o��@oK�@o
=@n��@nV@n@m��@m�-@mO�@m�@l��@l�@lj@l�@k��@kC�@k@j��@j�!@jn�@j-@j�@i��@i��@i�7@iG�@h��@h�@hb@g�w@g\)@g
=@f��@f�+@f�+@f�+@f�+@fv�@f5?@e�@e�@d��@d��@c�F@c�@cC�@b�H@b-@a��@a��@a��@a��@ax�@ahs@a%@`��@`�9@`1'@`b@_�@_�P@_l�@^��@^�y@^�@^��@^��@^ȴ@^��@^ff@]�-@]�-1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B!�B!�B �B"�B!�B �B�B�B�B�BhBbB\BVBPBPBVBDB
=B1B+BBBBBBBB  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�mB�mB�fB�fB�fB�fB�`B�`B�ZB�TB�NB�NB�HB�HB�BB�;B�;B�;B�;B�5B�5B�/B�/B�)B�)B�#B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�
B�
B�
B�
B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBȴBȴBƨBǮBǮBƨBƨBƨBƨBƨBƨBŢBŢBÖBÖBBB��B��B��B��B��B�}B�wB�qB�jB�^B�XB�RB�RB�LB�LB�LB�FB�?B�9B�3B�-B�-B�-B�'B�!B�B�B�!B�'B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�3B�9B�3B�3B�9B�9B�9B�?B�?B�LB�dB�qB�}B�}B��B��BĜBƨBƨBƨBŢBǮBȴBǮBƨBŢBĜBÖBÖBÖBĜBǮB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBȴBȴBȴBȴBǮBǮBǮBǮBȴBȴBȴBȴBɺBɺBȴBǮBǮBŢBŢBŢBŢBĜBĜBĜBĜBĜBĜBĜBĜBÖBĜBÖBÖBÖBÖBÖBĜBĜBĜBĜBŢBŢBŢBŢBŢBŢ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B}BwB}B}B}B}BwBwBqBqBkBwBqBqBqBlBqBlBlBqBqBlBrBlBlBrBxB}BxBxBxB~B~B~B~B~B~B~B~B~B~B~B�B�B�B�B�B�B�B�B�B�B~B~B~B~B~B~B~B~BxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBxBrBxBxBrBrBrBrBrBrBrBrBrBrBxBxBxBxBxBxBxBrBrBrBrBrBrBrBrBrB~B�B�B�B~B�B�B�B�B�B�B�BsBgBgBNB	6B0B*B$BBB$BBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�zB�tB�tB�oB�oB�oB�iB�iB�bB�bB�\B�VB�VB�PB�JB�DB�>B�>B�7B�7B�7B�7B�1B�1B�+B�&B� B� B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɿBɿBɿBȸBǲBƬBƬBƬBƬBƬBŦBŦBĠBÛBB��B��B��B��B��B��B��B�|B��B��B�|B�|B�|B�|B�|B�|B�vB�vB�jB�jB�cB�cB�]B�]B�]B�XB�XB�RB�LB�FB�?B�3B�-B�'B�'B�!B�!B�!B�B�B�B�	B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�	B�B�	B�	B�B�B�B�B�B�"B�:B�GB�RB�RB�XB�^B�qB�}B�}B�}B�wB��B��B��B�}B�wB�qB�kB�kB�kB�qB��BÜBÜBÜBĢBĢBĢBĢBÜBĢBĢBĢBÜBÜBĢBŨBƮBƮBƮBǳBǳBǳBǳBǳBǳBǴBƮBƮBƮBƮBƮBƮBŨBŨBŨBŨBŨBŨBĢBĢBĢBĢBĢBĢBĢBĢBĢBĢBĢBĢBĢBÜBÜBÜBÜBÜBÜBÜBÜBÜBBBBBBBBBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�xB�xB�xB�xB�rB�rB�rB�rB�rB�rB�rB�rB�lB�rB�lB�lB�lB�lB�lB�rB�rB�rB�rB�xB�xB�xB�xB�xB�x1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0079867                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904512018110709045120181107090451  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172718  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172718  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090451  IP  PSAL            @ffD�� G�O�                