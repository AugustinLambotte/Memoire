CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  ,   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
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
_FillValue                 ,  BL   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  Dx   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 ,  M(   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  OT   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  X   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 ,  `�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  b�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 ,  k�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  m�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  vl   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 ,     PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �H   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 ,  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �$   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �0   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �4   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �8   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �<   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �@   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090451  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               &A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�X��W1   @�XѢCe@P	��Z��@���DW�2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   AffA   A0  A>ffAP  A`  Ap  A�  A���A���A���A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A�33B   B��B  B��B��B  B  B  B   B$  B(  B,  B0  B3��B8  B<  B@  BD  BH  BL  BP  BT  BXffB\ffB`  Bd  Bh  BlffBp  Bt  Bx  B|ffB�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  C�fC� C�C	��C�C��C  C� C  C� C�C��C   C"� C%�C'��C*�C,��C/�C1��C4�C6� C8�fC;ffC>  C@� CC  CE� CH  CJ� CM  CO� CR  CT��CW  CY� C\  C^� Ca�Cc� Cf  Ch� Ck  Cm� Cp�Cr� Cu  Cw� Cy�fC|� C  C�� C�  C�@ C�� C�� C��3C�33C�s3C��3C��3C�33C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C���C��C�L�C�� C��3C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C��3C�  C�L�Cŀ C�� C��3C�33Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C��C�� C�  D � D  D@ D� D� DfD@ D	� D
� D  D@ Dy�D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D��D!@ D"� D#� D%  D&FfD'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJy�DK��DL��DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ�fD\fD]@ D^�fD_�fDafDb@ Dcy�Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D�  D�� D�` D�3D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D��D���D�` D�3D��3D�C3D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�c3D�3Dã3D�C3D�� Dŀ D�  D�� D�` D�3DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�c3D�  DҜ�D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� D�|�D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D�� D�` D�  D� D�@ D�� D�|�D��D��D�\�D�  D� D�@ D�� D�|�D��D�� D�` D�  D� D�@ D���D�|�D�  D�� D�` D�  D� D�C3D��3D�3D�  D��3D�` D�  D��3D�@ D�� D�� D�Ff1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@ff@@  @�  @�  @�  @�  A   AffA   A0  A>ffAP  A`  Ap  A�  A���A���A���A�  A�  A�  A���A�  A�  A�  A�  A�  A�  A�  A�33B   B��B  B��B��B  B  B  B   B$  B(  B,  B0  B3��B8  B<  B@  BD  BH  BL  BP  BT  BXffB\ffB`  Bd  Bh  BlffBp  Bt  Bx  B|ffB�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B���B�  B�33B�  B�  C�fC� C�C	��C�C��C  C� C  C� C�C��C   C"� C%�C'��C*�C,��C/�C1��C4�C6� C8�fC;ffC>  C@� CC  CE� CH  CJ� CM  CO� CR  CT��CW  CY� C\  C^� Ca�Cc� Cf  Ch� Ck  Cm� Cp�Cr� Cu  Cw� Cy�fC|� C  C�� C�  C�@ C�� C�� C��3C�33C�s3C��3C��3C�33C�s3C�� C�  C�@ C�� C�� C�  C�@ C�� C���C��C�L�C�� C��3C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C��3C�  C�L�Cŀ C�� C��3C�33Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C��C�� C�  D � D  D@ D� D� DfD@ D	� D
� D  D@ Dy�D� D  D@ D� D� D  D@ D� D� D  D@ D� D� D��D!@ D"� D#� D%  D&FfD'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJy�DK��DL��DN@ DO� DP� DR  DS@ DT� DU� DW  DX@ DY� DZ�fD\fD]@ D^�fD_�fDafDb@ Dcy�Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D�  D�� D�` D�3D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D��D���D�` D�3D��3D�C3D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�c3D�3Dã3D�C3D�� Dŀ D�  D�� D�` D�3DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�c3D�  DҜ�D�@ D�� DԀ D�  D�� D�` D�  Dנ D�@ D�� D�|�D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D�� D�` D�  D� D�@ D�� D�|�D��D��D�\�D�  D� D�@ D�� D�|�D��D�� D�` D�  D� D�@ D���D�|�D�  D�� D�` D�  D� D�C3D��3D�3D�  D��3D�` D�  D��3D�@ D�� D�� D�Ff1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@Դ9@Ԭ@ԣ�@ԋD@ԋD@ԛ�@���@��`@Ԭ@ԛ�@Ԭ@ԛ�@�z�@ԓu@���@��`@ԓu@�j@�(�@�(�@���@���@Դ9@ԣ�@Դ9@��/@���@ԛ�@Ԭ@�A�@җ�@�~�@�^5@�%@ϝ�@��@�^5@�C�@�J@���@Ĭ@ă@�1'@�  @���@��@ă@���@Ł@š�@�/@ģ�@Ĭ@ļj@���@���@Ĵ9@�Q�@öF@���@��@��^@�J@�~�@��@�;d@�\)@���@��m@��m@��;@��;@��
@���@��F@���@�|�@�S�@�ȴ@���@���@�$�@��H@���@�@�"�@�@�+@���@��@�t�@�l�@�S�@�;d@��@��!@���@�33@�o@��y@��@�\)@���@��!@���@�ff@�=q@�V@�~�@�^5@��@�?}@���@��@�7L@�x�@�@���@��+@�M�@��#@��7@���@���@��D@��@���@�|�@�"�@��!@��!@��H@�V@���@�p�@�p�@�7L@���@��D@�Q�@��@���@���@�~�@��@��7@�G�@��@�%@��`@��9@���@�r�@�Z@�b@���@�S�@��@�5?@��#@��#@��@��@��@��h@�Z@��@�ƨ@�9X@�A�@�1'@� �@���@�o@���@��\@�M�@�@�7L@��/@���@��@�bN@��m@�t�@���@��w@�33@��H@��\@��@��@���@��7@�O�@�/@�V@��/@���@��j@��u@�Q�@�1@��@�+@�ȴ@���@�n�@��@��T@���@�x�@�X@�7L@�V@��j@�1'@�b@��@�  @��@�ƨ@�t�@�\)@�K�@�o@��@���@��R@�M�@�5?@�-@��@�$�@��#@�hs@�%@���@��j@��u@�j@�9X@��@�1@��w@���@�|�@�\)@�33@�
=@��H@���@�v�@�V@�{@��^@�X@�?}@��@��`@�Z@��@��
@���@�l�@��y@�^5@�J@���@�p�@�/@��/@�1'@��@���@���@���@��+@�M�@��T@���@��@�/@��`@��@�1'@��@���@�+@��R@���@��+@�-@���@��-@��@�x�@�X@�G�@���@���@��D@�(�@���@�dZ@�K�@�
=@��@��H@��!@�~�@�n�@�5?@��^@��h@��7@�G�@��@��9@�bN@�b@��@��
@���@�dZ@�dZ@�l�@�S�@�K�@���@��!@�~�@�ff@��@���@���@�hs@�&�@��j@�bN@�Q�@� �@l�@
=@~�y@~E�@~5?@~E�@~@}��@}/@}V@|��@|Z@|�@|�@|�@{��@{�m@{��@{��@{C�@{o@z�@z��@z~�@z^5@z-@y�@y�^@y�7@yhs@y�@y�@x�`@xĜ@x�u@xbN@x  @w��@x  @xA�@x��@x��@x��@x�`@xĜ@x�`@x�`@x�`@x��@y�@y&�@y&�@y&�@y&�@y�@y�@y%@y%@y&�@x�`@xb@w��@v��@vv�@vȴ@v$�@tj@st�@sdZ@s@s@r�!@r��@r�@r��@r�!@r�!@r�\@q�@q�7@q��@r=q@r��@r~�@rn�@s@s��@tj@tz�@xĜ@y�^@x�u@w��@w�P@w�@w�@w��@w�@w|�@w
=@v�@v�y@x1'@w��@wK�@v�R@v�+@vff@vE�@u�@u?}@t�@tj@t9X@tZ@tz�@tz�@t�@t�D@sƨ@s��@s�@s�m@t(�@s�
@sC�@s"�@s�F@u�-@t�@tz�@t(�@tj@tj@t�@t1@s��@s�F@s��@sS�@r�@r�\@rM�@rJ@q�@q�7@q�@p��@p�u@pbN@p �@o�w@oK�@o
=@n�@n�+@n5?@m�@mO�@mV@l��@l(�@k�
@kdZ@k@j��@i�#@i&�@hĜ@hA�@g�;@g�w@g\)@g
=@fȴ@f�R@f��@fv�@f{@e�h@d�/@d�D@d�@cdZ@b�!@b��@b~�@b~�@b=q@bJ@a��@a��@a�@a��@ahs@a%@`��@`A�@_��@_K�@^ȴ@^5?@^@]�T@]@]p�@]V@\�@\Z@[ƨ@[��@[S�@Z�H@Z�\1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@Դ9@Ԭ@ԣ�@ԋD@ԋD@ԛ�@���@��`@Ԭ@ԛ�@Ԭ@ԛ�@�z�@ԓu@���@��`@ԓu@�j@�(�@�(�@���@���@Դ9@ԣ�@Դ9@��/@���@ԛ�@Ԭ@�A�@җ�@�~�@�^5@�%@ϝ�@��@�^5@�C�@�J@���@Ĭ@ă@�1'@�  @���@��@ă@���@Ł@š�@�/@ģ�@Ĭ@ļj@���@���@Ĵ9@�Q�@öF@���@��@��^@�J@�~�@��@�;d@�\)@���@��m@��m@��;@��;@��
@���@��F@���@�|�@�S�@�ȴ@���@���@�$�@��H@���@�@�"�@�@�+@���@��@�t�@�l�@�S�@�;d@��@��!@���@�33@�o@��y@��@�\)@���@��!@���@�ff@�=q@�V@�~�@�^5@��@�?}@���@��@�7L@�x�@�@���@��+@�M�@��#@��7@���@���@��D@��@���@�|�@�"�@��!@��!@��H@�V@���@�p�@�p�@�7L@���@��D@�Q�@��@���@���@�~�@��@��7@�G�@��@�%@��`@��9@���@�r�@�Z@�b@���@�S�@��@�5?@��#@��#@��@��@��@��h@�Z@��@�ƨ@�9X@�A�@�1'@� �@���@�o@���@��\@�M�@�@�7L@��/@���@��@�bN@��m@�t�@���@��w@�33@��H@��\@��@��@���@��7@�O�@�/@�V@��/@���@��j@��u@�Q�@�1@��@�+@�ȴ@���@�n�@��@��T@���@�x�@�X@�7L@�V@��j@�1'@�b@��@�  @��@�ƨ@�t�@�\)@�K�@�o@��@���@��R@�M�@�5?@�-@��@�$�@��#@�hs@�%@���@��j@��u@�j@�9X@��@�1@��w@���@�|�@�\)@�33@�
=@��H@���@�v�@�V@�{@��^@�X@�?}@��@��`@�Z@��@��
@���@�l�@��y@�^5@�J@���@�p�@�/@��/@�1'@��@���@���@���@��+@�M�@��T@���@��@�/@��`@��@�1'@��@���@�+@��R@���@��+@�-@���@��-@��@�x�@�X@�G�@���@���@��D@�(�@���@�dZ@�K�@�
=@��@��H@��!@�~�@�n�@�5?@��^@��h@��7@�G�@��@��9@�bN@�b@��@��
@���@�dZ@�dZ@�l�@�S�@�K�@���@��!@�~�@�ff@��@���@���@�hs@�&�@��j@�bN@�Q�@� �@l�@
=@~�y@~E�@~5?@~E�@~@}��@}/@}V@|��@|Z@|�@|�@|�@{��@{�m@{��@{��@{C�@{o@z�@z��@z~�@z^5@z-@y�@y�^@y�7@yhs@y�@y�@x�`@xĜ@x�u@xbN@x  @w��@x  @xA�@x��@x��@x��@x�`@xĜ@x�`@x�`@x�`@x��@y�@y&�@y&�@y&�@y&�@y�@y�@y%@y%@y&�@x�`@xb@w��@v��@vv�@vȴ@v$�@tj@st�@sdZ@s@s@r�!@r��@r�@r��@r�!@r�!@r�\@q�@q�7@q��@r=q@r��@r~�@rn�@s@s��@tj@tz�@xĜ@y�^@x�u@w��@w�P@w�@w�@w��@w�@w|�@w
=@v�@v�y@x1'@w��@wK�@v�R@v�+@vff@vE�@u�@u?}@t�@tj@t9X@tZ@tz�@tz�@t�@t�D@sƨ@s��@s�@s�m@t(�@s�
@sC�@s"�@s�F@u�-@t�@tz�@t(�@tj@tj@t�@t1@s��@s�F@s��@sS�@r�@r�\@rM�@rJ@q�@q�7@q�@p��@p�u@pbN@p �@o�w@oK�@o
=@n�@n�+@n5?@m�@mO�@mV@l��@l(�@k�
@kdZ@k@j��@i�#@i&�@hĜ@hA�@g�;@g�w@g\)@g
=@fȴ@f�R@f��@fv�@f{@e�h@d�/@d�D@d�@cdZ@b�!@b��@b~�@b~�@b=q@bJ@a��@a��@a�@a��@ahs@a%@`��@`A�@_��@_K�@^ȴ@^5?@^@]�T@]@]p�@]V@\�@\Z@[ƨ@[��@[S�@Z�H@Z�\1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B�B�B �B"�B%�B(�B'�B(�B(�B)�B)�B+B-B1'B5?B7LB8RB8RB6FB6FB7LB7LB7LB6FB33B0!B)�B%�B"�B&�B(�B.B0!B1'B49B6FB7LB8RB8RB8RB8RB8RB8RB8RB8RB:^B;dB<jB=qB?}BA�BB�BD�BD�BE�BI�BJ�BJ�BJ�BJ�BI�BH�BG�BF�B?}B>wB=qB;dB6FB49B49B5?B5?B6FB7LB:^B<jB<jB:^B:^B?}B=qBD�BF�BF�BJ�BJ�BH�BG�BG�BG�BC�BC�BA�B@�B?}B=qB>wB?}B=qB;dB:^B:^B:^B9XB8RB8RB6FB5?B49B2-B1'B0!B0!B/B/B/B/B/B/B/B.B-B,B+B(�B'�B(�B(�B(�B'�B%�B �B�B �B#�B#�B#�B"�B!�B�B�B�B�B�B�B�B�B�B�B{BuB�B�B�B{BuBoBhBbB\B\BVBVBVBVBPBPBJBDB
=B	7B1B+B+B%BBBBBBBBB  B  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�mB�fB�`B�`B�TB�HB�HB�BB�;B�5B�)B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBȴBǮBǮBƨBƨBƨBŢBŢBĜBÖBÖBB��B��B��B�}B�}B�}B�wB�qB�qB�jB�dB�dB�^B�^B�XB�RB�LB�FB�FB�?B�?B�9B�?B�9B�9B�9B�3B�3B�-B�'B�'B�!B�!B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�XB�jB�^B�RB�XB�XB�XB�^B�dB�dB�^B�^B�dB�}B�}B�wB�qB�qB�qB�wB�qB�jB�jB�jB�qB�wB�}B��B��B��B��B��B��BBÖBÖBBBŢB��BɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBǮBƨBŢBŢBĜBŢBƨBƨBƨBƨBƨBǮBǮBǮBǮBƨBƨBŢBŢBŢBŢBŢBƨBƨBƨBƨBƨBǮBȴBǮBǮBƨBŢBƨBŢBŢBŢBŢBŢBŢBŢBĜBĜBĜBĜBĜBĜBĜ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111BBBB%B%BBBBBBB%B%BBBBBB%B%BBBBBBBBBBBBB�BB BDBPBbB uBoB uB uB!{B!{B"�B$�B(�B,�B.�B/�B/�B-�B-�B.�B.�B.�B-�B*�B'�B!{BcBQBiB vB%�B'�B(�B+�B-�B.�B/�B/�B/�B/�B/�B/�B/�B/�B1�B2�B3�B4�B6�B9B:B<B<B= BA8BB?BB?BB?BB?BA8B@2B?,B>&B6�B5�B4�B2�B-�B+�B+�B,�B,�B-�B.�B1�B3�B3�B1�B1�B6�B4�B<B>'B>&BB?BB?B@2B?,B?,B?,B;B;B9B8B6�B4�B5�B6�B4�B2�B1�B1�B1�B0�B/�B/�B-�B,�B+�B)�B(�B'�B'�B&�B&�B&�B&�B&�B&�B&�B%�B$�B#�B"�B wBqB wB wB wBqBdBGB4BGBXBXBXBSBMB@B:B4B.B"BBB
BB
B�B
�BB
BB�B
�B	�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�yB�yB�yB�sB�sB�lB�fB�fB�fB�fB�`B�`B�TB�NB�NB�HB�HB�AB�AB�<B�<B�6B�6B�0B�0B�*B�*B�$B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��BռBӱBѥBПBϘB͌B̆B̆BˀB�zB�tB�tB�nB�hB�bB�\B�VB�PB�JB�CB�CB�=B�7B�7B�1B�1B�1B�,B�,B�&B� B� B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�|B��B�|B�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�jB�jB�jB�jB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�jB�vB�|B�|B�|B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�jB�^B�XB�XB�XB�XB�^B�^B�^B�^B�^B�^B�XB�XB�^B�jB�pB�jB�pB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B�B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�!B�!B�B�B�-B�KB�DB�DB�DB�QB�WB�WB�]B�WB�]B�]B�]B�]B�]B�]B�WB�WB�WB�WB�WB�WB�WB�QB�QB�WB�WB�WB�WB�WB�RB�RB�RB�RB�LB�LB�EB�EB�EB�9B�3B�-B�-B�'B�-B�3B�3B�3B�3B�3B�9B�9B�9B�9B�3B�3B�-B�-B�-B�-B�-B�3B�3B�3B�3B�3B�9B�?B�9B�9B�3B�-B�3B�-B�-B�-B�.B�.B�.B�.B�(B�(B�(B�(B�(B�(B�(1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0082766                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904512018110709045120181107090451  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172719  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172719  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090451  IP  PSAL            @ffD�FfG�O�                