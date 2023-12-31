CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:04Z last update (coriolis COQC software)   
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
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090437  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @��6�j1   @���5-@N�.;G���Aj����1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@Fff@�33@�  @�  @�  A��AffAffA0  A@  AP  A`  Ap  A�  A���A�  A�  A���A�  A�  A�  A�  A�33A�  A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BPffBT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C�C� C  C� C  C��C   C"ffC%  C'� C*  C,� C/�C1� C4  C6� C9�C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp  Cr� Cu�Cw� Cz  C|� C  C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C�� C�  C�@ C�� C���C�  C�@ Cŀ C�� C�  C�L�Cʌ�C���C��C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C��3C�@ C��C���C�  C�@ C� C�� C��C�@ C�� C�� C�  C�� C��3D � D  DFfD� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D�fDfDFfD�fD�fD fD!FfD"�fD#�fD%  D&9�D'� D(� D*  D+FfD,� D-�fD/fD0FfD1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?FfD@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK��DL��DN@ DO� DP� DR  DS@ DT� DU��DW  DX@ DY�fDZ� D\  D]@ D^�fD_�fDa  Db@ Dc�fDd� Df  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� DpfDq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D��3D�@ D���D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�3D�� D�C3D�� D�|�D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  Dà D�C3D��3Dŀ D��D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�\�D���DҠ D�@ D�� DԀ D�#3D�� D�\�D�  Dנ D�@ D��3Dـ D��D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�3D� D�C3D�� D� D��D�� D�c3D�  D� D�@ D�� D� D�#3D��3D�c3D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�#3D�� D�` D�3D��3D�@ D�� D�� D�  D���D�` D�  D�� D�FfD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@Fff@�33@�  @�  @�  A��AffAffA0  A@  AP  A`  Ap  A�  A���A�  A�  A���A�  A�  A�  A�  A�33A�  A���A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BPffBT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�33B�  B�  B�  B�  B�33B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C�C� C  C� C  C��C   C"ffC%  C'� C*  C,� C/�C1� C4  C6� C9�C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  ChffCk  Cm��Cp  Cr� Cu�Cw� Cz  C|� C  C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C���C�� C�  C�@ C�� C���C�  C�@ Cŀ C�� C�  C�L�Cʌ�C���C��C�@ C�s3C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C��C�@ C� C�� C�  C�@ C� C�� C��3C�@ C��C���C�  C�@ C� C�� C��C�@ C�� C�� C�  C�� C��3D � D  DFfD� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D�fDfDFfD�fD�fD fD!FfD"�fD#�fD%  D&9�D'� D(� D*  D+FfD,� D-�fD/fD0FfD1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?FfD@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK��DL��DN@ DO� DP� DR  DS@ DT� DU��DW  DX@ DY�fDZ� D\  D]@ D^�fD_�fDa  Db@ Dc�fDd� Df  Dg@ Dh�fDi� Dk  Dl@ Dm�fDn� DpfDq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D��3D�@ D���D�|�D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�3D�� D�C3D�� D�|�D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  Dà D�C3D��3Dŀ D��D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�\�D���DҠ D�@ D�� DԀ D�#3D�� D�\�D�  Dנ D�@ D��3Dـ D��D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�3D� D�C3D�� D� D��D�� D�c3D�  D� D�@ D�� D� D�#3D��3D�c3D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�#3D�� D�` D�3D��3D�@ D�� D�� D�  D���D�` D�  D�� D�FfD���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�#A�#A�#A�#A�#A�#A�#A�#A�#A�;A�;A�;A�;A�mA�A�mA�mA�mA�TA�mA�A�A�A�A�A�A�A�A�A�A�A�A�A�A�A�A��A��A��A��A�A�
A�
A��A��AƨA�FA�FA�FA��Al�A�A��A1'A��A(�A�H@��
@�E�@���@���@��@��H@���@�S�@�hs@�A�@���@��@��
@ۅ@��@��@�l�@Ұ!@�J@ѩ�@щ7@�?}@��#@���@���@��@�@�V@Ұ!@���@�v�@҇+@ҸR@җ�@�n�@�5?@�=q@���@ӝ�@�C�@�t�@�t�@�-@��#@д9@�Q�@���@��m@͡�@̬@��`@��@˶F@�K�@���@�?}@�/@���@�ƨ@� �@�hs@�
=@�n�@Ł@��@�G�@�V@Ų-@��@ũ�@ŉ7@�X@ă@�K�@�-@��T@�Ĝ@�bN@�r�@�|�@��@��@��T@��@�?}@��j@��F@��@�o@�"�@��!@�n�@�?}@���@�z�@��@�(�@�dZ@�$�@�X@�-@���@��9@���@�A�@��;@���@�o@�5?@�X@�&�@��/@�|�@�\)@�l�@�|�@�+@��y@�ff@��T@���@��@�?}@���@��D@�1'@� �@��@���@�=q@��F@���@�X@��@��@��R@�@���@�t�@�ff@��@���@���@��9@�1'@��@�@��@�?}@��@��`@���@���@�A�@��m@��@��H@�M�@��^@��^@�-@��+@�\)@�%@�(�@�;d@�ƨ@�"�@�n�@�J@���@��
@��\@��@�?}@�dZ@�V@�@�@��@�?}@��@��;@�S�@�33@��@�~�@�=q@�@�O�@�X@��@�Ĝ@�  @���@�K�@�o@��@���@�=q@���@��h@���@��@�l�@�ȴ@�-@�`B@��/@��u@��
@�33@��H@���@�V@���@�p�@��@���@�z�@���@��w@��@�"�@��@��+@�$�@��@��^@�hs@��@���@�r�@���@��m@�1'@�b@��@��F@��@�K�@�"�@��@���@��R@�n�@�=q@��@��h@�`B@�?}@�&�@�%@���@�%@���@���@���@��D@�j@�j@�(�@�bN@�I�@�bN@�(�@|�@~ȴ@~�R@~��@~��@~5?@}�@}p�@}�-@|�/@|(�@|j@|��@|1@{�@{��@|1@|1@{�m@{t�@z�@z��@z=q@y�^@y�#@y�7@yx�@y��@y��@y7L@xĜ@xr�@x �@w�@w�P@w�@v�R@vV@v{@v��@v��@vȴ@v�R@v@t�@tj@t�j@t�@t��@u�@up�@u`B@t�@s�F@so@r~�@r�H@s@r�!@r�H@s33@sC�@s�@s�
@t1@s�
@s��@s��@s�F@sƨ@s�
@t(�@tj@tz�@tz�@t�D@t�j@uV@u�h@u�-@u�-@u��@v$�@vV@vV@vE�@vff@v5?@v@u�-@u��@up�@u�@uV@u/@u/@uV@t�/@t�D@tZ@t9X@t�@t�@t�@t(�@s�m@s�m@s�
@sƨ@s��@st�@sS�@so@r��@r��@r��@r~�@r-@r^5@rM�@r�@q�@q�#@q�^@q��@q��@qX@q7L@p�`@p�9@p��@p�@p�u@p�@pbN@p1'@o�@o�@o�@o��@o�@o�P@o+@n�R@n�+@nv�@nV@nE�@n$�@m��@m�h@m?}@l�/@l�@l�D@lz�@mV@l�@l�D@lj@l(�@k�F@k�@kt�@kC�@j�!@j=q@i�#@i��@iX@i&�@h��@hĜ@h�9@h�9@hr�@hb@g�w@g�P@g\)@g
=@f�y@f�+@fff@fV@e�@e@eO�@eV@d�D@d�@c�m@c�m@c��@c��@ct�@c33@b��@bn�@b^5@b-@ahs@aX@a�@a%@`�9@`  @_+@^��@]��@]�-@]@^@]�h@]`B@]?}@\��@\�j@\I�@\(�@[�F@[��@[�@[o@Z��@Z�!@Z��@Zn�@ZJ@Y��@Yhs@Y&�@X��@XĜ@X�u@Xr�@Xr�@XQ�@X �1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�#A�#A�#A�#A�#A�#A�#A�#A�#A�;A�;A�;A�;A�mA�A�mA�mA�mA�TA�mA�A�A�A�A�A�A�A�A�A�A�A�A�A�A�A�A��A��A��A��A�A�
A�
A��A��AƨA�FA�FA�FA��Al�A�A��A1'A��A(�A�H@��
@�E�@���@���@��@��H@���@�S�@�hs@�A�@���@��@��
@ۅ@��@��@�l�@Ұ!@�J@ѩ�@щ7@�?}@��#@���@���@��@�@�V@Ұ!@���@�v�@҇+@ҸR@җ�@�n�@�5?@�=q@���@ӝ�@�C�@�t�@�t�@�-@��#@д9@�Q�@���@��m@͡�@̬@��`@��@˶F@�K�@���@�?}@�/@���@�ƨ@� �@�hs@�
=@�n�@Ł@��@�G�@�V@Ų-@��@ũ�@ŉ7@�X@ă@�K�@�-@��T@�Ĝ@�bN@�r�@�|�@��@��@��T@��@�?}@��j@��F@��@�o@�"�@��!@�n�@�?}@���@�z�@��@�(�@�dZ@�$�@�X@�-@���@��9@���@�A�@��;@���@�o@�5?@�X@�&�@��/@�|�@�\)@�l�@�|�@�+@��y@�ff@��T@���@��@�?}@���@��D@�1'@� �@��@���@�=q@��F@���@�X@��@��@��R@�@���@�t�@�ff@��@���@���@��9@�1'@��@�@��@�?}@��@��`@���@���@�A�@��m@��@��H@�M�@��^@��^@�-@��+@�\)@�%@�(�@�;d@�ƨ@�"�@�n�@�J@���@��
@��\@��@�?}@�dZ@�V@�@�@��@�?}@��@��;@�S�@�33@��@�~�@�=q@�@�O�@�X@��@�Ĝ@�  @���@�K�@�o@��@���@�=q@���@��h@���@��@�l�@�ȴ@�-@�`B@��/@��u@��
@�33@��H@���@�V@���@�p�@��@���@�z�@���@��w@��@�"�@��@��+@�$�@��@��^@�hs@��@���@�r�@���@��m@�1'@�b@��@��F@��@�K�@�"�@��@���@��R@�n�@�=q@��@��h@�`B@�?}@�&�@�%@���@�%@���@���@���@��D@�j@�j@�(�@�bN@�I�@�bN@�(�@|�@~ȴ@~�R@~��@~��@~5?@}�@}p�@}�-@|�/@|(�@|j@|��@|1@{�@{��@|1@|1@{�m@{t�@z�@z��@z=q@y�^@y�#@y�7@yx�@y��@y��@y7L@xĜ@xr�@x �@w�@w�P@w�@v�R@vV@v{@v��@v��@vȴ@v�R@v@t�@tj@t�j@t�@t��@u�@up�@u`B@t�@s�F@so@r~�@r�H@s@r�!@r�H@s33@sC�@s�@s�
@t1@s�
@s��@s��@s�F@sƨ@s�
@t(�@tj@tz�@tz�@t�D@t�j@uV@u�h@u�-@u�-@u��@v$�@vV@vV@vE�@vff@v5?@v@u�-@u��@up�@u�@uV@u/@u/@uV@t�/@t�D@tZ@t9X@t�@t�@t�@t(�@s�m@s�m@s�
@sƨ@s��@st�@sS�@so@r��@r��@r��@r~�@r-@r^5@rM�@r�@q�@q�#@q�^@q��@q��@qX@q7L@p�`@p�9@p��@p�@p�u@p�@pbN@p1'@o�@o�@o�@o��@o�@o�P@o+@n�R@n�+@nv�@nV@nE�@n$�@m��@m�h@m?}@l�/@l�@l�D@lz�@mV@l�@l�D@lj@l(�@k�F@k�@kt�@kC�@j�!@j=q@i�#@i��@iX@i&�@h��@hĜ@h�9@h�9@hr�@hb@g�w@g�P@g\)@g
=@f�y@f�+@fff@fV@e�@e@eO�@eV@d�D@d�@c�m@c�m@c��@c��@ct�@c33@b��@bn�@b^5@b-@ahs@aX@a�@a%@`�9@`  @_+@^��@]��@]�-@]@^@]�h@]`B@]?}@\��@\�j@\I�@\(�@[�F@[��@[�@[o@Z��@Z�!@Z��@Zn�@ZJ@Y��@Yhs@Y&�@X��@XĜ@X�u@Xr�@Xr�@XQ�@X �1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�1B�1B�1B�1B�+B�+B�1B�7B�7B�7B�=B�7B�7B�=B�7B�7B�7B�7B�=B�=B�=B�7B�7B�=B�JB�JB�DB�JB�PB�PB�PB�PB�PB�PB�DB�=B�7B�7B�1B�%B�B�B�VB�dB�RB�3BBĜB��B��B��B�fB�sB�B�B�B��B��B�B�yB�HB�;B�5B�5B�;B�BB�B�B�B�B��B��B��B��B��B��B��B��B��B��B��BB+B%B+B
=B+B%BBB+B%B��B��B��B  BB  B��B��B��B��B��B��BB��B��B��B��B��B��B  BBBBB��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�sB�mB�fB�`B�fB�ZB�BB�5B�TB�BB�HB�TB�NB�HB�HB�BB�/B�/B�/B�)B�#B�B�B�B�B�B�B�
B�
B�B�B�B�B��B�B�#B�BB�sB�B��B��B
=B	7B+BBB��B��B��B��B��B�B�B�B�B�B�B�B�B�yB�sB�sB�mB�mB�fB�`B�`B�fB�B�B��B��B��B��B��B��B��B��B��B�B�B�TB��B��BǮBǮBƨBƨBŢBÖB��B�}B�}B�wB�wB�qB�jB�jB�qB�qB�jB�dB�XB�XB�XB�XB�RB�LB�LB�FB�?B�9B�3B�-B�'B�!B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�{B��B�{B�uB�{B��B��B��B��B�{B�{B�uB�oB�uB�uB�uB�{B�{B�uB�uB�oB�oB�oB�hB�bB�bB�\B�\B�hB�oB�oB�oB�hB�\B�VB�\B�bB�bB�oB�oB�oB�bB�bB�VB�VB�\B�bB�bB�hB�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�'B�'B�-B�3B�3B�3B�9B�9B�9B�9B�9B�?B�?B�?B�LB�RB�RB�RB�RB�XB�XB�XB�^B�^B�^B�^B�^B�^B�dB�dB�dB�jB�jB�qB�qB�qB�qB�qB�qB�qB�qB�qB�wB�wB�}B�}B�}B�}B�}B�}B�}B��BBÖBBBÖBBB��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�qB�qB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�}B�wB�wB�wB�wB�wB�qB�jB�dB�^B�^B�^B�dB�jB�jB�jB�jB�jB�qB�qB�qB�jB�jB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�wB�}B�}B�}1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�LB�FB�FB�FB�FB�FB�FB�FB�FB�FB�FB�FB�FB�LB�LB�LB�LB�FB�FB�LB�RB�RB�RB�XB�RB�RB�XB�RB�RB�RB�RB�XB�XB�XB�RB�RB�XB�eB�eB�_B�eB�kB�kB�kB�kB�kB�kB�_B�XB�RB�RB�LB�@B�4B�:B�qB�B�mB�NB��B÷B��B��B��B�B�B�B�B��B��B��B�B�B�cB�VB�PB�PB�VB�]B�B�B�B��B��B��B��B��B��B��B��B��B��B��B�B  BEB?BEB	WBEB?B,B3BEB?B�B�B�B�B  B�B�B��B��B��B��B�B  B��B��B��B��B��B�B�B&B&B&B  B�B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B��B�B�B�B�B�B�{B�B�uB�]B�PB�oB�]B�cB�oB�iB�cB�cB�]B�JB�JB�JB�DB�>B�8B�8B�8B�2B�2B�+B�%B�%B�B�B�B�B�B�B�>B�]B�B��B��B�	B	WBQBEB3B !B�B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�{B�{B�B�B�B��B�B�B��B��B�B��B��B��B�B�B�oB�B��B��B��B��B��BĽB±B��B��B��B��B��B��B��B��B��B��B��B�B�sB�sB�sB�sB�mB�gB�gB�aB�ZB�TB�NB�HB�BB�<B�/B�/B�)B�$B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�xB�xB��B��B��B��B��B�xB�rB�xB�~B�~B��B��B��B�~B�~B�rB�rB�xB�~B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�$B�0B�0B�0B�0B�0B�0B�7B�7B�=B�=B�CB�CB�IB�OB�OB�OB�UB�UB�UB�UB�UB�[B�[B�[B�hB�nB�nB�nB�nB�tB�tB�tB�zB�zB�zB�zB�zB�zB�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B±B��B��B±B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB�zB�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.00087271                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904372018110709043720181107090437  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172704  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172704  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090437  IP  PSAL            @ffD���G�O�                