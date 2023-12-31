CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:07Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z        x  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    B   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        x  D4   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    L�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     x  N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  WD   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    _�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  a�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    jT   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  lt   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  t�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    }d   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     x  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �    HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �@   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �P   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �T   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �d   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �h   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �l   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �p   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090439  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�!����1   @�!�w�@N�^��P�BDV�9>1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A^ffAp  A�  A���A���A���A���A�  A�  A�  A�  A�  A�  A�  A���A���A�  A�  B   B  B  B  B  B  B��B��B��B#��B(  B,ffB0ffB4ffB8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B{��B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C��C  C� C�fC� C   C"� C%  C'� C*  C,� C.�fC1ffC4  C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  COffCR  CT� CW  CY��C\  C^� Ca  Cc� Cf  ChffCk  Cm� Co�fCr� Cu  CwffCz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�L�C���C���C��C�L�CŌ�C�� C�  C�@ Cʀ C�� C��3C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�33C� C���C�  C�@ C� C�� C�  C�@ C�� C���C��C�s3C�  D � D  D@ D� D� D  D@ D	� D
� D  DFfD�fD� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"�fD#� D%  D&@ D'y�D(��D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D59�D6y�D7��D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJy�DK� DM  DN@ DO� DP� DR  DSFfDT� DU�fDWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db9�Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dmy�Dn��Dp  Dq@ Dr� Ds� Du  Dv@ Dwy�Dx� Dz  D{@ D|� D}� D  D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D��3D�#3D�� D�c3D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�\�D���D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�<�D�� Dŀ D�  D��3D�` D�  Dȣ3D�@ D��3Dʀ D�  D�� D�c3D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�<�D���DԀ D�  D�� D�` D�  Dנ D�@ D���Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�<�D�� D�3D�#3D�� D�\�D�  D� D�@ D�� D� D�&fD��fD�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A^ffAp  A�  A���A���A���A���A�  A�  A�  A�  A�  A�  A�  A���A���A�  A�  B   B  B  B  B  B  B��B��B��B#��B(  B,ffB0ffB4ffB8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B{��B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B���B���B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C��C  C� C�fC� C   C"� C%  C'� C*  C,� C.�fC1ffC4  C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  COffCR  CT� CW  CY��C\  C^� Ca  Cc� Cf  ChffCk  Cm� Co�fCr� Cu  CwffCz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C��C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�L�C���C���C��C�L�CŌ�C�� C�  C�@ Cʀ C�� C��3C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�33C� C���C�  C�@ C� C�� C�  C�@ C�� C���C��C�s3C�  D � D  D@ D� D� D  D@ D	� D
� D  DFfD�fD� D  D@ D� D� D  D@ D� D� D  D@ D� D� D   D!@ D"�fD#� D%  D&@ D'y�D(��D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D59�D6y�D7��D9  D:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJy�DK� DM  DN@ DO� DP� DR  DSFfDT� DU�fDWfDX@ DY� DZ� D\  D]@ D^� D_� Da  Db9�Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dmy�Dn��Dp  Dq@ Dr� Ds� Du  Dv@ Dwy�Dx� Dz  D{@ D|� D}� D  D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D��3D�#3D�� D�c3D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D���D�� D�@ D�� D��3D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�\�D���D�� D�@ D�� D�� D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D��3D��3D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  Dà D�<�D�� Dŀ D�  D��3D�` D�  Dȣ3D�@ D��3Dʀ D�  D�� D�c3D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�<�D���DԀ D�  D�� D�` D�  Dנ D�@ D���Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�#3D�� D�` D�  D��D�@ D�� D� D�  D�� D�` D�  D� D�<�D�� D�3D�#3D�� D�\�D�  D� D�@ D�� D� D�&fD��fD�3311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��P@�|�@���@�b@�r�@���@���@���@��@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�X@�O�@�X@�X@�X@�O�@�O�@�O�@�G�@�G�@�G�@�G�@�?}@�7L@�/@�&�@�/@��@��@�%@���@��@�Ĝ@��j@��@��@��9@��9@��j@��j@��9@���@���@���@��u@��u@��u@��u@��@��@�z�@�r�@�r�@�j@�Q�@�Z@�I�@�I�@�I�@�A�@�A�@�9X@�1'@�(�@�(�@�(�@� �@� �@��@��@�  @��m@�ƨ@��@�33@��@��R@���@��\@�ff@�M�@�5?@�=q@�$�@�@�@�J@�@���@���@�p�@�&�@���@�r�@���@�ƨ@��@� �@��/@��j@��@�1@��;@��@�  @��F@���@�\)@��y@�J@��^@���@���@�hs@�V@�bN@�;d@��!@�v�@�ff@�^5@�-@��@��@��@�V@�V@��@��@��`@�1'@��;@�ƨ@�t�@�"�@��@�v�@���@��^@���@�x�@���@��7@���@��T@���@��^@�`B@��@��9@�9X@���@�33@��y@���@���@�M�@�~�@�@���@���@���@�Ĝ@�Z@��@��F@�\)@�33@���@�ȴ@��!@��\@��\@�=q@�@���@��7@�V@���@��u@��@�bN@�b@��w@�\)@�@��@���@��y@��P@��P@��+@�J@��T@��^@�x�@�O�@�O�@�G�@�?}@�/@�/@�`B@���@��@�`B@�7L@�%@��@�Ĝ@���@��u@��u@��u@��@�I�@�1'@�b@���@���@��@���@�|�@�l�@�\)@�;d@�33@�+@�
=@���@��@��H@���@���@���@��\@��\@��\@�~�@�n�@�V@�=q@�$�@��@�{@��@��^@��-@���@���@���@���@�hs@��@���@�Ĝ@���@�bN@�b@�@��@�w@|�@l�@l�@
=@~�R@~ff@~5?@}�@}�-@}?}@}/@|��@|�j@|(�@{��@{�@{�@{S�@{33@{o@z~�@zJ@yhs@yx�@yhs@y�@x�u@xb@xb@w��@w�@v�@v��@v��@w\)@w
=@v�y@vff@v5?@u�T@u@u?}@uV@t�@tI�@t�@s�F@s�@s33@r�H@r�\@r=q@q��@q��@q�^@q��@qx�@qG�@q7L@q&�@p��@pĜ@p�9@p��@p�9@q�@q7L@qX@qX@qx�@q�7@qx�@q�7@q��@q�^@q�^@q�#@q��@r�@r=q@r�\@r��@r�H@s@sC�@sdZ@s�@s��@s�F@s�m@t9X@t�j@t��@t��@t�@t��@u?}@uO�@u�@t�D@tz�@t�j@t�@t�@t�j@tz�@t�@t�D@t�j@t�/@tI�@t(�@s�
@sƨ@sS�@st�@sS�@r�@r��@r��@r�!@r=q@q��@q�@q�#@q��@q�^@qx�@qG�@q&�@qX@qX@p�`@pbN@pQ�@pA�@p1'@o�@o��@o��@o|�@o+@o�@nv�@n��@nȴ@n�R@nȴ@nff@n5?@n5?@n$�@nE�@n5?@m�-@l��@l��@l�@l��@lj@lZ@l��@l9X@k�
@k�F@kS�@ko@k@j��@j�\@jJ@i��@i�7@ihs@i�@h��@hr�@hA�@h �@hb@g|�@g��@gK�@f��@f�+@f5?@e�T@e�h@e�@e/@d�j@dj@d9X@c�
@c��@c�@c33@b�@b��@b=q@a��@a�#@a�^@a�7@`��@`Q�@`b@`  @`  @_�;@_��@_\)@_
=@^ff@^@]�-@]O�@\�j@\j@\�@[�m@[��@[S�@["�@Z�H@Z��@Z^5@Z=q@Z-@ZJ@Y��@Y&�@X�`@X�@XQ�@X1'@W�@Wl�@WK�@V��@V��@V5?@V{@U�T@U�-@U�h@U�h@Up�@U�@T�@T�@T�D@TI�@T�@S�F@S��@St�@S33@S@R��@R��@Rn�@R�@RJ@Q��@Q��@QG�@Q%@P��@PQ�@P1'@P  @O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��P@�|�@���@�b@�r�@���@���@���@��@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�X@�O�@�X@�X@�X@�O�@�O�@�O�@�G�@�G�@�G�@�G�@�?}@�7L@�/@�&�@�/@��@��@�%@���@��@�Ĝ@��j@��@��@��9@��9@��j@��j@��9@���@���@���@��u@��u@��u@��u@��@��@�z�@�r�@�r�@�j@�Q�@�Z@�I�@�I�@�I�@�A�@�A�@�9X@�1'@�(�@�(�@�(�@� �@� �@��@��@�  @��m@�ƨ@��@�33@��@��R@���@��\@�ff@�M�@�5?@�=q@�$�@�@�@�J@�@���@���@�p�@�&�@���@�r�@���@�ƨ@��@� �@��/@��j@��@�1@��;@��@�  @��F@���@�\)@��y@�J@��^@���@���@�hs@�V@�bN@�;d@��!@�v�@�ff@�^5@�-@��@��@��@�V@�V@��@��@��`@�1'@��;@�ƨ@�t�@�"�@��@�v�@���@��^@���@�x�@���@��7@���@��T@���@��^@�`B@��@��9@�9X@���@�33@��y@���@���@�M�@�~�@�@���@���@���@�Ĝ@�Z@��@��F@�\)@�33@���@�ȴ@��!@��\@��\@�=q@�@���@��7@�V@���@��u@��@�bN@�b@��w@�\)@�@��@���@��y@��P@��P@��+@�J@��T@��^@�x�@�O�@�O�@�G�@�?}@�/@�/@�`B@���@��@�`B@�7L@�%@��@�Ĝ@���@��u@��u@��u@��@�I�@�1'@�b@���@���@��@���@�|�@�l�@�\)@�;d@�33@�+@�
=@���@��@��H@���@���@���@��\@��\@��\@�~�@�n�@�V@�=q@�$�@��@�{@��@��^@��-@���@���@���@���@�hs@��@���@�Ĝ@���@�bN@�b@�@��@�w@|�@l�@l�@
=@~�R@~ff@~5?@}�@}�-@}?}@}/@|��@|�j@|(�@{��@{�@{�@{S�@{33@{o@z~�@zJ@yhs@yx�@yhs@y�@x�u@xb@xb@w��@w�@v�@v��@v��@w\)@w
=@v�y@vff@v5?@u�T@u@u?}@uV@t�@tI�@t�@s�F@s�@s33@r�H@r�\@r=q@q��@q��@q�^@q��@qx�@qG�@q7L@q&�@p��@pĜ@p�9@p��@p�9@q�@q7L@qX@qX@qx�@q�7@qx�@q�7@q��@q�^@q�^@q�#@q��@r�@r=q@r�\@r��@r�H@s@sC�@sdZ@s�@s��@s�F@s�m@t9X@t�j@t��@t��@t�@t��@u?}@uO�@u�@t�D@tz�@t�j@t�@t�@t�j@tz�@t�@t�D@t�j@t�/@tI�@t(�@s�
@sƨ@sS�@st�@sS�@r�@r��@r��@r�!@r=q@q��@q�@q�#@q��@q�^@qx�@qG�@q&�@qX@qX@p�`@pbN@pQ�@pA�@p1'@o�@o��@o��@o|�@o+@o�@nv�@n��@nȴ@n�R@nȴ@nff@n5?@n5?@n$�@nE�@n5?@m�-@l��@l��@l�@l��@lj@lZ@l��@l9X@k�
@k�F@kS�@ko@k@j��@j�\@jJ@i��@i�7@ihs@i�@h��@hr�@hA�@h �@hb@g|�@g��@gK�@f��@f�+@f5?@e�T@e�h@e�@e/@d�j@dj@d9X@c�
@c��@c�@c33@b�@b��@b=q@a��@a�#@a�^@a�7@`��@`Q�@`b@`  @`  @_�;@_��@_\)@_
=@^ff@^@]�-@]O�@\�j@\j@\�@[�m@[��@[S�@["�@Z�H@Z��@Z^5@Z=q@Z-@ZJ@Y��@Y&�@X�`@X�@XQ�@X1'@W�@Wl�@WK�@V��@V��@V5?@V{@U�T@U�-@U�h@U�h@Up�@U�@T�@T�@T�D@TI�@T�@S�F@S��@St�@S33@S@R��@R��@Rn�@R�@RJ@Q��@Q��@QG�@Q%@P��@PQ�@P1'@P  @O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�%B�=B�7B�=B�DB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B�'B�-B�-B�9B�LB�XB�XB�RB�RB�LB�LB�FB�FB�?B�?B�9B�-B�'B�'B�'B�'B�'B�!B�!B�!B�!B�'B�'B�-B�3B�3B�-B�3B�3B�3B�3B�-B�'B�!B�!B�!B�!B�-B�-B�3B�?B�?B�?B�9B�9B�-B�'B�!B�B�B�B�B�B�!B�3B�FB�!B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�uB�oB�oB�oB�oB�hB�bB�hB�bB�\B�\B�\B�bB�oB�hB�hB�bB�bB�bB�bB�\B�\B�\B�VB�VB�VB�VB�PB�PB�JB�JB�JB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�PB�PB�VB�\B�bB�bB�hB�oB�oB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�!B�!B�!B�!B�'B�'B�'B�-B�-B�-B�3B�9B�9B�?B�9B�?B�?B�FB�FB�LB�LB�LB�LB�RB�LB�RB�XB�XB�^B�^B�^B�^B�dB�dB�jB�dB�dB�dB�jB�jB�dB�jB�qB�jB�jB�qB�qB�jB�qB�qB�qB�qB�qB�qB�wB�wB�wB�wB�wB�wB�wB�qB�wB�wB�qB�qB�qB�qB�wB�wB�wB�qB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�qB�qB�qB�qB�jB�jB�jB�qB�qB�qB�qB�qB�jB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�jB�dB�dB�dB�dB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�^B�^11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B~�B~�B�B�B�B�B�/B�)B�/B�5B�YB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�fB�`B�`B�`B�`B�`B�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�lB�lB�lB�lB�lB�lB�lB�rB�rB�rB�rB�rB�xB�xB�rB�xB�xB�xB�xB�xB�xB�~B��B��B��B��B��B��B�lB��B��B�B�B�B�*B�=B�IB�IB�CB�CB�=B�=B�7B�7B�0B�0B�*B�B�B�B�B�B�B�B�B�B�B�B�B�B�$B�$B�B�$B�$B�$B�$B�B�B�B�B�B�B�B�B�$B�0B�0B�0B�*B�*B�B�B�B�B�B�B�B�B�B�$B�7B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�yB�B�B�yB�yB�sB�sB�sB�sB�mB�mB�gB�gB�aB�aB�aB�aB�ZB�TB�ZB�TB�NB�NB�NB�TB�aB�ZB�ZB�TB�TB�TB�TB�NB�NB�NB�HB�HB�HB�HB�CB�CB�=B�=B�=B�7B�7B�7B�=B�=B�=B�=B�=B�=B�=B�CB�CB�IB�NB�TB�TB�ZB�aB�aB�gB�mB�mB�sB�sB�yB�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�+B�+B�1B�+B�1B�1B�8B�8B�>B�>B�>B�>B�DB�>B�DB�JB�JB�PB�PB�PB�PB�VB�VB�\B�VB�VB�VB�\B�\B�VB�\B�cB�\B�\B�cB�cB�\B�cB�cB�cB�cB�cB�cB�iB�iB�iB�iB�iB�iB�iB�cB�iB�iB�cB�cB�cB�cB�iB�iB�iB�cB�cB�cB�cB�cB�cB�cB�cB�\B�\B�\B�cB�cB�cB�cB�\B�\B�\B�cB�cB�cB�cB�cB�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�VB�VB�VB�VB�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�P11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0020074                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904392018110709043920181107090439  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172707  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172707  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090439  IP  PSAL            @ffD�33G�O�                