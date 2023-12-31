CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:15Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090447  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�GQ�\(�1   @�GQ�\(�@O���O��@T{s9!h8   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@9��@�  @�  @�  @�  A��AffA   A0  AA��AP  A^ffAp  A�  A�  A�  A�  A�33A�  A���A���A�  A�  A�  A�  A�  A�  A�  A���B   B  BffB  B  B  BffB  B ffB$ffB(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\ffB`  Bd  Bg��Bl  Bp  Bt  BxffB|  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B���B�  B�33B�  B���B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33C  C� C�fC	ffC  C� C�C��C  C��C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CTffCW  CY��C\�C^��Ca  Cc� Cf�Ch� Cj�fCm� Co�fCr� Cu�Cw��Cz  C|� C  C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�L�Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��C�L�Cތ�C�� C��3C�33C� C���C��C�L�C� C�3C�  C�@ C� C�� C��3C�@ C� C�� C�  C�@ C�� C�� C�  C�� C��D �fDfD@ D� D�fD  D@ D	� D
� D  D@ D� D�fDfDFfD� D� D  D@ D�fD� D  D@ D�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0FfD1� D2� D4  D59�D6� D7� D9  D:@ D;� D<�fD>fD?@ D@� DA�fDCfDDFfDE�fDF� DG��DI@ DJ�fDK� DM  DN@ DO� DP� DQ��DS9�DTy�DU��DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  DlFfDm� Dn� Dp  Dq@ Dr� Ds� DufDv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D��D�� D�` D�  Dȣ3D�C3D��3Dʀ D�  D˼�D�\�D���D͠ D�@ D�� Dπ D�  D�� D�` D���DҠ D�C3D�� DԀ D��D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D���D� D�#3D��3D�c3D�3D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�#3D��3D�c3D�3D� D�@ D�� D�|�D��D��D�` D�3D�� D�@ D�� D�� D�  D�� D�\�D���D���D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@9��@�  @�  @�  @�  A��AffA   A0  AA��AP  A^ffAp  A�  A�  A�  A�  A�33A�  A���A���A�  A�  A�  A�  A�  A�  A�  A���B   B  BffB  B  B  BffB  B ffB$ffB(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\ffB`  Bd  Bg��Bl  Bp  Bt  BxffB|  B�  B�  B�  B�  B�33B�  B�  B�33B�  B�  B�  B�  B�  B���B�  B�33B�  B���B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33C  C� C�fC	ffC  C� C�C��C  C��C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C>  C@� CC  CE� CH  CJ� CM  CO� CR  CTffCW  CY��C\�C^��Ca  Cc� Cf�Ch� Cj�fCm� Co�fCr� Cu�Cw��Cz  C|� C  C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�33C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�L�Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��C�L�Cތ�C�� C��3C�33C� C���C��C�L�C� C�3C�  C�@ C� C�� C��3C�@ C� C�� C�  C�@ C�� C�� C�  C�� C��D �fDfD@ D� D�fD  D@ D	� D
� D  D@ D� D�fDfDFfD� D� D  D@ D�fD� D  D@ D�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0FfD1� D2� D4  D59�D6� D7� D9  D:@ D;� D<�fD>fD?@ D@� DA�fDCfDDFfDE�fDF� DG��DI@ DJ�fDK� DM  DN@ DO� DP� DQ��DS9�DTy�DU��DW  DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  DlFfDm� Dn� Dp  Dq@ Dr� Ds� DufDv@ Dw� Dx� Dz  D{@ D|� D}� D  D�  D�� D�` D�  D�� D�@ D�� D�� D�  D���D�\�D���D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D��D�� D�` D�  Dȣ3D�C3D��3Dʀ D�  D˼�D�\�D���D͠ D�@ D�� Dπ D�  D�� D�` D���DҠ D�C3D�� DԀ D��D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D���D� D�#3D��3D�c3D�3D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�#3D��3D�c3D�3D� D�@ D�� D�|�D��D��D�` D�3D�� D�@ D�� D�� D�  D�� D�\�D���D���D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��H@���@���@���@�@�@�@���@���@��y@��y@��H@��y@��H@��H@��y@���@���@���@�@�@�@�@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�o@�o@�o@�o@�o@��@��@��@��@��@��@��@��@��@�"�@�"�@�"�@��@��@��@��@�"�@�"�@�o@���@���@��H@��@��@��@��y@��@��y@��@��@���@�@���@���@�@�"�@�
=@�
=@�"�@�"�@�"�@��@�"�@��@���@���@��@��@��@��H@��@��@��@��@��@�ȴ@���@��R@��!@���@�n�@�M�@�=q@�-@�@��T@���@��^@�hs@�G�@�G�@��@�%@�A�@��@��@���@���@��@��@��@��@��@���@���@���@���@���@��P@��P@��P@���@��P@��P@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��@��F@��F@��F@��F@��F@��F@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@�ƨ@�ƨ@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@�ƨ@�ƨ@�ƨ@�ƨ@��w@��w@��w@��w@��w@��w@��F@��F@��F@��F@��@��F@��F@��F@��F@��F@��F@��F@��w@��w@��w@��F@��F@��F@��F@��F@��F@��F@��w@��w@��w@��F@��F@��F@��F@��F@��F@���@�|�@�t�@�33@���@�^5@���@��^@��h@�/@��j@��j@�z�@�1'@�  @��w@���@�C�@�"�@�
=@��H@��@���@��@�&�@���@���@�I�@��@�dZ@�
=@�~�@�5?@�n�@��w@��;@�  @� �@�1@���@���@���@�5?@���@�x�@�p�@���@�9X@���@�\)@�+@�@��!@�v�@���@�X@��@��/@�j@�1'@��m@�  @��
@��w@���@�t�@�+@��@���@�v�@��^@���@��/@���@��@�I�@���@���@�`B@���@�z�@��@�;d@��R@�V@�M�@��@��D@�z�@���@���@��@�o@��H@��!@�^5@��@��7@�p�@�X@���@��w@�C�@��!@��#@���@��@��@�j@�bN@�Z@��@�(�@�?}@��@��@���@���@�x�@�?}@�O�@�p�@�`B@�hs@��@�x�@�x�@�`B@�7L@�x�@�/@���@��u@�r�@�(�@��@��@��\@�M�@�@���@�O�@�Ĝ@���@���@�I�@�@\)@}��@|j@z�@zn�@x��@xQ�@w�@wl�@w�@v��@v$�@t��@t��@tj@t9X@t�@s��@sS�@s@r~�@rM�@r^5@st�@so@s"�@s"�@s��@vff@wl�@v�R@x��@y�#@z-@z=q@y��@y��@y�^@y�7@z~�@y��@zJ@z��@{�@{�m@|(�@|Z@|j@|j@|(�@{�
@{ƨ@{�F@{ƨ@{��@{o@{S�@{��@{@z�\@{@z��@yhs@y��@y��@y�#@y��@y��@y��@y��@y�^@y�@z�@z=q@z-@y��@x�9@wK�@vv�@u�@u?}@uV@t�@t�/@u�@u/@uO�@uO�@u@vff@vE�@vff@v5?@u�-@v$�@u�-@u��@up�@t�@t��@t��@tZ@t(�@s��@s��@s33@r��@r=q@q��@q�7@p��@q�@p�`@pr�@p �@o+@nE�@m�@m��@m�@l��@m?}@m`B@mO�@m�@l�/@l��@lZ@k��@kƨ@k�@kS�@k@j��@j�@i��@i7L@h��@hr�@hb@g��@g;d@g
=@f�@f�+@f$�@e��@ep�@e?}@d�@d�D@dz�@cƨ@ct�@c33@co@b��@b~�@a�^@a��@aG�@a%@`�9@`A�@_�@_�@_l�@_
=@_
=@^�y@^E�@^5?@^{@]�T@]�@]�@]O�@]/@\��@\��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��H@���@���@���@�@�@�@���@���@��y@��y@��H@��y@��H@��H@��y@���@���@���@�@�@�@�@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�o@�o@�o@�o@�o@��@��@��@��@��@��@��@��@��@�"�@�"�@�"�@��@��@��@��@�"�@�"�@�o@���@���@��H@��@��@��@��y@��@��y@��@��@���@�@���@���@�@�"�@�
=@�
=@�"�@�"�@�"�@��@�"�@��@���@���@��@��@��@��H@��@��@��@��@��@�ȴ@���@��R@��!@���@�n�@�M�@�=q@�-@�@��T@���@��^@�hs@�G�@�G�@��@�%@�A�@��@��@���@���@��@��@��@��@��@���@���@���@���@���@��P@��P@��P@���@��P@��P@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��@��F@��F@��F@��F@��F@��F@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@�ƨ@�ƨ@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@��w@�ƨ@�ƨ@�ƨ@�ƨ@��w@��w@��w@��w@��w@��w@��F@��F@��F@��F@��@��F@��F@��F@��F@��F@��F@��F@��w@��w@��w@��F@��F@��F@��F@��F@��F@��F@��w@��w@��w@��F@��F@��F@��F@��F@��F@���@�|�@�t�@�33@���@�^5@���@��^@��h@�/@��j@��j@�z�@�1'@�  @��w@���@�C�@�"�@�
=@��H@��@���@��@�&�@���@���@�I�@��@�dZ@�
=@�~�@�5?@�n�@��w@��;@�  @� �@�1@���@���@���@�5?@���@�x�@�p�@���@�9X@���@�\)@�+@�@��!@�v�@���@�X@��@��/@�j@�1'@��m@�  @��
@��w@���@�t�@�+@��@���@�v�@��^@���@��/@���@��@�I�@���@���@�`B@���@�z�@��@�;d@��R@�V@�M�@��@��D@�z�@���@���@��@�o@��H@��!@�^5@��@��7@�p�@�X@���@��w@�C�@��!@��#@���@��@��@�j@�bN@�Z@��@�(�@�?}@��@��@���@���@�x�@�?}@�O�@�p�@�`B@�hs@��@�x�@�x�@�`B@�7L@�x�@�/@���@��u@�r�@�(�@��@��@��\@�M�@�@���@�O�@�Ĝ@���@���@�I�@�@\)@}��@|j@z�@zn�@x��@xQ�@w�@wl�@w�@v��@v$�@t��@t��@tj@t9X@t�@s��@sS�@s@r~�@rM�@r^5@st�@so@s"�@s"�@s��@vff@wl�@v�R@x��@y�#@z-@z=q@y��@y��@y�^@y�7@z~�@y��@zJ@z��@{�@{�m@|(�@|Z@|j@|j@|(�@{�
@{ƨ@{�F@{ƨ@{��@{o@{S�@{��@{@z�\@{@z��@yhs@y��@y��@y�#@y��@y��@y��@y��@y�^@y�@z�@z=q@z-@y��@x�9@wK�@vv�@u�@u?}@uV@t�@t�/@u�@u/@uO�@uO�@u@vff@vE�@vff@v5?@u�-@v$�@u�-@u��@up�@t�@t��@t��@tZ@t(�@s��@s��@s33@r��@r=q@q��@q�7@p��@q�@p�`@pr�@p �@o+@nE�@m�@m��@m�@l��@m?}@m`B@mO�@m�@l�/@l��@lZ@k��@kƨ@k�@kS�@k@j��@j�@i��@i7L@h��@hr�@hb@g��@g;d@g
=@f�@f�+@f$�@e��@ep�@e?}@d�@d�D@dz�@cƨ@ct�@c33@co@b��@b~�@a�^@a��@aG�@a%@`�9@`A�@_�@_�@_l�@_
=@_
=@^�y@^E�@^5?@^{@]�T@]�@]�@]O�@]/@\��@\��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B6FB6FB6FB6FB6FB6FB5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B6FB5?B5?B5?B5?B5?B6FB6FB5?B5?B5?B5?B5?B6FB5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B49B49B5?B5?B5?B5?B5?B49B49B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B5?B49B5?B5?B49B49B49B49B49B49B33B33B2-B2-B2-B2-B1'B1'B1'B0!B0!B0!B/B/B.B.B.B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B,B,B,B+B+B)�B'�B%�B#�B#�B"�B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBhBbB\BVBJBJBVB{B�B�B�B�BoB+B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�sB�mB�fB�`B�TB�HB�HB�NB�NB�;B�#B�B��B��B��B��B��B��BɺBɺBȴBÖBĜBŢBƨBĜB��B��B��B��B�}B�wB�wB�qB�dB�LB�?B�3B�'B�B�B�B�B�B�B�B�B�9B�RB�RB�RB�LB�FB�FB�FB�RB�LB�RB�RB�XB�RB�RB�RB�XB�LB�FB�?B�?B�3B�'B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�uB�oB�hB�\B�\B�\B�\B�VB�VB�VB�PB�PB�PB�VB�oB�hB�oB�oB��B��B��B��B��B�B�B�B�B�B�B�B�3B�-B�9B�LB�XB�dB�jB�qB�wB�wB�wB�wB�wB�wB�wB�}B�}B��B��B��BBĜBĜBÖBƨBǮBǮBȴBȴBȴBȴBɺB��B��B��B��B��B��BȴBǮBƨBŢBŢBŢBǮBȴBȴBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBɺBȴBȴB��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBȴBȴBȴBǮBǮBǮBƨBƨBƨBƨBƨBŢBŢBŢBŢBŢBŢBŢBŢBŢBĜBĜBĜBĜBĜBĜBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBBBBBB11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B/�B/�B/�B/�B/�B/�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B/�B.�B.�B.�B.�B.�B/�B/�B.�B.�B.�B.�B.�B/�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B-�B-�B.�B.�B.�B.�B.�B-�B-�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B.�B-�B.�B.�B-�B-�B-�B-�B-�B-�B,�B,�B+�B+�B+�B+�B*�B*�B*�B)�B)�B)�B(�B(�B'�B'�B'�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B%�B%�B%�B$�B$�B#�B!}BpBdBdB^BXBKBKBEB?B9B9B4B.B.B(B(B(B!BBB	BB
�B	�B�B�B�B�B�B	BBBBB�B �B��B�wB�qB�qB�kB�YB�SB�LB�FB�FB�@B�:B�:B�.B�!B�B�B�B�B�B�B�B�B�
B�B�B��B��B��B��B��B��B��B��B��BԵBϖBΐB̄B�~B�qB�`B�TB�MB�MB�GB�)B�/B�5B�;B�/B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�xB�xB�rB�lB�gB�[B�HB�<B�)B�#B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�BB�TB�[B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�#B�0B�0B�*B�<B�BB�BB�HB�HB�HB�HB�NB�UB�[B�aB�gB�gB�UB�HB�BB�<B�6B�6B�6B�BB�HB�HB�NB�NB�[B�gB�gB�mB�mB�gB�sB�mB�mB�sB�mB�sB�sB�yB�yB�yB�sB�sB�sB�sB�sB�mB�mB�mB�mB�gB�aB�UB�NB�HB�NB�HB�HB�UB�UB�UB�UB�UB�UB�UB�UB�UB�NB�NB�NB�NB�HB�HB�HB�BB�BB�BB�<B�<B�<B�<B�<B�6B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�$B�$B�$B�$B�$B�$B�$11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0062856                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904482018110709044820181107090448  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172715  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172715  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090448  IP  PSAL            @ffD�� G�O�                