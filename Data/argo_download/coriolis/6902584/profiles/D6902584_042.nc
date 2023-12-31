CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  	   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:18Z creation; 2018-11-05T17:27:21Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z        $  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        $  C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     $  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  V    TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ^D   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  `P   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ht   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  j�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  r�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   z�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  |�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     $  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �(   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �X   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �X   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �X   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �X             ,  �XArgo profile    3.1 1.2 19500101000000  20181105172618  20181107090453  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               *A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�b�˄
�1   @�b�N��@O�j!
�C[F��@1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  @���A  A!��A0  A@  AP  A`  Aq��A�  A�33A�  A�  A�33A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�33B   B  B��B  B  B  B  B  B   B#��B(  B,  B/��B4  B8ffB<  B@  BD  BH  BLffBP  BT  BXffB\  B`  BdffBh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  B�  B�  B�33B�33B�  B���B�  B�33B�  B�  B�33B�  B�  B�  B�33B�  B���B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C  C� C�C� C   C"� C%  C'ffC*  C,� C/  C1� C3�fC6ffC9  C;� C>  C@� CC�CE� CG�fCJ� CM�CO��CR�CT� CW  CY��C\  C^ffCa  Cc��Cf  ChffCk  Cm� Cp  Cr��Cu  Cw� Cz�C|� C  C�� C��C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԌ�C�� C�  C�@ C�s3Cڳ3C�  C�@ Cހ C�� C�  C�@ C�s3C�3C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C�� C��D �fD  D@ D� D� DfD@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D��D@ D� D� D   D!@ D"� D#�fD%fD&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>fD?@ D@y�DA� DC  DDFfDE�fDF� DH  DI@ DJ� DK��DM  DN@ DO� DP� DR  DS@ DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^� D_� DafDb@ Dc�fDd� Df  Dg@ Dh� Di�fDkfDl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�#3D�� D�` D�3D��3D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D���D�� D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D���D�@ D�� D��3D�  D���D�` D�3D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�\�D�  D�� D�<�D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�<�D���Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�#3D�� D�c3D�3DҠ D�<�D�� DԀ D��Dռ�D�` D�3Dנ D�@ D�� Dـ D�  D��3D�c3D�  Dܠ D�@ D�� D�|�D�#3D��3D�c3D�  D��fD�f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @ff@@  @�  @�  @�  @�  @���A  A!��A0  A@  AP  A`  Aq��A�  A�33A�  A�  A�33A�  A�  A�  A�  A�  A���A�  A�  A�  A�  A�33B   B  B��B  B  B  B  B  B   B#��B(  B,  B/��B4  B8ffB<  B@  BD  BH  BLffBP  BT  BXffB\  B`  BdffBh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  B�  B�  B�33B�33B�  B���B�  B�33B�  B�  B�33B�  B�  B�  B�33B�  B���B�  B�33B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  C  C� C  C	� C  C� C  C� C  C� C�C� C   C"� C%  C'ffC*  C,� C/  C1� C3�fC6ffC9  C;� C>  C@� CC�CE� CG�fCJ� CM�CO��CR�CT� CW  CY��C\  C^ffCa  Cc��Cf  ChffCk  Cm� Cp  Cr��Cu  Cw� Cz�C|� C  C�� C��C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�s3C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԌ�C�� C�  C�@ C�s3Cڳ3C�  C�@ Cހ C�� C�  C�@ C�s3C�3C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�L�C�� C�� C�  C�� C��D �fD  D@ D� D� DfD@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D� D��D@ D� D� D   D!@ D"� D#�fD%fD&@ D'� D(� D*  D+FfD,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>fD?@ D@y�DA� DC  DDFfDE�fDF� DH  DI@ DJ� DK��DM  DN@ DO� DP� DR  DS@ DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^� D_� DafDb@ Dc�fDd� Df  Dg@ Dh� Di�fDkfDl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|y�D}� D  D�#3D�� D�` D�3D��3D�C3D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�<�D���D�� D�  D���D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D���D�� D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D���D�@ D�� D��3D�  D���D�` D�3D�� D�<�D���D�� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�\�D�  D�� D�<�D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�<�D���Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�#3D�� D�c3D�3DҠ D�<�D�� DԀ D��Dռ�D�` D�3Dנ D�@ D�� Dـ D�  D��3D�c3D�  Dܠ D�@ D�� D�|�D�#3D��3D�c3D�  D��fD�f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�@�D@�D@�D@��u@���@���@���@���@���@���@���@��u@�D@�D@�D@��u@��u@�Q�@�1'@� �@�b@�  @��
@��;@��m@��@��m@��;@�F@�dZ@�dZ@�;d@�+@�"�@�
=@��H@��@�O�@�&�@�V@�V@��@�V@��@�/@��
@�"�@�
=@���@@��#@�`B@�O�@�7L@��@���@��
@�K�@�33@�"�@�
=@�!@�~�@�ff@�ff@�V@�ff@�v�@�^5@�$�@��@�J@��@��@��@�@�j@�j@�A�@�A�@�A�@�A�@�1'@���@�ƨ@�ƨ@�b@�@�$�@�J@�-@�-@�~�@�\@�\@旍@�R@�R@��@�9X@�b@�R@��@߮@�
=@ް!@�E�@�@ݩ�@�x�@�X@�7L@�%@ܬ@�b@�V@�9X@�1'@��@��;@��;@�ƨ@׶F@ְ!@Չ7@���@���@�Ĝ@�bN@Ӆ@�+@ҟ�@ѡ�@���@Ь@��
@�ff@�9X@�\)@ʰ!@�=q@�p�@ȓu@�A�@Ǯ@��@š�@ēu@�z�@��;@�K�@+@�5?@���@��^@��h@�x�@�?}@��`@�b@�o@�-@��@�?}@�r�@��m@�C�@�@�@�`B@�%@�Z@�|�@�"�@��!@�V@�ff@�V@�5?@�@���@�V@��F@�dZ@�+@�C�@�^5@���@�(�@��m@��@�33@��H@��@��y@��@��+@�-@���@�@�O�@�/@��@���@�(�@�E�@���@�z�@�Z@�1'@���@��@�C�@��@��@���@�E�@���@��j@��j@���@��@�r�@�Z@�9X@�(�@�b@���@���@�  @�1@���@��@�33@��!@�{@�@�p�@��-@��h@�/@�z�@��D@��j@�S�@�V@��@�`B@�A�@��w@��@���@���@��!@��;@�{@��^@��@�n�@��@�\)@�K�@�C�@�C�@�o@�v�@���@��`@��u@���@�7L@�p�@��h@���@�O�@���@�9X@�b@�ƨ@���@��P@�dZ@�dZ@���@��@��y@���@���@�J@�bN@�t�@��!@�-@��@��F@���@���@��@��@��@�t�@�dZ@�S�@��@���@�5?@�O�@�&�@���@��/@���@��/@�Ĝ@���@��`@�G�@�V@��y@���@���@��@��y@�V@���@��T@��#@���@��#@���@���@���@��-@���@�/@���@�z�@�1'@�;@��@l�@+@+@+@~v�@~5?@}��@}�T@}O�@}�@|�@|��@|Z@|(�@|�/@}V@|�@|I�@{��@{��@{S�@{dZ@{��@|(�@{��@}�@~��@�w@~��@�w@~��@}��@|z�@{t�@|I�@z��@|�@}�-@}��@}O�@|I�@|��@|9X@{��@|I�@|��@|�@|Z@|j@|j@|�D@|j@{�@z^5@z-@z�@y�^@y�^@y��@y�7@yx�@yX@y%@xQ�@wl�@w|�@w�P@v��@u�T@v{@vff@v��@vff@u?}@uO�@t��@s��@sdZ@r��@q�^@pA�@n��@nV@m��@mV@lz�@l(�@k��@j�H@jn�@h��@g�@e�h@d��@d�@d�/@dz�@b��@` �@^�+@\��@\�@[�@["�@[@[@Z�@Z�H@Y��@X��@X�`@XĜ@X��@X�@W�w@V��@V@U?}@Tj@R�\@P��@P1'@O�;@O�;@N�y@N$�@M�-@M`B@M�@L��@K�
@K�@K�
@K�m@L(�@Lj@L(�@K�
@KdZ@Ko@Jn�@I��@H�`@G��@G+@F��@Fȴ@F$�@E��@Ep�@D�@C��@A��@>��@>ȴ@>��@>$�@=�@<�@<��@<j@;ƨ@:M�@4�@1�7@.�y@(1'@"^5@ �@K�@
=@��@
=@
=@
=@�@�@+@;d@
=@�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�@�D@�D@�D@��u@���@���@���@���@���@���@���@��u@�D@�D@�D@��u@��u@�Q�@�1'@� �@�b@�  @��
@��;@��m@��@��m@��;@�F@�dZ@�dZ@�;d@�+@�"�@�
=@��H@��@�O�@�&�@�V@�V@��@�V@��@�/@��
@�"�@�
=@���@@��#@�`B@�O�@�7L@��@���@��
@�K�@�33@�"�@�
=@�!@�~�@�ff@�ff@�V@�ff@�v�@�^5@�$�@��@�J@��@��@��@�@�j@�j@�A�@�A�@�A�@�A�@�1'@���@�ƨ@�ƨ@�b@�@�$�@�J@�-@�-@�~�@�\@�\@旍@�R@�R@��@�9X@�b@�R@��@߮@�
=@ް!@�E�@�@ݩ�@�x�@�X@�7L@�%@ܬ@�b@�V@�9X@�1'@��@��;@��;@�ƨ@׶F@ְ!@Չ7@���@���@�Ĝ@�bN@Ӆ@�+@ҟ�@ѡ�@���@Ь@��
@�ff@�9X@�\)@ʰ!@�=q@�p�@ȓu@�A�@Ǯ@��@š�@ēu@�z�@��;@�K�@+@�5?@���@��^@��h@�x�@�?}@��`@�b@�o@�-@��@�?}@�r�@��m@�C�@�@�@�`B@�%@�Z@�|�@�"�@��!@�V@�ff@�V@�5?@�@���@�V@��F@�dZ@�+@�C�@�^5@���@�(�@��m@��@�33@��H@��@��y@��@��+@�-@���@�@�O�@�/@��@���@�(�@�E�@���@�z�@�Z@�1'@���@��@�C�@��@��@���@�E�@���@��j@��j@���@��@�r�@�Z@�9X@�(�@�b@���@���@�  @�1@���@��@�33@��!@�{@�@�p�@��-@��h@�/@�z�@��D@��j@�S�@�V@��@�`B@�A�@��w@��@���@���@��!@��;@�{@��^@��@�n�@��@�\)@�K�@�C�@�C�@�o@�v�@���@��`@��u@���@�7L@�p�@��h@���@�O�@���@�9X@�b@�ƨ@���@��P@�dZ@�dZ@���@��@��y@���@���@�J@�bN@�t�@��!@�-@��@��F@���@���@��@��@��@�t�@�dZ@�S�@��@���@�5?@�O�@�&�@���@��/@���@��/@�Ĝ@���@��`@�G�@�V@��y@���@���@��@��y@�V@���@��T@��#@���@��#@���@���@���@��-@���@�/@���@�z�@�1'@�;@��@l�@+@+@+@~v�@~5?@}��@}�T@}O�@}�@|�@|��@|Z@|(�@|�/@}V@|�@|I�@{��@{��@{S�@{dZ@{��@|(�@{��@}�@~��@�w@~��@�w@~��@}��@|z�@{t�@|I�@z��@|�@}�-@}��@}O�@|I�@|��@|9X@{��@|I�@|��@|�@|Z@|j@|j@|�D@|j@{�@z^5@z-@z�@y�^@y�^@y��@y�7@yx�@yX@y%@xQ�@wl�@w|�@w�P@v��@u�T@v{@vff@v��@vff@u?}@uO�@t��@s��@sdZ@r��@q�^@pA�@n��@nV@m��@mV@lz�@l(�@k��@j�H@jn�@h��@g�@e�h@d��@d�@d�/@dz�@b��@` �@^�+@\��@\�@[�@["�@[@[@Z�@Z�H@Y��@X��@X�`@XĜ@X��@X�@W�w@V��@V@U?}@Tj@R�\@P��@P1'@O�;@O�;@N�y@N$�@M�-@M`B@M�@L��@K�
@K�@K�
@K�m@L(�@Lj@L(�@K�
@KdZ@Ko@Jn�@I��@H�`@G��@G+@F��@Fȴ@F$�@E��@Ep�@D�@C��@A��@>��@>ȴ@>��@>$�@=�@<�@<��@<j@;ƨ@:M�@4�@1�7@.�y@(1'@"^5@ �@K�@
=@��@
=@
=@
=@�@�@+@;d@
=@�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBhsBhsBhsBhsBhsBhsBgmBgmBgmBhsBgmBgmBhsBhsBhsBgmBgmBgmBiyBiyBhsBiyBiyBiyBiyBiyBiyBiyBhsBiyBiyBiyBiyBiyBhsBhsBiyBl�Bk�Bl�Bl�Bl�Bk�Bk�BjBiyBl�Bk�Bk�Bk�Bk�Bm�Bk�Bk�BjBhsBgmBgmBffBe`Be`BffBgmBhsBk�BjBjBjBjBl�Bl�Bm�Bm�Bl�Bn�Bp�Bo�Bp�Bo�Bo�Bo�Bo�Bn�Bn�Bo�Bn�Bo�Bl�Bn�Bq�Bp�Bp�Bq�Bo�Bo�Bo�Bn�Bm�Bm�Bp�Bo�Bn�Bn�Bp�Bt�Bt�Bt�Bt�Bt�Bu�Bu�Bu�Bt�Bt�Bs�Br�Bt�Bv�Bv�Bw�Bw�Bv�Bw�Bv�Bx�By�By�By�Bx�Bx�By�Bx�Bx�Bw�Bw�Bv�Bu�Bs�Bs�Br�Br�Bq�Bq�Bp�Bo�Bn�Bm�Bk�BiyBiyBhsBffBe`BdZBdZBcTBbNBbNBbNB`BB^5B\)BYBXBW
BVBS�BQ�BP�BO�BM�BL�BJ�BI�BH�BG�BF�BF�BF�BE�BE�BC�BB�B?}B>wB=qB=qB:^B7LB5?B5?B49B33B2-B2-B2-B1'B1'B0!B/B.B.B-B,B+B(�B$�B!�B �B�B�B�B�B�B�B�B�B�B�B{B{B{B{BuBuBuBuBoBoBoBoBoBhBhBbBVBJBJBDBJBDB
=B	7B1B+BB  B  B��B��B��B��B��B�B�B�mB�HB�HB�NB�`B�sB�B�B�B�B�B�yB�mB�TB�TB�`B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�sB�sB�sB�mB�ZB�BB�/B�#B�B�B��B��B��BȴBŢBĜBÖBÖBÖBB��B�}B�qB�jB�jB�jB�jB�qB�qB�qB�wB��B�wB�RB�RB�LB�LB�LB�?B�3B�3B�3B�3B�3B�-B�-B�3B�3B�3B�-B�'B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�'B�-B�-B�'B�'B�'B�-B�-B�3B�9B�?B�jB�wB��B��B��B�}B�qB�XB�XB�dB�RB�wBB��B��B�}B��B��B��BBÖBÖBĜBĜBŢBƨBƨBƨBǮBǮBǮBǮBǮBǮBȴBǮBǮBǮBǮBǮBǮBƨBƨBŢBǮBȴBɺBɺBɺBɺBȴBǮBǮBǮBƨBŢBĜBÖBÖBÖBÖBBBB��B��B�}B��B��B�}B�}B�wB�wB�qB�jB�jB�jB�jB�jB�jB�jB�jB�dB�jB�qB�qB�qB�jB�jB�jB�jB�jB�dB�dB�XB�LB�FB�?B�?B�9B�3B�-B�'B�'B�!B�!B�'B�-B�-B�3B�?B�?B�?B�?B�9B�3B�3B�-B�'B�'B�'B�!B�B�B�B�B�B�B�B�!B�!B�!B�'B�-B�'B�'B�!B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B^�B^�B^�B^�B^�B^�B]�B]�B]�B^�B]�B]�B^�B^�B^�B]�B]�B]�B_�B_�B^�B_�B_�B_�B_�B_�B_�B_�B^�B_�B_�B_�B_�B_�B^�B^�B_�Bb�Ba�Bb�Bb�Bb�Ba�Ba�B`�B_�Bb�Ba�Ba�Ba�Ba�Bc�Ba�Ba�B`�B^�B]�B]�B\�B[�B[�B\�B]�B^�Ba�B`�B`�B`�B`�Bb�Bb�Bc�Bc�Bb�Bd�Bf�Be�Bf�Be�Be�Be�Be�Bd�Bd�Be�Bd�Be�Bb�Bd�Bg�Bf�Bf�Bg�Be�Be�Be�Bd�Bc�Bc�Bf�Be�Bd�Bd�Bf�Bk	Bk	Bk	Bk	Bk	BlBlBlBk	Bk	BjBh�Bk	BmBmBnBnBmBnBmBo"Bp(Bp(Bp(Bo"Bo"Bp(Bo"Bo#BnBnBmBlBjBjBh�Bh�Bg�Bg�Bf�Be�Bd�Bc�Ba�B_�B_�B^�B\�B[�BZ�BZ�BY�BX�BX�BX�BV�BT�BRzBOhBNaBM[BLUBJJBH>BG7BF1BD%BCBAB@B?B>B<�B<�B<�B;�B;�B9�B8�B5�B4�B3�B3�B0�B-�B+�B+�B*�B)�B(�B(�B(�B'|B'|B&vB%pB$iB$iB#cB"]B!WBLB3B!BBBBB	BB�B�B�B�B�B
�B
�B
�B
�B	�B	�B	�B	�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B��B��B��B�kB�ZB�ZB�HB�/B�#B�B�B�B��B��BפBפBتBۼB��B��B��B��B��B��B��B��BٰBٰBۼB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BڷB֟BӌBрB�uB�bB�PB�>B�,B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�vB�vB�vB�vB�vB�vB�pB�pB�pB�pB�pB�jB�jB�pB�pB�vB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�B�B�B�B�B�B�B�B�	B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�xB�rB�mB�gB�sB��B��B��B��B��B��B��B��B�yB�gB�aB�OB�CB�CB�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>B�>11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0094263                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904542018110709045420181107090454  IF  ARFMCODA024c                                                                20181105172618                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172721  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172721  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090454  IP  PSAL            @ffD�fG�O�                