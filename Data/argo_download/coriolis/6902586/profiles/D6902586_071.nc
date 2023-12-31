CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T18:20:46Z creation; 2018-11-05T18:21:32Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105182046  20181107125409  6902586 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               GA   IF                                  2C  D   NOVA                            SN145                           n/a                             865 @�ĥR�1   @�ĥ�m�@N<�}l<h�CV��`t�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�33A   A  A!��A1��A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A���A�  A���A�  A�  A�  A���A�  A�  A�  B   B  B  B  B  BffB  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BC��BH  BL  BP  BT  BX  B\  B`  BdffBh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�33B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  Bə�B���B�  B�  B�  B�33B�  B�  B�  B�  B�33B�33C�C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C=�fC@� CC  CE� CH  CJ� CM  CO��CR  CTffCW  CY� C\  C^� C`�fCc� Cf  ChffCk  Cm� Cp  Cr��Cu  Cw� Cz�C|� C  C�� C��3C�33C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C��C�@ C�� C�� C�  C�@ Cŀ CƳ3C��3C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ Cڳ3C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��3C�33C�s3C�� C�  C�33C� C�� C�  C�@ C�� C�� C�  C��3C���C�  D � DfD@ D� D�fD  D@ D	� D
� D  D@ D� D� D  DFfD�fD� D  D@ D� D� D  D@ D� D��D��D!@ D"� D#� D%  D&@ D'� D(��D)��D+9�D,� D-� D/  D0@ D1� D2� D3��D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DI@ DJy�DK� DMfDN@ DOy�DP� DR  DS@ DT� DU�fDWfDX@ DY� DZ�fD\  D]9�D^� D_� Da  Db@ Dc� Dd� Df  DgFfDh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  DvFfDw� Dx� DzfD{@ D|� D}�fDfD�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�  D�� D�C3D�� D�� D�  D�� D�` D�3D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�<�D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  Dà D�@ D�� DŃ3D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D��D�� D�` D�  D͠ D�@ D�� Dπ D�  Dм�D�\�D�  DҠ D�C3D��3Dԃ3D�#3D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�c3D�3Dܣ3D�@ D�� Dހ D�  D�� D�` D�  D� D�<�D���D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D��D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D�3D�@ D�� D� D�  D�� D�` D�3D�� D�@ D�� D�|�D�  D�� D�` D�  D���D�	�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @ff@@  @�  @�  @�  @�33A   A  A!��A1��A@  AP  A`  Ap  A�  A�  A�  A�  A�  A�  A���A�  A���A�  A�  A�  A���A�  A�  A�  B   B  B  B  B  BffB  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BC��BH  BL  BP  BT  BX  B\  B`  BdffBh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�  B�  B�  B�  B���B���B�  B�33B�33B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  Bə�B���B�  B�  B�  B�33B�  B�  B�  B�  B�33B�33C�C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'� C*  C,� C/  C1� C4  C6� C9  C;� C=�fC@� CC  CE� CH  CJ� CM  CO��CR  CTffCW  CY� C\  C^� C`�fCc� Cf  ChffCk  Cm� Cp  Cr��Cu  Cw� Cz�C|� C  C�� C��3C�33C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C���C���C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�33C�� C���C��C�@ C�� C�� C�  C�@ Cŀ CƳ3C��3C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ Cڳ3C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��3C�33C�s3C�� C�  C�33C� C�� C�  C�@ C�� C�� C�  C��3C���C�  D � DfD@ D� D�fD  D@ D	� D
� D  D@ D� D� D  DFfD�fD� D  D@ D� D� D  D@ D� D��D��D!@ D"� D#� D%  D&@ D'� D(��D)��D+9�D,� D-� D/  D0@ D1� D2� D3��D5@ D6� D7� D8��D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DI@ DJy�DK� DMfDN@ DOy�DP� DR  DS@ DT� DU�fDWfDX@ DY� DZ�fD\  D]9�D^� D_� Da  Db@ Dc� Dd� Df  DgFfDh� Di� Dk  Dl@ Dm� Dn� Dp  Dq@ Dr� Ds� Du  DvFfDw� Dx� DzfD{@ D|� D}�fDfD�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�3D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D��3D�c3D�  D�� D�C3D�� D�� D�  D�� D�` D�3D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�<�D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�<�D���D�� D�  D�� D�` D�  Dà D�@ D�� DŃ3D�#3D��3D�` D�  DȠ D�@ D�� Dʀ D��D�� D�` D�  D͠ D�@ D�� Dπ D�  Dм�D�\�D�  DҠ D�C3D��3Dԃ3D�#3D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�c3D�3Dܣ3D�@ D�� Dހ D�  D�� D�` D�  D� D�<�D���D� D�  D�� D�` D�  D� D�@ D�� D�|�D�  D��D�` D�  D� D�@ D�� D� D�  D�� D�` D�3D�3D�@ D�� D� D�  D�� D�` D�3D�� D�@ D�� D�|�D�  D�� D�` D�  D���D�	�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A"-A"1'A"bA!�#A z�AG�A�HA�jA=qA1A�A�
A��AƨA�^A�A��A�hAp�AdZAXAS�AO�AG�A&�A9XA33@�9X@�G�@�@ٲ-@ԛ�@�1'@̼j@�@��@�33@�r�@�n�@�1'@���@�\)@�/@��D@�5?@�z�@��@���@���@���@��@���@�
=@�l�@���@�A�@�j@�A�@��;@�l�@�o@��T@��D@��@�33@���@���@�M�@���@���@�X@��@���@���@��@�E�@���@�/@��`@� �@�v�@�I�@���@�^5@�-@��@��7@�X@�G�@��^@��@��T@�x�@��@�1'@�(�@�  @�S�@��!@��!@�=q@�@��@�@�?}@�ƨ@�"�@�=q@�z�@���@���@��w@�~�@��@�X@��u@�1'@���@��@��h@�V@�Ĝ@��@��@�\)@�o@��@���@���@���@��+@�ff@�V@��\@�~�@�=q@�J@���@��@���@�r�@�b@�S�@���@���@�n�@�$�@��-@�hs@�7L@��@���@�Ĝ@���@�z�@�Q�@��;@���@�|�@�@���@���@�~�@�v�@�n�@�n�@�M�@�-@��#@��-@��7@�&�@���@��`@���@��9@��u@�z�@�Z@�A�@� �@�;@��@K�@~�y@~E�@}�@|z�@{�m@{��@{33@{"�@{dZ@{@y��@y%@xQ�@x �@w��@w�P@v�@v{@u�T@vff@v��@v�y@v��@vV@u�T@u�T@uO�@t��@t�D@tz�@t9X@s��@sƨ@sdZ@sC�@r�@r��@r~�@r-@q�#@q�7@q&�@p��@p�u@pr�@p1'@p1'@p1'@pb@o�;@o�P@oK�@n��@n��@nff@nE�@n$�@m�T@m�@mO�@mO�@mO�@mO�@l��@l�@l��@l��@lj@lI�@l1@k�m@kƨ@k�F@k��@k��@k�F@k�
@l1@k��@kƨ@k�F@k�F@kC�@j�H@j��@jn�@j�@i��@ix�@i�@hĜ@g|�@g+@g�@f�@f�@f�@f�y@fȴ@f�+@e�T@e��@ep�@d�@dj@d9X@d9X@d(�@c�
@c��@c��@c�@d1@dz�@dz�@d�@dj@dj@dI�@d�@cƨ@b�H@b��@b~�@b�\@bn�@b~�@b��@b��@b�\@b�@cdZ@c��@cƨ@c�m@c�
@c�F@c��@ct�@cC�@b�@b^5@b�@b�@a�@a7L@_�@^�y@^�R@_K�@`b@`bN@`��@`�u@a7L@a�@`�9@`Ĝ@`bN@`1'@^�y@^@]p�@\��@\�D@\�j@\�@\1@\�@\�@[�m@[��@[�@[ƨ@[ƨ@[�
@[��@[dZ@[C�@[��@[ƨ@[�m@[��@[�
@[S�@Z�H@Z��@Z��@Z��@Z�@Z�H@Z�!@Z�\@Z^5@Z~�@Z~�@Z�\@Zn�@ZM�@ZJ@Y��@Y�^@Y��@Y��@Y��@Y�7@YG�@Y7L@YG�@YG�@Y&�@X��@Xr�@W��@W�@W�;@W�;@W�w@W�;@W�;@W�w@W|�@W|�@W�P@W�@V�y@V�@Vȴ@V�@V�@V�@V�@V��@V��@V�+@Vv�@VE�@VE�@V{@U�@V@V{@V$�@V{@U�@U@U�-@U�-@U��@U��@U�h@U�@Up�@Up�@U�@U�@U�@U��@U�h@U�@U�@U�@U��@V@V@Vȴ@W��@XbN@X��@Y%@XĜ@XĜ@[�
@\�j@]��@^{@^��@_
=@_l�@_�;@aX@b�!@ct�@dZ@cƨ@d9X@eO�@e�-@f$�@fff@fv�@f��@f��@fȴ@fȴ@f�R@f�R@f�R@f�R@fv�@fE�@f5?@f5?@e�T@e`B@e�@d�@dz�@d(�@c�F@c�@cdZ@co@b~�@bn�@b�@a��@a�7@a&�@`��@`�9@`Q�@`b@`1'@_�@_|�@^��@^�y@^��@^V@]@]��@]`B@]�@]/@\�@\�j@\(�@\1@[��@[�F@[t�@[33@Z��@Z-@Z=q@Z�@Y��@Y��@Y&�@XbN@X�9@X�u@XA�@W��@V��@V@U��@U`B@T�/@T��@T��@S�m@S�@S��@S��@R��@R^5@R��@R^5@Q�#@R^5@Qx�@P�9@PQ�@PA�@P�u@P��@P��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 A"-A"1'A"bA!�#A z�AG�A�HA�jA=qA1A�A�
A��AƨA�^A�A��A�hAp�AdZAXAS�AO�AG�A&�A9XA33@�9X@�G�@�@ٲ-@ԛ�@�1'@̼j@�@��@�33@�r�@�n�@�1'@���@�\)@�/@��D@�5?@�z�@��@���@���@���@��@���@�
=@�l�@���@�A�@�j@�A�@��;@�l�@�o@��T@��D@��@�33@���@���@�M�@���@���@�X@��@���@���@��@�E�@���@�/@��`@� �@�v�@�I�@���@�^5@�-@��@��7@�X@�G�@��^@��@��T@�x�@��@�1'@�(�@�  @�S�@��!@��!@�=q@�@��@�@�?}@�ƨ@�"�@�=q@�z�@���@���@��w@�~�@��@�X@��u@�1'@���@��@��h@�V@�Ĝ@��@��@�\)@�o@��@���@���@���@��+@�ff@�V@��\@�~�@�=q@�J@���@��@���@�r�@�b@�S�@���@���@�n�@�$�@��-@�hs@�7L@��@���@�Ĝ@���@�z�@�Q�@��;@���@�|�@�@���@���@�~�@�v�@�n�@�n�@�M�@�-@��#@��-@��7@�&�@���@��`@���@��9@��u@�z�@�Z@�A�@� �@�;@��@K�@~�y@~E�@}�@|z�@{�m@{��@{33@{"�@{dZ@{@y��@y%@xQ�@x �@w��@w�P@v�@v{@u�T@vff@v��@v�y@v��@vV@u�T@u�T@uO�@t��@t�D@tz�@t9X@s��@sƨ@sdZ@sC�@r�@r��@r~�@r-@q�#@q�7@q&�@p��@p�u@pr�@p1'@p1'@p1'@pb@o�;@o�P@oK�@n��@n��@nff@nE�@n$�@m�T@m�@mO�@mO�@mO�@mO�@l��@l�@l��@l��@lj@lI�@l1@k�m@kƨ@k�F@k��@k��@k�F@k�
@l1@k��@kƨ@k�F@k�F@kC�@j�H@j��@jn�@j�@i��@ix�@i�@hĜ@g|�@g+@g�@f�@f�@f�@f�y@fȴ@f�+@e�T@e��@ep�@d�@dj@d9X@d9X@d(�@c�
@c��@c��@c�@d1@dz�@dz�@d�@dj@dj@dI�@d�@cƨ@b�H@b��@b~�@b�\@bn�@b~�@b��@b��@b�\@b�@cdZ@c��@cƨ@c�m@c�
@c�F@c��@ct�@cC�@b�@b^5@b�@b�@a�@a7L@_�@^�y@^�R@_K�@`b@`bN@`��@`�u@a7L@a�@`�9@`Ĝ@`bN@`1'@^�y@^@]p�@\��@\�D@\�j@\�@\1@\�@\�@[�m@[��@[�@[ƨ@[ƨ@[�
@[��@[dZ@[C�@[��@[ƨ@[�m@[��@[�
@[S�@Z�H@Z��@Z��@Z��@Z�@Z�H@Z�!@Z�\@Z^5@Z~�@Z~�@Z�\@Zn�@ZM�@ZJ@Y��@Y�^@Y��@Y��@Y��@Y�7@YG�@Y7L@YG�@YG�@Y&�@X��@Xr�@W��@W�@W�;@W�;@W�w@W�;@W�;@W�w@W|�@W|�@W�P@W�@V�y@V�@Vȴ@V�@V�@V�@V�@V��@V��@V�+@Vv�@VE�@VE�@V{@U�@V@V{@V$�@V{@U�@U@U�-@U�-@U��@U��@U�h@U�@Up�@Up�@U�@U�@U�@U��@U�h@U�@U�@U�@U��@V@V@Vȴ@W��@XbN@X��@Y%@XĜ@XĜ@[�
@\�j@]��@^{@^��@_
=@_l�@_�;@aX@b�!@ct�@dZ@cƨ@d9X@eO�@e�-@f$�@fff@fv�@f��@f��@fȴ@fȴ@f�R@f�R@f�R@f�R@fv�@fE�@f5?@f5?@e�T@e`B@e�@d�@dz�@d(�@c�F@c�@cdZ@co@b~�@bn�@b�@a��@a�7@a&�@`��@`�9@`Q�@`b@`1'@_�@_|�@^��@^�y@^��@^V@]@]��@]`B@]�@]/@\�@\�j@\(�@\1@[��@[�F@[t�@[33@Z��@Z-@Z=q@Z�@Y��@Y��@Y&�@XbN@X�9@X�u@XA�@W��@V��@V@U��@U`B@T�/@T��@T��@S�m@S�@S��@S��@R��@R^5@R��@R^5@Q�#@R^5@Qx�@P�9@PQ�@PA�@P�u@P��@P��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
�;B
�5B
�5B
�/B
�BB
�`B
�sB
�yB
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�`B
�#B
��B
��B/BbNBZB�B�7B�9B�!B�-B�qB�jB�XB�XB�3B�'B�B�B�9B��B�;B��B�yB�)BŢBȴB��B�B�`B��BB%B%B%BBBBB��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�sB�ZB�5B�)B�B�B�B�
B�
B�B�B�#B�B�B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��BƨBĜBŢBƨBÖB��B��B�}B�wB�}BƨBɺBɺBɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBǮBɺBɺBȴBȴBǮBƨBŢBĜBÖB��B��B��B�}B�}B�wB�wB�qB�qB�qB�qB�qB�jB�jB�dB�^B�^B�XB�RB�RB�RB�RB�RB�RB�LB�LB�FB�FB�?B�?B�9B�9B�9B�9B�3B�3B�3B�3B�-B�-B�'B�'B�!B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�oB�oB�hB�\B�VB�VB�VB�VB�VB�\B�VB�VB�PB�JB�JB�DB�=B�=B�=B�=B�=B�=B�=B�=B�DB�JB�JB�JB�JB�JB�JB�DB�=B�7B�1B�1B�1B�1B�1B�1B�1B�1B�1B�=B�DB�DB�DB�DB�DB�DB�DB�=B�=B�7B�1B�1B�+B�%B�B�B�B�B�B�B�B�%B�+B�+B�%B�%B�B�B�B� B}�B|�B|�B|�B|�Bz�B{�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�By�By�Bz�B{�B{�B{�Bz�By�By�By�Bx�Bx�By�By�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bv�Bu�Bu�Bu�Bu�Bu�Bt�Bs�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bq�Bq�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bo�Bo�Bo�Bo�Bn�Bo�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bq�Bt�Bv�Bw�Bw�Bw�Bx�B�B�B�%B�1B�=B�DB�JB�VB�uB��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�!B�'B�'B�'B�-B�-B�-B�-B�-B�-B�-B�-B�-B�-B�3B�-B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�9B�3B�3B�-B�3B�3B�3B�3B�3B�3B�9B�9B�9B�9B�3B�9B�9B�?B�?B�?B�9B�9B�?B�?B�?B�?B�?B�9B�?B�?B�9B�9B�3B�-B�-B�-B�'B�-B�-B�'B�!B�-B�-B�!B�!B�-B�-B�-B�3B�-B�!B�!B�'B�9B�?B�?11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B
�;B
�5B
�5B
�/B
�BB
�`B
�sB
�yB
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�B
�`B
�#B
��B
��B/BbNBZB�B�7B�9B�!B�-B�qB�jB�XB�XB�3B�'B�B�B�9B��B�;B��B�yB�)BŢBȴB��B�B�`B��BB%B%B%BBBBB��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�sB�ZB�5B�)B�B�B�B�
B�
B�B�B�#B�B�B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��BƨBĜBŢBƨBÖB��B��B�}B�wB�}BƨBɺBɺBɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBǮBɺBɺBȴBȴBǮBƨBŢBĜBÖB��B��B��B�}B�}B�wB�wB�qB�qB�qB�qB�qB�jB�jB�dB�^B�^B�XB�RB�RB�RB�RB�RB�RB�LB�LB�FB�FB�?B�?B�9B�9B�9B�9B�3B�3B�3B�3B�-B�-B�'B�'B�!B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�uB�oB�oB�hB�\B�VB�VB�VB�VB�VB�\B�VB�VB�PB�JB�JB�DB�=B�=B�=B�=B�=B�=B�=B�=B�DB�JB�JB�JB�JB�JB�JB�DB�=B�7B�1B�1B�1B�1B�1B�1B�1B�1B�1B�=B�DB�DB�DB�DB�DB�DB�DB�=B�=B�7B�1B�1B�+B�%B�B�B�B�B�B�B�B�%B�+B�+B�%B�%B�B�B�B� B}�B|�B|�B|�B|�Bz�B{�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�By�By�Bz�B{�B{�B{�Bz�By�By�By�Bx�Bx�By�By�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bv�Bu�Bu�Bu�Bu�Bu�Bt�Bs�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bq�Bq�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bo�Bo�Bo�Bo�Bn�Bo�Bo�Bo�Bo�Bn�Bn�Bn�Bn�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bn�Bo�Bo�Bo�Bq�Bt�Bv�Bw�Bw�Bw�Bx�B�B�B�%B�1B�=B�DB�JB�VB�uB��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�!B�'B�'B�'B�-B�-B�-B�-B�-B�-B�-B�-B�-B�-B�3B�-B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�9B�3B�3B�-B�3B�3B�3B�3B�3B�3B�9B�9B�9B�9B�3B�9B�9B�?B�?B�?B�9B�9B�?B�?B�?B�?B�?B�9B�?B�?B�9B�9B�3B�-B�-B�-B�'B�-B�-B�'B�!B�-B�-B�!B�!B�-B�-B�-B�3B�-B�!B�!B�'B�9B�?B�?11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected . OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                             201811071254102018110712541020181107125410  IF  ARFMCODA024c                                                                20181105182046                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105182132  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105182132  QCF$                G�O�G�O�G�O�0000000000002040GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107125410  IP  PSAL            @ffD�	�G�O�                