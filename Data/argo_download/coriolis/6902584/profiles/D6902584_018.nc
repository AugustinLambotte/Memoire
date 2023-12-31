CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:08Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090440  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�&�r� �1   @�&�/�`y@N���B�A�B5����K2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @�  @�  @�33@�  @���A  A   A0  A@  AQ��Aa��Ap  A�  A�  A�  A�33A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  BffBffB  B��B��B  B   B$  B(  B,  B0  B4  B8  B<  B@  BDffBHffBLffBPffBTffBXffB\ffB`ffBdffBhffBlffBpffBt  Bx  B|ffB�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B���B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	ffC�fC� C  C� C  CffC  C� C �C"� C$�fC'� C*�C,� C/  C1� C4  C6� C9�C;� C=�fC@� CC�CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cy�fC|� C�C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C���C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C��C�L�C���C���C�  C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C��C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C��C�L�C� C�3C�  C�@ C� C�� C�  C�@ C�s3C�� C�  C�@ C� C���C�  C�@ C�� C�� C�  C�� C�  D � DfD@ D� D� D  D@ D	� D
� D  DFfD� D� D  DFfD� D� D  D@ D� D� D  DFfD�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,�fD-�fD/  D0@ D1� D2� D4  D5@ D6�fD7� D9  D:@ D;y�D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DOy�DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� D`��Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn��Do��Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D��D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�C3D��3D��3D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D���D�� D�  D�� D�` D�  D�� D�<�D�� D��3D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D��3D�  D�� D�` D�  Dà D�@ D�� D�|�D�  D�� D�` D�  DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�  D�� D�` D�  DҠ D�@ D�� DԀ D��D�� D�` D���Dל�D�@ D�� Dـ D�#3D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D�3D�  D�� D�` D�  D�3D�@ D�� D� D�  D��3D�` D�  D� D�<�D�� D� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�` D�  D��3D�@ D�� D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@@  @�  @�  @�33@�  @���A  A   A0  A@  AQ��Aa��Ap  A�  A�  A�  A�33A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  BffBffB  B��B��B  B   B$  B(  B,  B0  B4  B8  B<  B@  BDffBHffBLffBPffBTffBXffB\ffB`ffBdffBhffBlffBpffBt  Bx  B|ffB�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B���B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  C  C� C  C	ffC�fC� C  C� C  CffC  C� C �C"� C$�fC'� C*�C,� C/  C1� C4  C6� C9�C;� C=�fC@� CC�CE� CH  CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cy�fC|� C�C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�33C�� C�� C�  C�@ C���C�� C�  C�L�C�� C�� C�  C�@ C�� C�� C��C�L�C���C���C�  C�L�C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ C�s3C�� C��C�@ CԀ C�� C�  C�@ Cـ C�� C�  C�@ Cހ C�� C��C�L�C� C�3C�  C�@ C� C�� C�  C�@ C�s3C�� C�  C�@ C� C���C�  C�@ C�� C�� C�  C�� C�  D � DfD@ D� D� D  D@ D	� D
� D  DFfD� D� D  DFfD� D� D  D@ D� D� D  DFfD�fD� D   D!@ D"� D#� D%  D&@ D'� D(� D*  D+@ D,�fD-�fD/  D0@ D1� D2� D4  D5@ D6�fD7� D9  D:@ D;y�D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DOy�DP� DR  DS@ DT� DU� DW  DX@ DY� DZ� D\  D]@ D^� D_� D`��Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn��Do��Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{FfD|� D}� D  D��D�� D�` D�  D�� D�C3D�� D��3D�  D�� D�` D�  D��3D�C3D��3D��3D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D���D�� D�  D�� D�` D�  D�� D�<�D�� D��3D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D��3D�  D�� D�` D�  Dà D�@ D�� D�|�D�  D�� D�` D�  DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dσ3D�  D�� D�` D�  DҠ D�@ D�� DԀ D��D�� D�` D���Dל�D�@ D�� Dـ D�#3D�� D�` D�  Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D�3D�  D�� D�` D�  D�3D�@ D�� D� D�  D��3D�` D�  D� D�<�D�� D� D�  D�� D�` D�  D�� D�@ D���D�� D�#3D�� D�` D�  D��3D�@ D�� D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�Ĝ@�%@�&�@�O�@�?}@��@��`@��`@���@�V@�G�@��h@���@���@���@���@���@��-@���@���@��7@��7@��7@��@�x�@��7@�x�@�`B@�O�@�7L@�/@�&�@�&�@�7L@�7L@�/@�/@�/@�&�@�&�@��@��@�&�@�7L@�7L@�7L@�7L@�7L@�7L@�/@�/@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�?}@�?}@�7L@�7L@�?}@�7L@�/@�/@�/@�/@�/@�/@�/@�/@�/@�&�@�V@�V@�V@�V@�%@�%@���@���@��@��@��@��@���@���@���@���@���@���@��@��/@��/@��/@��/@��/@��/@��/@��/@�A�@�@�-@���@�&�@� �@�x�@��@���@���@��P@�o@���@���@�Ĝ@��D@�bN@�b@�+@�@��H@���@��@��@���@���@�?}@��D@�r�@�j@��@�;d@�ȴ@��+@�v�@�n�@��@���@�p�@�G�@�&�@��`@�Ĝ@���@�Z@��;@���@�dZ@�+@�"�@��+@�M�@�{@��@��@��h@�&�@��/@�bN@��w@�;d@�+@�"�@�"�@�+@�
=@�v�@�E�@��@��#@��#@��@�p�@�r�@�  @��;@�|�@�;d@�33@�
=@�ȴ@���@���@�=q@�@��-@���@���@���@�`B@�Ĝ@��D@�9X@�1@��;@�ƨ@���@�ƨ@��
@��m@��F@���@��@��y@��+@�^5@�=q@�{@��#@���@��#@��#@��@��@�G�@���@���@���@�j@�I�@�A�@�1@��@��m@��w@���@�|�@�|�@�t�@�S�@�K�@�S�@�33@���@�ȴ@���@��@���@��R@���@�v�@�n�@�M�@�$�@�J@���@�@���@���@�p�@�`B@�X@�X@�G�@�/@��@���@���@��@�j@�bN@�I�@� �@��@~��@~��@~��@~V@}�T@}?}@|��@|��@|j@{�
@{t�@{"�@z�\@zM�@y�@y�#@y��@yx�@yG�@x��@y&�@x��@x1'@vȴ@v��@vff@u�T@u�-@u�@u@up�@t�D@tZ@s�F@r�!@r�@q��@qX@qx�@q��@qx�@q&�@p�@o�w@o�P@o��@o�@oK�@o
=@o�@n��@o
=@o�@nff@n$�@m�-@n{@m��@n$�@n�+@n�@nV@m�T@n�R@n�@o
=@n�R@nV@m�T@m@m@m��@m�-@m�-@m�@n��@o|�@p �@p��@p  @pr�@qhs@q��@r=q@r=q@r�\@rn�@r�!@s33@s33@sC�@sS�@sdZ@st�@sdZ@st�@s��@s��@s�F@s��@s��@so@r�H@r��@r^5@q�@q��@q��@q�^@q��@q�7@qhs@qX@qG�@q&�@q�@p��@p��@p�`@p��@pĜ@p�@pbN@pbN@p1'@o�@o�;@o�@o\)@o�@n��@n��@n�y@n�R@n�+@nv�@nV@n@m@m�@m`B@m?}@l��@l�j@lZ@k�
@k��@k"�@j�@j��@jM�@j�@i�^@ihs@i�@hĜ@h��@h�9@h�@h�@h�9@h  @g|�@g;d@g+@g+@g
=@fȴ@f�R@f�R@f��@fE�@e@e@e�-@ep�@e/@d�/@d�j@d�@d�@d�@eV@d�/@dj@dZ@c�F@cC�@b��@b��@b~�@b=q@b�@bJ@a�#@a�^@ahs@aG�@a&�@a%@`��@`��@`bN@_�;@_\)@^v�@]�@]@\��@\(�@[ƨ@[S�@Zn�@Y�^@YG�@Y7L@Y&�@X�`@X�`@X�`@X��@XbN@Xb@W��@W|�@W+@V�@Vȴ@W�w@XA�@XbN@XA�@X  @W��@W�@W�@W
=@V�y@Vv�@V5?@U�@U@Up�@U/@UO�@U?}@UO�@U?}@U/@T�@S��@Sƨ@Sƨ@S�F@S33@S@R�H@Rn�@RM�@R^5@R~�@R�!@R~�@R~�@RJ@Q��@Q&�@PbN@O��@N�y@Nȴ@N�+@M�@Mp�@M/@M�-@N@N�@N��@NV@N@M�-@Mp�@M�@L�j@Lj@L9X@L�@K��@Kƨ@K��@K�@K"�@Jn�@J�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@�Ĝ@�%@�&�@�O�@�?}@��@��`@��`@���@�V@�G�@��h@���@���@���@���@���@��-@���@���@��7@��7@��7@��@�x�@��7@�x�@�`B@�O�@�7L@�/@�&�@�&�@�7L@�7L@�/@�/@�/@�&�@�&�@��@��@�&�@�7L@�7L@�7L@�7L@�7L@�7L@�/@�/@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�7L@�?}@�?}@�7L@�7L@�?}@�7L@�/@�/@�/@�/@�/@�/@�/@�/@�/@�&�@�V@�V@�V@�V@�%@�%@���@���@��@��@��@��@���@���@���@���@���@���@��@��/@��/@��/@��/@��/@��/@��/@��/@�A�@�@�-@���@�&�@� �@�x�@��@���@���@��P@�o@���@���@�Ĝ@��D@�bN@�b@�+@�@��H@���@��@��@���@���@�?}@��D@�r�@�j@��@�;d@�ȴ@��+@�v�@�n�@��@���@�p�@�G�@�&�@��`@�Ĝ@���@�Z@��;@���@�dZ@�+@�"�@��+@�M�@�{@��@��@��h@�&�@��/@�bN@��w@�;d@�+@�"�@�"�@�+@�
=@�v�@�E�@��@��#@��#@��@�p�@�r�@�  @��;@�|�@�;d@�33@�
=@�ȴ@���@���@�=q@�@��-@���@���@���@�`B@�Ĝ@��D@�9X@�1@��;@�ƨ@���@�ƨ@��
@��m@��F@���@��@��y@��+@�^5@�=q@�{@��#@���@��#@��#@��@��@�G�@���@���@���@�j@�I�@�A�@�1@��@��m@��w@���@�|�@�|�@�t�@�S�@�K�@�S�@�33@���@�ȴ@���@��@���@��R@���@�v�@�n�@�M�@�$�@�J@���@�@���@���@�p�@�`B@�X@�X@�G�@�/@��@���@���@��@�j@�bN@�I�@� �@��@~��@~��@~��@~V@}�T@}?}@|��@|��@|j@{�
@{t�@{"�@z�\@zM�@y�@y�#@y��@yx�@yG�@x��@y&�@x��@x1'@vȴ@v��@vff@u�T@u�-@u�@u@up�@t�D@tZ@s�F@r�!@r�@q��@qX@qx�@q��@qx�@q&�@p�@o�w@o�P@o��@o�@oK�@o
=@o�@n��@o
=@o�@nff@n$�@m�-@n{@m��@n$�@n�+@n�@nV@m�T@n�R@n�@o
=@n�R@nV@m�T@m@m@m��@m�-@m�-@m�@n��@o|�@p �@p��@p  @pr�@qhs@q��@r=q@r=q@r�\@rn�@r�!@s33@s33@sC�@sS�@sdZ@st�@sdZ@st�@s��@s��@s�F@s��@s��@so@r�H@r��@r^5@q�@q��@q��@q�^@q��@q�7@qhs@qX@qG�@q&�@q�@p��@p��@p�`@p��@pĜ@p�@pbN@pbN@p1'@o�@o�;@o�@o\)@o�@n��@n��@n�y@n�R@n�+@nv�@nV@n@m@m�@m`B@m?}@l��@l�j@lZ@k�
@k��@k"�@j�@j��@jM�@j�@i�^@ihs@i�@hĜ@h��@h�9@h�@h�@h�9@h  @g|�@g;d@g+@g+@g
=@fȴ@f�R@f�R@f��@fE�@e@e@e�-@ep�@e/@d�/@d�j@d�@d�@d�@eV@d�/@dj@dZ@c�F@cC�@b��@b��@b~�@b=q@b�@bJ@a�#@a�^@ahs@aG�@a&�@a%@`��@`��@`bN@_�;@_\)@^v�@]�@]@\��@\(�@[ƨ@[S�@Zn�@Y�^@YG�@Y7L@Y&�@X�`@X�`@X�`@X��@XbN@Xb@W��@W|�@W+@V�@Vȴ@W�w@XA�@XbN@XA�@X  @W��@W�@W�@W
=@V�y@Vv�@V5?@U�@U@Up�@U/@UO�@U?}@UO�@U?}@U/@T�@S��@Sƨ@Sƨ@S�F@S33@S@R�H@Rn�@RM�@R^5@R~�@R�!@R~�@R~�@RJ@Q��@Q&�@PbN@O��@N�y@Nȴ@N�+@M�@Mp�@M/@M�-@N@N�@N��@NV@N@M�-@Mp�@M�@L�j@Lj@L9X@L�@K��@Kƨ@K��@K�@K"�@Jn�@J�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB� B�B�B�B�B�B� B� B� B� B�B�B�B�B�B�B�B�%B�B�B�%B�+B�+B�%B�%B�+B�%B�%B�+B�1B�DB�PB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB��B��B��B��B��B��B�'B�LB�XB�^B�dB�jB�qB�dB�dB�jB�jB�qB�jB�jB�jB�dB�^B�^B�dB�dB�dB�^B�^B�^B�XB�RB�LB�RB�RB�XB�^B�XB�XB�XB�XB�^B�dB�dB�dB�dB�dB�^B�^B�^B�XB�RB�RB�RB�RB�RB�LB�LB�FB�?B�9B�9B�9B�9B�9B�9B�3B�3B�-B�-B�-B�-B�'B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�oB�oB�oB�oB�oB�hB�VB�VB�VB�PB�PB�PB�VB�PB�JB�DB�=B�1B�+B�+B�%B�+B�+B�+B�+B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�+B�1B�+B�+B�7B�=B�=B�7B�7B�1B�1B�1B�1B�1B�7B�=B�PB�\B�hB�{B�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�!B�!B�'B�'B�'B�'B�-B�-B�-B�-B�3B�3B�3B�9B�9B�3B�3B�3B�3B�9B�9B�9B�9B�9B�9B�9B�?B�9B�9B�9B�9B�9B�9B�9B�?B�?B�?B�FB�FB�FB�FB�?B�?B�?B�?B�FB�FB�FB�FB�FB�FB�FB�FB�FB�FB�FB�RB�XB�^B�^B�jB�jB�dB�dB�dB�^B�^B�XB�XB�^B�^B�^B�^B�^B�dB�dB�dB�dB�dB�dB�dB�^B�XB�RB�LB�LB�FB�?B�?B�?B�9B�-B�-B�-B�3B�3B�9B�9B�9B�9B�9B�9B�9B�?B�?B�FB�XB�dB�jB�jB�jB�qB�jB�jB�jB�jB�jB�jB�qB�qB�qB�wB�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�wB�wB�wB�wB�}B��B��BBBÖBBB��B�}B�wB�qB�qB�wB�qB�jB�jB�}B��BĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBB��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B}]B~cB~cBiBiBiB}]B}]B}]B}]B~cBiB�vB�vB�|B�vB�|B��B�|B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�3B�?B�?B�RB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�xB�qB�qB�qB�kB�kB�kB�kB�kB�kB�kB�eB�eB�eB�kB�kB�eB�_B�YB�YB�YB�SB�YB�SB�SB�YB�_B�_B�YB�YB�MB�FB�@B�@B�@B�@B�@B�@B�@B�FB�FB�:B�4B�.B�.B�(B�.B�(B�(B�(B�(B�(B�(B�"B�"B�"B�"B�"B�"B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�wB�wB�}B�}B�}B�wB�}B�}B�}B�}B�}B�wB�wB�}B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�"B�(B�.B�.B�4B�:B�@B�@B�FB�MB�SB�SB�SB�SB�SB�YB�YB�YB�YB�_B�_B�_B�_B�_B�_B�eB�kB�kB�qB�xB�xB�xB�xB�~B�~B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0025752                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904402018110709044020181107090440  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172708  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172708  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090440  IP  PSAL            @��D��3G�O�                