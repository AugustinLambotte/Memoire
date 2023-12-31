CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:16Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090448  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL                A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�I�|�h1   @�I�|�h@O�w�4��@��47�8   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A���A�  A�33A�33A�33A�  A���A�  A�  A�  A���A���B ffBffB  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bw��B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�33C�C� C  C	� C  C��C  C� C�fC� C  C� C �C"��C%�C'� C*  C,� C.�fC1� C4�C6� C9  C;� C>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW  CYffC\  C^��Ca  Cc��Cf  ChffCk  Cm� Cp  Cr� Cu  Cw��Cz  C|ffC~�fC�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C��3C�  C�L�C�� C�� C�  C�33C�� C�� C�  C�@ C���C�� C�  C�33C�� C���C�  C�@ Cŀ C�� C�  C�@ Cʌ�C���C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�@ Cٌ�C�� C��3C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D �fD  D@ D�fD� D  D@ D	� D
��D  DFfD� D� D  D@ D� D�fDfD@ Dy�D� D  D@ D� D�fD   D!9�D"y�D#��D$��D&9�D'y�D(� D*fD+FfD,� D-�fD/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DIFfDJ�fDK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY� DZ� D\  D]@ D^�fD_�fDa  Db@ Dcy�Dd� Df  Dg@ Dh� Di� DkfDlFfDm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|�fD}�fDfD�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D��3D��3D�#3D��3D�c3D�  D���D�@ D�� D�|�D�  D�� D�` D�  Dà D�<�D�� Dŀ D�  D��3D�c3D�3DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�<�D�� Dπ D�  D�� D�` D�  Dң3D�C3D�� DԀ D��D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�3Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�#3D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D���D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@@  @�  @�  @�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A���A���A�  A�33A�33A�33A�  A���A�  A�  A�  A���A���B ffBffB  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bw��B|  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B���B�  B�33C�C� C  C	� C  C��C  C� C�fC� C  C� C �C"��C%�C'� C*  C,� C.�fC1� C4�C6� C9  C;� C>  C@� CC  CE� CH  CJ��CM  CO� CR  CT� CW  CYffC\  C^��Ca  Cc��Cf  ChffCk  Cm� Cp  Cr� Cu  Cw��Cz  C|ffC~�fC�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C�� C��C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C��3C�  C�L�C�� C�� C�  C�33C�� C�� C�  C�@ C���C�� C�  C�33C�� C���C�  C�@ Cŀ C�� C�  C�@ Cʌ�C���C�  C�@ Cπ C�� C��C�@ CԀ C�� C�  C�@ Cٌ�C�� C��3C�@ Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C�  C�� C�  D �fD  D@ D�fD� D  D@ D	� D
��D  DFfD� D� D  D@ D� D�fDfD@ Dy�D� D  D@ D� D�fD   D!9�D"y�D#��D$��D&9�D'y�D(� D*fD+FfD,� D-�fD/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DIFfDJ�fDK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY� DZ� D\  D]@ D^�fD_�fDa  Db@ Dcy�Dd� Df  Dg@ Dh� Di� DkfDlFfDm� Dn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx� Dz  D{@ D|�fD}�fDfD�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D�� D�@ D�� D�� D��D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D���D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�#3D�� D�` D�  D�� D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�c3D�  D�� D�@ D��3D��3D�#3D��3D�c3D�  D���D�@ D�� D�|�D�  D�� D�` D�  Dà D�<�D�� Dŀ D�  D��3D�c3D�3DȠ D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�<�D�� Dπ D�  D�� D�` D�  Dң3D�C3D�� DԀ D��D�� D�` D�  Dנ D�@ D�� Dـ D�  D�� D�` D�3Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D��3D� D�#3D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D���D�<�D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��H@��H@��@��R@��R@��@��y@��y@���@���@�V@�M�@�=q@�5?@��@��@�$�@�$�@�{@��@�{@�J@�J@�J@�@��@���@��@���@���@��^@�x�@�`B@�O�@��@�O�@�V@�%@��@��;@��w@��F@���@�K�@���@�-@�E�@�V@�M�@���@�@��^@���@��-@��-@��^@��-@���@��-@��^@��-@��-@��@�G�@�&�@���@�Ĝ@�Ĝ@�Ĝ@��9@��9@���@���@�j@�(�@�  @���@�dZ@�S�@�K�@�K�@�C�@�C�@�33@�+@�33@�+@�"�@�"�@�"�@�"�@�"�@�"�@��@��@��@��@��@���@���@���@��\@�V@�^5@�M�@�E�@�V@�M�@�M�@�-@��@��-@�G�@��@�V@�V@��@�&�@�&�@�&�@�/@���@��-@���@��w@�@��R@��R@��R@�ȴ@��@���@��@���@�ȴ@���@��H@��y@��y@�E�@�=q@�@��h@��@�G�@�7L@�7L@�/@���@��9@��@�A�@���@�ȴ@���@��T@�X@���@��w@�|�@�o@�@��@��R@��\@��+@���@���@��R@��!@��!@��!@�+@�r�@��@���@��@��@�Ĝ@��@��D@�Z@�(�@�1@��m@�dZ@�C�@�\)@�C�@�
=@��!@�v�@�-@�@���@�hs@�V@��`@���@�z�@�A�@�(�@�b@��@��@� �@��m@��@��@�hs@���@��@���@���@�1'@�1@���@���@�l�@�;d@�;d@�C�@�
=@��R@���@���@��H@��@��@��@��@��H@�ff@�5?@�{@���@��T@���@��@���@��@�$�@�{@�J@�@���@��@���@��-@���@��@�O�@�X@��@�V@�V@���@���@�1'@� �@�b@���@���@��F@���@���@���@��@�ƨ@���@��P@�|�@�dZ@�K�@���@�ȴ@���@�n�@��h@�O�@�7L@��/@��j@��@��D@��@�ff@��@��@��^@�z�@��
@�\)@�"�@�v�@��h@���@���@� �@�l�@�"�@��@��!@�=q@���@��h@��@�O�@�?}@��/@��@�ƨ@���@��@�
=@���@��+@�=q@��@��^@�@��7@�X@�/@��`@��D@�1'@��@�@�P@l�@�@�@��@~��@}�@}�@}p�@|I�@{�
@{�F@{"�@z~�@z=q@y�#@xĜ@x  @w�@w��@xb@x �@w��@wK�@wK�@w
=@w
=@v�R@vȴ@v�y@wK�@w��@x�@y�^@zJ@z=q@y�@yhs@yG�@x��@x��@y%@y�@y%@y��@y��@y��@z=q@z~�@z��@z��@z�H@z�H@z�H@z�H@z�H@z��@z��@z��@z��@z��@z=q@y�@y��@y�@y��@yX@yX@yhs@yx�@y��@y��@y%@x��@w�@w�P@w�@w�;@w�@xQ�@w��@xA�@xĜ@x�`@x�`@y&�@yX@yG�@y%@x��@xQ�@w�@w�@w�@w
=@v�@vȴ@v�+@v$�@v@u�@u�T@u�T@u�@u�@u�T@u��@u�@up�@u�@t�j@t�j@t�j@tj@tZ@t(�@s��@st�@sC�@r��@q�7@p��@pQ�@p��@o�;@o\)@o+@ol�@o��@o��@o��@ol�@n�y@nV@nE�@m@m`B@m�@l�/@l�@lj@k�m@k�F@kƨ@k�
@kt�@j�@j��@j=q@iG�@h��@h �@g��@g�@f�@fȴ@f�+@fff@fE�@fE�@f5?@e�@e�T@e��@e@e�h@eO�@e�@d�/@dZ@d�@c��@c��@c��@c"�@b�\@a��@a��@a�@a&�@`��@`�@` �@`  @_�P@_+@^��@^V@]@]�@]/@\�@\�j@\��@\z�@[dZ@Z��@Z�@Y��@Yhs@Y&�@X��@X�9@X�9@X�u@XA�@W�;@W��@W��@W��@Wl�@W�@V�@V��@U�T@U�h@Up�@UV@T��@T�@T�@T�@T�j@T9X@T�@T�@T�@T1@S�m@S��@St�@St�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��H@��H@��@��R@��R@��@��y@��y@���@���@�V@�M�@�=q@�5?@��@��@�$�@�$�@�{@��@�{@�J@�J@�J@�@��@���@��@���@���@��^@�x�@�`B@�O�@��@�O�@�V@�%@��@��;@��w@��F@���@�K�@���@�-@�E�@�V@�M�@���@�@��^@���@��-@��-@��^@��-@���@��-@��^@��-@��-@��@�G�@�&�@���@�Ĝ@�Ĝ@�Ĝ@��9@��9@���@���@�j@�(�@�  @���@�dZ@�S�@�K�@�K�@�C�@�C�@�33@�+@�33@�+@�"�@�"�@�"�@�"�@�"�@�"�@��@��@��@��@��@���@���@���@��\@�V@�^5@�M�@�E�@�V@�M�@�M�@�-@��@��-@�G�@��@�V@�V@��@�&�@�&�@�&�@�/@���@��-@���@��w@�@��R@��R@��R@�ȴ@��@���@��@���@�ȴ@���@��H@��y@��y@�E�@�=q@�@��h@��@�G�@�7L@�7L@�/@���@��9@��@�A�@���@�ȴ@���@��T@�X@���@��w@�|�@�o@�@��@��R@��\@��+@���@���@��R@��!@��!@��!@�+@�r�@��@���@��@��@�Ĝ@��@��D@�Z@�(�@�1@��m@�dZ@�C�@�\)@�C�@�
=@��!@�v�@�-@�@���@�hs@�V@��`@���@�z�@�A�@�(�@�b@��@��@� �@��m@��@��@�hs@���@��@���@���@�1'@�1@���@���@�l�@�;d@�;d@�C�@�
=@��R@���@���@��H@��@��@��@��@��H@�ff@�5?@�{@���@��T@���@��@���@��@�$�@�{@�J@�@���@��@���@��-@���@��@�O�@�X@��@�V@�V@���@���@�1'@� �@�b@���@���@��F@���@���@���@��@�ƨ@���@��P@�|�@�dZ@�K�@���@�ȴ@���@�n�@��h@�O�@�7L@��/@��j@��@��D@��@�ff@��@��@��^@�z�@��
@�\)@�"�@�v�@��h@���@���@� �@�l�@�"�@��@��!@�=q@���@��h@��@�O�@�?}@��/@��@�ƨ@���@��@�
=@���@��+@�=q@��@��^@�@��7@�X@�/@��`@��D@�1'@��@�@�P@l�@�@�@��@~��@}�@}�@}p�@|I�@{�
@{�F@{"�@z~�@z=q@y�#@xĜ@x  @w�@w��@xb@x �@w��@wK�@wK�@w
=@w
=@v�R@vȴ@v�y@wK�@w��@x�@y�^@zJ@z=q@y�@yhs@yG�@x��@x��@y%@y�@y%@y��@y��@y��@z=q@z~�@z��@z��@z�H@z�H@z�H@z�H@z�H@z��@z��@z��@z��@z��@z=q@y�@y��@y�@y��@yX@yX@yhs@yx�@y��@y��@y%@x��@w�@w�P@w�@w�;@w�@xQ�@w��@xA�@xĜ@x�`@x�`@y&�@yX@yG�@y%@x��@xQ�@w�@w�@w�@w
=@v�@vȴ@v�+@v$�@v@u�@u�T@u�T@u�@u�@u�T@u��@u�@up�@u�@t�j@t�j@t�j@tj@tZ@t(�@s��@st�@sC�@r��@q�7@p��@pQ�@p��@o�;@o\)@o+@ol�@o��@o��@o��@ol�@n�y@nV@nE�@m@m`B@m�@l�/@l�@lj@k�m@k�F@kƨ@k�
@kt�@j�@j��@j=q@iG�@h��@h �@g��@g�@f�@fȴ@f�+@fff@fE�@fE�@f5?@e�@e�T@e��@e@e�h@eO�@e�@d�/@dZ@d�@c��@c��@c��@c"�@b�\@a��@a��@a�@a&�@`��@`�@` �@`  @_�P@_+@^��@^V@]@]�@]/@\�@\�j@\��@\z�@[dZ@Z��@Z�@Y��@Yhs@Y&�@X��@X�9@X�9@X�u@XA�@W�;@W��@W��@W��@Wl�@W�@V�@V��@U�T@U�h@Up�@UV@T��@T�@T�@T�@T�j@T9X@T�@T�@T�@T1@S�m@S��@St�@St�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB7LB6FB6FB7LB6FB7LB6FB7LB8RB9XB8RB7LB8RB8RB8RB9XB8RB8RB8RB8RB8RB8RB8RB8RB8RB8RB9XB:^B9XB:^B9XB:^B;dB;dB;dB<jB:^B9XB:^B9XB:^B8RB8RB7LB7LB6FB6FB5?B5?B5?B49B5?B49B49B49B49B49B49B33B49B33B33B33B33B33B33B2-B2-B1'B1'B1'B1'B1'B0!B0!B0!B/B/B/B.B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B-B,B,B,B,B,B+B+B+B+B+B+B+B+B)�B(�B'�B'�B'�B'�B'�B'�B'�B'�B'�B(�B)�B%�B"�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BhBVBPB
=B+BBBB  B  B��B��B��B��B��B  B  BBBBBDBVBVBVBVBVBVBVBPBPBPBJBDBDBJBDBDB
=B	7B1B1B+B%BBBBBBBBBBB  B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�sB�mB�mB�mB�mB�fB�fB�fB�mB�yB�yB�sB�mB�sB�sB�mB�mB�fB�`B�ZB�TB�BB�;B�;B�5B�/B�/B�)B�B��B��B��B��B��BɺBǮBƨBĜB��B�}B�qB�jB�^B�XB�XB�RB�FB�?B�9B�9B�9B�9B�9B�-B�!B�B�B�B�B�B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�-B�-B�3B�?B�FB�LB�XB�^B�^B�dB�jB�jB�jB�qB�qB�wB�wB�wB�wB�}B��B��B��B��B��BBÖBĜBÖBB��B��BBÖBÖBŢBŢBǮBȴBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBǮBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBǮBǮBƨBƨBƨBƨBǮBƨBƨBŢBŢBŢBƨBƨBƨBƨBƨBƨBŢBŢBŢBŢBŢBŢBŢBŢBŢBƨBĜBÖBBBBBBBBBBBBBBBBBBBBBBBÖBÖBÖBÖBÖBĜBĜBÖBÖBĜBĜBĜBĜ11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B0�B/�B/�B0�B/�B0�B/�B0�B1�B2�B1�B0�B1�B1�B1�B2�B1�B1�B1�B1�B1�B1�B1�B1�B1�B1�B2�B3�B2�B3�B2�B3�B4�B4�B4�B5�B3�B2�B3�B2�B3�B1�B1�B0�B0�B/�B/�B.�B.�B.�B-zB.�B-zB-zB-zB-zB-zB-zB,tB-zB,tB,tB,tB,tB,tB,tB+nB+nB*hB*hB*hB*hB*hB)bB)bB)bB(]B(]B(]B'VB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB&PB%JB%JB%JB%JB%JB$DB$DB$DB$DB$DB$DB$DB$DB#>B"8B!2B!2B!2B!2B!2B!2B!2B!2B!2B"8B#>B%BB�B�B�B�B�B�B�B�B�B�BBBB�B�B�B�B�B�B�B�B�B�B�B�B�B�B
�B�B�B�B oB�cB�PB�JB�EB�EB�?B�?B�?B�?B�?B�EB�EB�JB�JB�PB�]B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBuB oB�iB�]B�]B�VB�WB�QB�KB�KB�KB�KB�KB�EB�,B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B߬B߬B߬B�B�B�B�B�B�B�B�B�B߬BާBݡBܛBىB؂B؂B�|B�vB�vB�pB�^B�@B�4B�-B�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�wB�kB�eB�eB�^B�XB�RB�RB�MB�GB�MB�GB�GB�AB�;B�4B�.B�.B�.B�.B�.B�4B�4B�4B�(B�"B�"B�"B�B�B�B�	B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�
B�B�"B�.B�GB�SB�YB�YB�YB�YB�YB�YB�_B�fB�fB�wB�wB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�"B�"B�(B�(B�(B�(B�(B�"B�"B�"B�"B�B�"B�"B�"B�"B�B�B�B�B�B�B�B�
B�
B�B�"B�(B�(B�(B�"B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0065667                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904482018110709044820181107090448  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172716  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172716  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090448  IP  PSAL            @��D�� G�O�                