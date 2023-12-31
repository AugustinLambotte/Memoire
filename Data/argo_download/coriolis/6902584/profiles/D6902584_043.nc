CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:18Z creation; 2018-11-05T17:27:22Z last update (coriolis COQC software)   
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
resolution        =���   axis      Z          9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ^   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       `   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   h,   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       j4   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       rL   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   zd   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       |l   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �    HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �P   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �`   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �d   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �t   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �x   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �|   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172618  20181107090454  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               +A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�eOz��1   @�eO��,@O#I'�>��C��/v�1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�33@�33A   A  A   A0  A@  AP  A`  Ap  A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\ffB`  Bd  Bh  Bl  BpffBt  Bx  B|ffB�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�33B�33B�  B�  B�  B���B���B�  B�33B�33B�  B�  B�  B���B�  B�  B���B�  B�  B�  B�  B�  B�  C�C� C  C	� C  CffC  C��C  CffC  C� C   C"� C%  C'� C*  C,� C.�fC1� C4  C6� C9  C;� C>  C@� CC  CE� CH�CJ��CM  CO� CR  CT� CW  CY� C\  C^ffCa  Cc��Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C��3C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ Cŀ C�� C�  C�33Cʀ C���C�  C�@ Cπ C�� C�  C�@ CԀ C�� C��3C�33Cـ C���C�  C�33C�s3C�� C��C�@ C� C�� C��3C�@ C��C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C��C�� C�  D �fD  D@ D� D� D  D@ D	y�D
� D  D9�D� D� D  D@ D� D� D��D@ D� D� D��D@ D� D��D   D!@ D"y�D#��D%  D&@ D'� D(� D*  D+9�D,� D-� D/  D0@ D1� D2� D4  D5@ D6�fD7� D9fD:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY�fDZ�fD\fD]@ D^� D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� DkfDlFfDm� Dn��Do��Dq9�Dr� Ds� Du  Dv9�Dwy�Dx��Dz  D{@ D|� D}� D  D��D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  DƼ�D�` D�  DȠ D�C3D�� Dʀ D�  D�� D�c3D�3Dͣ3D�@ D�� Dπ D�  Dм�D�\�D���DҜ�D�@ D�� DԀ D�  D�� D�\�D�  Dף3D�@ D�� Dـ D�  D��3D�` D�  Dܜ�D�<�D���Dހ D�  D�� D�Y�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@@  @�  @�  @�33@�33A   A  A   A0  A@  AP  A`  Ap  A�  A���A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BH  BL  BP  BT  BX  B\ffB`  Bd  Bh  Bl  BpffBt  Bx  B|ffB�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�33B�33B�  B�  B�  B���B���B�  B�33B�33B�  B�  B�  B���B�  B�  B���B�  B�  B�  B�  B�  B�  C�C� C  C	� C  CffC  C��C  CffC  C� C   C"� C%  C'� C*  C,� C.�fC1� C4  C6� C9  C;� C>  C@� CC  CE� CH�CJ��CM  CO� CR  CT� CW  CY� C\  C^ffCa  Cc��Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�L�C�� C��3C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ Cŀ C�� C�  C�33Cʀ C���C�  C�@ Cπ C�� C�  C�@ CԀ C�� C��3C�33Cـ C���C�  C�33C�s3C�� C��C�@ C� C�� C��3C�@ C��C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C�� C�� C��C�� C�  D �fD  D@ D� D� D  D@ D	y�D
� D  D9�D� D� D  D@ D� D� D��D@ D� D� D��D@ D� D��D   D!@ D"y�D#��D%  D&@ D'� D(� D*  D+9�D,� D-� D/  D0@ D1� D2� D4  D5@ D6�fD7� D9fD:@ D;� D<� D>  D?@ D@� DA� DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DW  DXFfDY�fDZ�fD\fD]@ D^� D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� DkfDlFfDm� Dn��Do��Dq9�Dr� Ds� Du  Dv9�Dwy�Dx��Dz  D{@ D|� D}� D  D��D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�C3D��3D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D���D�|�D�  D�� D�` D�  D���D�<�D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D���D�@ D�� D�|�D�  D�� D�c3D�3D�� D�@ D�� D�� D�  D�� D�` D���D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  Dà D�@ D�� Dŀ D�  DƼ�D�` D�  DȠ D�C3D�� Dʀ D�  D�� D�c3D�3Dͣ3D�@ D�� Dπ D�  Dм�D�\�D���DҜ�D�@ D�� DԀ D�  D�� D�\�D�  Dף3D�@ D�� Dـ D�  D��3D�` D�  Dܜ�D�<�D���Dހ D�  D�� D�Y�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�A��A��A	l�AG�A�HA��A�A$�A��AXAr�A;dA ffA  �@�C�@���@�E�@�+@���@���@�M�@�$�@���@�z�@�l�@�ff@�V@�K�@��H@��H@��@�$�@�hs@�@�l�@��H@ꗍ@�\@�v�@�5?@��#@�G�@���@�Z@�S�@�ff@�E�@���@�^@�O�@�%@�Z@�l�@���@�?}@��@��@�Ĝ@�@��@�ȴ@�n�@��@��#@�x�@�V@��@���@ܓu@�9X@�1@ۮ@�dZ@�o@ڇ+@�-@��@���@�@��`@�I�@��@���@׶F@׍P@�l�@�S�@�C�@�+@�
=@��@��y@ְ!@�v�@�-@���@ղ-@Ձ@��@�bN@���@ӝ�@ҸR@�X@мj@�Q�@�t�@�5?@�hs@���@���@�~�@�J@�O�@��
@��T@ř�@�x�@�hs@�?}@��@�Ĝ@Ĭ@ċD@��
@��H@�5?@��7@�7L@�%@�Ĝ@�z�@��@�C�@���@�{@���@��j@�(�@�ƨ@���@�C�@��#@�`B@�X@�hs@���@���@���@�hs@�j@�A�@�;d@���@���@��@�dZ@��!@�V@�=q@�p�@��j@��;@��y@�ff@�%@�I�@�t�@�33@��R@�A�@�l�@�S�@�33@���@��!@�E�@�$�@��@��@��T@��T@��@��@��@�/@��/@���@���@�9X@��m@��P@��H@��^@��@���@�7L@��@��@�V@�%@���@�r�@�1'@�1@���@���@�ƨ@��
@��@��@�\)@�@��@�V@�x�@�/@��9@�(�@�l�@�
=@�-@��7@��h@�hs@�7L@��@��@���@���@���@�p�@��@��@�@��-@��-@��@���@���@�p�@�?}@�X@�O�@�G�@��/@���@�t�@�;d@�;d@�33@�33@�K�@�l�@�t�@�dZ@�"�@��@��!@��!@���@��!@���@�ȴ@�
=@�33@�t�@���@��@���@��T@���@���@�7L@�&�@��@��@��@��9@�1'@���@�l�@�K�@�K�@�"�@��y@�v�@��#@��@�hs@��`@��@�9X@�b@��@�dZ@�C�@�+@�
=@��@���@�=q@�`B@��@��9@�(�@�b@��@��;@�ƨ@��w@�dZ@��R@�v�@�@��@���@�z�@�j@�bN@�Z@�b@�b@�1@��@��
@�t�@�o@�@�@��y@��y@��@���@���@��\@�v�@�ff@�5?@��#@�`B@�r�@��F@��@�ȴ@�~�@�n�@�J@���@��-@��7@�%@��@��/@��j@���@��u@�bN@�(�@�@�w@�@�@~��@~ȴ@~��@~�+@~ff@}@}O�@}/@|I�@{o@z�\@z�@y�#@y�^@yG�@x��@x�`@x�@xA�@xb@w�@w
=@u�h@uO�@uV@t��@tz�@tZ@s��@s�m@s��@s��@sS�@sC�@s@r�\@r^5@r-@qhs@p��@p�@p �@o��@o��@n��@n5?@m��@l9X@kC�@j��@j�@i��@iX@h��@hb@g;d@fV@f@e@e`B@eV@c��@co@b=q@a�^@`�9@`��@` �@_�w@^��@]@\�/@[�m@Z��@Z^5@Y�#@Y�@X�u@W�w@W+@V��@V{@U�@UO�@U/@T�j@T(�@S�F@Q��@P��@PA�@O\)@N�R@N@Lj@K��@KC�@I��@G�w@G�@FV@E�-@E�@E?}@EV@D��@C�
@C33@A��@?�@>{@=O�@<z�@;S�@9��@8  @7l�@5�T@3��@2�@0Ĝ@/�w@.@-�h@-/@,��@+33@)��@(bN@'\)@&ȴ@%?}@$��@$(�@$�@#dZ@"�@"��@"M�@!&�@ ��@ ��@�w@V@p�@��@��@  @K�@o@+@
��@
��@
��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�A��A��A	l�AG�A�HA��A�A$�A��AXAr�A;dA ffA  �@�C�@���@�E�@�+@���@���@�M�@�$�@���@�z�@�l�@�ff@�V@�K�@��H@��H@��@�$�@�hs@�@�l�@��H@ꗍ@�\@�v�@�5?@��#@�G�@���@�Z@�S�@�ff@�E�@���@�^@�O�@�%@�Z@�l�@���@�?}@��@��@�Ĝ@�@��@�ȴ@�n�@��@��#@�x�@�V@��@���@ܓu@�9X@�1@ۮ@�dZ@�o@ڇ+@�-@��@���@�@��`@�I�@��@���@׶F@׍P@�l�@�S�@�C�@�+@�
=@��@��y@ְ!@�v�@�-@���@ղ-@Ձ@��@�bN@���@ӝ�@ҸR@�X@мj@�Q�@�t�@�5?@�hs@���@���@�~�@�J@�O�@��
@��T@ř�@�x�@�hs@�?}@��@�Ĝ@Ĭ@ċD@��
@��H@�5?@��7@�7L@�%@�Ĝ@�z�@��@�C�@���@�{@���@��j@�(�@�ƨ@���@�C�@��#@�`B@�X@�hs@���@���@���@�hs@�j@�A�@�;d@���@���@��@�dZ@��!@�V@�=q@�p�@��j@��;@��y@�ff@�%@�I�@�t�@�33@��R@�A�@�l�@�S�@�33@���@��!@�E�@�$�@��@��@��T@��T@��@��@��@�/@��/@���@���@�9X@��m@��P@��H@��^@��@���@�7L@��@��@�V@�%@���@�r�@�1'@�1@���@���@�ƨ@��
@��@��@�\)@�@��@�V@�x�@�/@��9@�(�@�l�@�
=@�-@��7@��h@�hs@�7L@��@��@���@���@���@�p�@��@��@�@��-@��-@��@���@���@�p�@�?}@�X@�O�@�G�@��/@���@�t�@�;d@�;d@�33@�33@�K�@�l�@�t�@�dZ@�"�@��@��!@��!@���@��!@���@�ȴ@�
=@�33@�t�@���@��@���@��T@���@���@�7L@�&�@��@��@��@��9@�1'@���@�l�@�K�@�K�@�"�@��y@�v�@��#@��@�hs@��`@��@�9X@�b@��@�dZ@�C�@�+@�
=@��@���@�=q@�`B@��@��9@�(�@�b@��@��;@�ƨ@��w@�dZ@��R@�v�@�@��@���@�z�@�j@�bN@�Z@�b@�b@�1@��@��
@�t�@�o@�@�@��y@��y@��@���@���@��\@�v�@�ff@�5?@��#@�`B@�r�@��F@��@�ȴ@�~�@�n�@�J@���@��-@��7@�%@��@��/@��j@���@��u@�bN@�(�@�@�w@�@�@~��@~ȴ@~��@~�+@~ff@}@}O�@}/@|I�@{o@z�\@z�@y�#@y�^@yG�@x��@x�`@x�@xA�@xb@w�@w
=@u�h@uO�@uV@t��@tz�@tZ@s��@s�m@s��@s��@sS�@sC�@s@r�\@r^5@r-@qhs@p��@p�@p �@o��@o��@n��@n5?@m��@l9X@kC�@j��@j�@i��@iX@h��@hb@g;d@fV@f@e@e`B@eV@c��@co@b=q@a�^@`�9@`��@` �@_�w@^��@]@\�/@[�m@Z��@Z^5@Y�#@Y�@X�u@W�w@W+@V��@V{@U�@UO�@U/@T�j@T(�@S�F@Q��@P��@PA�@O\)@N�R@N@Lj@K��@KC�@I��@G�w@G�@FV@E�-@E�@E?}@EV@D��@C�
@C33@A��@?�@>{@=O�@<z�@;S�@9��@8  @7l�@5�T@3��@2�@0Ĝ@/�w@.@-�h@-/@,��@+33@)��@(bN@'\)@&ȴ@%?}@$��@$(�@$�@#dZ@"�@"��@"M�@!&�@ ��@ ��@�w@V@p�@��@��@  @K�@o@+@
��@
��@
��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB/B7LB:^BJ�BI�BK�BK�BJ�BJ�BS�BS�BVBYBXBW
BW
BVB]/B]/B\)B\)B\)B[#BZB]/B\)B\)B`BBbNBbNBaHBaHBcTBcTBffBffBffBffBe`Be`Be`BdZBdZBdZBdZBhsBjBk�Bl�Bm�Bl�Bk�Bk�Bl�Bo�Bo�Bo�Bn�Bn�Bm�Bn�Bp�Bp�Bp�Bo�Bp�Bq�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bq�Bp�Bp�Bp�Bo�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bq�Bp�Bo�Bo�Bn�Bm�Bn�Bm�Bm�Bl�Bl�Bk�BjBiyBiyBhsBffBdZBdZBdZBcTBcTBcTBcTBcTBbNBbNB`BB`BB_;B^5B^5B]/B]/B\)B[#BZBXBW
BVBS�BR�BQ�BP�BO�BM�BL�BL�BL�BM�BL�BL�BK�BJ�BH�BF�BB�B@�B>wB=qB<jB;dB:^B8RB5?B33B0!B.B+B(�B&�B%�B#�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BoBbB\BJB+B	7B	7B1B1B1B	7B	7B1B1B+B%B%BBBBBBBBB��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�fB�fB�fB�fB�fB�mB�mB�sB�mB�fB�fB�fB�fB�fB�mB�yB�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�sB�mB�`B�`B�`B�ZB�ZB�TB�NB�HB�HB�BB�;B�;B�5B�5B�/B�/B�/B�/B�)B�)B�#B�B�B�B�
B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBɺBɺBɺBɺBɺBɺBɺBɺBȴBǮBŢBĜBÖBÖBÖBBBBBBBBBBB��B��B��B��B��B��B�}B�}B�}B�}B��B��B��B��B��B��B�}B��B��B��B��B��B��B��B��BB��BBBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBÖBÖBÖBÖBBB��B��B��B��B��B�}B�}B�}B�wB�wB�wB�qB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�XB�XB�RB�RB�RB�RB�LB�LB�LB�LB�LB�FB�?B�9B�9B�9B�3B�-B�-B�-B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B%B-NB0`B@�B?�BA�BA�B@�B@�BI�BI�BLBOBNBMBMBLBS1BS1BR+BR+BR+BQ%BPBS1BR+BR+BVDBXPBXPBWJBWJBYVBYVB\hB\hB\hB\hB[bB[bB[bBZ\BZ\BZ\BZ\B^uB`�Ba�Bb�Bc�Bb�Ba�Ba�Bb�Be�Be�Be�Bd�Bd�Bc�Bd�Bf�Bf�Bf�Be�Bf�Bg�Bf�Bf�Bf�Bf�Bf�Bf�Bf�Bf�Bg�Bf�Bf�Bf�Be�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bg�Bf�Be�Be�Bd�Bc�Bd�Bc�Bc�Bb�Bb�Ba�B`�B_}B_}B^wB\jBZ_BZ_BZ_BYYBYYBYYBYYBYYBXSBXSBVGBVGBU@BT:BT:BS5BS5BR/BQ)BP#BNBMBL
BI�BH�BG�BF�BE�BC�BB�BB�BB�BC�BB�BB�BA�B@�B>�B<�B8�B6�B4�B3zB2sB1mB0gB.[B+IB)=B&+B$B!BB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B
�B|BoBiBXB�9B�EB�EB�?B�?B�?B�EB�EB�?B�?B�9B�3B�3B�-B�-B�-B�-B�'B�!B�B�B�	B��B��B��B��B��B�B�B�B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�wB�wB�wB�wB�wB�~B�~BބB�~B�wB�wB�wB�wB�wB�~BߊB��B�B�B�B�B�B�B�B�B��BߊBߊBߊBߊBߊBބB�~B�rB�rB�rB�lB�lB�fB�`B�ZB�ZB�TB�MB�MB�GB�HB�BB�BB�BB�BB�<B�<B�6B�0B�*B�#B�B�B�B�B�B�B�B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�|B�|B�|B�|B�|B�|B�|B�|B�vB�vB�vB�pB�pB�jB�jB�jB�jB�eB�eB�eB�eB�eB�_B�XB�RB�RB�RB�LB�FB�FB�FB�@B�;B�;B�5B�5B�5B�5B�.B�.B�.B�(B�(B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0097078                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904542018110709045420181107090454  IF  ARFMCODA024c                                                                20181105172618                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172722  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172722  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090454  IP  PSAL            @ffD�Y�G�O�                