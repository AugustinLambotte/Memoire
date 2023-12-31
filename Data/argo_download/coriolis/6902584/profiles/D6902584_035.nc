CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  3   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:17Z last update (coriolis COQC software)   
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
_FillValue                  ,  �0             ,  �0Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090449  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               #A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�QQ�~1   @�QQ�}�U@O��ve?�@q�4gz+1   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@@  @�  @�  @���@�  A��A  A   A0  A@  AP  A`  AnffA~ffA�33A�33A�  A���A�  A���A�  A�  A�  A�33A�33A�33A�33A�33A�  B   BffB  B��B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@ffBD  BH  BL  BPffBT  BX  B\  B`  Bd  Bh  BlffBp  Bt  Bx  B|  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C�fC	� C  C� C  C� C  C��C  C� C   C"� C%  C'� C*�C,� C/  C1� C4�C6��C9  C;� C>  C@� CC�CE��CH�CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Ce�fCh� Ck�Cm� Cp  Cr� Cu  CwffCz  C|� C  C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C��3C�@ C���C���C��C�L�C���C���C��C�@ C�s3C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԌ�C���C�  C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�3C��3C�@ C��C�� C�  C�@ C�� C�� C��C�� C�  D � D  D@ D� D� D  D@ D	� D
��D��D@ D� D� DfDFfD�fD�fDfDFfD� D� D  D@ D� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,y�D-� D/  D0FfD1�fD2� D3��D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@y�DA��DC  DD@ DE� DF� DH  DI@ DJy�DK��DM  DN@ DO� DP� DR  DS@ DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^y�D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDpfDqFfDr� Ds� Du  Dv@ Dw� Dx� DzfD{FfD|�fD}� D  D�  D�� D�` D���D���D�<�D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�c3D�  D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D��3D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D��3D�� D�  D���D�` D�3Dã3D�C3D�� Dŀ D�  D�� D�` D�  DȠ D�<�D�� Dʀ D��D�� D�` D�  D͠ D�<�D���Dπ D�#3D�� D�` D�3DҠ D�@ D�� Dԃ3D�#3D�� D�\�D�  Dנ D�@ D�� Dـ D�  D�� D�\�D�3Dܠ D�@ D�� Dހ D�  D�� D�c3D�  D� D�@ D�� D�3D�#3D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�\�D���D� D�@ D�� D� D�  D�� D�` D�  D���D�@ D�� D�� D��D�� D�c3D�  D�� D�@ D�� D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@@  @�  @�  @���@�  A��A  A   A0  A@  AP  A`  AnffA~ffA�33A�33A�  A���A�  A���A�  A�  A�  A�33A�33A�33A�33A�33A�  B   BffB  B��B  B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@ffBD  BH  BL  BPffBT  BX  B\  B`  Bd  Bh  BlffBp  Bt  Bx  B|  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  C  C� C�fC	� C  C� C  C� C  C��C  C� C   C"� C%  C'� C*�C,� C/  C1� C4�C6��C9  C;� C>  C@� CC�CE��CH�CJ� CM  CO� CR  CT� CW  CY� C\  C^� Ca  Cc� Ce�fCh� Ck�Cm� Cp  Cr� Cu  CwffCz  C|� C  C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C��3C�@ C���C���C��C�L�C���C���C��C�@ C�s3C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C�  C�@ Cʀ C�� C�  C�@ Cπ C�� C�  C�@ CԌ�C���C�  C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C��C�@ C� C�3C��3C�@ C��C�� C�  C�@ C�� C�� C��C�� C�  D � D  D@ D� D� D  D@ D	� D
��D��D@ D� D� DfDFfD�fD�fDfDFfD� D� D  D@ D� D� D   D!FfD"� D#� D%  D&@ D'� D(� D*  D+@ D,y�D-� D/  D0FfD1�fD2� D3��D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@y�DA��DC  DD@ DE� DF� DH  DI@ DJy�DK��DM  DN@ DO� DP� DR  DS@ DT� DU�fDW  DX@ DY� DZ� D\  D]@ D^y�D_� Da  Db@ Dc� Dd�fDf  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDpfDqFfDr� Ds� Du  Dv@ Dw� Dx� DzfD{FfD|�fD}� D  D�  D�� D�` D���D���D�<�D���D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�c3D�  D���D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�\�D�  D�� D�@ D�� D��3D�#3D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�|�D�  D��3D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D��3D�� D�  D���D�` D�3Dã3D�C3D�� Dŀ D�  D�� D�` D�  DȠ D�<�D�� Dʀ D��D�� D�` D�  D͠ D�<�D���Dπ D�#3D�� D�` D�3DҠ D�@ D�� Dԃ3D�#3D�� D�\�D�  Dנ D�@ D�� Dـ D�  D�� D�\�D�3Dܠ D�@ D�� Dހ D�  D�� D�c3D�  D� D�@ D�� D�3D�#3D�� D�` D�  D� D�@ D�� D� D�  D�� D�c3D�  D� D�@ D�� D� D�  D�� D�\�D���D� D�@ D�� D� D�  D�� D�` D�  D���D�@ D�� D�� D��D�� D�c3D�  D�� D�@ D�� D��311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�O�@�X@�hs@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�p�@�hs@�hs@�hs@�hs@�hs@�X@�O�@�G�@�?}@�7L@�7L@�7L@�?}@�?}@�?}@�?}@�X@�hs@�`B@�X@�O�@�G�@�O�@�hs@�`B@�p�@�x�@�hs@�`B@�&�@�X@���@��`@���@��@��@�%@���@���@�Ĝ@���@��@�K�@��7@���@��j@� �@�1@�1@��m@��@�\)@�33@�+@�o@�@��y@���@��R@���@�E�@��T@��h@�?}@���@��@�j@�1@�  @�  @��m@��P@�+@���@��R@�=q@��@�`B@�%@��`@��/@���@��9@���@�I�@��m@�\)@�@��H@��@��+@�M�@���@��@�/@��@��`@��u@�A�@�9X@�1'@��m@���@��@�\)@�@���@�~�@�^5@�=q@��@��@�@���@���@�O�@�/@��@�%@��@���@��/@���@��D@�r�@�Ĝ@���@��u@��D@�z�@�9X@�  @�ƨ@���@���@�\)@��y@��!@�^5@��T@���@�/@���@���@�Ĝ@�j@� �@��;@�C�@�+@�C�@�l�@��@�S�@�S�@�o@��@�n�@�-@�-@�E�@�{@���@��^@�@��T@���@���@���@���@�O�@��@��/@��@�  @��m@���@��@���@�;d@�v�@�@���@�X@��#@��7@�%@��D@�z�@�Z@��@� �@� �@��@�b@�1@���@�  @��w@��P@�|�@�t�@�K�@�dZ@���@���@��@�|�@�l�@�C�@�@���@�^5@�E�@�5?@�{@��@��#@���@��@���@�J@�{@�{@�@�@�@���@��@���@��^@�@���@��h@��7@�X@�G�@�&�@�V@�%@��/@��/@��/@���@��9@���@��u@�r�@�Z@��@���@���@��@��@�K�@�o@���@�o@�ff@�@���@��^@���@��^@��@���@�=q@��T@�?}@�z�@�I�@�  @���@��@�ƨ@��P@�dZ@�C�@�
=@��H@��\@�V@�=q@�J@��#@��-@���@��h@�hs@�&�@��/@�Ĝ@���@�j@�1'@�@�P@�@~ȴ@~5?@~@}@}�@|��@}�@|�@|9X@|(�@|�@{�F@{��@|j@|�D@{�
@{@zJ@y��@yX@y&�@y�@x��@xQ�@xb@w�w@w��@w|�@w�@v��@w
=@v�@v�+@u�-@u��@u�T@v{@v@v@u�@v@u�h@u`B@u�-@u�-@u�@t�@u@v��@w�@yG�@z^5@{�@{��@{t�@{C�@{�F@{ƨ@{��@{��@zJ@yhs@yx�@y��@y�#@y�^@zJ@y��@y�#@y&�@x �@xbN@yhs@yx�@v��@v��@vv�@vE�@v�R@wK�@w\)@w��@xr�@x�u@x  @v�+@vE�@w+@w�@x1'@xQ�@x  @w|�@w�@w�@w��@w��@w��@w�;@w��@x��@x��@x��@x��@x�`@y%@y7L@y&�@y�@x�`@x��@xbN@xA�@xb@w|�@vff@v5?@t��@t(�@t�@s�F@s�@st�@sC�@so@r��@r��@r��@r��@s@sC�@s�@s33@r�!@q��@q&�@q%@p�u@o�w@o�P@o�w@o��@o�@o|�@oK�@o
=@n��@n�@n��@nff@nV@nV@nff@n{@m��@mp�@m?}@l�/@lj@k��@k�@kS�@k"�@j��@j=q@i�#@ix�@jJ@j=q@iG�@hbN@h  @g�;@g�;@h1'@h �@g�@f��@f��@f$�@ep�@eO�@d�@d�j@d9X@c�m@cC�@b��@bn�@b=q@a�#@ahs@aX@a�^@a�^@ax�@ax�@`��@` �@_�w@_l�@_K�@_K�@_K�@_�@^�@^�+@^5?@^$�@]@]p�@]p�@]O�@\�/@\��@\��@\�@\��@\1@[�@[�
@[�m@[��@[�F@[��@[�F@[C�@[@Z��@Z�!@Zn�@Z=q@Y��@Y�#@Y�7@Yx�@YG�@Y7L@X��@X��@X��@XA�@W��@W��@V��@U@U�@T�@S�m@S��@T1@S��@S�
@S�F11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�O�@�X@�hs@�p�@�p�@�p�@�p�@�p�@�hs@�hs@�p�@�hs@�hs@�hs@�hs@�hs@�X@�O�@�G�@�?}@�7L@�7L@�7L@�?}@�?}@�?}@�?}@�X@�hs@�`B@�X@�O�@�G�@�O�@�hs@�`B@�p�@�x�@�hs@�`B@�&�@�X@���@��`@���@��@��@�%@���@���@�Ĝ@���@��@�K�@��7@���@��j@� �@�1@�1@��m@��@�\)@�33@�+@�o@�@��y@���@��R@���@�E�@��T@��h@�?}@���@��@�j@�1@�  @�  @��m@��P@�+@���@��R@�=q@��@�`B@�%@��`@��/@���@��9@���@�I�@��m@�\)@�@��H@��@��+@�M�@���@��@�/@��@��`@��u@�A�@�9X@�1'@��m@���@��@�\)@�@���@�~�@�^5@�=q@��@��@�@���@���@�O�@�/@��@�%@��@���@��/@���@��D@�r�@�Ĝ@���@��u@��D@�z�@�9X@�  @�ƨ@���@���@�\)@��y@��!@�^5@��T@���@�/@���@���@�Ĝ@�j@� �@��;@�C�@�+@�C�@�l�@��@�S�@�S�@�o@��@�n�@�-@�-@�E�@�{@���@��^@�@��T@���@���@���@���@�O�@��@��/@��@�  @��m@���@��@���@�;d@�v�@�@���@�X@��#@��7@�%@��D@�z�@�Z@��@� �@� �@��@�b@�1@���@�  @��w@��P@�|�@�t�@�K�@�dZ@���@���@��@�|�@�l�@�C�@�@���@�^5@�E�@�5?@�{@��@��#@���@��@���@�J@�{@�{@�@�@�@���@��@���@��^@�@���@��h@��7@�X@�G�@�&�@�V@�%@��/@��/@��/@���@��9@���@��u@�r�@�Z@��@���@���@��@��@�K�@�o@���@�o@�ff@�@���@��^@���@��^@��@���@�=q@��T@�?}@�z�@�I�@�  @���@��@�ƨ@��P@�dZ@�C�@�
=@��H@��\@�V@�=q@�J@��#@��-@���@��h@�hs@�&�@��/@�Ĝ@���@�j@�1'@�@�P@�@~ȴ@~5?@~@}@}�@|��@}�@|�@|9X@|(�@|�@{�F@{��@|j@|�D@{�
@{@zJ@y��@yX@y&�@y�@x��@xQ�@xb@w�w@w��@w|�@w�@v��@w
=@v�@v�+@u�-@u��@u�T@v{@v@v@u�@v@u�h@u`B@u�-@u�-@u�@t�@u@v��@w�@yG�@z^5@{�@{��@{t�@{C�@{�F@{ƨ@{��@{��@zJ@yhs@yx�@y��@y�#@y�^@zJ@y��@y�#@y&�@x �@xbN@yhs@yx�@v��@v��@vv�@vE�@v�R@wK�@w\)@w��@xr�@x�u@x  @v�+@vE�@w+@w�@x1'@xQ�@x  @w|�@w�@w�@w��@w��@w��@w�;@w��@x��@x��@x��@x��@x�`@y%@y7L@y&�@y�@x�`@x��@xbN@xA�@xb@w|�@vff@v5?@t��@t(�@t�@s�F@s�@st�@sC�@so@r��@r��@r��@r��@s@sC�@s�@s33@r�!@q��@q&�@q%@p�u@o�w@o�P@o�w@o��@o�@o|�@oK�@o
=@n��@n�@n��@nff@nV@nV@nff@n{@m��@mp�@m?}@l�/@lj@k��@k�@kS�@k"�@j��@j=q@i�#@ix�@jJ@j=q@iG�@hbN@h  @g�;@g�;@h1'@h �@g�@f��@f��@f$�@ep�@eO�@d�@d�j@d9X@c�m@cC�@b��@bn�@b=q@a�#@ahs@aX@a�^@a�^@ax�@ax�@`��@` �@_�w@_l�@_K�@_K�@_K�@_�@^�@^�+@^5?@^$�@]@]p�@]p�@]O�@\�/@\��@\��@\�@\��@\1@[�@[�
@[�m@[��@[�F@[��@[�F@[C�@[@Z��@Z�!@Zn�@Z=q@Y��@Y�#@Y�7@Yx�@YG�@Y7L@X��@X��@X��@XA�@W��@W��@V��@U@U�@T�@S�m@S��@T1@S��@S�
@S�F11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B"�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BuBuBoBoBhBhBbB\B\BPBJBDBDB
=B
=B
=B	7B1B+B%BBBBBBBB  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�fB�fB�fB�mB�mB�mB�mB�mB�fB�`B�ZB�ZB�ZB�TB�NB�NB�NB�TB�TB�TB�TB�NB�HB�HB�BB�;B�/B�/B�/B�5B�/B�#B�B�B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBȴBȴBɺB��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBȴBȴBȴBȴBȴBȴBǮBǮBǮBǮBƨBƨBƨBŢBŢBŢBĜBÖBÖBB��BB��B�wB�wB�qB�qB�wB�}BÖB��B��B�qB�XB�XB�RB�LB�LB�LB�FB�?B�?B�9B�3B�-B�-B�-B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�!B�9B�?B�?B�?B�FB�LB�RB�RB�?B�9B�?B�FB�FB�LB�XB�XB�XB�RB�LB�RB�dB�^B�FB�FB�FB�FB�RB�^B�dB�wB��B��B�}B�qB�qB��BÖBĜBŢBŢBĜBƨBƨBƨBǮBǮBǮBȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBɺBɺBɺB��B��B��B��B��B��B��B��BɺBɺBȴBȴBȴBɺBɺBɺBɺBɺBɺBɺBɺBɺBȴBȴBɺB��BɺBɺBɺBɺBȴBȴBȴBǮBǮBǮBǮBǮBǮBƨBȴBȴBƨBŢBĜBĜBŢBƨBƨBƨBŢBŢBĜBÖBÖBBBBB��B��B��B��B��B��B��BÖBÖBÖBĜBÖBBBBBBÖBBBBBÖBBBÖBÖBBBÖBÖBÖBÖBBÖBĜBĜBĜBĜBŢBŢBŢBŢBŢBŢBŢBƨBƨBŢBŢBŢBƨBƨBƨBƨBƨBƨBƨBŢBÖB��B��B�}B��B��BB��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB�BBBB3BBBBBBBBBBBBB�B�B�B�B�B�B�B�B�B�B�B
�B
�B	�B	�B�B�B�B�B�B�B�B�B�B�B�B �B��B��B�}B�}B�}B�vB�vB�pB�jB�dB�dB�eB�_B�YB�YB�YB�SB�SB�LB�FB�@B�:B�:B�4B�4B�4B�4B�4B�4B�.B�(B�!B�!B�!B�B�B�!B�(B�(B�(B�4B�4B�4B�4B�.B�.B�(B�(B�!B�!B�B�B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BۻBڵBڵBڵBۻBۻBۻBۻBڵBٯBٯBةBעBՖBՖBՖB֜BՖBӊB�xB�lB�fB�`B�lB�fB�ZB�MB�MB�GB�GB�GB�GB�GB�GB�AB�AB�AB�AB�;B�;B�;B�;B�;B�AB�AB�AB�AB�AB�;B�6B�0B�*B�*B�#B�#B�#B�B�B�#B�*B�*B�*B�*B�*B�*B�*B�*B�*B�*B�*B�#B�#B�#B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�yB�yB�sB�mB�mB�gB�aB�aB�\B�\B�UB�UB�UB�UB�UB�OB�IB�IB�IB�IB�UB�UB�OB�IB�=B�7B�7B�1B�1B�1B�*B�*B�*B�*B�*B�$B�$B�*B�$B�$B�B�$B�$B�$B�*B�*B�*B�*B�$B�$B�*B�*B�*B�*B�=B�OB�bB�yB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�B�B�B�B�B�B�B�B�B�B�B�+B�7B�7B�<B�<B�BB�HB�HB�NB�NB�NB�NB�NB�HB�BB�=B�=B�+B�$B�$B�$B�$B�$B�$B�$B�$B�$B�+B�+B�1B�7B�=B�=B�7B�+B�$B�$B�B�B�B�$B�$B�$B�$B�$B�$B�$B�$B�$B�B�B�$B�+B�$B�$B�$B�$B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B��B��B��B��B��B�B��B��B��B��B�B��B��B�B�B��B��B�B�B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0074152                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904502018110709045020181107090450  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172717  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172717  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090450  IP  PSAL            @��D��3G�O�                