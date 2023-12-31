CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:16Z creation; 2018-11-05T17:27:06Z last update (coriolis COQC software)   
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
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172616  20181107090438  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�P����1   @�P����@N���L�'�B�����8   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@Fff@�33@�33@�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A���A���A���A���A���A�  A�33A�33A�  A�  A�  B   B  B��B  B  B  B  B  B   B$ffB(  B,  B0  B3��B8  B<  B@  BD  BH  BL  BO��BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B{��B�  B�  B�  B�  B�33B�33B�  B���B���B�  B�  B�  B���B���B���B�  B�33B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B���B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  C  C� C  C	ffC�fC� C�C��C�C� C�fCffC�fC"� C%  C'� C*  C,� C.�fC1� C4  C6ffC9  C;��C>  C@� CC  CEffCH  CJ� CM  COffCR  CT� CW  CY� C\  C^� Ca  Cc��Cf�ChffCj�fCm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C��3C�@ Cό�C�� C�  C�@ CԀ C�� C�  C�@ Cٌ�C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C��3C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D�fD  D@ D� D� D   D!@ D"�fD#� D$��D&@ D'� D(� D*  D+FfD,�fD-�fD/fD0FfD1� D2� D4  D5@ D6y�D7� D9fD:@ D;� D<� D>  D?9�D@� DA�fDC  DD@ DE�fDF� DH  DI@ DJ� DK� DM  DNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DYy�DZ��D\  D]@ D^� D_��Da  Db@ Dc� Dd�fDffDgFfDh� Di� Dk  Dl@ Dm�fDn�fDp  Dq@ Dr� Ds��Du  DvFfDw� Dx� Dz  D{@ D|� D}��D~��D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D���D�|�D��D���D�\�D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D��D���D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�C3D��3D�� D��D�� D�` D�  D�� D�<�D�� D�� D�#3D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�\�D�  D͠ D�C3D�� Dπ D�  D�� D�` D�  DҠ D�@ D��3Dԃ3D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D��3D�` D���Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D��D�@ D�� D� D�#3D�� D�` D�  D�3D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D��3D�FfD�ٚ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@Fff@�33@�33@�  @�  A   A  A   A0  A@  AP  A`  Ap  A�  A�  A�  A�  A�  A���A���A���A���A���A�  A�33A�33A�  A�  A�  B   B  B��B  B  B  B  B  B   B$ffB(  B,  B0  B3��B8  B<  B@  BD  BH  BL  BO��BT  BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B{��B�  B�  B�  B�  B�33B�33B�  B���B���B�  B�  B�  B���B���B���B�  B�33B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B���B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  C  C� C  C	ffC�fC� C�C��C�C� C�fCffC�fC"� C%  C'� C*  C,� C.�fC1� C4  C6ffC9  C;��C>  C@� CC  CEffCH  CJ� CM  COffCR  CT� CW  CY� C\  C^� Ca  Cc��Cf�ChffCj�fCm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C���C�  C�@ C�� C���C�  C�@ C���C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�L�C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�@ Cŀ CƳ3C�  C�@ Cʀ C�� C��3C�@ Cό�C�� C�  C�@ CԀ C�� C�  C�@ Cٌ�C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�� C��3C�@ C�� C�� C�  C�� C�  D � D  D@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� D  D@ D� D�fD  D@ D� D� D   D!@ D"�fD#� D$��D&@ D'� D(� D*  D+FfD,�fD-�fD/fD0FfD1� D2� D4  D5@ D6y�D7� D9fD:@ D;� D<� D>  D?9�D@� DA�fDC  DD@ DE�fDF� DH  DI@ DJ� DK� DM  DNFfDO� DP� DR  DS@ DT� DU� DW  DX@ DYy�DZ��D\  D]@ D^� D_��Da  Db@ Dc� Dd�fDffDgFfDh� Di� Dk  Dl@ Dm�fDn�fDp  Dq@ Dr� Ds��Du  DvFfDw� Dx� Dz  D{@ D|� D}��D~��D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D���D�|�D��D���D�\�D�  D�� D�@ D�� D�� D�#3D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D��D���D�` D�  D�� D�@ D�� D�� D�  D���D�` D�  D�� D�C3D��3D�� D��D�� D�` D�  D�� D�<�D�� D�� D�#3D�� D�` D�3D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�  D���D�` D�  Dà D�@ D�� Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D�  D�� D�\�D�  D͠ D�C3D�� Dπ D�  D�� D�` D�  DҠ D�@ D��3Dԃ3D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  D��3D�` D���Dܠ D�@ D�� Dހ D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D��D�@ D�� D� D�#3D�� D�` D�  D�3D�@ D�� D� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�c3D�  D��3D�FfD�ٚ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�n�@�~�@Ƈ+@Ƈ+@Ƈ+@�n�@�n�@�v�@�~�@Ƈ+@�~�@�~�@�~�@Ə\@Ə\@Ə\@Ƈ+@Ɨ�@ư!@Ƨ�@Ƨ�@Ƨ�@ư!@Ɨ�@Ƨ�@ư!@Ɨ�@ư!@ƸR@���@���@���@���@���@���@���@���@�ȴ@���@���@���@���@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@���@��@��@��@��H@��H@��H@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@���@ƸR@���@ƸR@ƸR@Ƨ�@�-@öF@�J@�?}@��@�K�@�@���@��@�I�@�b@��m@��
@�ƨ@��
@�  @��;@� �@��@�5?@���@��j@��F@�"�@��@�x�@��@��;@���@�$�@�^5@��@��h@�X@���@�G�@�S�@��R@�v�@�E�@�@�p�@�Ĝ@�j@�(�@�ƨ@�K�@�S�@�"�@���@��@��`@��D@��@��P@�l�@�C�@���@�{@�p�@�V@�Ĝ@� �@��;@���@�\)@�C�@�33@�ȴ@��@���@�`B@�7L@���@�z�@�A�@���@���@���@�dZ@�@���@�ff@�J@�@�hs@�%@���@��@��@��D@�  @��m@���@�dZ@�C�@�33@�o@��!@�n�@�5?@���@���@���@�O�@�G�@�?}@��@�%@���@�j@�1'@���@���@�t�@�|�@�\)@�33@�
=@��@���@��+@�^5@�5?@��T@�x�@�?}@�/@�%@��j@�I�@�9X@�  @��w@���@�;d@�+@�K�@�33@�"�@��@���@�~�@��\@��+@��\@���@�E�@�$�@�-@�$�@�{@�@�@�p�@�O�@�hs@�p�@���@���@�{@�=q@�M�@�5?@�$�@��@��@�V@�n�@�5?@�{@��@��@�J@��@���@�x�@�?}@�Ĝ@���@��D@�r�@��@��@��;@��
@��@���@��H@�@���@�~�@�V@�-@��#@��h@��7@�O�@���@���@��D@�9X@��@~�R@}@}`B@}/@|�/@|j@{�m@{S�@{"�@z��@y�#@yX@xA�@w��@w
=@v�R@v�R@vff@u�@u�h@t��@u`B@uO�@u�@t�@t9X@s�
@sƨ@s��@sC�@r�@r~�@r^5@q�7@q&�@p��@q&�@qx�@qx�@qhs@qx�@qx�@qhs@qhs@q7L@q%@p��@q%@qG�@qX@qx�@qX@qG�@qX@q��@q�^@q��@q�@q��@rn�@rn�@rn�@r^5@r-@r��@r�\@r~�@r�!@s"�@s"�@s��@t(�@tz�@t��@t��@t�j@t��@t�/@u�@t��@t�@t��@uO�@up�@u�h@u��@u�@uO�@t��@t�D@t�@s��@sC�@s@r�!@r-@r^5@r=q@q��@q7L@q7L@qG�@p��@p�9@pĜ@p��@p�u@pbN@pA�@pb@p1'@p �@o�@ol�@o+@o
=@o
=@o
=@n�y@nff@m�T@m��@mp�@m?}@m�@l�/@l�@l��@l�D@lj@lZ@l�@k�m@kƨ@k��@kdZ@j�H@j��@j�!@j��@jM�@jn�@j�\@j�H@j��@j�\@j^5@j-@i�@i�^@i�7@h��@h�@h1'@g�@g�w@g|�@g;d@f��@f5?@f@e�-@e`B@e?}@d�/@d��@d�@c�
@c�
@c�m@c�m@cƨ@c��@cdZ@cC�@c33@b�H@b^5@a�@a�7@ahs@ahs@aG�@`�`@`�9@`�@`A�@`b@_��@_;d@^ȴ@^�R@^��@^E�@]@]p�@]?}@\��@\��@\�D@\(�@[�
@[�@[dZ@[C�@[o@[@Z��@Z�!@Z~�@Z^5@Z=q@ZJ@Y�#@Y��@Yhs@X��@X�9@X�u@Xb@W��@W��@WK�@W�@Vȴ@Vff@U�@T�/@T��@T��@T�D@T(�@S�
@S�F@S�
@S��@SdZ@R�@R��@Rn�@RM�@R�@RJ@Q�#@QX@P�9@P�u@PbN@Pb@O��@O�w@O�@Nȴ@NV@N$�@M@M�-@M`B@L�/@L��@LZ@K��@K�m@L(�@K�F@K�
@K�
@K33@K33@K@J�!@Jn�@J=q@J�@I��@I�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�n�@�~�@Ƈ+@Ƈ+@Ƈ+@�n�@�n�@�v�@�~�@Ƈ+@�~�@�~�@�~�@Ə\@Ə\@Ə\@Ƈ+@Ɨ�@ư!@Ƨ�@Ƨ�@Ƨ�@ư!@Ɨ�@Ƨ�@ư!@Ɨ�@ư!@ƸR@���@���@���@���@���@���@���@���@�ȴ@���@���@���@���@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@���@��@��@��@��H@��H@��H@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@��y@���@ƸR@���@ƸR@ƸR@Ƨ�@�-@öF@�J@�?}@��@�K�@�@���@��@�I�@�b@��m@��
@�ƨ@��
@�  @��;@� �@��@�5?@���@��j@��F@�"�@��@�x�@��@��;@���@�$�@�^5@��@��h@�X@���@�G�@�S�@��R@�v�@�E�@�@�p�@�Ĝ@�j@�(�@�ƨ@�K�@�S�@�"�@���@��@��`@��D@��@��P@�l�@�C�@���@�{@�p�@�V@�Ĝ@� �@��;@���@�\)@�C�@�33@�ȴ@��@���@�`B@�7L@���@�z�@�A�@���@���@���@�dZ@�@���@�ff@�J@�@�hs@�%@���@��@��@��D@�  @��m@���@�dZ@�C�@�33@�o@��!@�n�@�5?@���@���@���@�O�@�G�@�?}@��@�%@���@�j@�1'@���@���@�t�@�|�@�\)@�33@�
=@��@���@��+@�^5@�5?@��T@�x�@�?}@�/@�%@��j@�I�@�9X@�  @��w@���@�;d@�+@�K�@�33@�"�@��@���@�~�@��\@��+@��\@���@�E�@�$�@�-@�$�@�{@�@�@�p�@�O�@�hs@�p�@���@���@�{@�=q@�M�@�5?@�$�@��@��@�V@�n�@�5?@�{@��@��@�J@��@���@�x�@�?}@�Ĝ@���@��D@�r�@��@��@��;@��
@��@���@��H@�@���@�~�@�V@�-@��#@��h@��7@�O�@���@���@��D@�9X@��@~�R@}@}`B@}/@|�/@|j@{�m@{S�@{"�@z��@y�#@yX@xA�@w��@w
=@v�R@v�R@vff@u�@u�h@t��@u`B@uO�@u�@t�@t9X@s�
@sƨ@s��@sC�@r�@r~�@r^5@q�7@q&�@p��@q&�@qx�@qx�@qhs@qx�@qx�@qhs@qhs@q7L@q%@p��@q%@qG�@qX@qx�@qX@qG�@qX@q��@q�^@q��@q�@q��@rn�@rn�@rn�@r^5@r-@r��@r�\@r~�@r�!@s"�@s"�@s��@t(�@tz�@t��@t��@t�j@t��@t�/@u�@t��@t�@t��@uO�@up�@u�h@u��@u�@uO�@t��@t�D@t�@s��@sC�@s@r�!@r-@r^5@r=q@q��@q7L@q7L@qG�@p��@p�9@pĜ@p��@p�u@pbN@pA�@pb@p1'@p �@o�@ol�@o+@o
=@o
=@o
=@n�y@nff@m�T@m��@mp�@m?}@m�@l�/@l�@l��@l�D@lj@lZ@l�@k�m@kƨ@k��@kdZ@j�H@j��@j�!@j��@jM�@jn�@j�\@j�H@j��@j�\@j^5@j-@i�@i�^@i�7@h��@h�@h1'@g�@g�w@g|�@g;d@f��@f5?@f@e�-@e`B@e?}@d�/@d��@d�@c�
@c�
@c�m@c�m@cƨ@c��@cdZ@cC�@c33@b�H@b^5@a�@a�7@ahs@ahs@aG�@`�`@`�9@`�@`A�@`b@_��@_;d@^ȴ@^�R@^��@^E�@]@]p�@]?}@\��@\��@\�D@\(�@[�
@[�@[dZ@[C�@[o@[@Z��@Z�!@Z~�@Z^5@Z=q@ZJ@Y�#@Y��@Yhs@X��@X�9@X�u@Xb@W��@W��@WK�@W�@Vȴ@Vff@U�@T�/@T��@T��@T�D@T(�@S�
@S�F@S�
@S��@SdZ@R�@R��@Rn�@RM�@R�@RJ@Q�#@QX@P�9@P�u@PbN@Pb@O��@O�w@O�@Nȴ@NV@N$�@M@M�-@M`B@L�/@L��@LZ@K��@K�m@L(�@K�F@K�
@K�
@K33@K33@K@J�!@Jn�@J=q@J�@I��@I�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�bB�bB�bB�bB�bB�\B�\B�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�hB�hB�hB�bB�hB�hB�bB�hB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�\B�\B�VB�PB�DB�7B�+B�B�7B��B�jB��BĜBĜBÖBÖB��B��B��BBBBÖBĜBŢBÖBB��B��B��B��BĜB��B��B��BɺBȴB��BɺBɺBɺBɺB��BŢBĜBĜBĜBÖBÖBBB��B��B��B��B��B�}B�wB�wB�wB�wB�wB�}B��B��B�}B�qB�jB�dB�^B�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�RB�RB�RB�LB�LB�RB�RB�RB�LB�LB�FB�FB�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�9B�3B�3B�-B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�oB�hB�bB�bB�bB�bB�\B�VB�VB�\B�\B�\B�VB�VB�PB�PB�PB�JB�JB�DB�JB�DB�=B�=B�=B�DB�DB�JB�JB�PB�PB�PB�PB�PB�PB�VB�\B�\B�\B�bB�bB�hB�uB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�!B�'B�'B�'B�'B�'B�'B�'B�!B�'B�!B�!B�'B�-B�!B�!B�!B�'B�'B�-B�-B�-B�-B�-B�3B�3B�9B�?B�9B�9B�?B�?B�?B�FB�FB�?B�?B�?B�?B�?B�FB�FB�FB�LB�LB�LB�LB�LB�LB�LB�RB�LB�LB�RB�RB�RB�XB�XB�dB�jB�jB�jB�jB�jB�jB�qB�jB�jB�jB�jB�jB�jB�jB�jB�dB�dB�dB�dB�^B�^B�^B�^B�^B�^B�^B�^B�^B�dB�^B�^B�dB�dB�dB�dB�^B�^B�^B�^B�^B�^B�XB�XB�XB�^B�^B�^B�^B�^B�^B�XB�XB�XB�^B�^B�^B�^B�dB�dB�^B�dB�dB�dB�dB�dB�dB�dB�dB�dB�dB�^B�dB�dB�^B�dB�dB�dB�^B�^B�^B�^B�^B�XB�XB�XB�XB�XB�XB�XB�XB�XB�^B�^B�^B�^B�^B�^B�^B�^B�^B�^B�XB�XB�XB�XB�XB�XB�XB�XB�RB�RB�XB�XB�XB�XB�RB�XB�RB�RB�XB�^B�^B�^B�^B�XB�^B�^B�^B�^B�^B�^B�^B�d1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�sB�gB�UB�sB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�oB�oB�iB�iB�iB�iB�iB�cB�cB�cB�cB�]B�]B�WB�WB�PB�PB�PB�PB�PB�JB�JB�JB�JB�DB�DB�DB�>B�>B�8B�8B�2B�2B�,B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�%B�B�B�%B�%B�B�B�B�B�B�B�B�%B�2B�8B�>B�DB�>B�>B�>B�>B�DB�JB�DB�DB�JB�JB�JB�DB�DB�>B�>B�8B�2B�2B�2B�2B�,B�,B�2B�2B�,B�B�%B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB�zB�zB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�%B�,B�,B�8B�8B�>B�DB�PB�WB�]B�cB�cB�cB�cB�cB�cB�cB�]B�cB�]B�]B�cB�iB�]B�]B�]B�cB�cB�iB�iB�iB�iB�iB�oB�oB�uB�{B�uB�uB�{B�{B�{B��B��B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0017238                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904392018110709043920181107090439  IF  ARFMCODA024c                                                                20181105172616                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172706  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172706  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090439  IP  PSAL            @ffD�ٚG�O�                