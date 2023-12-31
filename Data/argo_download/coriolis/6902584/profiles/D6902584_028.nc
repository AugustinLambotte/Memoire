CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  1   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2018-11-05T17:26:17Z creation; 2018-11-05T17:27:14Z last update (coriolis COQC software)   
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
_FillValue                 4  B`   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  MX   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  O�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  XP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  a   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  cH   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  l   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  n@   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  w   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 4  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �$   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �d   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �t   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �x   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��             ,  ��Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090446  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�?�R�H]1   @�?��[�$@Ow�ڔi�A�����o2   IRIDIUM A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @��@Fff@�  @�33@�  @�33A   A��A   A0  A@  AQ��A`  Ap  A�  A�  A�  A�  A�  A�  A���A���A���A�  A�  A�  A���A���A�  A���B   B��B��B��B��B  B  B  B ffB$ffB(ffB,ffB0  B4  B7��B;��B@  BDffBHffBL  BP  BS��BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�33B�  B���B���B�  B�  B�  B�  B�  B���B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  C�C��C�C	� C  C� C�C� C  C� C  C��C �C"� C%  C'��C*  C,� C/  C1� C4  C6��C9  C;��C>�C@��CC  CEffCG�fCJffCM  CO� CR  CT� CV�fCYffC\  C^� Ca�Cc��Cf  Ch� Cj�fCm� Cp�Cr� Ct�fCwffCy�fC|� C  C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�33C�� C���C��C�@ C�� C�� C�  C�33C�s3C�� C��C�L�C���C���C��C�@ C���C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C��C�@ Cʀ C���C�  C�@ Cπ C�� C�  C�@ CԀ C�� C��C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C�� C��3C�  C�� C��3D � D  D@ D�fD� D  D@ D	y�D
� D  D@ D� D� DfD@ D� D� D  D@ D� D� D  D@ D�fD� D   D!FfD"� D#�fD%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DEy�DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU��DV��DX9�DY� DZ�fD\fD]FfD^�fD_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDpfDq@ Dr� Ds� Du  DvFfDw� Dx� DzfD{FfD|�fD}�fD  D�  D�� D�\�D�  D�� D�C3D��3D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�3D�� D�<�D�� D�� D�#3D�� D�` D�  D�� D�C3D�� D�� D�#3D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D��3D�c3D�  D���D�C3D�� D�� D�  D�� D�\�D�  D�� D�C3D��3D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D���D�� D�  D���D�` D�3D��3D�@ D���D�� D�  D��3D�` D�  Dà D�@ D���Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D��D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D���Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D��D�� D�` D���D� D�@ D�� D� D�#3D�� D�` D�  D� D�@ D�� D�3D�#3D��3D�` D���D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D���D�� D�C3D�� D�|�D��D�� D�\�D�3D��3D�c3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@Fff@�  @�33@�  @�33A   A��A   A0  A@  AQ��A`  Ap  A�  A�  A�  A�  A�  A�  A���A���A���A�  A�  A�  A���A���A�  A���B   B��B��B��B��B  B  B  B ffB$ffB(ffB,ffB0  B4  B7��B;��B@  BDffBHffBL  BP  BS��BX  B\  B`  Bd  Bh  Bl  Bp  Bt  Bx  B|  B�  B�  B�  B�  B�33B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�  B���B�  B�  B�  B�33B�  B���B���B�  B�  B�  B�  B�  B���B�  B�  B�  B�33B�33B�33B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  B�  C�C��C�C	� C  C� C�C� C  C� C  C��C �C"� C%  C'��C*  C,� C/  C1� C4  C6��C9  C;��C>�C@��CC  CEffCG�fCJffCM  CO� CR  CT� CV�fCYffC\  C^� Ca�Cc��Cf  Ch� Cj�fCm� Cp�Cr� Ct�fCwffCy�fC|� C  C���C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C��3C�33C�� C���C��C�@ C�� C�� C�  C�33C�s3C�� C��C�L�C���C���C��C�@ C���C�� C��3C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ Cŀ C�� C��C�@ Cʀ C���C�  C�@ Cπ C�� C�  C�@ CԀ C�� C��C�@ Cـ C�� C�  C�L�Cހ C�� C�  C�@ C� C�� C�  C�@ C� C�� C�  C�L�C� C�� C�  C�@ C� C�� C�  C�@ C�� C��3C�  C�� C��3D � D  D@ D�fD� D  D@ D	y�D
� D  D@ D� D� DfD@ D� D� D  D@ D� D� D  D@ D�fD� D   D!FfD"� D#�fD%  D&@ D'� D(� D*  D+@ D,� D-� D/  D0@ D1� D2� D4  D5@ D6� D7� D9  D:@ D;� D<� D>  D?@ D@� DA��DC  DD@ DEy�DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU��DV��DX9�DY� DZ�fD\fD]FfD^�fD_� Da  Db@ Dc� Dd� Df  Dg@ Dh� Di� Dk  Dl@ Dm� Dn�fDpfDq@ Dr� Ds� Du  DvFfDw� Dx� DzfD{FfD|�fD}�fD  D�  D�� D�\�D�  D�� D�C3D��3D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�3D�� D�<�D�� D�� D�#3D�� D�` D�  D�� D�C3D�� D�� D�#3D�� D�` D�  D�� D�C3D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D���D���D�@ D�� D�� D�  D��3D�c3D�  D���D�C3D�� D�� D�  D�� D�\�D�  D�� D�C3D��3D�� D�  D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D���D�� D�  D���D�` D�3D��3D�@ D���D�� D�  D��3D�` D�  Dà D�@ D���Dŀ D�  D�� D�` D�  DȠ D�@ D�� Dʀ D��D�� D�` D�  D͠ D�@ D�� Dπ D�  D�� D�` D�  DҠ D�@ D�� DԀ D�  D�� D�` D���Dנ D�@ D�� Dـ D�  D�� D�` D�  Dܠ D�@ D�� Dހ D��D�� D�` D���D� D�@ D�� D� D�#3D�� D�` D�  D� D�@ D�� D�3D�#3D��3D�` D���D� D�@ D�� D� D�#3D��3D�` D�  D� D�@ D�� D� D�  D�� D�` D���D�� D�C3D�� D�|�D��D�� D�\�D�3D��3D�c3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�X@�O�@��7@���@���@��T@��T@���@��@���@���@�@�=q@��@���@��T@��@�@�@��#@��^@���@�J@���@�5?@�ȴ@���@�
=@�o@�dZ@���@�;d@�dZ@�;d@�"�@���@��@�"�@���@�@��@�"�@�+@�;d@�33@�"�@��H@��@��!@��!@�;d@�K�@�|�@�t�@��@��F@�l�@�t�@�l�@�|�@��;@��;@���@�dZ@�33@�+@�o@��@�"�@��@�"�@�"�@�+@�"�@�"�@��@�o@��@��@�"�@�"�@�C�@�C�@�C�@�C�@�K�@�C�@�K�@�C�@�;d@�33@�;d@�C�@�;d@�;d@�C�@�C�@�K�@�K�@�+@��@�33@�+@���@��y@��H@��H@��@���@���@��H@�ȴ@���@��H@��y@�o@�@�"�@��@��y@���@�E�@��@�@��T@��@�@�{@�$�@�J@�@��-@�@���@��@���@�X@���@���@��j@��/@��`@��9@��j@��j@��j@��j@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@��`@��/@��/@��`@���@��/@��/@���@���@���@��@��@���@��u@���@���@���@��u@��u@��@��@��9@���@���@��@��@��@�r�@�r�@�I�@�A�@�I�@�Z@�Q�@�I�@�A�@�9X@�1'@�1'@�1'@�1'@�1'@�9X@�9X@�9X@�A�@�1'@�1'@�9X@�9X@�1'@�9X@�9X@�1'@� �@��@� �@��@��@��@�b@�b@��@��@��@�  @���@���@��@��m@��;@�ƨ@��@���@���@���@��@�|�@��@��P@��@�|�@�\)@�33@��@��@�o@�o@�o@��@��@���@��@��@��y@��y@��y@��@��!@���@��\@���@�V@�$�@�=q@�M�@�~�@���@��H@�ȴ@��R@��!@���@�V@�@���@�G�@���@���@�bN@�b@�  @��m@��@�$�@��@��T@��#@�@�x�@�V@�V@��@���@�Z@�A�@�@~�+@~�R@\)@�P@;d@;d@;d@+@~ȴ@~V@~5?@~{@}��@}V@|j@|Z@{�m@{C�@{@z-@y�@y�7@yx�@yx�@y&�@wl�@vȴ@v5?@up�@tZ@t1@t(�@t9X@s��@r�\@r�\@so@sC�@s@r��@r��@r~�@q��@q&�@q�@q�@p�`@p1'@ol�@o�@n�+@n$�@m��@m@m@m��@m��@m�T@n@n$�@n5?@nv�@n�y@n��@o+@o\)@o�@p  @pb@pb@pb@p �@p��@q%@qG�@rM�@s��@t�@s��@rM�@sdZ@sS�@s"�@so@r�\@r�H@st�@s��@t�@t1@tz�@t�D@sƨ@st�@sS�@sS�@r�@s@s"�@s��@s�m@s�m@t1@t9X@uV@u�@u��@t�@t��@t�j@t��@s@r�@r��@r�\@rn�@r�\@rn�@r~�@rn�@r=q@q��@q�@q��@qx�@qX@qG�@q&�@p��@q�@q&�@q&�@q7L@qX@qG�@p��@pr�@p �@p �@pb@o��@o�w@o|�@o�@n�+@nV@n{@m�T@m�@m�-@m�h@m@m�-@m�@mp�@m�-@m�T@n{@n5?@n5?@n5?@n5?@n5?@m��@m�h@l��@k@j��@jn�@j�@jJ@i�7@h��@hbN@g�;@g\)@g+@f��@fȴ@f��@f5?@e�@fE�@f��@fE�@e�@e��@e?}@d�@d1@c��@cƨ@c��@d1@cƨ@cS�@b�@b�\@b^5@b^5@bJ@a�^@a��@aX@`�`@`r�@`bN@`  @_�w@_��@_|�@_
=@^�@^ȴ@^E�@]@]V@\�j@\�@\9X@\1@[��@[33@[@Zn�@Z=q@Z�@Y�#@Yx�@YG�@X��@X  @W�;@W�w@W��@W|�@W+@W+@V�y@Vȴ@V��@Vv�@V$�@V@U��@U?}@U�@T��@T�D@T1@S��@S��@S��@S33@S"�@R��@R�H@R��@R�!@R�\@R-@Q�@Q��@Q�7@QX@Q�@PĜ@PbN@P1'@P �@O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�X@�O�@��7@���@���@��T@��T@���@��@���@���@�@�=q@��@���@��T@��@�@�@��#@��^@���@�J@���@�5?@�ȴ@���@�
=@�o@�dZ@���@�;d@�dZ@�;d@�"�@���@��@�"�@���@�@��@�"�@�+@�;d@�33@�"�@��H@��@��!@��!@�;d@�K�@�|�@�t�@��@��F@�l�@�t�@�l�@�|�@��;@��;@���@�dZ@�33@�+@�o@��@�"�@��@�"�@�"�@�+@�"�@�"�@��@�o@��@��@�"�@�"�@�C�@�C�@�C�@�C�@�K�@�C�@�K�@�C�@�;d@�33@�;d@�C�@�;d@�;d@�C�@�C�@�K�@�K�@�+@��@�33@�+@���@��y@��H@��H@��@���@���@��H@�ȴ@���@��H@��y@�o@�@�"�@��@��y@���@�E�@��@�@��T@��@�@�{@�$�@�J@�@��-@�@���@��@���@�X@���@���@��j@��/@��`@��9@��j@��j@��j@��j@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@�Ĝ@��`@��/@��/@��`@���@��/@��/@���@���@���@��@��@���@��u@���@���@���@��u@��u@��@��@��9@���@���@��@��@��@�r�@�r�@�I�@�A�@�I�@�Z@�Q�@�I�@�A�@�9X@�1'@�1'@�1'@�1'@�1'@�9X@�9X@�9X@�A�@�1'@�1'@�9X@�9X@�1'@�9X@�9X@�1'@� �@��@� �@��@��@��@�b@�b@��@��@��@�  @���@���@��@��m@��;@�ƨ@��@���@���@���@��@�|�@��@��P@��@�|�@�\)@�33@��@��@�o@�o@�o@��@��@���@��@��@��y@��y@��y@��@��!@���@��\@���@�V@�$�@�=q@�M�@�~�@���@��H@�ȴ@��R@��!@���@�V@�@���@�G�@���@���@�bN@�b@�  @��m@��@�$�@��@��T@��#@�@�x�@�V@�V@��@���@�Z@�A�@�@~�+@~�R@\)@�P@;d@;d@;d@+@~ȴ@~V@~5?@~{@}��@}V@|j@|Z@{�m@{C�@{@z-@y�@y�7@yx�@yx�@y&�@wl�@vȴ@v5?@up�@tZ@t1@t(�@t9X@s��@r�\@r�\@so@sC�@s@r��@r��@r~�@q��@q&�@q�@q�@p�`@p1'@ol�@o�@n�+@n$�@m��@m@m@m��@m��@m�T@n@n$�@n5?@nv�@n�y@n��@o+@o\)@o�@p  @pb@pb@pb@p �@p��@q%@qG�@rM�@s��@t�@s��@rM�@sdZ@sS�@s"�@so@r�\@r�H@st�@s��@t�@t1@tz�@t�D@sƨ@st�@sS�@sS�@r�@s@s"�@s��@s�m@s�m@t1@t9X@uV@u�@u��@t�@t��@t�j@t��@s@r�@r��@r�\@rn�@r�\@rn�@r~�@rn�@r=q@q��@q�@q��@qx�@qX@qG�@q&�@p��@q�@q&�@q&�@q7L@qX@qG�@p��@pr�@p �@p �@pb@o��@o�w@o|�@o�@n�+@nV@n{@m�T@m�@m�-@m�h@m@m�-@m�@mp�@m�-@m�T@n{@n5?@n5?@n5?@n5?@n5?@m��@m�h@l��@k@j��@jn�@j�@jJ@i�7@h��@hbN@g�;@g\)@g+@f��@fȴ@f��@f5?@e�@fE�@f��@fE�@e�@e��@e?}@d�@d1@c��@cƨ@c��@d1@cƨ@cS�@b�@b�\@b^5@b^5@bJ@a�^@a��@aX@`�`@`r�@`bN@`  @_�w@_��@_|�@_
=@^�@^ȴ@^E�@]@]V@\�j@\�@\9X@\1@[��@[33@[@Zn�@Z=q@Z�@Y�#@Yx�@YG�@X��@X  @W�;@W�w@W��@W|�@W+@W+@V�y@Vȴ@V��@Vv�@V$�@V@U��@U?}@U�@T��@T�D@T1@S��@S��@S��@S33@S"�@R��@R�H@R��@R�!@R�\@R-@Q�@Q��@Q�7@QX@Q�@PĜ@PbN@P1'@P �@O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B��B��B��B��B��BɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBȴBȴBȴBȴBȴBȴBȴBȴBȴBǮBǮBȴBǮBƨBŢBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBBBBBBBBBBBBBBBBBBBBBBBBBBBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�wB�wB�wB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�jB�dB�dB�dB�jB�qB�}B�}B�wB�wB�wB�qB�jB�dB�^B�RB�LB�FB�?B�9B�9B�3B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�oB�hB�\B�\B�\B�\B�\B�JB�PB�VB�\B�\B�\B�\B�VB�PB�JB�JB�JB�JB�DB�7B�7B�+B�+B�%B�%B�+B�+B�1B�1B�7B�7B�=B�DB�PB�VB�VB�\B�bB�oB�uB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�!B�'B�3B�FB�XB�XB�LB�LB�LB�RB�?B�?B�?B�?B�FB�LB�RB�XB�^B�^B�^B�^B�^B�^B�dB�dB�dB�dB�jB�qB�wB�wB��B��B��B�}B�}B�}B�}B��B��B��B��B��B��B��B��B��B��B��BBÖBBÖBŢBǮBǮBȴBɺB��BɺB��BɺBȴBǮBĜBĜBĜBÖBĜBÖBBBB��B��B��B��B��B��B��BÖBĜBĜBĜBĜBĜBÖBBBBÖBÖBÖBÖBBBÖBÖBÖBÖBÖBÖBÖBBBBBBBBBBB��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�wB�wB�wB�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�wB�wB�qB�wB�wB�wB�}B�}B�}B��B��B��B��B�}B��B��B��B��B�}B�}B�}B�}B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�4B�4B�4B�9B�9B�9B�9B�9B�9B�9B�9B�9B�-B�9B�3B�9B�3B�9B�9B�9B�3B�3B�?B�9B�?B�?B�?B�KB�EB�EB�KB�EB�EB�KB�EB�?B�EB�KB�EB�?B�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�EB�?B�QB�QB�EB�KB�KB�KB�KB�EB�KB�KB�EB�?B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�9B�?B�?B�EB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�9B�9B�9B�9B�3B�-B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�'B�!B�!B�'B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�	B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�vB�vB�vB�vB�jB�dB�dB�eB�^B�XB�RB�@B�@B�FB�LB�RB�RB�RB�RB�RB�LB�FB�FB�FB�@B�:B�-B�-B�'B�!B�!B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�	B�B�3B�FB�:B�-B�@B�@B�@B�@B�@B�LB�XB�^B�eB�eB�qB�qB�qB�qB�qB�qB�qB�wB�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�
B�B�
B�B�"B�"B�(B�.B�4B�.B�4B�.B�(B�"B�B�B�B�
B�B�
B�B�B�B��B��B��B��B��B��B��B�
B�B�B�B�B�B�
B�B�B�B�
B�
B�
B�
B�B�B�
B�
B�
B�
B�
B�
B�
B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.005418                                                                                                                                                                                                              No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904462018110709044620181107090446  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172714  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.5                                                                 20181105172714  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090446  IP  PSAL            @��D�c3G�O�                