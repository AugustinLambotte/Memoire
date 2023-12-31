CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  2   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
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
_FillValue                  ,  �             ,  �Argo profile    3.1 1.2 19500101000000  20181105172617  20181107090447  6902584 BSH                                                             Birgit KLEIN                                                    PRES            TEMP            PSAL               A   IF                                  2C  D   NOVA                            SN143                           n/a                             865 @�D�F)1   @�Dѩ�|X@O�ũ�܈�@m���h1   GPS     A   A   A   Primary sampling: averaged [10sec sampling;50cbar interval,50cbar average from 20000cbar to 5000cbar;25cbar interval,25cbar average from 5000cbar to 1000cbar;10cbar interval,10cbar average from 1000cbar to 20cbar;transition bins not included]                 @ff@@  @�  @�  @�  @�  A   A��A   A0  A@  AP  A`  Ap  A�  A�  A���A�  A�  A�  A�  A�33A�  A���A���A�  A�33A�33A�  A�  B   B  B  B  B��B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BW��B[��B_��Bd  Bh  BlffBpffBt  Bx  B|  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�33B�  B�  B�  B�  B���B�  B�  B���B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'ffC)�fC,ffC.�fC1� C4�C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  COffCR  CT��CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�s3C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�33Cŀ C�� C�  C�@ C�s3C�� C��C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��3C�@ Cތ�C�� C�  C�@ C� C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�3C�  C�@ C�� C��3C��3C�s3C�  D � DfD@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� DfD@ D� D� D  D@ D� D��D   D!@ D"� D#� D%  D&@ D'�fD(�fD*  D+@ D,� D-� D/  D09�D1� D2� D4  D5@ D6� D7� D9  D:@ D;�fD<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DV��DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  DgFfDh�fDi� Dk  Dl@ Dm�fDn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx�fDz  D{@ D|� D}� D  D��D���D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�  D���D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  Dà D�@ D��3Dŀ D�  D�� D�c3D�  DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dπ D�  D�� D�` D���DҜ�D�@ D�� D�|�D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dރ3D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D���D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D�� D�\�D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�I�D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ff@@  @�  @�  @�  @�  A   A��A   A0  A@  AP  A`  Ap  A�  A�  A���A�  A�  A�  A�  A�33A�  A���A���A�  A�33A�33A�  A�  B   B  B  B  B��B  B  B  B   B$  B(  B,  B0  B4  B8  B<  B@  BD  BHffBL  BP  BT  BW��B[��B_��Bd  Bh  BlffBpffBt  Bx  B|  B�33B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B���B���B�  B�  B�  B�  B�  B�33B�  B���B�  B�  B�33B�  B�  B�  B�  B���B�  B�  B���B�  B�  B�33B�  B�  B�  B�  B�  B�  B�33B�  B�  C  C� C  C	� C  C� C  C� C  C� C  C� C   C"� C%  C'ffC)�fC,ffC.�fC1� C4�C6��C9  C;� C>  C@� CC  CE� CH  CJ� CM  COffCR  CT��CW  CY� C\  C^� Ca  Cc� Cf  Ch� Ck  Cm� Cp  Cr� Cu  Cw� Cz  C|� C  C�� C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C�� C��3C�  C�@ C�� C�� C�  C�@ C�� C�� C��C�L�C�� C�� C�  C�@ C�s3C�� C�  C�@ C���C�� C�  C�@ C�s3C�� C�  C�@ C�� C���C�  C�@ C�� C�� C�  C�@ C�� C�� C�  C�@ C���C�� C�  C�33Cŀ C�� C�  C�@ C�s3C�� C��C�@ Cπ C�� C�  C�@ CԀ C�� C�  C�@ Cـ C�� C��3C�@ Cތ�C�� C�  C�@ C� C���C�  C�@ C� C�� C�  C�@ C� C�� C�  C�@ C� C�3C�  C�@ C�� C��3C��3C�s3C�  D � DfD@ D� D� D  D@ D	� D
� D  D@ D� D� D  D@ D� D� DfD@ D� D� D  D@ D� D��D   D!@ D"� D#� D%  D&@ D'�fD(�fD*  D+@ D,� D-� D/  D09�D1� D2� D4  D5@ D6� D7� D9  D:@ D;�fD<� D>  D?@ D@� DA��DC  DD@ DE� DF� DH  DI@ DJ� DK� DM  DN@ DO� DP� DR  DS@ DT� DU� DV��DX@ DY� DZ� D\  D]@ D^� D_� Da  Db@ Dc� Dd� Df  DgFfDh�fDi� Dk  Dl@ Dm�fDn� Dp  Dq@ Dr� Ds� Du  Dv@ Dw� Dx�fDz  D{@ D|� D}� D  D��D���D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D���D�� D�  D�� D�c3D�  D���D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D��3D�  D�� D�` D�  D�� D�@ D�� D�� D�  D��3D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D��3D�@ D�� D�� D�  D�� D�c3D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D��3D�#3D�� D�` D�  D���D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D�� D�� D�  D�� D�` D�  D�� D�@ D��3D�� D�  D�� D�` D�  D�� D�@ D��3D��3D�  D�� D�` D�  Dà D�@ D��3Dŀ D�  D�� D�c3D�  DȜ�D�@ D�� Dʀ D�  D�� D�` D�  D͠ D�@ D��3Dπ D�  D�� D�` D���DҜ�D�@ D�� D�|�D�  D�� D�` D�  Dנ D�@ D�� Dـ D�  Dڼ�D�` D�  Dܠ D�@ D�� Dރ3D�  D�� D�` D�  D� D�@ D��3D� D�  D�� D�` D���D� D�@ D�� D� D�  D�� D�` D�  D� D�@ D�� D� D�  D�� D�` D�  D� D�C3D�� D� D�  D�� D�\�D�  D��3D�@ D�� D�� D�  D�� D�` D�  D�� D�I�D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�z�@�j@�j@�j@�Z@�bN@�j@�z�@��@�r�@�r�@�j@�j@�bN@�r�@�r�@�Z@�A�@�A�@�I�@�A�@�I�@�I�@�A�@�1'@�1'@�1'@�9X@�1'@�1'@�I�@�9X@�A�@�1'@��@��@� �@�(�@�(�@�1'@�1'@�1'@�I�@�I�@�I�@�A�@�9X@�Q�@�9X@� �@�(�@�(�@� �@� �@��@� �@��@��@��@� �@�1'@�9X@�A�@�A�@�I�@�I�@�I�@�I�@�Q�@�Z@�Z@�bN@�bN@�j@�j@�j@�bN@�j@�(�@��F@�ƨ@���@��w@���@��@���@��@���@��w@���@�ƨ@��F@��w@��w@��F@��;@��
@���@��F@��F@��@�(�@� �@�1'@�1'@�9X@�(�@�1'@�  @��@��
@�I�@��@���@�|�@��@���@��@��m@��w@��P@�dZ@�"�@�-@���@���@��7@�`B@�O�@�O�@�&�@�%@��@���@��@��w@��@���@�n�@�M�@�5?@�{@�$�@��@�@��@�X@�V@�Ĝ@�z�@��
@���@�v�@�J@���@�bN@�  @���@�K�@�33@�
=@���@��R@��R@��!@���@��@�v�@��+@�V@�$�@���@��@��#@���@��#@�x�@�?}@��9@��D@��D@��@��@�r�@�Z@�I�@�9X@� �@��@�9X@�(�@��@��@�b@�1@�A�@�Q�@�bN@��/@�&�@�hs@�x�@�x�@��@�hs@�?}@�X@���@��@��@��@��@��@��@���@��@��T@��T@��#@���@�@��^@��^@��-@���@���@��h@��7@�x�@�`B@�`B@�?}@��@��@�%@�V@��@��@�%@��`@��/@���@��j@��9@��u@��D@��@�z�@�bN@�Z@�A�@�  @��F@��w@�|�@�\)@�o@��@��!@�ff@��\@���@���@���@���@�ff@�V@�V@�E�@���@�&�@��-@��D@�I�@�  @�+@���@�@��u@�1@��@�;d@��+@���@��T@���@���@���@��@��@��@�{@���@��@�%@�Ĝ@�I�@�ƨ@��@���@���@�dZ@�;d@���@���@�5?@�M�@�{@�&�@��j@��@�r�@�1'@��;@��F@��P@�|�@�dZ@�C�@��@���@�~�@�ff@�M�@�M�@�J@���@��@�X@�`B@�`B@�`B@�X@�G�@�%@��/@���@��9@��@�9X@�(�@�b@��@�P@~��@~$�@}@}?}@}�@}V@|�@|�/@|�j@|��@|(�@z�H@z�@y&�@w�w@w\)@w
=@w
=@vȴ@vE�@u�@u@vE�@v�@vV@uV@t�@tz�@t9X@tI�@u?}@v{@u��@v@v{@v$�@vV@vff@vV@v$�@v��@w;d@w�w@w�@xb@xb@xb@x  @x  @x  @x  @x �@xA�@y��@y�^@yG�@x��@xbN@x1'@xb@x �@w��@w�;@w�@xA�@xr�@xr�@xbN@w�;@w�@xb@wl�@w;d@v��@v��@v�y@v��@v��@v��@v��@v��@v��@v�+@v�+@v�+@v�+@v�+@vv�@v�+@v��@v�@v��@u��@vE�@vȴ@v5?@up�@uO�@t�@sƨ@r��@q�#@r-@q�#@r��@s��@r~�@q�@q��@q�#@r�@r=q@r=q@r�@q�#@q��@r�@q��@q�@q��@qhs@qX@p�`@p �@pb@p  @o�@o�@o|�@o+@n�y@nȴ@n�+@nE�@m�T@m`B@l�/@l9X@kƨ@k��@kS�@j�@j�!@j~�@i��@iG�@i7L@iG�@h��@h�@h1'@g��@g|�@g;d@g
=@f�R@fV@f$�@e@e�@e?}@e/@e�@eV@d�@d��@dj@d1@c�
@c��@ct�@b�@b��@bM�@a�@a�^@a&�@`Ĝ@`bN@` �@_��@_K�@_
=@^�y@^ȴ@^�@^v�@^E�@]�@]�@\��@\��@\�@\(�@[�m@[�F@[@Z�\@Z~�@Z�@Y�@Y��@Y�7@Yx�@YX@Y%@X�9@X�u@X �@W��@Wl�@W+@V��@Vȴ@V��@VV@V5?@U�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�z�@�j@�j@�j@�Z@�bN@�j@�z�@��@�r�@�r�@�j@�j@�bN@�r�@�r�@�Z@�A�@�A�@�I�@�A�@�I�@�I�@�A�@�1'@�1'@�1'@�9X@�1'@�1'@�I�@�9X@�A�@�1'@��@��@� �@�(�@�(�@�1'@�1'@�1'@�I�@�I�@�I�@�A�@�9X@�Q�@�9X@� �@�(�@�(�@� �@� �@��@� �@��@��@��@� �@�1'@�9X@�A�@�A�@�I�@�I�@�I�@�I�@�Q�@�Z@�Z@�bN@�bN@�j@�j@�j@�bN@�j@�(�@��F@�ƨ@���@��w@���@��@���@��@���@��w@���@�ƨ@��F@��w@��w@��F@��;@��
@���@��F@��F@��@�(�@� �@�1'@�1'@�9X@�(�@�1'@�  @��@��
@�I�@��@���@�|�@��@���@��@��m@��w@��P@�dZ@�"�@�-@���@���@��7@�`B@�O�@�O�@�&�@�%@��@���@��@��w@��@���@�n�@�M�@�5?@�{@�$�@��@�@��@�X@�V@�Ĝ@�z�@��
@���@�v�@�J@���@�bN@�  @���@�K�@�33@�
=@���@��R@��R@��!@���@��@�v�@��+@�V@�$�@���@��@��#@���@��#@�x�@�?}@��9@��D@��D@��@��@�r�@�Z@�I�@�9X@� �@��@�9X@�(�@��@��@�b@�1@�A�@�Q�@�bN@��/@�&�@�hs@�x�@�x�@��@�hs@�?}@�X@���@��@��@��@��@��@��@���@��@��T@��T@��#@���@�@��^@��^@��-@���@���@��h@��7@�x�@�`B@�`B@�?}@��@��@�%@�V@��@��@�%@��`@��/@���@��j@��9@��u@��D@��@�z�@�bN@�Z@�A�@�  @��F@��w@�|�@�\)@�o@��@��!@�ff@��\@���@���@���@���@�ff@�V@�V@�E�@���@�&�@��-@��D@�I�@�  @�+@���@�@��u@�1@��@�;d@��+@���@��T@���@���@���@��@��@��@�{@���@��@�%@�Ĝ@�I�@�ƨ@��@���@���@�dZ@�;d@���@���@�5?@�M�@�{@�&�@��j@��@�r�@�1'@��;@��F@��P@�|�@�dZ@�C�@��@���@�~�@�ff@�M�@�M�@�J@���@��@�X@�`B@�`B@�`B@�X@�G�@�%@��/@���@��9@��@�9X@�(�@�b@��@�P@~��@~$�@}@}?}@}�@}V@|�@|�/@|�j@|��@|(�@z�H@z�@y&�@w�w@w\)@w
=@w
=@vȴ@vE�@u�@u@vE�@v�@vV@uV@t�@tz�@t9X@tI�@u?}@v{@u��@v@v{@v$�@vV@vff@vV@v$�@v��@w;d@w�w@w�@xb@xb@xb@x  @x  @x  @x  @x �@xA�@y��@y�^@yG�@x��@xbN@x1'@xb@x �@w��@w�;@w�@xA�@xr�@xr�@xbN@w�;@w�@xb@wl�@w;d@v��@v��@v�y@v��@v��@v��@v��@v��@v��@v�+@v�+@v�+@v�+@v�+@vv�@v�+@v��@v�@v��@u��@vE�@vȴ@v5?@up�@uO�@t�@sƨ@r��@q�#@r-@q�#@r��@s��@r~�@q�@q��@q�#@r�@r=q@r=q@r�@q�#@q��@r�@q��@q�@q��@qhs@qX@p�`@p �@pb@p  @o�@o�@o|�@o+@n�y@nȴ@n�+@nE�@m�T@m`B@l�/@l9X@kƨ@k��@kS�@j�@j�!@j~�@i��@iG�@i7L@iG�@h��@h�@h1'@g��@g|�@g;d@g
=@f�R@fV@f$�@e@e�@e?}@e/@e�@eV@d�@d��@dj@d1@c�
@c��@ct�@b�@b��@bM�@a�@a�^@a&�@`Ĝ@`bN@` �@_��@_K�@_
=@^�y@^ȴ@^�@^v�@^E�@]�@]�@\��@\��@\�@\(�@[�m@[�F@[@Z�\@Z~�@Z�@Y�@Y��@Y�7@Yx�@YX@Y%@X�9@X�u@X �@W��@Wl�@W+@V��@Vȴ@V��@VV@V5?@U�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BuBoBoBhBhBhBhBbBbB\B\BPB
=B+B%BBBBBBBBBB  B��B��B��B��B��B�B�B�B�B�B�B�B�yB�yB�sB�sB�sB�yB�yB�yB�yB�sB�sB�mB�mB�mB�fB�fB�`B�ZB�NB�HB�HB�HB�HB�HB�BB�BB�HB�BB�BB�HB�HB�NB�NB�NB�NB�ZB�ZB�`B�mB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�sB�sB�mB�fB�fB�`B�ZB�TB�TB�NB�NB�TB�TB�ZB�ZB�TB�TB�NB�NB�NB�BB�/B�;B�B�B�
B��B��B��BɺBǮBƨBŢBB��B��B�}B�}B�jB�}BBÖBÖBB��B�}B�wB�jB�dB�dB�dB�dB�^B�XB�XB�RB�FB�FB�?B�-B�!B�!B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�'B�?B�FB�?B�?B�?B�?B�FB�LB�LB�RB�XB�^B�dB�jB�jB�dB�jB�wB�qB�jB�jB�jB�qB�qB�qB�wB�}B�}B��B��BBÖBÖBÖBĜBŢBƨBȴBȴBǮBɺB��B��BɺBɺBǮBŢBÖBBĜBĜBǮBɺBǮBŢBŢBƨBǮBȴBɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBǮBǮBȴBǮBǮBƨBƨBƨBƨBƨBŢBŢBĜBĜBŢBŢBŢBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBÖBÖBÖBÖBÖBÖBBBBBB��B��B��B��B��BB��B��B��B��B��B��B��B��B��B��B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B�}B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�B�B�B�B�B�B�B�B�B�B�B�B~B~B�B~B�B~B~B~B~B~B~B~B~B~BxB~B~B~B~B~B~B~B~B~B~B~B~B~BxB~B~B~B~B~BxBxB~BxBxBxBxBxBxBxBxB~BxBxBxB~B~BxBxB~B~B~B~B~B~B~B~BxBxBxBxBxBxBlBlBsBlBlBlBlBlBlBlBlBlBlBlBlBlBlBlBsBlBlBlBsBsBxBxBxBsBsBsBsBsBsBlBfBfBfBfBsBlBlBfB`B`BNBHBHBABABABAB
;B
;B	5B	5B)BBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�sB�lB�`B�ZB�ZB�TB�TB�NB�NB�NB�TB�TB�TB�TB�NB�NB�HB�HB�HB�BB�BB�<B�6B�*B�$B�$B�$B�$B�$B�B�B�$B�B�B�$B�$B�*B�*B�*B�*B�6B�6B�<B�IB�TB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�lB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�lB�lB�lB�lB�lB�lB�lB�fB�fB�fB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�ZB�ZB�TB�TB�TB�TB�OB�OB�IB�BB�BB�<B�6B�0B�0B�*B�*B�0B�0B�6B�6B�0B�0B�*B�*B�*B�B�B�B��B��B��B��B��BǰB×B��B��B��B�mB�aB�aB�[B�[B�HB�[B�mB�tB�tB�mB�aB�[B�UB�HB�BB�BB�BB�BB�<B�6B�6B�1B�%B�%B�B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�nB�nB�hB�nB�nB�hB�hB�hB�tB�{B�tB�hB�hB�hB�hB�hB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�%B�B�B�B�B�%B�+B�+B�1B�7B�=B�CB�IB�IB�CB�IB�VB�PB�IB�IB�IB�PB�PB�PB�VB�\B�\B�bB�hB�nB�uB�uB�uB�{B��B��BBB��BØBťBťBØBØB��B��B�uB�nB�{B�{B��BØB��B��B��B��B��BBØBØBØBØBğBťBťBťBťBťBťBğBğBğBğBğBğBğBğBğBğBğBğBĠBÙBBB��B��BB��B��B��B��B��B��B��B��B��B�{B�{B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�nB�nB�nB�nB�nB�hB�hB�hB�hB�hB�nB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�cB�cB�cB�cB�cB�c1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            OW : r=0.99998 , vertically averaged dS =-0.0059944                                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit). The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                     201811070904472018110709044720181107090447  IF  ARFMCODA024c                                                                20181105172617                      G�O�G�O�G�O�                IF  ARGQCOQC3.5                                                                 20181105172715  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.5                                                                 20181105172715  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20181107090447  IP  PSAL            @ffD���G�O�                