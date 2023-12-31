CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:24Z creation; 2020-11-17T12:19:03Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8(   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8h   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9,   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         90   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    98   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9<   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9D   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9L   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9T   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9X   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9`   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9d   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9h   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9l   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :l   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  B4   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D(   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]h   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _\   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  p�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  x�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  z�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �T   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �H   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �h   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �l   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �p   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �t   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �x   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �<   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �<   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �<   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �<Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170406  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               EA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @���q�1   @���q�@Rh�T,A2�*���=ŕ8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A���A�  A���A�33A�ffA���A�ffAٙ�A�33A�ffA�ffA�ffA�ffB��B��B��BffBffBffBffB33B#33B(  B+33B/33B2��B8  B<ffB@��BC33BG33BLffBPffBT  BX  B\  B`ffBd  BhffBlffBo33Bs33Bw��B|ffB���B�ffB�  B�ffB���B�ffB�33B���B���B�ffB�ffB�ffB���B���B�33B�33B�33B�33B�33B�33B�33B�  B���B���B�  B�33B�33B�  B���B���B�  B�33B�  B���B�  B���B���B���B���B���B�  B�33Bޙ�B�ffB晚B�33B�  B�  B�  B�  B�33C��C��C33CL�C	L�CL�CL�CffC��CL�C�CffC��C� CL�C33C!� C#��C%�3C'� C)ffC+L�C-�C/33C1�C3  C5  C7ffC9�fC;��C=��C?��CA��CC� CEffCGL�CI33CK  CMffCO��CQ��CS� CUffCW33CY�C[  C]ffC_�3CaffCc33CeffCgL�Ci33Ck33CmL�CoffCq33CsL�Cu� CwL�Cy33C{  C}ffC�fC��3C��3C��3C��3C��3C��3C��3C�� C�� C���C���C���C���C���C��fC��3C��3C��fC��3C���C���C�� C�� C�� C�� C�� C�� C�� C���C���C�ٚC��fC��3C�� C�� C���C��fC��3C�� C���C��fC��3C��3C���C��fC��3C�� C���C��fC��3C���C�ٚC��3C�� C���C��fC��3C�� C���C��fC��3C���C�ٚC��fC��3C�� C¦fCó3C�� CŦfCƳ3C���CȦfC�� C���C˦fC�� C���CΦfC�� C�ٚC�� CҦfC�� C��fC�� Cֳ3Cי�Cس3C���Cڳ3Cۙ�C܀ Cݳ3C�ٚC���C�� C�3C�fC㙚C䙚C��C��C��C� C� C� C� C� C� C� C� C�� C� C� C��C���C�� C���C���C���C���C�s3C�&fD 9�Dl�D��D�D@ Dy�D�3D�3D
33Ds3D��D��D@ Dl�D� D3DFfDy�D�3D��D,�DffD�fD� D  Dy�D ٚD"�D#9�D$` D%��D&��D(,�D)y�D*� D+�3D-,�D.� D/�3D1�D2L�D3��D4��D63D79�D8` D9�fD:�fD<,�D=s3D>��D@fDA@ DBs3DC�3DD��DFFfDG�fDH� DJfDK,�DLy�DM�fDN�3DP,�DQl�DR�3DS��DU,�DVl�DW��DX��DZ,�D[s3D\��D^fD_9�D`s3Da�3Db��Dd&fDeffDf��Dg�3Di33Djy�Dk� Dm�Dn9�Dol�Dp� Dq�3Ds33Dtl�Du�3Dv�3Dx@ Dy��Dz� D{�3D},�D~ffD� D�p D� D�� D�S3D�3D��fD�@ D���D�y�D�fD��3D�S3D��fD���D�C3D��3D��3D�#3D��fD�Y�D���D�� D�@ D��3D�y�D��D�� D�c3D���D��3D�9�D�� D�|�D�3D���D�i�D�fD��3D�C3D��3D��fD�&fD��fD�i�D��D���D�0 D��3D�y�D�  D���D�VfD���D��fD�C3D��3D��3D�&fD���D�S3D���D��fD�@ D���D�|�D��D��3D�Y�D�  D��3D�@ D�ٚD�vfD�3D��fD�Y�D���D��fD�@ D�ٚD�s3D��D��3D�c3D�  D�� D�<�D��3D��3D�#3D��3D�c3D�3D��3D�FfD�� D�vfD��D��fD�` D���D���D�33D�� D�� D�  D�� D�Y�D�  Dģ3D�<�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A���A�  A���A�33A�ffA���A�ffAٙ�A�33A�ffA�ffA�ffA�ffB��B��B��BffBffBffBffB33B#33B(  B+33B/33B2��B8  B<ffB@��BC33BG33BLffBPffBT  BX  B\  B`ffBd  BhffBlffBo33Bs33Bw��B|ffB���B�ffB�  B�ffB���B�ffB�33B���B���B�ffB�ffB�ffB���B���B�33B�33B�33B�33B�33B�33B�33B�  B���B���B�  B�33B�33B�  B���B���B�  B�33B�  B���B�  B���B���B���B���B���B�  B�33Bޙ�B�ffB晚B�33B�  B�  B�  B�  B�33C��C��C33CL�C	L�CL�CL�CffC��CL�C�CffC��C� CL�C33C!� C#��C%�3C'� C)ffC+L�C-�C/33C1�C3  C5  C7ffC9�fC;��C=��C?��CA��CC� CEffCGL�CI33CK  CMffCO��CQ��CS� CUffCW33CY�C[  C]ffC_�3CaffCc33CeffCgL�Ci33Ck33CmL�CoffCq33CsL�Cu� CwL�Cy33C{  C}ffC�fC��3C��3C��3C��3C��3C��3C��3C�� C�� C���C���C���C���C���C��fC��3C��3C��fC��3C���C���C�� C�� C�� C�� C�� C�� C�� C���C���C�ٚC��fC��3C�� C�� C���C��fC��3C�� C���C��fC��3C��3C���C��fC��3C�� C���C��fC��3C���C�ٚC��3C�� C���C��fC��3C�� C���C��fC��3C���C�ٚC��fC��3C�� C¦fCó3C�� CŦfCƳ3C���CȦfC�� C���C˦fC�� C���CΦfC�� C�ٚC�� CҦfC�� C��fC�� Cֳ3Cי�Cس3C���Cڳ3Cۙ�C܀ Cݳ3C�ٚC���C�� C�3C�fC㙚C䙚C��C��C��C� C� C� C� C� C� C� C� C�� C� C� C��C���C�� C���C���C���C���C�s3C�&fD 9�Dl�D��D�D@ Dy�D�3D�3D
33Ds3D��D��D@ Dl�D� D3DFfDy�D�3D��D,�DffD�fD� D  Dy�D ٚD"�D#9�D$` D%��D&��D(,�D)y�D*� D+�3D-,�D.� D/�3D1�D2L�D3��D4��D63D79�D8` D9�fD:�fD<,�D=s3D>��D@fDA@ DBs3DC�3DD��DFFfDG�fDH� DJfDK,�DLy�DM�fDN�3DP,�DQl�DR�3DS��DU,�DVl�DW��DX��DZ,�D[s3D\��D^fD_9�D`s3Da�3Db��Dd&fDeffDf��Dg�3Di33Djy�Dk� Dm�Dn9�Dol�Dp� Dq�3Ds33Dtl�Du�3Dv�3Dx@ Dy��Dz� D{�3D},�D~ffD� D�p D� D�� D�S3D�3D��fD�@ D���D�y�D�fD��3D�S3D��fD���D�C3D��3D��3D�#3D��fD�Y�D���D�� D�@ D��3D�y�D��D�� D�c3D���D��3D�9�D�� D�|�D�3D���D�i�D�fD��3D�C3D��3D��fD�&fD��fD�i�D��D���D�0 D��3D�y�D�  D���D�VfD���D��fD�C3D��3D��3D�&fD���D�S3D���D��fD�@ D���D�|�D��D��3D�Y�D�  D��3D�@ D�ٚD�vfD�3D��fD�Y�D���D��fD�@ D�ٚD�s3D��D��3D�c3D�  D�� D�<�D��3D��3D�#3D��3D�c3D�3D��3D�FfD�� D�vfD��D��fD�` D���D���D�33D�� D�� D�  D�� D�Y�D�  Dģ3D�<�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����
=��
=��������+����+��E���ff��ff��E���E���$ݿ�E���E���$ݿ��T����˿�˿�˿�����`B��`B��`B��?}��`B��`B��˿��T����$ݿ�ff��+���T��?}��z��9X��z���/�䛦��9X��t���Mӿ��������n����Ͽ��
��F�㕁��F��F���Ͽ�9X���Ͽ��Ͽ�9X��Z�䛦��j��;d�ۅ��V��A���!��C��ݑh��t���$ݿ���j����E�����$ݿ�������33��=q��?}�����������+����Ѓ��vɿ�;d���#��`B��������z῟���Q쿳S��������ff��$ݿ�$ݿ��P���9��r���7L��=q��-�������/�����33�����������R��V������h��xտ������7��/��(���^5��?}���忢Mӿ�O߿�Kǿ����R��$ݿ�  �q녿cS��\��U��Q녿C�
�8Q�*=q��Ͼ�E���J��7L����<j�C��ȴ9����=@�=��>}�>���>��>�>�ȴ>ݲ->ۥ�? Ĝ?&�?�? Ĝ>�E�>���>��
>�R=�w=�o=Ƨ�>$�>D��>j~�>��7>�(�>��j?��?	��?��?bN?��?K�? A�?p�?^5??}?�F?�w?"M�?!��?"J?"M�?"�\?%�?&�y?%�T?!%? �?!��?%`B?&$�?#��?)�^?+ƨ?(�9?�-?��?�!?�?	x�?V?��?{?	��?
��?r�?��?��? Ĝ?o?S�?��>��m>�>�->��>���>��/>޸R>ȴ9>��>���>�z�>�bN>�\)>�n�>�1'>�/>�Ĝ>��>�>��j>�Q�>���>�F>�1>���>�A�>�5?>�(�>�"�>�"�>ڟ�>ܬ>�/>ݲ->ݲ->�/>��>�hs>��>�dZ>�%>�+>���>�\)>�bN>�t�>�1'>��j>�33>�33>�->��D>���>�l�>�ff>�ff>��y>��>��/>�G�>�;d>�"�>�
=>��>�\)>�\)>���>�+>~��>cS�>_;d>\(�>Xb>Kƨ>;dZ>9X>9X>9X>:^5>;dZ><j>$�/> Ĝ>�u>
=q>o>bN>n�>bN>
=q>o=�=�x�=�/=���=��=� �=�O�=Y�=e`B=}�=e`B=@�<�/;�`B;o        ��o;�`B<T��<�o<e`B<t�;ě�;�o�D����o�o��o���
���
��j���<j��o�y�#�P�`�H�9�P�`�D���Y��e`B�u��hs���-��-��^5��vɽƧ����/��`B��xս�����1'�
=q�V�t�����,1�5?}�6E��;dZ�?|�A�7�E�˾P�`�V�["Ѿ`A��hr��q���t�j�vȴ�z�H��%���\����������+��C���O߾�\)���;��n�������P���������"Ѿ�����R��A���G����徢�徢MӾ�MӾ�Z��r���V��V��V������׾�����?}��E���KǾ�j��|�\����Ǯ�ȴ9�ɺ^������C���I�������;���`������
=��
=�����"Ѿ޸R��A���G���Z���1��������33��9X��?}���پ��پ�����^5�   �   �G��Mӿo�������/��T����ÿ	�^��ƨ��Ϳ��V����� ſ���녿녿-�-�n���!��!�n���!�t���j�9X��j���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ��
=��
=��������+����+��E���ff��ff��E���E���$ݿ�E���E���$ݿ��T����˿�˿�˿�����`B��`B��`B��?}��`B��`B��˿��T����$ݿ�ff��+���T��?}��z��9X��z���/�䛦��9X��t���Mӿ��������n����Ͽ��
��F�㕁��F��F���Ͽ�9X���Ͽ��Ͽ�9X��Z�䛦��j��;d�ۅ��V��A���!��C��ݑh��t���$ݿ���j����E�����$ݿ�������33��=q��?}�����������+����Ѓ��vɿ�;d���#��`B��������z῟���Q쿳S��������ff��$ݿ�$ݿ��P���9��r���7L��=q��-�������/�����33�����������R��V������h��xտ������7��/��(���^5��?}���忢Mӿ�O߿�Kǿ����R��$ݿ�  �q녿cS��\��U��Q녿C�
�8Q�*=q��Ͼ�E���J��7L����<j�C��ȴ9����=@�=��>}�>���>��>�>�ȴ>ݲ->ۥ�? Ĝ?&�?�? Ĝ>�E�>���>��
>�R=�w=�o=Ƨ�>$�>D��>j~�>��7>�(�>��j?��?	��?��?bN?��?K�? A�?p�?^5??}?�F?�w?"M�?!��?"J?"M�?"�\?%�?&�y?%�T?!%? �?!��?%`B?&$�?#��?)�^?+ƨ?(�9?�-?��?�!?�?	x�?V?��?{?	��?
��?r�?��?��? Ĝ?o?S�?��>��m>�>�->��>���>��/>޸R>ȴ9>��>���>�z�>�bN>�\)>�n�>�1'>�/>�Ĝ>��>�>��j>�Q�>���>�F>�1>���>�A�>�5?>�(�>�"�>�"�>ڟ�>ܬ>�/>ݲ->ݲ->�/>��>�hs>��>�dZ>�%>�+>���>�\)>�bN>�t�>�1'>��j>�33>�33>�->��D>���>�l�>�ff>�ff>��y>��>��/>�G�>�;d>�"�>�
=>��>�\)>�\)>���>�+>~��>cS�>_;d>\(�>Xb>Kƨ>;dZ>9X>9X>9X>:^5>;dZ><j>$�/> Ĝ>�u>
=q>o>bN>n�>bN>
=q>o=�=�x�=�/=���=��=� �=�O�=Y�=e`B=}�=e`B=@�<�/;�`B;o        ��o;�`B<T��<�o<e`B<t�;ě�;�o�D����o�o��o���
���
��j���<j��o�y�#�P�`�H�9�P�`�D���Y��e`B�u��hs���-��-��^5��vɽƧ����/��`B��xս�����1'�
=q�V�t�����,1�5?}�6E��;dZ�?|�A�7�E�˾P�`�V�["Ѿ`A��hr��q���t�j�vȴ�z�H��%���\����������+��C���O߾�\)���;��n�������P���������"Ѿ�����R��A���G����徢�徢MӾ�MӾ�Z��r���V��V��V������׾�����?}��E���KǾ�j��|�\����Ǯ�ȴ9�ɺ^������C���I�������;���`������
=��
=�����"Ѿ޸R��A���G���Z���1��������33��9X��?}���پ��پ�����^5�   �   �G��Mӿo�������/��T����ÿ	�^��ƨ��Ϳ��V����� ſ���녿녿-�-�n���!��!�n���!�t���j�9X��j���11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB1B
=BDBPBJBJBPBVBPBVBVBVB\B\BVBbBhBhB�B�B�B�B�B�B�B�B�B �B#�B&�B.B2-B2-B5?BT�Bv�B�bB��B��B�`B�B��B��B�B(�BM�B>wB@�BL�Bm�Bo�Br�Bt�Bv�Bw�Bx�Bz�B�VB�uB��B��B��B��B��BƨB��B��B�B�/B��B�B��B��B��B��B	  B	B	%B	1B	
=B	PB	�B	�B	$�B	.B	2-B	/B	=qB	C�B	I�B	L�B	N�B	XB	bNB	jB	v�B	�B	�B	�1B	�VB	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	�?B	�RB	��B	ŢB	ǮB	��B	�#B	�;B	�HB	�B	��B
1B
�B
"�B
,B
0!B
33B
?}B
D�B
G�B
I�B
[#B
bNB
n�B
�B
�JB
��B
��B
��B
�B
�LB
��B
��B
�B
�sB
�B
��B1B�BPB�B �B$�B33B8RBB�BH�BH�BB�BE�BT�B\)B]/BgmBffBbNB`BB\)B_;BD�BH�BP�BR�BVBaHBe`BjBx�B� B�JB�%B�7B�=B�=B�JB�bB�\B�\B�PB�\B�hB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�bB�\B�JB�PB�VB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BBBBB!B)BB,B >B >B#QB%_B(sB+�B2�B6�B6�B9�BY�B{fB��B�=B��B�B� B�[B�B(B-�BRtBCBE&BQqBr4BtFBwVByaB{qB|vB}{B�B��B�B�PB�eB��B�.B�1B�RB�tB֒B��B��B�gB�<B�fB�iB�tB��B	�B		�B	
�B	�B	�B	�B	5B	2B	)�B	2�B	6�B	3�B	BB	HFB	NiB	Q{B	S�B	\�B	f�B	o1B	{zB	��B	��B	��B	�
B	�KB	�hB	�vB	�|B	�{B	��B	��B	��B	��B	��B	��B	�B	�5B	�SB	�cB	ԔB	��B	��B	��B	�BB	��B
�B
VB
'�B
0�B
4�B
7�B
D5B
ITB
LfB
NrB
_�B
g	B
sSB
��B
�B
�:B
��B
��B
��B
�B
�EB
�B
��B
�1B
�bB�B�BFBB fB%�B)�B7�B=BGSBMxBMwBGQBJgBY�B`�Ba�Bl0Bk,BgBeB`�Bc�BI_BMvBU�BW�BZ�BfBj"BoBB}�B��B�B��B��B�B�B�B�&B�"B�"B�B�!B�/B�=B�GB�GB�OB�LB�UB�XB�`B�VB�_B�jB�kB�nB�gB��B��B�xB�kB�gB�XB�mB�TB�YB�YB�aB�kB�`B�gB�XB�YB�XB�`B�dB�tB�aB�gB�ZB�ZB�[B�FB�LB�LB�&B� B�B�B�B�4B�MB�vB�xB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�yB�B�xB�kB�}B��B��B��B��B��B��B��B��B��B��B�uB��B�~B��B�xB��B��B��B��B��B��B��B��B��B��B��B��B�|B��B�{B�B��B��B�~B��B��B��B��B��B��B��B��B�|B�B�}B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B�|B�uB�{B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�|B�B�~B�zB�wB��B�B�|B�}B�~B�~B�~B�~B��B�~B�}B�{B�B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0046275 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040620220127170406  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144428  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144428  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A���D�<�G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170406  IP  PSAL            A���D�<�G�O�                