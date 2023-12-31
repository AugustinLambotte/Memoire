CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T14:43:25Z creation; 2020-11-17T12:19:03Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z          :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   L�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       V�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ^�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       `�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   h�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       j�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       r�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   z�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       |�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   �   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �t   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �x   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �|   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �H   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �H   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �H   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �HArgo profile    3.1 1.2 19500101000000  20200828144325  20220127170407  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               JA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @�%XH+�1   @�%Y�A�@Q�h�����0��aK,�1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 9.5 dbar]                                                      A!��A1��AA��AQ��A^ffAq��A���A�  A�  A���A�  A�  A�33A�  A���A�  A���A���A�  A���A�  A�  B ��B  B��B  BffB33B  B��B   B$  B'33B,  B/��B333B733B;��B@ffBC��BG��BK33BP  BS��BXffB\ffB`ffBe33BhffBlffBpffBt  Bw33Bz��B33B�  B�  B���B���B�33B�33B�33B���B���B���B�33B���B�  B���B�33B�  B���B���B���B�ffB�33B�  B�  B�  B�ffB�ffB���B�  B�  B�  B�  B�  B�33B���B�ffBǙ�B�  B���B�  B�  Bڙ�B���B�  B晚B�ffB���B�  B���B���B���CL�CffC��CL�C	33CL�C��C� C� CffC� C� C� C� C33CL�C!ffC#ffC%� C'33C)L�C+ffC-ffC/� C1� C3��C533C7ffC9� C;L�C=�C?ffCA��CCffCE33CG� CI�3CK��CM� COL�CQ33CS  CU�CWffCY��C[��C]�3C_�3Ca��Cc��Ce�3Cg�3Ci�3Ck�3Cm��Co��Cq��Cs� Cu� Cw� Cy� C{� C}ffCffC��3C��3C��fC��fC��3C��3C��3C��3C�� C�� C�� C���C���C���C�ٚC��fC��fC��fC��3C��3C�� C�� C���C���C���C�� C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C���C��fC��fC��3C�� C���C��fC��3C�� C���C��3C�� C���C��3C�� C��fC��3C���C��fC�� C���C��3C�� C�ٚC�� C���C��3C���C��3C���C³3C���Cĳ3Cř�C�� C�ٚC�� Cɳ3Cʙ�Cˌ�C̳3C�ٚC���C�� Cг3CѦfCҙ�Cә�Cԙ�CՌ�C֌�C׌�C،�Cٌ�Cڌ�Cۙ�Cܙ�Cݙ�CަfCߦfC�3C�3C�� C���C�ٚC��fC��3C�3C� C��C� C�fC�fC��fC���C���C���C���C���C���C��fC��fC��3C�� C���C��fC�Y�C��fD @ D�fD��D��D9�D�fD��D��D
@ Dy�D�3D�3D,�Dl�D�3D  D&fDs3D� D��D,�Dl�D��D��D33Ds3D �3D!��D#FfD$s3D%� D&��D(9�D)l�D*��D+�3D-@ D.s3D/�fD0��D233D3ffD4�3D5��D733D8l�D9�fD;�D<@ D=� D>��D?��DA@ DB� DC��DD��DFFfDGy�DH��DI� DK9�DL��DM��DOfDP@ DQ� DR� DT  DU@ DV� DW�fDY�DZ9�D[ffD\��D^  D_S3D`� Da� DcfDd@ Dey�Df�3Dg�3Di,�DjffDk�fDl�fDn&fDoffDp��Dq�3Ds9�Dt�fDu� Dv��Dx9�Dy� Dz�fD{�3D}&fD~y�D�3D���D�)�D���D�\�D��3D�� D�@ D���D�y�D��D��3D�\�D��3D���D�<�D�ٚD��fD�#3D��3D�c3D�  D���D�<�D���D�|�D�  D�� D�c3D��fD���D�C3D�ٚD�� D�)�D�� D�Y�D��3D�� D�I�D��fD��3D�  D�� D�` D�  D�� D�9�D�� D��fD�#3D�� D�c3D�fD���D�33D�ٚD�|�D�  D��fD�` D���D��fD�C3D��fD�vfD��D���D�` D��fD���D�C3D���D�vfD��D�ɚD�c3D���D���D�@ D���D�y�D�fD��3D�VfD��fD���D�9�D�ٚD�y�D��D�� D�` D���D��3D�33D��fD�vfD�fD���D�Y�D�  D���D�33D���D�� D�&fD¼�D�S3D��fDē3D�9�D�ٚD�y�D�  Dǹ�D�VfD��3Dɓ3D�<�D�ٚ111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A!��A1��AA��AQ��A^ffAq��A���A�  A�  A���A�  A�  A�33A�  A���A�  A���A���A�  A���A�  A�  B ��B  B��B  BffB33B  B��B   B$  B'33B,  B/��B333B733B;��B@ffBC��BG��BK33BP  BS��BXffB\ffB`ffBe33BhffBlffBpffBt  Bw33Bz��B33B�  B�  B���B���B�33B�33B�33B���B���B���B�33B���B�  B���B�33B�  B���B���B���B�ffB�33B�  B�  B�  B�ffB�ffB���B�  B�  B�  B�  B�  B�33B���B�ffBǙ�B�  B���B�  B�  Bڙ�B���B�  B晚B�ffB���B�  B���B���B���CL�CffC��CL�C	33CL�C��C� C� CffC� C� C� C� C33CL�C!ffC#ffC%� C'33C)L�C+ffC-ffC/� C1� C3��C533C7ffC9� C;L�C=�C?ffCA��CCffCE33CG� CI�3CK��CM� COL�CQ33CS  CU�CWffCY��C[��C]�3C_�3Ca��Cc��Ce�3Cg�3Ci�3Ck�3Cm��Co��Cq��Cs� Cu� Cw� Cy� C{� C}ffCffC��3C��3C��fC��fC��3C��3C��3C��3C�� C�� C�� C���C���C���C�ٚC��fC��fC��fC��3C��3C�� C�� C���C���C���C�� C���C�� C�� C�� C�� C�� C�� C�� C�� C�� C���C��fC��fC��3C�� C���C��fC��3C�� C���C��3C�� C���C��3C�� C��fC��3C���C��fC�� C���C��3C�� C�ٚC�� C���C��3C���C��3C���C³3C���Cĳ3Cř�C�� C�ٚC�� Cɳ3Cʙ�Cˌ�C̳3C�ٚC���C�� Cг3CѦfCҙ�Cә�Cԙ�CՌ�C֌�C׌�C،�Cٌ�Cڌ�Cۙ�Cܙ�Cݙ�CަfCߦfC�3C�3C�� C���C�ٚC��fC��3C�3C� C��C� C�fC�fC��fC���C���C���C���C���C���C��fC��fC��3C�� C���C��fC�Y�C��fD @ D�fD��D��D9�D�fD��D��D
@ Dy�D�3D�3D,�Dl�D�3D  D&fDs3D� D��D,�Dl�D��D��D33Ds3D �3D!��D#FfD$s3D%� D&��D(9�D)l�D*��D+�3D-@ D.s3D/�fD0��D233D3ffD4�3D5��D733D8l�D9�fD;�D<@ D=� D>��D?��DA@ DB� DC��DD��DFFfDGy�DH��DI� DK9�DL��DM��DOfDP@ DQ� DR� DT  DU@ DV� DW�fDY�DZ9�D[ffD\��D^  D_S3D`� Da� DcfDd@ Dey�Df�3Dg�3Di,�DjffDk�fDl�fDn&fDoffDp��Dq�3Ds9�Dt�fDu� Dv��Dx9�Dy� Dz�fD{�3D}&fD~y�D�3D���D�)�D���D�\�D��3D�� D�@ D���D�y�D��D��3D�\�D��3D���D�<�D�ٚD��fD�#3D��3D�c3D�  D���D�<�D���D�|�D�  D�� D�c3D��fD���D�C3D�ٚD�� D�)�D�� D�Y�D��3D�� D�I�D��fD��3D�  D�� D�` D�  D�� D�9�D�� D��fD�#3D�� D�c3D�fD���D�33D�ٚD�|�D�  D��fD�` D���D��fD�C3D��fD�vfD��D���D�` D��fD���D�C3D���D�vfD��D�ɚD�c3D���D���D�@ D���D�y�D�fD��3D�VfD��fD���D�9�D�ٚD�y�D��D�� D�` D���D��3D�33D��fD�vfD�fD���D�Y�D�  D���D�33D���D�� D�&fD¼�D�S3D��fDē3D�9�D�ٚD�y�D�  Dǹ�D�VfD��3Dɓ3D�<�D�ٚ111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���&�y�0�`�Fff�bMӿnV�z�H�|푿}p�������翌���Kǿ��-��p���I���(����������������R������ٿ����T��?}���/����ȴ���xտ�"ѿ�I����R���;��&��녿��;���`��녿ļj�ě��Ł�Ƨ��+��ȴ��Kǿ��y��?}��33��-�����Ͽ��/��$ݿ��y��+�Ƈ+��`B��Kǿ��$ݿě�����=q��t������ÿ�1�ě����Ϳ�(���b��r������I���~����7��~���$ݿ��-�����������m��;d���Ϳ��;��33������ff���m��Mӿ�j���ٿ��^���������hs��V��J���˿��m���D�����bN�u��d���P�`�;��!�7��u���������;D��=8Q�=�+>k�?Z?z�?��?j?�?dZ?	�^>�^5>��>O�>��>��>�t�>���>?|�>bN>G�>���>��^>�I�>��9>cS�>S��>A�7>49X>+>��>,1>+>6E�>(��>=p�>9X>\)>!��>C��>�(�>�V>��;>ě�>�X>���>��j>�?}>�&�>�V>�9X>�>���>�1>���>�ƨ>vȴ>o��>hr�>aG�>W
=>@�>,1>�u=��#=���<D���'�%���T�����h����1&�(�þ$�/����;d��C���������{���m�o�����پ%�   ��l�������Ƨ�9X���罝�-�����hs��+�L�ͽ#�
�49X�o��o<t�<49X=C�=��=���=��=��w=��P<T��        ;�o<�C�<���=���=>D��>["�>���>�G�>��>�1>���>�(�>��u>�n�>�=q>r�!>\(�>I�^>:^5>,1>bN=�l�=�
==��`=�E�=�Q�=Ƨ�=�=��>o>8Q�>Q�>fff>��7>���>���>�z�>�t�>��>��>��>�\)>���>�O�>�I�>�=q>��^>�1'>���>~��>w��>s�F>k�>fff>fff>e`B>`A�>]/>Z�>T��>O�;>D��>A�7>=p�>5?}>+>$�/>�w>�>�>�>��>�>V>1'>�>o=��#=��=�l�=�S�=�"�=���=�v�=�9X=�O�=]/=0 �=�P<�<���<�1<��
<u    �D���u�e`B���㼴9X���ͽ+�t���w�',1�<j�<j�D���]/�aG��aG��ixսy�#�}󶽉7L��hs��C���C�������w����ě����`���ͽ����S����F�����m�o�	7L�O߾bN�n�����R�$�/�+�1&�2-�333�8Q�;dZ�>vɾA�7�F��L�;O�;�S�ϾW
=�["ѾaG��gl��n���s�F�vȴ�{�m�}�}󶾀����\������𾈴9������I���O߾�V��\)��hs��n���zᾕ���b���������-��5?���R������S����/���/��ff��ff��l���r����羫����D�� ž��!���j��E���KǾ������#��^5���H��p����۾��۾��7��J��������1'������O߾�\)���;��녾��Ͼ��׍P��b�ٙ���"Ѿݲ-�޸R�߾w��G���MӾ��/���y��~���1��{���-��?}��E���ȴ���پ�X��dZ��vɿ���o�o��
�`B��T�+�1'��ÿ	xտ	���D��h����;��׿�����E��ȴ�KǿQ�����#��H��m�dZ�j��R�   � ��!%�!%�!G��!�7�!�7�"Mӿ#���#S��"��#���&$ݿ&ff�'+�'��'*~��-O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111�&�y�0�`�Fff�bMӿnV�z�H�|푿}p�������翌���Kǿ��-��p���I���(����������������R������ٿ����T��?}���/����ȴ���xտ�"ѿ�I����R���;��&��녿��;���`��녿ļj�ě��Ł�Ƨ��+��ȴ��Kǿ��y��?}��33��-�����Ͽ��/��$ݿ��y��+�Ƈ+��`B��Kǿ��$ݿě�����=q��t������ÿ�1�ě����Ϳ�(���b��r������I���~����7��~���$ݿ��-�����������m��;d���Ϳ��;��33������ff���m��Mӿ�j���ٿ��^���������hs��V��J���˿��m���D�����bN�u��d���P�`�;��!�7��u���������;D��=8Q�G�O�G�O�?Z?z�?��?j?�?dZ?	�^>�^5>��>O�>��>��>�t�>���>?|�>bN>G�>���>��^>�I�>��9>cS�>S��>A�7>49X>+>��>,1>+>6E�>(��>=p�>9X>\)>!��>C��>�(�>�V>��;>ě�>�X>���>��j>�?}>�&�>�V>�9X>�>���>�1>���>�ƨ>vȴ>o��>hr�>aG�>W
=>@�>,1>�u=��#=���<D���'�%���T�����h����1&�(�þ$�/����;d��C���������{���m�o�����پ%�   ��l�������Ƨ�9X���罝�-�����hs��+�L�ͽ#�
�49X�o��o<t�<49X=C�=��=���=��=��w=��P<T��        ;�o<�C�<���=���=>D��>["�>���>�G�>��>�1>���>�(�>��u>�n�>�=q>r�!>\(�>I�^>:^5>,1>bN=�l�=�
==��`=�E�=�Q�=Ƨ�=�=��>o>8Q�>Q�>fff>��7>���>���>�z�>�t�>��>��>��>�\)>���>�O�>�I�>�=q>��^>�1'>���>~��>w��>s�F>k�>fff>fff>e`B>`A�>]/>Z�>T��>O�;>D��>A�7>=p�>5?}>+>$�/>�w>�>�>�>��>�>V>1'>�>o=��#=��=�l�=�S�=�"�=���=�v�=�9X=�O�=]/=0 �=�P<�<���<�1<��
<u    �D���u�e`B���㼴9X���ͽ+�t���w�',1�<j�<j�D���]/�aG��aG��ixսy�#�}󶽉7L��hs��C���C�������w����ě����`���ͽ����S����F�����m�o�	7L�O߾bN�n�����R�$�/�+�1&�2-�333�8Q�;dZ�>vɾA�7�F��L�;O�;�S�ϾW
=�["ѾaG��gl��n���s�F�vȴ�{�m�}�}󶾀����\������𾈴9������I���O߾�V��\)��hs��n���zᾕ���b���������-��5?���R������S����/���/��ff��ff��l���r����羫����D�� ž��!���j��E���KǾ������#��^5���H��p����۾��۾��7��J��������1'������O߾�\)���;��녾��Ͼ��׍P��b�ٙ���"Ѿݲ-�޸R�߾w��G���MӾ��/���y��~���1��{���-��?}��E���ȴ���پ�X��dZ��vɿ���o�o��
�`B��T�+�1'��ÿ	xտ	���D��h����;��׿�����E��ȴ�KǿQ�����#��H��m�dZ�j��R�   � ��!%�!%�!G��!�7�!�7�"Mӿ#���#S��"��#���&$ݿ&ff�'+�'��'*~��-O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB2-B�JB�)B;dB`BB��B��B�mB49Br�B��B�BBhB'�B6FBN�By�B�B�wB�;BB]/Bz�B�B�DB��B�`B	7B/BL�Bn�B�B��B��B	%B	�B	)�B	9XB	B�B	T�B	e`B	q�B	u�B	{�B	�B	�=B	�bB	��B	��B	��B	��B	��B	�B	�B	�B	�'B	�?B	�dB	�wB	��B	ÖB	ÖB	ÖB	ƨB	��B	��B	��B	��B	��B	��B	��B	��B	�B	�
B	�NB	�5B	�B	�yB	�B	��B	��B	��B	��B
  B
B
B
B
B
+B
hB
PB
�B
#�B
2-B
?}B
I�B
L�B
O�B
ZB
]/B
dZB
iyB
o�B
~�B
�B
�bB
��B
�B
�}B
ǮB
��B
�B
��B	7BPB�B
��B9XBJ�BP�BH�BI�BG�BW
B7LB8RB.B9XB>wB>wBD�B;dB<jBA�BI�BL�BN�BO�BN�BL�BM�BM�BM�BL�BM�BN�BN�BP�BQ�BYBZB^5B`BB\)Bp�Br�Br�Br�Br�Bq�Bq�Bs�Bp�Bs�Br�Br�Bs�Br�Bt�Bq�Bp�Bq�Bq�Bp�Bp�Bp�Bn�Bl�Bk�Bt�BaHBbNBbNB^5BVB\)B^5B^5B^5BbNBO�BhsBffBe`Be`Be`Be`BffBe`BffBgmBjBk�BjBjBk�Bl�Bl�Bm�Bm�Bn�Bp�Bo�Bs�Bs�Bv�Bv�By�B{�Bz�B{�B{�B{�Bz�B{�Bw�Bx�Bx�Bz�B�By�B�B�VB�bB�oB�PB��B��B��B��B��B�VB�oB�\B�PB�PB�PB�DB�DB�DB�=B�DB�%B�%B�%B�=B�DB�JB�\B�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B6wB��B�}B?�Bd�B��B�UB��B8�BwB�=B�B�B,RB:�BS>B~AB�xB��B�B�Ba�BOB�|B��B�MB��B�B3�BQBBsB��B�lB�>B	
�B	B	.yB	=�B	GB	Y{B	i�B	v+B	zDB	�iB	��B	��B	��B	�B	�B	�(B	�9B	�kB	��B	��B	��B	��B	��B	��B	��B	�B	�B	�B	�B	�*B	�GB	�lB	�iB	�uB	�qB	�zB	�|B	�tB	څB	یB	��B	�B	�B	��B	�/B	�[B	�PB
 iB
�B
�B
�B
�B
�B
�B
�B
�B
�B
B
(]B
6�B
DB
N@B
QSB
TgB
^�B
a�B
h�B
nB
t(B
��B
��B
��B
�JB
��B
�
B
�:B
يB
ܠB wB�B�G�O�G�O�B=�BOSBU{BMIBNNBLBB[�B;�B<�B2�B=�BC	BCBI0B?�B@�BFBNNBQaBSlBTsBSlBQaBRfBRbBRdBQ_BRgBSmBSmBUyBV~B]�B^�Bb�Bd�B`�Bu9BwFBwBBwEBwEBv@Bv@BxIBu9BxKBwCBwDBxLBwEByRBv>Bu:Bv=Bv?Bu7Bu5Bu8Bs+BqBpByLBe�Bf�Bf�Bb�BZ�B`�Bb�Bb�Bb�Bf�BTnBmBj�Bi�Bi�Bi�Bi�Bj�Bi�Bj�BlBoBpBoBoBpBqBqBr%Br%Bs)Bu8Bt2BxJBxIB{\B{\B~lB�zBtB�xB�zB�xBuB�zB|dB}iB}jBwB��B~pB��B��B��B�B��B�'B�$B�B�<B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�#B�B�.B�;B�3B�5B�7B�3B�7B�7B�;B�BB�=B�@B�@B�;B�6B�-B�4B�-B�3B�5B�;B�;B�.B�0B�5B�4B�4B�1B�4B�-B�/B�/B�/B�/B�/B�5B�=B�9B�;B�<B�@B�9B�2B�:B�@B�<B�9B�<B�:B�AB�;B�6B�4B�/B�-B�-B�-B�7B�.B�+B�'B�.B�1B�1B�/B�3B�5B�9B�1B�/B�,B�6B�-B�-B�-B�4B�4B�4B�:B�7B�:B�BB�BB�CB�GB�HB�>B�GB�HB�EB�DB�FB�IB�>B�HB�EB�FB�MB�MB�LB�MB�LB�FB�FB�NB�MB�MB�IB�GB�LB�IB�LB�MB�NB�OB�MB�JB�JB�JB�NB�QB�NB�KB�KB�VB�SB�RB�KB�SB�SB�UB�UB�YB�TB�[B�TB�TB�RB�YB�SB�TB�ZB�RB�SB�SB�TB�TB�XB�SB�TB�XB�TB�YB�[B�SB�XB�YB�YB�YB�XB�[B�YB�ZB�YB�ZB�ZB�XB�XB�YB�XB�XB�XB�WB�]B�`B�`B�_B�`B�`B�aB�aB�aB�^B�^B�`B�^B�bB�aB�^B�`B�aB�]B�]B�_B�_B�_B�_B�dB�aB�^B�`B�dB�eB�dB�dB�`B�gB�fB�fB�hB�fB�eB�fB�kB�jB�fB�kB�kB�dB�kB�oB�lB�oB�lB�kB�lB�kB�kB�lB�mB�jB�jB�eB�jB�lB�jB�kB�kB�kB�qB�qB�qB�sB�iB�kB�qB�eB�kB�kB�kB�kB�jB�kB�k111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.004451 (+/- 0.01)                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040720220127170407  IF  ARFMCODA035h                                                                20200828144325                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144431  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200828144431  QCF$                G�O�G�O�G�O�0000000000004000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A!��D�ٚG�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170407  IP  PSAL            A!��D�ٚG�O�                