CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
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
resolution        =���   axis      Z          :p   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   B�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          D�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   L�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       V�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ^�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       `�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   i    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       k   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       s    PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   {8   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       }@   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   �X   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       �`   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �$   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �4   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �8   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �H   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �L   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �P   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �T   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �x   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170407  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               GA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @��I��J1   @��I��J@R,�U����/߇�j�p8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A1��A@  AL��A\��Ap  A���A�33A���A�ffA�ffA���A�  A�  A���A�33A�33A�ffA�ffA�33A���A�33B   B��B��BffBffB��B  B��B ��B#��B'��B,ffB/��B3��B8ffB;��B>��BC��BHffBL  BP  BS��BW��B[��B_33Bc��Bg33Bk��Bo��Bs��Bw��B|  B~��B���B���B�  B�ffB�  B���B�  B�33B���B�ffB���B�33B�  B���B�  B�ffB�33B���B�ffB�33B���B�ffB�33B���B���B�33B�33B���B���B�33B�  B���B�ffB�33B���Bș�B�33B���B�ffB�33B���B���B㙚B癚B���B���B���B�  B�  B�33C33C33CL�CL�C	ffC� C�3C��CffC  C�C33C33CffC� C��C!�3C#��C%�fC'� C)  C+�C-33C/� C1��C3��C5��C7��C9��C;L�C=L�C?ffCA� CC��CE33CGL�CIffCK33CMffCO��CQL�CS�CUL�CW� CYL�C[� C]��C_��Ca� CcL�Ce33Cg�CiffCk��Cm�3Co��Cq� CsffCuL�Cw33Cy�C{  C}ffC�fC��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C���C�� C�� C���C���C���C�ٚC�ٚC��fC��3C�� C�� C���C���C��fC��3C��fC��fC��fC��fC���C���C��3C��fC��fC��fC��fC��fC��fC��fC��3C��3C��3C�� C���C��fC��3C�� C���C��fC�� C���C��fC�� C���C��fC�� C���C��3C�� C��fC��3C���C��3C�� C�ٚC�� CæfC�� C�ٚCƳ3CǙ�Cȳ3C���Cʳ3Cˌ�C̳3C���Cγ3Cό�Cг3C�ٚC�� CӦfCԌ�Cճ3C�ٚC�� CئfCٙ�C�� C��fC���Cݳ3CަfCߙ�C�� C��fC�ٚC���C�3C�fC��C�3C��fC�ٚC�� C�3C�fC��C� C�3C��fC�ٚC���C�3C��fC���C���C���C��3C��fC�� C��fD @ D�3D�fD��D,�D` D��D	3D
S3D�3D�3D3D@ DffD��D�3D,�Ds3D��DfD@ Dy�D��D�3D9�D� D ��D!��D#,�D$y�D%�3D'�D(@ D)� D*�fD,  D-FfD.�3D/� D0�fD233D3�fD4� D5�3D7,�D8y�D9�fD;fD<FfD=�fD>�fD@�DA@ DBl�DC��DD�3DF,�DG�fDH��DJ  DK33DLy�DM��DO  DP,�DQffDR��DT�DU@ DVy�DW��DY  DZL�D[� D\�3D]�fD_33D`�fDa��Db�fDd9�De��Df��DhfDi@ Dj� Dk� Dm  Dn@ Do�fDp��Dq�3Ds  DtffDu��Dv��Dx@ Dy�3Dz��D|�D}9�D~` D��D�|�D�#3D��3D�VfD���D���D�9�D�ٚD�y�D��D��fD�\�D�fD��3D�@ D�� D�� D�  D���D�` D�fD��3D�@ D�� D��3D�&fD�ɚD�` D��fD���D�9�D��3D�|�D��D���D�Y�D���D���D�@ D��fD�|�D� D���D�c3D�  D���D�9�D�ٚD�y�D��D��3D�\�D�  D��fD�C3D�� D�� D�#3D��3D�ffD�  D��fD�@ D�ٚD�vfD��D���D�Y�D�3D�� D�<�D���D�y�D��D��fD�Y�D�  D��fD�<�D��3D�|�D�&fD��3D�` D���D���D�<�D���D�� D�&fD�ɚD�i�D�  D��3D�9�D�� D���D�#3D���D�c3D�  D��fD�C3D���D�y�D�#3D�� D�Y�D��fDē3D�C3D�� D�|�D��DǶfD�ffD�3Dɣ3D�C3D���D�y�D��D̀ 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A1��A@  AL��A\��Ap  A���A�33A���A�ffA�ffA���A�  A�  A���A�33A�33A�ffA�ffA�33A���A�33B   B��B��BffBffB��B  B��B ��B#��B'��B,ffB/��B3��B8ffB;��B>��BC��BHffBL  BP  BS��BW��B[��B_33Bc��Bg33Bk��Bo��Bs��Bw��B|  B~��B���B���B�  B�ffB�  B���B�  B�33B���B�ffB���B�33B�  B���B�  B�ffB�33B���B�ffB�33B���B�ffB�33B���B���B�33B�33B���B���B�33B�  B���B�ffB�33B���Bș�B�33B���B�ffB�33B���B���B㙚B癚B���B���B���B�  B�  B�33C33C33CL�CL�C	ffC� C�3C��CffC  C�C33C33CffC� C��C!�3C#��C%�fC'� C)  C+�C-33C/� C1��C3��C5��C7��C9��C;L�C=L�C?ffCA� CC��CE33CGL�CIffCK33CMffCO��CQL�CS�CUL�CW� CYL�C[� C]��C_��Ca� CcL�Ce33Cg�CiffCk��Cm�3Co��Cq� CsffCuL�Cw33Cy�C{  C}ffC�fC��fC�ٚC�ٚC���C���C�� C�� C�� C�� C�� C���C�� C�� C���C���C���C�ٚC�ٚC��fC��3C�� C�� C���C���C��fC��3C��fC��fC��fC��fC���C���C��3C��fC��fC��fC��fC��fC��fC��fC��3C��3C��3C�� C���C��fC��3C�� C���C��fC�� C���C��fC�� C���C��fC�� C���C��3C�� C��fC��3C���C��3C�� C�ٚC�� CæfC�� C�ٚCƳ3CǙ�Cȳ3C���Cʳ3Cˌ�C̳3C���Cγ3Cό�Cг3C�ٚC�� CӦfCԌ�Cճ3C�ٚC�� CئfCٙ�C�� C��fC���Cݳ3CަfCߙ�C�� C��fC�ٚC���C�3C�fC��C�3C��fC�ٚC�� C�3C�fC��C� C�3C��fC�ٚC���C�3C��fC���C���C���C��3C��fC�� C��fD @ D�3D�fD��D,�D` D��D	3D
S3D�3D�3D3D@ DffD��D�3D,�Ds3D��DfD@ Dy�D��D�3D9�D� D ��D!��D#,�D$y�D%�3D'�D(@ D)� D*�fD,  D-FfD.�3D/� D0�fD233D3�fD4� D5�3D7,�D8y�D9�fD;fD<FfD=�fD>�fD@�DA@ DBl�DC��DD�3DF,�DG�fDH��DJ  DK33DLy�DM��DO  DP,�DQffDR��DT�DU@ DVy�DW��DY  DZL�D[� D\�3D]�fD_33D`�fDa��Db�fDd9�De��Df��DhfDi@ Dj� Dk� Dm  Dn@ Do�fDp��Dq�3Ds  DtffDu��Dv��Dx@ Dy�3Dz��D|�D}9�D~` D��D�|�D�#3D��3D�VfD���D���D�9�D�ٚD�y�D��D��fD�\�D�fD��3D�@ D�� D�� D�  D���D�` D�fD��3D�@ D�� D��3D�&fD�ɚD�` D��fD���D�9�D��3D�|�D��D���D�Y�D���D���D�@ D��fD�|�D� D���D�c3D�  D���D�9�D�ٚD�y�D��D��3D�\�D�  D��fD�C3D�� D�� D�#3D��3D�ffD�  D��fD�@ D�ٚD�vfD��D���D�Y�D�3D�� D�<�D���D�y�D��D��fD�Y�D�  D��fD�<�D��3D�|�D�&fD��3D�` D���D���D�<�D���D�� D�&fD�ɚD�i�D�  D��3D�9�D�� D���D�#3D���D�c3D�  D��fD�C3D���D�y�D�#3D�� D�Y�D��fDē3D�C3D�� D�|�D��DǶfD�ffD�3Dɣ3D�C3D���D�y�D��D̀ 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���̬���Ϳ��� ſ�G���������ϿԼj�ա˿�E���ff��E���`B�Լj���
�щ7��|������Ϳ�I��̬��V��/��V���Ϳ�~���=q�ɺ^���ٿ��ٿ�
=���ٿǍP�ɺ^��"ѿ˥��dZ�������Ϳ͑h�θR��J��-��hs��V���ٿŁ����������/��Ĝ������~���푿z��s33�vE��t9X�t9X�^vɿF�y�>5?�5��$���
~����
=q���1'�hs��ٿ7
=���T�����S�F�š˾�xվ����G��S�Ͼ�!��l��t���h��O�<�/>J��>��\>�33>�ff>ٙ�>���-V��v�=��>�  ?"�\>�z�?"�\?e`B?]p�?��y?��?�(�?��D?��H?gl�?&ff>�`B?NV?�9X?��R?�p�?�|�?�hs?�Ĝ?��P?���?��?�&�?�E�?���?��?�^5?�Z?�ȴ?���?�I�?{"�?hr�?lI�?|�?��?�Ĝ?��?� �?�?}?�hs?��-?��?{dZ?u?}?n�?c�
?^�R?_�w?^��?`�?d�?^v�?]�-?]�-?f$�?r�!?tz�?j=q?g+?l��?tz�?r�?p��?m�h?n{?j=q?f$�?\�?R-?Tz�?VE�?M��?K?I�^?PbN?S��?S��?S��?S�F?St�?R�?R�!?R-?O�;?O\)?O�?N�?N��?N{?M�h?J~�?F$�?C�
?C��?E�T?DZ?A%?8Q�?:�?G+?D�/?E�?F��?C��?D�/?E�T?Fff?G�?Fff?C��?@�??;d?:^5?7K�?5?}?6�+?6ȴ?6ȴ?6?333?1hs?2n�?0bN?/�?,��?-O�?/�;?.�?,I�?+C�?)7L?*=q?+ƨ?)7L?(1'?'+?%��?$Z?#�
?#o?"M�?!��?�R?/?"�?^5?^5?^5?�?�u?E�?t�?��?�?bN?V?	7L?	7L?��?�9?�y?��?�T?$�?$�?��?�?o?��?M�?�7>���>�?}>��>��>>�V>�~�>��y>��
>���>�/>�
=>���>�bN>���>ɺ^>Ƨ�>�%>��H>�E�>���>�~�>�M�>�;d>���>�ƨ>��>�J>x��>k�>dZ>Q�>H�9>=p�>6E�>,1>�->�+>	7L>1'>J=�x�=�
==���=���=���=�o=q��=T��=D��=\)<�9X<T��;ě�;��
    ��o�t��49X�T����C����㼣�
��1������`B��h��h��w�,1�49X�8Q�<j�L�ͽy�#��+���P��������ȴ9���`��
=��
=��S���l���h���#���
=q�O߾hs��+��R�!���'.{�333�<j�?|�D���L�;Q녾["ѾaG��bMӾdZ�e`B�l�D�s�F�z�H�~�۾�%��������������˾��˾����ƨ��O߾�bN������
=�����(���Ĝ�������/����羬1��V��V�������F����X���H��vɾ�J�����7L��ƨ���`��n������Ͼ�
=�ٙ��ڟ��ܬ��;d��Ĝ������`B�������{��׾�-��F��ȴ��X���H��j��vɿ   � A���7�J��\�o��
����`B��T�ff�l��1'�	��ƨ�V��h�V���� ſhs�������녿n���z��j��
=�b����H�dZ���j�/�5?�|�   �!���"Mӿ#���$�/�%`B�&ff�&�y�'(�ÿ)��*~��*���+��.{�/��/\)�/�;�1&�2-�2�!�3t��4�j�4���6E��8b�8���9X�:^5�;��<(��<푿=/�=p��>5?�?|�@Ĝ�B�\�C�
�D���E`B�E�˿E��11111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  �̬���Ϳ��� ſ�G���������ϿԼj�ա˿�E���ff��E���`B�Լj���
�щ7��|������Ϳ�I��̬��V��/��V���Ϳ�~���=q�ɺ^���ٿ��ٿ�
=���ٿǍP�ɺ^��"ѿ˥��dZ�������Ϳ͑h�θR��J��-��hs��V���ٿŁ����������/��Ĝ������~���푿z��s33�vE��t9X�t9X�^vɿF�y�>5?�5��$���
~����
=q���1'�hs��ٿ7
=���T�����S�F�š˾�xվ����G��S��G�O�G�O��t���h��O�<�/>J��>��\>�33>�ff>ٙ�>���-V��v�=��>�  ?"�\>�z�?"�\?e`BG�O�G�O�?��?�(�?��D?��H?gl�?&ff>�`B?NV?�9X?��R?�p�?�|�?�hs?�Ĝ?��P?���?��?�&�?�E�?���?��?�^5?�Z?�ȴ?���?�I�?{"�?hr�?lI�?|�?��?�Ĝ?��?� �?�?}?�hs?��-?��?{dZ?u?}?n�?c�
?^�R?_�w?^��?`�?d�?^v�?]�-?]�-?f$�?r�!?tz�?j=q?g+?l��?tz�?r�?p��?m�h?n{?j=q?f$�?\�?R-?Tz�?VE�?M��?K?I�^?PbN?S��?S��?S��?S�F?St�?R�?R�!?R-?O�;?O\)?O�?N�?N��?N{?M�h?J~�?F$�?C�
?C��?E�T?DZ?A%?8Q�?:�?G+?D�/?E�?F��?C��?D�/?E�T?Fff?G�?Fff?C��?@�??;d?:^5?7K�?5?}?6�+?6ȴ?6ȴ?6?333?1hs?2n�?0bN?/�?,��?-O�?/�;?.�?,I�?+C�?)7L?*=q?+ƨ?)7L?(1'?'+?%��?$Z?#�
?#o?"M�?!��?�R?/?"�?^5?^5?^5?�?�u?E�?t�?��?�?bN?V?	7L?	7L?��?�9?�y?��?�T?$�?$�?��?�?o?��?M�?�7>���>�?}>��>��>>�V>�~�>��y>��
>���>�/>�
=>���>�bN>���>ɺ^>Ƨ�>�%>��H>�E�>���>�~�>�M�>�;d>���>�ƨ>��>�J>x��>k�>dZ>Q�>H�9>=p�>6E�>,1>�->�+>	7L>1'>J=�x�=�
==���=���=���=�o=q��=T��=D��=\)<�9X<T��;ě�;��
    ��o�t��49X�T����C����㼣�
��1������`B��h��h��w�,1�49X�8Q�<j�L�ͽy�#��+���P��������ȴ9���`��
=��
=��S���l���h���#���
=q�O߾hs��+��R�!���'.{�333�<j�?|�D���L�;Q녾["ѾaG��bMӾdZ�e`B�l�D�s�F�z�H�~�۾�%��������������˾��˾����ƨ��O߾�bN������
=�����(���Ĝ�������/����羬1��V��V�������F����X���H��vɾ�J�����7L��ƨ���`��n������Ͼ�
=�ٙ��ڟ��ܬ��;d��Ĝ������`B�������{��׾�-��F��ȴ��X���H��j��vɿ   � A���7�J��\�o��
����`B��T�ff�l��1'�	��ƨ�V��h�V���� ſhs�������녿n���z��j��
=�b����H�dZ���j�/�5?�|�   �!���"Mӿ#���$�/�%`B�&ff�&�y�'(�ÿ)��*~��*���+��.{�/��/\)�/�;�1&�2-�2�!�3t��4�j�4���6E��8b�8���9X�:^5�;��<(��<푿=/�=p��>5?�?|�@Ĝ�B�\�C�
�D���E`B�E�˿E��11111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�PB�bB��B��B�TB��BB
=BbB�B0!B;dBG�Bl�B�7B�bB��B��B�B�^B�}BȴB��B��B�#B�B��B��B��B	B	1B	PB	{B	�B	�B	$�B	&�B	'�B	(�B	0!B	33B	49B	<jB	G�B	N�B	M�B	XB	\)B	_;B	W
B	l�B	w�B	|�B	q�B	�bB	�PB	��B	�!B	�-B	�FB	�'B	��B	ǮB	��B	��B	�#B	�5B	�NB	�ZB	�sB	�B	�sB	�BB	�TB	�B	�B

=B
�B
�B
{B
#�B
1'B	�HB
:^B
E�B
W
B
ffB
n�B
}�B
�+B
�DB
�\B
��B
�oB
��B
�B
��B
�NB
��B1B%B%B
�)B1'B33B2-B8RB{B	7B�BC�BN�BZBK�BZBZBI�BJ�BL�BN�BL�BT�BYBM�BZB]/BYB_;BcTB\)BZB\)Bk�Bm�Bs�B{�Bu�By�Bt�Bs�Br�Bm�Bl�Bm�BiyBk�Bo�Bs�Bt�Bx�Bz�Bz�B{�Bk�B�%B�B�B�%B�7B�DB�DB�JB�DB�PB�7B�DB�1B�+B�DB�1B�B�+B�+B�1B�PB�VB�PB�VB�\B�\B�\B�\B�hB�hB�hB�hB�oB�oB�oB�hB�hB�oB�bB�uB�oB�oB�oB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B��B��B��B�bB��B|B	�B�B�B.B4�B?�BL?BqB��B��B�:B�gB��B��B�B�LB�uB�zB߸B�0B�yB	�B	 �B	�B	�B	�B	B	B	#PB	)vB	+�B	,�B	-�B	4�B	7�B	8�B	AB	LGB	StB	RnB	\�B	`�B	c�B	[�B	q)B	|lB	��B	vHB	�B	��B	��B	��B	��B	��B	��B	�*B	�MB	�iB	�sB	��B	��B	��B	��B	�B	��B	�B	��B	��B	�YB	�YB
�B
1B
0B
B
(zG�O�G�O�B
?B
JIB
[�B
kB
s?B
��B
��B
��B
�B
�=B
�B
�6B
��B
�4B
��B
�zB�B
�G�O�G�O�B5�B7�B6�B=B*B�B$kBHHBS�B^�BPxB^�B^�BNkBOrBQ{BS�BQ{BY�B]�BR�B^�Ba�B]�Bc�BhB`�B^�B`�Bp7BrEBxjB��BzvB~�BynBxiBw`BrBBq;BrBBn+Bp4BtQBxgByoB}�B�B�B��Bp3B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�#B� B�B�B�#B�B�&B�!B�!B� B� B�AB�;B�<B�9B�RB�FB�XB�TB�ZB�ZB�_B�SB�ZB�AB�OB�EB�RB�VB�TB�]B�dB�WB�_B�eB�WB�]B�^B�_B�mB�qB�eB�kB�qB�mB�qB�qB�qB�qB�rB�wB�wB�xB�pB�rB�pB�sB�sB�pB�kB�pB�sB�pB�kB�pB�kB�rB�jB�xB�pB�rB�qB�xB�kB�sB�sB�rB�yB�yB�qB�jB�yB�dB�sB�iB�^B�dB�dB�dB�cB�dB�]B�dB�dB�eB�dB�dB�^B�^B�^B�_B�aB�_B�]B�XB�WB�WB�YB�XB�SB�SB�SB�TB�RB�NB�QB�QB�QB�RB�OB�LB�MB�MB�QB�QB�QB�SB�TB�RB�MB�QB�VB�SB�QB�RB�TB�TB�WB�WB�_B�WB�]B�]B�_B�_B�]B�_B�_B�eB�cB�cB�jB�jB�iB�jB�jB�jB�qB�kB�oB�kB�iB�jB�mB�hB�jB�iB�dB�dB�cB�dB�eB�jB�iB�cB�kB�iB�jB�kB�jB�jB�kB�kB�kB�kB�nB�oB�oB�qB�qB�oB�oB�rB�oB�oB�pB�pB�pB�pB�pB�pB�pB�pB�pB�vB�vB�vB�tB�xB�uB�wB�xB�uB�uB�uB�sB�}B�vB�xB�}B�xB�vB�{B�}B�}B�zB�~B�|B�~B�|B�zB�|B�|B�{B�}B�B�}B�wB�{B�}B�|B�|B��B�~B�|B�|B�~B�|B�~B�}B�}B�}B�}B�{B�}B�}B�{B��B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111144111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.0045694 (+/- 0.01)                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040720220127170407  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144429  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144429  QCF$                G�O�G�O�G�O�0000000000004000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            A1��D̀ G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170407  IP  PSAL            A1��D̀ G�O�                