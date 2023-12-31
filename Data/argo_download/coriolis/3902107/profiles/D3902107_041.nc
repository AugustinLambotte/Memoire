CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     	N_HISTORY          N_CALIB          	   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-28T15:12:16Z creation; 2020-11-23T11:33:26Z last update (BSH ARSQ software)    
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
_FillValue                  ,  �HArgo profile    3.1 1.2 19500101000000  20200828151216  20220127170814  3902107 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               )A   IF                                  2C  D   ARVOR                           AI2600-18EU007                  5900A04                         844 @�,8�q�1   @�,8�q�@S�ƣ<�@zM~��p8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 8.6 dbar]                                                      A��A!��A1��A@  AQ��Aa��Aq��A~ffA�  A���A���A���A���A���A���A�  A���A���A���A���A���A�33A�  B ffB  BffBffBffBffBffBffB ffB$  B(  B,  B0��B3��B7��B<  B@ffBDffBHffBL  BO��BS��BX��B\  B`ffBd  Bg��Bk33Bo33Bs33Bx  B{��B33B�  B�  B�  B�  B�  B�33B�33B�33B���B���B�  B�  B�33B���B���B���B�  B�33B���B���B�33B���B�  B�33B���B�  B�ffB�  B�ffB���B�33B���B�  B�ffB�  BǙ�B���BΙ�Bҙ�B���B���B���B♚B晚BꙚBB�B���B���B���CffC� C� C�3C	L�CffC� C� C��CL�CffCffC� C��C33CL�C!ffC#ffC%��C'L�C)ffC+��C-L�C/� C1�3C3� C5L�C7�C9ffC;��C=� C?L�CA33CC  CEffCG�3CI�3CK��CM� CO� CQffCSffCUffCWffCYffC[ffC]ffC_ffCaffCc� Ce� Cg��Ci�3Ck��Cm�fCoffCq  Cs�CuL�CwffCy� C{�3C}�fC� C���C��fC��3C���C�ٚC��3C�� C���C��fC��3C���C�ٚC��3C�� C���C��3C�� C��fC��3C���C��fC�� C�ٚC��3C�� C���C��3C�ٚC��3C�� C���C��3C���C��fC�� C���C��fC�� C���C��fC��3C�� C���C��fC�� C���C�ٚC��fC��fC��3C�� C�� C���C���C���C���C��3C��fC���C���C���C��fC��fC��3C�� C���C�CæfC�� Cř�CƦfC�� Cș�CɦfC�� C˦fC̳3C���CΦfCό�Cг3C���Cҳ3Cә�C�� C�ٚC���C׳3Cؙ�C�� C�ٚC۳3Cܙ�Cݳ3C�� CߦfC�3C���C�fC�3C���C�3C�� C�ٚC�3C陚C�� C�ٚC�3C홚C�� C�ٚC���C�3C�fC��C��3C��fC�ٚC���C�� C��3C�ffC�ٚD &fDffD�fD�fD,�Ds3D��D��D
9�D�3D��D��D&fDs3D� D�D9�Dl�D�3D  D9�Dl�D�fD��D,�D�fD � D!��D#33D$s3D%�3D&�3D(33D)l�D*�fD,fD-@ D.y�D/��D0��D2@ D3�fD4�3D5�fD733D8�fD9� D;  D<@ D=� D>�fD@�DA@ DBl�DC��DEfDFFfDG�fDH�fDJ�DK33DLffDM�3DO  DP9�DQs3DR� DT  DU9�DVs3DW��DX�3DZ9�D[�fD\��D]�fD_33D`� Da�3Db�fDd33De��Df�fDg��Di33Djl�Dk�fDl�fDn  DoffDp�fDq� Ds9�Dty�Du� Dv�fDx&fDy� Dz�fD|  D}@ D~�fD�fD��3D�#3D��3D�\�D���D���D�6fD���D��fD�  D���D�\�D���D��3D�<�D��fD�vfD�fD���D�\�D�  D��3D�C3D��D�� D�3D��fD�ffD�3D���D�@ D��fD�� D��D�� D�c3D�fD���D�33D��fD�y�D��D��3D�ffD�  D��fD�9�D���D��3D�)�D�� D�VfD���D���D�C3D��D�� D��D��fD�ffD��fD��fD�9�D�� D��3D��D��3D�\�D�	�D��3D�@ D��fD�� D��D�� D�ffD�  D�� D�C3D��3D�|�D�fD���D�c3D���D���D�C3D���D�y�D�fD���D�\�D���D�� D�33D��fD�|�D�3D���D�c3D���D��fD�0 D�� D�|�D��D° D�\�D�	�DĦfD�FfD��fDƆfD��Dǰ D�\�D�3Dɠ D�#3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A��A!��A1��A@  AQ��Aa��Aq��A~ffA�  A���A���A���A���A���A���A�  A���A���A���A���A���A�33A�  B ffB  BffBffBffBffBffBffB ffB$  B(  B,  B0��B3��B7��B<  B@ffBDffBHffBL  BO��BS��BX��B\  B`ffBd  Bg��Bk33Bo33Bs33Bx  B{��B33B�  B�  B�  B�  B�  B�33B�33B�33B���B���B�  B�  B�33B���B���B���B�  B�33B���B���B�33B���B�  B�33B���B�  B�ffB�  B�ffB���B�33B���B�  B�ffB�  BǙ�B���BΙ�Bҙ�B���B���B���B♚B晚BꙚBB�B���B���B���CffC� C� C�3C	L�CffC� C� C��CL�CffCffC� C��C33CL�C!ffC#ffC%��C'L�C)ffC+��C-L�C/� C1�3C3� C5L�C7�C9ffC;��C=� C?L�CA33CC  CEffCG�3CI�3CK��CM� CO� CQffCSffCUffCWffCYffC[ffC]ffC_ffCaffCc� Ce� Cg��Ci�3Ck��Cm�fCoffCq  Cs�CuL�CwffCy� C{�3C}�fC� C���C��fC��3C���C�ٚC��3C�� C���C��fC��3C���C�ٚC��3C�� C���C��3C�� C��fC��3C���C��fC�� C�ٚC��3C�� C���C��3C�ٚC��3C�� C���C��3C���C��fC�� C���C��fC�� C���C��fC��3C�� C���C��fC�� C���C�ٚC��fC��fC��3C�� C�� C���C���C���C���C��3C��fC���C���C���C��fC��fC��3C�� C���C�CæfC�� Cř�CƦfC�� Cș�CɦfC�� C˦fC̳3C���CΦfCό�Cг3C���Cҳ3Cә�C�� C�ٚC���C׳3Cؙ�C�� C�ٚC۳3Cܙ�Cݳ3C�� CߦfC�3C���C�fC�3C���C�3C�� C�ٚC�3C陚C�� C�ٚC�3C홚C�� C�ٚC���C�3C�fC��C��3C��fC�ٚC���C�� C��3C�ffC�ٚD &fDffD�fD�fD,�Ds3D��D��D
9�D�3D��D��D&fDs3D� D�D9�Dl�D�3D  D9�Dl�D�fD��D,�D�fD � D!��D#33D$s3D%�3D&�3D(33D)l�D*�fD,fD-@ D.y�D/��D0��D2@ D3�fD4�3D5�fD733D8�fD9� D;  D<@ D=� D>�fD@�DA@ DBl�DC��DEfDFFfDG�fDH�fDJ�DK33DLffDM�3DO  DP9�DQs3DR� DT  DU9�DVs3DW��DX�3DZ9�D[�fD\��D]�fD_33D`� Da�3Db�fDd33De��Df�fDg��Di33Djl�Dk�fDl�fDn  DoffDp�fDq� Ds9�Dty�Du� Dv�fDx&fDy� Dz�fD|  D}@ D~�fD�fD��3D�#3D��3D�\�D���D���D�6fD���D��fD�  D���D�\�D���D��3D�<�D��fD�vfD�fD���D�\�D�  D��3D�C3D��D�� D�3D��fD�ffD�3D���D�@ D��fD�� D��D�� D�c3D�fD���D�33D��fD�y�D��D��3D�ffD�  D��fD�9�D���D��3D�)�D�� D�VfD���D���D�C3D��D�� D��D��fD�ffD��fD��fD�9�D�� D��3D��D��3D�\�D�	�D��3D�@ D��fD�� D��D�� D�ffD�  D�� D�C3D��3D�|�D�fD���D�c3D���D���D�C3D���D�y�D�fD���D�\�D���D�� D�33D��fD�|�D�3D���D�c3D���D��fD�0 D�� D�|�D��D° D�\�D�	�DĦfD�FfD��fDƆfD��Dǰ D�\�D�3Dɠ D�#3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���o��o�6E��(r��%�T��˾��y��J���`Ĝ�k�h�ÿpbN�b���Qhs�U?}�^�R�h�9��^5���
���˿�hs��n���o��񪿑��������7L�:���t�����=���=�=ě�=,1<�9X<�C�<�C���1�Ƨ��/<���<�C���j=�=��m>��>���?ȴ?A�7?YX?lI�?�?�v�?ȓu?��#?�Q�@J@�/@
�\@(�@�P@��@\)@7L@��@$I�@*�@+"�@-�@17L@1�@2�\@3"�@3��@3C�@2�@/��@.ff@,��@,9X@*�!@(�@)&�@(A�@&�y@&ff@%�h@$9X@#33@"J@"�@"J@!�^@#@"�@!&�@   @�P@�w@ r�@ �u@!G�@"=q@"�@!��@!�@!%@!�@"�@!��@"��@#"�@!��@ ��@��@l�@��@��@�D@(�@��@��@?}@?}@?}@/@V@�@��@z�@j@9X@��@�@C�@"�@j@�h@!%@$I�@$Z@$1@#��@#�m@$(�@$��@%/@%/@%�@&$�@%�@$�/@%/@$�j@$�@$�@$��@$��@%�@%?}@%?}@%/@$��@%V@$��@$z�@#��@#�@"�@"-@!x�@!�@ ��@ Q�@�w@�P@��@�-@�T@V@�j@I�@"�@��@��@~�@^5@��@�u@Q�@��@E�@�-@/@�D@�j@�j@z�@dZ@��@%@1'@|�@��@?}@l�@J@ A�?���?�C�?���?�?��?�;d?�O�?��?�?��?陚?�r�?���?�ff?�J?� �?��?��?ݑh?ܬ?��?���?׮?�+?�+?�+?�$�?ա�?Ձ?Ձ?��?��/?�9X?��
?Ұ!?���?�G�?�G�?�%?�%?�  ?�|�?��;?Ͼw?�\)?Ͳ-?�I�?�j?�(�?�(�?�dZ?�"�?ʟ�?�ȴ?�?�z�?��?�J?���?�hs?��7?�hs?�Ĝ?�Ĝ?��?��;?�\)?�{?���?�b?���?�z�?�S�?��!?�M�?�n�?�J?��`?��`?���?�/?�/?�O�?��?�X?��?��j?�J?�  ?�Ĝ?�|�?���?���?��?�$�?�&�?��?���?���?��+?��
?�n�?��`?�n�?�$�?��+?�V?��-?��?�r�?�?��j?�n�?|j?w
=?p�`?j~�?f$�?`Ĝ?XQ�?Qhs?D�/?A%?;�m?0��?&��?#o?v�?�+?��?��>��>�1>߾w>Ǯ>�%>�^5>��!>���>�{>���>���>��>��^>��\>|�>q��>F��>6E�>&�y>V=ě�=��T=L��=+<�/<���<���<��
<T��:�o���
�49X�o���',1�D����%��Ƨ��"ѽ���%�1'�\)�t����)��0 žB�\�L�;Z��aG��q���vȴ�|푾�J��$ݾ��^��\)��bN��t���zᾘb������(���;d��G���S���S����
���/��ff��~���V��-���پ��#���H���H��  ����Ƨ��7L������Ͼ�b�ݲ-��G�������MӾ�S����/���xվ�xվ���h��{��!��E����#����vɿ A��%����/�ff��9�	�^�C���Ϳ��V������녿9X�����ȴ��ٿ������^5�푿�R�   � A��   � Ĝ�!G��"Mӿ#S��$Z�%`B�&�y�(1'�)7L�)��*���*���*���+ƨ�,�Ϳ.V�.���/\)�0bN�2-�333�49X�4�j�5�6�+�7
=�7�P�8Q�8���9��9X�9���:��:�H�;"ѿ;dZ�;dZ�;��<j�<푿=/�=�-�=�>v�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111�o��o�6E��(r��%�T��˾��y��J���`Ĝ�k�h�ÿpbN�b���Qhs�U?}�^�R�h�9��^5���
���˿�hs��n���o��񪿑��������7L�:���t�����=���=�=ě�=,1<�9X<�C�<�C���1�Ƨ��/<���<�C���j=�=��m>��>���?ȴ?A�7?YX?lI�?�?�v�?ȓu?��#?�Q�@J@�/@
�\@(�@�P@��@\)@7L@��@$I�@*�@+"�@-�@17L@1�@2�\@3"�@3��@3C�@2�@/��@.ff@,��@,9X@*�!@(�@)&�@(A�@&�y@&ff@%�h@$9X@#33@"J@"�@"J@!�^@#@"�@!&�@   @�P@�w@ r�@ �u@!G�@"=q@"�@!��@!�@!%@!�@"�@!��@"��@#"�@!��@ ��@��@l�@��@��@�D@(�@��@��@?}@?}@?}@/@V@�@��@z�@j@9X@��@�@C�@"�@j@�h@!%@$I�@$Z@$1@#��@#�m@$(�@$��@%/@%/@%�@&$�@%�@$�/@%/@$�j@$�@$�@$��@$��@%�@%?}@%?}@%/@$��@%V@$��@$z�@#��@#�@"�@"-@!x�@!�@ ��@ Q�@�w@�P@��@�-@�T@V@�j@I�@"�@��@��@~�@^5@��@�u@Q�@��@E�@�-@/@�D@�j@�j@z�@dZ@��@%@1'@|�@��@?}@l�@J@ A�?���?�C�?���?�?��?�;d?�O�?��?�?��?陚?�r�?���?�ff?�J?� �?��?��?ݑh?ܬ?��?���?׮?�+?�+?�+?�$�?ա�?Ձ?Ձ?��?��/?�9X?��
?Ұ!?���?�G�?�G�?�%?�%?�  ?�|�?��;?Ͼw?�\)?Ͳ-?�I�?�j?�(�?�(�?�dZ?�"�?ʟ�?�ȴ?�?�z�?��?�J?���?�hs?��7?�hs?�Ĝ?�Ĝ?��?��;?�\)?�{?���?�b?���?�z�?�S�?��!?�M�?�n�?�J?��`?��`?���?�/?�/?�O�?��?�X?��?��j?�J?�  ?�Ĝ?�|�?���?���?��?�$�?�&�?��?���?���?��+?��
?�n�?��`?�n�?�$�?��+?�V?��-?��?�r�?�?��j?�n�?|j?w
=?p�`?j~�?f$�?`Ĝ?XQ�?Qhs?D�/?A%?;�m?0��?&��?#o?v�?�+?��?��>��>�1>߾w>Ǯ>�%>�^5>��!>���>�{>���>���>��>��^>��\>|�>q��>F��>6E�>&�y>V=ě�=��T=L��=+<�/<���<���<��
<T��:�o���
�49X�o���',1�D����%��Ƨ��"ѽ���%�1'�\)�t����)��0 žB�\�L�;Z��aG��q���vȴ�|푾�J��$ݾ��^��\)��bN��t���zᾘb������(���;d��G���S���S����
���/��ff��~���V��-���پ��#���H���H��  ����Ƨ��7L������Ͼ�b�ݲ-��G�������MӾ�S����/���xվ�xվ���h��{��!��E����#����vɿ A��%����/�ff��9�	�^�C���Ϳ��V������녿9X�����ȴ��ٿ������^5�푿�R�   � A��   � Ĝ�!G��"Mӿ#S��$Z�%`B�&�y�(1'�)7L�)��*���*���*���+ƨ�,�Ϳ.V�.���/\)�0bN�2-�333�49X�4�j�5�6�+�7
=�7�P�8Q�8���9��9X�9���:��:�H�;"ѿ;dZ�;dZ�;��<j�<푿=/�=�-�=�>v�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B�`Bq�B�B�B��BBQ�B�?B	7BcTB��B�;B9XBdZBu�B��B�/B�5B9XB^5Bv�B��B�B�'B��B��B�BBW
B��B�3B�HB�B��B		7B	DB	�B	$�B	2-B	+B	33B	E�B	M�B	W
B	y�B	�DB	��B	�B	�qB	�TB	��B
\B
>wB
YB
s�B
��B
�RB
��B
�B
�`B
��BBuB�B!�B2-BC�BT�BW
BaHBk�Bn�Bp�Bt�Bx�By�B~�B|�Bz�B{�B� B}�B}�B� B�B�B�B�B�%B�B�B�B�%B�1B�DB�DB�=B�7B�DB�JB�PB�VB�hB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�!B�'B�-B�-B�?B�XB��BŢBŢBƨBƨBǮBǮB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BƨBǮBȴBɺBɺB��BɺBǮBǮBĜBÖBB��B�wB�jB�-B�3B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B�B�B�FB�FB�LB�FB�FB�LB�LB�RB�RB�LB�FB�?B�?B�?B�9B�9B�3B�-B�'B�!B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B��B�`Bq�B�B�B��BBQ�B�?B	7BcTB��B�;B9XBdZBu�B��B�/B�5B9XB^5Bv�B��B�B�'B��B��B�BBW
B��B�3B�HB�B��B		7B	DB	�B	$�B	2-B	+B	33B	E�B	M�B	W
B	y�B	�DB	��B	�B	�qB	�TB	��B
\B
>wB
YB
s�B
��B
�RB
��B
�B
�`B
��BBuB�B!�B2-BC�BT�BW
BaHBk�Bn�Bp�Bt�Bx�By�B~�B|�Bz�B{�B� B}�B}�B� B�B�B�B�B�%B�B�B�B�%B�1B�DB�DB�=B�7B�DB�JB�PB�VB�hB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�!B�'B�-B�-B�?B�XB��BŢBŢBƨBƨBǮBǮB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BƨBǮBȴBɺBɺB��BɺBǮBǮBĜBÖBB��B�wB�jB�-B�3B�B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B�B�B�FB�FB�LB�FB�FB�LB�LB�RB�RB�LB�FB�?B�?B�?B�9B�9B�3B�-B�'B�!B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               No adjustment was necessary. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                                                      202011231133262022012717081420220127170814  IF  ARFMCODA035h                                                                20200828151216                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828151254  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828151254  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201123113326  IP  PSAL            A��D�#3G�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170814  IP  PSAL            A��D�#3G�O�                