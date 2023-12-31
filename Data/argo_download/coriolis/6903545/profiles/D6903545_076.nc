CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:30Z creation; 2023-08-05T07:55:35Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_035h      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8    FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8(   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    88   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8H   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8P   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9    	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9,   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    90   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     94   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9T   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9t   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Bl   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  DT   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Ud   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  \�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ft   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  h\   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  o�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  w�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  yl   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �(   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �8   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �<   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �L   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �P   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �T   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �X   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �|   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200829092230  20230805075535  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               LA   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @�"x��>�1   @�"y^o��@P�Mc�w�8ݪۙ��1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 6.7 dbar]                                                      @�33A��A��AffA0  A>ffAQ��A`  Aq��A���A�  A���A���A�  A���A���A�  A�  A���A���A�  A���A�  A�  A�  A�33B33BffB  B��BffBffBffB ffB$  B(ffB,  B0ffB4  B8  B;��B@  BD  BH  BL��BO33BS��BW��B[��B_��Bc��BhffBk��BpffBs��Bx  B|ffB�33B���B���B���B���B���B�  B�  B���B���B�ffB�33B�33B�33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B�  B�33B�  B���B�  B���B�ffB���B�33B�33B���B���B���B�  B�  B�33B�ffB���B�ffB���B���B�ffB�  CffCL�CffCL�C	L�CL�CL�C33C33CL�CffCffC� C33CffC��C!L�C#�C%ffC'��C)ffC+33C-�C/  C1ffC3��C5�3C7� C9��C;� C=� C?� CAffCCffCE� CG��CIL�CKL�CM� CO33CQL�CSffCU�CWL�CY� C[33C]ffC_��CaffCc�CeffCg��CiL�Ck�CmffCo��CqffCs33CuffCw� CyL�C{�C}ffC� C��fC���C��fC��3C���C��3C���C��fC���C��fC�� C��fC���C��3C���C��3C���C�� C�ٚC�� C��3C���C���C��3C���C�� C��fC���C��3C�ٚC�� C��3C���C���C��3C�ٚC���C��3C��fC���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��3C��fC�ٚC���C�� C��3C��3C��fC��fC��fC���C���C���C���C�CÀ CĀ Cų3C��3C��fC��fC�ٚC�ٚC�ٚC�ٚC��fC��fC��fC��fC��fC��fC��3C�� CՀ C֌�Cי�CئfCٳ3C�� C�ٚC��fC��3C�� Cߌ�C���C�� C���C�fC�fC�3C�3C�� C�� C���CꙚC�fC�3C��3C�� C�� C���C�fC�fC�3C�� C���C��fC��fC�� C���C�s3C���D &fDl�D�3D  DL�Dy�D��D	fD
@ D� D� D  DFfDs3D�fD��DL�D�fD� D��D33Dl�D��D��D@ Dy�D ��D!��D#L�D$��D%��D&��D(33D)� D*��D+�3D-33D.l�D/�fD1fD2,�D3s3D4��D6  D7FfD8y�D9��D:�fD<&fD=` D>�fD?��DAFfDB�fDC�3DD�fDF9�DG��DH��DJ�DKS3DL� DM�3DOfDP9�DQs3DR��DS��DU@ DV�fDW��DX�3DZ9�D[l�D\��D^fD_9�D`l�Da�fDb�3DdFfDel�Df��Dh  Di,�Djy�Dk�fDl�3Dn,�Do� Dp��DrfDs9�Dts3Du��Dv�fDx  Dy` Dz�fD{��D}L�D~�fD� D��3D�&fD���D�VfD���D��3D�@ D�ٚD��3D�#3D�� D�` D�  D��3D�9�D���D��fD�  D���D�S3D��3D�� D�33D��3D�vfD��D��fD�c3D���D���D�6fD��3D�s3D�fD���D�\�D�  D��fD�@ D��fD�� D��D��fD�c3D�3D��3D�FfD�ٚD�� D�&fD���D�S3D���D��3D�@ D�ٚD�vfD�3D��3D�S3D��3D��fD�9�D��fD��3D��D���D�c3D�  D���D�<�D�� D�y�D�3D���D�Y�D���D���D�@ D��3D��fD��D��fD�VfD���D�c311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @�33A��A��AffA0  A>ffAQ��A`  Aq��A���A�  A���A���A�  A���A���A�  A�  A���A���A�  A���A�  A�  A�  A�33B33BffB  B��BffBffBffB ffB$  B(ffB,  B0ffB4  B8  B;��B@  BD  BH  BL��BO33BS��BW��B[��B_��Bc��BhffBk��BpffBs��Bx  B|ffB�33B���B���B���B���B���B�  B�  B���B���B�ffB�33B�33B�33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B�  B�33B�  B���B�  B���B�ffB���B�33B�33B���B���B���B�  B�  B�33B�ffB���B�ffB���B���B�ffB�  CffCL�CffCL�C	L�CL�CL�C33C33CL�CffCffC� C33CffC��C!L�C#�C%ffC'��C)ffC+33C-�C/  C1ffC3��C5�3C7� C9��C;� C=� C?� CAffCCffCE� CG��CIL�CKL�CM� CO33CQL�CSffCU�CWL�CY� C[33C]ffC_��CaffCc�CeffCg��CiL�Ck�CmffCo��CqffCs33CuffCw� CyL�C{�C}ffC� C��fC���C��fC��3C���C��3C���C��fC���C��fC�� C��fC���C��3C���C��3C���C�� C�ٚC�� C��3C���C���C��3C���C�� C��fC���C��3C�ٚC�� C��3C���C���C��3C�ٚC���C��3C��fC���C���C��3C��fC�ٚC���C�� C��3C��fC���C���C�� C��3C��3C��fC�ٚC���C�� C��3C��3C��fC��fC��fC���C���C���C���C�CÀ CĀ Cų3C��3C��fC��fC�ٚC�ٚC�ٚC�ٚC��fC��fC��fC��fC��fC��fC��3C�� CՀ C֌�Cי�CئfCٳ3C�� C�ٚC��fC��3C�� Cߌ�C���C�� C���C�fC�fC�3C�3C�� C�� C���CꙚC�fC�3C��3C�� C�� C���C�fC�fC�3C�� C���C��fC��fC�� C���C�s3C���D &fDl�D�3D  DL�Dy�D��D	fD
@ D� D� D  DFfDs3D�fD��DL�D�fD� D��D33Dl�D��D��D@ Dy�D ��D!��D#L�D$��D%��D&��D(33D)� D*��D+�3D-33D.l�D/�fD1fD2,�D3s3D4��D6  D7FfD8y�D9��D:�fD<&fD=` D>�fD?��DAFfDB�fDC�3DD�fDF9�DG��DH��DJ�DKS3DL� DM�3DOfDP9�DQs3DR��DS��DU@ DV�fDW��DX�3DZ9�D[l�D\��D^fD_9�D`l�Da�fDb�3DdFfDel�Df��Dh  Di,�Djy�Dk�fDl�3Dn,�Do� Dp��DrfDs9�Dts3Du��Dv�fDx  Dy` Dz�fD{��D}L�D~�fD� D��3D�&fD���D�VfD���D��3D�@ D�ٚD��3D�#3D�� D�` D�  D��3D�9�D���D��fD�  D���D�S3D��3D�� D�33D��3D�vfD��D��fD�c3D���D���D�6fD��3D�s3D�fD���D�\�D�  D��fD�@ D��fD�� D��D��fD�c3D�3D��3D�FfD�ٚD�� D�&fD���D�S3D���D��3D�@ D�ٚD�vfD�3D��3D�S3D��3D��fD�9�D��fD��3D��D���D�c3D�  D���D�<�D�� D�y�D�3D���D�Y�D���D���D�@ D��3D��fD��D��fD�VfD���D�c311111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��>߾w>��7�Z��hr��/��R񪾒n����𾄛��J�%�T�M��N��7KǾj~��׍P���T��Q�p�����ÿ
���	�^��p���-�.���o��1'�bM�>>v�>�{?E�?*~�>.{=� ſt��PbN��?}��=q��hs��G���t���ff���������Ϳ�t���󶿢�\��|��?}��G��͑h��{�Η���녿�t��ӕ����Ͽա˿���E�������ۿ͑h���������A���^5��7L���#��C���������dZ�����C����翹7L��Q쿸Q쿽�h�������D�������������|/���^�����\���w��"ѿ�A���X�~�R�y���h�ÿ]�-�Y���N{�@A��5��,�Ϳ%��X��˾�S���=q���|푾m�h�dZ�Q녾1&�ě��C�<�h=�+=��#>)��>R�>vȴ>�>��>�?}>�33>�|�>�t�>���?
~�?&�?33??�H?,1?(1'?(r�?*~�?*��?-�h?;dZ?I��?R-?R�?VE�?_;d?co?\(�?G+?4��?3t�?2�?49X?6ȴ?:��?>�R?Gl�?["�?g�?e�T?r-?���?��?�v�?��7?�33?�n�?��y?���?���?��?�~�?��?�`B?�z�?���?�`B?���?�t�?��j?��?�(�?��?�7L?z��?vȴ?q��?q�?s��?{��?{"�?z�?xQ�?r�?g+?]/?O�;?F��?E�?Q&�?T��?V?T�j?H�9?H��?A�7?B��?D�?A��?>��?>�R?>�R??�w?@A�?@Ĝ?@Ĝ?A%?AG�?AG�?A%?@Ĝ?@  ?@Ĝ?=�-?�P?�?�!?�u?�F?ƨ?��?l�?��?��?�?Z?��?�7>�j>�Q�>���?M�>�v�>���>��>�l�>��
>���>���>�G�>߾w>ܬ>�"�>�"�>ۥ�>ۥ�>޸R>�G�>���>�S�>�S�>��
>��y>�l�>��>��>�ff>ڟ�>�n�>��`>��>�|�>�r�>��R>���>���>��u>�t�>�hs>�C�>���>��>��>��>��>�$�>���>��>��>���>��>��>�J>e`B>A�7>1&�>%�T>��>��>1&�>;dZ>1&�>+>:^5>B�\>D��>@�>0 �>#�
>#�
>!��>C��>J��>cS�>~��>��7>�J>x��>r�!>m�h>hr�>hr�>j~�>j~�>ix�>m�h>n��>cS�>["�>bM�>z�H>r�!>["�>�C�>��/>�r�>��>���>��P>��>��>��`>�\)>��>�+>��^>��>���>��\>z�H>s�F>n��>ix�>`A�>T��>P�`>I�^>E��>A�7>7K�>5?}>1&�>"��>�->�+>I�>+>�=��=�/=ě�=���=�C�=��=�o=u=Y�=49X<�h<��
<T��;��
��o�49X���
���ͽ\)�\)�'P�`��\)��{�\����G���h�����#�����+�
=q��+��R�&�y�1&�<j�D���H�9�O�;�V�Y��bMӾk��p�׾y�#�~�۾��7��+��=q��O߾�I���O߾�hs��t����+�����(����R��MӾ�ff��xվ��D�����&龳33���F����Q쾺^5��j��p���%����ȴ9��=q��ƨ��O߾��`��z�����׍P�ڟ��ܬ��/�ݲ-��;d�߾w��G���Z���y��~����׾�?}��X��j���۾�|�   ��|��|��|��|��|��|��|���۾��۾�vɾ�vɾ�vɾ�v�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   >߾w>��7�Z��hr��/��R񪾒n����𾄛��J�%�T�M��N��7KǾj~��׍P���T��Q�p�����ÿ
���	�^��p���-�.���o��1'�bM�>>v�>�{?E�?*~�>.{=� ſt��PbN��?}��=q��hs��G���t���ff���������Ϳ�t���󶿢�\��|��?}��G��͑h��{�Η���녿�t��ӕ����Ͽա˿���E�������ۿ͑h���������A���^5��7L���#��C���������dZ�����C����翹7L��Q쿸Q쿽�h�������D�������������|/���^�����\���w��"ѿ�A���X�~�R�y���h�ÿ]�-�Y���N{�@A��5��,�Ϳ%��X��˾�S���=q���|푾m�h�dZ�Q녾1&�ě��C�<�h=�+=��#>)��>R�>vȴ>�>��>�?}>�33>�|�>�t�>���?
~�?&�?33??�H?,1?(1'?(r�?*~�?*��?-�h?;dZ?I��?R-?R�?VE�?_;d?co?\(�?G+?4��?3t�?2�?49X?6ȴ?:��?>�R?Gl�?["�?g�?e�T?r-?���?��?�v�?��7?�33?�n�?��y?���?���?��?�~�?��?�`B?�z�?���?�`B?���?�t�?��j?��?�(�?��?�7L?z��?vȴ?q��?q�?s��?{��?{"�?z�?xQ�?r�?g+?]/?O�;?F��?E�?Q&�?T��?V?T�j?H�9?H��?A�7?B��?D�?A��?>��?>�R?>�R??�w?@A�?@Ĝ?@Ĝ?A%?AG�?AG�?A%?@Ĝ?@  ?@Ĝ?=�-?�P?�?�!?�u?�F?ƨ?��?l�?��?��?�?Z?��?�7>�j>�Q�>���?M�>�v�>���>��>�l�>��
>���>���>�G�>߾w>ܬ>�"�>�"�>ۥ�>ۥ�>޸R>�G�>���>�S�>�S�>��
>��y>�l�>��>��>�ff>ڟ�>�n�>��`>��>�|�>�r�>��R>���>���>��u>�t�>�hs>�C�>���>��>��>��>��>�$�>���>��>��>���>��>��>�J>e`B>A�7>1&�>%�T>��>��>1&�>;dZ>1&�>+>:^5>B�\>D��>@�>0 �>#�
>#�
>!��>C��>J��>cS�>~��>��7>�J>x��>r�!>m�h>hr�>hr�>j~�>j~�>ix�>m�h>n��>cS�>["�>bM�>z�H>r�!>["�>�C�>��/>�r�>��>���>��P>��>��>��`>�\)>��>�+>��^>��>���>��\>z�H>s�F>n��>ix�>`A�>T��>P�`>I�^>E��>A�7>7K�>5?}>1&�>"��>�->�+>I�>+>�=��=�/=ě�=���=�C�=��=�o=u=Y�=49X<�h<��
<T��;��
��o�49X���
���ͽ\)�\)�'P�`��\)��{�\����G���h�����#�����+�
=q��+��R�&�y�1&�<j�D���H�9�O�;�V�Y��bMӾk��p�׾y�#�~�۾��7��+��=q��O߾�I���O߾�hs��t����+�����(����R��MӾ�ff��xվ��D�����&龳33���F����Q쾺^5��j��p���%����ȴ9��=q��ƨ��O߾��`��z�����׍P�ڟ��ܬ��/�ݲ-��;d�߾w��G���Z���y��~����׾�?}��X��j���۾�|�   ��|��|��|��|��|��|��|���۾��۾�vɾ�vɾ�vɾ�v�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB;dB<jB��B��B��BVB%�BJ�B�B��BB9XB^5BhsB�1B��B�B0!B��B�XB��B,B1'BK�BQ�BffB�oB�#BB�B�NB�B
=BuBhB�BDB�B+B8RBJ�BF�B@�BM�Br�B�7B~�B�B��B��B��B�3B��B��B��B��B�B�B�#B�)B�5B�NB�mB�B��B��B��B	%B	DB	JB	DB	oB	{B	�B	�B	�B	�B	�B	�B	�B	�B	%�B	-B	8RB	>wB	C�B	J�B	O�B	S�B	[#B	`BB	iyB	p�B	y�B	�DB	�uB	��B	�B	�dB	ƨB	ɺB	��B	�)B	�yB	��B
B
hB
�B
/B
9XB
G�B
K�B
M�B
O�B
Q�B
ZB
e`B
r�B
y�B
�1B
�oB
��B
��B
��B
��B
�3B
�LB
�jB
B
ǮB
��B
�B
�/B
�;B
�BB
�ZB
�B
�B
�B
�B
��B
��BBB
=BJBhBuB{BuBuBPB\BhBoB�B�B�B"�B-B/B33B8RBB�BF�BJ�BL�BQ�BN�BN�BW
BW
BXBT�BT�BL�BL�BM�BQ�BN�BQ�BS�BYBYBZBYBR�BT�BR�BP�BH�BZBYBXBXB\)BXBS�BO�BL�BN�BS�BT�BT�BR�BS�BT�BVBS�BT�BVBXBYBZB]/B]/B\)B^5B^5B]/B^5B_;B_;B^5B]/B_;B[#BS�BR�BW
BZBZBXBW
BXBXBXBXBW
BW
BYBZB\)B^5B_;B`BB]/B^5B]/B^5B]/B]/B]/B_;B`BB`BBaHBaHBcTBdZBe`BffBgmBgmBgmBhsBhsBhsBiyBiyBiyBhsBe`BdZBdZBdZBdZBe`BdZBdZBdZBdZBe`BffBgmBgmBgmBhsBiyBiyBiyBiyBjBjBjBiyBhsBiyBiyBiyBjBk�Bo�Bp�Bp�Br�Bt�Bt�Bt�Bt�Bt�Bu�Bw�Bz�B{�B� B�B�B�B�B�B�%B�B�%B�+B�+B�7B�=B�=B�=B�=B�DB�bB�\B�\B�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B;dB<jB��B��B��BVB%�BJ�B�B��BB9XB^5BhsB�1B��B�B0!B��B�XB��B,B1'BK�BQ�BffB�oB�#BB�B�NB�B
=BuBhB�BDB�B+B8RBJ�BF�B@�BM�Br�B�7B~�B�B��B��B��B�3B��B��B��B��B�B�B�#B�)B�5B�NB�mB�B��B��B��B	%B	DB	JB	DB	oB	{B	�B	�B	�B	�B	�B	�B	�B	�B	%�B	-B	8RB	>wB	C�B	J�B	O�B	S�B	[#B	`BB	iyB	p�B	y�B	�DB	�uB	��B	�B	�dB	ƨB	ɺB	��B	�)B	�yB	��B
B
hB
�B
/B
9XB
G�B
K�B
M�B
O�B
Q�B
ZB
e`B
r�B
y�B
�1B
�oB
��B
��B
��B
��B
�3B
�LB
�jB
B
ǮB
��B
�B
�/B
�;B
�BB
�ZB
�B
�B
�B
�B
��B
��BBB
=BJBhBuB{BuBuBPB\BhBoB�B�B�B"�B-B/B33B8RBB�BF�BJ�BL�BQ�BN�BN�BW
BW
BXBT�BT�BL�BL�BM�BQ�BN�BQ�BS�BYBYBZBYBR�BT�BR�BP�BH�BZBYBXBXB\)BXBS�BO�BL�BN�BS�BT�BT�BR�BS�BT�BVBS�BT�BVBXBYBZB]/B]/B\)B^5B^5B]/B^5B_;B_;B^5B]/B_;B[#BS�BR�BW
BZBZBXBW
BXBXBXBXBW
BW
BYBZB\)B^5B_;B`BB]/B^5B]/B^5B]/B]/B]/B_;B`BB`BBaHBaHBcTBdZBe`BffBgmBgmBgmBhsBhsBhsBiyBiyBiyBhsBe`BdZBdZBdZBdZBe`BdZBdZBdZBdZBe`BffBgmBgmBgmBhsBiyBiyBiyBiyBjBjBjBiyBhsBiyBiyBiyBjBk�Bo�Bp�Bp�Br�Bt�Bt�Bt�Bt�Bt�Bu�Bw�Bz�B{�B� B�B�B�B�B�B�%B�B�%B�+B�+B�7B�=B�=B�=B�=B�DB�bB�\B�\B�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755352023080507553520230805075535  IF  ARFMCODA035h                                                                20200829092230                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092341  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC4.6                                                                 20200829092341  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134655  IP  PSAL            @�33D�c3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172541  IP  PSAL            @�33D�c3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075535  IP  PSAL            @�33D�c3G�O�                