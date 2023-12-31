CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-08-29T09:22:29Z creation; 2023-08-05T07:55:34Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z           :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z           C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  J�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���        L�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        S�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Z�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        \�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  c�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        e�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        l�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        u�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  |�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o        ~�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �T   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �d   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �h   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �x   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �|   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200829092229  20230805075534  6903545 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               :A   IF                                  2C  D   ARVOR                           AI2600-18EU001                  5900A04                         844 @��wl�l1   @��wl�l@R1H��4��.������8   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      A1��A<��ANffAa��AnffA���A�33A���A�33A�ffA���A���A���A�ffA�ffA�ffA�ffA�ffA�33A�  A�  A�33B��B��B��B��B��B��B  B   B$ffB(  B,ffB0  B4  B7��B<  B?��BD  BH  BLffBP  BT  BX  B\  B`  Bd  Bh  BlffBo33Bs��Bx  B|  B�33B���B���B�33B���B�  B�33B���B�ffB���B�33B���B�ffB���B�33B���B�ffB���B�ffB�33B���B���B�33B���B�33B�  B���B�ffB���B�ffB�33B���B���B�ffB���Bƙ�B�ffB���B�ffB�  B�33B���B�ffB♚B�  BꙚB���B�33B���B�ffB���C��C�3C�3C�3C	��C�fC��C�3C�3C��C��C��C��C��C��C�3C!�3C#��C%ffC'  C)�C+L�C-ffC/� C1�3C3��C5ffC7  C9�C;33C=ffC?� CA��CC��CEffCG  CIL�CKffCMffCO� CQ��CS��CUL�CWffCY� C[��C]��C_L�CaffCcffCeffCg� Ci� Ck��CmL�CoffCq��CsL�Cu� Cw�3CyffC{33C}ffC��C��3C���C��3C�� C��fC��3C���C��fC�� C�ٚC�� C��fC��3C���C��3C���C��3C���C��3C���C��3C���C��fC���C��3C���C��3C��fC���C��3C���C��3C���C��3C�ٚC��3C��fC���C��3C�ٚC�� C��fC���C��3C�ٚC���C��3C���C���C��3C���C��3C���C�� C�ٚC�� C��fC���C��3C�ٚC�� C��3C���C�� C�ٚC�� C¦fCÌ�Cĳ3C���CƳ3CǙ�CȌ�Cɳ3C���C�� C̦fC͌�Cγ3C�ٚC���Cѳ3CҦfCә�CԌ�Cճ3C��fC�ٚC���C���C�� C�� Cܳ3Cݳ3C޳3C߳3C�3C�3C�3C�3C�� C�� C�� C���C���C�ٚC��fC��fC��3C��3C�� C� C�� C�3C�� C�� C�� C�� C���C���C�ٚC��3C�Y�C�ٚD 33Dy�D�fD�3D,�D� D��D	  D
@ Dl�D��D�D@ Dy�D��D  DL�Dy�D��D� D9�D�3D�3D3DS3D��D ��D!��D#L�D$� D%��D&�3D(,�D)l�D*��D+��D-33D.� D/��D1  D29�D3s3D4�fD6fD7FfD8l�D9�3D;  D<,�D=ffD>��D@3DAL�DB�fDC��DD�3DF,�DGffDH� DI� DK  DL` DM� DN�fDP33DQ� DR�fDT3DU9�DVl�DW��DX� DZ9�D[�3D\� D]��D_9�D`�fDa��Db�3Dd@ Dey�Df�3Dg��Di,�Djl�Dk�3Dl��Dn9�Do� Dp��Dr  Ds@ Dt�fDu�3Dv��DxFfDyl�Dz��D|fD}@ D~s3D�fD�|�D�&fD��fD�ffD�fD�� D�9�D�ٚD�|�D�#3D�� D�\�D���D�� D�C3D��fD�|�D�fD��3D�c3D�  D��3D�33D��fD�y�D�  D��fD�\�D��fD�� D�I�D��3D�� D�  D�� D�` D�3D��fD�9�D�� D�s3D�fD��fD�Y�D���D�� D�<�D�ٚD�vfD�3D��3D�VfD���D�� D�9�D��3D�|�D��D���D�` D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A1��A<��ANffAa��AnffA���A�33A���A�33A�ffA���A���A���A�ffA�ffA�ffA�ffA�ffA�33A�  A�  A�33B��B��B��B��B��B��B  B   B$ffB(  B,ffB0  B4  B7��B<  B?��BD  BH  BLffBP  BT  BX  B\  B`  Bd  Bh  BlffBo33Bs��Bx  B|  B�33B���B���B�33B���B�  B�33B���B�ffB���B�33B���B�ffB���B�33B���B�ffB���B�ffB�33B���B���B�33B���B�33B�  B���B�ffB���B�ffB�33B���B���B�ffB���Bƙ�B�ffB���B�ffB�  B�33B���B�ffB♚B�  BꙚB���B�33B���B�ffB���C��C�3C�3C�3C	��C�fC��C�3C�3C��C��C��C��C��C��C�3C!�3C#��C%ffC'  C)�C+L�C-ffC/� C1�3C3��C5ffC7  C9�C;33C=ffC?� CA��CC��CEffCG  CIL�CKffCMffCO� CQ��CS��CUL�CWffCY� C[��C]��C_L�CaffCcffCeffCg� Ci� Ck��CmL�CoffCq��CsL�Cu� Cw�3CyffC{33C}ffC��C��3C���C��3C�� C��fC��3C���C��fC�� C�ٚC�� C��fC��3C���C��3C���C��3C���C��3C���C��3C���C��fC���C��3C���C��3C��fC���C��3C���C��3C���C��3C�ٚC��3C��fC���C��3C�ٚC�� C��fC���C��3C�ٚC���C��3C���C���C��3C���C��3C���C�� C�ٚC�� C��fC���C��3C�ٚC�� C��3C���C�� C�ٚC�� C¦fCÌ�Cĳ3C���CƳ3CǙ�CȌ�Cɳ3C���C�� C̦fC͌�Cγ3C�ٚC���Cѳ3CҦfCә�CԌ�Cճ3C��fC�ٚC���C���C�� C�� Cܳ3Cݳ3C޳3C߳3C�3C�3C�3C�3C�� C�� C�� C���C���C�ٚC��fC��fC��3C��3C�� C� C�� C�3C�� C�� C�� C�� C���C���C�ٚC��3C�Y�C�ٚD 33Dy�D�fD�3D,�D� D��D	  D
@ Dl�D��D�D@ Dy�D��D  DL�Dy�D��D� D9�D�3D�3D3DS3D��D ��D!��D#L�D$� D%��D&�3D(,�D)l�D*��D+��D-33D.� D/��D1  D29�D3s3D4�fD6fD7FfD8l�D9�3D;  D<,�D=ffD>��D@3DAL�DB�fDC��DD�3DF,�DGffDH� DI� DK  DL` DM� DN�fDP33DQ� DR�fDT3DU9�DVl�DW��DX� DZ9�D[�3D\� D]��D_9�D`�fDa��Db�3Dd@ Dey�Df�3Dg��Di,�Djl�Dk�3Dl��Dn9�Do� Dp��Dr  Ds@ Dt�fDu�3Dv��DxFfDyl�Dz��D|fD}@ D~s3D�fD�|�D�&fD��fD�ffD�fD�� D�9�D�ٚD�|�D�#3D�� D�\�D���D�� D�C3D��fD�|�D�fD��3D�c3D�  D��3D�33D��fD�y�D�  D��fD�\�D��fD�� D�I�D��3D�� D�  D�� D�` D�3D��fD�9�D�� D�s3D�fD��fD�Y�D���D�� D�<�D�ٚD�vfD�3D��3D�VfD���D�� D�9�D��3D�|�D��D���D�` D��3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����n���o��V��(����^�������D���H������1'���忌�D�e`B�333�+ƨ��#�&��ƨ����>vɾs�F��녾����t�������n����^����~�۾cS��["ѾN��:^5�\=q��;ě��q��<D��<���=���=��=�
=>A�7>W
=>���>�n�>���>���>���>�?}>ۥ�>�/>�l�>�M�>�G�>�G�>ݲ->�$�>š�>Ǯ>��>ܬ>��T>��?J?r�?n�?�?:^5?H1'?d�/?y�#?�G�?�K�?�/?�C�?�C�?�^5?���?�=q?�^5?�K�?�^5?��+?�?�?��?�$�?�ff?��?�\)?���?�9X?�J?��?�(�?���?z^5?v?|(�?��?�C�?��?��?�hs?��H?�hs?�I�@��@A�@5?@?}@/@�-@v�@�w@bN@7L@��@)�7@*�!@(��@'�w@!%@�7@�@ȴ@Z@
=@\)@A�@	��@dZ@�@v�@�u@�-@
��@
�@��@ƨ@%?�v�?��?��H?���?��R@t�@��@�@j@Z@=q?��;?�o?�u?��y?��?ߝ�?׮?�J?ӶF?��
?��/?Լj?ԛ�?��?ԛ�?�9X?�9X?��
?ҏ\?У�?ϝ�?θR?�V?�^5?�9X?�bN?�5??���?��?�^5?�X?��?��+?�K�?��P?�+?��y?�?���?��j?�Z?��
?�S�?���?�v�?� �?�Ĝ?���?�v�?�O�?�V?��?��?�dZ?���?���?�E�?��?�I�?�E�?�Z?�-?���?��?�V?�/?�;d?���?��R?�^5?�?���?���?{dZ?nV?i�^?h�9?g�?l1?q��?y�?xQ�?t9X?m�h?g�?d��?]�-?O�?H�9?Gl�?F$�?@A�?@  ?=�-?;��?8��?7�P?7�P?6ȴ?4z�? A�?��?�?	x�?
��?�?33?bN?�?��?r�?{?|�?!%?!��?#��?%��?&��?'+?)7L?,�D?&��?&�y?"M�? Ĝ?|�?/?�w?"��?#o?!��? Ĝ?v�?dZ?"�?dZ?"�?j?v�?(�?�P?hs?��?9X?�`?
��?�
>�^5>�ȴ>�9X>�Z>Ձ>���>�dZ>���>�ff>�b>�n�>�bN>�bN>�=q>���>�  >vȴ>k�>cS�>Y�>W
=>R�>I�^>8Q�>-V>(��>!��>�P>
=q>1'>�>J=��#=�G�=���=�Q�=� �=��
=�t�=Y�=8Q�=m�h=�t�=�7L=]/;�`B�o�H�9�q����%�ixսL�ͽ8Q�49X�,1�#�
�,1�L�ͽixս�o��o��7L���T��E���vɽ����vɽ�vɽ�vɽě��ȴ9�ȴ9�������ͽ����
=��G���`B���ٽ��m�o�o�
=q�n���+���#�
�$�/�+�2-�:^5�C���J���Kƨ�V�^5?�dZ�fff�hr��gl��hr��j~��m�h�s�F�z�H��  ���7��o������$ݾ����=q������O߾����bN��녾�n���񪾓t����Ͼ��������+��
=���u���㾜����R��A���G���`B���þ�1����������33��9X���j���j��?}��E���KǾ��پ�Q쾸�����H��푾�p���󶾿|��%���7��+������111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111��n���o��V��(����^�������D���H������1'���忌�D�e`B�333�+ƨ��#�&��ƨ����>vɾs�F��녾����t�������n����^����~�۾cS��["ѾN��:^5�\G�O�G�O�G�O�G�O�<���=���=��=�
=>A�7>W
=>���>�n�>���>���>���>�?}>ۥ�>�/>�l�>�M�>�G�>�G�>ݲ->�$�>š�>Ǯ>��>ܬ>��T>��?J?r�?n�?�?:^5?H1'?d�/?y�#?�G�?�K�?�/?�C�?�C�?�^5?���?�=q?�^5?�K�?�^5?��+?�?�?��?�$�?�ff?��?�\)?���?�9X?�J?��?�(�?���?z^5?v?|(�?��?�C�?��?��?�hs?��H?�hs?�I�@��@A�@5?@?}@/@�-@v�@�w@bN@7L@��@)�7@*�!@(��@'�w@!%@�7@�@ȴ@Z@
=@\)@A�@	��@dZ@�@v�@�u@�-@
��@
�@��@ƨ@%?�v�?��?��H?���?��R@t�@��@�@j@Z@=q?��;?�o?�u?��y?��?ߝ�?׮?�J?ӶF?��
?��/?Լj?ԛ�?��?ԛ�?�9X?�9X?��
?ҏ\?У�?ϝ�?θR?�V?�^5?�9X?�bN?�5??���?��?�^5?�X?��?��+?�K�?��P?�+?��y?�?���?��j?�Z?��
?�S�?���?�v�?� �?�Ĝ?���?�v�?�O�?�V?��?��?�dZ?���?���?�E�?��?�I�?�E�?�Z?�-?���?��?�V?�/?�;d?���?��R?�^5?�?���?���?{dZ?nV?i�^?h�9?g�?l1?q��?y�?xQ�?t9X?m�h?g�?d��?]�-?O�?H�9?Gl�?F$�?@A�?@  ?=�-?;��?8��?7�P?7�P?6ȴ?4z�? A�?��?�?	x�?
��?�?33?bN?�?��?r�?{?|�?!%?!��?#��?%��?&��?'+?)7L?,�D?&��?&�y?"M�? Ĝ?|�?/?�w?"��?#o?!��? Ĝ?v�?dZ?"�?dZ?"�?j?v�?(�?�P?hs?��?9X?�`?
��?�
>�^5>�ȴ>�9X>�Z>Ձ>���>�dZ>���>�ff>�b>�n�>�bN>�bN>�=q>���>�  >vȴ>k�>cS�>Y�>W
=>R�>I�^>8Q�>-V>(��>!��>�P>
=q>1'>�>J=��#=�G�=���=�Q�=� �=��
=�t�=Y�=8Q�=m�h=�t�=�7L=]/;�`B�o�H�9�q����%�ixսL�ͽ8Q�49X�,1�#�
�,1�L�ͽixս�o��o��7L���T��E���vɽ����vɽ�vɽ�vɽě��ȴ9�ȴ9�������ͽ����
=��G���`B���ٽ��m�o�o�
=q�n���+���#�
�$�/�+�2-�:^5�C���J���Kƨ�V�^5?�dZ�fff�hr��gl��hr��j~��m�h�s�F�z�H��  ���7��o������$ݾ����=q������O߾����bN��녾�n���񪾓t����Ͼ��������+��
=���u���㾜����R��A���G���`B���þ�1����������33��9X���j���j��?}��E���KǾ��پ�Q쾸�����H��푾�p���󶾿|��%���7��+������111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	�!B	�-B	�3B	��B	��B	�B	�#B	�)B	�/B	�;B	�B	�
B
B	��B	��B
VB
�B
-B
A�B
G�B
S�B
[#B
]/B
]/B
^5B
_;B
bNB
bNB
k�B
iyB
l�B
p�B
s�B
|�B
��B
q�B
ɺB
y�B
�B
�%B
�B
�1B
��B
�B
��B
��B
��B
�B
�B
�!B
�B
�9B
�9B
�qB
�XB
�RB
�'B
�^B
�qB
�}B
��B
ŢB
ȴB
��B
��B
��B
�B
�TB
�B
�B
�B+B
=B\BhB�B�B�B�B�B�B{BoB\BJBVBPBbBhBbB�B�B�B�B�B�B�BbBoB�B�BD�B/B'�B/B^5BS�Bv�B��B��B�-B�B�-B�'B�FB�^B��B�}B��B�ZB�`B�TB�ZB�;B��BƨB�9B�-B�RB�RB�3B�qBÖBɺBȴB��B��BȴBǮBĜB��B�}B�jB�qB�qB�}B�qB��B��B��B��B��B��BɺB��B�qB�jBɺB�RB�-B�?B�RB�RB�jB�XB�^B�^B�^B�^B��B�qB�jB�qB�dB�dB�^B�LB�LB�3B�3B�-B�'B�!B�!B�'B�!B�-B�?B�9B�9B�?B�?B�?B�FB�LB�LB�XB�LB�RB�^B�^B�dB�dB�dB�dB�dB�dB�jB�dB�^B�RB�9B�9B�-B�3B�-B�3B�3B�3B�?B�FB�?B�RB�-B�'B�!B�B��B��B��B��B��B�B�-B�'B�!B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B��B�B�B�B�B�!B�B�B�!B�!B�!B�!B�'B�3B�9B�9B�3B�3B�FB�FB�?B�3B�9B�3B�3B�'B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B	�!B	�-B	�3B	��B	��B	�B	�#B	�)B	�/B	�;B	�B	�
B
B	��B	��B
VB
�B
-B
A�B
G�B
S�B
[#B
]/B
]/B
^5B
_;B
bNB
bNB
k�B
iyB
l�B
p�B
s�B
|�G�O�G�O�G�O�G�O�B
�B
�%B
�B
�1B
��B
�B
��B
��B
��B
�B
�B
�!B
�B
�9B
�9B
�qB
�XB
�RB
�'B
�^B
�qB
�}B
��B
ŢB
ȴB
��B
��B
��B
�B
�TB
�B
�B
�B+B
=B\BhB�B�B�B�B�B�B{BoB\BJBVBPBbBhBbB�B�B�B�B�B�B�BbBoB�B�BD�B/B'�B/B^5BS�Bv�B��B��B�-B�B�-B�'B�FB�^B��B�}B��B�ZB�`B�TB�ZB�;B��BƨB�9B�-B�RB�RB�3B�qBÖBɺBȴB��B��BȴBǮBĜB��B�}B�jB�qB�qB�}B�qB��B��B��B��B��B��BɺB��B�qB�jBɺB�RB�-B�?B�RB�RB�jB�XB�^B�^B�^B�^B��B�qB�jB�qB�dB�dB�^B�LB�LB�3B�3B�-B�'B�!B�!B�'B�!B�-B�?B�9B�9B�?B�?B�?B�FB�LB�LB�XB�LB�RB�^B�^B�dB�dB�dB�dB�dB�dB�jB�dB�^B�RB�9B�9B�-B�3B�-B�3B�3B�3B�?B�FB�?B�RB�-B�'B�!B�B��B��B��B��B��B�B�-B�'B�!B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B��B�B�B�B�B�!B�B�B�!B�!B�!B�!B�'B�3B�9B�9B�3B�3B�FB�FB�?B�3B�9B�3B�3B�'B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111144441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202308050755342023080507553420230805075534  IF  ARFMCODA035h                                                                20200829092229                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200829092324  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200829092324  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200915134654  IP  PSAL            A1��D��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607172540  IP  PSAL            A1��D��3G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20230805075534  IP  PSAL            A1��D��3G�O�                