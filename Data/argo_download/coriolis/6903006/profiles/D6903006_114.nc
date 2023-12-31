CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2023-01-04T00:20:23Z creation; 2023-01-04T00:23:01Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_054b      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    :�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    :�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  ;   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  ;D   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  ;�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        ;�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     <   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     <8   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    <X   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         <\   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    <d   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            <h   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           <p   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           <x   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    <�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    <�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    <�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    <�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        =�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�        =�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  D�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  E�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  I4   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  J   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  M�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  R   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  U�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  V�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  Z    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  [   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ^�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  b$   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  c   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  f�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  g|   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    w�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    w�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    w�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    w�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  w�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    w�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    w�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    w�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         w�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         w�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        x    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    x   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  k   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    kL   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    oL   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    sL   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  wLArgo profile    3.1 1.2 19500101000000  20230104002023  20230131085301  6903006 NARVAL                                                          Camille DAUBORD                                                 MTIME           PRES            TEMP            PSAL               rA   IF                                  2C  D   ARVOR                           AI2600-19FR101                  5900A04                         844 @�Zq�r1   @�Zq�r@S*@
�Li����w�8   GPS         A   A   A   Primary sampling: averaged [10 sec sampling, 10 dbar average from 1000 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                     A.�~    A.�~    �8EȠ   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �r�m�@  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���l�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��"�9P  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����(  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���.�   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��ò��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���#H  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��l��  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����l  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���,   A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����ô  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ���m�8  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    �� �.F  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ��J�͐  A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    A.�~    ����  990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990999999999999990A!��A,��A<��AK33A^ffAnffA|��A�33A�ffA���A�  A���A���A�  A�  A�  A�33A�ffA���A�  A�33A�33B   B��B��B��B��B  B  BffB   B$  B(  B,ffB/33B4  B8  B<ffB@ffBD  BHffBS��Bg��B|��B���B�33B�33B���B�  B�33B�ffB�  B�  B���B�ffB�  C 33C�fC	��C�fC��C�fC�C"�fC'�3C-  C1�fC6�3C<  C@�fCE�3CK  CP  CU  CZ  C^�3Cc�fCi33Cn  Cr�fCw�3C|��C��C���C��C�� C��fC�s3C�  C�s3C�ٚC�� C��3C�Y�C�  C���C�  C�s3C��fC�ffC��fC�ffC��fC�Y�C�  C�� C�  C�s3C��fCČ�C�  C�s3C��fC�ffC��fC�s3C�  C؀ C��C݀ C��3C�ffC��fC�ffC��3C�ffC�  C� C��fC�Y�C�ٚC�ffC�  D 9�Ds3D� D��D,�Dl�D��D��D
,�Ds3D��DfD33DffD��D�3D,�D�fD� D��D9�D� D�fD��D33Dy�D � D"�D#9�D$ffD%��D'�D(@ D)y�D*�fD,fD-9�D.s3D/��D0��D2,�D3s3D4��D5��D7@ D8�fD9��D;�D<9�D=ffD>�3D?�3DA,�DB�fDC�fDEfDFFfDGl�DH�3DJ  DK33DLffDM�3DO  DPFfDQ�fDR��DS��DU&fDVffDW�fDX�fDZ33D[�fD\� D]��D_33D`l�Da��Db�3Dd33Des3Df��Dh  DiFfDjl�Dk�3Dl��Dn&fDos3Dp� Dq�3Ds33Dt��Du�3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A!��A,��A<��AK33A^ffAnffA|��A�33A�ffA���A�  A���A���A�  A�  A�  A�33A�ffA���A�  A�33A�33B   B��B��B��B��B  B  BffB   B$  B(  B,ffB/33B4  B8  B<ffB@ffBD  BHffBS��Bg��B|��B���B�33B�33B���B�  B�33B�ffB�  B�  B���B�ffB�  C 33C�fC	��C�fC��C�fC�C"�fC'�3C-  C1�fC6�3C<  C@�fCE�3CK  CP  CU  CZ  C^�3Cc�fCi33Cn  Cr�fCw�3C|��C��C���C��C�� C��fC�s3C�  C�s3C�ٚC�� C��3C�Y�C�  C���C�  C�s3C��fC�ffC��fC�ffC��fC�Y�C�  C�� C�  C�s3C��fCČ�C�  C�s3C��fC�ffC��fC�s3C�  C؀ C��C݀ C��3C�ffC��fC�ffC��3C�ffC�  C� C��fC�Y�C�ٚC�ffC�  D 9�Ds3D� D��D,�Dl�D��D��D
,�Ds3D��DfD33DffD��D�3D,�D�fD� D��D9�D� D�fD��D33Dy�D � D"�D#9�D$ffD%��D'�D(@ D)y�D*�fD,fD-9�D.s3D/��D0��D2,�D3s3D4��D5��D7@ D8�fD9��D;�D<9�D=ffD>�3D?�3DA,�DB�fDC�fDEfDFFfDGl�DH�3DJ  DK33DLffDM�3DO  DPFfDQ�fDR��DS��DU&fDVffDW�fDX�fDZ33D[�fD\� D]��D_33D`l�Da��Db�3Dd33Des3Df��Dh  DiFfDjl�Dk�3Dl��Dn&fDos3Dp� Dq�3Ds33Dt��Du�3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�����ٿ�r���b��r��ؓu��Kǿ�l����y��ff��ȴ��ȴ��ff��ȴ��ȴ��
=��E���33�щ7�Ұ!��S���bN��Q��9X������/���j���j�������Ͽ��𿻥㿾V���h��푿��D��1��dZ�����������������^���y����I����u����d���V�1&�<j�(�þ�����t�>��?PbN?~5??���?��`?�
=?�bN?�l�?Ұ!?�$�?�|�?�?��?�Z?��`?���?���?���?���?���?�S�@�w@+@��?��H?�"�?ش9?Ѓ?�/?�X?�+?���?�;d?�o?�I�?�"�?׮?�-?Ƈ+?��?�bN?���?��/?�l�?�?�  ?���?���?��?�t�?���?�v�?ؓu?�o?�V?��?���?�$�?��h?�33?�=q?��?���?�?�O�?���?���?�E�?�Z?~5??kƨ?^�R?q&�?l1?j��?]p�?U?S33?\j?`  ?{"�?���?|�?]�?O�?L1?T��?T�j?L�D?P�`?rn�?k?d�?L�D?:^5?0bN?2�!?WK�?_;d?[��?U�?X��?[dZ?@  ?49X?%�T?�?ȴ?hs?O�?	x�?��?�>�dZ?�7>��H>��j>�V>�>�>ݲ->�b>���>�33>�Z>��>�(�>��>�I�>���>�ƨ>���>r�!>]/>I�^>:^5>333>�->�u>\)>�+>�+>bN=�`B=ě�=�v�=��-=�hs=Y�=,1=�P<��<�<�9X<T��;�`B:�o:�o;o<t�    ��/�P�`�D���49X��w�49X�@��T���ixս#�
�t��#�
111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111���ٿ�r���b��r��ؓu��Kǿ�l����y��ff��ȴ��ȴ��ff��ȴ��ȴ��
=��E���33�щ7�Ұ!��S���bN��Q��9X������/���j���j�������Ͽ��𿻥㿾V���h��푿��D��1��dZ�����������������^���y����I����u����d���V�1&�<j�(�þ�����t�>��?PbN?~5??���?��`?�
=?�bN?�l�?Ұ!?�$�?�|�?�?��?�Z?��`?���?���?���?���?���?�S�@�w@+@��?��H?�"�?ش9?Ѓ?�/?�X?�+?���?�;d?�o?�I�?�"�?׮?�-?Ƈ+?��?�bN?���?��/?�l�?�?�  ?���?���?��?�t�?���?�v�?ؓu?�o?�V?��?���?�$�?��h?�33?�=q?��?���?�?�O�?���?���?�E�?�Z?~5??kƨ?^�R?q&�?l1?j��?]p�?U?S33?\j?`  ?{"�?���?|�?]�?O�?L1?T��?T�j?L�D?P�`?rn�?k?d�?L�D?:^5?0bN?2�!?WK�?_;d?[��?U�?X��?[dZ?@  ?49X?%�T?�?ȴ?hs?O�?	x�?��?�>�dZ?�7>��H>��j>�V>�>�>ݲ->�b>���>�33>�Z>��>�(�>��>�I�>���>�ƨ>���>r�!>]/>I�^>:^5>333>�->�u>\)>�+>�+>bN=�`B=ě�=�v�=��-=�hs=Y�=,1=�P<��<�<�9X<T��;�`B:�o:�o;o<t�    ��/�P�`�D���49X��w�49X�@��T���ixս#�
�t��#�
111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB	.B	0!B	.B	/B	0!B	2-B	1'B	2-B	49B	33B	1'B	33B	33B	1'B	2-B	33B	5?B	:^B	5?B	:^B	L�B	G�B	J�B	e`B	iyB	k�B	jB	m�B	n�B	p�B	q�B	s�B	t�B	u�B	t�B	u�B	u�B	v�B	v�B	u�B	u�B	w�B	y�B	{�B	�B	�uB	��B	�qB	ǮB	�HB	�yB	��B
-B
W
B
dZB
�JB
��B
�B
��BDB{B,B49BD�BW
BcTBm�Bn�Bo�Bp�Bt�Bz�Bz�B{�B�B�1B��B�3B�B��B�B}�Bu�Bn�Bu�B�7B�1B�1B�VB��B��B��B�oB�7B�7B�DB�VB��B��B��B�B�'B�RB�FB�?B�3B�3B�B�B��B�hB�PB�7B�%B�PB�uB��B��B�{B�DB�B�B�B~�B�B{�Bw�B|�B� B}�By�Bx�Bz�B� B�B�JB�oB�VB�B~�B� B�B�7B�%B�DB��B��B�{B�bB�PB�+B�1B��B��B��B��B��B��B��B��B�hB�hB�hB�hB�oB�oB�hB�uB�hB�uB�uB�uB�uB�{B�{B��B�{B�uB�hB�bB�bB�hB�bB�bB�bB�{B��B�{B�{B�oB�uB�uB�oB�uB�oB�{B��B��B�uB�oB�hB�hB�hB�bB�bB�hB�oB�oB�oB�oB�oB�oB�{B�{B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B	.B	0!B	.B	/B	0!B	2-B	1'B	2-B	49B	33B	1'B	33B	33B	1'B	2-B	33B	5?B	:^B	5?B	:^B	L�B	G�B	J�B	e`B	iyB	k�B	jB	m�B	n�B	p�B	q�B	s�B	t�B	u�B	t�B	u�B	u�B	v�B	v�B	u�B	u�B	w�B	y�B	{�B	�B	�uB	��B	�qB	ǮB	�HB	�yB	��B
-B
W
B
dZB
�JB
��B
�B
��BDB{B,B49BD�BW
BcTBm�Bn�Bo�Bp�Bt�Bz�Bz�B{�B�B�1B��B�3B�B��B�B}�Bu�Bn�Bu�B�7B�1B�1B�VB��B��B��B�oB�7B�7B�DB�VB��B��B��B�B�'B�RB�FB�?B�3B�3B�B�B��B�hB�PB�7B�%B�PB�uB��B��B�{B�DB�B�B�B~�B�B{�Bw�B|�B� B}�By�Bx�Bz�B� B�B�JB�oB�VB�B~�B� B�B�7B�%B�DB��B��B�{B�bB�PB�+B�1B��B��B��B��B��B��B��B��B�hB�hB�hB�hB�oB�oB�hB�uB�hB�uB�uB�uB�uB�{B�{B��B�{B�uB�hB�bB�bB�hB�bB�bB�bB�{B��B�{B�{B�oB�uB�uB�oB�uB�oB�{B��B��B�uB�oB�hB�hB�hB�bB�bB�hB�oB�oB�oB�oB�oB�oB�{B�{B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
MTIME           PRES            TEMP            PSAL            not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            not applicable                                                                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement necessary until cycle 289. ASD observed form cycle 290 to the end.Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V02 & ARGO2021V03 -                                                                                 20230131085301202301310853012023013108530120230131085301IF  ARFMCODA054b                                                                20230104002023                      G�O�G�O�G�O�                IF  ARGQCOQC6.0                                                                 20230104002301  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC6.0                                                                 20230104002301  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  2.0 CTD2021V02 & ARGO2021V03                                        20230131085301  IP  PSAL            A!��Du�3G�O�                