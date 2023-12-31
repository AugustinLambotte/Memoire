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
_FillValue                 �  BL   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  DD   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  L    PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  N   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  qx   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  yT   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  {L   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �    HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �X   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �\   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �`   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �d   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �h   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �,   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �,   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �,   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �,Argo profile    3.1 1.2 19500101000000  20200828144324  20220127170406  3902102 ARGO POLAND                                                     Waldemar Walczowski                                             PRES            TEMP            PSAL               CA   IF                                  2C  D   ARVOR                           AI2600-17EU026                  5900A04                         844 @���q�1   @���q�@R�5(�|�&ɰb�k8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 2.5 dbar]                                                      @�33@�33@�33A��A  A��A0  A@  AQ��Aa��Aq��A�  A�33A�33A�  A���A���A�  A���A�33A�ffA���A���A�33A噚A�ffA�  A�ffB��B��B  B��B33B��B��B ��B#��B(ffB+��B/��B4ffB8ffB<ffB@ffBC33BG��BL  BP��BT  BX  B[��B_��Bc33Bg��Bk��Bp  Bs��Bw33B|  B��B���B�  B���B���B�  B�ffB���B���B�  B���B�  B�33B�  B���B���B���B�  B�33B���B���B���B�ffB�ffB�ffB���B���B���B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B���Bҙ�B֙�Bڙ�B�  B♚B晚B���BB�ffB���B�  B���C� C��C�3CffC	33C� C� C� C� C� C��CL�C� C�3C� CL�C!33C#ffC%��C'�3C)� C+ffC-ffC/ffC1L�C3L�C533C733C933C;�C=�C?�CA�CC�CE�CG�CI�CK�CM  CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc33Ce33Cg33Ci33CkL�CmL�CoffCqffCsffCuffCwffCyffC{ffC}L�CL�C���C���C��3C��fC�ٚC���C���C�� C��3C��3C��fC���C���C���C���C���C���C��fC��fC��3C��3C��3C�� C�� C�� C�� C��3C��3C��3C��fC��fC��fC��fC��fC��fC��3C�� C���C��3C�� C���C��fC�� C��fC�� C�ٚC�� C��fC���C��3C�ٚC�� C��fC���C��3C��fC�ٚC���C���C�� C��3C��3C��3C��3C��3C��3C³3Có3C�� C�� C�� C���C���C�ٚC��fC��3C�� C̀ CΌ�Cό�CЙ�CѦfCҦfCӳ3CԳ3C�� C�� C���C���C�ٚC��fC��3C�� C݀ Cތ�Cߙ�C�fC�3C�3C�� C���C�ٚC��fC�3C� C��CꙚC�fC�3C�� C���C�ٚC��3C�� C��C��C���C��fC��fC��3C��3C�� C�� C��D 9�DffD�3D  DL�Ds3D� D	  D
L�D�fD��D�3D@ D��D��D3D9�DffD��D��D@ D�fD��D��DFfD� D � D"  D#,�D$y�D%�fD&��D(33D)l�D*�fD+�fD-9�D.l�D/��D1fD2@ D3�fD4� D5��D7@ D8� D9�fD:�3D<9�D=� D>��D?�3DA&fDBffDC�fDD��DF33DGy�DH��DJ  DK33DLl�DM�fDN� DP9�DQy�DR�fDS��DU33DV�fDW�3DX��DZ@ D[� D\��D^  D_@ D`�fDa��Db�3Dd33De� Df�fDg��Di33Dj� Dk� Dm�Dn9�DoffDp��Dr�Ds@ Dty�Du��Dv��Dx33Dys3Dz��D{��D}FfD~��D�3D�� D�  D��fD�` D�  D��3D�C3D��3D��fD��D��3D�Y�D�  D��fD�@ D�ٚD�vfD�3D��3D�c3D�3D��3D�6fD���D��3D��D���D�VfD���D��fD�C3D���D�|�D��D���D�Y�D���D�� D�C3D��fD���D��D�� D�VfD��fD���D�C3D�� D�y�D�3D���D�ffD�	�D���D�33D��fD�y�D�#3D�� D�` D�  D�� D�@ D��3D��3D�fD���D�\�D��3D���D�C3D���D�vfD��D��fD�c3D���D���D�33D��3D�|�D�)�D��fD�c3D�  D���D�9�D�ٚD�s3D� D�� D�\�D��D���D�FfD�� D�|�D�fD���D�c3D�fD��f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�33@�33@�33A��A  A��A0  A@  AQ��Aa��Aq��A�  A�33A�33A�  A���A���A�  A���A�33A�ffA���A���A�33A噚A�ffA�  A�ffB��B��B  B��B33B��B��B ��B#��B(ffB+��B/��B4ffB8ffB<ffB@ffBC33BG��BL  BP��BT  BX  B[��B_��Bc33Bg��Bk��Bp  Bs��Bw33B|  B��B���B�  B���B���B�  B�ffB���B���B�  B���B�  B�33B�  B���B���B���B�  B�33B���B���B���B�ffB�ffB�ffB���B���B���B�  B�  B�  B�  B�  B�  B�33B�33B�  B���B���Bҙ�B֙�Bڙ�B�  B♚B晚B���BB�ffB���B�  B���C� C��C�3CffC	33C� C� C� C� C� C��CL�C� C�3C� CL�C!33C#ffC%��C'�3C)� C+ffC-ffC/ffC1L�C3L�C533C733C933C;�C=�C?�CA�CC�CE�CG�CI�CK�CM  CO�CQ�CS�CU�CW�CY�C[�C]�C_�Ca�Cc33Ce33Cg33Ci33CkL�CmL�CoffCqffCsffCuffCwffCyffC{ffC}L�CL�C���C���C��3C��fC�ٚC���C���C�� C��3C��3C��fC���C���C���C���C���C���C��fC��fC��3C��3C��3C�� C�� C�� C�� C��3C��3C��3C��fC��fC��fC��fC��fC��fC��3C�� C���C��3C�� C���C��fC�� C��fC�� C�ٚC�� C��fC���C��3C�ٚC�� C��fC���C��3C��fC�ٚC���C���C�� C��3C��3C��3C��3C��3C��3C³3Có3C�� C�� C�� C���C���C�ٚC��fC��3C�� C̀ CΌ�Cό�CЙ�CѦfCҦfCӳ3CԳ3C�� C�� C���C���C�ٚC��fC��3C�� C݀ Cތ�Cߙ�C�fC�3C�3C�� C���C�ٚC��fC�3C� C��CꙚC�fC�3C�� C���C�ٚC��3C�� C��C��C���C��fC��fC��3C��3C�� C�� C��D 9�DffD�3D  DL�Ds3D� D	  D
L�D�fD��D�3D@ D��D��D3D9�DffD��D��D@ D�fD��D��DFfD� D � D"  D#,�D$y�D%�fD&��D(33D)l�D*�fD+�fD-9�D.l�D/��D1fD2@ D3�fD4� D5��D7@ D8� D9�fD:�3D<9�D=� D>��D?�3DA&fDBffDC�fDD��DF33DGy�DH��DJ  DK33DLl�DM�fDN� DP9�DQy�DR�fDS��DU33DV�fDW�3DX��DZ@ D[� D\��D^  D_@ D`�fDa��Db�3Dd33De� Df�fDg��Di33Dj� Dk� Dm�Dn9�DoffDp��Dr�Ds@ Dty�Du��Dv��Dx33Dys3Dz��D{��D}FfD~��D�3D�� D�  D��fD�` D�  D��3D�C3D��3D��fD��D��3D�Y�D�  D��fD�@ D�ٚD�vfD�3D��3D�c3D�3D��3D�6fD���D��3D��D���D�VfD���D��fD�C3D���D�|�D��D���D�Y�D���D�� D�C3D��fD���D��D�� D�VfD��fD���D�C3D�� D�y�D�3D���D�ffD�	�D���D�33D��fD�y�D�#3D�� D�` D�  D�� D�@ D��3D��3D�fD���D�\�D��3D���D�C3D���D�vfD��D��fD�c3D���D���D�33D��3D�|�D�)�D��fD�c3D�  D���D�9�D�ٚD�s3D� D�� D�\�D��D���D�FfD�� D�|�D�fD���D�c3D�fD��f11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����o��33��33��o�����\��n���Mӿ�-��J��녿�G���7������;��|��;d���ۿ����O߿�/��/��/��V��/��V��V��/��V��V��/��V��V��O߿ݑh�ݲ-�ݲ-��5?��5?��V��5?��V��vɿ�V��V�ޗ��޸R�޸R�޸R�ޗ��޸R����޸R�޸R�޸R�޸R��vɿ�V��{��{����{��V��V��5?��V��5?��5?��5?��{��{��{���ݑh��p���O߿��Ϳ�ƨ���H�ش9��Kǿ�$ݿ���?}�Լj��S����ѩ����`��Ĝ�ϝ����ۿ����ۿȴ9��xտ�X��1'�ǍP��+�Ƈ+�öF�°!��bN��O߿�J�T���5��D��A����\�^5?�=p��8Q�+>6E�>�K�>���?`B?!%?0bN?R�?u?}?�J?��?�;d?� �?�hs?�Ĝ?�A�?� �?�\)?�{?�{?��^?���?�"�?�"�?�?�J?�-?��?�33?��`?�  ?��w?�  ?�Ĝ?��;?�-?��h?���?�O�?�S�?�33?��?�?���?��?�O�?�V?��?���?���?��?��D?��?���?��^?�~�?���?���?�x�?�b?�ff?�?�?}?�S�?�bN?��?�j?�X?���?��`?��?�1?��#?�ff?;d?u?v?s�F?s33?st�?st�?r�?r-?r-?q&�?pbN?o�?n��?mO�?m��?mV?l�D?kƨ?k�?j~�?j=q?p�`?r-?l��?j~�?i�^?j��?u?y�?{�m?}/?�A�?�bN?��7?�hs?~5??wK�?p��?p �?lI�?l1?k�?j��?i��?g�?b��?b��?a��?a�7?bM�?bJ?a��?e`B?d�?c�
?bJ?a��?`  ?\�?X��?S��?Q��?Q&�?O�;?N{?NV?NV?N{?<�?%`B?dZ?z�?�;?{? �?,��?9�?A%?F�y?BM�?9�?7K�?/�;?1�?4z�?5?}?.��?+C�?.{?.��?.��?/�?/�;?0bN?0�`?1&�?1&�?1&�?1&�?)7L?"M�?|�?!%?"J?   ?�m?��?��?�+?9X?z�?�+?�h?
~�?�T?��?��?G�>�|�>��#>�!>�&�>��>�x�>�ff>ۥ�>�bN>ȴ9>��7>�ȴ>�1>�ff>��>��P>�
=>��+>��9>{�m>|�>�J>y�#>^5?>N�>;dZ>49X>=p�>@�>\)>�=��=�
==��`=�l�=��`=�\)=aG�=�+=�1=y�#=,1=��=o<�`B=Y�=0 �=�w<�=+=�P<�/<��
<u;�`B:�o    ��o��o�T�����
���ͽ49X��o���������T�����
=��;d����	7L�hs�z�t������!���''"��&�y�0 ž8Q�8Q�C���O�;�V�Xb�`A��fff�k��l�D�r�!�z�H�~�۾�J���˾��^�������;�hs��t���zᾖ���u�����������-��Z���y��l����羬1���h�������������׾��F���j��?}����ȴ��KǾ�X��푾�p���vɾ�  �Õ��š˾Ƨ��7L������C���ƨ���;�����`��n����ϾՁ��b�ؓu����������ڟ��ڟ��ܬ�ݲ-�ݲ-��5?�߾w�����Z���y����þ�~���~�����h�� ž�&��!��F��33��F��ȴ���پ��پ�X��� ���7�����\�`B�r��
=q�\)��ϿX��m� ��&$ݿ'+�'+�(�9�)�^�.��0�׿3��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��o��33��33��o�����\��n���Mӿ�-��J��녿�G���7������;��|��;d���ۿ����O߿�/��/��/��V��/��V��V��/��V��V��/��V��V��O߿ݑh�ݲ-�ݲ-��5?��5?��V��5?��V��vɿ�V��V�ޗ��޸R�޸R�޸R�ޗ��޸R����޸R�޸R�޸R�޸R��vɿ�V��{��{����{��V��V��5?��V��5?��5?��5?��{��{��{���ݑh��p���O߿��Ϳ�ƨ���H�ش9��Kǿ�$ݿ���?}�Լj��S����ѩ����`��Ĝ�ϝ����ۿ����ۿȴ9��xտ�X��1'�ǍP��+�Ƈ+�öF�°!��bN��O߿�J�T���5��D��A����\�^5?�=p��8Q�+>6E�>�K�>���?`B?!%?0bN?R�?u?}?�J?��?�;d?� �?�hs?�Ĝ?�A�?� �?�\)?�{?�{?��^?���?�"�?�"�?�?�J?�-?��?�33?��`?�  ?��w?�  ?�Ĝ?��;?�-?��h?���?�O�?�S�?�33?��?�?���?��?�O�?�V?��?���?���?��?��D?��?���?��^?�~�?���?���?�x�?�b?�ff?�?�?}?�S�?�bN?��?�j?�X?���?��`?��?�1?��#?�ff?;d?u?v?s�F?s33?st�?st�?r�?r-?r-?q&�?pbN?o�?n��?mO�?m��?mV?l�D?kƨ?k�?j~�?j=q?p�`?r-?l��?j~�?i�^?j��?u?y�?{�m?}/?�A�?�bN?��7?�hs?~5??wK�?p��?p �?lI�?l1?k�?j��?i��?g�?b��?b��?a��?a�7?bM�?bJ?a��?e`B?d�?c�
?bJ?a��?`  ?\�?X��?S��?Q��?Q&�?O�;?N{?NV?NV?N{?<�?%`B?dZ?z�?�;?{? �?,��?9�?A%?F�y?BM�?9�?7K�?/�;?1�?4z�?5?}?.��?+C�?.{?.��?.��?/�?/�;?0bN?0�`?1&�?1&�?1&�?1&�?)7L?"M�?|�?!%?"J?   ?�m?��?��?�+?9X?z�?�+?�h?
~�?�T?��?��?G�>�|�>��#>�!>�&�>��>�x�>�ff>ۥ�>�bN>ȴ9>��7>�ȴ>�1>�ff>��>��P>�
=>��+>��9>{�m>|�>�J>y�#>^5?>N�>;dZ>49X>=p�>@�>\)>�=��=�
==��`=�l�=��`=�\)=aG�=�+=�1=y�#=,1=��=o<�`B=Y�=0 �=�w<�=+=�P<�/<��
<u;�`B:�o    ��o��o�T�����
���ͽ49X��o���������T�����
=��;d����	7L�hs�z�t������!���''"��&�y�0 ž8Q�8Q�C���O�;�V�Xb�`A��fff�k��l�D�r�!�z�H�~�۾�J���˾��^�������;�hs��t���zᾖ���u�����������-��Z���y��l����羬1���h�������������׾��F���j��?}����ȴ��KǾ�X��푾�p���vɾ�  �Õ��š˾Ƨ��7L������C���ƨ���;�����`��n����ϾՁ��b�ؓu����������ڟ��ڟ��ܬ�ݲ-�ݲ-��5?�߾w�����Z���y����þ�~���~�����h�� ž�&��!��F��33��F��ȴ���پ��پ�X��� ���7�����\�`B�r��
=q�\)��ϿX��m� ��&$ݿ'+�'+�(�9�)�^�.��0�׿3��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�dB�qB�XB��B�BoB:^BJ�BS�BZB`BB�B��B�sB�B&�B7LBN�Bl�B}�B�B�JB�oB��B��B��B��B��B��B�B�-B��BBƨBɺB��B��B��B��B��B��B��B��B��B��B�
B�
B�B�B�#B�5B�BB�HB�TB�fB�B�B��B��B	  B	B	\B	�B	�B	�B	&�B	'�B	'�B	)�B	+B	-B	1'B	33B	=qB	A�B	G�B	P�B	bNB	iyB	w�B	�B	�B	�B	�B	�B	�1B	�7B	�DB	�VB	�bB	�oB	�oB	�oB	��B	��B	��B	��B	��B	��B	��B	��B	��B	��B	�B	ȴB	�B
�B
5?B
G�B
N�B
l�B
|�B
�B
�7B
�{B
��B
�5B
��BB�B'�B8RB@�BK�B]/B`BBaHBcTBe`BdZBcTBdZBe`BgmBgmBhsBjBl�Bm�Bk�Bm�Bn�Bq�Bq�Bs�Bt�Bu�Bx�By�Bw�B�\B�hB�hB�PB�JB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B�uB�\B�\B�VB�VB�VB�PB�VB�VB�bB�\B�VB�VB�\B�VB�PB�PB�VB�VB�VB�PB�PB�JB�JB�{B�oB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�hB�PB�PB�\B�VB�\B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�B�B�B�{B�BBB?BOnBX�B^�Bd�B��B�WB�&B5B+�B<BS�BqEB��B��B�B�*B�NB�B��B��B��B��B��B��B�>B�MB�eB�wB�}BъBׯB֪BصB׮BصBٻBصBصB��B��B��B��B��B��B��B�B�B�!B�BB�oB��B	�B	�B	�B	B	YB	^B	#wB	+�B	,�B	,�B	.�B	/�B	1�B	5�B	7�B	B2B	FJB	LoB	U�B	gB	n;B	|�B	��B	��B	��B	��B	��B	��B	��B	�B	�B	�'B	�3B	�3B	�2B	�EB	�jB	�vB	�^B	�xB	�{B	�{B	��B	��B	��B	��B	�yB	�EB
 pB
:B
LxB
S�B
q[B
��B
��B
�B
�IB
�QB
�B
��B�BbB,�B=*BE\BP�BbBeBf Bh*Bj8Bi3Bh,Bi4Bj7BlFBlDBmJBoWBqdBrjBp^BrhBsqBv�Bv�Bx�By�Bz�B}�B~�B|�B�6B�CB�BB�+B�$B�aB�pB�jB�sB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�vB��B�xB�UB�\B�PB�7B�6B�.B�0B�.B�*B�0B�1B�;B�6B�/B�2B�4B�1B�(B�+B�/B�/B�1B�*B�*B�%B�$B�XB�JB�NB�OB�UB�nB�rB��B��B��B��B��B��B��B�~B�|B��B�{B�uB�sB�~B�tB�kB�tB�iB�tB�{B�~B��B�{B��B��B��B��B��B��B��B��B��B��B��B�~B�B�B��B�yB��B�?B�(B�)B�6B�0B�2B�mB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B��B��B��B��B��B�B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001 (+/- 1e-05) , vertically averaged dS =0.004711 (+/- 0.01)                                                                                                                                                                                             No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No adjustment was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                               Salinity drift or offset detected - OW fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OW Method, 1.1,  -CTD2021V02 & ARGO2021V03 -                                                                                                            202011171219032022012717040620220127170406  IF  ARFMCODA035h                                                                20200828144324                      G�O�G�O�G�O�                IF  ARGQCOQC4.6                                                                 20200828144426  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.6                                                                 20200828144426  QCF$                G�O�G�O�G�O�0000000000000000PL  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology 20201117121903  IP  PSAL            @�33D��fG�O�                PL  ARSQOW  1.1 CTD2021V02 & ARGO2021V03                                        20220127170406  IP  PSAL            @�33D��fG�O�                