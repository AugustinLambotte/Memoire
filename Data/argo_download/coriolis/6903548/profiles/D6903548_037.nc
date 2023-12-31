CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  S   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-06-20T11:56:51Z creation; 2021-06-07T15:42:40Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_029d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        	L  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  D$   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        	L  Fx   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  O�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     	L  R   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  [d   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  d�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  g   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  pP   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  r�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  {�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  �<   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 T  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     	L  �0   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
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
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20190620115651  20210607154240  6903548 NorArgo                                                         Kjell Arne Mork                                                 PRES            TEMP            PSAL               %A   IF                                  2C  D   ARVOR                           AI2600-18EU002                  5900A04                         844 @،w�[1   @،w�C� @TXth���@,��Ss1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 2 dbar average from 500 dbar to 100 dbar; 10 sec sampling, 1 dbar average from 100 dbar to 4.5 dbar]                                                      @�33@�33@�33A��A  A!��A.ffA>ffAP  Ac33Aq��A�  A���A�33A���A�33A���A�  A���A�33A�ffA���A�  A�  A���A�  A�33A�ffB33B��B  BffB33BffB  B��B#33B'33B+��B0ffB4ffB8  B;33B@  BC��BG33BLffBPffBT  BX  B[��B_��Bb��Bg��BlffBp  Bs33Bw33B{��B�  B���B���B���B�  B���B���B�  B���B���B���B�33B�  B���B�33B�  B���B�33B�  B���B���B���B�33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B�ffB���B�ffB���B�  B���B♚B�33B�33B�33B�33B�ffB���B���C�CffC��CffC	33C� C�3CffC33CffC��CffC33CffC��CL�C!� C#�3C%ffC'� C)��C+ffC-� C/�3C1� C333C5ffC7��C9ffC;�C=ffC?��CAffCC�CEffCG� CIL�CK� CM�3CO� CQL�CS� CU��CW��CYffC[33C]  C_ffCa�3Cc��Ce��Cg� CiffCkffCmffCoffCqffCsffCuffCwffCyffC{ffC}ffCffC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C���C���C���C�ٚC��fC��3C�� C�� C���C���C��fC�� C�� C�� C��fC���C�� C�� C��3C��3C��3C�� C��3C�� C�� C�� C���C���C��fC��fC�� C�� C���C��fC��fC��3C�� C���C��fC��3C�� C���C��3C�� C���C��fC�� C�� C�ٚC��3C�� C���C��fC��3C���C��fC³3C�� Cę�CŦfC�� Cǌ�CȦfCɳ3C�� C���C̦fCͦfCγ3Cϳ3C�� Cљ�CҦfCӳ3CԳ3C�� C֙�CצfCس3C�� C�ٚC۳3C�� C���Cޙ�CߦfC�� CᙚC�fC�3C�� C�� C���C�ٚC�fC�3C�� C���C�ٚC��3CC�fC�� C�C�fC�� C���C��3C���C��3C���C�� C�� C�ٚD 33D� D��D��D9�D�fD�3D	  D
,�Ds3D��D  DL�D�fD�fDfDFfD��D�3D�DL�D�fD�fD��D&fD` D ��D"�D#FfD$�fD%�fD'�D(9�D)l�D*��D+�3D-&fD.y�D/��D0��D2@ D3� D4� D6  D7@ D8��D9�3D:�fD<33D=��D>��D@�DAL�DB�3DC�3DE3DFS3DG�3DH�3DJ3DKS3DL� DM�fDN��DP@ DQ��DR�3DT  DU&fDVy�DW�fDY3DZ@ D[s3D\� D^�D_33D`` Da�fDcfDdFfDe� Df��Dg��Di@ Dj�fDk��Dl��Dn&fDos3Dp� Dr�Ds9�Dtl�Du� Dv��Dx,�Dy� Dz�fD{��D}33D~��D� D�y�D�3D�� D�\�D�	�D���D�9�D�� D�s3D�fD���D�` D���D���D�6fD��fD�s3D�3D��3D�VfD���D���D�<�D��3D��fD��D��fD�\�D�fD���D�6fD��fD�vfD�fD���D�` D��fD���D�C3D�� D�|�D��D���D�Y�D���D���D�@ D��fD�� D��D���D�Y�D���D�� D�@ D�� D��3D�&fD�� D�\�D���D���D�<�D��3D�y�D�#3D�� D�\�D���D��fD�0 D���D��fD�&fD��fD�c3D�3D��3D�FfD��D�|�D� D���D�VfD���D��3D�@ D���D�|�D��D���D�\�D�  D��3D�C3D��3D�|�D�3D���D�` D�fD�� D�9�D�� D�|�D�fD��3D�c3D�3Dģ3D�C3D��fD�y�D��D��3D�VfD���Dɠ D�33D�ٚD�|�D�fD̰ D�\�D�fDΣ3D�@ D�� D�|�D�#3DѼ�D�VfD���Dӣ3D�@ D���D�y�D��Dּ�D�` D�  Dأ3D�FfD�� D�vfD�  Dۼ�D�VfD�  Dݠ D�<�D�� D߀ D�  D��3D�ffD���D� D�C3D��fD�|�D�#3D幚D�` D�	�D�3D�@ D���D�|�D��D�� D�c3D�fD쩚D�<�D��3D�y�D�#3D�ɚD�l�D�3D� D�6fD�ٚD�|�D�  D���D�S3D���D���D�6fD�ٚD��3D�#31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�33@�33@�33A��A  A!��A.ffA>ffAP  Ac33Aq��A�  A���A�33A���A�33A���A�  A���A�33A�ffA���A�  A�  A���A�  A�33A�ffB33B��B  BffB33BffB  B��B#33B'33B+��B0ffB4ffB8  B;33B@  BC��BG33BLffBPffBT  BX  B[��B_��Bb��Bg��BlffBp  Bs33Bw33B{��B�  B���B���B���B�  B���B���B�  B���B���B���B�33B�  B���B�33B�  B���B�33B�  B���B���B���B�33B�  B�  B���B���B���B���B���B���B���B���B���B���B���B���B�ffB���B�ffB���B�  B���B♚B�33B�33B�33B�33B�ffB���B���C�CffC��CffC	33C� C�3CffC33CffC��CffC33CffC��CL�C!� C#�3C%ffC'� C)��C+ffC-� C/�3C1� C333C5ffC7��C9ffC;�C=ffC?��CAffCC�CEffCG� CIL�CK� CM�3CO� CQL�CS� CU��CW��CYffC[33C]  C_ffCa�3Cc��Ce��Cg� CiffCkffCmffCoffCqffCsffCuffCwffCyffC{ffC}ffCffC��3C��3C��3C��3C��3C��3C��3C�� C�� C�� C���C���C���C�ٚC��fC��3C�� C�� C���C���C��fC�� C�� C�� C��fC���C�� C�� C��3C��3C��3C�� C��3C�� C�� C�� C���C���C��fC��fC�� C�� C���C��fC��fC��3C�� C���C��fC��3C�� C���C��3C�� C���C��fC�� C�� C�ٚC��3C�� C���C��fC��3C���C��fC³3C�� Cę�CŦfC�� Cǌ�CȦfCɳ3C�� C���C̦fCͦfCγ3Cϳ3C�� Cљ�CҦfCӳ3CԳ3C�� C֙�CצfCس3C�� C�ٚC۳3C�� C���Cޙ�CߦfC�� CᙚC�fC�3C�� C�� C���C�ٚC�fC�3C�� C���C�ٚC��3CC�fC�� C�C�fC�� C���C��3C���C��3C���C�� C�� C�ٚD 33D� D��D��D9�D�fD�3D	  D
,�Ds3D��D  DL�D�fD�fDfDFfD��D�3D�DL�D�fD�fD��D&fD` D ��D"�D#FfD$�fD%�fD'�D(9�D)l�D*��D+�3D-&fD.y�D/��D0��D2@ D3� D4� D6  D7@ D8��D9�3D:�fD<33D=��D>��D@�DAL�DB�3DC�3DE3DFS3DG�3DH�3DJ3DKS3DL� DM�fDN��DP@ DQ��DR�3DT  DU&fDVy�DW�fDY3DZ@ D[s3D\� D^�D_33D`` Da�fDcfDdFfDe� Df��Dg��Di@ Dj�fDk��Dl��Dn&fDos3Dp� Dr�Ds9�Dtl�Du� Dv��Dx,�Dy� Dz�fD{��D}33D~��D� D�y�D�3D�� D�\�D�	�D���D�9�D�� D�s3D�fD���D�` D���D���D�6fD��fD�s3D�3D��3D�VfD���D���D�<�D��3D��fD��D��fD�\�D�fD���D�6fD��fD�vfD�fD���D�` D��fD���D�C3D�� D�|�D��D���D�Y�D���D���D�@ D��fD�� D��D���D�Y�D���D�� D�@ D�� D��3D�&fD�� D�\�D���D���D�<�D��3D�y�D�#3D�� D�\�D���D��fD�0 D���D��fD�&fD��fD�c3D�3D��3D�FfD��D�|�D� D���D�VfD���D��3D�@ D���D�|�D��D���D�\�D�  D��3D�C3D��3D�|�D�3D���D�` D�fD�� D�9�D�� D�|�D�fD��3D�c3D�3Dģ3D�C3D��fD�y�D��D��3D�VfD���Dɠ D�33D�ٚD�|�D�fD̰ D�\�D�fDΣ3D�@ D�� D�|�D�#3DѼ�D�VfD���Dӣ3D�@ D���D�y�D��Dּ�D�` D�  Dأ3D�FfD�� D�vfD�  Dۼ�D�VfD�  Dݠ D�<�D�� D߀ D�  D��3D�ffD���D� D�C3D��fD�|�D�#3D幚D�` D�	�D�3D�@ D���D�|�D��D�� D�c3D�fD쩚D�<�D��3D�y�D�#3D�ɚD�l�D�3D� D�6fD�ٚD�|�D�  D���D�S3D���D���D�6fD�ٚD��3D�#31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@����7�G��G��   ��푾�푾�p���j��^5��X��Q��X���پ�E����m���#���ۿG���7��� ���X���>5?}?:^5?h��?�M�?���?�%?���?�X?�V?�bN?��/?���?�S�?�r�?��`?ؓu?�S�@X@
J@�@3��@=��@>�+@L�D@R�H@R�H@X  @[t�@[ƨ@Z�@\�@UO�@S@TI�@V�y@A�7@8��@NE�@J�!@C�m@>�y@81'@@r�@LZ@A&�@@1'@97L@4�@>E�@F�@P1'@P  @P1'@O�P@O�P@Tj@T��@R��@P �@M�@H �@F{@F@F5?@E�@E�-@E��@G�@G;d@F@E?}@Qx�@Q�^@<�/@(Ĝ@1�@2�!@.v�@$1@!G�@�@�H@hs@�`@r�@�@�@z�@!�7@!��@\)@ Q�@!��@"^5@ ��@!��@"~�@"-@'�w@*M�@*��@*=q@(�@%�-@%��@&v�@"��@!��@"M�@#��@$�D@$��@$�@#�
@$I�@$��@$�@$��@$��@%O�@%�-@%�@&ȴ@'�@'
=@'l�@'�;@(  @)&�@)�^@)��@*-@)�#@)�^@)�^@)��@)�7@(Ĝ@(��@(�`@( �@'��@'
=@&��@&V@%�@%O�@$�j@$�@$j@$Z@$1@#ƨ@#��@#C�@"�!@"^5@"n�@"=q@"�@"J@"J@"�@"�@"J@!��@!��@!%@ �u@  �@l�@�y@��@��@�+@@�-@O�@�@z�@33@^5@=q@��@G�@�9@r�@Q�@r�@�@Q�@ �@\)@�y@��@�+@�@V@�@@n�@�#@x�@J@M�@�@��@��@�^@��@�@�^@~�@~�@��@�\@M�@��@�7@G�@��@��@�^@��@G�@%@��@��@Ĝ@bN@1'@��@��@|�@�@��@�@��@?}@O�@?}@�@��@�h@�h@�h@O�@�D@z�@(�@9X@I�@I�@�D@�@��@I�@�@�
@�m@�F@33@o@
n�@	��@	�@�u@�;@K�@ȴ@�+@?}@��@9X@�@7L@   ?�ƨ?��+?�33?�  ?��?�"�?�P?�bN?�^5?�K�?�Z?��`?���?�9X?��;?�I�?���?���?�j?��#?�$�?�33?��!?�|�?�+?�o?�  ?�p�?�b?���?��7?{�m?r�!?k�?cS�?W
=?PbN?I7L?=�?:��?1&�?,�D?$�?dZ?9X?&�?
��?��>�v�>���>�h>�Ĝ>�(�>��>��>��j>�`B>�/>���>�I�>�$�>�%>n��>\(�>E��>'�>z�==ȴ9=�-=��=�\)=y�#=@�=�w<�`B<�C�;�`B    �e`B��9X��w�u�����hs���P��1������S���h���#�J�I���u�(�þ:^5�>vɾC���P�`�fff�l�D�m�h�o���p�׾y�#��  ���\������$ݾ�7L���;�\)��n������������㾜���5?��G����T��l����羭V���׾�-��Q쾺�H����o��$ݾ�=q��V���;��녾�����ۥ��(���;d�߾w��Ĝ���
��ff��~������33��Q���H��dZ���m�   ��7�o�Z��˿���y�1'��ÿ
=q�	��
=q�
����ƨ��Ϳ��������;�bN��`����n��t��9X��j����?}��+��+�ȴ��P��u�X�����#��#��������/��-�5?�vɿ;d��w� Ĝ�!G��!���"J�"��#S��#���#���$Z�%��%�T�&ff�&��&��'l��'(r��(r��(�ÿ)�^�)�^�*���+C��+��+ƨ�,1�,�D�-V�-O߿-��-��-��.��/��/�;�0 ſ0�`�1&�1���2-�2-�2�!�333�3t��3t��49X�4�j�4���5?}�5��5?}�5��5�5�5��6E��7
=�7Kǿ7Kǿ7Kǿ8b�8b�8�u�8�u�8�u�8�u�9��9�#�:^5�;"ѿ;�m�;�m�<(��<��<푿=p��=�-�=�-�=�>vɿ>�R�>�ۿ>�ۿ?;d�?;d�?|�?�w�@  �@  �@  �@A��@��@��@Ĝ�A%�A%�AG��AG�1111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ��7�G��G��   ��푾�푾�p���j��^5��X��Q��X���پ�E����m���#���ۿG���7��� ���X���>5?}?:^5?h��?�M�?���?�%?���?�X?�V?�bN?��/?���?�S�?�r�?��`?ؓu?�S�@X@
J@�@3��@=��@>�+@L�D@R�H@R�H@X  @[t�@[ƨ@Z�@\�@UO�@S@TI�@V�y@A�7@8��@NE�@J�!@C�m@>�yG�O�G�O�@LZ@A&�@@1'@97L@4�@>E�@F�@P1'@P  @P1'@O�P@O�P@Tj@T��@R��@P �@M�@H �@F{@F@F5?@E�@E�-@E��@G�@G;d@F@E?}@Qx�@Q�^@<�/@(Ĝ@1�@2�!@.v�@$1@!G�@�@�H@hs@�`@r�@�@�@z�@!�7@!��@\)@ Q�@!��@"^5@ ��@!��@"~�@"-@'�w@*M�@*��@*=q@(�@%�-@%��@&v�@"��@!��@"M�@#��@$�D@$��@$�@#�
@$I�@$��@$�@$��@$��@%O�@%�-@%�@&ȴ@'�@'
=@'l�@'�;@(  @)&�@)�^@)��@*-@)�#@)�^@)�^@)��@)�7@(Ĝ@(��@(�`@( �@'��@'
=@&��@&V@%�@%O�@$�j@$�@$j@$Z@$1@#ƨ@#��@#C�@"�!@"^5@"n�@"=q@"�@"J@"J@"�@"�@"J@!��@!��@!%@ �u@  �@l�@�y@��@��@�+@@�-@O�@�@z�@33@^5@=q@��@G�@�9@r�@Q�@r�@�@Q�@ �@\)@�y@��@�+@�@V@�@@n�@�#@x�@J@M�@�@��@��@�^@��@�@�^@~�@~�@��@�\@M�@��@�7@G�@��@��@�^@��@G�@%@��@��@Ĝ@bN@1'@��@��@|�@�@��@�@��@?}@O�@?}@�@��@�h@�h@�h@O�@�D@z�@(�@9X@I�@I�@�D@�@��@I�@�@�
@�m@�F@33@o@
n�@	��@	�@�u@�;@K�@ȴ@�+@?}@��@9X@�@7L@   ?�ƨ?��+?�33?�  ?��?�"�?�P?�bN?�^5?�K�?�Z?��`?���?�9X?��;?�I�?���?���?�j?��#?�$�?�33?��!?�|�?�+?�o?�  ?�p�?�b?���?��7?{�m?r�!?k�?cS�?W
=?PbN?I7L?=�?:��?1&�?,�D?$�?dZ?9X?&�?
��?��>�v�>���>�h>�Ĝ>�(�>��>��>��j>�`B>�/>���>�I�>�$�>�%>n��>\(�>E��>'�>z�==ȴ9=�-=��=�\)=y�#=@�=�w<�`B<�C�;�`B    �e`B��9X��w�u�����hs���P��1������S���h���#�J�I���u�(�þ:^5�>vɾC���P�`�fff�l�D�m�h�o���p�׾y�#��  ���\������$ݾ�7L���;�\)��n������������㾜���5?��G����T��l����羭V���׾�-��Q쾺�H����o��$ݾ�=q��V���;��녾�����ۥ��(���;d�߾w��Ĝ���
��ff��~������33��Q���H��dZ���m�   ��7�o�Z��˿���y�1'��ÿ
=q�	��
=q�
����ƨ��Ϳ��������;�bN��`����n��t��9X��j����?}��+��+�ȴ��P��u�X�����#��#��������/��-�5?�vɿ;d��w� Ĝ�!G��!���"J�"��#S��#���#���$Z�%��%�T�&ff�&��&��'l��'(r��(r��(�ÿ)�^�)�^�*���+C��+��+ƨ�,1�,�D�-V�-O߿-��-��-��.��/��/�;�0 ſ0�`�1&�1���2-�2-�2�!�333�3t��3t��49X�4�j�4���5?}�5��5?}�5��5�5�5��6E��7
=�7Kǿ7Kǿ7Kǿ8b�8b�8�u�8�u�8�u�8�u�9��9�#�:^5�;"ѿ;�m�;�m�<(��<��<푿=p��=�-�=�-�=�>vɿ>�R�>�ۿ>�ۿ?;d�?;d�?|�?�w�@  �@  �@  �@A��@��@��@Ĝ�A%�A%�AG��AG�1111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��BBBBBBBÖBÖBŢBƨBĜBÖB��BĜBƨBBBÖBĜBƨB�BB  BdZB�qB�BB�B	{B	0!B	49B	;dB	@�B	F�B	M�B	ZB	jB	w�B	�+B	�B	ÖB	�#B	��B
uB
{�B
�!B
�!B
��B
�BoB%�B7LB;dB@�B=qBG�BD�BE�BN�BC�B>wBt�Bm�BjBffB[#B�B�VB�B~�Bs�Bk�B�DB�bB�B�'B�'B�9B�!B��B�}B�dB�XB�LB�dB�'B�'B�3B�'B�3B�-B�-B�LB�^B�LB��BȴB��B�{B��B��B��B�bB�\B�+B}�B�B�%B�1B�DB�\B�JB��B��B��B��B�B�B�B�B�!B�9B��BBĜBŢBB�qB�jB��B�qB�dB�qB��BŢBƨBǮBȴB��B��B��B��B��B��B��B��B��B�
B�B�B�)B�5B�HB�TB�TB�ZB�ZB�TB�ZB�ZB�`B�fB�fB�fB�`B�`B�`B�ZB�ZB�TB�NB�NB�NB�NB�NB�NB�NB�NB�TB�HB�BB�HB�HB�NB�HB�NB�NB�NB�TB�TB�NB�NB�HB�BB�BB�;B�;B�;B�5B�5B�5B�/B�/B�5B�#B�B�B�B�B�B�B�
B�B�B�B�B�B�B��B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBƨBŢBƨBĜBĜBB��B�}B�}B�dB�jB�^B�LB�LB�9B�9B�3B�3B�-B�?B�9B�9B�9B�3B�-B�3B�-B�'B�!B�!B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�1111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 B�B�B�B�B�B�B�B�"B�"B�.B�4B�(B�"B�B�(B�4B�B�B�"B�(B�4B��B�Bg�B��B��B�BB	B	3�B	7�B	>�B	DB	J4B	Q_B	]�B	nB	{[B	��B	��B	�"B	ޯB	�HB
B
sB
��B
��B
ׄB
�B�B)oB:�B>�BDB@�BK:BH(BI.BReBG"BBBxHBqBnBi�G�O�G�O�B��B��B��BwBBoB��B��B��B��B��B��B��B�B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�_B�@B�uB�B�cB�WB�QB��B��B��B��B��B��B��B��B��B��B�DB�JB�QB��B��B��B��B��B��B��B�B�B�(B�.B�B��B��B�B��B��B��B�B�.B�4B�:B�@B�MB�YB�YB�YB�_B�eB�qB�xB؊BږBۜBܣBߵB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��BޯBܣBܣBݩBܣBܣBܣBږBۜBܣBۜBۜBِBِB؊BږBׄB�~B�~B�xB�kB�qB�eB�qB�xB�qB�xB�qB�qB�qB�xB�qB�~B�~B�~B�~B�xB�xB�xB�qB�qB�qB�qB�xB�xB�qB�qB�qB�qB�kB�qB�kB�kB�eB�kB�eB�YB�YB�YB�YB�_B�YB�_B�_B�_B�_B�_B�YB�_B�_B�_B�_B�eB�kB�kB�kB�kB�kB�qB�kB�kB�kB�kB�kB�eB�eB�_B�_B�_B�_B�_B�YB�YB�YB�SB�MB�FB�FB�4B�.B�4B�(B�(B�B�B�	B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�|B�uB�uB�oB�iB�iB�cB�cB�cB�]B�]B�]B�]B�WB�WB�WB�QB�JB�JB�DB�DB�DB�DB�DB�DB�DB�>B�>B�8B�8B�2B�8B�8B�2B�2B�8B�8B�8B�8B�8B�8B�8B�8B�2B�8B�8B�>B�>B�8B�>B�>B�>B�>B�>B�DB�DB�DB�DB�DB�DB�DB�JB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�QB�QB�QB�QB�WB�QB�QB�WB�WB�QB�WB�WB�WB�WB�WB�WB�WB�WB�WB�]B�WB�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�]B�cB�cB�cB�cB�cB�cB�cB�iB�iB�cB�iB�cB�cB�iB�iB�iB�iB�iB�iB�iB�iB�iB�oB�iB�oB�iB�oB�oB�oB�iB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�oB�oB�uB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B�|B��B��B�|B��B�|B��B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                                                none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1.0001, vertically averaged dS= 0.0034623                                                                                                                                                                                                                    No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Significant salinity drift present  - correction applied using OW method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                      202106071542402021060715424020210607154240  IF  ARFMCODA029d                                                                20190620115651                      G�O�G�O�G�O�                IF  ARGQCOQC4.2                                                                 20190620115722  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.2                                                                 20190620115722  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210607154240  IP  PSAL            @�33D�#3G�O�                