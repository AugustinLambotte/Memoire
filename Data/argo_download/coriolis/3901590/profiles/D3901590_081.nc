CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   P   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2019-01-16T13:02:03Z creation; 2019-01-16T13:02:37Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_027d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7@   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7P   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7T   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7X   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7h   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7x   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8@   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8D   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8H   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8L   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8l   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8p   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8t   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           8�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            8�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9    PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        @  :   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  ;X   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        @  ;�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  <�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     @  =8   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  >x   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  ?�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  @   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  AH   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  A�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  B�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  D   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  Dh   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  P  E�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     @  E�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    P�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    P�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    P�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    P�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  P�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    P�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    P�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    P�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Q   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Q   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Q   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Q   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  G8   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Gh   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Jh   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Mh   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  PhArgo profile    3.1 1.2 19500101000000  20190116130203  20210522160233  3901590 ARGO-BSH                                                        Birgit KLEIN                                                    PRES            TEMP            PSAL               QA   IF                                  2C  D   APEX                            7615                            102015                          846 @؍m�W��1   @؍n+�`@RO|�hs��/���1   ARGOS   A   A   F   Primary sampling: discrete []                                                                                                                                                                                                                                      @�33@�ffAffAy��A���A�  A噚B��B ffB4ffBJ  Br��B�ffB���B�  B�  B�33B�33CffCL�CffC�fC)��C4L�C>  CHL�CQ�3C\�Ce�fCo�3Cy��C��C�33C���C��3C�L�C��3C�@ C�&fC�  C��3C�� C�L�Cԙ�C�33C�ffC�&fDFfD	S3D�fD3DFfD"�fD(�fD.� D;� DH  DTffDa&fDm�fDz�D�<�D�6fD���D���D�I�D���D�ɚD�3D�@ D�y�D��3D��D�0 DԆfDڼ�D���D�FfD�|�D�11111111111111111111111111111111111111111111111111111111111111111111111111111111@���@�  A#33A~fgA�33A�ffA�  B��B!��B5��BK33Bt  B�  B�34B���BǙ�B���B���C�3C��C�3C 33C*�C4��C>L�CH��CR  C\fgCf33Cp  Cz�C�33C�Y�C��3C��C�s3C��C�ffC�L�C�&fC��C��fC�s3C�� C�Y�C��C�L�DY�D	ffDٙD&fDY�D"��D(ٙD.�3D;�3DH3DTy�Da9�Dm��Dz  D�FgD�@ D��gD�4D�S4D��4D��4D��D�I�D��4D���D�gD�9�DԐ D��gD�4D�P D�gD�411111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@5?@5?@5?@5?@5?@5?@$�@5?@V@ff@�+@�@��?��R?�1?��w?���?�o?{��?u?o\)?k�?b��?T�j?St�?Rn�?KC�?I7L?2�!?(1'?�H?�?I�?�? Ĝ>�ff>��
>��>�J>�  >��h>��;>��>Kƨ>'�=��=���=@�<�h�o�+��%���ͽ���z�G��l�D��J��\)������徬�D��KǾǮ��hs�߾w������H����	7L����!��+�����ۿ"J�%�˿*=q�,�D�/\)11111111111111111111111111111111111111111111111111111111111111111111111111111111@5?@5?@5?@5?@5?@5?@$�@5?@V@ff@�+@�@��?��R?�1?��w?���?�o?{��?u?o\)?k�?b��?T�j?St�?Rn�?KC�?I7L?2�!?(1'?�H?�?I�?�? Ĝ>�ff>��
>��>�J>�  >��h>��;>��>Kƨ>'�=��=���=@�<�h�o�+��%���ͽ���z�G��l�D��J��\)������徬�D��KǾǮ��hs�߾w������H����	7L����!��+�����ۿ"J�%�˿*=q�,�D�/\)11111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
�HB
�BB
�HB
�BB
�BB
�HB
�BB
�BB
�HB
�HB
�TB
��B@�B�!B�!B�B�B��B��B�B�!B�!B�!B�B�9B�9B�3B�?B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�'B�B�!B�'B�'B�'B�-B�'B�-B�'B�!B�!B�'B�!B�'B�-B�-B�-B�3B�3B�9B�3B�3B�-B�-B�'B�'B�'B�'B�'B�-B�-B�-B�3B�3B�3B�3B�333333333333333333333333333333333333333333333333333333333333333333333333333333333G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�44444444444444444444444444444444444444444444444444444444444444444444444444444444G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED=PSAL (re-calculated by using PRES_ADJUSTED) + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                           Surface pressure = -0.3 dbar                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            OW : r=1 , vertically averaged dS =-0.011338                                                                                                                                                                                                                    Pressure adjusted by using pressure offset at the sea surface. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                     No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Severe sensor malfunction, measurements are bad and not adjustable                                                                                                                                                                                              202105221602332021052216023320210522160233  IF  ARFMCODA027d                                                                20190116130203                      G�O�G�O�G�O�                IF  ARGQCOQC3.8                                                                 20190116130237  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC3.8                                                                 20190116130237  QCF$                G�O�G�O�G�O�0000000000000000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology 20190605092715  IP  PSAL            @�33D�G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20210522160233  IP  PSAL            @�33D�G�O�                