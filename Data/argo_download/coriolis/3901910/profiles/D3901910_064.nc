CDF   #   
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2018-09-15T20:43:41Z creation; 2019-06-15T17:32:26Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  �  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  �  8`   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  `  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9@   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9H   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9L   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                  @  9P   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                  @  9�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                  @  9�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                  @  :   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    :\   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         :d   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    :t   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            :x   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           :�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           :�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    :�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    :�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    :�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    :�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    :�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  J`   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  M�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  [T   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  ^�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  lH   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  y�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  }<   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �0   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  �P   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 d  �D   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  �  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ތ   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ި   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                     ް   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                     ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  �8   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ˘   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ј   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ט   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T  ݘ   T  ݘ          T  ݘArgo profile    3.1 1.2 19500101000000  20180915204341  20190615173226  3901910 3901910 MOCCA-POLAND                                                    MOCCA-POLAND                                                    Waldemar Walczowski                                             Waldemar Walczowski                                             PRES            TEMP            PSAL            PRES            TEMP            PSAL               @   @AA  IFIF                                                                2C  2B  DA  ARVOR                           ARVOR                           AI2600-16FR073                  AI2600-16FR073                  5900A01                         5900A01                         844 844 @؁���6�@؁���6�11  @؁�e�8 @؁�e�8 @T#���K@T#���K@ETu�L@ETu�L11  GPS     GPS     AA  BC  BF  Primary sampling: averaged [10 sec sampling, 5 dbar average from 2000 dbar to 1000 dbar; 10 sec sampling, 2 dbar average from 1000 dbar to 200 dbar; 10 sec sampling, 1 dbar average from 200 dbar to 2.3 dbar]                                                 Near-surface sampling: averaged, unpumped [10 sec sampling, 1 dbar average from 2.3 dbar to surface]                                                                                                                                                                  @9��@s33@���@���@���@���A  A$��A333AA��AP  A\��Al��A{33A�  A�33A�ffA���A���A�  A�33A�ffAř�A���A�  A�33A�ffA�A���B ffBffBffB  B  B  B  B��B��B#��B'33B+33B/33B2��B6��B:��B>��BBffBFffBK��BQ33BU33BX��B\ffB`ffBd  Bh  Bl  Bp  Bt  Bx  B{��B��B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B���B���B�33B���B�ffB���B�33B���B�ffB�  B�ffB�  B���B�33B���B�33B�  B˙�B�33B���B�ffB�  Bՙ�B�ffB�  Bܙ�B�ffB�  B���B㙚B�33B���BꙚB�ffB�33B���B���B�B�33B�  B���B���B�ffC �C  C�fC��C�3C��C� C�fCL�C	33C
�C�C  C�fC��C�3C��C�fCL�CL�C33C�C  C�fC��C�3C��C� C�fCL�C33C�C   C �fC!��C"�3C$  C%L�C&33C'  C'�fC(�3C*  C+L�C,�C-  C-�fC.�3C/�fC133C2  C2�fC3�3C4�fC633C7�C7�fC8�3C:  C;33C<�C=  C=�fC>�3C?�fCAL�CB�CC  CC�fCD��CE�3CF��CHffCJ��CLffCN33CP� CR�3CT� CV33CX� CZ��C\��C^ffC`33Cb� Cd�3Cf� ChL�Cj�Cl� Cn��Cp�3Cr� Ct� CvL�Cx33Cz�C|  C~ffC�ffC�Y�C�Y�C�L�C�@ C�&fC��C��C��C�  C�  C��3C�  C�  C�  C�  C�  C�  C��C��C��C��C��C��C��C�&fC�&fC�&fC�33C�33C�@ C�L�C�L�C�Y�C�ffC�33C�  C��C��C�&fC�33C�33C�@ C�L�C�Y�C�ffC�33C�  C��C�33C�L�C�Y�C�33C�  C��C�&fC�33C�@ C�L�C�Y�C�ffC�s3C�@ C��C��C�&fC�33C�@ C�L�C�ffC�@ C��C��C�33C�L�C�ffC�33C�  C��C�&fC�33C�@ C�L�C�Y�C�ffC�33C�  C��C��C�33C�@ C�L�C�Y�C�s3C�@ C��C�&fC�@ C�L�C�ffC�@ C��C�&fC�@ C�Y�C�@ C��C�33C�Y�C�33C��C�33C�L�C�33C�  C�&fC�@ C�Y�C�33C��C�&fC�@ C�Y�C�@ C��C�33C�Y�C�33D fD ��D,�D��D�D� D33D�fD3D��D  D�3D33D��D,�D�fD  D��D	3D	��D
&fD
��D3D��D&fD� D�D�3DfD��D  D��D3D�fD&fD��D3D��D&fD��D33D�fD  D�3DfD� D3D�fD  D�3DfD� D�D�3D33D�fD�D�3DfD�fD�D�3D,�D�fD   D ��D!3D!�3D"�D"��D#�D#��D$fD$��D%,�D%��D%��D&�3D'  D'ffD(L�D(�3>���?���@   G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                               @9��@s33@���@���@���@���A  A$��A333AA��AP  A\��Al��A{33A�  A�33A�ffA���A���A�  A�33A�ffAř�A���A�  A�33A�ffA�A���B ffBffBffB  B  B  B  B��B��B#��B'33B+33B/33B2��B6��B:��B>��BBffBFffBK��BQ33BU33BX��B\ffB`ffBd  Bh  Bl  Bp  Bt  Bx  B{��B��B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B�  B���B���B�33B���B�ffB���B�33B���B�ffB�  B�ffB�  B���B�33B���B�33B�  B˙�B�33B���B�ffB�  Bՙ�B�ffB�  Bܙ�B�ffB�  B���B㙚B�33B���BꙚB�ffB�33B���B���B�B�33B�  B���B���B�ffC �C  C�fC��C�3C��C� C�fCL�C	33C
�C�C  C�fC��C�3C��C�fCL�CL�C33C�C  C�fC��C�3C��C� C�fCL�C33C�C   C �fC!��C"�3C$  C%L�C&33C'  C'�fC(�3C*  C+L�C,�C-  C-�fC.�3C/�fC133C2  C2�fC3�3C4�fC633C7�C7�fC8�3C:  C;33C<�C=  C=�fC>�3C?�fCAL�CB�CC  CC�fCD��CE�3CF��CHffCJ��CLffCN33CP� CR�3CT� CV33CX� CZ��C\��C^ffC`33Cb� Cd�3Cf� ChL�Cj�Cl� Cn��Cp�3Cr� Ct� CvL�Cx33Cz�C|  C~ffC�ffC�Y�C�Y�C�L�C�@ C�&fC��C��C��C�  C�  C��3C�  C�  C�  C�  C�  C�  C��C��C��C��C��C��C��C�&fC�&fC�&fC�33C�33C�@ C�L�C�L�C�Y�C�ffC�33C�  C��C��C�&fC�33C�33C�@ C�L�C�Y�C�ffC�33C�  C��C�33C�L�C�Y�C�33C�  C��C�&fC�33C�@ C�L�C�Y�C�ffC�s3C�@ C��C��C�&fC�33C�@ C�L�C�ffC�@ C��C��C�33C�L�C�ffC�33C�  C��C�&fC�33C�@ C�L�C�Y�C�ffC�33C�  C��C��C�33C�@ C�L�C�Y�C�s3C�@ C��C�&fC�@ C�L�C�ffC�@ C��C�&fC�@ C�Y�C�@ C��C�33C�Y�C�33C��C�33C�L�C�33C�  C�&fC�@ C�Y�C�33C��C�&fC�@ C�Y�C�@ C��C�33C�Y�C�33D fD ��D,�D��D�D� D33D�fD3D��D  D�3D33D��D,�D�fD  D��D	3D	��D
&fD
��D3D��D&fD� D�D�3DfD��D  D��D3D�fD&fD��D3D��D&fD��D33D�fD  D�3DfD� D3D�fD  D�3DfD� D�D�3D33D�fD�D�3DfD�fD�D�3D,�D�fD   D ��D!3D!�3D"�D"��D#�D#��D$fD$��D%,�D%��D%��D&�3D'  D'ffD(L�D(�3>���?���@   G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                               @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@uO�@�ȴ@���@�Ĝ@�{@��@�Ĝ@���@��F@�Q�@�j@�r�@���@���@��/@��@�p�@���@��@�E�@�J@�$�@�~�@���@���@�ȴ@���@��R@���@�ff@�v�@��+@��\@���@���@��R@���@��R@���@��y@�@�;d@�l�@�t�@��@��@��
@��m@��w@���@��\@���@��\@�V@�{@��7@�bN@�|�@�;d@�;d@���@�M�@�7L@��@�l�@�-@�`B@�hs@�j@�l�@���@�7L@�&�@��@�V@���@���@��7@��@���@���@��@�X@���@�j@��@�Ĝ@�O�@��@�r�@�r�@�j@�r�@�bN@�I�@���@�ƨ@�\)@�ff@��^@��@�?}@��@�%@��`@��u@�A�@���@�V@��^@�7L@��/@�"�@�j@�C�@�M�@��w@�&�@�V@�;d@��
@�O�@�(�@�V@���@�7L@�G�@��\@��!@��H@���@�G�@�p�@�j@��@�1@��;@�;d@�X@�/@�&�@��D@�z�@�|�@��y@���@�n�@�=q@���@��\@�@�
=@�
=@��;@�;d@�=q@��@�V@���@�x�@���@���@��`@��R@��!@��\@��
@}O�@p�@k�
@kdZ@j��@hĜ@fE�@f��@g�w@eV@b�\@`A�@_|�@_l�@^ff@Xr�@V��@V�y@W�@W�P@Y�#@^V@^V@g+@jM�@k@p  @p�u@p�`@q%@q%@p�9@p�u@o�@o
=@n�y@nv�@n$�@mO�@l1@i%@h��@e�@`Q�@Yhs@Q�^@J�H@F{@C�m@B�H@B�\@A��@?�@>�y@>ff@=�T@=p�@=�@=��@=�-@=��@=�-@=�-@=��@=�T@=�@=`B@<��@<z�@9�#@3t�@2~�@2M�@1�^@1�^@2=q@2n�@1�@1��@1�7@1�@2��@6ȴ@6�@5/@4I�@5�@5�@4(�@3"�@1��@1��@1��@1hs@1X@1X@1G�@1X@1x�@1�7@1�7@1�7@17L@0�9@0��@1&�@1�#@2=q@2��@2�H@2�H@2��@2n�@1��@1�@1��@1�#@1��@0��@0 �@/�w@.��@,Z@+�F@+�F@+dZ@,1@,9X@,�@,(�@+ƨ@+@)��@';d@%?}@$9X@#�@#�
@&ff@%O�@"��@"��@#@"�@"��@"J@�;@l�@��@@C�@��@�j@@�^@hs@�9@;d@�R@E�@�h@C�@
�!@
~�@	�@�R@��@I�@@�@ �u?�|�?��?��?�S�?� �?�{?�ff?�5??���?�7L?�l�?Լj?�V?��y?�`B?��/?�o?�J?���?���?�V?��H?�`B?�;d?�I�?��F?���?�\)?�C�?��#?�K�?��j?���?��\?���?���?�?}?���?�A�?|�?|j?t9X?p �?e�?^v�?[�m?YX?XQ�?Rn�?NV?H�9?C�
?@A�?>v�?3��?ƨ?�y?S�>��m>�Ĝ>�I�>��>�  >��H>�\)>dZ>.{=���=�C�=D����󶽙����E������xվO߾+�cS������hs���Ͼ��+�������������������������������������������������㾟;d@<��@;C�@HĜG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111311                                                                                                                                                                                                                                                                                                                                                                                                                                               @uO�@�ȴ@���@�Ĝ@�{G�O�G�O�@���@��F@�Q�@�j@�r�@���@���@��/@��@�p�@���@��@�E�@�J@�$�@�~�@���@���@�ȴ@���@��R@���@�ff@�v�@��+@��\@���@���@��R@���@��R@���@��y@�@�;d@�l�@�t�@��@��@��
@��m@��w@���@��\@���@��\@�V@�{@��7@�bN@�|�@�;d@�;d@���@�M�@�7L@��@�l�@�-@�`B@�hs@�j@�l�@���@�7L@�&�@��@�V@���@���@��7@��@���@���@��@�X@���@�j@��@�Ĝ@�O�@��@�r�@�r�@�j@�r�@�bN@�I�@���@�ƨ@�\)@�ff@��^@��@�?}@��@�%@��`@��u@�A�@���@�V@��^@�7L@��/@�"�@�j@�C�@�M�@��w@�&�@�V@�;d@��
@�O�@�(�@�VG�O�G�O�@�G�@��\@��!@��H@���@�G�@�p�@�j@��@�1@��;@�;d@�X@�/@�&�@��D@�z�@�|�@��y@���@�n�@�=q@���@��\@�@�
=@�
=@��;@�;d@�=q@��@�V@���@�x�@���@���@��`@��R@��!@��\@��
@}O�G�O�G�O�@kdZ@j��@hĜ@fE�@f��@g�w@eV@b�\@`A�@_|�@_l�@^ff@Xr�@V��@V�y@W�@W�P@Y�#@^V@^V@g+@jM�@k@p  @p�u@p�`@q%@q%@p�9@p�u@o�@o
=@n�y@nv�@n$�@mO�@l1@i%@h��@e�@`Q�@Yhs@Q�^@J�H@F{@C�m@B�H@B�\@A��@?�@>�y@>ff@=�T@=p�@=�@=��@=�-@=��@=�-@=�-@=��@=�T@=�@=`B@<��@<z�@9�#@3t�@2~�@2M�@1�^@1�^@2=q@2n�@1�@1��@1�7@1�@2��@6ȴ@6�@5/@4I�@5�@5�@4(�@3"�@1��@1��@1��@1hs@1X@1X@1G�@1X@1x�@1�7@1�7@1�7@17L@0�9@0��@1&�@1�#@2=q@2��@2�H@2�H@2��@2n�@1��@1�@1��@1�#@1��@0��@0 �@/�w@.��@,Z@+�F@+�F@+dZ@,1@,9X@,�@,(�@+ƨ@+@)��@';d@%?}@$9X@#�@#�
@&ff@%O�@"��@"��@#@"�@"��@"J@�;@l�@��@@C�@��@�j@@�^@hs@�9@;d@�R@E�@�h@C�@
�!@
~�@	�@�R@��@I�@@�@ �u?�|�?��?��?�S�?� �?�{?�ff?�5??���?�7L?�l�?Լj?�V?��y?�`B?��/?�o?�J?���?���?�V?��H?�`B?�;d?�I�?��F?���?�\)?�C�?��#?�K�?��j?���?��\?���?���?�?}?���?�A�?|�?|j?t9X?p �?e�?^v�?[�m?YX?XQ�?Rn�?NV?H�9?C�
?@A�?>v�?3��?ƨ?�y?S�>��m>�Ĝ>�I�>��>�  >��H>�\)>dZ>.{=���=�C�=D����󶽙����E������xվO߾+�cS������hs���Ͼ��+�������������������������������������������������㾟;d@<��@;C�@HĜG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114411111111111111111111111111111111111111111144111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111311                                                                                                                                                                                                                                                                                                                                                                                                                                               ;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�BR�B	�^B	��B
B
��B,B�B"�BM�BVBYB\)B^5BbNBaHBe`BjBk�Bo�Bl�Bl�Bo�BhsBs�Bv�Bu�Bt�Bt�Bv�Bs�Bs�Bs�Bt�Bu�Bu�Bu�Bv�Bw�Bv�Bz�B}�B|�B��B��B��B�B��B��B�/B�B�B�B�B��B��BBDBB+B+B
=BPB�B�B�B�B�B!�B7LB-B49B/B33B2-B49B33B6FB2-B(�B'�B(�B �BoB%B�B49B7LBD�BH�BK�BF�BH�BI�BJ�BK�BK�BK�BS�BT�BO�BO�BR�BO�BP�BP�BO�BO�B[#BQ�BQ�BO�BM�B?}BJ�B7LB2-B.B@�B:^BI�B=qB+B�B!�BJ�BPBB)�B'�B,B-B1'BA�BD�BA�B@�BA�BB�BB�BQ�BA�B=qBC�BA�B;dB;dB:^B8RB9XB;dB:^B<jB?}BE�BI�BE�BJ�B9XB33B49B.B49B#�B7LB7LB49B5?B$�B\)B�B��B��B��B��B�B��BVB��B��B��B�B��B�yB�mB�mB�B�mB�BBDB��B�B�B1'B.B,B.B/B2-B/B/B/B33B1'B/B1'B2-B)�B(�B,B�B�B\B  B�B�B�B�yB�B�yB�yB�B�sB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�HB�5B�5B�5B�/B�HB�ZB�ZB�TB�ZB�`B�yB�B��B��B�B��B�B�B��B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B  B  B  B  B��B��B��B  B  BBB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B��B��B�B�yB�B�mB�`B�fB�sB�ZB�`B�TB�TB�HB�BB�HB�BB�5B�;B�/B�5B�ZB�)B�5B�B�B�)B�B�B�B��B��B��B��B��B��BƨBŢBƨBǮBǮBĜBŢBǮB��BȴB��B��BĜBƨB��B��B��BB��BĜB�FB��B�}BÖBB��BBB��BÖBÖBĜBBBÖBÖB�}B�wB�qB�}B�}B�B�XB�dB�dB�LB�XB�9B�3B�XB�9B�FB�9B�-B�B�!B�!B��B�B�B��B��B��B�qB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B	��B	$�Bu�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111441111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111411111411441114414441111111111114111111111111111141111141111144111141411114111114114111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111333                                                                                                                                                                                                                                                                                                                                                                                                                                               BR�B	�^B	��B
B
��G�O�G�O�B"�BM�BVBYB\)B^5BbNBaHBe`BjBk�Bo�Bl�Bl�Bo�BhsBs�Bv�Bu�Bt�Bt�Bv�Bs�Bs�Bs�Bt�Bu�Bu�Bu�Bv�Bw�Bv�Bz�B}�B|�B��B��B��B�B��B��B�/B�B�B�B�B��B��BBDBB+B+B
=BPB�B�B�B�B�B!�B7LB-B49B/B33B2-B49B33B6FB2-B(�B'�B(�B �BoG�O�G�O�B49B7LBD�BH�BK�BF�BH�BI�BJ�BK�BK�BK�BS�BT�BO�BO�BR�BO�BP�BP�BO�BO�G�O�BQ�BQ�BO�BM�B?}G�O�B7LB2-G�O�G�O�B:^BI�B=qG�O�G�O�B!�G�O�G�O�G�O�B)�B'�B,B-B1'BA�BD�BA�B@�BA�BB�BB�G�O�BA�B=qBC�BA�B;dB;dB:^B8RB9XB;dB:^B<jB?}BE�BI�BE�G�O�B9XB33B49B.B49G�O�B7LB7LB49B5?B$�G�O�G�O�B��B��B��B��G�O�B��G�O�B��B��B��B�G�O�B�yB�mB�mB�B�mG�O�BBDG�O�B�B�B1'B.B,B.B/B2-B/B/B/B33B1'B/B1'B2-B)�B(�B,B�B�B\B  B�B�B�B�yB�B�yB�yB�B�sB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�HB�5B�5B�5B�/B�HB�ZB�ZB�TB�ZB�`B�yB�B��B��B�B��B�B�B��B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B  B  B  B  B��B��B��B  B  BBB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B��B��B�B�yB�B�mB�`B�fB�sB�ZB�`B�TB�TB�HB�BB�HB�BB�5B�;B�/B�5B�ZB�)B�5B�B�B�)B�B�B�B��B��B��B��B��B��BƨBŢBƨBǮBǮBĜBŢBǮB��BȴB��B��BĜBƨB��B��B��BB��BĜB�FB��B�}BÖBB��BBB��BÖBÖBĜBBBÖBÖB�}B�wB�qB�}B�}B�B�XB�dB�dB�LB�XB�9B�3B�XB�9B�FB�9B�-B�B�!B�!B��B�B�B��B��B��B�qB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B	��B	$�Bu�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111441111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111411111411441114414441111111111114111111111111111141111141111144111141411114111114114111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111333                                                                                                                                                                                                                                                                                                                                                                                                                                               <#�
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
G�O�<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
G�O�G�O�<#�
<#�
<#�
G�O�G�O�<#�
G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
G�O�G�O�<#�
<#�
<#�
<#�
G�O�<#�
G�O�<#�
<#�
<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = a1 * PSAL + a0                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Not applicable                                                                                                                                                                                                                                                  Not applicable                                                                                                                                                                                                                                                  a1=1, a0=0.0000                                                                                                                                                                                                                                                 No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              No adjustment performed (values duplicated)                                                                                                                                                                                                                     No adjustment performed (values duplicated)                                                                                                                                                                                                                     Real-time salinity adjustment                                                                                                                                                                                                                                   201906151732262019061517322620190615173226201809152043412018091520434120180129132202IF  IF  ARFMARFMCODACODA024c024c                                                                                                                                2018091520434120180915204341                                        G�O�G�O�G�O�G�O�G�O�G�O�                                IF  IF  ARGQARGQCOQCCOQC3.5 3.5                                                                                                                                 2018091520435320180915204353QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�000000000008FB7E0000000000689B7EIF  IF  ARGQARGQCOQCCOQC3.5 3.5                                                                                                                                 2018091520435320180915204353QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�00000000000040000000000000600000IF      CORT    COOA    6.2     RTQCGL01                                                                                                                        20180917160030              QCF$    TEMP                                G�O�B�  G�O�G�O�G�O�5                               IF      CORT    COOA    6.2     RTQCGL01                                                                                                                        20180917162923              QCF$    PSAL                                G�O�B�  G�O�G�O�G�O�5                               IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C>�3G�O�C>�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            CB�G�O�CB�G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   B�33G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   B�ffG�O�B�ffG�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   BꙚG�O�BꙚG�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   B���G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            B�33G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C3�3G�O�C3�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            BꙚG�O�BꙚG�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            B���G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            B���G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C�fG�O�C�fG�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   B���G�O�B���G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C�fG�O�C�fG�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C�3G�O�C�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C   G�O�C   G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C8�3G�O�C8�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C>�3G�O�C>�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   CB�G�O�CB�G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C&33G�O�C&33G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL_ADJUSTED                   C2  G�O�C2  G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C8�3G�O�C8�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C3�3G�O�C3�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C2  G�O�C2  G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C&33G�O�C&33G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C   G�O�C   G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            C�3G�O�C�3G�O�?�  G�O�                                IF              SCOO    0.39                                                                                                                                    20180918100607              CF      PSAL                            B�ffG�O�B�ffG�O�?�  G�O�                                IF              COCQ    3.2                                                                                                                                     20180918084122                                                      G�O�G�O�G�O�G�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2018V02 + ARGO climatology                                                                 20190615173226              IP      PSAL                            @9��G�O�D(�3G�O�G�O�G�O�                                