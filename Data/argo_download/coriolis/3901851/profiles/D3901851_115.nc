CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  6   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2020-01-27T11:41:33Z creation; 2022-10-09T12:10:06Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_032h      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
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
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  C�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  E�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  N�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  P�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  Y�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  b�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  d�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  m�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  o�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  x�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 8  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �    HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �$   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �(   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �,   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �l   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �|   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��Argo profile    3.1 1.2 19500101000000  20200127114133  20221009121006  3901851 MOCCA-POLAND                                                    Waldemar Walczowski                                             PRES            TEMP            PSAL               sA   IF                                  2C  D   ARVOR                           AI2600-16FR014                  5900A00                         844 @�������1   @������@Q%��.��5uF�k�w1   GPS     A   B   B   Primary sampling: averaged [10 sec sampling, 4 dbar average from 2000 dbar to 1400 dbar; 10 sec sampling, 2 dbar average from 1400 dbar to 250 dbar; 10 sec sampling, 1 dbar average from 250 dbar to 2.1 dbar]                                                    @��@9��@s33@�33@�33@�33A   A��A!��A.ffAA��AS33A^ffAq��A���A�ffA�ffA�  A�  A�ffA�33A���A�ffA���A�33A�  A�33A�ffA���A�33B��B33B33B33B��B��B��B ��B$ffB'��B+��B.��B2ffB6ffB:  B?��BE��BI��BM33BQ33BU33BX��B\��B`ffBg��Bk��Bp  Bs33Bw��B|��B�  B���B�  B���B�ffB�33B�  B���B�ffB�33B���B���B�ffB�33B�  B�  B���B���B�ffB�ffB�33B�33B�  B���B���B���B���B���B�ffB�ffB�ffB�ffB�ffB�ffB�33B�33B�33B�33B�33B�33B�33B�33B�ffB�ffB�ffB�ffBܙ�B���B���B���B�  B�33B�33B�ffB홚BB�B���B�  B�  B�  B�  B�33C �C�C�C�C�C�C�C33C33C
�C
��C  C33C  C�3C�fC�C��C  C33C  C�fC��C��C  C��C  C33C  C��C��C�fC!�C"  C"��C#��C$�fC&�C'  C'��C(��C)�fC+33C,�C-�C-�fC.��C/��C0��C1� C2�fC4L�C533C633C7�C8  C8�fC9��C:��C;�3C<��C=�fC?ffC@L�CA33CB33CC  CD  CE  CE�fCF��CG��CH�3CI�3CJ�3CK��CL��CM� CN��CO� CP� CQ� CR� CS� CT� CU� CV� CW��CX� CY� CZ��C[� C\��C]��C^� C_��C`��Ca��Cb�3Cc�3Cd�3Ce�3Cf�3Cg�3Ch��Ci�3Cj��Ck��Cl��Cm��Cn�fCo�fCp�fCq�fCr�fCt  Cu  Cv�Cw�Cx�Cy33Cz�3C|��C~��C�33C�  C�  C��C��C��C�&fC�33C�33C�33C�@ C�@ C�@ C�L�C�ffC�ffC�33C�  C��C��C��C��C��C�&fC�33C�33C�33C�33C�@ C�@ C�L�C�Y�C�ffC�ffC�s3C�33C�  C�  C�  C�s3C�L�C��C��C�&fC�33C�33C�@ C�L�C�&fC�33C�@ C��C�&fC�@ C��C�33C�L�C�&fC��C�&fC�@ C�&fC��C�33C�L�C�@ C�&fC��C�&fC�L�C�33C��C�@ C�Y�C�L�C�33C��C��C�33C�Y�C�L�C�@ C�33C�&fC�&fC��C��C�33C�ffC�Y�C�L�C�@ C�33C�&fC�&fC��C��C��C��C��C��C��C��C��C���C�&fC�L�C�L�C�L�C�L�C�L�C�L�C�L�C�L�C�&fC�33C�L�C��C�&fC�33C�@ C�Y�C�33C��C�&fC�@ C��C�@ D ,�D � D�D� D,�D� D3D�fD�D�fD�D��DfD�3D&fD� D3D��D	  D	��D
33D
�fD&fD� D�D�3D3D�3D3D��D�D�fDfD� D  D� D�D��D9�D��D9�D��D33D�3D,�D��D,�D��D&fD�fD  D� D�D�3D3D��DfD� D�D�3D,�D�fD  D� D   D ��D!�D!�3D"�D"��D#�D#�fD$fD$�fD%  D%y�D&�D&�3D',�D'��D(&fD(� D)9�D)��D*3D*�3D+3D+��D,  D,� D-  D-� D.  D.� D/&fD/�fD0�D0�3D1�D1� D2&fD2�3D3�D3�fD43D4��D5  D5��D6�D6�fD73D7� D8&fD8��D9fD9�3D:  D:�3D;  D;�3D<  D<�3D=  D=��D>  D>�3D?fD?�3D@&fD@��DAfDA��DB,�DB� DC�DC�3DD�DD�fDE�DE�3DF,�DF�fDG  DG��DH3DH��DIfDI�fDJ  DJ��DK9�DK��DL33DL��DM,�DM��DN,�DN�fDO&fDO�fDP  DP� DQ  DQ� DR  DR�fDS&fDS�fDT&fDT�fDU,�DU�3DV33DV� DW  DW� DXfDX��DY3DY�3DZ�DZ� D[&fD[�fD\&fD\��D]�D]� D^,�D^� D^��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@9��@s33@�33@�33@�33A   A��A!��A.ffAA��AS33A^ffAq��A���A�ffA�ffA�  A�  A�ffA�33A���A�ffA���A�33A�  A�33A�ffA���A�33B��B33B33B33B��B��B��B ��B$ffB'��B+��B.��B2ffB6ffB:  B?��BE��BI��BM33BQ33BU33BX��B\��B`ffBg��Bk��Bp  Bs33Bw��B|��B�  B���B�  B���B�ffB�33B�  B���B�ffB�33B���B���B�ffB�33B�  B�  B���B���B�ffB�ffB�33B�33B�  B���B���B���B���B���B�ffB�ffB�ffB�ffB�ffB�ffB�33B�33B�33B�33B�33B�33B�33B�33B�ffB�ffB�ffB�ffBܙ�B���B���B���B�  B�33B�33B�ffB홚BB�B���B�  B�  B�  B�  B�33C �C�C�C�C�C�C�C33C33C
�C
��C  C33C  C�3C�fC�C��C  C33C  C�fC��C��C  C��C  C33C  C��C��C�fC!�C"  C"��C#��C$�fC&�C'  C'��C(��C)�fC+33C,�C-�C-�fC.��C/��C0��C1� C2�fC4L�C533C633C7�C8  C8�fC9��C:��C;�3C<��C=�fC?ffC@L�CA33CB33CC  CD  CE  CE�fCF��CG��CH�3CI�3CJ�3CK��CL��CM� CN��CO� CP� CQ� CR� CS� CT� CU� CV� CW��CX� CY� CZ��C[� C\��C]��C^� C_��C`��Ca��Cb�3Cc�3Cd�3Ce�3Cf�3Cg�3Ch��Ci�3Cj��Ck��Cl��Cm��Cn�fCo�fCp�fCq�fCr�fCt  Cu  Cv�Cw�Cx�Cy33Cz�3C|��C~��C�33C�  C�  C��C��C��C�&fC�33C�33C�33C�@ C�@ C�@ C�L�C�ffC�ffC�33C�  C��C��C��C��C��C�&fC�33C�33C�33C�33C�@ C�@ C�L�C�Y�C�ffC�ffC�s3C�33C�  C�  C�  C�s3C�L�C��C��C�&fC�33C�33C�@ C�L�C�&fC�33C�@ C��C�&fC�@ C��C�33C�L�C�&fC��C�&fC�@ C�&fC��C�33C�L�C�@ C�&fC��C�&fC�L�C�33C��C�@ C�Y�C�L�C�33C��C��C�33C�Y�C�L�C�@ C�33C�&fC�&fC��C��C�33C�ffC�Y�C�L�C�@ C�33C�&fC�&fC��C��C��C��C��C��C��C��C��C���C�&fC�L�C�L�C�L�C�L�C�L�C�L�C�L�C�L�C�&fC�33C�L�C��C�&fC�33C�@ C�Y�C�33C��C�&fC�@ C��C�@ D ,�D � D�D� D,�D� D3D�fD�D�fD�D��DfD�3D&fD� D3D��D	  D	��D
33D
�fD&fD� D�D�3D3D�3D3D��D�D�fDfD� D  D� D�D��D9�D��D9�D��D33D�3D,�D��D,�D��D&fD�fD  D� D�D�3D3D��DfD� D�D�3D,�D�fD  D� D   D ��D!�D!�3D"�D"��D#�D#�fD$fD$�fD%  D%y�D&�D&�3D',�D'��D(&fD(� D)9�D)��D*3D*�3D+3D+��D,  D,� D-  D-� D.  D.� D/&fD/�fD0�D0�3D1�D1� D2&fD2�3D3�D3�fD43D4��D5  D5��D6�D6�fD73D7� D8&fD8��D9fD9�3D:  D:�3D;  D;�3D<  D<�3D=  D=��D>  D>�3D?fD?�3D@&fD@��DAfDA��DB,�DB� DC�DC�3DD�DD�fDE�DE�3DF,�DF�fDG  DG��DH3DH��DIfDI�fDJ  DJ��DK9�DK��DL33DL��DM,�DM��DN,�DN�fDO&fDO�fDP  DP� DQ  DQ� DR  DR�fDS&fDS�fDT&fDT�fDU,�DU�3DV33DV� DW  DW� DXfDX��DY3DY�3DZ�DZ� D[&fD[�fD\&fD\��D]�D]� D^,�D^� D^��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@N��@N�R@N��@N�+@N��@Nv�@Nff@Nff@O�P@S"�@X�`@b��@d(�@Qhs@5��@0r�@:�@��?�R?-��?��#@;d@l�@E�@1@v�@9�#@^ȴ@7;d@
^5?�b>|푽��w��|�:^5���D��Q쿔9X�4z�U�]�-�i7L�v�+��%�]p��I�^�]p����F����{��푿�Z�[�m�R�!���5?}�'�=�S�>���>��<�h>��`>���?�dZ?��R?�t�?�O�?��H?��?�-?ٙ�@�@��@ȴ@S�@��@��@1@ff@�R@ȴ@x�@�D@�-@�@-@�@ Ĝ?�dZ?�7L?�j?��?�b?���?��?���?��?�?�~�?׮?��m?š�?�J?ļj?ě�?��
?��?�E�?�=q?���?�?}?�?}?�?}?�Z?ě�?ļj?�Z?��?�K�?�$�?��`?��?���?�`B?��?��?��D?��T?�ȴ?ӶF?ȓu?�?�=q?�`B?��w?��?�?�+?���?��?��?�?�5??��?�x�?��P?�
=?�ff?�E�?�?���?��?��?��7?~5??q&�?j=q?i7L?g�?\(�?R�?Kƨ?I�^?H1'?E`B?6?2-?0�`?0�`?-O�?+ƨ?+C�?)��?)�^?(�9?)��?(�9?"��?!��?!G�? �? �?�w?�w?dZ?��?9X?�?�?��?C�?��?��?$�?`B?�/?Z?�/?�?�\>���>�?}>�->� �>��>>�V>�h>>�X>��>�X>��j>�F>��>��>���>� �>�->�>��y>��/>ܬ>�(�>ڟ�>ؓu>ۥ�>�"�>�"�>�b>ؓu>�/>ܬ>�"�>��;>���>��`>�ƨ>�=q>ɺ^>�ƨ>���>׍P>��>ٙ�>�\)>�7L>��>�  >�j>��m>��m>�X>�K�>�?}>� �>�~�>�r�>�A�>�/>���>��>�b>�bN>�O�>�I�>�I�>�I�>�ƨ>�C�>��9>�1'>���>���>���>��>�  >�%>�%>��>|�>~��>�  >�  >}�>u>q��>o��>n��>n��>p��>o��>q��>s�F>s�F>u>vȴ>|�>�  >�%>�o>�o>��>���>�+>�+>�+>�1'>��^>�7L>�C�>�O�>�\)>�bN>��;>��`>�z�>��u>���>�"�>��->��R>�5?>��->��->�/>�/>��->��->��>�/>�A�>�Ĝ>�G�>�M�>�M�>�Z>�ff>��T>�`B>�`B>�l�>�r�>�r�>���>���>�>�~�>���>�r�>��y>�ff>�ff>�`B>�S�>�G�>�A�>�A�>��R>��->�/>�/>�/>��->��->�5?>�5?>��>��>���>���>��>��>��u>��+>�>�>�>�>�>��>���>��>��`>��;>��>�O�>�I�>��>���>�$�>��>�o>�%>~��>|�>z�H>vȴ>t�j>q��>l�D>ix�>fff>\(�>Xb>W
=>W
=>W
=>V>R�>O�;>M��>H�9>D��>A�7>>v�><j>9X>7K�>6E�>5?}>49X>1&�>,1>%�T>!��>�w>��>�>hs>bN>I�>
=q>+>�>%=��=��=�=�h=�l�=�"�=�"�=��`=ȴ9=ě�=� �=��
=��w=��w=��P=�O�=�%=u=��=y�#=}�=u=ix�=Y�=H�9=H�9=@�=49X=#�
=�P=\)<�<ě�<�j<�j<�j<u<T��<o;��
;��
;��
;�o;D��    �D����t���1��9X�ě����t��'8Q�H�9�P�`�m�h��C���hs��hs��hs���w�� Ž�j�\����
=��"ѽ�G���+��P�����R�&�y�'-V�0 ž49X�6E��9X�;dZ�<j�>vɾ@��D���G��L�;O�;�R�V�W
=�Xb�Z��]/�]/�]/�]/�]/�`A��e`B�e`B�k��r�!�w�پz�H���������\�����+���^���������녾��������A����徧�xվ�xվ�~���~���~���~���~�������11111111111111111111111111111111111111111111111111111111441111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @N��@N�R@N��@N�+@N��@Nv�@Nff@Nff@O�P@S"�@X�`@b��@d(�@Qhs@5��@0r�@:�@��?�R?-��?��#@;d@l�@E�@1@v�@9�#@^ȴ@7;d@
^5?�b>|푽��w��|�:^5���D��Q쿔9X�4z�U�]�-�i7L�v�+��%�]p��I�^�]p����F����{��푿�Z�[�m�R�!���5?}G�O�G�O�>���>��<�h>��`G�O�G�O�?��R?�t�?�O�?��H?��?�-?ٙ�@�@��@ȴ@S�@��@��@1@ff@�R@ȴ@x�@�D@�-@�@-@�@ Ĝ?�dZ?�7L?�j?��?�b?���?��?���?��?�?�~�?׮?��m?š�?�J?ļj?ě�?��
?��?�E�?�=q?���?�?}?�?}?�?}?�Z?ě�?ļj?�Z?��?�K�?�$�?��`?��?���?�`B?��?��?��D?��T?�ȴ?ӶF?ȓu?�?�=q?�`B?��w?��?�?�+?���?��?��?�?�5??��?�x�?��P?�
=?�ff?�E�?�?���?��?��?��7?~5??q&�?j=q?i7L?g�?\(�?R�?Kƨ?I�^?H1'?E`B?6?2-?0�`?0�`?-O�?+ƨ?+C�?)��?)�^?(�9?)��?(�9?"��?!��?!G�? �? �?�w?�w?dZ?��?9X?�?�?��?C�?��?��?$�?`B?�/?Z?�/?�?�\>���>�?}>�->� �>��>>�V>�h>>�X>��>�X>��j>�F>��>��>���>� �>�->�>��y>��/>ܬ>�(�>ڟ�>ؓu>ۥ�>�"�>�"�>�b>ؓu>�/>ܬ>�"�>��;>���>��`>�ƨ>�=q>ɺ^>�ƨ>���>׍P>��>ٙ�>�\)>�7L>��>�  >�j>��m>��m>�X>�K�>�?}>� �>�~�>�r�>�A�>�/>���>��>�b>�bN>�O�>�I�>�I�>�I�>�ƨ>�C�>��9>�1'>���>���>���>��>�  >�%>�%>��>|�>~��>�  >�  >}�>u>q��>o��>n��>n��>p��>o��>q��>s�F>s�F>u>vȴ>|�>�  >�%>�o>�o>��>���>�+>�+>�+>�1'>��^>�7L>�C�>�O�>�\)>�bN>��;>��`>�z�>��u>���>�"�>��->��R>�5?>��->��->�/>�/>��->��->��>�/>�A�>�Ĝ>�G�>�M�>�M�>�Z>�ff>��T>�`B>�`B>�l�>�r�>�r�>���>���>�>�~�>���>�r�>��y>�ff>�ff>�`B>�S�>�G�>�A�>�A�>��R>��->�/>�/>�/>��->��->�5?>�5?>��>��>���>���>��>��>��u>��+>�>�>�>�>�>��>���>��>��`>��;>��>�O�>�I�>��>���>�$�>��>�o>�%>~��>|�>z�H>vȴ>t�j>q��>l�D>ix�>fff>\(�>Xb>W
=>W
=>W
=>V>R�>O�;>M��>H�9>D��>A�7>>v�><j>9X>7K�>6E�>5?}>49X>1&�>,1>%�T>!��>�w>��>�>hs>bN>I�>
=q>+>�>%=��=��=�=�h=�l�=�"�=�"�=��`=ȴ9=ě�=� �=��
=��w=��w=��P=�O�=�%=u=��=y�#=}�=u=ix�=Y�=H�9=H�9=@�=49X=#�
=�P=\)<�<ě�<�j<�j<�j<u<T��<o;��
;��
;��
;�o;D��    �D����t���1��9X�ě����t��'8Q�H�9�P�`�m�h��C���hs��hs��hs���w�� Ž�j�\����
=��"ѽ�G���+��P�����R�&�y�'-V�0 ž49X�6E��9X�;dZ�<j�>vɾ@��D���G��L�;O�;�R�V�W
=�Xb�Z��]/�]/�]/�]/�]/�`A��e`B�e`B�k��r�!�w�پz�H���������\�����+���^���������녾��������A����徧�xվ�xվ�~���~���~���~���~�������11111111111111111111111111111111111111111111111111111111441111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB+B,B+B(�B(�B'�B'�B)�BE�Be`B��B"�BF�B�{B{�B�JB�fB{�B� Bt�B+B�B�#B�BH�Bx�BÖB{�B(�BB�B�/B��B��BɺB�BB�BR�B�=B�bB��B��B�LB��BǮB��B��B��B	B	JB��B	�B	�B	<jB	VB	iyB	9XB	t�B	��B	��B	��B	�
B	�uB
B
W
B
W
B
hsB
^5B
{�B
�B
��B
�#B
��B
��B
�/B
�/B
�B
�B
�B
��B
��BB
�NB
�B
�B
�B
�B
��B
�B
��B
�B
�B
�B
�B
�fB  B
�B
��B
��B
�B
�B
�B
�B
��B
�B
�B
��B
��B
��B
��B
��B  B
��B
��BB1B%BBDB+BJB%B%B	7B1BbBbB�B/B-B&�B�B �B �B#�B$�B$�B&�B&�B33B)�B.B/B1'B0!B.B0!B1'B0!B33B0!B1'B1'B-B2-B-B.B.B.B1'B49B/B0!B0!B49B1'B2-B1'B49B33B33B49B49B49B33B5?B7LB5?B7LB7LB6FB7LB6FB8RB;dB:^B:^B;dB;dB<jB=qB=qB=qB>wB>wB>wB>wB?}B?}BC�B@�BA�BA�BA�BA�BC�BE�BC�BA�BG�BI�BG�BG�BH�BF�BF�BJ�BH�BI�BJ�BK�BL�BJ�BK�BM�BN�BL�BM�BL�BM�BQ�BM�BM�BO�BO�BM�BR�BP�BQ�BQ�BR�BT�BS�BT�BS�BR�BYBR�BT�BT�BT�BVBW
BW
BW
BW
BW
BXBZBZBZBYBYBZB[#B[#B[#B[#B[#B[#B[#B]/B\)B]/B]/B^5B_;B^5B^5B_;B`BBaHBaHBaHBbNBbNBbNBbNBbNBcTBdZBe`BffBffBffBffBgmBhsBhsBjBjBk�Bl�Bm�Bn�Bn�Bo�Bo�Bp�Bp�Bq�Bq�Bs�Bs�Bt�Bt�Bv�Bw�Bw�By�By�Bz�Bz�Bz�Bz�B{�B{�B{�B|�B|�B~�B~�B~�B~�B� B�B�B�B�B�B�B�B�%B�B�+B�+B�1B�1B�7B�7B�=B�=B�=B�DB�=B�DB�DB�DB�DB�JB�JB�JB�JB�JB�PB�PB�\B�\B�bB�bB�bB�hB�hB�oB�oB�oB�oB�oB�oB�oB�uB�uB�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111441111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B+B,B+B(�B(�B'�B'�B)�BE�Be`B��B"�BF�B�{B{�B�JB�fB{�B� Bt�B+B�B�#B�BH�Bx�BÖB{�B(�BB�B�/B��B��BɺB�BB�BR�B�=B�bB��B��B�LB��BǮB��B��B��B	B	JB��B	�B	�B	<jB	VG�O�G�O�B	t�B	��B	��B	��G�O�G�O�B
B
W
B
W
B
hsB
^5B
{�B
�B
��B
�#B
��B
��B
�/B
�/B
�B
�B
�B
��B
��BB
�NB
�B
�B
�B
�B
��B
�B
��B
�B
�B
�B
�B
�fB  B
�B
��B
��B
�B
�B
�B
�B
��B
�B
�B
��B
��B
��B
��B
��B  B
��B
��BB1B%BBDB+BJB%B%B	7B1BbBbB�B/B-B&�B�B �B �B#�B$�B$�B&�B&�B33B)�B.B/B1'B0!B.B0!B1'B0!B33B0!B1'B1'B-B2-B-B.B.B.B1'B49B/B0!B0!B49B1'B2-B1'B49B33B33B49B49B49B33B5?B7LB5?B7LB7LB6FB7LB6FB8RB;dB:^B:^B;dB;dB<jB=qB=qB=qB>wB>wB>wB>wB?}B?}BC�B@�BA�BA�BA�BA�BC�BE�BC�BA�BG�BI�BG�BG�BH�BF�BF�BJ�BH�BI�BJ�BK�BL�BJ�BK�BM�BN�BL�BM�BL�BM�BQ�BM�BM�BO�BO�BM�BR�BP�BQ�BQ�BR�BT�BS�BT�BS�BR�BYBR�BT�BT�BT�BVBW
BW
BW
BW
BW
BXBZBZBZBYBYBZB[#B[#B[#B[#B[#B[#B[#B]/B\)B]/B]/B^5B_;B^5B^5B_;B`BBaHBaHBaHBbNBbNBbNBbNBbNBcTBdZBe`BffBffBffBffBgmBhsBhsBjBjBk�Bl�Bm�Bn�Bn�Bo�Bo�Bp�Bp�Bq�Bq�Bs�Bs�Bt�Bt�Bv�Bw�Bw�By�By�Bz�Bz�Bz�Bz�B{�B{�B{�B|�B|�B~�B~�B~�B~�B� B�B�B�B�B�B�B�B�%B�B�+B�+B�1B�1B�7B�7B�=B�=B�=B�DB�=B�DB�DB�DB�DB�JB�JB�JB�JB�JB�PB�PB�\B�\B�bB�bB�bB�hB�hB�oB�oB�oB�oB�oB�oB�oB�uB�uB�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��11111111111111111111111111111111111111111111111111111111441111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              202210091210062022100912100620221009121006  IF  ARFMCODA032h                                                                20200127114133                      G�O�G�O�G�O�                IF  ARGQCOQC4.4                                                                 20200127114233  QCP$                G�O�G�O�G�O�000000000008FB7EIF  ARGQCOQC4.4                                                                 20200127114233  QCF$                G�O�G�O�G�O�0000000000004000GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200212132919  IP  PSAL            @��D^��G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20200830134206  IP  PSAL            @��D^��G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology 20210308172402  IP  PSAL            @��D^��G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V01 + ARGO climatology 20211025140934  IP  PSAL            @��D^��G�O�                GE  ARSQOW  2.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20221009121006  IP  PSAL            @��D^��G�O�                