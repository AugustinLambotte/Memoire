CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2017-04-26T09:43:49Z creation; 2018-06-23T14:07:29Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8X   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     8�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     8�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9    JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9(   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >��	4E�   
_FillValue        A.�~            9,   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           94   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9<   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9D   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9H   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9P   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9T   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9X   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9\   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :\   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  :`   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  B(   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  D   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ]h   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _\   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g$   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  i   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  p�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  x�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  z�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �d   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  �X   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �|   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �    SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �P   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �P   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �P   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �P  �P             ,  �PArgo profile    3.1 1.2 19500101000000  20170426094349  20180623140729  3901874 MOCCA-EU                                                        Romain Cancouet                                                 PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR                           AR2600-16FR037                  5900A00                         844 @��!�q�1   @����J @R��Q�?��G�z�1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 250 dbar; 10 sec sampling, 2 dbar average from 250 dbar to 1400 dbar; 10 sec sampling, 4 dbar average from 1400 dbar to 1000 dbar]                                                     A(  AI��AfffA�33A�33A�  A���A���Aՙ�A�ffA�ffB  B
  B33B��B!��B)��B133B8��B@��BI33BQ33BZ  BbffBj  Bq��BzffB�33B�  B�  B�  B���B���B���B���B�33B�  B���B�33B�33B�  B���B�ffB���B�33B�ffB���BЙ�B���B�33B�  B���B�  B�ffB���B�ffB���B���B���C ��C�3CffCL�CffC
��CffC33CffC��CffCL�CffC��CL�C33C L�C"ffC$��C&ffC(33C*ffC,��C.�3C0ffC233C4ffC6ffC8��C:�3C<ffC>33C@L�CB� CD��CF�3CHffCJ33CLL�CN� CP��CRffCT33CVL�CX� CZ��C\� C^L�C`� Cb��CdffCf�Ch33CjL�ClL�CnffCpffCrffCtffCv� Cx� Cz� C|��C~��C�Y�C�33C��C��C�&fC�&fC�33C�@ C�L�C�33C��C�&fC�33C�@ C�L�C�33C��C�33C�L�C�Y�C�@ C�&fC�33C�L�C�Y�C�33C��C�33C�@ C�L�C�Y�C�@ C�&fC�33C�@ C�@ C�L�C�Y�C�33C��C��C�&fC�33C�33C�33C�33C�33C�@ C�@ C�L�C�L�C�L�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�33C�33C�Y�C�Y�C�Y�C�L�C�L�C�L�C�L�C�@ C�@ C�33C�&fC��C�33C�L�C�@ C�33C�&fC�33C�L�C�33C�&fC�33C�@ C�&fC�33C�@ C�@ C�L�C�33C�33C�@ C�@ C�L�C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�33C�33C�&fC�@ C�33C�&fC�33C�@ C�33C��C�&fC�33C�@ C�@ C�@ C�L�C�L�C�L�C�L�C�Y�C�L�C�L�C�@ D   D �3D�D��D�D��D  D� D&fD�3D�D�fD�D�3D�D�3D  D��D	fD	��D
�D
� D�D� D&fD��D  D��D3D�3D  D��D3D�fD�D�3D�D�3D�D� D3D� D3D��D�D��D3D�3D  D��D�D��D�D�3D�D�3D�D��D�D�3D  D�3D3D�3D �D �3D!�D!�3D"�D"��D#�D#��D$3D$��D%�D%�3D&�D&��D'�D'��D(3D(��D)3D)�3D*3D*��D+�D+� D,3D,�3D-�D-��D.3D.� D/3D/��D0�D0�3D1�D1��D23D2�3D3�D3��D4�D4��D53D5��D6�D6� D7  D7� D8fD8��D93D9��D:  D:�3D;�D;��D<�D<��D=3D=�3D>3D>��D?�D?�fD@�D@��DA�DA��DB�DB�3DC3DC��DD  DD� DE�DE��DF�DF��DGfDG��DH�DH� DI  DI��DJ3DJ��DK�DK��DL3DL�3DM�DM��DN3DN��DO  DO��DP�DP�3DQ�DQ��DR  DR�fDS�DS��DT3DT��DU�DU� DV3DV�3DW  DW�3DX�DX��DY�DY�3DZ�DZ��D[3D[��D\3D\�3D]�D]��D^�D^�3D_3D_�3D`3D`�3Da3Da�3Db3Db�3Dc3Dc��Dd�Dd�3De�De�3Df�Df��Dg�Dg�3Dh�Dh�3Di3Di��Dj3Dj��Dk�Dk�3Dl�Dl��Dm3Dm��Dn�Dn��Do3Do�3Dp3Dp��Dq3Dq��Dr�Dr��Ds�Ds�3Dt�Dt�3Du3Du�3Dv�Dv��Dw3Dw��Dx�Dx�3Dy�Dy��Dz�Dz��D{�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A(  AI��AfffA�33A�33A�  A���A���Aՙ�A�ffA�ffB  B
  B33B��B!��B)��B133B8��B@��BI33BQ33BZ  BbffBj  Bq��BzffB�33B�  B�  B�  B���B���B���B���B�33B�  B���B�33B�33B�  B���B�ffB���B�33B�ffB���BЙ�B���B�33B�  B���B�  B�ffB���B�ffB���B���B���C ��C�3CffCL�CffC
��CffC33CffC��CffCL�CffC��CL�C33C L�C"ffC$��C&ffC(33C*ffC,��C.�3C0ffC233C4ffC6ffC8��C:�3C<ffC>33C@L�CB� CD��CF�3CHffCJ33CLL�CN� CP��CRffCT33CVL�CX� CZ��C\� C^L�C`� Cb��CdffCf�Ch33CjL�ClL�CnffCpffCrffCtffCv� Cx� Cz� C|��C~��C�Y�C�33C��C��C�&fC�&fC�33C�@ C�L�C�33C��C�&fC�33C�@ C�L�C�33C��C�33C�L�C�Y�C�@ C�&fC�33C�L�C�Y�C�33C��C�33C�@ C�L�C�Y�C�@ C�&fC�33C�@ C�@ C�L�C�Y�C�33C��C��C�&fC�33C�33C�33C�33C�33C�@ C�@ C�L�C�L�C�L�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�33C�33C�Y�C�Y�C�Y�C�L�C�L�C�L�C�L�C�@ C�@ C�33C�&fC��C�33C�L�C�@ C�33C�&fC�33C�L�C�33C�&fC�33C�@ C�&fC�33C�@ C�@ C�L�C�33C�33C�@ C�@ C�L�C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�@ C�33C�33C�&fC�@ C�33C�&fC�33C�@ C�33C��C�&fC�33C�@ C�@ C�@ C�L�C�L�C�L�C�L�C�Y�C�L�C�L�C�@ D   D �3D�D��D�D��D  D� D&fD�3D�D�fD�D�3D�D�3D  D��D	fD	��D
�D
� D�D� D&fD��D  D��D3D�3D  D��D3D�fD�D�3D�D�3D�D� D3D� D3D��D�D��D3D�3D  D��D�D��D�D�3D�D�3D�D��D�D�3D  D�3D3D�3D �D �3D!�D!�3D"�D"��D#�D#��D$3D$��D%�D%�3D&�D&��D'�D'��D(3D(��D)3D)�3D*3D*��D+�D+� D,3D,�3D-�D-��D.3D.� D/3D/��D0�D0�3D1�D1��D23D2�3D3�D3��D4�D4��D53D5��D6�D6� D7  D7� D8fD8��D93D9��D:  D:�3D;�D;��D<�D<��D=3D=�3D>3D>��D?�D?�fD@�D@��DA�DA��DB�DB�3DC3DC��DD  DD� DE�DE��DF�DF��DGfDG��DH�DH� DI  DI��DJ3DJ��DK�DK��DL3DL�3DM�DM��DN3DN��DO  DO��DP�DP�3DQ�DQ��DR  DR�fDS�DS��DT3DT��DU�DU� DV3DV�3DW  DW�3DX�DX��DY�DY�3DZ�DZ��D[3D[��D\3D\�3D]�D]��D^�D^�3D_3D_�3D`3D`�3Da3Da�3Db3Db�3Dc3Dc��Dd�Dd�3De�De�3Df�Df��Dg�Dg�3Dh�Dh�3Di3Di��Dj3Dj��Dk�Dk�3Dl�Dl��Dm3Dm��Dn�Dn��Do3Do�3Dp3Dp��Dq3Dq��Dr�Dr��Ds�Ds�3Dt�Dt�3Du3Du�3Dv�Dv��Dw3Dw��Dx�Dx�3Dy�Dy��Dz�Dz��D{�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@ǶF@ǥ�@ǝ�@ǅ@ǅ@�S�@�@��@Ł@��@�bN@��@�n�@�hs@��@�p�@t(�@lI�@f��@`�u@W|�@P�`@H�9@>{@6ȴ@/�P@+ƨ@+"�@*�\@&��@$1@\)@�m@r�@�@�@V@Q�@��@��@�
@n�@   ?���?�t�?�?��T?�M�?ݑh?��?��/?��`?�~�?�9X?�p�?���?��\?�\)?��?�b?�n�?���?�j?���?��
?�Ĝ?��;?��?�&�?�M�?��7?�v�?�=q?�1'?���?��?~5??yX?tz�?o��?j=q?g+?c��?`A�?]/?Z^5?W��?St�?Q&�?NV?M�h?I�^?I7L?E�T?E��?Fff?F$�?E`B?D�?@Ĝ?=�?:��?7K�?5?}?2n�?0��?.{?*=q?%��?"�\? Ĝ?�R?dZ?�?�u?�+?�?�j?��?9X?hs?bN?\)?V?{?�h?I�?	��?�y?��?Z?�y?Z?o?%? �? �>�v�>��m>��#>�K�>�>�F>�&�>>�h>�D>�>�~�>�x�>�l�>�S�>�G�>�A�>޸R>�/>ۥ�>��>�>�hs>�I�>�C�>ɺ^>�+>��>�$�>Ƨ�>ě�>�%>�p�>���>�>�33>���>���>�{>�V>���>��T>�Ĝ>���>�
=>�bN>�C�>��>���>���>���>��>��\>�J>��7>}�>t�j>m�h>fff>cS�>`A�>V>G�>7K�>&�y>��>�>t�>V>o=�F=�`B=�l�=���=ȴ9=��=�v�=�j=�^5=�1=��=�C�=�O�=�+=ix�=L��=D��=@�=8Q�=�w=C�=o=o=t�=�w=#�
='�=,1='�=#�
=�P=o<�1<�o<�C�<�C�<�o<o    ��o��`B�t��u��o��t���j�ě���j��j��j��j�u�T���u���
���ͼ��C��<j�P�`�]/�aG��ixսu��o�����7L��\)��hs��hs��hs��\)��\)��t�������t���+��7L��+��C�������w�������Ƨ�ě��Ƨ�����������`��������;d��;d��G���S���`B��l�����F�����#�����$ݾ1'�O߾bN�O߾I��bN�t��z�z�z�z������P��P���������w�#�
�(�þ+�0 ž8Q�=p��>vɾ?|�@��B�\�B�\�@��@��C���F��I�^�H�9�J���L�;N��N��O�;�O�;�O�;�Q녾Q녾S�ϾT���V�Xb�Y��Y��Z��Z��["Ѿ\(��]/�\(��\(��\(��["ѾZ��Z��\(��cS��fff�fff�gl��gl��dZ�cS��dZ�fff�hr��ixվj~��j~��j~��j~��k��k��m�h�n���p�׾p�׾p�׾q���q���r�!�s�F�t�j�t�j�u�u�vȴ�w�پx���x���z�H�z�H�{�m�|푾}󶾀  ���������7��%���7���7���7��J��o����������������������˾��˾���+��+��1'���9��7L��=q�������;�O߾���V��V��\)��\)���;��bN���`���`���`��녾�񪾓t���zᾕ������P��b���u�����"Ѿ��㾜(���/���-���R���w������S���`B���T��ff���y���xվ����D���h������������������������������� ž��׾�&龱&龱&龱���������!��33���j��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ǶF@ǥ�@ǝ�@ǅ@ǅ@�S�@�@��@Ł@��@�bN@��@�n�@�hs@��@�p�@t(�@lI�@f��@`�u@W|�@P�`@H�9@>{@6ȴ@/�P@+ƨ@+"�@*�\@&��@$1@\)@�m@r�@�@�@V@Q�@��@��@�
@n�@   ?���?�t�?�?��T?�M�?ݑh?��?��/?��`?�~�?�9X?�p�?���?��\?�\)?��?�b?�n�?���?�j?���?��
?�Ĝ?��;?��?�&�?�M�?��7?�v�?�=q?�1'?���?��?~5??yX?tz�?o��?j=q?g+?c��?`A�?]/?Z^5?W��?St�?Q&�?NV?M�h?I�^?I7L?E�T?E��?Fff?F$�?E`B?D�?@Ĝ?=�?:��?7K�?5?}?2n�?0��?.{?*=q?%��?"�\? Ĝ?�R?dZ?�?�u?�+?�?�j?��?9X?hs?bN?\)?V?{?�h?I�?	��?�y?��?Z?�y?Z?o?%? �? �>�v�>��m>��#>�K�>�>�F>�&�>>�h>�D>�>�~�>�x�>�l�>�S�>�G�>�A�>޸R>�/>ۥ�>��>�>�hs>�I�>�C�>ɺ^>�+>��>�$�>Ƨ�>ě�>�%>�p�>���>�>�33>���>���>�{>�V>���>��T>�Ĝ>���>�
=>�bN>�C�>��>���>���>���>��>��\>�J>��7>}�>t�j>m�h>fff>cS�>`A�>V>G�>7K�>&�y>��>�>t�>V>o=�F=�`B=�l�=���=ȴ9=��=�v�=�j=�^5=�1=��=�C�=�O�=�+=ix�=L��=D��=@�=8Q�=�w=C�=o=o=t�=�w=#�
='�=,1='�=#�
=�P=o<�1<�o<�C�<�C�<�o<o    ��o��`B�t��u��o��t���j�ě���j��j��j��j�u�T���u���
���ͼ��C��<j�P�`�]/�aG��ixսu��o�����7L��\)��hs��hs��hs��\)��\)��t�������t���+��7L��+��C�������w�������Ƨ�ě��Ƨ�����������`��������;d��;d��G���S���`B��l�����F�����#�����$ݾ1'�O߾bN�O߾I��bN�t��z�z�z�z������P��P���������w�#�
�(�þ+�0 ž8Q�=p��>vɾ?|�@��B�\�B�\�@��@��C���F��I�^�H�9�J���L�;N��N��O�;�O�;�O�;�Q녾Q녾S�ϾT���V�Xb�Y��Y��Z��Z��["Ѿ\(��]/�\(��\(��\(��["ѾZ��Z��\(��cS��fff�fff�gl��gl��dZ�cS��dZ�fff�hr��ixվj~��j~��j~��j~��k��k��m�h�n���p�׾p�׾p�׾q���q���r�!�s�F�t�j�t�j�u�u�vȴ�w�پx���x���z�H�z�H�{�m�|푾}󶾀  ���������7��%���7���7���7��J��o����������������������˾��˾���+��+��1'���9��7L��=q�������;�O߾���V��V��\)��\)���;��bN���`���`���`��녾�񪾓t���zᾕ������P��b���u�����"Ѿ��㾜(���/���-���R���w������S���`B���T��ff���y���xվ����D���h������������������������������� ž��׾�&龱&龱&龱���������!��33���j��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�3B�3B�3B�3B�3B�3B�-B�RB�XB�jB�qB�wB��B�BbB�BhB{B�BbB�B�BoBJBuB	7BB1B
=BJB	7B+B+BB��B��B��B��B�B�B��B��B�B�B�B�B�B�yB�yB�fB�sB�NB�NB�/B�B�#B�B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBȴBǮBǮBǮBǮBŢBƨBƨBǮBŢBƨBŢBǮBƨBǮBƨBǮBŢBŢBŢBÖBÖBÖBÖBB��B��B�}B�}B��B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�}B�}B�}B�}B�wB�qB�qB�}B�}B�wB�wB�qB�wB�wB�}B�wB�wB�wB�qB�qB�qB�qB�qB�qB�wB�qB�jB�qB�qB�jB�jB�jB�jB�jB�dB�dB�dB�^B�^B�^B�dB�^B�dB�dB�^B�dB�^B�XB�XB�XB�RB�RB�RB�RB�RB�RB�LB�?B�FB�?B�?B�?B�?B�9B�9B�9B�9B�9B�9B�9B�3B�3B�-B�-B�3B�-B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  B�3B�3B�3B�3B�3B�3B�-B�RB�XB�jB�qB�wB��B�BbB�BhB{B�BbB�B�BoBJBuB	7BB1B
=BJB	7B+B+BB��B��B��B��B�B�B��B��B�B�B�B�B�B�yB�yB�fB�sB�NB�NB�/B�B�#B�B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBǮBȴBǮBǮBǮBǮBŢBƨBƨBǮBŢBƨBŢBǮBƨBǮBƨBǮBŢBŢBŢBÖBÖBÖBÖBB��B��B�}B�}B��B�}B�}B�}B�}B�}B�}B�wB�wB�wB�wB�wB�}B�}B�}B�}B�wB�qB�qB�}B�}B�wB�wB�qB�wB�wB�}B�wB�wB�wB�qB�qB�qB�qB�qB�qB�wB�qB�jB�qB�qB�jB�jB�jB�jB�jB�dB�dB�dB�^B�^B�^B�dB�^B�dB�dB�^B�dB�^B�XB�XB�XB�RB�RB�RB�RB�RB�RB�LB�?B�FB�?B�?B�?B�?B�9B�9B�9B�9B�9B�9B�9B�3B�3B�-B�-B�3B�-B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              201806231407292018062314072920180623140729  IF  ARFMCODA011b                                                                20170426094349                      G�O�G�O�G�O�                IF  ARGQCOQC3.1                                                                 20170426094356  QCP$                G�O�G�O�G�O�000000000008FB5EIF  ARGQCOQC3.1                                                                 20170426094356  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.0 CTD2017V1                                                       20171114122626  IP  PSAL            A(  D{�G�O�                GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology 20180623140729  IP  PSAL            A(  D{�G�O�                