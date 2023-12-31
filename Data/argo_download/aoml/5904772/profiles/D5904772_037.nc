CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:21Z creation      
references        (http://www.argodatamgt.org/Documentation   comment       	free text      user_manual_version       3.2    Conventions       Argo-3.2 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         (Matthew Alkire, University of Washington      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7,   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  74   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7t   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     88   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8X   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8x   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~       axis      T           8|   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~            8�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    A�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    ]�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  _�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    g�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  i�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  q�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    y�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  {�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    �|   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �|   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �t   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �    HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �0   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �4   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �D   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �H   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �L   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004021  20190405134355  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               %A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�ZJ��z1   @�Z�$@N*�G�{�B��Q�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    %A   B   B   @���@���A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C5�fC8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@�fDAfDA�fDB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt� Dy� D��3D�0 D���D��3D��fD�@ D�� D���D�fD�33D��3D���D�fD�0 D�6fD�ɚD�  D�P D�p D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�  @�33A33A+33AK33Ak33A���A���A���A���Ař�Aՙ�A噚A���B��B
��B��B��B"��B*��B2��B:��BB��BJ��BR��BZ��Bb��Bj��Br��Bz��B�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffC �3C�3C�3C�3C�3C
�3C�3C�3C�3C�3C�3C�3C�3C�3C�3C�3C �3C"�3C$�3C&�3C(�3C*�3C,�3C.�3C0�3C2�3C4�3C6��C8�3C:�3C<�3C>�3C@�3CB�3CD�3CF�3CH�3CJ�3CL�3CN�3CP�3CR�3CT�3CV�3CX�3CZ�3C\�3C^�3C`�3Cb�3Cd�3Cf�3Ch�3Cj�3Cl�3Cn�3Cp�3Cr�3Ct�3Cv�3Cx�3Cz�3C|�3C~�3C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�D ,�D ��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D	,�D	��D
,�D
��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D ,�D ��D!,�D!��D",�D"��D#,�D#��D$,�D$��D%,�D%��D&,�D&��D',�D'��D(,�D(��D),�D)��D*,�D*��D+,�D+��D,,�D,��D-,�D-��D.,�D.��D/,�D/��D0,�D0��D1,�D1��D2,�D2��D3,�D3��D4,�D4��D5,�D5��D6,�D6��D7,�D7��D8,�D8��D9,�D9��D:,�D:��D;,�D;��D<,�D<��D=,�D=��D>,�D>��D?,�D?��D@,�D@�3DA33DA�3DB,�DB��DC,�DC��DD,�DD��DE,�DE��DF,�DF��DG,�DG��DH,�DH��DI,�DI��DJ,�DJ��DK,�DK��DL,�DL��DM,�DM��DN,�DN��DO,�DO��DP,�DP��DQ,�DQ��DR,�DR��DS,�DS��DT,�DT��DU,�DU��DV,�DV��DW,�DW��DX,�DX��DY,�DY��DZ,�DZ��D[,�D[��D\,�D\��D],�D]��D^,�D^��D_,�D_��D`,�D`��Da,�Da��Db,�Db��Dc,�Dc��Dd,�Dd��De,�De��Df,�Df��Dg,�Dg��Dh,�Dh��Di,�Di��Dj,�Dj��Dk,�Dk��Dl,�Dl��Dm,�Dm��Dn,�Dn��Do,�Do��Dp,�Dp��Dq,�Dq��Dr,�Dr��Ds,�Ds��Dt,�Dt��Du�Dy��D�	�D�FfD�� D�ٙD��D�VfD��fD��3D�,�D�I�D���D��3D�,�D�FfD�L�D�� D�fD�ffD�fD�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@��P@�t�@�l�@�dZ@�l�@�l�@�l�@�l�@�t�@�l�@�l�@�l�@�t�@�t�@�l�@�l�@�t�@�t�@�|�@�t�@�t�@�t�@�t�@�|�@��@��@�|�@�|�@�|�@�t�@�|�@�|�@��@��@��@��@��@��@��@��P@��@��@��@�dZ@�\)@�S�@�o@�=q@�$�@���@���@�O�@�V@��j@�1'@~5?@}�@z~�@yX@v�+@u@up�@t(�@st�@s@q�^@q�#@q�7@p�9@p�u@pA�@pb@o�@o�@o��@o��@o�w@o��@o��@o�P@o�P@o�P@o��@o��@o��@o�P@o|�@ol�@ol�@ol�@ol�@o|�@o|�@o��@o�P@o��@oK�@n��@nv�@nff@nV@nV@m��@m��@m�@m�@l��@l�/@lj@l(�@l1@k�m@k�F@kt�@kC�@k33@k"�@j��@j^5@jM�@j=q@jJ@i��@i��@i&�@h�u@hA�@g�@g�P@g�P@g
=@fv�@fV@fV@fȴ@g�w@h  @h  @g�;@g�;@g�w@g�w@g�w@g��@hb@g�;@g��@g�@g�@g�@g�w@g�@h  @h  @hb@h �@h �@h �@h �@h1'@hA�@hbN@h�@hr�@hbN@hbN@hr�@hr�@hbN@hQ�@hA�@h �@h1'@h �@hA�@hr�@h�u@h�u@h�u@hbN@hQ�@hQ�@h �@h  @g��@g��@g|�@g\)@gK�@gK�@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@g\)@gl�@g\)@g\)@g\)@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@g\)@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gK�@gK�@gK�@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@g\)@g\)@g\)@g\)@g\)@gK�@g+@f�y@fȴ@f��@fv�@fV@fv�@f�+@f�+@f�+@f�+@f�+@fff@fff@fV@fV@fff@fff@fv�@fv�@fff@fff@fff@fV@fV@f5?@f5?@f$�@f{@f@e�T@e��@e�-@e��@e��@e�h@e�h@e�h@ep�@eO�@e?}@e/@e�@e�@e�@d��@d��@d�/@d��@d�j@d��@d��@d��@d�D@d�D@d�D@dz�@dZ@dj@dZ@dI�@d9X@d(�@d�@d1@c�
@c�
@c�
@c�
@c�
@c�
@c�
@c�F@c��@c��@ct�@cdZ@cS�@cC�@c"�@co@co@co@co@co@c@c@b�@b��@b��@b��@b�!@b��@b~�@b^5@bM�@bM�@bM�@b=q@b=q@b=q@b-@b-@b�@bJ@a�@a�@a�@a�@a�@a��@a�^@a�^@a��@a��@ahs@aX@a7L@a&�@a�@a�@a%@`��@`�9@`�@`r�@`bN@`bN@`bN@`A�@`  @_�@_�;@_��@_�w@_�@_�@_�@_�@_�P@_|�@_|�@_l�@_;d@_+@_
=@_
=@^�y@^�@^�@^�R@^�+@^E�@^5?@^@]��@]@]�-@]�-@]�-@]�-@]�-@]�-@]��@]�h@]�@]`B@]O�@]?}@]�@]V@\��@\��@\��@\�@\�@\�/@\�/@\�/@\�/@\�/@\��@\��@\��@\�j@\�j@\�j@\�@\�D@\z�@\j@\j@\Z@\I�@\I�@\9X@\9X@\9X@\(�@\(�@\9X@\9X@\9X@\9X@\(�@\�@\1@[��@[��@[��@[�m@[�
@[�m@[�
@[ƨ@[�F@[��@[��@[�@[t�@[t�@[dZ@[33@Z��@Y�@X�`@Xr�@XA�@W�;@Wl�@W�@Wl�@XA�@Z�@`�`@fff@fv�@fV@d9X@a�7@^��@\�@Y�#111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @���@��P@�t�@�l�@�dZ@�l�@�l�@�l�@�l�@�t�@�l�@�l�@�l�@�t�@�t�@�l�@�l�@�t�@�t�@�|�@�t�@�t�@�t�@�t�@�|�@��@��@�|�@�|�@�|�@�t�@�|�@�|�@��@��@��@��@��@��@��@��P@��@��@��@�dZ@�\)@�S�@�o@�=q@�$�@���@���@�O�@�V@��j@�1'@~5?@}�@z~�@yX@v�+@u@up�@t(�@st�@s@q�^@q�#@q�7@p�9@p�u@pA�@pb@o�@o�@o��@o��@o�w@o��@o��@o�P@o�P@o�P@o��@o��@o��@o�P@o|�@ol�@ol�@ol�@ol�@o|�@o|�@o��@o�P@o��@oK�@n��@nv�@nff@nV@nV@m��@m��@m�@m�@l��@l�/@lj@l(�@l1@k�m@k�F@kt�@kC�@k33@k"�@j��@j^5@jM�@j=q@jJ@i��@i��@i&�@h�u@hA�@g�@g�P@g�P@g
=@fv�@fV@fV@fȴ@g�w@h  @h  @g�;@g�;@g�w@g�w@g�w@g��@hb@g�;@g��@g�@g�@g�@g�w@g�@h  @h  @hb@h �@h �@h �@h �@h1'@hA�@hbN@h�@hr�@hbN@hbN@hr�@hr�@hbN@hQ�@hA�@h �@h1'@h �@hA�@hr�@h�u@h�u@h�u@hbN@hQ�@hQ�@h �@h  @g��@g��@g|�@g\)@gK�@gK�@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@g\)@gl�@g\)@g\)@g\)@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@g\)@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gK�@gK�@gK�@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@gl�@g\)@g\)@g\)@g\)@g\)@gK�@g+@f�y@fȴ@f��@fv�@fV@fv�@f�+@f�+@f�+@f�+@f�+@fff@fff@fV@fV@fff@fff@fv�@fv�@fff@fff@fff@fV@fV@f5?@f5?@f$�@f{@f@e�T@e��@e�-@e��@e��@e�h@e�h@e�h@ep�@eO�@e?}@e/@e�@e�@e�@d��@d��@d�/@d��@d�j@d��@d��@d��@d�D@d�D@d�D@dz�@dZ@dj@dZ@dI�@d9X@d(�@d�@d1@c�
@c�
@c�
@c�
@c�
@c�
@c�
@c�F@c��@c��@ct�@cdZ@cS�@cC�@c"�@co@co@co@co@co@c@c@b�@b��@b��@b��@b�!@b��@b~�@b^5@bM�@bM�@bM�@b=q@b=q@b=q@b-@b-@b�@bJ@a�@a�@a�@a�@a�@a��@a�^@a�^@a��@a��@ahs@aX@a7L@a&�@a�@a�@a%@`��@`�9@`�@`r�@`bN@`bN@`bN@`A�@`  @_�@_�;@_��@_�w@_�@_�@_�@_�@_�P@_|�@_|�@_l�@_;d@_+@_
=@_
=@^�y@^�@^�@^�R@^�+@^E�@^5?@^@]��@]@]�-@]�-@]�-@]�-@]�-@]�-@]��@]�h@]�@]`B@]O�@]?}@]�@]V@\��@\��@\��@\�@\�@\�/@\�/@\�/@\�/@\�/@\��@\��@\��@\�j@\�j@\�j@\�@\�D@\z�@\j@\j@\Z@\I�@\I�@\9X@\9X@\9X@\(�@\(�@\9X@\9X@\9X@\9X@\(�@\�@\1@[��@[��@[��@[�m@[�
@[�m@[�
@[ƨ@[�F@[��@[��@[�@[t�@[t�G�O�@[33@Z��@Y�@X�`@Xr�@XA�@W�;@Wl�@W�@Wl�@XA�@Z�@`�`@fff@fv�@fV@d9X@a�7@^��@\�@Y�#111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�oB�hB�oB�hB�oB�hB�hB�hB�hB�bB�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�JB�PB�JB�JB�JB�JB�DB�7B�7B�1B�1B�1B�+B�%B�%B�%B�+B�7B�=B�=B�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B� B� B� B~�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B|�B|�B|�B|�B{�B{�B{�B{�B{�B{�B{�B{�B{�Bz�Bz�Bz�Bz�Bz�By�By�By�By�By�Bx�Bx�Bx�Bw�Bw�Bw�Bw�Bw�Bw�Bw�Bw�Bw�Bw�Bw�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bt�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Bs�Br�Br�Br�Br�Br�Bq�Bp�Bn�Bm�Bl�Bl�Bk�Bk�BjBk�Bn�Bw�B�JB��B��B�B�B�B�B�'B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B�0B�-B�-B�-B�.B�.B�0B�-B�.B�.B�.B�.B�.B�.B�-B�-B�.B�-B�-B�0B�.B�-B�.B�.B�,B�.B�.B�-B�0B�-B�.B�-B�0B�0B�-B�,B�.B�1B�.B�.B�.B�0B�0B�.B�.B�-B�.B�+B�)B�*B�)B�#B�"B�B�B�B�
B�B�B�B�
B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B~�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B|�B|�B|�B|�B{�B{�B{�B{�B{�B{�B{�B{�B{�Bz�Bz�Bz�Bz�Bz�By|By|By{ByzBy{BxuBxuBxvBwnBwoBwoBwqBwnBwnBwnBwoBwnBwmBwoBviBviBviBvfBvhBviBvhBvhBucBucBucBudBudBubBucBueBueBudBueBucBudBucBucBt]Bt]Bt\Bt_Bt]Bt_Bt_Bt^Bt]Bt\Bt]Bt\Bt]Bt\Bt]Bt^Bt^BsVBsVBsWBsTBsWBsWBsVBsVBsUBsZBsXBrPBrOBrOBrQG�O�BqJBpCBn9Bm4Bl-Bl.Bk(Bk(Bj#Bk(Bn<BwtB��B�eB��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.7 dbar.                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051343552019040513435520190405134355  AO  ARCAADJP                                                                    20181006004021    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004021  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004021  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134355  IP                  G�O�G�O�G�O�                