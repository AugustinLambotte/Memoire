CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:26Z creation      
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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004026  20190405134359  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               6A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�1��_�1   @�1Ǜ�@NLj~��#�A�r� Ĝ1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    6A   B   B   @�33@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BHffBP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.y�D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� DtٚDy��D� D�9�D���D��fD�3D�FfD���D���D�3D�)�D��3D�� D���D�6fDڜ�D��3D�	�D�6fD�p D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��G@׮A�
A+�
AK�
Ak�
A��A��A��A��A��A��A��A��B��B
��B��B��B"��B*��B2��B:��BB��BK\)BR��BZ��Bb��Bj��Br��Bz��B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�z�B�B�z�B�z�B�z�B�z�B�z�C �qC�qC�qC�qC�qC
�qC�qC�qC�qC�qC�qC�qC�qC�qC�qC�qC �qC"�qC$�qC&�qC(�qC*�qC,�qC.�qC0�qC2�qC4�qC6�qC8�qC:�qC<�qC>�qC@�qCB�qCD�qCF�qCH�qCJ�qCL�qCN�qCP�qCR�qCT�qCV�qCX�qCZ�qC\�qC^�qC`�qCb�qCd�qCf�qCh�qCj�qCl�qCn�qCp�qCr�qCt�qCv�qCx�qCz�qC|�qC~�qC�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�C�^�D /\D �\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D	/\D	�\D
/\D
�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D/\D�\D /\D �\D!/\D!�\D"/\D"�\D#/\D#�\D$/\D$�\D%/\D%�\D&/\D&�\D'/\D'�\D(/\D(�\D)/\D)�\D*/\D*�\D+/\D+�\D,/\D,�\D-/\D-�\D./\D.��D//\D/�\D0/\D0�\D1/\D1�\D2/\D2�\D3/\D3�\D4/\D4�\D5/\D5�\D6/\D6�\D7/\D7�\D8/\D8�\D9/\D9�\D:/\D:�\D;/\D;�\D</\D<�\D=/\D=�\D>/\D>�\D?/\D?�\D@/\D@�\DA/\DA�\DB/\DB�\DC/\DC�\DD/\DD�\DE/\DE�\DF/\DF�\DG/\DG�\DH/\DH�\DI/\DI�\DJ/\DJ�\DK/\DK�\DL/\DL�\DM/\DM�\DN/\DN�\DO/\DO�\DP/\DP�\DQ/\DQ�\DR/\DR�\DS/\DS�\DT/\DT�\DU/\DU�\DV/\DV�\DW/\DW�\DX/\DX�\DY/\DY�\DZ/\DZ�\D[/\D[�\D\/\D\�\D]/\D]�\D^/\D^�\D_/\D_�\D`/\D`�\Da/\Da�\Db/\Db�\Dc/\Dc�\Dd/\Dd�\De/\De�\Df/\Df�\Dg/\Dg�\Dh/\Dh�\Di/\Di�\Dj/\Dj�\Dk/\Dk�\Dl/\Dl�\Dm/\Dm�\Dn/\Dn�\Do/\Do�\Dp/\Dp�\Dq/\Dq�\Dr/\Dr�\Ds/\Ds�\Dt/\Dt�\Du�Dy�)D�'�D�QHD��{D��D��D�^D��{D��HD��D�AHD���D��D�{D�NDڴ{D���D�!HD�ND�D��H111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�Ĝ@�@�@� �@�33@�33@�+@�
=@���@���@���@��y@��H@���@�ȴ@���@�R@�+@�=q@��#@�1'@�v�@�  @��@�?}@��!@���@�hs@���@���@�M�@�
=@�9X@�K�@�^5@��@��T@��@�33@�V@�G�@��m@��^@���@���@� �@�ƨ@��@�;d@��@�t�@��
@��m@��m@��m@�  @�1@�b@� �@��@�9X@���@��/@���@��`@���@��w@�;d@��H@���@��y@���@��@��/@�(�@�A�@��
@�hs@���@��7@��D@�r�@�z�@��m@�
=@�^5@�J@���@��-@���@��@�\)@�l�@��@���@��T@��@��@��-@�7L@��`@��`@���@���@��u@� �@�|�@���@��
@���@��R@�"�@��@��+@�$�@��#@���@��T@��@��^@�O�@�O�@�`B@��h@��@�/@���@���@���@���@��^@�`B@�G�@��@���@� �@��;@�+@�~�@�J@�p�@�%@�Ĝ@��@�{@��@��/@���@��@�G�@�X@�O�@�?}@�/@��@�%@���@��@�Ĝ@���@�I�@��m@���@�o@��`@��m@�\)@���@���@�&�@��@�r�@�1@�\)@�o@���@��R@��R@���@��+@�v�@�ff@�n�@�M�@�-@���@���@��-@���@���@���@���@�p�@�  @~E�@~5?@|�@�9X@�1'@� �@� �@�A�@�w@
=@}�T@|�@{�m@{�F@z��@y�@y�7@y%@x��@x�@xQ�@x1'@x �@xb@w�@x  @x1'@yhs@y��@y�@y��@y�@y�7@y&�@x��@x��@x�9@xbN@xb@w��@w|�@w�@vE�@u@t��@t�/@u�@uO�@uO�@t��@u�@u�@u`B@u��@up�@u/@t��@t�/@t�j@t��@tj@sƨ@s33@r�H@r��@r�!@r�!@r�!@s@s�@s��@st�@so@q��@q�7@q��@q�@q�#@q��@q�7@qG�@pr�@pb@pb@o�@pb@p �@pb@o�@o�w@o�w@o�@ol�@n�y@m�T@mp�@m?}@m?}@mO�@mO�@mO�@mO�@m?}@m�@l�j@l��@l��@l��@l��@l��@m?}@l��@l�j@l�/@l�/@l�@l��@l��@l�j@l�j@l�j@l�@l�D@lz�@lz�@lj@lj@lI�@k��@k��@k��@k��@k��@k�@kt�@kS�@ko@j�@k@k@k@k@j�H@j��@j�!@j��@j~�@j^5@j-@jJ@jJ@jJ@i�#@i��@i��@ix�@i�7@iX@iG�@i7L@i&�@i%@h��@h��@i%@i�@i%@h��@h�9@h�9@h��@h��@h��@hr�@h �@hb@hb@h �@h �@h �@h �@h �@h �@h �@g�@g�;@g�;@g��@g��@g��@g�w@g��@g|�@gl�@g\)@gK�@g;d@g;d@g+@f�y@fȴ@f�@f��@f�+@fv�@fff@f�+@f��@f�+@fv�@fV@fE�@fE�@f5?@f$�@e�@e�@e��@e@e@e�-@e�-@e�h@ep�@eO�@e?}@e/@eV@d��@d��@d�j@d��@d�D@dz�@dz�@dj@dZ@dI�@dI�@dI�@dj@dj@dz�@d�D@d��@d��@d��@d�D@dZ@d9X@d(�@d�@c�m@cƨ@c�F@c�@cdZ@cS�@c33@c33@c33@c"�@co@c@b�@b�@co@co@c@b�H@b��@b�\@b^5@a�@a��@a�7@ahs@ahs@ax�@ax�@ahs@aG�@a&�@`��@`��@`�u@`�u@`r�@`A�@`1'@`A�@`A�@`A�@_K�@]@\��@[��@\�j@[@ZM�@Y��@WK�@W;d@Z=q@_�@d��@g�;@i��@j��@ix�@f�@ct�@_�@^v�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @��@��@�Ĝ@�@�@� �@�33@�33@�+@�
=@���@���@���@��y@��H@���@�ȴ@���@�R@�+@�=q@��#@�1'@�v�@�  @��@�?}@��!@���@�hs@���@���@�M�@�
=@�9X@�K�@�^5@��@��T@��@�33@�V@�G�@��m@��^@���@���@� �@�ƨ@��@�;d@��@�t�@��
@��m@��m@��m@�  @�1@�b@� �@��@�9X@���@��/@���@��`@���@��w@�;d@��H@���@��y@���@��@��/@�(�@�A�@��
@�hs@���@��7@��D@�r�@�z�@��m@�
=@�^5@�J@���@��-@���@��@�\)@�l�@��@���@��T@��@��@��-@�7L@��`@��`@���@���@��u@� �@�|�@���@��
@���@��R@�"�@��@��+@�$�@��#@���@��T@��@��^@�O�@�O�@�`B@��h@��@�/@���@���@���@���@��^@�`B@�G�@��@���@� �@��;@�+@�~�@�J@�p�@�%@�Ĝ@��@�{@��@��/@���@��@�G�@�X@�O�@�?}@�/@��@�%@���@��@�Ĝ@���@�I�@��m@���@�o@��`@��m@�\)@���@���@�&�@��@�r�@�1@�\)@�o@���@��R@��R@���@��+@�v�@�ff@�n�@�M�@�-@���@���@��-@���@���@���@���@�p�@�  @~E�@~5?@|�@�9X@�1'@� �@� �@�A�@�w@
=@}�T@|�@{�m@{�F@z��@y�@y�7@y%@x��@x�@xQ�@x1'@x �@xb@w�@x  @x1'@yhs@y��@y�@y��@y�@y�7@y&�@x��@x��@x�9@xbN@xb@w��@w|�@w�@vE�@u@t��@t�/@u�@uO�@uO�@t��@u�@u�@u`B@u��@up�@u/@t��@t�/@t�j@t��@tj@sƨ@s33@r�H@r��@r�!@r�!@r�!@s@s�@s��@st�@so@q��@q�7@q��@q�@q�#@q��@q�7@qG�@pr�@pb@pb@o�@pb@p �@pb@o�@o�w@o�w@o�@ol�@n�y@m�T@mp�@m?}@m?}@mO�@mO�@mO�@mO�@m?}@m�@l�j@l��@l��@l��@l��@l��@m?}@l��@l�j@l�/@l�/@l�@l��@l��@l�j@l�j@l�j@l�@l�D@lz�@lz�@lj@lj@lI�@k��@k��@k��@k��@k��@k�@kt�@kS�@ko@j�@k@k@k@k@j�H@j��@j�!@j��@j~�@j^5@j-@jJ@jJ@jJ@i�#@i��@i��@ix�@i�7@iX@iG�@i7L@i&�@i%@h��@h��@i%@i�@i%@h��@h�9@h�9@h��@h��@h��@hr�@h �@hb@hb@h �@h �@h �@h �@h �@h �@h �@g�@g�;@g�;@g��@g��@g��@g�w@g��@g|�@gl�@g\)@gK�@g;d@g;d@g+@f�y@fȴ@f�@f��@f�+@fv�@fff@f�+@f��@f�+@fv�@fV@fE�@fE�@f5?@f$�@e�@e�@e��@e@e@e�-@e�-@e�h@ep�@eO�@e?}@e/@eV@d��@d��@d�j@d��@d�D@dz�@dz�@dj@dZ@dI�@dI�@dI�@dj@dj@dz�@d�D@d��@d��@d��@d�D@dZ@d9X@d(�@d�@c�m@cƨ@c�F@c�@cdZ@cS�@c33@c33@c33@c"�@co@c@b�@b�@co@co@c@b�H@b��@b�\@b^5@a�@a��@a�7@ahs@ahs@ax�@ax�@ahs@aG�@a&�@`��@`��@`�u@`�u@`r�@`A�@`1'@`A�@`A�G�O�@_K�@]@\��@[��@\�j@[@ZM�@Y��@WK�@W;d@Z=q@_�@d��@g�;@i��@j��@ix�@f�@ct�@_�@^v�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
��B
��B
B
ÖB
ÖB
ǮB
��B
��B
��B
��B
��B
��B
��B
��B
��B
�B
��B
��B
��B
��B
ȴB
ĜB
��B,B)�B.B6FBP�BhsBo�Br�By�Bx�B�B�DB�DB�=B�7B�1B�+B�B�B�B�B� B~�B�B�B�B�B�B�B�=B�VB�\B�bB�hB�uB��B��B��B��B��B��B��B��B�B�B��B��B��B��B�B�B��B��B��B��B��B�9B�?B�9B�3B�9B�?B�9B�-B�!B�B�B�B�B��B�B�B�B�B��B�B�B�B�B�B�B�!B�B�!B�B�B�B�!B�B�B�'B�'B�!B�!B�!B�'B�-B�3B�3B�3B�3B�9B�?B�FB�FB�^B�jB�jB�jB�jB�jB�dB�LB�!B��B��B��B��B��B��B��B��B��B��B�'B�FB�FB�LB�XB�dB�dB�jB�jB�jB�jB�jB�jB�qB�wB�wB�qB�dB�XB�-B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�oB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�hB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�DB�DB�DB�JB�JB�JB�JB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B� B� B� B� B� B}�Bz�Bx�Bv�Bz�Bx�Bu�Bt�Bn�Bn�Bw�B�%B��B��B�-B�jB�dB�dB�dB�FB�^111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B
�$B
�$B
�-B
�8B
�7B
�NB
ЃB
ЅB
ыB
ғB
ӗB
ӕB
ӗB
әB
ԛB
գB
ԝB
ԝB
ӖB
ЄB
�RB
�8B
�qB+�B)�B-�B5�BP~BhBo7BrIByuBxoB��B��B��B��B��B��B��B��B��B��B��B�B~�B��B��B��B��B��B��B��B��B��B��B�B�B�$B�:B�GB�DB�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B��B��B��B��B�}B�hB�]B�XB�NB�OB�WB�cB��B��B��B��B��B��B��B��B�B�B�B�B�B�B�
B�B�B�B��B��B��B��B��B��B�|B�uB�pB�hB�^B�PB�KB�BB�DB�BB�KB�QB�PB�QB�VB�VB�WB�ZB�UB�_B�_B�\B�UB�VB�QB�?B�.B�3B�DB�WB�VB�VB�WB�UB�NB�JB�CB�@B�:B�8B�(B�#B�%B�!B�!B�B�B�B�"B�"B�"B�'B�-B�?B�DB�DB�DB�DB�EB�CB�EB�DB�DB�=B�=B�@B�9B�3B�.B�&B� B�%B�.B�.B�+B�.B�0B�2B�:B�@B�AB�CB�:B�7B�;B�9B�9B�4B�+B�,B�,B�+B�+B�+B�1B�?B�?B�?B�9B�,B�'B�*B�3B�1B�*B�+B�&B� B�B�B�B� B� B�%B�'B� B�#B�B�B�B�	B�B� B�
B�B�B�B�B�B�
B�B�B�
B�
B�
B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�B�	B�	B�B�	B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�G�O�B}�BzBxsBveBz~BxtBu]BtZBn8Bn7BwqB��B�?B��B��B�B�B�B�B��B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.74 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051343592019040513435920190405134359  AO  ARCAADJP                                                                    20181006004026    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004026  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004026  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134359  IP                  G�O�G�O�G�O�                