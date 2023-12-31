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
_FillValue                    ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �D   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �T   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �X   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �h   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �l   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �p   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004026  20190405134400  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               :A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�;����1   @�;�wwy�@N�/��w�A��G�{1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    :A   B   B   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B7��B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�3Dy��D��fD�33D���D�� D���D�0 D�� D�� D��D�0 D��3DǶfD��3D�I�DږfD��D��fD�9�D�fD��f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @j�H@���@أ�AQ�A,Q�ALQ�AlQ�A�(�A�(�A�(�A�(�A�(�A�(�A�(�A�(�B{B{B{B{B#{B+{B3{B:�BC{BK{BS{B[{Bc{Bk{Bs{B{{B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=BŊ=BɊ=B͊=Bъ=BՊ=Bي=B݊=B�=B�=B�=B�=B�=B��=B��=B��=C ��C�C�C�C�C
�C�C�C�C�C�C�C�C�C�C�C �C"�C$�C&�C(�C*�C,�C.�C0�C2�C4�C6�C8�C:�C<�C>�C@�CB�CD�CF�CH�CJ�CL�CN�CP�CR�CT�CV�CX�CZ�C\�C^�C`�Cb�Cd�Cf�Ch�Cj�Cl�Cn�Cp�Cr�Ct�Cv�Cx�Cz�C|�C~�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�D 1HD �HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD	1HD	�HD
1HD
�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD 1HD �HD!1HD!�HD"1HD"�HD#1HD#�HD$1HD$�HD%1HD%�HD&1HD&�HD'1HD'�HD(1HD(�HD)1HD)�HD*1HD*�HD+1HD+�HD,1HD,�HD-1HD-�HD.1HD.�HD/1HD/�HD01HD0�HD11HD1�HD21HD2�HD31HD3�HD41HD4�HD51HD5�HD61HD6�HD71HD7�HD81HD8�HD91HD9�HD:1HD:�HD;1HD;�HD<1HD<�HD=1HD=�HD>1HD>�HD?1HD?�HD@1HD@�HDA1HDA�HDB1HDB�HDC1HDC�HDD1HDD�HDE1HDE�HDF1HDF�HDG1HDG�HDH1HDH�HDI1HDI�HDJ1HDJ�HDK1HDK�HDL1HDL�HDM1HDM�HDN1HDN�HDO1HDO�HDP1HDP�HDQ1HDQ�HDR1HDR�HDS1HDS�HDT1HDT�HDU1HDU�HDV1HDV�HDW1HDW�HDX1HDX�HDY1HDY�HDZ1HDZ�HD[1HD[�HD\1HD\�HD]1HD]�HD^1HD^�HD_1HD_�HD`1HD`�HDa1HDa�HDb1HDb�HDc1HDc�HDd1HDd�HDe1HDe�HDf1HDf�HDg1HDg�HDh1HDh�HDi1HDi�HDj1HDj�HDk1HDk�HDl1HDl�HDm1HDm�HDn1HDn�HDo1HDo�HDp1HDp�HDq1HDq�HDr1HDr�HDs1HDs�HDt1HDt�HDu{Dy�D�
D�K�D��qD��D�qD�H�D���D��D�%qD�H�D���D��
D���D�b>Dگ
D��qD�
D�R>D�
D��
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�\)@�\)@�\)@�dZ@�dZ@�dZ@�dZ@�dZ@�l�@�l�@�l�@�l�@�l�@�l�@�l�@�t�@�t�@�t�@�t�@�t�@�t�@�t�@�|�@�|�@�|�@�|�@��@��@��@��@��@��@��P@��P@��P@��P@��P@��P@��P@��P@���@���@��P@��P@���@���@���@���@���@���@���@���@���@���@���@���@��@�^5@�J@��@���@���@�V@��9@��u@�t�@�o@�M�@�/@���@�1'@�1'@�ƨ@�dZ@���@�@���@�p�@�`B@�`B@�O�@�O�@�&�@��@��u@�Q�@�A�@� �@�  @�ƨ@���@�l�@�33@��@�~�@�@���@�x�@��7@�`B@���@���@��j@��@���@�j@�  @�b@��F@��@��R@���@�J@��h@��@���@���@�j@�1'@�(�@�ƨ@�33@���@���@�v�@�M�@�J@���@�p�@�X@��@�Ĝ@���@�z�@�1'@� �@�b@�  @��@�b@�1@�  @�@�w@�P@K�@~�@~ȴ@~��@~��@~�+@~$�@~@}��@}��@}p�@|��@}?}@}V@|��@|��@|��@|�j@|��@|z�@|I�@|(�@|�@|1@{�m@{��@{C�@z�H@z��@z�!@z�!@z�H@z��@z�!@z�\@z=q@y�^@y�@y%@x�`@xbN@xA�@x �@w�@w|�@w\)@w
=@v�@v��@v��@v�+@vff@v5?@u��@u@u@u@u@u@u@u��@u�-@u�@up�@t�/@tj@s�m@s�m@sƨ@s�F@s�@sS�@sC�@sC�@s"�@s@r�!@rn�@rn�@r^5@r-@q��@q�#@qhs@qG�@q7L@qG�@q7L@q&�@q&�@q7L@q&�@q7L@q7L@q�@p�`@pr�@p1'@p �@p  @o�@o�@o�@o��@o��@oK�@n�y@nv�@nV@n@m�h@l��@l��@l��@l��@l��@l��@m?}@m�@m@m�h@mO�@l�/@l��@l�j@l�@l�j@l�j@l�/@l�/@l��@lj@lI�@l9X@l(�@l(�@l(�@l9X@lI�@l9X@l9X@l9X@l9X@l9X@k�
@kdZ@j�@k@ko@k"�@ko@k@k@j��@j��@j~�@j~�@j~�@j^5@jM�@jM�@j=q@j-@i��@i�^@i��@i�7@ix�@ihs@iG�@i7L@i�@h�`@h��@h��@hQ�@g�@g�@g�@g�P@g;d@g+@g
=@f��@g+@g�w@h  @g��@g�@g�P@gK�@g
=@g�@gK�@gl�@g�P@g�;@h  @h  @g�w@g+@f�y@f�@fȴ@f�@fȴ@f�R@fȴ@fff@fV@fE�@fff@fff@f5?@f$�@f{@f{@fE�@f$�@f{@f{@f$�@f5?@fE�@f5?@f{@e��@e@e�T@f$�@fE�@fE�@f$�@e�@f@f@e�T@e��@e@e�-@e��@e�-@e��@e�@eO�@e/@e/@e?}@e?}@e/@e�@d�@d��@d�j@d�j@d�j@d�@d�@d�@d�@d�@d��@d�D@dj@dZ@dI�@d9X@d9X@d(�@d�@d�@d�@d�@d1@c��@c�
@c�
@c�
@cƨ@c��@cS�@cS�@cS�@cC�@c33@co@b�@b�H@b��@b��@b�!@b�!@b�!@b��@b�!@b��@b��@b��@b��@b��@b�@b�@b�H@b��@b�!@b�!@b=q@a�^@a��@a�7@ax�@ax�@b~�@bn�@a��@a�^@a�^@bM�@b~�@b�!@b�H@co@c"�@c"�@c@b�H@b�!@b�\@bM�@b=q@b=q@bM�@b�\@c"�@c�@cC�@c33@ct�@ct�@cdZ@c�@ct�@ct�@ct�@cdZ@cdZ@`Ĝ@bJ@`A�@^v�@_+@_�@]@a7L@a��@d�D@gK�@hbN@g�;@g��@e��@b�@`�`@\��@Y�#@Xr�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 @�\)@�\)@�\)@�dZ@�dZ@�dZ@�dZ@�dZ@�l�@�l�@�l�@�l�@�l�@�l�@�l�@�t�@�t�@�t�@�t�@�t�@�t�@�t�@�|�@�|�@�|�@�|�@��@��@��@��@��@��@��P@��P@��P@��P@��P@��P@��P@��P@���@���@��P@��P@���@���@���@���@���@���@���@���@���@���@���@���@��@�^5@�J@��@���@���@�V@��9@��u@�t�@�o@�M�@�/@���@�1'@�1'@�ƨ@�dZ@���@�@���@�p�@�`B@�`B@�O�@�O�@�&�@��@��u@�Q�@�A�@� �@�  @�ƨ@���@�l�@�33@��@�~�@�@���@�x�@��7@�`B@���@���@��j@��@���@�j@�  @�b@��F@��@��R@���@�J@��h@��@���@���@�j@�1'@�(�@�ƨ@�33@���@���@�v�@�M�@�J@���@�p�@�X@��@�Ĝ@���@�z�@�1'@� �@�b@�  @��@�b@�1@�  @�@�w@�P@K�@~�@~ȴ@~��@~��@~�+@~$�@~@}��@}��@}p�@|��@}?}@}V@|��@|��@|��@|�j@|��@|z�@|I�@|(�@|�@|1@{�m@{��@{C�@z�H@z��@z�!@z�!@z�H@z��@z�!@z�\@z=q@y�^@y�@y%@x�`@xbN@xA�@x �@w�@w|�@w\)@w
=@v�@v��@v��@v�+@vff@v5?@u��@u@u@u@u@u@u@u��@u�-@u�@up�@t�/@tj@s�m@s�m@sƨ@s�F@s�@sS�@sC�@sC�@s"�@s@r�!@rn�@rn�@r^5@r-@q��@q�#@qhs@qG�@q7L@qG�@q7L@q&�@q&�@q7L@q&�@q7L@q7L@q�@p�`@pr�@p1'@p �@p  @o�@o�@o�@o��@o��@oK�@n�y@nv�@nV@n@m�h@l��@l��@l��@l��@l��@l��@m?}@m�@m@m�h@mO�@l�/@l��@l�j@l�@l�j@l�j@l�/@l�/@l��@lj@lI�@l9X@l(�@l(�@l(�@l9X@lI�@l9X@l9X@l9X@l9X@l9X@k�
@kdZ@j�@k@ko@k"�@ko@k@k@j��@j��@j~�@j~�@j~�@j^5@jM�@jM�@j=q@j-@i��@i�^@i��@i�7@ix�@ihs@iG�@i7L@i�@h�`@h��@h��@hQ�@g�@g�@g�@g�P@g;d@g+@g
=@f��@g+@g�w@h  @g��@g�@g�P@gK�@g
=@g�@gK�@gl�@g�P@g�;@h  @h  @g�w@g+@f�y@f�@fȴ@f�@fȴ@f�R@fȴ@fff@fV@fE�@fff@fff@f5?@f$�@f{@f{@fE�@f$�@f{@f{@f$�@f5?@fE�@f5?@f{@e��@e@e�T@f$�@fE�@fE�@f$�@e�@f@f@e�T@e��@e@e�-@e��@e�-@e��@e�@eO�@e/@e/@e?}@e?}@e/@e�@d�@d��@d�j@d�j@d�j@d�@d�@d�@d�@d�@d��@d�D@dj@dZ@dI�@d9X@d9X@d(�@d�@d�@d�@d�@d1@c��@c�
@c�
@c�
@cƨ@c��@cS�@cS�@cS�@cC�@c33@co@b�@b�H@b��@b��@b�!@b�!@b�!@b��@b�!@b��@b��@b��@b��@b��@b�@b�@b�H@b��@b�!@b�!@b=q@a�^@a��@a�7@ax�@ax�@b~�@bn�@a��@a�^@a�^@bM�@b~�@b�!@b�H@co@c"�@c"�@c@b�H@b�!@b�\@bM�@b=q@b=q@bM�@b�\@c"�@c�@cC�@c33@ct�@ct�@cdZ@c�@ct�@ct�@ct�G�O�@cdZ@`Ĝ@bJ@`A�@^v�@_+@_�@]@a7L@a��@d�D@gK�@hbN@g�;@g��@e��@b�@`�`@\��@Y�#@Xr�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB<jB=qBP�BffBm�Bp�Bs�Bv�B}�B�B�B�VB�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�oB�oB�uB�uB�uB�uB�{B�{B��B�{B�{B�{B�uB�uB�uB�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�uB�oB�oB�oB�oB�oB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�\B�\B�VB�VB�PB�PB�PB�JB�JB�JB�JB�PB�VB�\B�VB�VB�VB�PB�PB�PB�PB�VB�VB�\B�\B�\B�VB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�+B�+B�+B�+B�+B�+B�%B�%B�B�B�B�B�B�B�%B�%B�B�B�B�+B�+B�+B�1B�1B�7B�1B�1B�1B�1B�+B�+B�+B�+B�+B�1B�7B�=B�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�%B�=B�+B�B�+B�+B�+B�oB��B��B�B�?B�RB�dB�dB�dB�XB�RB�XB�^1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B< B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B;�B< B;�B;�B;�B< B;�B;�B;�B;�B;�B<B;�B;�B;�B;�B;�B;�B;�B;�B=BPyBe�Bm%Bp8BsLBv[B}�B��B��B��B��B� B�JB�JB�^B�fB�nB�nB�tB�}B�{B�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�wB�wB�wB�{B��B��B��B�{B�zB�yB�{B�{B�zB��B�{B�xB�oB�oB�jB�dB�_B�_B�^B�_B�_B�XB�XB�YB�YB�YB�YB�QB�PB�RB�JB�IB�EB�AB�>B�AB�AB�AB�@B�=B�?B�KB�JB�KB�JB�JB�JB�JB�IB�DB�DB�JB�GB�IB�KB�KB�MB�EB�EB�WB�XB�VB�XB�WB�WB�VB�XB�WB�VB�WB�WB�UB�SB�QB�TB�TB�PB�TB�RB�RB�PB�TB�TB�LB�FB�FB�FB�FB�?B�AB�@B�8B�9B�9B�3B�3B�3B�3B�6B�3B�3B�6B�3B�3B�3B�3B�3B�1B�3B�3B�3B�4B�-B�-B�2B�4B�4B�4B�2B�+B�.B�,B�.B�,B�.B�'B�'B�)B�&B�(B�'B�'B�'B�)B�"B�(B�(B�'B�'B�1B�.B�+B�)B�'B�)B�)B�'B�$B�'B�!B�!B�!B� B�B�B�B�B�	B�B�B�
B�B�
B�
B�B�B�B�B�B�B�	B�	B�
B�B�B�B�B�B�B�
B�
B�
B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B��B��B��B�	B�B�iB��B��B��B�B�B�B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.77 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344002019040513440020190405134400  AO  ARCAADJP                                                                    20181006004026    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004026  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004026  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134400  IP                  G�O�G�O�G�O�                