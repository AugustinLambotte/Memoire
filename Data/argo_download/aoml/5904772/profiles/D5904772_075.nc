CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:30Z creation      
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
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004030  20190405134403  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               KA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�f/7���1   @�f0u��@O6fffff�A�n��P1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    KA   B   B   @333@�  @�  A   A   A@  A`  A�  A�  A���A�33A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BW��B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D��D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2y�D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�3Dy��D�  D�@ D���D��fD�fD�33D�p D��fD�3D�@ D�fD�ɚD��D�6fDڐ D��3D���D�P D�fD�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @dz�@���@أ�AQ�A,Q�ALQ�AlQ�A�(�A�(�A���A�\)A�(�A�(�A�(�A�(�B{B{B{B{B#{B+{B3{B;{BC{BK{BS{BZ�Bc{Bk{Bs{B{{B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=BŊ=BɊ=B͊=Bъ=BՊ=Bي=B݊=B�=B�=B�=B�=B�=B��=B��=B��=C �C�C�C�C�C
�C�C�C�C�C�C�C�C�C�C�C �C"�C$�C&�C(�C*�C,�C.�C0�C2�C4�C6�C8�C:�C<�C>�C@�CB�CD�CF�CH�CJ�CL�CN�CP�CR�CT�CV�CX�CZ�C\�C^�C`�Cb�Cd�Cf�Ch�Cj�Cl�Cn�Cp�Cr�Ct�Cv�Cx�Cz�C|�C~�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�D 1HD �HD1HD�HD1HD�HD*�D�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD	1HD	�HD
1HD
�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD 1HD �HD!1HD!�HD"1HD"�HD#1HD#�HD$1HD$�HD%1HD%�HD&1HD&�HD'1HD'�HD(1HD(�HD)1HD)�HD*1HD*�HD+1HD+�HD,1HD,�HD-1HD-�HD.1HD.�HD/1HD/�HD01HD0�HD11HD1�HD21HD2��D31HD3�HD41HD4�HD51HD5�HD61HD6�HD71HD7�HD81HD8�HD91HD9�HD:1HD:�HD;1HD;�HD<1HD<�HD=1HD=�HD>1HD>�HD?1HD?�HD@1HD@�HDA1HDA�HDB1HDB�HDC1HDC�HDD1HDD�HDE1HDE�HDF1HDF�HDG1HDG�HDH1HDH�HDI1HDI�HDJ1HDJ�HDK1HDK�HDL1HDL�HDM1HDM�HDN1HDN�HDO1HDO�HDP1HDP�HDQ1HDQ�HDR1HDR�HDS1HDS�HDT1HDT�HDU1HDU�HDV1HDV�HDW1HDW�HDX1HDX�HDY1HDY�HDZ1HDZ�HD[1HD[�HD\1HD\�HD]1HD]�HD^1HD^�HD_1HD_�HD`1HD`�HDa1HDa�HDb1HDb�HDc1HDc�HDd1HDd�HDe1HDe�HDf1HDf�HDg1HDg�HDh1HDh�HDi1HDi�HDj1HDj�HDk1HDk�HDl1HDl�HDm1HDm�HDn1HDn�HDo1HDo�HDp1HDp�HDq1HDq�HDr1HDr�HDs1HDs�HDt1HDt�HDu{Dy�D��D�X�D��qD��
D�
D�K�D���D��
D��D�X�D�/
D��>D�%qD�O
Dڨ�D���D�qD�h�D�
D�Ȥ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�V@��@��@��@��@��@��@�&�@�&�@�&�@�&�@��@��@�&�@��@�&�@��@��@�V@��@���@�Ĝ@��9@�z�@�Z@�Z@�r�@��P@��@�|�@�t�@�l�@�S�@�33@�o@��@��@��R@�v�@�V@�-@�G�@�Ĝ@��u@�z�@�bN@�9X@�b@���@�C�@��@�hs@��@�/@�l�@�M�@�7L@��@�~�@��@�I�@���@���@���@��@�?}@���@�  @�ƨ@�;d@��H@���@�{@���@�G�@��/@��D@��@�l�@�l�@�"�@��y@���@�^5@��^@�hs@�G�@�V@��@�j@�A�@���@��@�"�@��@�ȴ@���@��+@�M�@��@��@���@��h@�X@��@���@��9@���@��u@��D@�j@�A�@�1@�S�@��@��H@�ȴ@���@���@�ff@�^5@�v�@�5?@��@��@��T@���@��^@��@�`B@�G�@�7L@�V@��`@��j@��u@�bN@�9X@�(�@�1@�@l�@K�@;d@~��@~ȴ@~ȴ@~ȴ@~ȴ@~E�@}�@}��@}@}�-@}�h@}p�@}O�@}?}@}�@}V@|�@|�/@|��@|�j@|�D@|z�@|�D@|z�@|�@|1@{��@{�
@{��@{�@{dZ@{C�@{S�@{"�@z�\@zM�@z�@y�#@y��@y��@y�7@y�^@z=q@z=q@zJ@y�@y��@y��@y�^@y��@y��@y��@y��@y�7@yX@x�`@x��@x�@xbN@xA�@w��@w|�@wl�@w\)@wl�@wl�@wl�@w;d@v�y@vȴ@vȴ@v��@vv�@vv�@vff@vE�@v@u�T@u@up�@u?}@uO�@up�@up�@up�@u`B@u/@u�@t��@tz�@tI�@t1@sƨ@st�@sdZ@s33@s"�@r�H@r�!@r�\@r^5@r-@rJ@rJ@q��@q��@q��@qx�@q7L@q�@p��@p��@p�u@pbN@pbN@pQ�@pQ�@pQ�@p1'@o�@o�@o�@o��@o�w@o�P@ol�@oK�@o+@n�R@nV@n$�@n@m��@m/@mV@mV@l��@l�/@l�@lj@lj@l�@k"�@j��@j�@j�H@j��@j�@k"�@kC�@kS�@kdZ@kdZ@kS�@k33@ko@kC�@kC�@k"�@j�@j��@j�H@k@k@k@k@ko@k@j�\@j^5@j^5@j~�@j��@j��@j�!@j��@j~�@j=q@j=q@j�@i��@i��@ix�@ihs@ix�@i��@j�@j=q@jM�@j�@i��@i��@i��@i�7@iG�@i7L@i�@i�@iX@ihs@i�^@i�@i��@ix�@ihs@i7L@ix�@ihs@iX@iX@iG�@i7L@i7L@i&�@i7L@i�@i%@h��@h��@h��@hĜ@h��@h�u@hr�@hQ�@h1'@h1'@h1'@h1'@h �@hb@h  @g�;@g�;@g��@g��@g��@g�P@gl�@g\)@gK�@g+@g�@f��@f�y@f�y@fȴ@f��@f�+@f�+@f�+@fv�@fff@fff@fE�@f5?@f5?@f5?@fV@fff@f�+@f�+@f�+@f�+@fv�@fv�@fv�@fv�@fv�@fE�@f$�@f$�@f$�@e�@e�T@e�T@e��@e@e@e@e��@e�@e�@e�@ep�@eO�@e/@e?}@e?}@e?}@d�/@d��@d�D@d�D@d�D@dj@dZ@d9X@dI�@d9X@d(�@d�@d1@c��@c�
@cƨ@cƨ@cƨ@c�
@c�
@cƨ@c��@c�@cS�@c"�@c@b�H@b��@b��@b�\@b~�@b~�@b~�@b�\@b��@b��@b��@b�\@b�\@bM�@bJ@a�^@a�^@a�^@a�^@a�^@a�^@a��@a�^@a��@a��@a�7@ax�@aX@`�`@_�w@^�y@]�h@\��@\9X@\j@]��@a��@e�T@g+@e@hQ�@f�y@d��@a��@^�y@\�j@Z�\@Vff@R��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 @��@��@�V@��@��@��@��@��@��@�&�@�&�@�&�@�&�@��@��@�&�@��@�&�@��@��@�V@��@���@�Ĝ@��9@�z�@�Z@�Z@�r�@��P@��@�|�@�t�@�l�@�S�@�33@�o@��@��@��R@�v�@�V@�-@�G�@�Ĝ@��u@�z�@�bN@�9X@�b@���@�C�@��@�hs@��@�/@�l�@�M�@�7L@��@�~�@��@�I�@���@���@���@��@�?}@���@�  @�ƨ@�;d@��H@���@�{@���@�G�@��/@��D@��@�l�@�l�@�"�@��y@���@�^5@��^@�hs@�G�@�V@��@�j@�A�@���@��@�"�@��@�ȴ@���@��+@�M�@��@��@���@��h@�X@��@���@��9@���@��u@��D@�j@�A�@�1@�S�@��@��H@�ȴ@���@���@�ff@�^5@�v�@�5?@��@��@��T@���@��^@��@�`B@�G�@�7L@�V@��`@��j@��u@�bN@�9X@�(�@�1@�@l�@K�@;d@~��@~ȴ@~ȴ@~ȴ@~ȴ@~E�@}�@}��@}@}�-@}�h@}p�@}O�@}?}@}�@}V@|�@|�/@|��@|�j@|�D@|z�@|�D@|z�@|�@|1@{��@{�
@{��@{�@{dZ@{C�@{S�@{"�@z�\@zM�@z�@y�#@y��@y��@y�7@y�^@z=q@z=q@zJ@y�@y��@y��@y�^@y��@y��@y��@y��@y�7@yX@x�`@x��@x�@xbN@xA�@w��@w|�@wl�@w\)@wl�@wl�@wl�@w;d@v�y@vȴ@vȴ@v��@vv�@vv�@vff@vE�@v@u�T@u@up�@u?}@uO�@up�@up�@up�@u`B@u/@u�@t��@tz�@tI�@t1@sƨ@st�@sdZ@s33@s"�@r�H@r�!@r�\@r^5@r-@rJ@rJ@q��@q��@q��@qx�@q7L@q�@p��@p��@p�u@pbN@pbN@pQ�@pQ�@pQ�@p1'@o�@o�@o�@o��@o�w@o�P@ol�@oK�@o+@n�R@nV@n$�@n@m��@m/@mV@mV@l��@l�/@l�@lj@lj@l�@k"�@j��@j�@j�H@j��@j�@k"�@kC�@kS�@kdZ@kdZ@kS�@k33@ko@kC�@kC�@k"�@j�@j��@j�H@k@k@k@k@ko@k@j�\@j^5@j^5@j~�@j��@j��@j�!@j��@j~�@j=q@j=q@j�@i��@i��@ix�@ihs@ix�@i��@j�@j=q@jM�@j�@i��@i��@i��@i�7@iG�@i7L@i�@i�@iX@ihs@i�^@i�@i��@ix�@ihs@i7L@ix�@ihs@iX@iX@iG�@i7L@i7L@i&�@i7L@i�@i%@h��@h��@h��@hĜ@h��@h�u@hr�@hQ�@h1'@h1'@h1'@h1'@h �@hb@h  @g�;@g�;@g��@g��@g��@g�P@gl�@g\)@gK�@g+@g�@f��@f�y@f�y@fȴ@f��@f�+@f�+@f�+@fv�@fff@fff@fE�@f5?@f5?@f5?@fV@fff@f�+@f�+@f�+@f�+@fv�@fv�@fv�@fv�@fv�@fE�@f$�@f$�@f$�@e�@e�T@e�T@e��@e@e@e@e��@e�@e�@e�@ep�@eO�@e/@e?}@e?}@e?}@d�/@d��@d�D@d�D@d�D@dj@dZ@d9X@dI�@d9X@d(�@d�@d1@c��@c�
@cƨ@cƨ@cƨ@c�
@c�
@cƨ@c��@c�@cS�@c"�@c@b�H@b��@b��@b�\@b~�@b~�@b~�@b�\@b��@b��@b��@b�\@b�\@bM�@bJ@a�^@a�^@a�^@a�^@a�^@a�^@a��@a�^@a��@a��@a�7@ax�G�O�@`�`@_�w@^�y@]�h@\��@\9X@\j@]��@a��@e�T@g+@e@hQ�@f�y@d��@a��@^�y@\�j@Z�\@Vff@R��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�uB�uB�uB�uB�{B�{B�{B��B��B��B��B��B��B��B��B��B�{B�{B�{B��B��B��B��B��B��B�{B�{B�{B�{B��B��B��B��B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�{B��B��B��B�{B�{B�{B�{B�uB�uB�uB�uB�{B�{B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�VB�VB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�DB�DB�DB�DB�DB�JB�JB�JB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�DB�7B�%B�%B�B�B�B�B�1B�{B��B�B�B�XB�dB�jB�qB�qB�}B��B�qB�d1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}B�xB�xB�pB�lB�lB�dB�gB�_B�aB�_B�UB�`B�^B�`B�`B�^B�aB�^B�_B�^B�\B�WB�QB�RB�SB�RB�RB�NB�QB�PB�PB�YB�YB�YB�YB�WB�SB�ZB�RB�SB�RB�SB�RB�RB�KB�GB�FB�DB�FB�FB�DB�GB�FB�DB�JB�GB�EB�GB�GB�EB�GB�EB�EB�EB�@B�AB�@B�@B�7B�8B�8B�:B�8B�8B�8B�8B�8B�8B�8B�:B�8B�<B�;B�;B�9B�9B�;B�9B�9B�;B�9B�;B�9B�9B�;B�;B�BB�?B�@B�BB�?B�?B�@B�@B�>B�?B�EB�HB�LB�LB�EB�EB�NB�LB�GB�FB�EB�GB�UB�WB�WB�VB�VB�WB�_B�_B�eB�cB�aB�eB�eB�`B�`B�`B�dB�`B�`B�`B�`B�`B�`B�`B�`B�YB�[B�ZB�ZB�YB�YB�[B�[B�ZB�[B�ZB�YB�RB�RB�RB�ZB�ZB�ZB�\B�ZB�ZB�RB�RB�TB�MB�MB�OB�MB�MB�GB�HB�GB�GB�FB�DB�DB�DB�AB�BB�BB�AB�?B�?B�;B�:B�:B�;B�;B�:B�<B�<B�:B�:B�:B�;B�:B�;B�>B�;B�4B�6B�6B�-B�0B�.B�)B�"B� B�$B�"B�$B�B�B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�	B�	B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B�B� B� B�B� B�B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B� B��B�B� B�B�B� B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B��B��B��B�B�sB��B��B��B�B�B�B�B�B�)B�B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.77 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344032019040513440320190405134403  AO  ARCAADJP                                                                    20181006004030    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004030  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004030  QCF$                G�O�G�O�G�O�48000           UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134403  IP                  G�O�G�O�G�O�                