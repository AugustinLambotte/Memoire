CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:32Z creation      
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
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004032  20190405134404  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               QA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�u(�y�G1   @�u)Q���@O$�j~���C;��S��1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    QA   B   B   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A���A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�33B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CS�fCV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�3Dy�3D�fD�<�D�vfD��3D��D�6fD�|�D�� D�fD�C3D�y�D�� D���D�L�Dډ�D�ɚD��D�@ D�3D���1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @j�H@���@أ�AQ�A,Q�ALQ�AlQ�A�(�A�(�A�(�A�(�A�(�A�(�A���A�(�B{B{B{B{B#{B+{B3{B;{BC{BK{BS{B[{Bc{Bk{Bs{B{{B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=B��=BŊ=BɊ=B͊=Bъ=BՊ=Bي=B݊=B�=B�=B�=B�=B�=B��=B��pB��=C �C�C�C�C�C
�C�C�C�C�C�C�C�C�C�C�C �C"�C$�C&�C(�C*�C,�C.�C0�C2�C4�C6�C8�C:�C<�C>�C@�CB�CD�CF�CH�CJ�CL�CN�CP�CR�CT��CV�CX�CZ�C\�C^�C`�Cb�Cd�Cf�Ch�Cj�Cl�Cn�Cp�Cr�Ct�Cv�Cx�Cz�C|�C~�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�C�b�D 1HD �HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD	1HD	�HD
1HD
�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD1HD�HD 1HD �HD!1HD!�HD"1HD"�HD#1HD#�HD$1HD$�HD%1HD%�HD&1HD&�HD'1HD'�HD(1HD(�HD)1HD)�HD*1HD*�HD+1HD+�HD,1HD,�HD-1HD-�HD.1HD.�HD/1HD/�HD01HD0�HD11HD1�HD21HD2�HD31HD3�HD41HD4�HD51HD5�HD61HD6�HD71HD7�HD81HD8�HD91HD9�HD:1HD:�HD;1HD;�HD<1HD<�HD=1HD=�HD>1HD>�HD?1HD?�HD@1HD@�HDA1HDA�HDB1HDB�HDC1HDC�HDD1HDD�HDE1HDE�HDF1HDF�HDG1HDG�HDH1HDH�HDI1HDI�HDJ1HDJ�HDK1HDK�HDL1HDL�HDM1HDM�HDN1HDN�HDO1HDO�HDP1HDP�HDQ1HDQ�HDR1HDR�HDS1HDS�HDT1HDT�HDU1HDU�HDV1HDV�HDW1HDW�HDX1HDX�HDY1HDY�HDZ1HDZ�HD[1HD[�HD\1HD\�HD]1HD]�HD^1HD^�HD_1HD_�HD`1HD`�HDa1HDa�HDb1HDb�HDc1HDc�HDd1HDd�HDe1HDe�HDf1HDf�HDg1HDg�HDh1HDh�HDi1HDi�HDj1HDj�HDk1HDk�HDl1HDl�HDm1HDm�HDn1HDn�HDo1HDo�HDp1HDp�HDq1HDq�HDr1HDr�HDs1HDs�HDt1HDt�HDu${Dy�{D�
D�UqD��
D���D�%qD�O
D��qD��D�
D�[�D��>D�ؤD�qD�eqDڢ>D��>D�%qD�X�D��D��>1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A�7A7LAȴAƨA|�A`BA�A�/A�9A�+AffAI�A(�AJA�;@���@��`@�z�@�33@��;@�^5@�33@�I�@�K�@��@�ȴ@�~�@��j@�ƨ@�b@��
@�o@�{@�O�@��/@�r�@�1'@��@�|�@�"�@��y@��\@�5?@��-@�hs@�V@��@���@���@��-@��@�1'@��@�l�@�
=@�-@��^@��7@�G�@�/@�&�@��`@��D@�Z@��;@�\)@��!@�{@�hs@�Ĝ@�Q�@��m@�33@��@���@���@�J@���@�/@�Ĝ@�j@��m@��F@��@�C�@�;d@�S�@��R@��\@���@��H@�v�@���@��7@��@���@�bN@���@��F@���@�K�@�v�@�-@�J@��@��#@���@�hs@�G�@��@��@�bN@�  @��@�t�@�\)@�;d@�"�@�o@��\@�-@��@�@��h@�hs@�%@��`@��j@�I�@��m@�"�@���@��+@��@���@�bN@�9X@�  @��@��P@�\)@�|�@��@�S�@�"�@�@�$�@��#@��-@�x�@�/@��@���@��@�Z@��;@�|�@�+@��R@�ff@��@���@�X@�A�@��F@��P@�K�@�@��H@���@�E�@�5?@�$�@�@�J@�J@���@��#@���@��-@�x�@��@��D@�j@�I�@��m@��w@���@�S�@�
=@���@�V@�@��-@�p�@�7L@���@�Ĝ@��j@��@�(�@�1@��@���@���@���@��@�t�@�l�@�l�@�|�@��@�C�@�+@��@��@���@���@�ff@�M�@�=q@�$�@��@�{@���@���@��^@��@�p�@�7L@���@���@�Z@���@��F@�t�@�K�@�"�@�
=@���@��R@��!@�~�@�~�@�M�@�=q@�$�@�@��@���@�X@�?}@�7L@�7L@�&�@��/@��9@��@�z�@�z�@�Z@�9X@��@�  @��@�@�P@\)@+@~��@~�y@~�@~��@~E�@~@}��@}�@}/@|�/@|�j@|�@|j@|9X@|�@{ƨ@{��@{t�@{C�@z�@z�H@z�H@z�H@z��@z�!@z~�@zn�@zn�@z^5@z-@y��@yx�@yX@y7L@yG�@y&�@x�`@x��@x�@x�@x�@xr�@xb@w��@w|�@w;d@v�y@v�@v��@v5?@u�-@u�-@u�-@u�h@u�@up�@up�@u`B@uO�@u?}@u�@u�@uV@t�/@t�j@t��@tz�@tj@tI�@t�@s�m@s�F@s��@s�@so@so@r��@r~�@rn�@rM�@r-@r�@rJ@rJ@rJ@q��@q�@q��@q��@q��@qX@q7L@p�`@p�u@pr�@p��@p�@pQ�@p  @o�@o�;@o��@o��@o�w@o�w@o�P@o;d@o+@o;d@o�@nv�@nE�@nE�@nV@n@m�T@m��@m��@m�@m�@m�-@m�h@m�@mp�@m`B@m`B@mO�@m?}@m�@l�@l��@l�@lz�@lZ@l9X@l�@l�@l1@k��@k��@k�m@kƨ@kƨ@k��@k�@kdZ@kS�@kS�@kC�@kC�@kC�@k33@ko@j�@j�H@j�H@j��@j�!@j�!@j�!@j�!@j�!@j��@j��@j~�@jn�@jn�@jn�@j^5@j^5@j=q@j�@i��@i�#@i�^@ix�@iG�@i&�@i&�@i&�@i&�@i&�@i�@i%@h��@h�9@h��@h��@h��@h��@h��@h��@h�u@h�u@hr�@hA�@h �@h  @g�@g�;@g��@g��@g�P@gl�@gK�@g;d@g;d@g;d@g;d@g
=@f��@f�+@f�+@f�+@fv�@fV@fE�@f{@e�T@e�T@e@d�/@c�F@a�^@_�@]�@]/@]�@`r�@b�@d�D@d�@ep�@`Ĝ@_
=@Y��@U��@R^5@NE�@EV@?��@;��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 A�7A7LAȴAƨA|�A`BA�A�/A�9A�+AffAI�A(�AJA�;@���@��`@�z�@�33@��;@�^5@�33@�I�@�K�@��@�ȴ@�~�@��j@�ƨ@�b@��
@�o@�{@�O�@��/@�r�@�1'@��@�|�@�"�@��y@��\@�5?@��-@�hs@�V@��@���@���@��-@��@�1'@��@�l�@�
=@�-@��^@��7@�G�@�/@�&�@��`@��D@�Z@��;@�\)@��!@�{@�hs@�Ĝ@�Q�@��m@�33@��@���@���@�J@���@�/@�Ĝ@�j@��m@��F@��@�C�@�;d@�S�@��R@��\@���@��H@�v�@���@��7@��@���@�bN@���@��F@���@�K�@�v�@�-@�J@��@��#@���@�hs@�G�@��@��@�bN@�  @��@�t�@�\)@�;d@�"�@�o@��\@�-@��@�@��h@�hs@�%@��`@��j@�I�@��m@�"�@���@��+@��@���@�bN@�9X@�  @��@��P@�\)@�|�@��@�S�@�"�@�@�$�@��#@��-@�x�@�/@��@���@��@�Z@��;@�|�@�+@��R@�ff@��@���@�X@�A�@��F@��P@�K�@�@��H@���@�E�@�5?@�$�@�@�J@�J@���@��#@���@��-@�x�@��@��D@�j@�I�@��m@��w@���@�S�@�
=@���@�V@�@��-@�p�@�7L@���@�Ĝ@��j@��@�(�@�1@��@���@���@���@��@�t�@�l�@�l�@�|�@��@�C�@�+@��@��@���@���@�ff@�M�@�=q@�$�@��@�{@���@���@��^@��@�p�@�7L@���@���@�Z@���@��F@�t�@�K�@�"�@�
=@���@��R@��!@�~�@�~�@�M�@�=q@�$�@�@��@���@�X@�?}@�7L@�7L@�&�@��/@��9@��@�z�@�z�@�Z@�9X@��@�  @��@�@�P@\)@+@~��@~�y@~�@~��@~E�@~@}��@}�@}/@|�/@|�j@|�@|j@|9X@|�@{ƨ@{��@{t�@{C�@z�@z�H@z�H@z�H@z��@z�!@z~�@zn�@zn�@z^5@z-@y��@yx�@yX@y7L@yG�@y&�@x�`@x��@x�@x�@x�@xr�@xb@w��@w|�@w;d@v�y@v�@v��@v5?@u�-@u�-@u�-@u�h@u�@up�@up�@u`B@uO�@u?}@u�@u�@uV@t�/@t�j@t��@tz�@tj@tI�@t�@s�m@s�F@s��@s�@so@so@r��@r~�@rn�@rM�@r-@r�@rJ@rJ@rJ@q��@q�@q��@q��@q��@qX@q7L@p�`@p�u@pr�@p��@p�@pQ�@p  @o�@o�;@o��@o��@o�w@o�w@o�P@o;d@o+@o;d@o�@nv�@nE�@nE�@nV@n@m�T@m��@m��@m�@m�@m�-@m�h@m�@mp�@m`B@m`B@mO�@m?}@m�@l�@l��@l�@lz�@lZ@l9X@l�@l�@l1@k��@k��@k�m@kƨ@kƨ@k��@k�@kdZ@kS�@kS�@kC�@kC�@kC�@k33@ko@j�@j�H@j�H@j��@j�!@j�!@j�!@j�!@j�!@j��@j��@j~�@jn�@jn�@jn�@j^5@j^5@j=q@j�@i��@i�#@i�^@ix�@iG�@i&�@i&�@i&�@i&�@i&�@i�@i%@h��@h�9@h��@h��@h��@h��@h��@h��@h�u@h�u@hr�@hA�@h �@h  @g�@g�;@g��@g��@g�P@gl�@gK�@g;d@g;d@g;d@g;d@g
=@f��@f�+@f�+@f�+@fv�@fV@fE�@f{@e�T@e�TG�O�@d�/@c�F@a�^@_�@]�@]/@]�@`r�@b�@d�D@d�@ep�@`Ĝ@_
=@Y��@U��@R^5@NE�@EV@?��@;��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB@�B>wB<jB;dB:^B:^B9XB9XB9XB9XB9XB9XB9XB8RB33BE�BiyB�\B��B��B��B��B�-B�?B�XB�}BĜBĜBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBǮBƨBŢBĜBĜBÖBÖBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBŢBĜBĜBŢBǮBƨBŢBŢBĜBĜBĜBĜBÖBÖBB��B��B��B��B��B��B��B��B��B��B��B��B�}B�}B��B��B��B��B��B��B��B�}B�}B�}B�wB�wB�wB�wB�qB�jB�dB�dB�dB�XB�XB�XB�RB�RB�RB�XB�dB�dB�jB�qB�jB�^B�^B�XB�XB�RB�RB�RB�RB�RB�LB�FB�?B�FB�FB�?B�?B�9B�-B�'B�'B�'B�'B�'B�!B�!B�!B�!B�'B�'B�'B�-B�-B�-B�-B�-B�-B�-B�'B�'B�'B�'B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�VB�VB�=B�+B�B�B�7B�oB��B��B�B�FB�3B�9B�-B�-B�3B�9B��B��B�1111111111111141111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 B@!B>B<	B;B9�B9�B8�B8�B8�B8�B8�B8�B8�B7�G�O�BE>BiB��B�'B�DB�PB��B��B��B��B�B�3B�3B�TB�uBуB�{B�jB�pB�vB�vB�uB�uB�}BуBуB҈B҇BҊB҈B҉BсB�|B�jB�eB�\B�VB�YB�VB�VB�[B�^B�]B�^B�`B�nB�tB�nB�jB�dB�dB�QB�AB�?B�8B�2B�4B�-B�.B�2B�-B�+B�,B�,B�,B�/B�,B�,B�/B�-B�1B�9B�/B�3B�7B�DB�@B�9B�7B�2B�2B�1B�1B�,B�,B�$B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�~B��B�xB��B�B�~B�vB�xB�vB�yB�xB�xB�xB�xB�zB�yB�yB�zB�|B�yB�yB�yB�wB�yB�yB�yB�yB�sB�rB�rB�rB�rB�rB�rB�sB�sB�rB�pB�pB�qB�sB�sB�sB�sB�qB�qB�qB�pB�sB�sB�sB�sB�mB�lB�lB�nB�nB�nB�qB�nB�nB�pB�nB�lB�nB�nB�dB�gB�fB�gB�dB�bB�aB�aB�aB�cB�^B�aB�^B�aB�aB�aB�^B�^B�aB�\B�\B�[B�ZB�ZB�\B�[B�ZB�\B�ZB�[B�TB�VB�TB�UB�NB�OB�OB�NB�KB�OB�KB�NB�OB�KB�IB�IB�GB�JB�JB�IB�CB�BB�IB�IB�BB�BB�BB�BB�?B�BB�>B�>B�<B�>B�>B�>B�>B�9B�6B�6B�6B�6B�0B�0B�0B�0B�0B�6B�6B�4B�5B�7B�5B�5B�5B�0B�3B�0B�3B�.B�*B�*B�*B�*B�*B�*B�*B�*B�)B�*B�*B�*B�*B�%B�%B�%B�%B�'B�%B�#B�#B�%B�%B�%B�%B�#B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�B�
B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B��B��G�O�B��B��B��B��B��B��B��B�	B�EB�B��B��B��B��B��B��B��B��B��B��B��1111111111111141111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.77 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344042019040513440420190405134404  AO  ARCAADJP                                                                    20181006004032    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004032  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004032  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134404  IP                  G�O�G�O�G�O�                