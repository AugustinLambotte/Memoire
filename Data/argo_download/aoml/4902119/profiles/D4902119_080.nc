CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92018-11-06T01:28:57Z creation; 2023-02-09T14:06:17Z DMQC;      
references        (http://www.argodatamgt.org/Documentation   comment       	free text      user_manual_version       3.2    Conventions       Argo-3.2 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         iPRIMARY | https://orcid.org/0000-0001-5113-1068 | Deborah West-Mack, Woods Hole Oceanographic Institution         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7d   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7t   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7x   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7|   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  �  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  �  8<   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  `  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9$   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9(   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                  @  9,   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9l   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9t   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                  @  9x   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                  @  9�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                  @  9�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    :8   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~       axis      T           :@   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    :P   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~            :T   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           :d   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           :t   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    :�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    :�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    :�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    <�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    <�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  SL   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  X�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  u<   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �|   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �$   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �l   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �T   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 � �   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     � �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` *<   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   *�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   0�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   6�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T <�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   <�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   <�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   =    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   =   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � =   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   =�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   =�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    =�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        =�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        =�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       =�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    =�Argo profile    3.1 1.2 19500101000000  20181106012857  20230209090617  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               P   PAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�y)T2$t@�y)T2$t11  @�y)`��@�y)`��@N�Nz�VC@N�Nz�VC�;�f�P��;�f�P�11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�z�@�\@@  @}p�@��R@��R@޸R@�p�A�RA#33A@��A`  A�  A�Q�A��A�  A�Q�A�Q�A��A�Q�B (�B  B�
B  B   B(  B0(�B8  B@  BG�BP(�BX  B`  Bh  Bo�
Bx  B�{B�(�B�(�B�  B��
B��B��B��B�{B�{B��B��B�  B��B��
B�B��
B��B�{B�  B��
B��B�{B�  B��B��B��B�(�B�{B�{B�  B�  C   C��C��C��C��C	��C{C
=C��C  C{C
=C  C
=C{C
=C��C"  C$
=C&{C({C*{C,
=C.
=C0
=C2  C3�C5�HC7��C:
=C<  C>  C?��CB  CD{CF  CH  CJ  CL
=CN
=CO��CR  CT
=CU��CW�HCY��C\
=C]��C_�Cb
=Cd{Cf  Cg�Cj  Cl
=Cn  Cp  Cq��Cs��Cu��Cx  Cz{C|  C~  C�\C�\C�  C���C���C���C���C��C���C�  C�  C���C�  C���C�  C�\C�C�  C���C�  C�
=C���C��C�  C�
=C�C�  C���C�  C�  C���C���C�C�  C���C���C�C�  C�  C�
=C�
=C�C�C�
=C�
=C�
=C���C��C���C���C�  C�  C�  C�C�
=C�  C���C�  C�C���C��C��C���C�C�\C�C���C���C�  C�C�
=C�  C��C���C�  C�C�C�C�C�  C���C���C���C���C���C��C��C��C��C��C��C�  C�  C���C�
=C�
=C�C�  C���C�  C�C�  C���C���C�  C�
=C���C���C�  C���C�  C�C���C��C���C�
=C�  C���C���C���C�C�C�C�  C���C���C���C���D �D � D  D}qD  D}qD��DxRD�RDxRD  D��D�D�DD��D  D}qD	  D	��D
D
� D
�qD� DD��D  D� D��D� DD� D�D��D  D� D�D� D  D� D  D��D�qD}qD�qDz�D�qD}qD  D��D  D}qD  D� D�qD� D�D��D�D� D�D� D�qD}qD   D ��D!  D!}qD!��D"xRD"�qD#�D$�D$� D$�qD%z�D%��D&}qD'�D'}qD'�RD(}qD)�D)��D*�D*z�D*�qD+��D+�qD,z�D,��D-� D.�D.�D/D/� D0  D0� D1  D1��D2�D2� D3  D3� D4  D4� D4�qD5}qD5��D6}qD7  D7z�D7�qD8}qD8��D9}qD:�D:� D;�D;��D<�D<�D<�qD=z�D>�D>�D?�D?}qD@�D@��D@�qDA}qDB  DBxRDC  DC��DD  DDz�DD�qDE��DF  DF� DF��DGz�DG�RDH}qDI  DI� DI�qDJz�DJ�qDK�DL�DL}qDM�DM�DN  DN� DODO� DP�DP�DQ�DQ� DQ��DR}qDS�DS}qDT  DT�DU�DU��DV  DV� DW  DW}qDW�RDX� DYDY� DZ  DZ� D[  D[�D[�qD\xRD]  D]��D^�D^� D^�qD_z�D_��D`z�D`�qDaz�Da�qDb� Dc�Dc�Dd�Dd}qDe  De� Df  Df��DgDg}qDg�RDh}qDh�qDi� Dj�Dj�Dk  Dkz�Dl  Dl� Dm  Dm}qDn  Dn��Do  Do}qDp  Dp��Dq  Dq� DrDr��Ds�Ds� Dt  Dt}qDt��Du��Dv  Dvz�Dv�RDwz�Dx�Dx}qDy�Dy� Dy�qDz� Dz�qD{� D|  D|� D}�D}��D~�D~z�D~�qD}qD�qD�>�D�~�D���D��qD�@ D�� D�� D��D�@ D�}qD��qD��qD�>�D�� D��HD��D�C�D��HD���D��qD�@ D��HD�D�  D�>�D�~�D��qD�HD�>�D�|)D���D�  D�B�D�� D���D���D�>�D�}qD��qD��)D�<)D�� D�D�HD�@ D�� D�� D�  D�>�D�� D�� D���D�@ D�� D��HD�  D�>�D�}qD��qD�HD�B�D���D�D���D�>�D��HD�� D�  D�@ D�~�D��HD�  D�=qD�~�D�� D���D�@ D��HD���D���D�AHD�� D��qD�  D�@ D�~�D��HD���D�>�D�}qD��qD���D�AHD��HD��HD��D�B�D���D�� D��qD�=qD�|)D��qD��qD�>�D�~�D���D��qD�>�D�~�D�� D�  D�>�D�� D��HD�HD�B�D��HD���D�HD�>�D�~�D���D��qD�AHD��HD���D�HD�@ D�� D�� D�  D�>�D�}qD���D���D�@ D�}qD���D�HD�@ D�� D��HD��qD�>�D�� D���D���D�@ D�� D�� D���D�>�D�~�D���D���D�@ D�� D�� D�  D�=qD�~�D�� D�  D�B�D���D�� D��qD�>�D��HD�� D��qD�>�D�� D�� D�  D�@ D�� D�� D�HD�AHD��HD���D���D�AHD�� D�� D�  D�@ D�~�D��qD���D�=qD�� D�� D���D�>�D�~�D���D���D�@ D�� D��HD�HD�>�D�� D�� ?L��?u?�=q?��R?�{?�p�?���?�(�?��@�\@��@z�@(�@#�
@.{@5@@  @J=q@Q�@^�R@c�
@k�@s33@}p�@��
@���@�{@��@�@��H@��R@��
@��@��@��@�@���@��R@��
@Ǯ@���@У�@�@��H@޸R@�\@�@�@��@�@���@�p�A ��A33AAQ�A
�HA(�A�RA��A33AAQ�A=qA(�A{A ��A#33A%A(Q�A*=qA,(�A.�RA0��A333A5A8Q�A:�HA=p�A?\)AA�AC�
AEAH��AJ�HAN{AP  AQ�ATz�AVffAX��AZ�HA^{A`��Ab�\Ae�Ag
=Ah��Ak�An{Ap��As33Au�Aw
=Ax��A|(�A~{A�Q�A�G�A�=qA��A���A�{A�
=A�  A�G�A��\A��A�z�A�A�\)A�Q�A�G�A��\A��A��A�{A�
=A�  A�G�A��\A��A���A�{A�
=A�Q�A���A��\A��A���A�{A�
=A�  A�G�A��\A��
A���A�{A�
=A�  A�G�A��\A��A�z�A�A�
=A�  A���A�=qA�33A�z�A�p�A��RA�  A���A��A�33A�z�A�A�ffAǮA���A��A��HA�(�A�p�A�ffA�\)A�Q�A��Aҏ\AӅA���A�A�
=A�  A�G�A�=qA�33A�(�A�p�A�ffA�\)A��AᙚA�\A�A���A�A�RA�A���A��A��HA��
A��A�{A�
=A�  A�G�A�=qA��HA��
A��A�{A��RA�  A���A��A��\A��A�z�A�p�A�ffA�
=B   B z�B ��BG�BB=qB�\B�HB\)B�
B(�Bz�B��Bp�BB=qBffB�HB\)B�B  Bz�B��B	�B	p�B	�B
ffB
�RB
=B�B�
BQ�Bz�B��Bp�BB{BffB
=B\)B�B  Bz�B��BG�B��B{B�\B�RB33B�B(�Bz�B��BG�BB�BffB�HB33B�B  Bz�B��B�B��B{BffB�RB33B�B  Bz�B��BG�B��B{B�\B�HB33B�B (�B z�B ��B!G�B!B"{B"ffB"�HB#\)B#�B$  B$z�B$��B%G�B%��B&{B&ffB&�RB'33B'�B(  B(Q�B(��B)G�B)��B)�B*ffB*�HB+33B+\)B+�
B,Q�B,��B,��B-p�B-�B.{B.ffB.�HB/\)B/�B/�
B0Q�B0��B0��B1p�B1B2{B2ffB2�HB333B3\)B3�
B4Q�B4��B4��B5G�B5B5�B6=qB6�RB7
=B7�B7�B8(�B8z�B8��B8��B9p�B9�B:=qB:ffB:�RB;33B;\)B;�
B<(�B<z�B<��B=G�B=��B=B>{B>�\B>�HB?33B?�B@  B@(�B@z�B@��BAG�BA��BA�BBffBB�\BC
=BC\)BC�BD  BDQ�BD��BD��BEp�BEBF{BFffBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    ?�z�@�\@@  @}p�@��R@��R@޸R@�p�A�RA#33A@��A`  A�  A�Q�A��A�  A�Q�A�Q�A��A�Q�B (�B  B�
B  B   B(  B0(�B8  B@  BG�BP(�BX  B`  Bh  Bo�
Bx  B�{B�(�B�(�B�  B��
B��B��B��B�{B�{B��B��B�  B��B��
B�B��
B��B�{B�  B��
B��B�{B�  B��B��B��B�(�B�{B�{B�  B�  C   C��C��C��C��C	��C{C
=C��C  C{C
=C  C
=C{C
=C��C"  C$
=C&{C({C*{C,
=C.
=C0
=C2  C3�C5�HC7��C:
=C<  C>  C?��CB  CD{CF  CH  CJ  CL
=CN
=CO��CR  CT
=CU��CW�HCY��C\
=C]��C_�Cb
=Cd{Cf  Cg�Cj  Cl
=Cn  Cp  Cq��Cs��Cu��Cx  Cz{C|  C~  C�\C�\C�  C���C���C���C���C��C���C�  C�  C���C�  C���C�  C�\C�C�  C���C�  C�
=C���C��C�  C�
=C�C�  C���C�  C�  C���C���C�C�  C���C���C�C�  C�  C�
=C�
=C�C�C�
=C�
=C�
=C���C��C���C���C�  C�  C�  C�C�
=C�  C���C�  C�C���C��C��C���C�C�\C�C���C���C�  C�C�
=C�  C��C���C�  C�C�C�C�C�  C���C���C���C���C���C��C��C��C��C��C��C�  C�  C���C�
=C�
=C�C�  C���C�  C�C�  C���C���C�  C�
=C���C���C�  C���C�  C�C���C��C���C�
=C�  C���C���C���C�C�C�C�  C���C���C���C���D �D � D  D}qD  D}qD��DxRD�RDxRD  D��D�D�DD��D  D}qD	  D	��D
D
� D
�qD� DD��D  D� D��D� DD� D�D��D  D� D�D� D  D� D  D��D�qD}qD�qDz�D�qD}qD  D��D  D}qD  D� D�qD� D�D��D�D� D�D� D�qD}qD   D ��D!  D!}qD!��D"xRD"�qD#�D$�D$� D$�qD%z�D%��D&}qD'�D'}qD'�RD(}qD)�D)��D*�D*z�D*�qD+��D+�qD,z�D,��D-� D.�D.�D/D/� D0  D0� D1  D1��D2�D2� D3  D3� D4  D4� D4�qD5}qD5��D6}qD7  D7z�D7�qD8}qD8��D9}qD:�D:� D;�D;��D<�D<�D<�qD=z�D>�D>�D?�D?}qD@�D@��D@�qDA}qDB  DBxRDC  DC��DD  DDz�DD�qDE��DF  DF� DF��DGz�DG�RDH}qDI  DI� DI�qDJz�DJ�qDK�DL�DL}qDM�DM�DN  DN� DODO� DP�DP�DQ�DQ� DQ��DR}qDS�DS}qDT  DT�DU�DU��DV  DV� DW  DW}qDW�RDX� DYDY� DZ  DZ� D[  D[�D[�qD\xRD]  D]��D^�D^� D^�qD_z�D_��D`z�D`�qDaz�Da�qDb� Dc�Dc�Dd�Dd}qDe  De� Df  Df��DgDg}qDg�RDh}qDh�qDi� Dj�Dj�Dk  Dkz�Dl  Dl� Dm  Dm}qDn  Dn��Do  Do}qDp  Dp��Dq  Dq� DrDr��Ds�Ds� Dt  Dt}qDt��Du��Dv  Dvz�Dv�RDwz�Dx�Dx}qDy�Dy� Dy�qDz� Dz�qD{� D|  D|� D}�D}��D~�D~z�D~�qD}qD�qD�>�D�~�D���D��qD�@ D�� D�� D��D�@ D�}qD��qD��qD�>�D�� D��HD��D�C�D��HD���D��qD�@ D��HD�D�  D�>�D�~�D��qD�HD�>�D�|)D���D�  D�B�D�� D���D���D�>�D�}qD��qD��)D�<)D�� D�D�HD�@ D�� D�� D�  D�>�D�� D�� D���D�@ D�� D��HD�  D�>�D�}qD��qD�HD�B�D���D�D���D�>�D��HD�� D�  D�@ D�~�D��HD�  D�=qD�~�D�� D���D�@ D��HD���D���D�AHD�� D��qD�  D�@ D�~�D��HD���D�>�D�}qD��qD���D�AHD��HD��HD��D�B�D���D�� D��qD�=qD�|)D��qD��qD�>�D�~�D���D��qD�>�D�~�D�� D�  D�>�D�� D��HD�HD�B�D��HD���D�HD�>�D�~�D���D��qD�AHD��HD���D�HD�@ D�� D�� D�  D�>�D�}qD���D���D�@ D�}qD���D�HD�@ D�� D��HD��qD�>�D�� D���D���D�@ D�� D�� D���D�>�D�~�D���D���D�@ D�� D�� D�  D�=qD�~�D�� D�  D�B�D���D�� D��qD�>�D��HD�� D��qD�>�D�� D�� D�  D�@ D�� D�� D�HD�AHD��HD���D���D�AHD�� D�� D�  D�@ D�~�D��qD���D�=qD�� D�� D���D�>�D�~�D���D���D�@ D�� D��HD�HD�>�D�� D�� ?L��?u?�=q?��R?�{?�p�?���?�(�?��@�\@��@z�@(�@#�
@.{@5@@  @J=q@Q�@^�R@c�
@k�@s33@}p�@��
@���@�{@��@�@��H@��R@��
@��@��@��@�@���@��R@��
@Ǯ@���@У�@�@��H@޸R@�\@�@�@��@�@���@�p�A ��A33AAQ�A
�HA(�A�RA��A33AAQ�A=qA(�A{A ��A#33A%A(Q�A*=qA,(�A.�RA0��A333A5A8Q�A:�HA=p�A?\)AA�AC�
AEAH��AJ�HAN{AP  AQ�ATz�AVffAX��AZ�HA^{A`��Ab�\Ae�Ag
=Ah��Ak�An{Ap��As33Au�Aw
=Ax��A|(�A~{A�Q�A�G�A�=qA��A���A�{A�
=A�  A�G�A��\A��A�z�A�A�\)A�Q�A�G�A��\A��A��A�{A�
=A�  A�G�A��\A��A���A�{A�
=A�Q�A���A��\A��A���A�{A�
=A�  A�G�A��\A��
A���A�{A�
=A�  A�G�A��\A��A�z�A�A�
=A�  A���A�=qA�33A�z�A�p�A��RA�  A���A��A�33A�z�A�A�ffAǮA���A��A��HA�(�A�p�A�ffA�\)A�Q�A��Aҏ\AӅA���A�A�
=A�  A�G�A�=qA�33A�(�A�p�A�ffA�\)A��AᙚA�\A�A���A�A�RA�A���A��A��HA��
A��A�{A�
=A�  A�G�A�=qA��HA��
A��A�{A��RA�  A���A��A��\A��A�z�A�p�A�ffA�
=B   B z�B ��BG�BB=qB�\B�HB\)B�
B(�Bz�B��Bp�BB=qBffB�HB\)B�B  Bz�B��B	�B	p�B	�B
ffB
�RB
=B�B�
BQ�Bz�B��Bp�BB{BffB
=B\)B�B  Bz�B��BG�B��B{B�\B�RB33B�B(�Bz�B��BG�BB�BffB�HB33B�B  Bz�B��B�B��B{BffB�RB33B�B  Bz�B��BG�B��B{B�\B�HB33B�B (�B z�B ��B!G�B!B"{B"ffB"�HB#\)B#�B$  B$z�B$��B%G�B%��B&{B&ffB&�RB'33B'�B(  B(Q�B(��B)G�B)��B)�B*ffB*�HB+33B+\)B+�
B,Q�B,��B,��B-p�B-�B.{B.ffB.�HB/\)B/�B/�
B0Q�B0��B0��B1p�B1B2{B2ffB2�HB333B3\)B3�
B4Q�B4��B4��B5G�B5B5�B6=qB6�RB7
=B7�B7�B8(�B8z�B8��B8��B9p�B9�B:=qB:ffB:�RB;33B;\)B;�
B<(�B<z�B<��B=G�B=��B=B>{B>�\B>�HB?33B?�B@  B@(�B@z�B@��BAG�BA��BA�BBffBB�\BC
=BC\)BC�BD  BDQ�BD��BD��BEp�BEBF{BFffBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A%XA%A$�A#dZA"��A"Q�A"-A"JA!�mA!x�A!VA �yA ȴA �!A ��A ��A ��A �DA �A�FAG�A;dAA�FA1'A��A�hA33A�A
�HA�A�\A��@�r�@���@��@��w@��9@�G�@�@�?}@�P@�E�@�ff@���@홚@�7L@�Q�@��@�^@�w@�v�@��@�J@���@��@�V@�@�&�@�9X@�1@�;d@���@��@�\@�P@�|�@���@��@�^5@��@��@��@��m@��;@���@�@�dZ@�o@���@�ȴ@�E�@���@ᙚ@��#@�V@�~�@�\@��@⟾@�~�@�-@��@��#@�^@�X@�O�@���@��@�Z@�A�@�  @��m@���@ߕ�@��@�Q�@��u@��D@�(�@���@�A�@��@��@��@�Z@� �@� �@�I�@�Z@�bN@�j@��@��@��@��u@��@�@�9@��D@�r�@�bN@�Q�@�I�@�Z@�Z@�Z@�A�@�b@��
@߶F@ߥ�@ߕ�@�t�@�33@�o@�o@�33@�+@�33@�33@�"�@�o@�
=@��H@ާ�@ާ�@�v�@�E�@�5?@�-@�$�@��@�@݉7@݉7@ݑh@ݑh@�X@�X@�O�@�O�@�?}@��@��`@ܼj@ܴ9@ܣ�@�j@�Z@�I�@�1'@��@�b@ۮ@�S�@��@��@ڗ�@�$�@�J@�@ٲ-@�x�@�p�@�G�@�/@�7L@ش9@؛�@؛�@أ�@ج@���@��`@؛�@�r�@�j@�bN@�Z@�A�@� �@׾w@ו�@�K�@��y@�n�@�J@��T@�x�@��@ԣ�@�A�@�1@ӶF@�t�@�
=@�v�@�ff@�~�@җ�@җ�@ѡ�@�V@�V@��`@��@��`@Гu@�  @϶F@υ@��@��@�ȴ@Χ�@�~�@��@��#@�p�@��`@�b@��@��m@˶F@�dZ@�"�@��H@�ȴ@ʸR@ʟ�@�@Ɂ@ɉ7@ɑh@�x�@�`B@�?}@�Ĝ@ȃ@��@� �@ǶF@�K�@�
=@Ɨ�@�=q@���@���@�O�@��/@���@Ĭ@ċD@���@���@��m@��@���@�  @�  @���@��@�t�@��@�n�@�E�@���@��7@�X@��@�Ĝ@��;@��
@���@���@���@��P@���@�l�@�"�@�@�ȴ@���@���@��`@���@��@�E�@���@�X@�7L@�V@�bN@�  @��m@�dZ@�V@��h@�G�@�&�@��@��@�&�@�/@�7L@�/@�&�@���@��D@�j@�I�@�  @���@��F@��@�o@��\@�V@���@��^@�x�@�?}@��@�Q�@�1@�;d@��y@��@��H@���@�^5@�E�@��#@���@���@���@�hs@�O�@�7L@�&�@��@�z�@�1'@���@��y@�$�@�x�@���@��@���@��D@�j@�9X@���@�@��+@�E�@�$�@���@�x�@��@���@�Q�@��@��
@���@��@�t�@�
=@��!@���@���@�M�@��^@�X@�%@�(�@��w@�\)@�;d@��@�~�@��@��T@��@��#@���@��^@�x�@��@���@���@�r�@�j@�j@�Q�@�I�@�9X@�b@��m@��F@�+@���@���@�ff@���@�/@�V@��@�j@���@�;d@���@�E�@���@�@��7@�`B@�X@��@��9@�j@��@���@��w@��P@�"�@��!@��+@�v�@�ff@�V@�E�@�-@�-@��@�{@�{@��@��@��@��@��T@���@���@�x�@�`B@�`B@�`B@�G�@�?}@�&�@��`@��j@���@��u@�bN@�9X@�b@��w@��@�\)@�\)@�;d@�33@�+@��@�o@��@�@���@��R@��!@���@���@���@���@���@�v�@�E�@�M�@�E�@�5?@�5?@�5?@�$�@�{@��@��T@���@��7@�p�@�/@�&�@���@��`@���@�Ĝ@���@��u@��u@��@��@�
=@��y@���@�=q@��#@��-@���@�x�@�hs@�X@�X@�G�@��@��@��@���@���@��/@���@��j@��9@��@�1'@���@��y@���@��\@�~�@�E�@��@���@��@��#@���@��^@��7@��7@�x�@�p�@�X@�G�@�&�@�V@�V@���@��@���@��`@�Ĝ@��9@��@���@��u@��D@�z�@�z�@�z�@�j@�I�@�1'@�(�@��@�b@���@��@���@��@��;@��m@�ƨ@��F@��@���@���@��@�C�@�33@�33@�+@�"�@�o@�
=@�@��H@��!@�V@�J@��@���@��-@���@��h@��7@�x�@�G�@�/@�&�@�&�@��@��@�V@�V@�V@��/@���@���@��9@���@��D@�j@�9X@�  @�P@|�@|�@l�@l�@l�@\)@K�@K�@;d@;d@;d@+@
=@
=@~�@~�@~�@~ȴ@~�R@~��@~��@~��@~��@~��@~��@~��@~��@~��@~��@~�+@~��@~��@~��@~��@~��@~��@~��@~�R@~�R@~��@~�R@~�R@~ȴ@~ȴ@~ȴ@~ȴ@~�+@~ff@~v�@~�+@~�+@~�+@~�+@~�+@~�+@~��@~��@~�+@~�+@~�+@~�+@~v�@~ff@~E�@~$�@}�@}��@}`B@}V@}V@}/@}/@}�@}�@}V@}V@|��A%33A%t�A%�A%p�A%?}A%?}A%?}A%K�A%hsA%"�A$�/A$�uA$~�A$~�A$�DA$�A$JA#��A#�
A#A#��A#�PA#�A#`BA#O�A#7LA#&�A#"�A#oA"��A"�HA"��A"�\A"jA"^5A"ZA"ZA"Q�A"Q�A"M�A"I�A"=qA"9XA"1'A"(�A"(�A"$�A"�A"{A"bA"JA"JA"JA"1A"1A"JA"  A!�A!�mA!�TA!��A!�^A!�-A!��A!��A!��A!��A!�hA!�A!l�A!7LA!"�A!�A!�A!oA!oA!oA!VA!VA!VA!oA!VA!VA!VA!VA!
=A!
=A!%A ��A ��A �A �A �A �yA �`A �HA �/A �A ��A ��A ��A ��A ��A ��A ��A ��A ȴA ȴA ȴA ȴA ĜA ĜA ĜA ��A ��A �jA �jA �RA �!A �!A �A �!A �A �A �A �A �A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A �\A �\A �uA ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A �uA �DA �+A �A z�A v�A n�A bNA ^5A Q�A Q�A E�A 9XA 1'A (�A �A {A {A {A   A��A�TA�-A�^AAA�wA�wA�wA�FA�FA�FA�FA�FA�-A�-A�A�A��A��A��A��A��A��A��A�hA�PA�A`BAG�A&�A
=A�A�A��A�RA��A�+An�AE�A  A��A�PA�7At�AXA��A�/A��A��A��A��A��A��A�jA�RA�9A�9A�A�A�A��A��A�AVA1'AG�A�AZA�A�;A�^A�PAhsA;dAA�RAr�AE�AA�A5?A�AA��A��A�A�#A�^AA��A�^A�hA�7A|�A?}AC�A;dA�A��A��A�Ar�A�uA��A��An�AVA9XA�A  A��A��A�
A�mA�;A�A�A�A�#A�wA�^A��A�A��A��A��A��A��A��A��A��A��A��A��A��A��A�PA��A��A�PA�PA�7A�PA�PA�hA�hA�PA�hA�PA��A��A��A��A��A�PA�hA��A�A�A�PA�A�7A�7Al�A`BAO�AS�AXA?}A7LA33A7LA;dA"�A"�A&�A�A�A%A�`A��A��A��A��A�`A�/AĜA�DA9XA$�AAƨAdZA�RA��A�+A�+A^5A�A��A�Al�A/A�jA�uAQ�AJA�^A��Ax�AO�AA�A�RAr�AJA;dA
��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    A%XA%A$�A#dZA"��A"Q�A"-A"JA!�mA!x�A!VA �yA ȴA �!A ��A ��A ��A �DA �A�FAG�A;dAA�FA1'A��A�hA33A�A
�HA�A�\A��@�r�@���@��@��w@��9@�G�@�@�?}@�P@�E�@�ff@���@홚@�7L@�Q�@��@�^@�w@�v�@��@�J@���@��@�V@�@�&�@�9X@�1@�;d@���@��@�\@�P@�|�@���@��@�^5@��@��@��@��m@��;@���@�@�dZ@�o@���@�ȴ@�E�@���@ᙚ@��#@�V@�~�@�\@��@⟾@�~�@�-@��@��#@�^@�X@�O�@���@��@�Z@�A�@�  @��m@���@ߕ�@��@�Q�@��u@��D@�(�@���@�A�@��@��@��@�Z@� �@� �@�I�@�Z@�bN@�j@��@��@��@��u@��@�@�9@��D@�r�@�bN@�Q�@�I�@�Z@�Z@�Z@�A�@�b@��
@߶F@ߥ�@ߕ�@�t�@�33@�o@�o@�33@�+@�33@�33@�"�@�o@�
=@��H@ާ�@ާ�@�v�@�E�@�5?@�-@�$�@��@�@݉7@݉7@ݑh@ݑh@�X@�X@�O�@�O�@�?}@��@��`@ܼj@ܴ9@ܣ�@�j@�Z@�I�@�1'@��@�b@ۮ@�S�@��@��@ڗ�@�$�@�J@�@ٲ-@�x�@�p�@�G�@�/@�7L@ش9@؛�@؛�@أ�@ج@���@��`@؛�@�r�@�j@�bN@�Z@�A�@� �@׾w@ו�@�K�@��y@�n�@�J@��T@�x�@��@ԣ�@�A�@�1@ӶF@�t�@�
=@�v�@�ff@�~�@җ�@җ�@ѡ�@�V@�V@��`@��@��`@Гu@�  @϶F@υ@��@��@�ȴ@Χ�@�~�@��@��#@�p�@��`@�b@��@��m@˶F@�dZ@�"�@��H@�ȴ@ʸR@ʟ�@�@Ɂ@ɉ7@ɑh@�x�@�`B@�?}@�Ĝ@ȃ@��@� �@ǶF@�K�@�
=@Ɨ�@�=q@���@���@�O�@��/@���@Ĭ@ċD@���@���@��m@��@���@�  @�  @���@��@�t�@��@�n�@�E�@���@��7@�X@��@�Ĝ@��;@��
@���@���@���@��P@���@�l�@�"�@�@�ȴ@���@���@��`@���@��@�E�@���@�X@�7L@�V@�bN@�  @��m@�dZ@�V@��h@�G�@�&�@��@��@�&�@�/@�7L@�/@�&�@���@��D@�j@�I�@�  @���@��F@��@�o@��\@�V@���@��^@�x�@�?}@��@�Q�@�1@�;d@��y@��@��H@���@�^5@�E�@��#@���@���@���@�hs@�O�@�7L@�&�@��@�z�@�1'@���@��y@�$�@�x�@���@��@���@��D@�j@�9X@���@�@��+@�E�@�$�@���@�x�@��@���@�Q�@��@��
@���@��@�t�@�
=@��!@���@���@�M�@��^@�X@�%@�(�@��w@�\)@�;d@��@�~�@��@��T@��@��#@���@��^@�x�@��@���@���@�r�@�j@�j@�Q�@�I�@�9X@�b@��m@��F@�+@���@���@�ff@���@�/@�V@��@�j@���@�;d@���@�E�@���@�@��7@�`B@�X@��@��9@�j@��@���@��w@��P@�"�@��!@��+@�v�@�ff@�V@�E�@�-@�-@��@�{@�{@��@��@��@��@��T@���@���@�x�@�`B@�`B@�`B@�G�@�?}@�&�@��`@��j@���@��u@�bN@�9X@�b@��w@��@�\)@�\)@�;d@�33@�+@��@�o@��@�@���@��R@��!@���@���@���@���@���@�v�@�E�@�M�@�E�@�5?@�5?@�5?@�$�@�{@��@��T@���@��7@�p�@�/@�&�@���@��`@���@�Ĝ@���@��u@��u@��@��@�
=@��y@���@�=q@��#@��-@���@�x�@�hs@�X@�X@�G�@��@��@��@���@���@��/@���@��j@��9@��@�1'@���@��y@���@��\@�~�@�E�@��@���@��@��#@���@��^@��7@��7@�x�@�p�@�X@�G�@�&�@�V@�V@���@��@���@��`@�Ĝ@��9@��@���@��u@��D@�z�@�z�@�z�@�j@�I�@�1'@�(�@��@�b@���@��@���@��@��;@��m@�ƨ@��F@��@���@���@��@�C�@�33@�33@�+@�"�@�o@�
=@�@��H@��!@�V@�J@��@���@��-@���@��h@��7@�x�@�G�@�/@�&�@�&�@��@��@�V@�V@�V@��/@���@���@��9@���@��D@�j@�9X@�  @�P@|�@|�@l�@l�@l�@\)@K�@K�@;d@;d@;d@+@
=@
=@~�@~�@~�@~ȴ@~�R@~��@~��@~��@~��@~��@~��@~��@~��@~��@~��@~�+@~��@~��@~��@~��@~��@~��@~��@~�R@~�R@~��@~�R@~�R@~ȴ@~ȴ@~ȴ@~ȴ@~�+@~ff@~v�@~�+@~�+@~�+@~�+@~�+@~�+@~��@~��@~�+@~�+@~�+@~�+@~v�@~ff@~E�@~$�@}�@}��@}`B@}V@}V@}/@}/@}�@}�@}V@}V@|��A%33A%t�A%�A%p�A%?}A%?}A%?}A%K�A%hsA%"�A$�/A$�uA$~�A$~�A$�DA$�A$JA#��A#�
A#A#��A#�PA#�A#`BA#O�A#7LA#&�A#"�A#oA"��A"�HA"��A"�\A"jA"^5A"ZA"ZA"Q�A"Q�A"M�A"I�A"=qA"9XA"1'A"(�A"(�A"$�A"�A"{A"bA"JA"JA"JA"1A"1A"JA"  A!�A!�mA!�TA!��A!�^A!�-A!��A!��A!��A!��A!�hA!�A!l�A!7LA!"�A!�A!�A!oA!oA!oA!VA!VA!VA!oA!VA!VA!VA!VA!
=A!
=A!%A ��A ��A �A �A �A �yA �`A �HA �/A �A ��A ��A ��A ��A ��A ��A ��A ��A ȴA ȴA ȴA ȴA ĜA ĜA ĜA ��A ��A �jA �jA �RA �!A �!A �A �!A �A �A �A �A �A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A �\A �\A �uA ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A ��A �uA �DA �+A �A z�A v�A n�A bNA ^5A Q�A Q�A E�A 9XA 1'A (�A �A {A {A {A   A��A�TA�-A�^AAA�wA�wA�wA�FA�FA�FA�FA�FA�-A�-A�A�A��A��A��A��A��A��A��A�hA�PA�A`BAG�A&�A
=A�A�A��A�RA��A�+An�AE�A  A��A�PA�7At�AXA��A�/A��A��A��A��A��A��A�jA�RA�9A�9A�A�A�A��A��A�AVA1'AG�A�AZA�A�;A�^A�PAhsA;dAA�RAr�AE�AA�A5?A�AA��A��A�A�#A�^AA��A�^A�hA�7A|�A?}AC�A;dA�A��A��A�Ar�A�uA��A��An�AVA9XA�A  A��A��A�
A�mA�;A�A�A�A�#A�wA�^A��A�A��A��A��A��A��A��A��A��A��A��A��A��A��A�PA��A��A�PA�PA�7A�PA�PA�hA�hA�PA�hA�PA��A��A��A��A��A�PA�hA��A�A�A�PA�A�7A�7Al�A`BAO�AS�AXA?}A7LA33A7LA;dA"�A"�A&�A�A�A%A�`A��A��A��A��A�`A�/AĜA�DA9XA$�AAƨAdZA�RA��A�+A�+A^5A�A��A�Al�A/A�jA�uAQ�AJA�^A��Ax�AO�AA�A�RAr�AJA;dA
��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�BjBjBl�Bm�Bo�Bo�Bo�Bo�Bp�Bs�Bu�Bu�Bv�Bx�Bz�B|�B~�B� B�B�=B�1B�bB��B�B��B��B��B�uB�Bm�BdZB_;B`BBhsB� B��B�XB�wBB��B�RB�-B�LB�}BÖBƨBǮBŢBŢB�}B�RB�RB�^B�qB�qB��BƨB�
B��BB+BB��B��B�B�mB�mB�yB�mB�`B�ZB�ZB�TB�ZB�`B�`B�fB�`B�ZB�TB�TB�NB�TB�TB�fB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B  BBBB1B1B+B+B+B+B	7BDBJBPBVB\BhBhBoBoBuB{B{B{B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B{BuBoBhBbBhBuB{B�BoBhBhBhBuBoBuBuBuBoBhBhBhBbB\B\B\BVBPBVB\B\B\BbBbBbBbBbBbBuB{B{B{B{BuBuBuB{B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B �B �B �B�B �B �B �B �B �B �B �B �B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B{BuBoBhBhBhBbBbB\B\BPBJBJBJBDB
=B
=B
=B
=B
=B
=BDBDBDBDBDBDBDB
=B1B+B%BBBBBBBBB  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�sB�sB�mB�mB�fB�`B�`B�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�NB�NB�NB�NB�NB�NB�HB�HB�HB�HB�HB�HB�HB�HB�NB�NB�NB�NB�NB�NB�NB�NB�TB�TB�TB�TB�ZB�`B�`B�`B�`B�`B�`B�`B�ZB�ZB�ZB�ZB�TB�TB�TB�TB�TB�TB�TB�NB�TB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�HB�NB�NB�HB�HB�HB�HB�HB�HB�HB�HB�HB�BB�BB�BB�;B�;B�;B�;B�;B�;B�;B�5B�/B�)B�)B�#B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�
B�
B�
B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B�B��B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BjBgmBgmBo�BiyBjBiyBffBhsBl�Bo�Bk�Bk�BffBm�Bl�Bk�Bm�Bl�Bn�Bm�Bm�Bn�Bm�Bo�Bn�Bm�Bm�Bn�Bn�Bn�Bp�Bp�Bp�Bo�Bo�Bo�Bo�Bn�Bo�Bo�Bp�Bo�Bp�Bo�Bo�Bo�Bp�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bp�Bp�Bo�Bp�Bq�Bq�Bq�Br�Bq�Bq�Bq�Br�Br�Bu�Bv�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bt�Bt�Bt�Bu�Bt�Bu�Bt�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bv�Bu�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�Bx�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�B{�Bz�B{�B{�Bz�Bz�B{�Bz�B{�B{�B{�B{�B{�B{�B|�B|�B|�B}�B}�B~�B~�B~�B~�B~�B~�B� B~�B~�B� B� B~�B� B� B� B~�B~�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B� B� B�B�B�B�B�B�B�B�%B�%B�%B�7B�%B�1B�=B�=B�=B�=B�=B�DB�=B�7B�=B�=B�7B�7B�=B�7B�7B�1B�1B�+B�%B�+B�+B�+B�%B�+B�+B�=B�+B�1B�1B�1B�1B�1B�1B�+B�1B�7B�\B�hB�\B�bB�hB�uB��B�uB��B�{B�uB�oB�oB�hB�hB�hB�bB�\B�\B�JB�JB�DB�1B�7B�JB�7B��B�B��B��B��B��B��B��B��B�B�B�!B�B�B�B�B�B�B�B�B�B�!B�B��B�B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�oB��B�uB�uB�oB�uB�hB�{B�hB�\B�VB�PB�VB�PB�=B�DB�bB�1B�%B�B�+B�JB�B� B~�B~�B~�B�Bz�Bx�Bv�B}�B{�By�Bw�Bv�Bq�Br�Bp�Bo�BjBk�Bn�BjBn�Br�BffG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    Bm�Bo�Bp�Bp�Br	Bp�BpyBp�Br�Bt�Bv0Bv$BwBx�Bz�B|�B4B�<B�XB��B��B��B�B�]B�{B��B�B��B�hBxbBl8Bg�Bj�Bl�B��B��B�BÉB�0B�<B��B�B�5B�B�,B�OB�B��BȤBB�XB�B�HB��B�KB��B�hB�sB��B?BpBzB��B�B��B�B�B�~B�B�B�2B��B�B�hB�yB�B��B��B��B�TB�B�B�B��B�B�4B�lB�xB�B��B�'B�B��B��B�=B��B�lB�&B��B��B�B��B��B�B�PB�HB��B�B �BFB �B�B"B`BnBB&B�B	B4B>B*B%B]BzBSB_BgB�B�B�B�B�BhB�B�B�B�B�B�B�B�B�B�B�B�BbB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B0B"B�B�B�B.B�B�B�B�B�B�B�B�BHB�B�B�B�BYB�BB�B�B�B�B�B�B2B�BB0BJB(B�B,BB.BB�B�B�BB@BzBEBOB�B�BWBtB�BaB�B�BPB�B�B
B�B�B�B�B�B�BB6B�B�BjB�B�B�B�B�BzB�B�B�BoBrB�B�B�B,B�BBvB%B0B�BEB&BB�BXBTB�B�B�BB�B�B�B�B�B�B�B�B|B�BXB B oB!8B!B!hB!B!B �B!B �B �B �B �B!B!6B �B B B 'B�BBaB�B�BB�B�B�B2B�BlB:B�B�B�B�B�B�B�B�B�B�B�B<B�B�B�BBgB�B0B;B�B�B�B�B�B�BQB�B�B�BcBCB�B�B
lB
�B
]B
nB
ZB
�BlBmB`BmB*B�B!B`B	YB<B�B�B5B4BOBfB�B@B�B iB�2B��B�B��B�eB��B�=B�BB�8B�B��B�wB�^B��B��B�DB��B�\B�CB�B�YB�BB��B��B�B��B�B�B�B�B�B�B�+B�QB�B��B�B�B��B��B��B��B��B��B�aB�B��B��B�tB�TB�B�B�8B�B�B�RB��B��B�B�B�B�gB�B��B��B��B�B�B�B��B��B�B�dB�bB�bB�dB�nB�OB�hB�]B�OB�DB�MB�RB�B�mB�B�sB�B�~B�`B�aB�B�oB�B��B�B�B�kB�B�B�B��B�B�B�XB�B�bB�`B�gB�^B�FB�vB�B�vB�[B�]B�fB�QB�PB�PB�|B�B�BB�WB�aB�HB�GB�`B�cB�{B�bB�eB�B�oB�B�QBߌB�YB�]B�KB�kB�PB�DB�*B�iB�HB�]BےB��BڶB�]B�8B�VB�3B�2B�B�8B�QB�B�B�5B�B�5B�$B�$B�B�XB֐B��B�B�vB�B�B�WB�>B�3B�B�B�B�(B�FB�B�B�B�)B�B�4B�)B� B�B�B��B�B�0B�B�B�B�B�B�B�B�B�B�5B�&B�B�B�B� B�
B��B�	B�B��B�/B�B�
B�B�B�&B�^B�B��B�B�B�B�B�B�0B�KBԂB�lB�-B�*B�B�B�B�B�B�:B�B�B��B��B��B��B��B��B�?B� B�B�B�B�B�(B�8B�HB�EB��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B��B�B�(B�B�#B��B��B��B��B��B��B��B��B��BjBgmBgmBo�BiyBjBiyBffBhsBl�Bo�Bk�Bk�BffBm�Bl�Bk�Bm�Bl�Bn�Bm�Bm�Bn�Bm�Bo�Bn�Bm�Bm�Bn�Bn�Bn�Bp�Bp�Bp�Bo�Bo�Bo�Bo�Bn�Bo�Bo�Bp�Bo�Bp�Bo�Bo�Bo�Bp�Bo�Bo�Bo�Bo�Bo�Bo�Bo�Bn�Bp�Bp�Bo�Bp�Bq�Bq�Bq�Br�Bq�Bq�Bq�Br�Br�Bu�Bv�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bt�Bt�Bt�Bu�Bt�Bu�Bt�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bu�Bv�Bu�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bv�Bw�Bw�Bw�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�Bx�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�B{�Bz�B{�B{�Bz�Bz�B{�Bz�B{�B{�B{�B{�B{�B{�B|�B|�B|�B}�B}�B~�B~�B~�B~�B~�B~�B� B~�B~�B� B� B~�B� B� B� B~�B~�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B� B� B�B�B�B�B�B�B�B�%B�%B�%B�7B�%B�1B�=B�=B�=B�=B�=B�DB�=B�7B�=B�=B�7B�7B�=B�7B�7B�1B�1B�+B�%B�+B�+B�+B�%B�+B�+B�=B�+B�1B�1B�1B�1B�1B�1B�+B�1B�7B�\B�hB�\B�bB�hB�uB��B�uB��B�{B�uB�oB�oB�hB�hB�hB�bB�\B�\B�JB�JB�DB�1B�7B�JB�7B��B�B��B��B��B��B��B��B��B�B�B�!B�B�B�B�B�B�B�B�B�B�!B�B��B�B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�oB��B�uB�uB�oB�uB�hB�{B�hB�\B�VB�PB�VB�PB�=B�DB�bB�1B�%B�B�+B�JB�B� B~�B~�B~�B�Bz�Bx�Bv�B}�B{�By�Bw�Bv�Bq�Br�Bp�Bo�BjBk�Bn�BjBn�Br�BffG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                    <+L�<7L<0~�<+��<(X~<$z�<$i&<$�;<&�}<%F<#�H<#�<#�l<#�8<#�c<#�<#�N<%�<% <%��<@�	<i�T<u�<2m<&y<#�5<%�<R�|<�e�<nP�<N�c<T@�<ka�<1Q�<(�,<#�Q<%04<6�<)'7<.T�<(�a<&�<#ا<$G<$�<$,<%MY<$ح<*��<*ٜ<'<$Sa<#؄<#�<#�r<$MO<+��<Bշ<4�g<#�*<%'<,�u<%�<'��<1�B<#�D<#�m<$�j<$�.<$�(<$e.<%��<#�<#ף<#��<#�N<#�N<$v<$/<#�
<$@|<$><<#��<#�W<$5w<#�<#��<#�*<#�<#�N<$�<#�(<#�*<#�E<$�<#�D<$B�<$�<#�U<#�+<#��<#�+<#�r<#�M<$�<$�<#�W<#׎<$�<#�<#��<$�<#׺<#��<#�<#�<#�<#�<#�o<#��<#�{<#�r<#�^<#�<#�<#�c<#��<#�<<#�<#��<#�<#��<#׺<#��<#�<#�<#�+<#�M<#�5<#�<#�c<#ٛ<#��<#��<#ޫ<#�<#ܯ<#�0<#׎<#�<#؄<#��<#�<#�<#�<#�<<#�<#�<#�o<#׎<#�<#�<#�<#�<#�
<#�X<#�<#�<#�
<#ף<#�<#�]<#�<#�4<#�<#ף<#�<#�5<#��<#��<#ۮ<#��<#�C<$F<$r<#�<#��<#��<$&<#�<#��<#�]<#�5<#��<#�E<#��<#�0<$A�<#�l<#�<#�I<#�<#�N<#�C<#��<#��<#��<#׎<#׎<#�r<#��<$}<#��<$�<$"2<$=<$ K<#��<$*<$!><$2G<$�<#�N<$v<#��<$#(<$e.<#��<#��<#�r<#�C<%MY<$z�<#�{<#�N<#�0<#��<$f<$i&<#�H<#�<$ K<#�<#��<#��<#�<$f<#�N<$-<$o�<$�<#ޫ<#ף<#�!<$a<#��<#�W<#�r<#��<#ޫ<%2?<#��<#�{<#�I<#ۮ<#�r<#�<$=<#��<$�<#�<$)
<$.<#�H<$7�<$�<#��<#�!<$><<$3U<#�D<#ߜ<#�E<$^�<#�<#ڑ<#�i<#�{<#�X<#�<#��<#��<$L<$� <$+<#�<$6�<#�a<#�<$(<#�!<%G<#�<#�<#�<#׺<#�i<#�&<#�&<#��<#��<#��<#�<%p<$��<%v�<&�<$��<$��<#��<#ߜ<#�<$�<$}<#ޫ<$f�<&�<%�<#�a<#�J<#�X<#�<#�{<#�{<#�{<#׺<#�<#�&<$4e<#�N<#��<#��<$�<#�<#�<$:�<$Gd<#�W<$�<#��<#�W<#�<$r<$� <$r<$��<$�<#��<#�0<#��<#��<#��<$#(<#�*<#�^<#ٛ<#�<#��<#�+<#�o<#�+<$x+<$�<$k�<$ح<$�<$��<$R'<#��<#�o<#�D<#��<#�<$&<$��<$L<#��<#��<$�<$�<$}<$<<$N�<#��<#�<#�"<#��<#�]<$'<$G<#�<#�]<$ �<$|d<$k<$�<%04<$.<$#(<#�E<#�<$�Q<$x+<#�i<#�X<#ٛ<#�o<#�o<#�H<$�<$?[<#��<#�&<#ף<#�<#��<#׺<#�o<#�<#�<#�<$j|<$�<#�N<#�5<$�<$q@<#��<#�e<$]h<$��<$5w<$�Q<$<$�<#�<#�<#�U<#�$<#��<$$<$<<$�<#�N<#��<#�<$-<$3U<#�<#�o<#�<#�<#�o<#�r<#�<#�<#׺<#�<#�X<#�<#�<#�<#��<#��<#��<#ޫ<#��<#�
<#�<#��<#׺<#�8<#��<#��<#��<#��<#�<#�<#��<$<<#��<#�&<#�<#ޫ<#ף<#�<#�$<#��<#�<<#��<#�<#��<#׎<#׺<#��<#�<#�<#�<#�J<#�<#�{<#׺<#��<#�
<#�<#��<#�D<#��<#�<#ٛ<#�N<#�8<#��<#؄<#�<#��<#ڑ<#��<#�<#�c<#׺<$��<%b<#��<#�N<#��<$&<$k<#�<#�]<#�&<#�o<#�D<#�
<#�]<#��<#ף<#�<#�8<#ף<#�8<#�<#�<#��<#�<$�<$m,<$�<$�<#��<#��<#�5<#�<#ߜ<#��<#�<#��<#�l<#��<#�&<#��<#׺<#ܯ<#�*<#��<#ܯ<#�<#��<#�<#�i<#؄<#ޫ<#�<#׎<#��<#�{<#�{<#��<#�
<#�<#�<#�^<#��<#׎<#�D<#ף<#ڑ<#�{<#�{<#�i<#ا<#�<<#�^<#�D<#�{<#׎<#׎<#��<#�(<#ٛ<#�<#׎<#׺<#�<#ף<#ף<#��<#�<$<$ <#ߜ<#ޫ<#ۮ<#ۮ<#ף<#�{<#�]<#��<#��<#��<#�
<#�{<#�X<#׎<#�
<#�<#�!<#ף<#׺<#��<#�<#�r<#��<#�<#��<#�5<#ף<#�
<#�{<#�<#�<#�{<#�{<#�
<#�{<#�<#�<#�{<#؄<#�<#ڑ<#�<#�<#׎<#׎<#�{<#�
<#�<#�<#�<#�X<#�<#�<#�
<#�<#�X<#ף<#�<#�<#�<#�{<#�X<#��<#׺<#�<#�<<#׎<#�<#��<#�<#�<#�<#�^<#�<#�I<#׎<#�<#�<#�<#�<#�<#�i<#�<#��<#�<#�<#�<#׎<#׺<#�D<#��<#�8<#�<#�J<#��<#�<#��<#�<#�i<#�
<#�X<#�
<#׎<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285720181106012857QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285720181106012857QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                