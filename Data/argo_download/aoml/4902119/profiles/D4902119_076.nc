CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  v   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92018-11-06T01:28:56Z creation; 2023-02-09T14:06:16Z DMQC;      
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
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  P\   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  UH   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  h�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  m�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �D   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �0   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �|   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �,   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` d   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                      HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                       HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   (   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   0   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � 8   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar           HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�          HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Argo profile    3.1 1.2 19500101000000  20181106012856  20230209090616  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               L   LAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�o&�2#1@�o&�2#111  @�o&�L0@�o&�L0@Nb�5��@Nb�5���<Ӓ:)�z�<Ӓ:)�z11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AB  ?�\)@�\@B�\@}p�@�p�@��R@�  A   A\)A#33A?\)A`  A�  A�  A�  A�  A�  A�Q�A߮A�\)B   B(�B  B�
B (�B(  B0  B8(�B?�
BG\)BPz�BX(�B_�Bh  Bp(�Bxz�B�B�B�B��
B�B��B��B�B��B��
B��B��B�  B��B��
B��B�{B�  B��
B��B�{B�{B�(�B�{B�(�B�(�B��
B��
B�{B�=qB��B�  C �C{C
=C{C{C

=C
=C
=C{C{C
=C�C(�C  C
=C{C �C"
=C#�HC%�C(  C*  C,{C.  C/�HC1��C3�C5�C8  C9��C;��C>
=C@{CB
=CD{CF{CH
=CI�CL  CN
=CO�CR  CT{CV  CW��CY�C[�HC]�HC_��Cb
=Cd
=Cf  Cg�Cj  Ck��Cn  Cp
=Cr  Ct  Cv  Cx  Cz  C|
=C~
=C�
=C�
=C�C�  C�  C�  C�  C���C�C�\C�
=C�C�C�  C�
=C�  C���C���C���C��C���C�{C�C���C�  C�C���C��C���C�  C�  C�C�\C�C���C�C��C���C�  C���C�C���C��fC��C���C���C�C�C�\C�C���C���C�  C�  C�  C���C���C�  C�C�C�
=C�C�
=C�C���C���C��C��C�  C�\C�\C�C�
=C�
=C�C�  C���C���C���C���C���C���C��C���C�  C�
=C�  C���C�  C�C��C���C�
=C�C�C�  C�  C�  C���C�C���C���C�  C�C�
=C�
=C�  C��C���C�  C�C�C�  C�
=C�\C���C��C���C�  C�
=C�  C���C���C���C��C���C���C�C���D u�D ��D}qD  D��D�D��D�D�DD��D�D� D�RDz�D��Dz�D��D	xRD	�RD
z�D�D�D�D��D�D}qD��D��D  D�D�D��D  D}qD  D�D�qD� D  Dz�D  D��D��D}qDD��D�D�D�D� D�D}qD�qD}qD�qD�DD� D  D� D�qD��D �D ��D!�D!z�D!�qD"�D#�D#z�D#�qD$� D$�qD%� D&D&��D&�qD'� D(�D(z�D(��D)� D)�qD*}qD+�D+� D,  D,� D,��D-z�D-��D.xRD.�RD/xRD/�RD0xRD1  D1�D2�D2�D2�RD3}qD4D4z�D4�qD5��D6D6� D6�qD7z�D7�RD8z�D8��D9z�D9��D:xRD:��D;z�D;�qD<��D=  D=��D>  D>� D>��D?}qD@  D@��DADA� DA�qDBz�DB�RDCz�DC�qDDz�DD��DE� DF�DF� DF�qDG}qDG�qDHz�DI�DI}qDI��DJz�DJ��DK}qDL�DL�DM�DM� DN  DN� DO  DO� DO�qDP}qDP�qDQ�DR  DRxRDR��DSz�DT  DT� DU  DU�DV�DV}qDV��DW� DW��DXz�DY  DY��DZ  DZz�DZ��D[z�D[�qD\� D]  D]��D^
=D^}qD^��D_� D`�D`�DaDa��Db  Db� Db��Dc}qDd  Ddu�Dd�RDez�De�qDfu�Df�qDg� Dh  Dh�Di
=Di��Dj�Dj� Dj�qDk� Dl�Dl��Dm�Dm��Dn�Dn��DoDo��Do�qDp��DqDq}qDq��Dr� Ds  Ds}qDt  Dt}qDu�Du��Dv�Dv� Dw  Dw� Dx�Dx�=DyDy� Dy��Dz}qD{�D{}qD{�RD|� D}�D}z�D}�RD~xRD~��D� D��D�C�D�~�D��qD�HD�AHD���D��HD��)D�=qD�~�D�� D�  D�@ D���D���D��D�AHD���D�� D�  D�AHD��HD���D��)D�>�D��HD�D�HD�>�D�� D�D�HD�=qD�� D��HD�  D�>�D�~�D���D��qD�>�D�~�D��qD���D�AHD��HD��qD�  D�AHD�� D�� D���D�@ D���D�D�HD�@ D�}qD��qD���D�>�D�� D�� D��D�@ D�� D��HD���D�AHD�� D��qD��qD�=qD�� D��HD�HD�>�D�� D�� D��D�ED��HD�� D�HD�@ D�}qD���D���D�=qD�� D��qD���D�@ D�~�D��HD�HD�@ D�~�D��HD�  D�AHD��HD�� D�  D�>�D��HD��HD��qD�>�D�� D��HD���D�C�?B�\?aG�?�  ?�\)?��R?�33?�p�?���?�G�?��?��H@   @��@z�@��@�R@(��@5@8Q�@@  @E�@L��@Tz�@\(�@aG�@h��@s33@z�H@�G�@��@���@���@���@��@�@���@�p�@�G�@��@�=q@�{@��@�
=@�
=@��H@��R@\@�ff@�=q@�{@��@�z�@�Q�@�(�@޸R@��
@�@�@�{@�33@�@�Q�A   A�\A�
A�AffAQ�A	��A(�A  A  A�\A�
AffAQ�A��A�A�RA�RA!G�A#33A%�A'
=A(��A*�HA,��A.�RA0��A1�A4z�A7
=A8Q�A:=qA<(�A>{A@  AB�\ADz�AFffAH��AJ�HAJ=qAN{AP��AQG�AU�AVffAX��AZ�HA\(�A^�RAaG�Adz�Adz�Ag
=Ag�Al��Amp�An�RAp��Ar�\Au�Aw
=AxQ�Az�HA}p�A~�RA�Q�A���A��\A��A���A�A�ffA�\)A���A���A��\A��
A���A�ffA�{A��A���A���A�=qA��
A���A�A�
=A�  A���A��A��HA��
A��A�{A�A�  A���A�=qA�33A�z�A�p�A�ffA�\)A�Q�A���A�=qA�33A�z�A�p�A�ffA�\)A�G�A��A�=qA�33A�z�A�p�A�ffA�\)A���A���A��\A��A���A�A��RA�
=A���A�G�A\A��
A��A�{AƸRAǮA���A�=qA��HA��
A��A�{A�\)A�Q�A�G�A�z�AӅA��HA�p�A�ffA׮A�Q�A�G�Aڏ\A��A�z�A�p�A޸RA�\)A��A��A��HA�A���A�A�RA�  A���A��A�\A��
A���A�{A�
=A�  A�G�A��A��HA�(�A��A�{A�\)A�Q�A���A�=qA�33A�z�A�{A�ffA�\)B ��BBG�BBffB�\B33B�
BQ�B��Bp�B�BffB
=B�B  Bz�B��B	p�B	�B
ffB
�HB\)B  Bz�B��B��B�BffB33B�B  B��B�B��B=qB�\B
=B�BQ�B��BG�BB=qB�HB33B�
B��B��B��B{B�\B�HB�B(�B��B�BB=qB�RB\)B   B Q�B ��B"ffB"{B"�\B#
=B#�B$  B$��B$��B$��B&{B%B'
=B'�B(Q�B(��B)G�B)B*=qB*�RB+\)B+�
B,Q�B,��B-�B-��B.{B.�RB/33B/�B0(�B0��B1G�B1�B2=qB2�HB3�B4  B4��B5�B5B6ffB7
=B7�B8(�B8��B9p�B:ffB:ffB:�RB;�B<(�B<��B=�B=�B>=qB>�RB?�B@(�B@��BAG�BABBffBC33BC�BD(�BD��BEG�BE�BF�\BG
=G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   ?�\)@�\@B�\@}p�@�p�@��R@�  A   A\)A#33A?\)A`  A�  A�  A�  A�  A�  A�Q�A߮A�\)B   B(�B  B�
B (�B(  B0  B8(�B?�
BG\)BPz�BX(�B_�Bh  Bp(�Bxz�B�B�B�B��
B�B��B��B�B��B��
B��B��B�  B��B��
B��B�{B�  B��
B��B�{B�{B�(�B�{B�(�B�(�B��
B��
B�{B�=qB��B�  C �C{C
=C{C{C

=C
=C
=C{C{C
=C�C(�C  C
=C{C �C"
=C#�HC%�C(  C*  C,{C.  C/�HC1��C3�C5�C8  C9��C;��C>
=C@{CB
=CD{CF{CH
=CI�CL  CN
=CO�CR  CT{CV  CW��CY�C[�HC]�HC_��Cb
=Cd
=Cf  Cg�Cj  Ck��Cn  Cp
=Cr  Ct  Cv  Cx  Cz  C|
=C~
=C�
=C�
=C�C�  C�  C�  C�  C���C�C�\C�
=C�C�C�  C�
=C�  C���C���C���C��C���C�{C�C���C�  C�C���C��C���C�  C�  C�C�\C�C���C�C��C���C�  C���C�C���C��fC��C���C���C�C�C�\C�C���C���C�  C�  C�  C���C���C�  C�C�C�
=C�C�
=C�C���C���C��C��C�  C�\C�\C�C�
=C�
=C�C�  C���C���C���C���C���C���C��C���C�  C�
=C�  C���C�  C�C��C���C�
=C�C�C�  C�  C�  C���C�C���C���C�  C�C�
=C�
=C�  C��C���C�  C�C�C�  C�
=C�\C���C��C���C�  C�
=C�  C���C���C���C��C���C���C�C���D u�D ��D}qD  D��D�D��D�D�DD��D�D� D�RDz�D��Dz�D��D	xRD	�RD
z�D�D�D�D��D�D}qD��D��D  D�D�D��D  D}qD  D�D�qD� D  Dz�D  D��D��D}qDD��D�D�D�D� D�D}qD�qD}qD�qD�DD� D  D� D�qD��D �D ��D!�D!z�D!�qD"�D#�D#z�D#�qD$� D$�qD%� D&D&��D&�qD'� D(�D(z�D(��D)� D)�qD*}qD+�D+� D,  D,� D,��D-z�D-��D.xRD.�RD/xRD/�RD0xRD1  D1�D2�D2�D2�RD3}qD4D4z�D4�qD5��D6D6� D6�qD7z�D7�RD8z�D8��D9z�D9��D:xRD:��D;z�D;�qD<��D=  D=��D>  D>� D>��D?}qD@  D@��DADA� DA�qDBz�DB�RDCz�DC�qDDz�DD��DE� DF�DF� DF�qDG}qDG�qDHz�DI�DI}qDI��DJz�DJ��DK}qDL�DL�DM�DM� DN  DN� DO  DO� DO�qDP}qDP�qDQ�DR  DRxRDR��DSz�DT  DT� DU  DU�DV�DV}qDV��DW� DW��DXz�DY  DY��DZ  DZz�DZ��D[z�D[�qD\� D]  D]��D^
=D^}qD^��D_� D`�D`�DaDa��Db  Db� Db��Dc}qDd  Ddu�Dd�RDez�De�qDfu�Df�qDg� Dh  Dh�Di
=Di��Dj�Dj� Dj�qDk� Dl�Dl��Dm�Dm��Dn�Dn��DoDo��Do�qDp��DqDq}qDq��Dr� Ds  Ds}qDt  Dt}qDu�Du��Dv�Dv� Dw  Dw� Dx�Dx�=DyDy� Dy��Dz}qD{�D{}qD{�RD|� D}�D}z�D}�RD~xRD~��D� D��D�C�D�~�D��qD�HD�AHD���D��HD��)D�=qD�~�D�� D�  D�@ D���D���D��D�AHD���D�� D�  D�AHD��HD���D��)D�>�D��HD�D�HD�>�D�� D�D�HD�=qD�� D��HD�  D�>�D�~�D���D��qD�>�D�~�D��qD���D�AHD��HD��qD�  D�AHD�� D�� D���D�@ D���D�D�HD�@ D�}qD��qD���D�>�D�� D�� D��D�@ D�� D��HD���D�AHD�� D��qD��qD�=qD�� D��HD�HD�>�D�� D�� D��D�ED��HD�� D�HD�@ D�}qD���D���D�=qD�� D��qD���D�@ D�~�D��HD�HD�@ D�~�D��HD�  D�AHD��HD�� D�  D�>�D��HD��HD��qD�>�D�� D��HD���D�C�?B�\?aG�?�  ?�\)?��R?�33?�p�?���?�G�?��?��H@   @��@z�@��@�R@(��@5@8Q�@@  @E�@L��@Tz�@\(�@aG�@h��@s33@z�H@�G�@��@���@���@���@��@�@���@�p�@�G�@��@�=q@�{@��@�
=@�
=@��H@��R@\@�ff@�=q@�{@��@�z�@�Q�@�(�@޸R@��
@�@�@�{@�33@�@�Q�A   A�\A�
A�AffAQ�A	��A(�A  A  A�\A�
AffAQ�A��A�A�RA�RA!G�A#33A%�A'
=A(��A*�HA,��A.�RA0��A1�A4z�A7
=A8Q�A:=qA<(�A>{A@  AB�\ADz�AFffAH��AJ�HAJ=qAN{AP��AQG�AU�AVffAX��AZ�HA\(�A^�RAaG�Adz�Adz�Ag
=Ag�Al��Amp�An�RAp��Ar�\Au�Aw
=AxQ�Az�HA}p�A~�RA�Q�A���A��\A��A���A�A�ffA�\)A���A���A��\A��
A���A�ffA�{A��A���A���A�=qA��
A���A�A�
=A�  A���A��A��HA��
A��A�{A�A�  A���A�=qA�33A�z�A�p�A�ffA�\)A�Q�A���A�=qA�33A�z�A�p�A�ffA�\)A�G�A��A�=qA�33A�z�A�p�A�ffA�\)A���A���A��\A��A���A�A��RA�
=A���A�G�A\A��
A��A�{AƸRAǮA���A�=qA��HA��
A��A�{A�\)A�Q�A�G�A�z�AӅA��HA�p�A�ffA׮A�Q�A�G�Aڏ\A��A�z�A�p�A޸RA�\)A��A��A��HA�A���A�A�RA�  A���A��A�\A��
A���A�{A�
=A�  A�G�A��A��HA�(�A��A�{A�\)A�Q�A���A�=qA�33A�z�A�{A�ffA�\)B ��BBG�BBffB�\B33B�
BQ�B��Bp�B�BffB
=B�B  Bz�B��B	p�B	�B
ffB
�HB\)B  Bz�B��B��B�BffB33B�B  B��B�B��B=qB�\B
=B�BQ�B��BG�BB=qB�HB33B�
B��B��B��B{B�\B�HB�B(�B��B�BB=qB�RB\)B   B Q�B ��B"ffB"{B"�\B#
=B#�B$  B$��B$��B$��B&{B%B'
=B'�B(Q�B(��B)G�B)B*=qB*�RB+\)B+�
B,Q�B,��B-�B-��B.{B.�RB/33B/�B0(�B0��B1G�B1�B2=qB2�HB3�B4  B4��B5�B5B6ffB7
=B7�B8(�B8��B9p�B:ffB:ffB:�RB;�B<(�B<��B=�B=�B>=qB>�RB?�B@(�B@��BAG�BABBffBC33BC�BD(�BD��BEG�BE�BF�\BG
=G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A
��A
z�A
v�A
~�A
r�A
^5A
ffA
JA	�TA	�mA	�-A	t�A	dZA	XA	S�A	S�A	S�A	;dA	+A��A^5A��A�AffA\)A��A�TA ��@��@�ƨ@�@���@� �@�ƨ@�ff@��@�`B@�X@�?}@���@��@�1@�E�@���@��@��@���@�{@�M�@�~�@�+@@�~�@�=q@�J@��T@���@���@�{@�J@�$�@�-@�J@��^@�hs@�%@�@��m@�C�@�C�@�"�@�Ĝ@�F@�1@�h@�@�\@���@畁@�5?@��@���@�S�@���@◍@�@��T@�x�@�O�@�O�@�O�@�7L@�/@��@��@�%@�X@�h@�7@�h@�@�@��@��@�E�@�M�@�M�@�M�@�M�@�M�@�M�@�V@�V@�V@�M�@�E�@�5?@�-@�J@�@��@��#@���@���@�^@�-@�h@�h@��@���@��T@��T@��#@�^@��@�@�O�@�G�@�?}@��@��@��/@���@��/@��/@�Ĝ@��@���@�Z@�9X@��@�b@��@�(�@�A�@�A�@�I�@�A�@� �@��@��@��m@��;@�ƨ@߶F@ߥ�@ߝ�@ߝ�@ߕ�@�t�@�C�@��@�@��y@���@ް!@�v�@�V@�-@�5?@�$�@�{@���@ݑh@�G�@�p�@���@�J@��@�E�@ޗ�@ް!@�{@�@�X@���@�Ĝ@ܛ�@�Q�@�b@�b@���@��;@��;@��;@��m@��@�I�@�I�@�9X@��@��;@۾w@۾w@�ƨ@�ƨ@���@��m@��;@��
@��
@��
@��;@��;@��
@��
@��
@��;@��
@��m@��m@���@��m@���@۾w@�|�@�"�@��y@���@ڟ�@ڗ�@��@١�@ٺ^@ٲ-@١�@ف@�G�@�%@�Ĝ@�bN@׶F@�dZ@�S�@�C�@�+@�@ָR@�@պ^@Չ7@�`B@�X@�O�@�/@��`@�A�@�l�@�;d@�
=@��y@җ�@�V@ёh@�I�@Ϯ@���@�-@��@�z�@�;d@�?}@�+@Ĭ@���@�O�@�/@��`@��@�ƨ@��\@��+@�E�@��@�/@�%@���@��@��@�ƨ@���@�\)@��y@��\@�ff@�^5@�=q@�-@���@��7@�p�@�hs@�G�@�%@�Ĝ@��j@���@�%@�V@��@��@��@���@��/@�%@�/@�hs@���@���@�-@��@�hs@���@�r�@�9X@�@���@��7@���@��#@��h@�G�@�r�@�V@�bN@��`@�%@��/@�(�@�(�@� �@�O�@��@��@��@�b@�b@� �@�1@���@�ƨ@��P@��@�dZ@��y@�ȴ@��!@��@�r�@�V@�G�@�&�@��@�Ĝ@��u@�(�@��@�C�@�M�@��T@���@�x�@�G�@��@��`@�Ĝ@���@�1@��
@��w@�|�@�C�@�33@�o@�
=@���@���@��@���@�ff@�V@��@���@�r�@�1'@�ƨ@���@��@�\)@�S�@�S�@�C�@�33@��@�o@�o@���@�ff@���@���@��7@�/@���@��@��@�Z@�1'@�  @�S�@��@��\@�@�p�@���@���@��D@� �@��@��m@��
@���@�t�@�33@�o@���@���@��H@��R@�^5@���@�@�G�@�/@�&�@���@���@���@��u@�A�@��
@��P@�ȴ@�@�/@�r�@�A�@� �@�  @��@��
@��@�;d@��@�n�@��@��@�p�@��/@��@�z�@�r�@�Q�@�9X@�1@�t�@�\)@��@��@�
=@�
=@��@��@���@���@�~�@�M�@�@���@�hs@�O�@�7L@��@��@��@�%@�%@�%@�%@���@���@���@���@��@��/@��`@��/@���@���@���@���@�Ĝ@��@���@��@��D@�z�@��@�bN@�I�@�b@�1@�(�@�  @���@��
@��
@���@��@��P@�S�@�K�@�+@��\@�M�@�M�@�5?@�$�@�{@�J@��@��T@���@��^@��@�p�@�hs@�&�@�V@�V@�V@�V@�%@�Ĝ@��9@��@���@��D@��@�r�@�r�@�r�@�j@�bN@�Z@�A�@�9X@�(�@�b@���@��@��@���@���@�|�@�l�@�\)@�K�@�"�@�
=@�
=@�
=@�@���@��y@��y@��H@��@�ȴ@���@��!@��\@��\@�~�@�v�@�^5@�^5@�^5@�V@�E�@�E�@�M�@�E�@�5?@�-@�-@�$�@���@���@���@�@�@�@��@��@��#@���@�@�@���@�`BA
��A
�A
�A
�!A
��A
��A
�A
�DA
�+A
~�A
~�A
z�A
z�A
v�A
v�A
z�A
r�A
z�A
v�A
v�A
r�A
v�A
z�A
z�A
~�A
z�A
z�A
~�A
~�A
~�A
~�A
z�A
v�A
�A
~�A
�A
�A
v�A
jA
ZA
^5A
^5A
ZA
VA
A�A
Q�A
ffA
bNA
n�A
n�A
n�A
n�A
jA
n�A
n�A
ffA
bNA
ZA
ZA
^5A
M�A
1A	�A	�A	�A	�A	�A	�A	�TA	�TA	�TA	�;A	�;A	�TA	�TA	�TA	�A	�A	�A	�A	�A	�A	�A	�A	�A	�
A	�
A	��A	ƨA	A	A	�wA	�wA	A	A	A	�wA	�wA	A	�^A	��A	�hA	�A	�A	|�A	|�A	x�A	x�A	x�A	x�A	t�A	p�A	t�A	p�A	p�A	l�A	l�A	l�A	l�A	l�A	hsA	dZA	hsA	hsA	dZA	dZA	dZA	dZA	`BA	`BA	`BA	`BA	dZA	`BA	`BA	`BA	\)A	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	S�A	XA	S�A	S�A	S�A	XA	S�A	S�A	XA	S�A	XA	S�A	O�A	S�A	S�A	S�A	S�A	S�A	S�A	O�A	S�A	O�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	O�A	S�A	S�A	S�A	O�A	O�A	O�A	O�A	K�A	C�A	C�A	?}A	;dA	;dA	;dA	?}A	;dA	7LA	7LA	33A	33A	/A	/A	/A	/A	/A	/A	+A	+A	&�A	"�A	"�A	&�A	&�A	"�A	"�A	"�A	"�A	�A	�A	oA	VA	
=A	A��A��A�A�A�A�yA�HA�A��A��A��A�9A��A��A��A�uAz�A^5AM�A(�A{A  A��A�A�mA�TA�TA�;A�#A�
A��A�^A��A��A�7A�A�Ap�Ap�AdZA`BA`BA`BA`BA`BAS�A?}A?}A33A�A��A��A�9A�9AĜAĜA��AĜAĜA��A�jA��A~�An�AbNAM�A(�A �AJA�mA��A��A�^A��A�Ax�Ap�Al�AdZAS�AG�A?}A"�A�AoAoA
=A
=A%A%AA��A�A��A��A��A��A�9A��A��A�DA�\A~�AjAM�A5?A-A-A$�AbA��A�A�TA�
A�wA�hA"�A9XA��A+A �A&�A�A ��A33A33A ��A =q@��
@��@�t�@�;d@�+@�"�@�"�@��@�o@�@�@�@�@��y@���@���@�M�@��@��7@�G�@��`@���@���G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   A
��A
z�A
v�A
~�A
r�A
^5A
ffA
JA	�TA	�mA	�-A	t�A	dZA	XA	S�A	S�A	S�A	;dA	+A��A^5A��A�AffA\)A��A�TA ��@��@�ƨ@�@���@� �@�ƨ@�ff@��@�`B@�X@�?}@���@��@�1@�E�@���@��@��@���@�{@�M�@�~�@�+@@�~�@�=q@�J@��T@���@���@�{@�J@�$�@�-@�J@��^@�hs@�%@�@��m@�C�@�C�@�"�@�Ĝ@�F@�1@�h@�@�\@���@畁@�5?@��@���@�S�@���@◍@�@��T@�x�@�O�@�O�@�O�@�7L@�/@��@��@�%@�X@�h@�7@�h@�@�@��@��@�E�@�M�@�M�@�M�@�M�@�M�@�M�@�V@�V@�V@�M�@�E�@�5?@�-@�J@�@��@��#@���@���@�^@�-@�h@�h@��@���@��T@��T@��#@�^@��@�@�O�@�G�@�?}@��@��@��/@���@��/@��/@�Ĝ@��@���@�Z@�9X@��@�b@��@�(�@�A�@�A�@�I�@�A�@� �@��@��@��m@��;@�ƨ@߶F@ߥ�@ߝ�@ߝ�@ߕ�@�t�@�C�@��@�@��y@���@ް!@�v�@�V@�-@�5?@�$�@�{@���@ݑh@�G�@�p�@���@�J@��@�E�@ޗ�@ް!@�{@�@�X@���@�Ĝ@ܛ�@�Q�@�b@�b@���@��;@��;@��;@��m@��@�I�@�I�@�9X@��@��;@۾w@۾w@�ƨ@�ƨ@���@��m@��;@��
@��
@��
@��;@��;@��
@��
@��
@��;@��
@��m@��m@���@��m@���@۾w@�|�@�"�@��y@���@ڟ�@ڗ�@��@١�@ٺ^@ٲ-@١�@ف@�G�@�%@�Ĝ@�bN@׶F@�dZ@�S�@�C�@�+@�@ָR@�@պ^@Չ7@�`B@�X@�O�@�/@��`@�A�@�l�@�;d@�
=@��y@җ�@�V@ёh@�I�@Ϯ@���@�-@��@�z�@�;d@�?}@�+@Ĭ@���@�O�@�/@��`@��@�ƨ@��\@��+@�E�@��@�/@�%@���@��@��@�ƨ@���@�\)@��y@��\@�ff@�^5@�=q@�-@���@��7@�p�@�hs@�G�@�%@�Ĝ@��j@���@�%@�V@��@��@��@���@��/@�%@�/@�hs@���@���@�-@��@�hs@���@�r�@�9X@�@���@��7@���@��#@��h@�G�@�r�@�V@�bN@��`@�%@��/@�(�@�(�@� �@�O�@��@��@��@�b@�b@� �@�1@���@�ƨ@��P@��@�dZ@��y@�ȴ@��!@��@�r�@�V@�G�@�&�@��@�Ĝ@��u@�(�@��@�C�@�M�@��T@���@�x�@�G�@��@��`@�Ĝ@���@�1@��
@��w@�|�@�C�@�33@�o@�
=@���@���@��@���@�ff@�V@��@���@�r�@�1'@�ƨ@���@��@�\)@�S�@�S�@�C�@�33@��@�o@�o@���@�ff@���@���@��7@�/@���@��@��@�Z@�1'@�  @�S�@��@��\@�@�p�@���@���@��D@� �@��@��m@��
@���@�t�@�33@�o@���@���@��H@��R@�^5@���@�@�G�@�/@�&�@���@���@���@��u@�A�@��
@��P@�ȴ@�@�/@�r�@�A�@� �@�  @��@��
@��@�;d@��@�n�@��@��@�p�@��/@��@�z�@�r�@�Q�@�9X@�1@�t�@�\)@��@��@�
=@�
=@��@��@���@���@�~�@�M�@�@���@�hs@�O�@�7L@��@��@��@�%@�%@�%@�%@���@���@���@���@��@��/@��`@��/@���@���@���@���@�Ĝ@��@���@��@��D@�z�@��@�bN@�I�@�b@�1@�(�@�  @���@��
@��
@���@��@��P@�S�@�K�@�+@��\@�M�@�M�@�5?@�$�@�{@�J@��@��T@���@��^@��@�p�@�hs@�&�@�V@�V@�V@�V@�%@�Ĝ@��9@��@���@��D@��@�r�@�r�@�r�@�j@�bN@�Z@�A�@�9X@�(�@�b@���@��@��@���@���@�|�@�l�@�\)@�K�@�"�@�
=@�
=@�
=@�@���@��y@��y@��H@��@�ȴ@���@��!@��\@��\@�~�@�v�@�^5@�^5@�^5@�V@�E�@�E�@�M�@�E�@�5?@�-@�-@�$�@���@���@���@�@�@�@��@��@��#@���@�@�@���@�`BA
��A
�A
�A
�!A
��A
��A
�A
�DA
�+A
~�A
~�A
z�A
z�A
v�A
v�A
z�A
r�A
z�A
v�A
v�A
r�A
v�A
z�A
z�A
~�A
z�A
z�A
~�A
~�A
~�A
~�A
z�A
v�A
�A
~�A
�A
�A
v�A
jA
ZA
^5A
^5A
ZA
VA
A�A
Q�A
ffA
bNA
n�A
n�A
n�A
n�A
jA
n�A
n�A
ffA
bNA
ZA
ZA
^5A
M�A
1A	�A	�A	�A	�A	�A	�A	�TA	�TA	�TA	�;A	�;A	�TA	�TA	�TA	�A	�A	�A	�A	�A	�A	�A	�A	�A	�
A	�
A	��A	ƨA	A	A	�wA	�wA	A	A	A	�wA	�wA	A	�^A	��A	�hA	�A	�A	|�A	|�A	x�A	x�A	x�A	x�A	t�A	p�A	t�A	p�A	p�A	l�A	l�A	l�A	l�A	l�A	hsA	dZA	hsA	hsA	dZA	dZA	dZA	dZA	`BA	`BA	`BA	`BA	dZA	`BA	`BA	`BA	\)A	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	XA	S�A	XA	S�A	S�A	S�A	XA	S�A	S�A	XA	S�A	XA	S�A	O�A	S�A	S�A	S�A	S�A	S�A	S�A	O�A	S�A	O�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	S�A	O�A	S�A	S�A	S�A	O�A	O�A	O�A	O�A	K�A	C�A	C�A	?}A	;dA	;dA	;dA	?}A	;dA	7LA	7LA	33A	33A	/A	/A	/A	/A	/A	/A	+A	+A	&�A	"�A	"�A	&�A	&�A	"�A	"�A	"�A	"�A	�A	�A	oA	VA	
=A	A��A��A�A�A�A�yA�HA�A��A��A��A�9A��A��A��A�uAz�A^5AM�A(�A{A  A��A�A�mA�TA�TA�;A�#A�
A��A�^A��A��A�7A�A�Ap�Ap�AdZA`BA`BA`BA`BA`BAS�A?}A?}A33A�A��A��A�9A�9AĜAĜA��AĜAĜA��A�jA��A~�An�AbNAM�A(�A �AJA�mA��A��A�^A��A�Ax�Ap�Al�AdZAS�AG�A?}A"�A�AoAoA
=A
=A%A%AA��A�A��A��A��A��A�9A��A��A�DA�\A~�AjAM�A5?A-A-A$�AbA��A�A�TA�
A�wA�hA"�A9XA��A+A �A&�A�A ��A33A33A ��A =q@��
@��@�t�@�;d@�+@�"�@�"�@��@�o@�@�@�@�@��y@���@���@�M�@��@��7@�G�@��`@���@���G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�
B�B�B�
B�B�B�5B�;B�NB�mB�mB�sB��B��B��B+B	7B
=B
=BVBhB{B�B�B�B{BoBoBuBoBuB{B�B�B�B�B�B�B!�B �B �B!�B#�B(�B)�B+B+B+B+B)�B'�B&�B"�B!�B$�B'�B33BB�BF�B6FB-B(�B"�B�B�B\B	7B+BBBBBBBBBBBBBB%B	7B
=BJBVBVB\BbBoBoBoBoBuBuBuBuB{B{B{B{B{B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B!�B!�B"�B!�B"�B"�B#�B#�B#�B$�B%�B%�B%�B&�B'�B'�B(�B)�B)�B+B,B,B-B-B/B/B/B/B/B.B.B.B.B.B/B2-B49B5?B6FB7LB7LB7LB8RB8RB8RB8RB8RB8RB8RB7LB7LB7LB7LB7LB7LB7LB7LB6FB6FB7LB8RB8RB8RB8RB8RB7LB6FB6FB49B1'B.B'�B$�B�B�BJBB��B�B�B�B�B�B�B�B�sB�fB�`B�`B�ZB�ZB�TB�NB�NB�HB�HB�HB�HB�HB�HB�BB�HB�NB�NB�TB�ZB�TB�ZB�`B�fB�mB�sB�sB�sB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�ZB�TB�ZB�ZB�`B�`B�TB�;B�
B��B�B�
B�B��B��B��B�#B�`B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��BB+B+B+B1B1B+B+B+B	7B	7B	7B1B	7B	7B	7B	7B	7B	7B	7B	7B1B1B1B1B1B1B+B+B+B+B+B%B%B%B%B%B%B%B%B%B%B%B%B%B%BBBBBBBBBBBBBBBBBBBB  B  B  B  B  B��B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�sB�sB�sB�sB�sB�sB�sB�sB�sB�mB�mB�mB�mB�mB�mB�mB�mB�fB�mB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�TB�NB�TB�TB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�HB�HB�HB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B��B��B�B��B�B��B�B�
B�
B�
B�
B�B�
B�B�
B�B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�B�
B�
B�B�
B�
B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�
B�B�B�B�B�
B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�#B�B�#B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�/B�/B�5B�5B�;B�;B�;B�5B�;B�;B�BB�;B�;B�;B�;B�5B�;B�BBŢB�BB�HB�TB�ZB�TB�NB�NB�NB�NB�HB�HB�BB�BB�NB�NB�TB�NB�ZB�TB�TB�ZB�`B�`B�`B�`B�mB�mB�fB�fB�fB�mB�mB�mB�sB�yB�sB�sB�sB�sB�mB�mB�mB�mB�mB�sB�mB�mB�mB�mB�sB�mB�mB�sB�mB�sB�yB�sB�yB�mB�mB�mB�sB�mB�mB�fB�`B�`B�fB�B��B��B��B�B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B%B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111111111111111111111111111111111111111111111111141111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   BНB��B��B�8B�NB�$BѬBѰB��BєBђB�$B�B�
B�B�B�TB�HB״B��B�<B��B�fB�BB�1B�B�uB��B�
BBB
�B
B
�B@BiB�B�B�B�BBoBB�BBB�BB@BaB�B�B�B(B"B! B!B!�B#�B(�B)�B*�B+;B+|B+vB*�B(�B(B#�B!�B$�B%iB.�BB=BJB:B-�B+]B$�B�B.BvB

B BmB�BOB�BFBBB,BB,BBB�B�B	8B
0BBJBBB'B^BlBlBoBqBsBfBvBxB�B�B�B�B�B�B�B�B�B�B�B�B�B�BdBRBwB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B~BtB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B=BBZBBAB�BwBDB �B"�B JB \BJB�B�BB�B�B�B�B�B�B�BjBuB!�B!�B#B"#B#B"�B#�B#�B#�B$�B%�B%�B%�B&�B'�B'�B) B)�B)�B*�B,B+�B-B,�B/3B/DB/9B/B/�B.oB.WB.KB.-B.�B/�B2B4FB5\B6yB7�B7�B7�B8�B9KB8�B8mB8iB8zB8�B7�B8SB7�B7�B7�B7YB7YB7~B6�B7GB8�B8�B8�B8�B8�B8�B8�B8LB77B5PB2lB/�B)B&�B!�B�BBB�B��B� B��B�YB�pB�B��B�B��B�B�B��B�B��B�B�B��B��B�B�VB�tB�dB��B�B�rB�`B�B�B�B�iB�8B�$B�_B�XB�}B�tB�B�B�DB�CB�5B�BB�QB�#B��B��B�B�>B�B�B�B�rB�;B�B��B��B�B�SB��B�dB��B�RB�B��B��B�+BتB��B�dB�cB�B�B��B�.B�B�
B��B��B�`B��B��B�B��B��B�BTBKB�B�B�B�B�B�B	�B	�B	~B�B	�B	�B	kB	�B
B	�B	cB	�B�BLBaB>BIB2B`B�BBJB�B�B�B�B�BdBNBfB2B(B?B<BMB2B)B�B�B�B�B\B�B{B�B]B_BZBkBBsB�B6B�B�B@B tB �B HB B B�fB 5B�`B�0B�B��B�B�4B��B��B�OB��B�B��B�*B�5B�'B�B�iB��B�WB�#B�B�B��B�B��B��B��B��B�8B�1B�4B��B�B�NB��B�B��B��B�B��B��B��B��B��B��B�B�B�B��B�B�B��B��B��B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�|B�B�B�zB�B�B��B�B�PB�B�B�B�B�B�B�B��B�B�B�YB��B�zB�B�B�B�B�B�B�B�B��B�B�B��B�B�mB�kB�mB�wB��B�B�sB�uB�B�vB�B�gB�fB�pB�lB�mB�B�oB�zB�B�B��B�\B�rB�hB�B�vB�sB�vB�B�B�XB�WB�`B�aB�lB�WB�bB�aB�lB�cB�oB�B�VB�mB�bB�vB�NB�PB�^B�gB�MB�AB�[B�hB�ZB�NB�\B�B�PB�LB�AB�NB�NB�fB�\B�hB�YB�fB�JB�GB��B�aB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B��B��B�B��B�B��B�B�
B�
B�
B�
B�B�
B�B�
B�B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�
B�B�
B�
B�B�
B�
B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�G�O�B�B�B�B�B�
B�B�B�B�B�
B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�#B�B�#B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�/B�/B�5B�5B�;B�;B�;B�5B�;B�;B�BB�;B�;B�;B�;B�5B�;B�BG�O�B�BB�HB�TB�ZB�TB�NB�NB�NB�NB�HB�HB�BB�BB�NB�NB�TB�NB�ZB�TB�TB�ZB�`B�`B�`B�`B�mB�mB�fB�fB�fB�mB�mB�mB�sB�yB�sB�sB�sB�sB�mB�mB�mB�mB�mB�sB�mB�mB�mB�mB�sB�mB�mB�sB�mB�sB�yB�sB�yB�mB�mB�mB�sB�mB�mB�fB�`B�`B�fB�B��B��B��B�B�B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B%B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111111111111111111111111111111111111111111111111141111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                   <$E<#��<#�{<#�<#��<#�&<&Z�<$T�<#�$<$4e<$2G<#ޫ<#��<#�{<#�<#�<#�&<#��<$/%<&Gi<'Qf<&)�<'x�<*�f<&A�<+H<Z5�<1�<7��<@5+<-͝<$ub<$�<&��<$�2<$�<#�i<#��<#�H<$�<$�V<(��<$	<#�X<#��<#�o<#��<#�<#�4<#�C<#�*<#�^<#�N<#�<#�<#�<#�<#�<#�<#��<#�i<#��<$p<$ <$f<$'<$�<$� <#ף<#׎<(��<1K�<#�<,��<.�Y<$Y�<(I�<'|<'�<&�<'A><$^�<$aD<#�<$m,<#�e<$�<#�<#�
<#�<#�r<#�<#�*<#�<#�<#�N<#�<#�<#׎<#��<#�{<#�l<#�<#�<#��<#�<#�<#�
<#�<#�<#׺<#�<#�<#ף<#��<#��<#��<#�J<#�i<#�o<#�]<#�i<#׎<#��<#��<#ޫ<#�<#��<#�<#�o<#�<#׺<#ޫ<#��<#��<#�<#��<#��<#ޫ<#�<#�<#�{<#�<<#�
<#�8<#�8<#�$<#��<#ߜ<#�<#�{<#��<#�c<#�r<#�<#�<<#�X<#ޫ<#�<#�<#ף<#�i<#ۮ<#�<#��<#׎<#�<#؄<#�<#��<#�<#�r<#��<#�l<#ޫ<#�<#�^<#�<#�{<#��<#�D<#��<$�<#�Q<#��<$�<#�)<#��<#�<$ <#ٛ<$t <$�<$#(<$
<#�<#�<#��<#��<#�<<#��<#��<#�<#�<#׺<#��<#�M<#�<#�D<#�^<#�<#�^<#�<#�I<#�<#�I<#��<#�0<#�{<#�<#�<#��<#�<#�X<#�<#�<#�i<#�&<#ا<#�<#��<#��<#�+<#��<#��<$<#�N<#�<#�E<#��<$B�<$6�<#��<#׎<#ٛ<#��<#��<#��<#�	<$�<$��<$.<#�D<#ا<#��<#�&<$�<$��<#�N<#�<#�e<#׎<#׎<#ޫ<$p<$�<%�<#�<#�M<#��<$<<#��<%s<'<$��<$��<%'<%�<$�h<&�J<*�<+:<.p<*ٜ<'��<#�<#�N<#؄<%��<&y�<#��<#��<$��<$	�<#��<#�<$<<$/%<$v<#��<#��<$1:<$r<#��<#ף<#��<#�o<$%<#��<#��<#�<#��<#��<#�	<#׺<#��<#�U<#ף<#�D<#�X<#�<#��<#��<#�<#�U<#�<#�<#�!<$�<#�<$�.<$o�<$Z<$�<+�)<$C�<#��<#��<#�5<#��<$a<%X
<+'�<*w<$	�<#��<#��<$��<#�<#�<<&J$<(��<%��<#ۮ<#�<#�<#�<#ۮ<$.<#�r<#�<#�o<#��<$/%<#�&<#��<$0.<&�2<$�	<#�<#�+<#�*<$�<#�<$'<$J�<$:�<%m�<$*<#�"<#�l<#�<#�M<#�4<#�J<#�<$c�<#�<#��<#�<#�N<#�D<#�<#׎<#��<#�<#ߜ<#�<#�<#��<#�H<&L�<$aD<#�H<$*<#�&<#�+<#��<#׎<#�<#�<#ا<#��<#׎<#�<#�Q<$*<$+<$	�<#�e<$�<#��<#��<#�&<#��<#��<#�<$�e<#�<$m,<$�!<$<$L<#�N<$ <$0.<#��<#׎<#�<#��<#ߜ<#��<#��<#�8<#�<#��<#�<$�<$.<#�<$7�<#ܯ<#׺<#�<#�<#�<#ٛ<$r<$&<$�<%*<$�k<%$<$ʾ<#�<#ߜ<#�<#�D<#��<$�<$ �<$v<$j|<$e.<$#(<#�8<$o�<#�<#�l<#�c<#ޫ<#�l<#�<$z�<#�<#�<#�<#��<#�<#��<#ۮ<#�$<#�<#�J<#�4<$G<$�<#�Q<#�l<#��<#��<#�<#��<#��<#�<#�
<#�
<#�i<#�<#�
<#�<#�D<#��<#�{<#�i<#�i<#�
<#�X<#�0<#ף<#�r<#�i<#�I<#��<#��<#�i<#�^<#�+<#�"<#��<#��<#�<#��<#��<#�<#׎<#ޫ<#��<#�m<#�C<#�N<$x+<#�)<#�0<#�8<#��<#��<#׺<#�r<#��<#��<#�l<#��<#ٛ<#�$<#�)<#�8<#�
<#�<#�
<#��<#�<#��<#׎<#׺<#ۮ<#��<#�<#�<#�
<#��<#�{<#׎<#�r<#׺<#�<#�r<#�+<#��<#�<#��<#ף<#ۮ<#�o<#��<#�o<#�<#�r<#�<#�<#�{<#׎<#��<#�<#ף<#׎<#��<#׺<#�D<#�^<#�<#��<#ף<#ڑ<#�
<#�<#�X<#��<#�<#׎<#׎<#�<#�{<#�
<#ף<#�<#�<#�<#׎<#�
<#�
<#��<#ף<#�<#�i<#��<#�<#�<$�<#��<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                    WHOI    ARSQ    WHQC    V0.5                                                                                                                                              20200527000000    CF                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                