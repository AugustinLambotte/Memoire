CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  B   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92018-11-06T01:28:56Z creation; 2023-02-09T14:06:17Z DMQC;      
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
resolution        =���   axis      Z          <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  N�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       S@   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  eP   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       i�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       {�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       �x   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       �   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       �   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �,   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       Ѱ   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       �D   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  �T   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                      HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                      HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                      HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                       HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � (   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Argo profile    3.1 1.2 19500101000000  20181106012856  20230209090617  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               M   MAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�q�I�@�q�I�11  @�q�O�`@�q�O�`@N���	k�@N���	k��<LD|0��<LD|0�11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�\)?��H@@  @�  @�  @�G�@�  A   A��A#33A@  A`  A\)A�  A�  A��A��A�  A�  A�Q�B (�B�
B  B(�B (�B(  B/�
B8  B@(�BG\)BPQ�BX  B_�Bh(�Bp  Bw�
B�  B�{B�(�B�  B�{B�(�B��B�{B�{B�{B�{B��B�  B�{B�{B�(�B�{B�  B��B��
B��B��B��
B��B�  B�  B�(�B�{B�  B�  B��
B�  C {C  C  C��C  C

=C  C  C  C
=C  C��C��C�C�C  C {C"  C#��C%�C(  C*{C,{C.{C0�C2�C4�C6
=C8  C:  C<  C=��C?��CA��CD  CF
=CH
=CJ
=CL
=CN{CP  CQ�HCS�HCU�HCW��CZ
=C\�C^  C_�HCa��Cd
=Ce��Cg�HCi��Cl  Cn
=Cp{Cr  Cs�Cv  Cx{Cz  C{��C~
=C��C��C���C�
=C�C�  C�
=C�  C���C�  C�C�\C�C��C���C���C���C�  C�C�
=C���C��C��C��C��C��C�  C�
=C�  C���C���C�  C�C�  C�C�C���C���C�  C�  C�C�  C�  C���C�  C�C���C�  C���C���C�C�C�
=C�  C���C���C���C���C�  C�C�C�  C���C���C���C���C�C�  C���C�  C�C�  C�C�C�  C�  C�C�C�  C�  C���C���C���C���C���C�  C�  C�C�  C�  C�  C�  C�  C�  C�  C�  C�C���C�C�  C���C�C�C���C�  C�
=C�C�C�C�
=C�  C���C�
=C�C�  C�
=C�
=C�  C��C���C�C���C���C���C�  C�  C��C��C���D }qD  D��D�qDxRD�qD��D�qDz�D  D�D�D� D�qDz�D�D� D�qD	z�D	�RD
� DD� D��D}qD�D�D�D� D�RDz�D��DxRD  D�D�D� D  D� D  D� D�D� D  D}qD��D� D�D��D  Dz�D�RD� DD��D�D��D  Dz�D  D�D  D}qD   D }qD �qD!}qD!�qD"� D"�qD#�D$  D$� D%�D%��D%�qD&}qD'  D'}qD'�qD(}qD)�D)�D*  D*� D+�D+}qD+�qD,}qD,�qD-��D.�D.��D/�D/� D0  D0}qD1�D1��D1�qD2xRD2�qD3� D3��D4}qD5�D5� D6�D6�D7D7��D8�D8� D8�qD9��D:D:�D;�D;�D<  D<z�D<�qD=}qD=��D>� D?  D?z�D@  D@�DA�DA��DB  DB}qDB��DC}qDDDD��DE  DEz�DF  DF��DF�qDG� DH  DH� DI�DI��DJ  DJz�DJ�qDK� DK�qDL}qDM  DM� DM�qDN��DO�DOz�DO��DP� DP�qDQz�DQ��DR}qDS  DSz�DS�qDT� DU  DU��DV�DV��DV�qDW}qDW��DXxRDX�qDY��DZDZ� DZ��D[z�D\  D\�D]  D]��D^  D^}qD_  D_}qD_�qD`�Da  Da� Db�Db�DcDc��Dd�Dd� De  De�De�qDfz�Df�qDg� Dg�qDh}qDh�qDi� Dj  Dj�Dk  Dk� Dl  Dl}qDm�Dm�Dn�Dn�DoDo� Do�qDp}qDp�qDq��Dr  Dr}qDr�qDs��Dt�Dt� DuDu��Dv  Dv}qDv�qDw��Dx�Dx� Dx�qDy}qDz  Dz��D{  D{� D|�D|}qD}  D}�D~  D~}qD~�qD��D�HD�AHD�� D���D��qD�>�D��HD���D��)D�@ D��HD�� D�  D�=qD�� D�D�  D�=qD�~�D�� D�  D�@ D�~�D�� D�HD�@ D�~�D���D���D�@ D�� D���D���D�AHD�� D�� D�  D�@ D��HD�D�HD�@ D�� D��HD�  D�>�D��HD�D�HD�@ D��HD�D�  D�>�D��HD�� D�  D�@ D��HD��HD��D�+�?W
=?u?�=q?���?���?�Q�?���?�(�?��?��H@�@\)@��@!G�@(��@.{@8Q�@@  @G�@O\)@W
=@^�R@fff@n{@u@}p�@��\@�ff@�=q@�{@��@�z�@���@�p�@�G�@��
@��@��@�{@��@�@���@�p�@\@��@Ǯ@˅@�\)@�33@�Q�@��H@�p�@�G�@��@���@���@��@�z�@�Q�@�(�@��RAG�A33A�AffA�A��A�Ap�A�RA  A�A�
AA
=A��A=qA��A{A   A!G�A#�
A%�A'
=A(��A*�HA,(�A.{A0  A1�A3�
A5A7
=A8��A:�HA<��A>�RA@  AA�AC�
AEAG
=AH��AJ�HAL��AN{AP  AQ�AS�
AU�AW
=AX��AZ�HA\��A^{A`  Aa�Ac�
Ae�Ag
=Ai��Aj�HAl��An{Ap��Ar�\As�
AuAw�Ay��Az�HA|��A\)A�Q�A���A�=qA�33A�(�A���A�p�A��RA��A�Q�A�G�A��\A��A��
A��A�{A�\)A��A���A��A��HA��A�z�A�A��RA�\)A�Q�A���A��\A�33A�(�A��A�{A��RA�  A���A��A��\A��A���A�A�ffA�\)A�Q�A�G�A�=qA�33A��
A��A�A��RA��A���A���A��\A��A�z�A�p�A�{A�
=A�  A���A���A��HA��
A���A�p�A�ffA�\)A�Q�A�G�A��A��HA��
A���A�AƸRA�\)A�Q�A���A�=qA�33A��
A���A�{AθRA�\)A�Q�A�G�A�=qA�33A�(�A��A�{AָRA׮Aأ�Aٙ�Aڏ\AۅA�z�A�p�A�{A�
=A�Q�A�G�A��A��HA��
A���A�p�A�ffA�A�Q�A�G�A�=qA�33A�(�A���A�{A�
=A�  A��A�A�\A�A�(�A�p�A�ffA�
=A�  A���A��A��HA��
A���A�A�ffA�\)B (�B ��B ��Bp�B�B=qB�RB\)B�
BQ�B��B�B��B{B�\B
=B�B�
BQ�B��B	p�B	B
=qB
�RB33B�B(�B��B�Bp�B{B�\B�HB\)B�
BQ�B��B�BB{B�\B
=B�B�
BQ�B��B�B��B{B�\B�HB33B�
B(�Bz�B��Bp�B�B=qB�RB33B�B  Bz�B��BG�BB=qB�RB
=B�B   B Q�B ��B!G�B!B"{B"�\B#
=B#\)B#�
B$Q�B$��B%�B%��B&{B&�\B&�HB'\)B'�
B((�B(��B)G�B)��B)�B*ffB*�HB+33B+�B,(�B,��B,��B-p�B-�B.ffB.�RB/\)B/�
B0  B0��B1�B1p�B1B2ffB2�HB333B3�B4(�B4z�B5�B5��B5�B6ffB6�HB7\)B7�
B8Q�B8��B9�B9��B:{B:�\B;
=B;\)B;�
B<z�B<��B=G�B=B>=qB>�\B?33B?�B@  B@z�BA�BAp�BA�BB�\BB�HBC\)BD  BDQ�BD��BEp�BEBF=qBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  ?�\)?��H@@  @�  @�  @�G�@�  A   A��A#33A@  A`  A\)A�  A�  A��A��A�  A�  A�Q�B (�B�
B  B(�B (�B(  B/�
B8  B@(�BG\)BPQ�BX  B_�Bh(�Bp  Bw�
B�  B�{B�(�B�  B�{B�(�B��B�{B�{B�{B�{B��B�  B�{B�{B�(�B�{B�  B��B��
B��B��B��
B��B�  B�  B�(�B�{B�  B�  B��
B�  C {C  C  C��C  C

=C  C  C  C
=C  C��C��C�C�C  C {C"  C#��C%�C(  C*{C,{C.{C0�C2�C4�C6
=C8  C:  C<  C=��C?��CA��CD  CF
=CH
=CJ
=CL
=CN{CP  CQ�HCS�HCU�HCW��CZ
=C\�C^  C_�HCa��Cd
=Ce��Cg�HCi��Cl  Cn
=Cp{Cr  Cs�Cv  Cx{Cz  C{��C~
=C��C��C���C�
=C�C�  C�
=C�  C���C�  C�C�\C�C��C���C���C���C�  C�C�
=C���C��C��C��C��C��C�  C�
=C�  C���C���C�  C�C�  C�C�C���C���C�  C�  C�C�  C�  C���C�  C�C���C�  C���C���C�C�C�
=C�  C���C���C���C���C�  C�C�C�  C���C���C���C���C�C�  C���C�  C�C�  C�C�C�  C�  C�C�C�  C�  C���C���C���C���C���C�  C�  C�C�  C�  C�  C�  C�  C�  C�  C�  C�C���C�C�  C���C�C�C���C�  C�
=C�C�C�C�
=C�  C���C�
=C�C�  C�
=C�
=C�  C��C���C�C���C���C���C�  C�  C��C��C���D }qD  D��D�qDxRD�qD��D�qDz�D  D�D�D� D�qDz�D�D� D�qD	z�D	�RD
� DD� D��D}qD�D�D�D� D�RDz�D��DxRD  D�D�D� D  D� D  D� D�D� D  D}qD��D� D�D��D  Dz�D�RD� DD��D�D��D  Dz�D  D�D  D}qD   D }qD �qD!}qD!�qD"� D"�qD#�D$  D$� D%�D%��D%�qD&}qD'  D'}qD'�qD(}qD)�D)�D*  D*� D+�D+}qD+�qD,}qD,�qD-��D.�D.��D/�D/� D0  D0}qD1�D1��D1�qD2xRD2�qD3� D3��D4}qD5�D5� D6�D6�D7D7��D8�D8� D8�qD9��D:D:�D;�D;�D<  D<z�D<�qD=}qD=��D>� D?  D?z�D@  D@�DA�DA��DB  DB}qDB��DC}qDDDD��DE  DEz�DF  DF��DF�qDG� DH  DH� DI�DI��DJ  DJz�DJ�qDK� DK�qDL}qDM  DM� DM�qDN��DO�DOz�DO��DP� DP�qDQz�DQ��DR}qDS  DSz�DS�qDT� DU  DU��DV�DV��DV�qDW}qDW��DXxRDX�qDY��DZDZ� DZ��D[z�D\  D\�D]  D]��D^  D^}qD_  D_}qD_�qD`�Da  Da� Db�Db�DcDc��Dd�Dd� De  De�De�qDfz�Df�qDg� Dg�qDh}qDh�qDi� Dj  Dj�Dk  Dk� Dl  Dl}qDm�Dm�Dn�Dn�DoDo� Do�qDp}qDp�qDq��Dr  Dr}qDr�qDs��Dt�Dt� DuDu��Dv  Dv}qDv�qDw��Dx�Dx� Dx�qDy}qDz  Dz��D{  D{� D|�D|}qD}  D}�D~  D~}qD~�qD��D�HD�AHD�� D���D��qD�>�D��HD���D��)D�@ D��HD�� D�  D�=qD�� D�D�  D�=qD�~�D�� D�  D�@ D�~�D�� D�HD�@ D�~�D���D���D�@ D�� D���D���D�AHD�� D�� D�  D�@ D��HD�D�HD�@ D�� D��HD�  D�>�D��HD�D�HD�@ D��HD�D�  D�>�D��HD�� D�  D�@ D��HD��HD��D�+�?W
=?u?�=q?���?���?�Q�?���?�(�?��?��H@�@\)@��@!G�@(��@.{@8Q�@@  @G�@O\)@W
=@^�R@fff@n{@u@}p�@��\@�ff@�=q@�{@��@�z�@���@�p�@�G�@��
@��@��@�{@��@�@���@�p�@\@��@Ǯ@˅@�\)@�33@�Q�@��H@�p�@�G�@��@���@���@��@�z�@�Q�@�(�@��RAG�A33A�AffA�A��A�Ap�A�RA  A�A�
AA
=A��A=qA��A{A   A!G�A#�
A%�A'
=A(��A*�HA,(�A.{A0  A1�A3�
A5A7
=A8��A:�HA<��A>�RA@  AA�AC�
AEAG
=AH��AJ�HAL��AN{AP  AQ�AS�
AU�AW
=AX��AZ�HA\��A^{A`  Aa�Ac�
Ae�Ag
=Ai��Aj�HAl��An{Ap��Ar�\As�
AuAw�Ay��Az�HA|��A\)A�Q�A���A�=qA�33A�(�A���A�p�A��RA��A�Q�A�G�A��\A��A��
A��A�{A�\)A��A���A��A��HA��A�z�A�A��RA�\)A�Q�A���A��\A�33A�(�A��A�{A��RA�  A���A��A��\A��A���A�A�ffA�\)A�Q�A�G�A�=qA�33A��
A��A�A��RA��A���A���A��\A��A�z�A�p�A�{A�
=A�  A���A���A��HA��
A���A�p�A�ffA�\)A�Q�A�G�A��A��HA��
A���A�AƸRA�\)A�Q�A���A�=qA�33A��
A���A�{AθRA�\)A�Q�A�G�A�=qA�33A�(�A��A�{AָRA׮Aأ�Aٙ�Aڏ\AۅA�z�A�p�A�{A�
=A�Q�A�G�A��A��HA��
A���A�p�A�ffA�A�Q�A�G�A�=qA�33A�(�A���A�{A�
=A�  A��A�A�\A�A�(�A�p�A�ffA�
=A�  A���A��A��HA��
A���A�A�ffA�\)B (�B ��B ��Bp�B�B=qB�RB\)B�
BQ�B��B�B��B{B�\B
=B�B�
BQ�B��B	p�B	B
=qB
�RB33B�B(�B��B�Bp�B{B�\B�HB\)B�
BQ�B��B�BB{B�\B
=B�B�
BQ�B��B�B��B{B�\B�HB33B�
B(�Bz�B��Bp�B�B=qB�RB33B�B  Bz�B��BG�BB=qB�RB
=B�B   B Q�B ��B!G�B!B"{B"�\B#
=B#\)B#�
B$Q�B$��B%�B%��B&{B&�\B&�HB'\)B'�
B((�B(��B)G�B)��B)�B*ffB*�HB+33B+�B,(�B,��B,��B-p�B-�B.ffB.�RB/\)B/�
B0  B0��B1�B1p�B1B2ffB2�HB333B3�B4(�B4z�B5�B5��B5�B6ffB6�HB7\)B7�
B8Q�B8��B9�B9��B:{B:�\B;
=B;\)B;�
B<z�B<��B=G�B=B>=qB>�\B?33B?�B@  B@z�BA�BAp�BA�BB�\BB�HBC\)BD  BDQ�BD��BEp�BEBF=qBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�AZAVAQ�AQ�AVAM�AI�A5?A�
A�AoA��A�7A��A(�AA��A�hA�AdZAVA
�uA
I�A
I�A
1A	��A�HAQ�A1AdZAn�A�#A��A �@��-@�;d@���@��@�M�@���@��u@���@�\@�p�@�u@�=q@�  @�!@�=q@�^@�O�@�&�@�V@���@�D@�F@�\)@�v�@�@���@���@��#@�&�@�9@�j@�Q�@�9X@��@�t�@��@�5?@�-@�hs@��@�j@�j@�9X@ߕ�@�dZ@�33@��@ޏ\@�@�&�@�r�@��m@���@ۥ�@ە�@�C�@ڸR@ڟ�@���@��@��@���@�`B@��@ؓu@�Z@�9X@�b@��
@׶F@ׅ@ׅ@ׅ@�\)@�C�@�\)@ם�@��@أ�@��@��@��`@ؼj@�r�@�9X@��@׶F@�z�@���@��@ٙ�@�E�@ڗ�@�
=@�33@�S�@���@܃@ܓu@ܬ@ܬ@ܓu@�z�@�z�@܃@܋D@ܓu@ܛ�@ܓu@܋D@܋D@܃@�bN@�I�@�b@��@��@��;@�C�@���@ڟ�@��#@�x�@٩�@ٺ^@ٺ^@���@ٺ^@ٺ^@���@��#@��@�{@�J@��@���@ف@�%@ج@�l�@�~�@�@�`B@�~�@��
@أ�@�7L@���@�V@ڗ�@���@�o@��H@ڰ!@ڏ\@��@��T@�x�@��@�r�@׍P@�~�@պ^@�O�@�V@ԛ�@�r�@�bN@�A�@�(�@��;@�K�@ҏ\@�M�@�{@�/@�;d@��@Ώ\@�J@���@�@�@ͩ�@͉7@́@�p�@�X@��@��@���@�Ĝ@�Ĝ@�Ĝ@���@̴9@�Q�@��@�ƨ@˅@�;d@�"�@�@�-@�&�@��@��/@���@ȴ9@ȴ9@Ȭ@ȣ�@ȋD@�z�@�j@�Z@� �@��m@Ǯ@�\)@ƸR@�@�ƨ@ț�@��`@��`@ȼj@Ȭ@ȣ�@ȴ9@ȼj@���@�&�@���@ȣ�@ȓu@ȋD@�I�@�1@��
@ǍP@�S�@�33@�@��@Ƈ+@�V@�$�@��@���@�x�@��@�(�@�t�@�@�{@���@�(�@���@�S�@�$�@��@��^@��@�dZ@�^5@��@��#@���@�%@�A�@�"�@��y@�ȴ@�v�@�n�@�-@�@���@�J@��@�E�@��\@��!@��@�\)@�dZ@�\)@�C�@��@���@���@�n�@�J@��@��@���@�?}@���@�A�@�  @�ƨ@�l�@��\@�^5@�E�@�E�@�-@�5?@�-@�=q@�5?@�J@��T@���@��@��!@���@�@��-@���@�1'@�|�@��!@�J@��@�r�@�1'@��;@��F@��P@�S�@�33@�"�@���@���@���@�v�@�$�@��@�{@�@��T@���@���@���@��#@���@�@��-@���@�p�@���@��j@� �@���@��H@���@��\@�~�@�~�@�~�@�n�@�ff@�ff@�ff@�M�@�-@�@�p�@��@��`@���@�z�@�j@�I�@� �@��@�b@�1@�1@�1@���@��m@��;@��
@��w@���@�|�@�;d@�33@�"�@�o@�
=@�
=@�@��H@��R@�~�@�=q@�-@�-@�-@�$�@��@��@�{@�J@���@��#@��@�/@���@��9@��@��@��@���@��@��@��9@��9@��9@��j@�Ĝ@�Ĝ@��9@�9X@��@�|�@�;d@�ff@�@���@���@���@��@�x�@�x�@�hs@�p�@�hs@�`B@�`B@�`B@�X@�O�@�O�@�G�@�O�@�/@�&�@�&�@��@�%@�V@��9@��u@��@�I�@�1@���@���@��;@���@��w@��@��@��@���@���@�|�@�l�@�\)@�l�@�l�@�l�@�dZ@�dZ@�\)@�dZ@�dZ@�\)@�K�@�;d@�S�@�33@��y@���@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@��!@���@�v�@�V@�E�@�{@�@��@���@��^@��7@�`B@�%@��D@�Q�@�1'@��@��
@�+@�o@�o@���@��H@��R@���@��\@�n�@�^5@�^5@�V@�V@�J@���@���@��7@�`B@�X@�X@�X@�X@�`B@�`B@�`B@�`B@�`B@�`B@�`B@�`B@�hs@�hsAZAZAZAVAVAVAVAQ�AVAVAVAVA^5AVAVAQ�AZAVAM�AI�AM�AM�AM�AM�AM�AQ�AQ�AVAVAQ�AVAVAVAQ�AZAVAZAZAZAVAZAVAM�AI�AM�AI�AE�AI�AI�AM�AI�AI�AE�A5?AQ�AVA(�AE�AZAbNAQ�A(�A�A1'A(�AJA�A�AƨAƨA��A�wA��A�A�^A�
AAA�#AƨA|�Ap�A`BAXAO�AK�AG�AG�AC�A?}A7LA7LA33A/A+A"�AoAVA%AA��A��A�A�A�A�`A�HA�A�A��A��A��AȴA�jA�jA�jA�RA�RA�A�AbNAE�A1'A�AA�mA�#A�
A��AA�wA�-A��A�Ap�AXA7LA/A+A&�A&�A"�A"�A"�A�AoA
=AA��A�yA�/A��A��AĜA��Az�Ar�AjAVAI�AE�AA�A9XA5?A1'A-A(�A$�A �A�A �A �A�A�A�A�A�A{AbAbAbA1A1A1AA  A  A��A��A��A�A�A�A�A�A�mA�mA�#A�
A�
A��A��AƨAA�wA�wA�^A�FA�FA�FA�-A��A��A��A�hA�hA�PA�PA�PA�PA�PA�PA�PA�PA�hA�PA�PA�PA�PA�7A�7A�A�7A�A�A�A�A�A|�A|�A|�A�A�A�A�A|�A|�A|�At�At�At�Ap�At�Ap�AhsAhsA`BA`BA\)AXAO�AC�A;dA/A+A"�A�A�A�A�A�AoA
=A%AAAAA
��A
�A
�A
�`A
�A
��A
ȴA
��A
�!A
��A
��A
��A
��A
�A
bNA
ZA
I�A
E�A
I�A
I�A
I�A
I�A
I�A
M�A
M�A
M�A
M�A
I�A
I�A
M�A
I�A
M�A
I�A
M�A
M�A
M�A
I�A
M�A
M�A
M�A
M�A
M�A
M�A
M�A
I�A
M�A
M�A
E�A
E�A
M�A
I�A
I�A
E�A
A�A
A�A
=qA
A�A
A�A
A�A
A�A
A�A
=qA
=qA
5?A
1'A
 �A	�A	��A	�
A	�A	��A	�FA	A	�
A	�mA	�
A	�;A	�A	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	�A	XA	K�A	C�A	G�A	G�A	;dA	?}A	7LA	7LA	/A	�A	"�A	�A	VA��A��A�uAz�Az�AffA5?A��A�A�A|�AXAC�A?}A�A��A��A�uAjA5?AƨA�HA�9A�A��An�Av�Ar�AI�A=qA5?A-AbA  A��A�mA�TA�mA�TA�
AƨAƨAƨAA�wA�wAA��A`BA\)A\)G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  AZAVAQ�AQ�AVAM�AI�A5?A�
A�AoA��A�7A��A(�AA��A�hA�AdZAVA
�uA
I�A
I�A
1A	��A�HAQ�A1AdZAn�A�#A��A �@��-@�;d@���@��@�M�@���@��u@���@�\@�p�@�u@�=q@�  @�!@�=q@�^@�O�@�&�@�V@���@�D@�F@�\)@�v�@�@���@���@��#@�&�@�9@�j@�Q�@�9X@��@�t�@��@�5?@�-@�hs@��@�j@�j@�9X@ߕ�@�dZ@�33@��@ޏ\@�@�&�@�r�@��m@���@ۥ�@ە�@�C�@ڸR@ڟ�@���@��@��@���@�`B@��@ؓu@�Z@�9X@�b@��
@׶F@ׅ@ׅ@ׅ@�\)@�C�@�\)@ם�@��@أ�@��@��@��`@ؼj@�r�@�9X@��@׶F@�z�@���@��@ٙ�@�E�@ڗ�@�
=@�33@�S�@���@܃@ܓu@ܬ@ܬ@ܓu@�z�@�z�@܃@܋D@ܓu@ܛ�@ܓu@܋D@܋D@܃@�bN@�I�@�b@��@��@��;@�C�@���@ڟ�@��#@�x�@٩�@ٺ^@ٺ^@���@ٺ^@ٺ^@���@��#@��@�{@�J@��@���@ف@�%@ج@�l�@�~�@�@�`B@�~�@��
@أ�@�7L@���@�V@ڗ�@���@�o@��H@ڰ!@ڏ\@��@��T@�x�@��@�r�@׍P@�~�@պ^@�O�@�V@ԛ�@�r�@�bN@�A�@�(�@��;@�K�@ҏ\@�M�@�{@�/@�;d@��@Ώ\@�J@���@�@�@ͩ�@͉7@́@�p�@�X@��@��@���@�Ĝ@�Ĝ@�Ĝ@���@̴9@�Q�@��@�ƨ@˅@�;d@�"�@�@�-@�&�@��@��/@���@ȴ9@ȴ9@Ȭ@ȣ�@ȋD@�z�@�j@�Z@� �@��m@Ǯ@�\)@ƸR@�@�ƨ@ț�@��`@��`@ȼj@Ȭ@ȣ�@ȴ9@ȼj@���@�&�@���@ȣ�@ȓu@ȋD@�I�@�1@��
@ǍP@�S�@�33@�@��@Ƈ+@�V@�$�@��@���@�x�@��@�(�@�t�@�@�{@���@�(�@���@�S�@�$�@��@��^@��@�dZ@�^5@��@��#@���@�%@�A�@�"�@��y@�ȴ@�v�@�n�@�-@�@���@�J@��@�E�@��\@��!@��@�\)@�dZ@�\)@�C�@��@���@���@�n�@�J@��@��@���@�?}@���@�A�@�  @�ƨ@�l�@��\@�^5@�E�@�E�@�-@�5?@�-@�=q@�5?@�J@��T@���@��@��!@���@�@��-@���@�1'@�|�@��!@�J@��@�r�@�1'@��;@��F@��P@�S�@�33@�"�@���@���@���@�v�@�$�@��@�{@�@��T@���@���@���@��#@���@�@��-@���@�p�@���@��j@� �@���@��H@���@��\@�~�@�~�@�~�@�n�@�ff@�ff@�ff@�M�@�-@�@�p�@��@��`@���@�z�@�j@�I�@� �@��@�b@�1@�1@�1@���@��m@��;@��
@��w@���@�|�@�;d@�33@�"�@�o@�
=@�
=@�@��H@��R@�~�@�=q@�-@�-@�-@�$�@��@��@�{@�J@���@��#@��@�/@���@��9@��@��@��@���@��@��@��9@��9@��9@��j@�Ĝ@�Ĝ@��9@�9X@��@�|�@�;d@�ff@�@���@���@���@��@�x�@�x�@�hs@�p�@�hs@�`B@�`B@�`B@�X@�O�@�O�@�G�@�O�@�/@�&�@�&�@��@�%@�V@��9@��u@��@�I�@�1@���@���@��;@���@��w@��@��@��@���@���@�|�@�l�@�\)@�l�@�l�@�l�@�dZ@�dZ@�\)@�dZ@�dZ@�\)@�K�@�;d@�S�@�33@��y@���@�ȴ@�ȴ@�ȴ@�ȴ@�ȴ@��!@���@�v�@�V@�E�@�{@�@��@���@��^@��7@�`B@�%@��D@�Q�@�1'@��@��
@�+@�o@�o@���@��H@��R@���@��\@�n�@�^5@�^5@�V@�V@�J@���@���@��7@�`B@�X@�X@�X@�X@�`B@�`B@�`B@�`B@�`B@�`B@�`B@�`B@�hs@�hsAZAZAZAVAVAVAVAQ�AVAVAVAVA^5AVAVAQ�AZAVAM�AI�AM�AM�AM�AM�AM�AQ�AQ�AVAVAQ�AVAVAVAQ�AZAVAZAZAZAVAZAVAM�AI�AM�AI�AE�AI�AI�AM�AI�AI�AE�A5?AQ�AVA(�AE�AZAbNAQ�A(�A�A1'A(�AJA�A�AƨAƨA��A�wA��A�A�^A�
AAA�#AƨA|�Ap�A`BAXAO�AK�AG�AG�AC�A?}A7LA7LA33A/A+A"�AoAVA%AA��A��A�A�A�A�`A�HA�A�A��A��A��AȴA�jA�jA�jA�RA�RA�A�AbNAE�A1'A�AA�mA�#A�
A��AA�wA�-A��A�Ap�AXA7LA/A+A&�A&�A"�A"�A"�A�AoA
=AA��A�yA�/A��A��AĜA��Az�Ar�AjAVAI�AE�AA�A9XA5?A1'A-A(�A$�A �A�A �A �A�A�A�A�A�A{AbAbAbA1A1A1AA  A  A��A��A��A�A�A�A�A�A�mA�mA�#A�
A�
A��A��AƨAA�wA�wA�^A�FA�FA�FA�-A��A��A��A�hA�hA�PA�PA�PA�PA�PA�PA�PA�PA�hA�PA�PA�PA�PA�7A�7A�A�7A�A�A�A�A�A|�A|�A|�A�A�A�A�A|�A|�A|�At�At�At�Ap�At�Ap�AhsAhsA`BA`BA\)AXAO�AC�A;dA/A+A"�A�A�A�A�A�AoA
=A%AAAAA
��A
�A
�A
�`A
�A
��A
ȴA
��A
�!A
��A
��A
��A
��A
�A
bNA
ZA
I�A
E�A
I�A
I�A
I�A
I�A
I�A
M�A
M�A
M�A
M�A
I�A
I�A
M�A
I�A
M�A
I�A
M�A
M�A
M�A
I�A
M�A
M�A
M�A
M�A
M�A
M�A
M�A
I�A
M�A
M�A
E�A
E�A
M�A
I�A
I�A
E�A
A�A
A�A
=qA
A�A
A�A
A�A
A�A
A�A
=qA
=qA
5?A
1'A
 �A	�A	��A	�
A	�A	��A	�FA	A	�
A	�mA	�
A	�;A	�A	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	�A	XA	K�A	C�A	G�A	G�A	;dA	?}A	7LA	7LA	/A	�A	"�A	�A	VA��A��A�uAz�Az�AffA5?A��A�A�A|�AXAC�A?}A�A��A��A�uAjA5?AƨA�HA�9A�A��An�Av�Ar�AI�A=qA5?A-AbA  A��A�mA�TA�mA�TA�
AƨAƨAƨAA�wA�wAA��A`BA\)A\)G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�BÖBÖBÖBÖBÖBÖBÖBBĜBŢBƨBŢBĜBÖBĜBĜBŢBƨBƨBƨBǮB��B��B��B��B��B��B�B�ZB�`B�fB�fB�ZB�ZB�fB�ZB�`B�TB�NB�HB�HB�BB�BB�BB�;B�;B�;B�;B�;B�;B�;B�5B�5B�/B�/B�)B�)B�)B�)B�)B�)B�/B�5B�5B�/B�/B�)B�)B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBȴBȴBȴBǮBƨBŢBƨBĜBŢBŢBǮBƨBŢBĜBĜBĜBÖBÖBÖBŢBƨBƨBȴB��B��B��B��B��B�B�B�B�B�B�B�B�B�BB�TB�mB�B�B�B��B��B��B��BBBBBBB%B+B1B1B1B	7B1B1B1B1B	7B
=B
=B
=B	7B	7B+B+BB+B	7B
=B
=BDBDBDBJBPBVB\BbB\BVBPBDB	7BBB  BB
=BoB�B�B#�B%�B'�B+B-B.B-B-B,B+B)�B'�B%�B �B�B�B�B�B�B�B�B�B�B�B{BuBoBoBVB+BBBBBBBB%B%B+B+B1B	7B	7B	7B	7B	7B	7B	7B
=B
=B	7B	7B	7B	7B1B+B+B1B1B	7B	7B	7B
=B
=B
=B
=B
=B
=B
=B	7B	7B	7BDBbB�B�B"�B#�B$�B$�B$�B%�B%�B'�B,B-B.B0!B5?B6FB6FB7LB7LB7LB6FB6FB6FB6FB6FB6FB6FB5?B49B33B33B2-B0!B.B,B'�B&�B$�B�B�B�B�B�BuBoBhB\BPBDB1B+B+B%BBBB%B+B1B	7BDBJBPB\BhBhBhBoBoBoBoBuBuBuBoBhBhBhBhBbB\BVBPBPBPBPBPBPBPBJBJBDB	7B%BBBBB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B��B��B��B��B��B��B�B�B�B�B��B��B��B��B��B��B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�sB�sB�sB�sB�sB�sB�sB�sB�sB�mB�mB�mB�fB�fB�fB�fB�mB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�mBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBBÖBÖBÖBÖBÖBÖBÖBÖBBÖBBBÖBÖBÖBÖBÖBÖBÖBBBÖBÖBÖBBĜBBBÖBÖBBBBBBBBBŢBÖB��B��BBBŢBB��BÖBBBɺBƨBBƨBĜB��BÖBĜBŢBĜBBBȴBŢBƨBƨBƨBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBǮBƨBǮBƨBƨBƨBƨBƨBƨBƨBŢBƨBŢBŢBŢBŢBŢBŢBŢBĜBĜBĜBŢBƨBŢBŢBĜBŢBŢBŢBĜBĜBĜBĜBĜBÖBĜBĜBÖBŢBĜBĜBĜBÖBĜBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBÖBBBƨBÖBÖBĜBĜBĜBÖBĜBĜBÖBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBŢBĜBŢBĜBĜBĜBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBƨBƨBƨBƨBƨBǮBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBǮBƨBǮBǮBǮBǮBǮBǮBǮBǮBǮBȴBȴBȴBȴBǮBǮBǮBȴBȴBǮBȴB��BȴBɺBȴB��BɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺB��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B�B�
B�B�B�HB�TB�NB�NB�ZB�NB�HB�ZB�ZB�TB�TB�ZB�`B�ZB�`B�`B�ZB�ZB�ZB�`B�ZB�ZB�ZB�ZB�ZB�TB�NB�yB�`B�`B�ZG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  BñBçBØBÒB��B��B�fBėB�8B�B�B��B��B�~B�B�GB�EB��B�BǰB�'BˏB��BːB�3B�;B�TB�B�B��B�0B�LB�%B��B�9B��B�B��B�rB��B�DB�}B��B�B��B�B�HB��B� B��B�xB�]B�TB��B�jB��B݁B��B�8B�+B�gB�.B��BޢB�XB�RB�eB�B�ZB��B��B�nB�}B��B�B�OB��B�GB�7B�BлB�BͺB��BʑB��B��B��B�+B�tB��BǑBķBűB��B�VB�B�tB��B��B��B��B��B��BšBƨB��B��BʓB�bB�B�B�LB�B�qB�^BهB�mB؃B�SB��BߪB��B�B�vB�B�
B��B��B��B�*B�B�BB6B<BBBB B$B?B	=B4B?BbB\B	�B
kB
<B
[B
)B	�B�BFB�B�B	 B
8B
(BXBCB0B1B6BBhB�B�B�BB�BB�B6B �B�PB B2B�B�B#>B%{B'cB*�B-\B.^B-KB-�B,hB+�B*�B)B'LB"`B�BAB�B5B�B�B�B�B�BhB�B�B�B�BDB�B�B�BkB&BB<BMB-BBBRB�BnB	YB	LB	6B	6B	,B	bB	�B
�B
B	�B	�B	bB	tB	lB�B�BKBPB	ZB	7B	BB
HB
bB
XB
UB
ZB
�B
�B	�B	�B
'B
�B:BLB3B"�B$B$�B$�B$�B%�B%�B'�B,9B-�B.4B0:B5�B6�B6�B7�B7�B7}B6�B6�B6�B6�B6�B6�B6tB5�B4�B4�B4IB3lB1B.�B.kB(�B'�B&B "B)B�B_BB�B�B�BaB�B�B�BdB�B4B�BbB+BBB�B�B	B�B�BJBtB�B�B�B�B�BB�B�B�B?BbBB�B�B�B�B�BzBRBtBIB\B7B\B�B�B�B
cBNB*BmB2BdB  B�B�-B��B�VB��B�AB�IB�B�B�%B�B��B�
B�?B��B�B�CB��B��B��B��B��B��B��B��B��B��B��B��B�B��B�,B��B��B��B�0B��B��B��B��B��B��B��B��B��B��B�B��B�FB�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B�<B�4B�B�*B��B��B��B��B��B��B��B��B��B��B��B��B��B�xB�/B�cB�2B��B�LB� B��B�B��B�B�B��B�B�B�B�B�B�B�B�B�B�B��B�B�B�B�B�B�&B��B�B��B��B�B�B��B�B�B�B�B�B�B�B�B��B�B�B�B�B��B�B��B�B�B��B��B��B�rB��B�B��B�B�B�B�B�B�B�B�B�B�B��B�B�B��B�B��B��B�B�@B��B�B�B��B�B�B�}B�B�B�B�B�B�B�B�uB�B�{B��B�B�B�B�B�uB�fB�hB�kB�YB�gB�fB�gB�gB�fB�eB�fB�]B�kB�kBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBBÖBÖBÖBÖBÖBÖBÖBÖBBÖBBBÖBÖBÖBÖBÖBÖBÖBBBÖBÖBÖBBĜBBBÖBÖBBBBBBBBBŢBÖB��B��BBBŢBB��BÖBBBɺBƨBBƨBĜB��BÖBĜBŢBĜBBBȴBŢBƨBƨBƨBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBǮBƨBǮBƨBƨBƨBƨBƨBƨBƨBŢBƨBŢBŢBŢBŢBŢBŢBŢBĜBĜBĜBŢBƨBŢBŢBĜBŢBŢBŢBĜBĜBĜBĜBĜBÖBĜBĜBÖBŢBĜBĜBĜBÖBĜBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBÖBBBƨBÖBÖBĜBĜBĜBÖBĜBĜBÖBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBŢBĜBŢBĜBĜBĜBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBŢBƨBƨBƨBƨBƨBǮBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBǮBƨBǮBǮBǮBǮBǮBǮBǮBǮBǮBȴBȴBȴBȴBǮBǮBǮBȴBȴBǮBȴB��BȴBɺBȴB��BɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺB��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBɺBɺBɺB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B�B�
B�B�B�HB�TB�NB�NB�ZB�NB�HB�ZB�ZB�TB�TB�ZB�`B�ZB�`B�`B�ZB�ZB�ZB�`B�ZB�ZB�ZB�ZB�ZB�TB�NB�yB�`B�`B�ZG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                  <#�D<#��<#�<#�<#�J<#��<$Z�<'F<%��<%��<%k�<+��<'hA<&��<$v<$0.<$(<#��<#�m<$�3<%��<$XX<#�0<$Y�<%it<(c�<K;<B��<'��<(�D<&R`<$x+<)�1<9�|<.�<-Yc<$�<(�)<$�q<%��<$�e<'��<%�!<%rN<-c�<,sq<'�<$=<$MO<$"2<#�e<#��<#��<$@|<%�<$k<%>�<$5w<#׺<#�<#��<$�<$7�<#�H<#�+<#��<#�<$� <%�<$><<$W<#�N<$�<$��<#�{<#�<$x+<#�<#�4<#��<$j|<$��<$�<$�3<$c�<#ޫ<#�<#�o<$�<$U�<#�<$|d<#�D<#׺<#��<$-<#��<$]h<#��<#��<#�&<#�<#��<#�<#�<#�
<#�U<#��<#�<#�<$E<$XX<$1:<#�<#��<#��<#�N<#�<#�a<#�<$�w<$}<#�	<$R'<$��<$�<$+<#�<#�<$�Q<$Sa<#��<#�+<#�
<#��<#��<#�<#׎<#׎<#��<#׎<#ף<#�&<#�<#ף<#�^<#ܯ<#�<#�<#�<#��<$�b<$<<$�<$ʾ<$%<#�<#ا<#�<#�c<#�C<#�<#�C<#��<#�<#��<#�{<#�8<#ߜ<$ <$><<$"2<&�
<%�<$�q<$G<&y<'J�<%�<$�t<$�q<$k<#�	<$�<#��<#�<#�<#�e<$*<#�(<$&<$�<$��<%b�<%�j<%�<$,<#��<$3U<#�<#�D<#�<#ܯ<$ <$y�<$�(<#�g<#��<%b�<*i�<$<<$�<$H�<#�5<#�<#�<#�l<#�J<#�<<#ٛ<#ۮ<#��<#�e<#ڑ<#�c<#�<#�<#�i<#ܯ<$�<$�<#�U<#��<#�N<#ܯ<#�e<%�<%�<#��<#�<#��<#��<#�
<#�i<#�i<#�8<#�D<#��<#ٛ<#�m<#�<#�5<$�<$�t<#��<$�<%�<$<<#�<#�N<#�<#�{<#ا<#׺<#��<#��<#�^<$<#�*<#��<#�5<#�g<#�<#�H<#�<#�^<#�<#�l<$	<#�<#�<#�"<#�<$�<$(<%S�<$�;<%s<$y�<$#(<(;B<$)
<$Z�<%��<#��<#��<&�}<&D�<%�<#�N<#�g<#�<$�j<%�<&�<#�g<#��<#��<#׺<#�Q<#�<#�{<#��<#�o<#��<#�H<#�<#��<$%<#�<#�{<#�<$
�<#�4<#ܯ<#�N<$�<#��<#��<#�<$Z�<$�X<$"2<#�)<#�<$!><%<#�!<#�l<#�<#��<#�0<#�{<#��<#�{<#�<#�<%��<$��<'r#<$��<#�<#��<%K:<$�<$��<% <$�<%�!<$��<#��<$p<#��<#�<#�<#ޫ<#�]<#��<$�<#�<#�<$	<#׎<#׺<#��<#��<#�<#�i<#ף<#�i<#ף<#�<#��<#�o<#�<$E<#��<$<<$U�<$��<$ <#ף<#��<#�<#�
<#�c<#�i<#�<#�<#�8<#�J<#�l<$t <$<#�<#�<#��<#�D<#��<#�<#׎<#�{<#�{<#�<#�<#�C<#�<#�{<#׺<#�r<#�^<#��<#�)<#ף<#�<#ا<#�{<#�<#�X<#��<#�&<#��<#��<#��<#�<#�<#�{<#׎<#�<#׺<#׺<#��<#��<$�<$v<#�<#��<#ף<#�<#�<#�i<#�{<#�<#׎<#�<#�<#׎<#ף<#�<#��<$I�<#�a<$2G<$�<%"<$!><#�<#�4<#�&<#�<#ף<#�
<#��<#�i<#�{<#�{<#�<#�<#�{<#ף<#�<#�X<#׎<#ޫ<#�{<#�<#��<#��<#�<$
�<#�J<#�o<#��<#�<#�<#�
<#��<#��<#�<#�o<#�<#�<#��<#�D<#��<#��<#��<#؄<#�
<#�
<#׺<#�<#�$<#ף<#�<#��<#�<#�<#��<#��<#��<#ۮ<#ף<#�<#�<#�<#�<#��<#�l<#�^<#�J<#�*<#�<#�o<#�D<#�<#��<#�<#��<$r<$H�<#��<#ߜ<#�l<#��<$�<#��<#�<#�8<#�+<#�&<#�r<#�D<#�<#��<#�<#ף<#�<<$�<$,<#�<#ا<#��<#׺<#�
<#�<#�<#׎<#�<#�
<#�<#�<#�
<#�<#�
<#�I<#�<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                