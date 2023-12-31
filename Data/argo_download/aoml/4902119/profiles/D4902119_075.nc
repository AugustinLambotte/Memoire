CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  d   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92018-11-06T01:28:55Z creation; 2023-02-09T14:06:16Z DMQC;      
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
resolution        =���   axis      Z           <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  O�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���        T�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���        l|   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o        �   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o        ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o        �l   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o           PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  լ   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o        �t   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o        �\   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` |   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   0   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   8   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   @   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   H   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � P   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar           HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar           HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       $   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ,Argo profile    3.1 1.2 19500101000000  20181106012855  20230209090616  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               K   KAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�l�	{U�@�l�	{U�11  @�l�qڀ@�l�qڀ@N4�Ck��@N4�Ck���=/a���=/a��11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�=q@   @B�\@�G�@�G�@�G�@�  @��RA\)A#33A?\)A`  A\)A�\)A�  A�Q�A��AϮA�  A�Q�B   B�
B�
B  B   B(  B0  B8(�B@(�BG�BPQ�BW�B`  BhQ�BpQ�BxQ�B�(�B�  B�  B�(�B�=qB�(�B�{B�(�B�  B�{B�{B�  B�B��
B�B�B��B�  B��B�  B�{B�{B�=qB�{B��B�  B��B��B�  B�(�B�  B�(�C {C
=C  C��C�HC	�C�C��C
=C{C��C  C�C�C(�C{C �C!��C#�HC%�HC'��C*
=C+��C-��C0{C2  C4  C5��C8  C9��C;��C=�C?�HCB  CC��CE��CH{CJ
=CL  CM��CO��CR  CS��CV
=CX
=CZ
=C\
=C^  C`(�Cb
=Cd  Cf
=Cg��Cj  Cl
=Cm��Co�Cq�Cs��Cv
=Cx{Cz  C{�C~  C�C�C�  C���C���C�{C�  C���C�
=C�C���C���C���C�C�C��C�  C�C�  C���C��C�C�\C�
=C�  C�  C�C�  C�
=C���C���C�C�  C�C�
=C�
=C�
=C�
=C�C�
=C���C��C�  C�  C���C�\C�
=C�
=C�C���C���C�  C���C���C�  C�  C�  C�  C�  C���C���C��C�  C���C�C�C�  C�C�
=C�  C���C���C�C���C��C�  C�
=C���C���C��C��C�  C�
=C�
=C�  C���C���C���C�  C�  C�C�  C���C�  C�C���C���C�  C�
=C�C�C���C�  C�  C�C���C���C�
=C���C���C�C�C�  C���C���C���C���C�  C�C�  C�  C�  C�  C�C���C���C���C�  C���D ��D�D}qD�RDu�D�qD��D�D�D  D}qDD�=D�qDz�D�qD}qD�qD	xRD	�RD
}qD�D}qD�RDu�D�RDxRD  D� D�qD�DD�D
=D� D�RDz�D  D}qD��D}qD�D��D��D}qD�D}qD��Dz�D�qD��D  D��DD�D  Du�D��Dz�D��D� D�qD}qD �D � D!  D!}qD!�qD"}qD#  D#��D$  D$z�D%  D%� D&�D&�D'�D'�D(  D(}qD)  D)� D)�RD*z�D*�qD+� D,  D,��D-�D-� D.  D.��D/�D/�D0  D0z�D0��D1� D2D2��D3  D3}qD4  D4��D5D5u�D5�RD6� D7  D7� D7�qD8}qD8�qD9� D:�D:��D:�qD;z�D;�RD<xRD<��D=z�D=��D>}qD?  D?��D@�D@��D@�qDAxRDB�DB��DC�DC�DDDD�DE
=DE��DE��DF}qDF��DGxRDG�qDH}qDH��DIs3DI��DJz�DJ��DKxRDL  DL��DM  DM� DM��DN}qDO�DO� DP  DP�DQ�DQ� DQ�qDRz�DS  DS��DT�DT��DUDU�DVDV��DW  DWxRDW�qDX}qDY  DY��DZDZ� DZ��D[}qD[�qD\}qD]�D]�D^  D^z�D_  D_�D`�D`��D`�qDa� Db�Db��Db�qDc� Dd�Dd� Dd��Dez�De��DfxRDf��Dgz�Dg�qDh}qDi�Di� Di�qDj��Dk�Dk}qDl  Dl}qDl��Dm}qDm�qDn}qDn��Do}qDp  Dp� Dq�Dq}qDq�qDr� Dr�qDs� DtDt� Du  Du}qDu�qDv}qDv��Dwz�Dx  Dx}qDy  Dy��Dz  Dz}qD{  D{��D|�D|��D}  D}� D~  D~��D  Dz�D��D�=qD�}qD���D�HD�AHD�� D���D���D�@ D��HD�D��qD�=qD��HD��HD���D�>�D�� D��HD��D�@ D�� D��HD���D�@ D��HD��qD���D�AHD�� D�� D�  D�>�D�� D��HD���D�@ D�~�D��qD�  D�@ D�|)D��)D�HD�=qD�~�D��HD�  D�=qD�~�D�� D��qD�=qD�� D�� D��D�B�D��HD���D���D�AHD�� D�� D�  D�@ D��HD��HD�  D�>�D�� D��qD��qD�@ D�� D���D�  D�@ D��HD��HD�  D�AHD�� D��)D��qD�>�D�� D�� D��qD�=qD�}qD��HD�HD�AHD��D���?B�\?�  ?�  ?�\)?��
?�Q�?\?�
=?�G�?��@   @
=q@\)@
=@�R@#�
@(��@0��@5@J=q@E�@O\)@Q�@Y��@^�R@h��@u@p��@u@�  @��\@��@�ff@���@�\)@��@�Q�@�Q�@�
=@�p�@��\@��
@�ff@���@���@�\)@�33@�
=@���@�(�@��R@�G�@��
@�ff@���@���@�\)@��@�z�@�
=@ٙ�@��H@�p�@�G�@�\@��
@�=q@�=q@�\)@��@�33@�@���@�(�@��RA ��A�A�\A�
AAffAQ�A
�HA
=qA�A��A{A  A�
A33A�
AAffA��A�HA�RA�RA ��A!�A#�
A%A'�A*�HA+�A,��A0  A1�A3�
A5�A7
=A9��A;�A<��A?\)AAG�AC33AE�AH��AG�AK�AMp�AN{AQG�AS�
AUAXQ�AY��A\(�A^{A`  Ab�\Adz�AfffAhQ�Aj�HAl��Al��Ap��Ar�\Au�AvffAy��A{�A|(�A~{A�Q�A�G�A��A��HA�(�A��A���A��RA�  A�G�A�G�A��HA��A���A�{A��RA�A���A�G�A��\A��A�z�A��A�{A�\)A�Q�A�G�A��A�33A��HA���A�p�A��RA�Q�A���A���A��\A��A�(�A���A�ffA�
=A�  A���A�=qA��HA��
A��A�{A��RA��A���A��A��\A��A���A�A�ffA�\)A���A���A�=qA�33A���A��A�{A�\)A�Q�A���A��A�33A�(�A��AǮA�
=A�Q�A���A��A��HA��
A���A�AθRA�  AУ�Aљ�A��HAҏ\A�z�A�p�AָRAָRA�Q�A�G�A�=qAۅA�(�A��A�ffA�\)A�  A���A�=qA�33A��
A���A�{A�RA�A��A��A�\A�A�z�A�A�ffA�
=A�Q�A�G�A�=qA��HA�(�A��A�{A��RA�  A�G�A��A��\A��
A���A�A�ffA�\)B Q�B ��BG�B��B=qB�\B
=B\)Bz�Bz�B��B�B�B=qB�RB\)B�
B(�B��B	G�B	B
{B
ffB33B�B  Bz�B�Bp�B��B�\B�HB33B  BQ�B��BG�B�BffB�RB33B�BQ�B��B�B��B=qB�\B
=B�B�
Bz�B��Bp�BBffB�HB\)B�
BQ�B��BG�BB{B
=B
=B�B   B Q�B ��B!p�B!�B"=qB"ffB#\)B#�
B$z�B$��B%G�B&=qB%�B&�RB'
=B'�B(  B(z�B(��B)p�B)�B*=qB*=qB+\)B+�
B,��B,��B-��B-G�B-��B/
=B/
=B/
=B0(�B0��B0��B1p�B2{B2�\B2�HB3�B4  B4z�B4��B5B5�B6ffB6�HB7\)B7�
B8Q�B8��B9p�B9B:ffB:�HB;\)B;�
B<Q�B<��B=G�B=�B>=qB>�RB?33B?33B@  B@��BA�BA��BB{BC33BC
=BC�BD  BDz�BD��BEp�BF{BFffBF�HBG\)G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       ?�=q@   @B�\@�G�@�G�@�G�@�  @��RA\)A#33A?\)A`  A\)A�\)A�  A�Q�A��AϮA�  A�Q�B   B�
B�
B  B   B(  B0  B8(�B@(�BG�BPQ�BW�B`  BhQ�BpQ�BxQ�B�(�B�  B�  B�(�B�=qB�(�B�{B�(�B�  B�{B�{B�  B�B��
B�B�B��B�  B��B�  B�{B�{B�=qB�{B��B�  B��B��B�  B�(�B�  B�(�C {C
=C  C��C�HC	�C�C��C
=C{C��C  C�C�C(�C{C �C!��C#�HC%�HC'��C*
=C+��C-��C0{C2  C4  C5��C8  C9��C;��C=�C?�HCB  CC��CE��CH{CJ
=CL  CM��CO��CR  CS��CV
=CX
=CZ
=C\
=C^  C`(�Cb
=Cd  Cf
=Cg��Cj  Cl
=Cm��Co�Cq�Cs��Cv
=Cx{Cz  C{�C~  C�C�C�  C���C���C�{C�  C���C�
=C�C���C���C���C�C�C��C�  C�C�  C���C��C�C�\C�
=C�  C�  C�C�  C�
=C���C���C�C�  C�C�
=C�
=C�
=C�
=C�C�
=C���C��C�  C�  C���C�\C�
=C�
=C�C���C���C�  C���C���C�  C�  C�  C�  C�  C���C���C��C�  C���C�C�C�  C�C�
=C�  C���C���C�C���C��C�  C�
=C���C���C��C��C�  C�
=C�
=C�  C���C���C���C�  C�  C�C�  C���C�  C�C���C���C�  C�
=C�C�C���C�  C�  C�C���C���C�
=C���C���C�C�C�  C���C���C���C���C�  C�C�  C�  C�  C�  C�C���C���C���C�  C���D ��D�D}qD�RDu�D�qD��D�D�D  D}qDD�=D�qDz�D�qD}qD�qD	xRD	�RD
}qD�D}qD�RDu�D�RDxRD  D� D�qD�DD�D
=D� D�RDz�D  D}qD��D}qD�D��D��D}qD�D}qD��Dz�D�qD��D  D��DD�D  Du�D��Dz�D��D� D�qD}qD �D � D!  D!}qD!�qD"}qD#  D#��D$  D$z�D%  D%� D&�D&�D'�D'�D(  D(}qD)  D)� D)�RD*z�D*�qD+� D,  D,��D-�D-� D.  D.��D/�D/�D0  D0z�D0��D1� D2D2��D3  D3}qD4  D4��D5D5u�D5�RD6� D7  D7� D7�qD8}qD8�qD9� D:�D:��D:�qD;z�D;�RD<xRD<��D=z�D=��D>}qD?  D?��D@�D@��D@�qDAxRDB�DB��DC�DC�DDDD�DE
=DE��DE��DF}qDF��DGxRDG�qDH}qDH��DIs3DI��DJz�DJ��DKxRDL  DL��DM  DM� DM��DN}qDO�DO� DP  DP�DQ�DQ� DQ�qDRz�DS  DS��DT�DT��DUDU�DVDV��DW  DWxRDW�qDX}qDY  DY��DZDZ� DZ��D[}qD[�qD\}qD]�D]�D^  D^z�D_  D_�D`�D`��D`�qDa� Db�Db��Db�qDc� Dd�Dd� Dd��Dez�De��DfxRDf��Dgz�Dg�qDh}qDi�Di� Di�qDj��Dk�Dk}qDl  Dl}qDl��Dm}qDm�qDn}qDn��Do}qDp  Dp� Dq�Dq}qDq�qDr� Dr�qDs� DtDt� Du  Du}qDu�qDv}qDv��Dwz�Dx  Dx}qDy  Dy��Dz  Dz}qD{  D{��D|�D|��D}  D}� D~  D~��D  Dz�D��D�=qD�}qD���D�HD�AHD�� D���D���D�@ D��HD�D��qD�=qD��HD��HD���D�>�D�� D��HD��D�@ D�� D��HD���D�@ D��HD��qD���D�AHD�� D�� D�  D�>�D�� D��HD���D�@ D�~�D��qD�  D�@ D�|)D��)D�HD�=qD�~�D��HD�  D�=qD�~�D�� D��qD�=qD�� D�� D��D�B�D��HD���D���D�AHD�� D�� D�  D�@ D��HD��HD�  D�>�D�� D��qD��qD�@ D�� D���D�  D�@ D��HD��HD�  D�AHD�� D��)D��qD�>�D�� D�� D��qD�=qD�}qD��HD�HD�AHD��D���?B�\?�  ?�  ?�\)?��
?�Q�?\?�
=?�G�?��@   @
=q@\)@
=@�R@#�
@(��@0��@5@J=q@E�@O\)@Q�@Y��@^�R@h��@u@p��@u@�  @��\@��@�ff@���@�\)@��@�Q�@�Q�@�
=@�p�@��\@��
@�ff@���@���@�\)@�33@�
=@���@�(�@��R@�G�@��
@�ff@���@���@�\)@��@�z�@�
=@ٙ�@��H@�p�@�G�@�\@��
@�=q@�=q@�\)@��@�33@�@���@�(�@��RA ��A�A�\A�
AAffAQ�A
�HA
=qA�A��A{A  A�
A33A�
AAffA��A�HA�RA�RA ��A!�A#�
A%A'�A*�HA+�A,��A0  A1�A3�
A5�A7
=A9��A;�A<��A?\)AAG�AC33AE�AH��AG�AK�AMp�AN{AQG�AS�
AUAXQ�AY��A\(�A^{A`  Ab�\Adz�AfffAhQ�Aj�HAl��Al��Ap��Ar�\Au�AvffAy��A{�A|(�A~{A�Q�A�G�A��A��HA�(�A��A���A��RA�  A�G�A�G�A��HA��A���A�{A��RA�A���A�G�A��\A��A�z�A��A�{A�\)A�Q�A�G�A��A�33A��HA���A�p�A��RA�Q�A���A���A��\A��A�(�A���A�ffA�
=A�  A���A�=qA��HA��
A��A�{A��RA��A���A��A��\A��A���A�A�ffA�\)A���A���A�=qA�33A���A��A�{A�\)A�Q�A���A��A�33A�(�A��AǮA�
=A�Q�A���A��A��HA��
A���A�AθRA�  AУ�Aљ�A��HAҏ\A�z�A�p�AָRAָRA�Q�A�G�A�=qAۅA�(�A��A�ffA�\)A�  A���A�=qA�33A��
A���A�{A�RA�A��A��A�\A�A�z�A�A�ffA�
=A�Q�A�G�A�=qA��HA�(�A��A�{A��RA�  A�G�A��A��\A��
A���A�A�ffA�\)B Q�B ��BG�B��B=qB�\B
=B\)Bz�Bz�B��B�B�B=qB�RB\)B�
B(�B��B	G�B	B
{B
ffB33B�B  Bz�B�Bp�B��B�\B�HB33B  BQ�B��BG�B�BffB�RB33B�BQ�B��B�B��B=qB�\B
=B�B�
Bz�B��Bp�BBffB�HB\)B�
BQ�B��BG�BB{B
=B
=B�B   B Q�B ��B!p�B!�B"=qB"ffB#\)B#�
B$z�B$��B%G�B&=qB%�B&�RB'
=B'�B(  B(z�B(��B)p�B)�B*=qB*=qB+\)B+�
B,��B,��B-��B-G�B-��B/
=B/
=B/
=B0(�B0��B0��B1p�B2{B2�\B2�HB3�B4  B4z�B4��B5B5�B6ffB6�HB7\)B7�
B8Q�B8��B9p�B9B:ffB:�HB;\)B;�
B<Q�B<��B=G�B=�B>=qB>�RB?33B?33B@  B@��BA�BA��BB{BC33BC
=BC�BD  BDz�BD��BEp�BF{BFffBF�HBG\)G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�Av�An�A�A�DA��A�\A�+A�\A�DAM�A��A��A��A��A�uAz�A5?A1A�-A��A��A��A��A�\A�AQ�A7L@��H@��T@�?}@���@��@�A�@��@�
=@�9@�r�@웦@�u@�D@�V@�p�@�x�@���@�V@�V@��`@�Ĝ@�j@�9@�u@�@�z�@�I�@�1@�
=@�V@�E�@�+@�C�@ꗍ@�9@��@�ƨ@�P@���@��H@�K�@��m@�(�@�1'@� �@��@�F@�33@���@��@���@柾@�n�@�V@�@�h@�7@�p�@�X@�G�@�V@���@�j@䛦@��@�+@��@���@�%@��@߮@߅@�dZ@�S�@�dZ@�|�@�|�@��
@�@�ff@�!@���@�@�;d@��@�ȴ@���@�h@�7@�hs@��@�=q@�E�@�&�@��/@�/@�-@���@�@�^@�h@�`B@�V@��@�Z@�j@�z�@��@�9@�9@�@��D@�A�@߾w@�t�@�K�@�S�@�+@�@�
=@���@���@���@��@���@��@���@ݩ�@݁@�`B@�X@�`B@�p�@�hs@�`B@ݑh@�J@ް!@��y@��y@���@ާ�@ް!@�^5@�-@�5?@�@�Ĝ@�j@�r�@�b@�9X@���@��@�b@�9X@�Q�@�A�@�9X@�9X@�1'@�(�@�b@� �@�Z@ܓu@���@ܣ�@�bN@ܛ�@� �@�b@�(�@�Q�@�j@�j@�Z@�I�@�A�@�A�@�A�@�A�@�A�@�9X@��@�dZ@�+@��@ڏ\@�=q@�@���@��@��#@���@ٺ^@ى7@�G�@�V@أ�@�Q�@��
@ו�@��@��y@֏\@�^5@�$�@��@�@��@��#@ՙ�@�hs@�G�@�V@�%@�V@�?}@���@�l�@�$�@��T@�@�-@�V@�~�@�x�@�7L@��`@���@�O�@�p�@�{@���@� �@Լj@�V@�?}@��`@��@�Ĝ@ԓu@�9X@��@���@��
@��;@�  @ӶF@�o@�M�@��@�V@���@��@�J@�J@Ѳ-@Ь@�1@��@��@��@θR@�ff@͡�@�Z@�I�@�I�@�Z@�I�@�+@�V@��@���@Ǯ@�|�@ǅ@�~�@�V@��@Ĭ@ě�@ă@�bN@å�@\@���@��-@��7@�`B@�p�@�p�@���@�&�@�^5@�=q@�-@�V@�E�@��@�p�@�&�@�%@�%@���@���@��@�1'@�1@�l�@��H@�~�@��@��@��-@�?}@���@���@�Ĝ@��9@��@��@���@�z�@�9X@�t�@��@�X@�%@���@�bN@�Z@�Z@�Q�@�b@��@���@��P@�C�@��+@�^5@�J@���@�@�&�@��`@��@��
@�S�@��\@�$�@�p�@�bN@��@��R@���@�p�@�V@�Ĝ@��9@�bN@��w@�t�@��@��@�?}@��@��`@���@���@�"�@�ȴ@�V@���@��@�7L@�/@�O�@�G�@�O�@�?}@�/@�&�@��@�&�@�7L@�/@��j@�  @���@�l�@�C�@���@���@�V@���@�p�@�/@�V@��@�  @�@�v�@�=q@��@�7L@��D@�1'@�(�@��;@�\)@�
=@���@��\@�n�@�V@�M�@�=q@��@�{@�@���@�@��h@�O�@�/@�&�@�&�@�&�@��@��@��@���@��@��D@�r�@�r�@�j@�Z@�Q�@� �@�ƨ@��F@���@���@�l�@�K�@�;d@��@��@��H@���@��!@��!@���@���@�n�@�^5@�M�@�$�@�J@��@��@��@��@��@��T@�@��7@�x�@���@���@�X@�O�@�G�@��/@�Ĝ@�Ĝ@��j@��j@��j@�Ĝ@�Ĝ@��j@��j@�Ĝ@��j@�Ĝ@��j@��j@��j@��j@��j@�Ĝ@��j@�z�@�Q�@�I�@�A�@�9X@��@�1@���@��;@��;@��;@��
@��
@���@�ƨ@���@���@��P@��P@��P@��@��@�|�@�|�@�|�@�t�@�l�@�K�@�+@�
=@��y@�ȴ@���@���@���@���@���@��\@��+@��\@��\@��\@��\@��\@��\@��\@��\@���@��\@��\@��\@��+@�~�@�n�@�V@�M�@�=q@�-@�{@�@���@��#@�@�@��^@��-@���@���@���@��h@��@�O�@��@��@��/@�Ĝ@���@���@���@��u@��@�z�@�j@�j@�Q�@� �@���@�\)An�An�Ar�Av�A~�A~�A~�Az�Az�An�An�AjA^5AffAffAr�Az�Ar�Av�A~�A�DA�uA�uA�\A�\A�DA�DA�DA�\A�DA�DA�DA�DA�DA�\A�DA�DA�DA�uA��A��A��A��A��A��A��A��A�uA�DA�DA�\A�DA�\A�\A�DA�uA�uA�\A�\A�DA�+A�+A�+A�A�A�A�A�+A�DA�\A�\A�\A�\A�\A�\A�uA��A�uA�\A�\A�DA�DA�DA�\A�uA�uA��A�uA�uA�uA�+AjAz�A�+A�+A�DA�DA�DAz�AQ�AbAAAA  A��A��A  A�A�A�A��AA  AJAJAA��A�A�A�;A�mA5?A^5AA�A5?AJA��A�^A��A�FA��A=qA1'A�A��A�TA�;A�A��A�^A�A��A��A��A��A��A��A��A|�AdZAXA\)AO�A`BAG�A+A&�A�A�A�AVAA��A�HA�jA�A��A��A��A��A��A��A��A��A�uA�\A�\A�\A�\A�DA�DA�DA�DA�DA�DA�DA�+A�+A�+A�+A�+A�A�Az�Av�Az�Av�Az�Av�Ar�AjAbNAbNAVAM�AM�AE�AA�A9XA5?A(�A(�A(�A(�A-A(�A(�A(�A$�A(�A(�A$�A�A�A{AJA1AAA  A  A  A��A�A�A�A�A�mA�TA�;A��AA�wA�^A�^A�^A�-A��A��A��A�PAt�A;dA33A33A"�A
=A��A��A��A��A�A�A�A�A�yA�yA�`A�HA�HA�HA�/A�A��A��A��A��A��A��AȴA��AȴAĜAĜA��A��A��A��A��A��A��A��A��A��AĜAĜAĜAĜA��A��A��A�jA�9A�!A�A��A��A��A��A��A��A��A��A��A��A��A��A��A��A�uA�uA�uA�uA�uA�uA�\A�\A�\A�\A�uA�\A�DA�DA�+A�+A�DA�DA�+A�+A�+A�+A�A�+A�+A�+A�DA�DA�+A�A�A~�Az�Az�Az�Ar�An�AffAbNAZAVAM�AM�AQ�AM�AM�AM�AM�AE�AA�A5?A-A$�A �AbAA�A�^A��At�AG�A�A �A ȴA �!A n�A 5?A $�A J@�ƨ@���@�l�@�S�@�C�@�"�@���@��@��H@���@���@�ff@�E�@�5?@�5?@�5?@�=q@�=q@�5?@�5?@�$�@�@���@���@���@���@���@���@��^@���@���@��7@��@�x�@�x�@�p�@�p�@�hs@�hs@�XG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       Av�An�A�A�DA��A�\A�+A�\A�DAM�A��A��A��A��A�uAz�A5?A1A�-A��A��A��A��A�\A�AQ�A7L@��H@��T@�?}@���@��@�A�@��@�
=@�9@�r�@웦@�u@�D@�V@�p�@�x�@���@�V@�V@��`@�Ĝ@�j@�9@�u@�@�z�@�I�@�1@�
=@�V@�E�@�+@�C�@ꗍ@�9@��@�ƨ@�P@���@��H@�K�@��m@�(�@�1'@� �@��@�F@�33@���@��@���@柾@�n�@�V@�@�h@�7@�p�@�X@�G�@�V@���@�j@䛦@��@�+@��@���@�%@��@߮@߅@�dZ@�S�@�dZ@�|�@�|�@��
@�@�ff@�!@���@�@�;d@��@�ȴ@���@�h@�7@�hs@��@�=q@�E�@�&�@��/@�/@�-@���@�@�^@�h@�`B@�V@��@�Z@�j@�z�@��@�9@�9@�@��D@�A�@߾w@�t�@�K�@�S�@�+@�@�
=@���@���@���@��@���@��@���@ݩ�@݁@�`B@�X@�`B@�p�@�hs@�`B@ݑh@�J@ް!@��y@��y@���@ާ�@ް!@�^5@�-@�5?@�@�Ĝ@�j@�r�@�b@�9X@���@��@�b@�9X@�Q�@�A�@�9X@�9X@�1'@�(�@�b@� �@�Z@ܓu@���@ܣ�@�bN@ܛ�@� �@�b@�(�@�Q�@�j@�j@�Z@�I�@�A�@�A�@�A�@�A�@�A�@�9X@��@�dZ@�+@��@ڏ\@�=q@�@���@��@��#@���@ٺ^@ى7@�G�@�V@أ�@�Q�@��
@ו�@��@��y@֏\@�^5@�$�@��@�@��@��#@ՙ�@�hs@�G�@�V@�%@�V@�?}@���@�l�@�$�@��T@�@�-@�V@�~�@�x�@�7L@��`@���@�O�@�p�@�{@���@� �@Լj@�V@�?}@��`@��@�Ĝ@ԓu@�9X@��@���@��
@��;@�  @ӶF@�o@�M�@��@�V@���@��@�J@�J@Ѳ-@Ь@�1@��@��@��@θR@�ff@͡�@�Z@�I�@�I�@�Z@�I�@�+@�V@��@���@Ǯ@�|�@ǅ@�~�@�V@��@Ĭ@ě�@ă@�bN@å�@\@���@��-@��7@�`B@�p�@�p�@���@�&�@�^5@�=q@�-@�V@�E�@��@�p�@�&�@�%@�%@���@���@��@�1'@�1@�l�@��H@�~�@��@��@��-@�?}@���@���@�Ĝ@��9@��@��@���@�z�@�9X@�t�@��@�X@�%@���@�bN@�Z@�Z@�Q�@�b@��@���@��P@�C�@��+@�^5@�J@���@�@�&�@��`@��@��
@�S�@��\@�$�@�p�@�bN@��@��R@���@�p�@�V@�Ĝ@��9@�bN@��w@�t�@��@��@�?}@��@��`@���@���@�"�@�ȴ@�V@���@��@�7L@�/@�O�@�G�@�O�@�?}@�/@�&�@��@�&�@�7L@�/@��j@�  @���@�l�@�C�@���@���@�V@���@�p�@�/@�V@��@�  @�@�v�@�=q@��@�7L@��D@�1'@�(�@��;@�\)@�
=@���@��\@�n�@�V@�M�@�=q@��@�{@�@���@�@��h@�O�@�/@�&�@�&�@�&�@��@��@��@���@��@��D@�r�@�r�@�j@�Z@�Q�@� �@�ƨ@��F@���@���@�l�@�K�@�;d@��@��@��H@���@��!@��!@���@���@�n�@�^5@�M�@�$�@�J@��@��@��@��@��@��T@�@��7@�x�@���@���@�X@�O�@�G�@��/@�Ĝ@�Ĝ@��j@��j@��j@�Ĝ@�Ĝ@��j@��j@�Ĝ@��j@�Ĝ@��j@��j@��j@��j@��j@�Ĝ@��j@�z�@�Q�@�I�@�A�@�9X@��@�1@���@��;@��;@��;@��
@��
@���@�ƨ@���@���@��P@��P@��P@��@��@�|�@�|�@�|�@�t�@�l�@�K�@�+@�
=@��y@�ȴ@���@���@���@���@���@��\@��+@��\@��\@��\@��\@��\@��\@��\@��\@���@��\@��\@��\@��+@�~�@�n�@�V@�M�@�=q@�-@�{@�@���@��#@�@�@��^@��-@���@���@���@��h@��@�O�@��@��@��/@�Ĝ@���@���@���@��u@��@�z�@�j@�j@�Q�@� �@���@�\)An�An�Ar�Av�A~�A~�A~�Az�Az�An�An�AjA^5AffAffAr�Az�Ar�Av�A~�A�DA�uA�uA�\A�\A�DA�DA�DA�\A�DA�DA�DA�DA�DA�\A�DA�DA�DA�uA��A��A��A��A��A��A��A��A�uA�DA�DA�\A�DA�\A�\A�DA�uA�uA�\A�\A�DA�+A�+A�+A�A�A�A�A�+A�DA�\A�\A�\A�\A�\A�\A�uA��A�uA�\A�\A�DA�DA�DA�\A�uA�uA��A�uA�uA�uA�+AjAz�A�+A�+A�DA�DA�DAz�AQ�AbAAAA  A��A��A  A�A�A�A��AA  AJAJAA��A�A�A�;A�mA5?A^5AA�A5?AJA��A�^A��A�FA��A=qA1'A�A��A�TA�;A�A��A�^A�A��A��A��A��A��A��A��A|�AdZAXA\)AO�A`BAG�A+A&�A�A�A�AVAA��A�HA�jA�A��A��A��A��A��A��A��A��A�uA�\A�\A�\A�\A�DA�DA�DA�DA�DA�DA�DA�+A�+A�+A�+A�+A�A�Az�Av�Az�Av�Az�Av�Ar�AjAbNAbNAVAM�AM�AE�AA�A9XA5?A(�A(�A(�A(�A-A(�A(�A(�A$�A(�A(�A$�A�A�A{AJA1AAA  A  A  A��A�A�A�A�A�mA�TA�;A��AA�wA�^A�^A�^A�-A��A��A��A�PAt�A;dA33A33A"�A
=A��A��A��A��A�A�A�A�A�yA�yA�`A�HA�HA�HA�/A�A��A��A��A��A��A��AȴA��AȴAĜAĜA��A��A��A��A��A��A��A��A��A��AĜAĜAĜAĜA��A��A��A�jA�9A�!A�A��A��A��A��A��A��A��A��A��A��A��A��A��A��A�uA�uA�uA�uA�uA�uA�\A�\A�\A�\A�uA�\A�DA�DA�+A�+A�DA�DA�+A�+A�+A�+A�A�+A�+A�+A�DA�DA�+A�A�A~�Az�Az�Az�Ar�An�AffAbNAZAVAM�AM�AQ�AM�AM�AM�AM�AE�AA�A5?A-A$�A �AbAA�A�^A��At�AG�A�A �A ȴA �!A n�A 5?A $�A J@�ƨ@���@�l�@�S�@�C�@�"�@���@��@��H@���@���@�ff@�E�@�5?@�5?@�5?@�=q@�=q@�5?@�5?@�$�@�@���@���@���@���@���@���@��^@���@���@��7@��@�x�@�x�@�p�@�p�@�hs@�hs@�XG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�ZB�ZB�NB�NB�NB�NB�HB�NB�NB�HB�HB�BB�BB�;B�;B�;B�HB�yB��B��B��B�B�B��B��B%BJBPBbB�B�B�B#�B%�B&�B'�B)�B)�B+B)�B)�B)�B)�B)�B(�B'�B$�B!�B!�B&�B'�B%�B�B�B�B �B�B"�B%�B(�B,B,B,B,B+B)�B(�B'�B'�B&�B%�B%�B#�B"�B"�B"�B"�B!�B �B�B�B�B�B�B�BhBPB1B+B+B+B+B1B
=BDBVB�B�B�B�B �B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B{B{B�B�B�B�B�B�B�B�B�B �B �B!�B!�B �B �B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B �B!�B"�B"�B#�B#�B#�B"�B#�B#�B#�B#�B#�B#�B"�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B{BuBoBoBhBhBhBbBbB\B\BVBVBVBVBhBbBDB+B+B1B	7B
=BDB1B+B+B1BJBPBhB�B�B"�B%�B'�B(�B+B+B+B,B+B+B+B+B,B,B+B(�B(�B+B+B-B-B-B.B-B+B+B)�B)�B(�B(�B&�B%�B&�B&�B%�B$�B"�B�B�B�B�B{BuBhB\BJBDB
=B
=B	7B	7B1B1B+B1B
=BDBPBVBoB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuBuB{B{B{B{B{B{B{B{BuBuB{BuBuBuBuBoBoBhBhBbBbBbBbBVBPBDB
=B
=B
=B
=B
=B	7B	7B1B1B1B1B1B1B1B+B%BBBBBBBBBBBBBBBBBBBBB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�yB�sB�yB�yB�yB�yB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�sB�mB�mB�mB�mB�mB�mB�mB�mB�mB�mB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�fB�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�`B�`B�`B�fB�`B�fB�`B�`B�fB�ZB�ZB�ZB�ZB�ZB�fB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�ZB�`B�`B�`B�fB�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�ZB�ZB�ZB�`B�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�fB�ZB�TB�TB�TB�TB�TB�NB�ZB�sB�mB�ZB�ZB�ZB�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�`B�TB�ZB�TB�ZB�ZB�ZB�ZB�ZB�NB�/B�TB�NB�TB�`B�HB�ZB�ZB�TB�;B�HB�ZB�fB�NB�NB�HB�NB�ZB�`B�TB�HB�TB�NB�TB�NB�HB�NB�TB�NB�HB�HB�HB�NB�TB�NB�NB�HB�NB�NB�HB�NB�NB�`B�NB�NB�NB�NB�NB�NB�HB�NB�HB�NB�NB�NB�NB�NB�HB�NB�NB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�HB�HB�HB�HB�HB�HB�HB�NB�NB�NB�HB�NB�HB�NB�NB�NB�HB�NB�NB�HB�HB�HB�NB�HB�HB�NB�NB�NB�NB�NB�NB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�HB�HB�HB�NB�HB�HB�BB�HB�HB�HB�BB�BB�BB�BB�HB�HB�BB�BB�HB�HB�HB�HB�HB�HB�HB�BB�HB�HB�HB�BB�HB�BB�BB�BB�HB�BB�HB�BB�BB�BB�BB�BB�HB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�;B�BB�BB�BB�;B�;B�;B�5B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�5B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�BB�;B�;B�;B�;B�;B�;B�;B�;B�;B�BB�BB�BB�BB�BB�BB�BB�BB�HB�HB�HB�NB�HB�HB�HB�HB�HB�HB�HB�HB�HB�HB�HB�NB�NB�ZB�fB�`B�fB�yB�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       B�hB��B�,B�0B�B�B�NB�B�B�lB�[B�mB�6B�|B�B� B��B�TB�VB��B�rB��B�pB�iB��B�B�B�<B��B��B�B�sB��B�B	�B�BBnB�B�B"B#�B&�B&�B'�B*9B*0B+B*B*-B*B*B*FB)lB)XB%�B!�B �B&�B)B(�B!aB�BB!�B�B".B$�B(�B+�B,B,B,�B+�B*RB)#B( B(>B'4B&B&�B$'B"�B"�B"�B"�B"#B!$B�B�BBOBgB�B�B�B�BjBZBCBBB
.B
�B�B�B8B�B`B nB"B!FB�BPB�B�B5B�B�B^B
BB�B�B�B�B�BB<B�B�B�B�B�BhB�B�B�B$BfBB�B�B�B�B�B�B�B�BB�B�B�B�B�B�B�BtBlB�B�B4B�B�BXB�B!B �B!�B">B!B �B!zBB.B�BBWB�B�BXBNBtB�B�B�B�B�B�B�BEBWBWB  BBzB!vB �B!�B"�B"�B#�B#�B#�B"�B#�B#�B#�B#�B#�B$B#�B" B!B RB7B
B�B�B�B�B�B�BB�BBBBEB�B.B�B�B�B�BvB�B�B�B�B�B�B�BcBGB#B+BnB8B�B�B�B�B
B�B�B�BGBrBBLBB�B�B"JB%�B(iB(�B+~B+VB+�B,qB+ B+1B*�B*�B,xB- B,$B)>B(�B+�B*�B,�B-B-�B/�B.B+8B,?B*gB*0B)~B*+B(�B&B&�B&�B&B&�B$B %BBB�B�B�B�BFB�B`B
cB
yB
YB
�B	^B^BqBmB
%BQB�B�B�B�B�BjB�BCBiBB�B�B�B�B_B)B�B�BzB@B,B�B�BEB�B�B�B�B�B�B�B�BB�B�B|B�BB�B�B}B�B�BB�B�BB�B�B�B�B�B\B�BBoB<B�BB�B B�B�B�B
�B
�B
�B
]B
�B
;B	�B�B	�B	KB�BDB�B	lB:B�B�B�B�B�B$B�B#B	B+B+B"B"B
B�B(B�B*B�BdBFB pB �B�{B��B��B�WB�,B��B��B�`B��B�6B�SB�B��B�LB��B�:B��B�AB�4B�	B��B��B��B��B��B��B��B��B��B��B�B��B��B�B�B�B��B�B��B��B��B��B�B�B��B�B��B�&B��B��B�B��B��B�B��B��B�B��B�B�B�B�B��B�B��B��B��B��B��B�B�B�B��B��B��B�B�YB�B��B�B�B�*B�B�B�B�B�B�zB�B�B�B�{B�B�zB�B�B�B�B�B�{B�B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�~B�B�B�B�B�B�~B�B�B�B�B�B�B�B�B�|B�yB�B�B�B�B�nB�xB�yB�xB�yB�xB�zB�uB�gB�~B�zB�yB�B�B�B�B�B�B�B�B�B�B�B�B�vB�~B�}B�|B�{B�|B�zB�B�B��B�B�B�B�B�rB�fB�tB�~B�uB�B�jB�B�B�.B��B�{B�`B�`B�`B�ZB�`B�`B�`B�`B�`B�`B�`B�fB�`B�fB�`B�`B�fB�ZB�ZB�ZB�ZB�ZB�fB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�ZB�`B�`B�`B�fB�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�`B�ZB�`B�`B�`B�`B�ZB�ZB�ZB�`B�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�ZB�fB�ZB�TB�TB�TB�TB�TB�NB�ZB�sB�mB�ZB�ZB�ZB�`B�ZB�ZB�ZB�ZB�ZB�ZB�ZB�`B�TB�ZB�TB�ZB�ZB�ZB�ZB�ZB�NB�/B�TB�NB�TB�`B�HB�ZB�ZB�TB�;B�HB�ZB�fB�NB�NB�HB�NB�ZB�`B�TB�HB�TB�NB�TB�NB�HB�NB�TB�NB�HB�HB�HB�NB�TB�NB�NB�HB�NB�NB�HB�NB�NB�`B�NB�NB�NB�NB�NB�NB�HB�NB�HB�NB�NB�NB�NB�NB�HB�NB�NB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�HB�HB�HB�HB�HB�HB�HB�NB�NB�NB�HB�NB�HB�NB�NB�NB�HB�NB�NB�HB�HB�HB�NB�HB�HB�NB�NB�NB�NB�NB�NB�NB�NB�NB�HB�HB�NB�NB�NB�HB�HB�HB�HB�HB�NB�HB�HB�BB�HB�HB�HB�BB�BB�BB�BB�HB�HB�BB�BB�HB�HB�HB�HB�HB�HB�HB�BB�HB�HB�HB�BB�HB�BB�BB�BB�HB�BB�HB�BB�BB�BB�BB�BB�HB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�;B�BB�BB�BB�;B�;B�;B�5B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�5B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�BB�;B�;B�;B�;B�;B�;B�;B�;B�;B�BB�BB�BB�BB�BB�BB�BB�BB�HB�HB�HB�NB�HB�HB�HB�HB�HB�HB�HB�HB�HB�HB�HB�NB�NB�ZB�fB�`B�fB�yB�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                       <#�<<#��<#�J<#�<#�+<#ڑ<#�<#�<$��<$��<#�<$ѩ<&��<$�t<#�<$e.<$f<$�<'<$r<#�<#�m<#ߜ<#�<$C�<,��<6��<%��<$�k<$�b<&�}<5\<6��<VjL<-��<#��<#��<#�{<#�<$4e<$!><#�<<$)
<#�$<#�&<#�e<#�J<#�I<#��<#�^<#�<#׺<#�<$�<%`�<$�3<#�$<$�.<#ڑ<$ѩ<)�L<%�n<#�<#��<$�e<#�&<$(<$|d<#��<#�$<#ا<#؄<$�<$I�<#�<#�8<#��<#�<#�4<#�^<$_�<#�<#�$<#��<#��<#ٛ<#�<#�<#�o<#�<%^�<$B�<$g�<%,#<%b<%�6<$�<#�&<#��<#��<#ا<#�8<#׺<$�<)�L<$��<$p<#ڑ<#�<#�"<#��<$	�<$�!<$)
<#�<#�8<#�)<$]h<#�<<&e<#�N<$p<$MO<#�<#�X<#�$<#�N<#�<$.<$XX<#�<#��<#��<#׺<#�&<#�<#�i<#�<#��<$@|<$�<#�<#�<#�<#�e<#�<<#�D<#�<#��<$��<#��<#�D<#��<#�^<#�U<#�J<#׺<#׎<#�c<#�i<#�i<#�<$8�<$�V<#�(<#�<#��<#�r<#�0<#�a<#�4<#�<$:�<%MY<$�<#�<$/<#�<#�<#ף<#��<#�&<#��<#؄<#׎<#�<#��<#��<#��<#؄<#�N<#�<#�<#�U<#�N<#��<$6�<#�D<#�8<#��<#�]<#�<#�D<#ا<#�<#�<#�<#�<#�<#��<#��<$�J<#�<#�<$�<$<<#�<#ף<#��<#��<#��<#�<#�!<#��<#�<$'<$�<$E<#�Q<$8�<#�<$G<#�<#�<#ף<#�r<#ۮ<#�*<#�"<#�<#�^<#�<#׎<#׺<#��<$J�<'�<&�*<#�g<#�<#��<#��<#ۮ<%e<#��<$<<#�o<$F9<#�<$��<%>�<&/<$�<$�<#�l<$�<#�<$�<#�<$
<#��<#�<#��<#ף<#�l<#�N<$�Q<$��<#��<#�M<$3U<#��<#��<#�<<$#(<%�<$��<#��<%�<#��<#�J<$r<$�k<&e�<#�l<#�<#�$<#�8<&
(<%<#�g<#�<,��<#��<#�{<%gB<#��<&�z<$3U<#�o<#�r<#�<$��<%�b<$�<#�8<#��<#�<#��<#׎<$'<#�W<&'<#��<#�C<#�<#��<$'<$C�<#��<#�J<#�<#׎<#�&<$?[<$�<#��<$��<$_�<$%<$k<#�<#�(<$E<#��<#�e<#�<#�o<#�{<#�<#��<#�e<#��<%�<'��<$�<$<$(<#�<#ף<#�<#��<#��<$F<#�r<#�c<$'<$��<#�<$<<#�]<#��<$y�<#��<$!><$��<$_�<%�<$7�<$�7<%��<%B�<%'<%��<#�<$�<#��<#�*<$�<$��<$�<$
<&�<$�<$
�<#�$<#��<%�<$��<$'<$:�<$g�<$�<#��<#�i<#��<#�X<#��<#��<#��<#��<#�I<#׺<#�c<#؄<$B�<$��<$�<#��<#�<#�N<$�<$	�<$_�<$
<#��<#��<$!><$��<%�V<$o�<#�N<$v<%�<$�w<$a<#�C<$<<$H�<$�<$�<#�!<#��<#�r<#׺<#�<#�<#ף<#��<#�4<#�o<#�<#�	<#�E<#��<#�
<#�
<#ף<#�<#�<#��<#�<#�J<#��<#�<#�{<#ا<#�C<#�<$
�<#�o<#�+<#��<#��<#��<#�D<#ߜ<#�&<#��<#ߜ<#�o<#�<#ף<#ٛ<#��<#��<#ڑ<#��<#�+<#��<#׺<#�<#�<#�<#��<#��<#�"<#��<#ޫ<#�0<#��<#��<#�$<$*<#��<#�<#�{<#�<#�
<#�i<#�<#�i<#�
<#�X<#�{<#�i<#׎<#�<#�<#�<#�<#�X<#�<#��<#�<#ף<#׎<#׺<#�<#؄<#�D<#�*<#�<#�<#�<#�<#׎<#ף<#�^<#ף<#�<#�<#�<#�{<#�<#�i<#�
<#�<#�i<#ף<#�^<#�E<#��<#��<#ߜ<#�^<#�<#�
<#�X<#׎<#ף<#�i<#�i<#�<#�
<#�<#�
<#�<#�<#�<#�{<#�<#�<#�
<#�{<#��<#��<#��<#ף<#��<#�o<#�8<#�<#׺<#ޫ<#ڑ<#�<#�i<#��<#׺<#ף<#׺<#׎<#ٛ<#�<#�N<#�<#�D<#�8<#��<#�{<#�
<#ף<#��<#׺<#�D<#�<#�l<#�<$XX<$�<#�D<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285520181106012855QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285520181106012855QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                