CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS     N_CALIB       	N_HISTORY             
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
resolution        =���   axis      Z          <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   T�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       Z�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   r�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       x�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       ��   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       ��   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o       ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                   ��   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      �   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                      PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o      !   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` 9   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   9l   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   ?l   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   El   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T Kl   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   K�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   K�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   K�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   K�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � K�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   L`   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   L|   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    L�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        L�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        L�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       L�   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    L�Argo profile    3.1 1.2 19500101000000  20181106012856  20230209090617  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               O   OAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�v�'��@�v�'��11  @�v�8�P@�v�8�P@N�QC�rq@N�QC�rq�;�EM��`�;�EM��`11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�\)@   @@  @�  @�G�@�G�@�  @��RA  A$z�A?\)A`  A�Q�A�Q�A���A�  A�  AУ�A�Q�A�  B   B  B  B  B   B(  B0  B7�
B@  BG33BO�BX  B_�
Bg�Bp  Bx(�B�(�B�(�B�(�B�(�B�  B��B�  B�  B�{B�  B��B��B��B��
B�B��B�  B��B��
B�  B�=qB�(�B�(�B�{B�{B�{B��B��
B�  B��B��B�  B��C�HC��C  C��C
  C{C�C  C�C  C{C  C�C�C�C��C"
=C$  C%��C'�HC)�HC,  C.  C0  C2
=C4  C5��C7�C9�HC<  C>{C@{CB
=CD  CF
=CH�CJ
=CK��CM�CO��CR  CS�CU��CX  CY��C[��C^  C`
=Cb  Cd
=Cf
=Ch  Cj
=Cl{Cn
=Cp  Cq��Cs�Cu��Cx  Cz
=C|{C~  C�C�  C�
=C�C�  C���C���C�  C���C�  C�C���C���C���C���C���C���C���C�  C�  C���C���C�  C�C�C���C���C�  C���C���C���C���C��C���C���C���C���C�  C�  C�C�C�  C�C�  C���C���C���C�  C�
=C�C���C�  C�
=C�  C���C�  C���C���C�C�
=C���C���C�  C�  C���C���C�  C�C�  C�C�C�  C���C���C�C�  C���C���C���C�C�  C�  C���C���C�  C���C���C���C���C���C���C�  C�C�  C���C���C���C�  C���C���C���C���C�  C�
=C�C�  C���C�  C�  C���C�C�\C�C���C���C���C���C���C���C���C�  C�  C�C�  C���C���C���C���C��D xRD ��D}qD�D� D��D� D�D}qD��D}qD  D� D�D��D�qD}qD�qD	}qD	�qD
}qD
�qDz�D�D� D  D�DD}qD�qD� D�qD� D�D}qD  D�D�D��DD}qD�RDxRD�RDz�D�qD� D�D�DD}qD�RDxRD�qD� DD� D�qD� D�D}qD��D� D �D }qD �RD!z�D!�qD"� D#  D#��D$�D$��D%�D%� D%��D&� D'D'��D(  D(� D(�qD)��D*D*� D+�D+��D,�D,� D,�qD-��D.  D.� D.�qD/}qD0  D0}qD1  D1��D2  D2}qD2�qD3}qD4  D4��D5�D5}qD6  D6�D7  D7z�D7�qD8� D9  D9� D9�qD:� D;  D;z�D;�qD<� D=D=}qD=��D>z�D>�qD?��D@D@��DA  DA� DB  DB��DC�DCz�DC�qDD� DD�qDE}qDE�qDFz�DG  DG��DG�qDH}qDH�qDIz�DI�qDJ��DJ�qDKz�DL  DL�DM  DM� DN�DN��DO�DO� DP  DP��DQDQ� DQ�qDR��DS�DS� DS�qDTz�DU  DU}qDV  DV��DW�DW��DW��DXz�DY  DY}qDZ  DZ��D[�D[�D[�qD\z�D\�qD]� D]�qD^� D_  D_��D_�qD`z�D`�qDa� Da��DbxRDb�RDc� DdDd��De�De��Df�Df��DgDg��Dh  Dhz�Dh��Di� Dj  DjxRDj��Dkz�Dk�qDl}qDm  Dm��Dn�Dn}qDo�Do�Dp  Dpz�Dq�Dq� Dq�qDr��Dr�qDsz�Dt�Dt� Dt�qDu�Dv  Dv� DwDw� Dx  Dx}qDx�qDy� Dy�RDz}qD{  D{��D|�D|}qD}  D}� D}�qD~xRD~�qD� D�HD�AHD��HD�D���D�@ D�� D��HD��D�=qD�}qD��qD��)D�@ D�� D�� D��D�AHD�~�D�D�HD�AHD��HD�� D���D�=qD�~�D��HD�HD�@ D���D��HD��D�AHD�~�D�� D��qD�@ D�~�D�D�HD�@ D�� D��qD��qD�=qD�}qD��qD�HD�B�D�~�D���D�HD�AHD��HD��HD���D�@ D��HD�D��D�AHD�~�D���D��D�AHD�~�D��qD��qD�=qD�� D��HD���D�>�D�� D�� D�  D�@ D�~�D��qD���D�AHD�� D�� D�  D�=qD�}qD�� D�HD�AHD��HD��HD�  D�=qD�� D�D��D�AHD�~�D��HD�HD�>�D�� D��HD�HD�B�D��HD��HD�  D�@ D�~�D��)D�  D�AHD�~�D��)D��qD�>�D�� D��HD�  D�>�D��HD��HD�HD�@ D�}qD���D���D�=qD�� D�� D���D�B�D���D��HD��qD�<)D�z�D�� D��D�B�D��HD�� D���D�>�D�}qD��HD��D�>�D�~�D��HD�  D�=qD�� D���D��qD�AHD���D���D��qD�>�D��HD�D�HD�B�D�� D���D�  D�>�D�~�D���D���D�>�D��HD�D�HD�>�D�~�D���D���D�=qD�~�D�� D�HD�AHD�� D�� D�  D�AHD���D��HD�  D�@ D��HD�D���D�<)D�}qD�� D��D�@ D�}qD���D���D�=qD�� D�D�HD�@ D�� D�� D�  D�>�D�� D��HD�  D�>�D�� D���D���D�>�D�~�D���D���D�>�D�� D�� D�  D�B�D�� D���D�  D�@ D��HD�� D��qD�@ D��HD��HD�  D�>�D�}qD��qD���D�@ D�� D���D���D�>�D�� D�� D��?aG�?u?�\)?���?���?�Q�?\?�
=?�ff?��@�@��@z�@��@!G�@+�@333@:�H@E�@L��@Tz�@^�R@c�
@n{@s33@z�H@�G�@��@���@�{@��@�@���@�p�@�G�@��@���@��@�\)@�z�@�Q�@�(�@�  @��
@Ǯ@˅@�\)@�33@�Q�@�(�@�  @��@���@���@��@�z�@���@��RA�A�
AffAQ�A
=qA(�A�RAG�A�
AA�A��A(�A{A ��A#33A%A(Q�A*=qA,��A.�RA0��A333A5A8Q�A:=qA<(�A>{A@��AC33AE�AG�AI��AL(�AN{AP��AR�\AU�AW
=AX��A[�A^{A`��Ab�\Adz�Ag
=Ai��Ak�Amp�Ao\)Aq�Atz�AvffAx��Az�HA|��A\)A���A��A��HA�(�A��A�ffA�\)A�Q�A���A��HA��
A��A�{A�\)A���A���A�=qA��A���A�{A�
=A�  A�G�A��\A��A�z�A�A��RA�  A���A��A��A�z�A�p�A��RA��A���A�=qA�33A�(�A�p�A�ffA�\)A���A��A��HA��
A���A�{A�\)A�Q�A�G�A��\A��
A���A�A�
=A�  A�G�A�=qAÅA�z�A�A�
=AǮA���Aʏ\A˅A�(�A�A�
=A�  A���A�=qAӅA���A�AָRA�  A�G�A�=qA�33A�z�A�A޸RA߮A���A�=qA�33A�(�A�p�A�RA�\)A�Q�A陚A��HA��
A���A�A�
=A�Q�A�G�A��A�33A�z�A�p�A�ffA�\)A���A��A��HA��A���A�{A�
=B   B ��B�B��B{B�RB\)B�B(�B��Bp�B�BffB�HB�B  Bz�B	�B	��B
=qB
�RB33B�
BQ�B��BG�B�B�\B
=B\)B  B��B�B��B{B�RB\)B�B(�B��Bp�BB=qB�HB\)B�
BQ�B��Bp�BB=qB�RB\)B�B  Bz�B�Bp�BBffB�HB33B�B   B z�B ��B!�B!B"=qB"�\B"�HB#\)B#�
B$(�B$��B$��B%p�B%�B&=qB&�\B'
=B'�B'�
B((�B(��B)�B)��B)�B*=qB*�HB+33B+�B+�
B,Q�B,��B-�B-p�B-�B.=qB.�RB/
=B/�B0  B0(�B0��B1�B1G�B1B2=qB2�\B2�HB333B3�B4(�B4z�B4��B5�B5��B6{B6ffB6�RB733B7�B7�
B8(�B8��B8��B9p�B9B:=qB:�\B:�HB;33B;�B<  B<Q�B<��B=�B=p�B=B>=qB>�RB?
=B?\)B?�
B@(�B@z�B@��BAp�BABB{BBffBB�HBC33BC�BC�
BDQ�BD��BD��BEG�BEBF=qBFffBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       ?�\)@   @@  @�  @�G�@�G�@�  @��RA  A$z�A?\)A`  A�Q�A�Q�A���A�  A�  AУ�A�Q�A�  B   B  B  B  B   B(  B0  B7�
B@  BG33BO�BX  B_�
Bg�Bp  Bx(�B�(�B�(�B�(�B�(�B�  B��B�  B�  B�{B�  B��B��B��B��
B�B��B�  B��B��
B�  B�=qB�(�B�(�B�{B�{B�{B��B��
B�  B��B��B�  B��C�HC��C  C��C
  C{C�C  C�C  C{C  C�C�C�C��C"
=C$  C%��C'�HC)�HC,  C.  C0  C2
=C4  C5��C7�C9�HC<  C>{C@{CB
=CD  CF
=CH�CJ
=CK��CM�CO��CR  CS�CU��CX  CY��C[��C^  C`
=Cb  Cd
=Cf
=Ch  Cj
=Cl{Cn
=Cp  Cq��Cs�Cu��Cx  Cz
=C|{C~  C�C�  C�
=C�C�  C���C���C�  C���C�  C�C���C���C���C���C���C���C���C�  C�  C���C���C�  C�C�C���C���C�  C���C���C���C���C��C���C���C���C���C�  C�  C�C�C�  C�C�  C���C���C���C�  C�
=C�C���C�  C�
=C�  C���C�  C���C���C�C�
=C���C���C�  C�  C���C���C�  C�C�  C�C�C�  C���C���C�C�  C���C���C���C�C�  C�  C���C���C�  C���C���C���C���C���C���C�  C�C�  C���C���C���C�  C���C���C���C���C�  C�
=C�C�  C���C�  C�  C���C�C�\C�C���C���C���C���C���C���C���C�  C�  C�C�  C���C���C���C���C��D xRD ��D}qD�D� D��D� D�D}qD��D}qD  D� D�D��D�qD}qD�qD	}qD	�qD
}qD
�qDz�D�D� D  D�DD}qD�qD� D�qD� D�D}qD  D�D�D��DD}qD�RDxRD�RDz�D�qD� D�D�DD}qD�RDxRD�qD� DD� D�qD� D�D}qD��D� D �D }qD �RD!z�D!�qD"� D#  D#��D$�D$��D%�D%� D%��D&� D'D'��D(  D(� D(�qD)��D*D*� D+�D+��D,�D,� D,�qD-��D.  D.� D.�qD/}qD0  D0}qD1  D1��D2  D2}qD2�qD3}qD4  D4��D5�D5}qD6  D6�D7  D7z�D7�qD8� D9  D9� D9�qD:� D;  D;z�D;�qD<� D=D=}qD=��D>z�D>�qD?��D@D@��DA  DA� DB  DB��DC�DCz�DC�qDD� DD�qDE}qDE�qDFz�DG  DG��DG�qDH}qDH�qDIz�DI�qDJ��DJ�qDKz�DL  DL�DM  DM� DN�DN��DO�DO� DP  DP��DQDQ� DQ�qDR��DS�DS� DS�qDTz�DU  DU}qDV  DV��DW�DW��DW��DXz�DY  DY}qDZ  DZ��D[�D[�D[�qD\z�D\�qD]� D]�qD^� D_  D_��D_�qD`z�D`�qDa� Da��DbxRDb�RDc� DdDd��De�De��Df�Df��DgDg��Dh  Dhz�Dh��Di� Dj  DjxRDj��Dkz�Dk�qDl}qDm  Dm��Dn�Dn}qDo�Do�Dp  Dpz�Dq�Dq� Dq�qDr��Dr�qDsz�Dt�Dt� Dt�qDu�Dv  Dv� DwDw� Dx  Dx}qDx�qDy� Dy�RDz}qD{  D{��D|�D|}qD}  D}� D}�qD~xRD~�qD� D�HD�AHD��HD�D���D�@ D�� D��HD��D�=qD�}qD��qD��)D�@ D�� D�� D��D�AHD�~�D�D�HD�AHD��HD�� D���D�=qD�~�D��HD�HD�@ D���D��HD��D�AHD�~�D�� D��qD�@ D�~�D�D�HD�@ D�� D��qD��qD�=qD�}qD��qD�HD�B�D�~�D���D�HD�AHD��HD��HD���D�@ D��HD�D��D�AHD�~�D���D��D�AHD�~�D��qD��qD�=qD�� D��HD���D�>�D�� D�� D�  D�@ D�~�D��qD���D�AHD�� D�� D�  D�=qD�}qD�� D�HD�AHD��HD��HD�  D�=qD�� D�D��D�AHD�~�D��HD�HD�>�D�� D��HD�HD�B�D��HD��HD�  D�@ D�~�D��)D�  D�AHD�~�D��)D��qD�>�D�� D��HD�  D�>�D��HD��HD�HD�@ D�}qD���D���D�=qD�� D�� D���D�B�D���D��HD��qD�<)D�z�D�� D��D�B�D��HD�� D���D�>�D�}qD��HD��D�>�D�~�D��HD�  D�=qD�� D���D��qD�AHD���D���D��qD�>�D��HD�D�HD�B�D�� D���D�  D�>�D�~�D���D���D�>�D��HD�D�HD�>�D�~�D���D���D�=qD�~�D�� D�HD�AHD�� D�� D�  D�AHD���D��HD�  D�@ D��HD�D���D�<)D�}qD�� D��D�@ D�}qD���D���D�=qD�� D�D�HD�@ D�� D�� D�  D�>�D�� D��HD�  D�>�D�� D���D���D�>�D�~�D���D���D�>�D�� D�� D�  D�B�D�� D���D�  D�@ D��HD�� D��qD�@ D��HD��HD�  D�>�D�}qD��qD���D�@ D�� D���D���D�>�D�� D�� D��?aG�?u?�\)?���?���?�Q�?\?�
=?�ff?��@�@��@z�@��@!G�@+�@333@:�H@E�@L��@Tz�@^�R@c�
@n{@s33@z�H@�G�@��@���@�{@��@�@���@�p�@�G�@��@���@��@�\)@�z�@�Q�@�(�@�  @��
@Ǯ@˅@�\)@�33@�Q�@�(�@�  @��@���@���@��@�z�@���@��RA�A�
AffAQ�A
=qA(�A�RAG�A�
AA�A��A(�A{A ��A#33A%A(Q�A*=qA,��A.�RA0��A333A5A8Q�A:=qA<(�A>{A@��AC33AE�AG�AI��AL(�AN{AP��AR�\AU�AW
=AX��A[�A^{A`��Ab�\Adz�Ag
=Ai��Ak�Amp�Ao\)Aq�Atz�AvffAx��Az�HA|��A\)A���A��A��HA�(�A��A�ffA�\)A�Q�A���A��HA��
A��A�{A�\)A���A���A�=qA��A���A�{A�
=A�  A�G�A��\A��A�z�A�A��RA�  A���A��A��A�z�A�p�A��RA��A���A�=qA�33A�(�A�p�A�ffA�\)A���A��A��HA��
A���A�{A�\)A�Q�A�G�A��\A��
A���A�A�
=A�  A�G�A�=qAÅA�z�A�A�
=AǮA���Aʏ\A˅A�(�A�A�
=A�  A���A�=qAӅA���A�AָRA�  A�G�A�=qA�33A�z�A�A޸RA߮A���A�=qA�33A�(�A�p�A�RA�\)A�Q�A陚A��HA��
A���A�A�
=A�Q�A�G�A��A�33A�z�A�p�A�ffA�\)A���A��A��HA��A���A�{A�
=B   B ��B�B��B{B�RB\)B�B(�B��Bp�B�BffB�HB�B  Bz�B	�B	��B
=qB
�RB33B�
BQ�B��BG�B�B�\B
=B\)B  B��B�B��B{B�RB\)B�B(�B��Bp�BB=qB�HB\)B�
BQ�B��Bp�BB=qB�RB\)B�B  Bz�B�Bp�BBffB�HB33B�B   B z�B ��B!�B!B"=qB"�\B"�HB#\)B#�
B$(�B$��B$��B%p�B%�B&=qB&�\B'
=B'�B'�
B((�B(��B)�B)��B)�B*=qB*�HB+33B+�B+�
B,Q�B,��B-�B-p�B-�B.=qB.�RB/
=B/�B0  B0(�B0��B1�B1G�B1B2=qB2�\B2�HB333B3�B4(�B4z�B4��B5�B5��B6{B6ffB6�RB733B7�B7�
B8(�B8��B8��B9p�B9B:=qB:�\B:�HB;33B;�B<  B<Q�B<��B=�B=p�B=B>=qB>�RB?
=B?\)B?�
B@(�B@z�B@��BAp�BABB{BBffBB�HBC33BC�BC�
BDQ�BD��BD��BEG�BEBF=qBFffBF�HBG\)BG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A%"�A%"�A%�A%"�A%"�A%"�A%&�A%?}A%S�A%`BA%&�A#�
A!33A��A��A"�AQ�A�AC�A�DA��A�
A7LAz�A�^AoA�mA�yA�A	�^A�A�A�PA�/A�FA�#@�+@��@�o@�@��`@�9X@�t�@��@�@��@���@��#@�1@�
=@�j@��@���@띲@�P@�@�F@��@�|�@��@�ȴ@�!@��@�R@�=q@���@��@��@��@���@��@��@���@���@旍@���@�7L@�V@�Q�@��m@��@��;@��m@��;@�@�+@���@�!@�=q@��@��@��T@�@��@�X@���@��D@�Q�@� �@��m@�1@��/@�G�@�p�@�j@��@�z�@�z�@���@���@���@�G�@��@��@���@�^@�^@���@���@�{@��T@�x�@�V@�Ĝ@�V@���@��/@��u@�/@�V@�ff@ᙚ@�&�@��`@��u@���@�/@�hs@�7@�-@��T@��@�n�@�+@�|�@�M�@�p�@��@��`@��@�9X@�(�@�1@���@��@��@��m@�1@��u@�9X@�I�@��u@���@��@�Ĝ@��@�  @ߕ�@�33@�"�@��@�o@���@���@އ+@�v�@�V@��@���@��T@݉7@�/@�&�@�O�@�hs@�-@�K�@�j@��@���@���@���@��D@�j@�A�@�9X@�  @ߥ�@�t�@�dZ@�S�@�"�@���@�=q@��@��T@��#@���@�@ݺ^@ݑh@�G�@�V@���@ܴ9@�Z@��@۾w@�l�@�@���@١�@���@׶F@ׅ@��H@ָR@���@��H@��H@��@��@��H@��y@���@�@�33@�|�@׾w@�ƨ@���@��;@��m@׍P@��y@֗�@�%@�Z@Ӿw@�dZ@��@�"�@��@�=q@с@�7L@�&�@���@���@д9@Гu@�r�@�I�@�9X@� �@�1@��@�j@�1@ϕ�@�C�@θR@Ώ\@�x�@�Ĝ@��;@˾w@˕�@�+@�
=@��@ʰ!@�{@�O�@��/@�
=@Ų-@�&�@��`@�I�@�ƨ@�+@�o@�ȴ@\@�n�@�ff@�@��#@��@���@�x�@�p�@�hs@�7L@�Ĝ@�Q�@�  @���@�1@�V@�hs@��#@���@��@�I�@��w@��@���@�$�@���@�O�@���@�r�@�9X@�b@���@���@�t�@�33@�~�@�=q@�5?@�5?@�-@���@���@���@�`B@�V@��u@�A�@���@�C�@��@���@�V@��#@��@��/@��@�r�@�Z@��
@��F@���@��@�K�@��H@���@���@���@�V@�J@���@��@�&�@��/@�r�@���@��H@��R@��!@��+@�ff@��@���@���@���@���@��7@��@��@�`B@�O�@�/@�V@��`@��@��@��@�z�@�r�@�Q�@���@�l�@�C�@�S�@��@�{@���@��@�`B@�X@�?}@��@���@��`@��/@���@���@��@�1'@��@�+@��!@�^5@�M�@�M�@�=q@��@��@���@�@���@�p�@�`B@�/@�V@�%@���@��9@��u@�j@�A�@� �@��;@��P@�C�@�
=@�ȴ@���@�~�@�V@�J@�J@���@��#@��^@���@�x�@��@��@��u@�r�@�bN@�A�@�1'@�1@�  @��@��
@��w@���@�t�@�\)@�C�@�"�@��@�ȴ@��!@��\@�^5@�-@��@�@�@��^@�p�@��@���@��j@��@�bN@�Q�@�9X@�1'@� �@�b@�1@��@���@���@�l�@�"�@��y@��\@�~�@�^5@�^5@�{@��@���@�@��^@��-@��h@�x�@�`B@�G�@�&�@��@��/@���@��/@��/@���@��j@��j@��@���@�r�@�I�@�A�@�9X@�9X@�1'@��@�b@�1@��w@�l�@�@��y@���@�{@�J@���@��#@��^@���@���@��@�p�@�O�@�O�@�/@���@���@���@�Ĝ@��u@�(�@��w@��F@��P@�v�@�M�@���@���@�?}@�V@��`@��@�1'@�  @��
@�ƨ@�ƨ@�ƨ@���@�
=@���@���@��+@��+@��+@�~�@�n�@�v�@�n�@�ff@�M�@�-@��@��T@�@���@���@��h@�/@���@�A�@��@��m@�\)@�"�@���@�ȴ@���@�n�@���@��7@�7L@���@�Z@�A�@�b@���@��P@�l�@�K�@��@���@��!@���@�~�@�^5@�V@�V@�V@�M�@�E�@�M�@�M�@�M�@�M�@�M�@�=q@�=q@�=q@�E�@�M�@�M�@�V@�V@�V@�V@�^5@�^5@�^5@�n�@�v�@�v�@�v�@�v�@�n�@�n�@�ff@�^5@�M�@�E�@�5?@�J@���@��^@��^@��-@��-@���@��h@�x�@�`B@�O�@�&�@�V@�A�@�1'@�1'@�1@��F@���@��;@��;@��m@��
@��w@��F@��w@�K�@�33@�
=@��@���@��!@��\@�V@�J@�J@�@�J@�J@�{@�{@�{@��@��@��@�$�@��@��@��@�-@�-@�-@�-@�5?@�5?@�5?@�-@��@��T@���@�@��^@���@���@���@���@���@���@��@�`B@�?}@�G�@�O�@�G�@�/@�&�@��@�%@���@�%@�V@���@��@��@��/@��/@��/@���@��9@��9@��9@��9@��9@��j@��j@��9@��j@��j@��j@��j@��j@�Ĝ@��@��@�Ĝ@�Ĝ@���@�Ĝ@�Ĝ@���@���@���@��/@���@���@��/A%"�A%"�A%"�A%�A%"�A%"�A%"�A%&�A%"�A%&�A%&�A%&�A%"�A%"�A%�A%�A%�A%�A%�A%"�A%"�A%�A%"�A%�A%�A%"�A%"�A%"�A%&�A%"�A%&�A%&�A%&�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%�A%"�A%&�A%"�A%"�A%"�A%/A%+A%"�A%&�A%&�A%33A%C�A%O�A%S�A%S�A%S�A%S�A%O�A%S�A%O�A%O�A%O�A%XA%dZA%dZA%dZA%`BA%dZA%`BA%S�A%XA%XA%`BA%\)A%O�A%S�A%K�A%?}A%?}A%C�A%+A%VA%
=A$��A$��A$��A$�`A$�HA$�`A$��A$��A$��A$r�A$ZA$-A#��A#C�A"�HA"��A"��A"��A"~�A"E�A"$�A"1A!�TA!ƨA!��A!t�A!;dA �/A ~�A jA ffA =qA�A�-AdZAC�A33A�A�A��A��A~�A^5A5?AbA  A�A��A��A�7AXAC�A
=A��A�!A��Az�A1'A�A�^A��Ax�A\)AXAK�AG�AC�AC�A7LA+A�A%A��A��A�/AĜA�9A�AbNAZAZAM�AM�AM�AE�AA�A=qA5?A5?A1'A-A-A$�A�A{AJAA�A�A�;A��A�^A�A��A��A��A��A��A�PAl�AK�A?}A33A/A+A"�A%A�A��A�9A��A��A��A��A��A�uA�\A�DA�Az�Av�An�A^5AI�A(�A�AbA1A1A��A�TA��A��A�wA�FA�A��A�AG�A
=A�jA�\Az�Ar�AZA5?A�mA��At�AS�A?}A��A�!A�+AZA5?A��A��A�PA�Ax�A?}A+A�`A�!A�DAbNAVAE�A �AbAƨA��AC�A�DA�A�#A�#A�;A�#A�
A��AƨA��A��A��A��A��A��A��AƨA�wA�FA�FA�FA�-A�FA�wAAA�-A��A�PA�AhsAC�A&�A�AVAVAVAoAoAoAoAoAoAVA
=AA��A��A��A�A�HA��AĜA��AjA1'A��A��A�wA�wA�-A�A��A�7A�A�At�A`BAS�A;dA/A+A+A"�A�A"�A�A�A%A��A�A�A�HA�HA�/A��A�!A��A��A�\A�+A~�Av�AI�A�A��A�#AƨA�FA��A�AS�A33A
=A
��A
�DA
=qA	��A	�
A	A	��A	�#A	�TA	�;A	�;A	�#A	�TA	��A
JA
bA
{A
�A
�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       A%"�A%"�A%�A%"�A%"�A%"�A%&�A%?}A%S�A%`BA%&�A#�
A!33A��A��A"�AQ�A�AC�A�DA��A�
A7LAz�A�^AoA�mA�yA�A	�^A�A�A�PA�/A�FA�#@�+@��@�o@�@��`@�9X@�t�@��@�@��@���@��#@�1@�
=@�j@��@���@띲@�P@�@�F@��@�|�@��@�ȴ@�!@��@�R@�=q@���@��@��@��@���@��@��@���@���@旍@���@�7L@�V@�Q�@��m@��@��;@��m@��;@�@�+@���@�!@�=q@��@��@��T@�@��@�X@���@��D@�Q�@� �@��m@�1@��/@�G�@�p�@�j@��@�z�@�z�@���@���@���@�G�@��@��@���@�^@�^@���@���@�{@��T@�x�@�V@�Ĝ@�V@���@��/@��u@�/@�V@�ff@ᙚ@�&�@��`@��u@���@�/@�hs@�7@�-@��T@��@�n�@�+@�|�@�M�@�p�@��@��`@��@�9X@�(�@�1@���@��@��@��m@�1@��u@�9X@�I�@��u@���@��@�Ĝ@��@�  @ߕ�@�33@�"�@��@�o@���@���@އ+@�v�@�V@��@���@��T@݉7@�/@�&�@�O�@�hs@�-@�K�@�j@��@���@���@���@��D@�j@�A�@�9X@�  @ߥ�@�t�@�dZ@�S�@�"�@���@�=q@��@��T@��#@���@�@ݺ^@ݑh@�G�@�V@���@ܴ9@�Z@��@۾w@�l�@�@���@١�@���@׶F@ׅ@��H@ָR@���@��H@��H@��@��@��H@��y@���@�@�33@�|�@׾w@�ƨ@���@��;@��m@׍P@��y@֗�@�%@�Z@Ӿw@�dZ@��@�"�@��@�=q@с@�7L@�&�@���@���@д9@Гu@�r�@�I�@�9X@� �@�1@��@�j@�1@ϕ�@�C�@θR@Ώ\@�x�@�Ĝ@��;@˾w@˕�@�+@�
=@��@ʰ!@�{@�O�@��/@�
=@Ų-@�&�@��`@�I�@�ƨ@�+@�o@�ȴ@\@�n�@�ff@�@��#@��@���@�x�@�p�@�hs@�7L@�Ĝ@�Q�@�  @���@�1@�V@�hs@��#@���@��@�I�@��w@��@���@�$�@���@�O�@���@�r�@�9X@�b@���@���@�t�@�33@�~�@�=q@�5?@�5?@�-@���@���@���@�`B@�V@��u@�A�@���@�C�@��@���@�V@��#@��@��/@��@�r�@�Z@��
@��F@���@��@�K�@��H@���@���@���@�V@�J@���@��@�&�@��/@�r�@���@��H@��R@��!@��+@�ff@��@���@���@���@���@��7@��@��@�`B@�O�@�/@�V@��`@��@��@��@�z�@�r�@�Q�@���@�l�@�C�@�S�@��@�{@���@��@�`B@�X@�?}@��@���@��`@��/@���@���@��@�1'@��@�+@��!@�^5@�M�@�M�@�=q@��@��@���@�@���@�p�@�`B@�/@�V@�%@���@��9@��u@�j@�A�@� �@��;@��P@�C�@�
=@�ȴ@���@�~�@�V@�J@�J@���@��#@��^@���@�x�@��@��@��u@�r�@�bN@�A�@�1'@�1@�  @��@��
@��w@���@�t�@�\)@�C�@�"�@��@�ȴ@��!@��\@�^5@�-@��@�@�@��^@�p�@��@���@��j@��@�bN@�Q�@�9X@�1'@� �@�b@�1@��@���@���@�l�@�"�@��y@��\@�~�@�^5@�^5@�{@��@���@�@��^@��-@��h@�x�@�`B@�G�@�&�@��@��/@���@��/@��/@���@��j@��j@��@���@�r�@�I�@�A�@�9X@�9X@�1'@��@�b@�1@��w@�l�@�@��y@���@�{@�J@���@��#@��^@���@���@��@�p�@�O�@�O�@�/@���@���@���@�Ĝ@��u@�(�@��w@��F@��P@�v�@�M�@���@���@�?}@�V@��`@��@�1'@�  @��
@�ƨ@�ƨ@�ƨ@���@�
=@���@���@��+@��+@��+@�~�@�n�@�v�@�n�@�ff@�M�@�-@��@��T@�@���@���@��h@�/@���@�A�@��@��m@�\)@�"�@���@�ȴ@���@�n�@���@��7@�7L@���@�Z@�A�@�b@���@��P@�l�@�K�@��@���@��!@���@�~�@�^5@�V@�V@�V@�M�@�E�@�M�@�M�@�M�@�M�@�M�@�=q@�=q@�=q@�E�@�M�@�M�@�V@�V@�V@�V@�^5@�^5@�^5@�n�@�v�@�v�@�v�@�v�@�n�@�n�@�ff@�^5@�M�@�E�@�5?@�J@���@��^@��^@��-@��-@���@��h@�x�@�`B@�O�@�&�@�V@�A�@�1'@�1'@�1@��F@���@��;@��;@��m@��
@��w@��F@��w@�K�@�33@�
=@��@���@��!@��\@�V@�J@�J@�@�J@�J@�{@�{@�{@��@��@��@�$�@��@��@��@�-@�-@�-@�-@�5?@�5?@�5?@�-@��@��T@���@�@��^@���@���@���@���@���@���@��@�`B@�?}@�G�@�O�@�G�@�/@�&�@��@�%@���@�%@�V@���@��@��@��/@��/@��/@���@��9@��9@��9@��9@��9@��j@��j@��9@��j@��j@��j@��j@��j@�Ĝ@��@��@�Ĝ@�Ĝ@���@�Ĝ@�Ĝ@���@���@���@��/@���@���@��/A%"�A%"�A%"�A%�A%"�A%"�A%"�A%&�A%"�A%&�A%&�A%&�A%"�A%"�A%�A%�A%�A%�A%�A%"�A%"�A%�A%"�A%�A%�A%"�A%"�A%"�A%&�A%"�A%&�A%&�A%&�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%"�A%�A%"�A%&�A%"�A%"�A%"�A%/A%+A%"�A%&�A%&�A%33A%C�A%O�A%S�A%S�A%S�A%S�A%O�A%S�A%O�A%O�A%O�A%XA%dZA%dZA%dZA%`BA%dZA%`BA%S�A%XA%XA%`BA%\)A%O�A%S�A%K�A%?}A%?}A%C�A%+A%VA%
=A$��A$��A$��A$�`A$�HA$�`A$��A$��A$��A$r�A$ZA$-A#��A#C�A"�HA"��A"��A"��A"~�A"E�A"$�A"1A!�TA!ƨA!��A!t�A!;dA �/A ~�A jA ffA =qA�A�-AdZAC�A33A�A�A��A��A~�A^5A5?AbA  A�A��A��A�7AXAC�A
=A��A�!A��Az�A1'A�A�^A��Ax�A\)AXAK�AG�AC�AC�A7LA+A�A%A��A��A�/AĜA�9A�AbNAZAZAM�AM�AM�AE�AA�A=qA5?A5?A1'A-A-A$�A�A{AJAA�A�A�;A��A�^A�A��A��A��A��A��A�PAl�AK�A?}A33A/A+A"�A%A�A��A�9A��A��A��A��A��A�uA�\A�DA�Az�Av�An�A^5AI�A(�A�AbA1A1A��A�TA��A��A�wA�FA�A��A�AG�A
=A�jA�\Az�Ar�AZA5?A�mA��At�AS�A?}A��A�!A�+AZA5?A��A��A�PA�Ax�A?}A+A�`A�!A�DAbNAVAE�A �AbAƨA��AC�A�DA�A�#A�#A�;A�#A�
A��AƨA��A��A��A��A��A��A��AƨA�wA�FA�FA�FA�-A�FA�wAAA�-A��A�PA�AhsAC�A&�A�AVAVAVAoAoAoAoAoAoAVA
=AA��A��A��A�A�HA��AĜA��AjA1'A��A��A�wA�wA�-A�A��A�7A�A�At�A`BAS�A;dA/A+A+A"�A�A"�A�A�A%A��A�A�A�HA�HA�/A��A�!A��A��A�\A�+A~�Av�AI�A�A��A�#AƨA�FA��A�AS�A33A
=A
��A
�DA
=qA	��A	�
A	A	��A	�#A	�TA	�;A	�;A	�#A	�TA	��A
JA
bA
{A
�A
�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B��B��B��B��B��B��B��B��B��B��B��B��B�1Bt�BhsBbNBaHB`BB_;B_;BbNBaHB\)BiyBo�Bt�B�B�DB�=B�uB�-BǮB��B��B��B�
B�;B�TB�TB�NB�HB�sB�yB�yB�yB�yB�sB�sB�sB�sB�yB�yB�yB�yB�sB�yB�yB�yB�mB�mB�mB�fB�fB�fB�TB�NB�;B�#B�#B�)B�)B�/B�/B�/B�5B�;B�BB�BB�BB�;B�;B�;B�;B�;B�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�;B�5B�5B�;B�;B�NB�`B�yB�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BBB  B��B��B��BBB%B+B1B
=BbBhB�B�BuB\BPBPBPBPBPBVBVBVBVBbBhB�B�B�B�B�B�B�B�B�B�B{BuBuBuBoBoBhBhBhBbBbBbBbBhBoB{B�B�B"�B)�B,B.B/B0!B0!B0!B0!B0!B0!B0!B0!B0!B/B/B.B.B.B/B.B.B/B/B/B/B/B/B/B/B/B.B,B+B'�B&�B$�B"�B"�B!�B!�B#�B$�B%�B%�B%�B&�B'�B(�B)�B,B/B1'B2-B2-B33B49B6FB6FB5?B2-B1'B0!B0!B0!B0!B/B,B)�B(�B(�B(�B(�B(�B(�B+B-B-B-B-B-B1'B2-B2-B1'B0!B/B/B/B1'B1'B1'B2-B1'B0!B/B-B)�B'�B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B#�B%�B(�B(�B(�B(�B'�B'�B'�B'�B&�B&�B%�B%�B%�B%�B%�B$�B$�B#�B#�B"�B#�B"�B"�B"�B"�B!�B!�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BoBbBVBVBPBPBPBPBPBPBPBPBPBPBPBPBPBPBPBJBPBPBPBPBPBPBDBDBDBJBJBDB
=B
=B
=B
=B
=B	7B	7B	7B	7B1B1B1B%BBBBBBBBBBBBBBB  BBB  BB  B  B  B  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�sB�sB�sB�mB�mB�fB�fB�`B�`B�`B�ZB�ZB�ZB�ZB�ZB�ZB�TB�TB�TB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�NB�HB�HB�HB�HB�HB�HB�BB�BB�;B�;B�5B�5B�5B�/B�/B�/B�/B�)B�)B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
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
B�B�
B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B�uB�{B�bB�\B�\B�PB�JB�JB�7B�1B�=B�B}�B� B�B}�B� Bz�Bw�Bv�Bw�Bu�Bu�Bs�Br�Br�Bq�Bo�Bn�Bn�Bn�Bm�Bl�BjBjBk�BhsBffBhsBhsBhsBffBe`BffBffBcTBdZBcTBcTBbNBbNBcTBcTBbNBbNB`BBaHBbNB`BBaHBcTBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BB`BB`BBaHBaHBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B`BB`BB`BB_;B_;B_;B^5B_;BaHB`BB`BBaHB`BB_;B_;B_;B_;B_;B_;B_;B_;B^5B^5B_;B`BBaHBbNBbNBbNBcTBbNBcTBcTBcTBbNBbNBaHB`BB`BBbNB`BBcTBcTBaHB`BB_;B_;BaHBcTBaHBaHB_;B_;BbNB_;B_;B_;B^5B`BB^5B]/B[#B_;B[#B[#B_;B[#B\)BZBZB[#BZB\)BbNBaHBiyBo�Bm�Bl�Bm�Bl�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bm�Bn�Bn�Bo�Bo�Bm�Bn�Bn�Bn�Bn�Bp�Bp�Br�Br�Bp�Bp�Bq�Bs�Bs�Bt�Bs�Bs�Bs�Bs�Bs�Bs�Br�Bs�Bs�Bs�Bs�Bt�Bu�Bu�Bu�Bu�Bv�Bw�Bw�By�B|�B�B�B�B�%B�+B�1B�7B�7B�=B�=B�=B�=B�JB�DB�JB�PB�PB�JB�PB�JB�JB�DB�JB�PB�DB�DB�PB�DB�7B�=B�1B�JB�JB�=B�=B�=B�7B�=B�=B�PB�7B�=B�DB�7B�7B�1B�=B�1B�=B�7B�JB�DB�=B�+B�7B�1B�DB�DB�PB�PB�VB�VB�VB�bB�hB�oB�oB�oB�{G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       B��B��B��B��B��B��B�:B�bB��B��B��B��B��Bz�Bl�Bd�BbBb;BakBa�Bh5Bh�Bd:Bk�Bq�BxBB�0B��B�B�AB��B�tB��B�FB�zBݪB�B�B��B�`B�5B�B��B�B�B�B��B�+B�B��B�qB��B��B�B�LB�mB�B�B�?B�B�B�1B�B�B�B�hB�XB�B�YB�5B�MB�7B�2BݎB�YB�.B��B�LB��B�2B�TB�0B�JBߍB��B��B�B��B�B�FB�PB�sB�uB�B�BޮBގB߃BߌB�B�*B��B�RB�B��B�B�}B�B�ZB�B�9B�$B��B��B��B��B��B��B��B�/B�~B�{B�DB�tB��B�B�7B��B�RBBJB �B�_B�dB�dB�B�B�B�B�B	4B�BQBB.B�B&BjB�B�BkB�BoB`BWBfB&B�B�BeB"B�B�BuB�B�B3BB�B�B�B�B�B�B�B�B�B�B�B�B�BuB9BHBNB�B!,B)�B+ZB.#B/�B0>B0TB0WB01B0~B0�B0jB0<B0:B/jB/�B.�B.�B.3B/$B.$B..B/+B/^B/�B/xB/GB/�B/�B/�B/�B.�B,�B,�B(LB(&B&zB#7B#�B"B!�B#�B$�B%�B%�B%�B&�B'�B(�B)�B+�B.�B1B2B2B3/B4�B7GB6�B7�B3GB2B0�B0�B0B0nB0%B-$B*uB)B)@B)3B)B)-B)2B+;B-"B-3B-0B,�B,�B1�B2�B2�B1�B0tB0�B06B0rB1bB1kB1�B2^B1SB0�B0
B.AB*�B*�B#�B�BB�BlByB�B�B�B�B�BB�B�B B�B�B�B�BLBOBB�BFB9B#/B%:B)@B*B)�B)�B(�B(tB(�B(�B'sB'nB&�B&=B&)B&MB&.B%%B%EB$�B$:B"�B#�B"�B#B#B#B"4B"SB"�B!GB xB �B^BB7BtB�BB.B�B�BaB�B�B�B�B.B�B�B�B�BB)B�BB�B B:BB�BbB�B�B�B�B�B�B_BhB[BQB�BmB�B�B�B�BUB�B\B^B�BHB�B�B5B�B�B�B
sB
nB
MB
fB
nB	zB	QB	DB	GB�B_B�B�B�B�B�B-BB)B@BTB=BB1BXB'BOB 3BBB jB6B BB ?B 6B mB B�pB�YB�ZB�9B�'B�7B�fB��B�
B�B�$B�#B�&B�uB��B�B�B��B�B��B�B��B��B�B�B�B�%B�B�B�B�JB��B��B�B�B�B��B��B��B�=B�;B�EB�6B��B��B�*B��B��B��B��B��B��B��B��B��B�B�)B�B�<B��B��B�B�B��B��B�B�B�B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B��B��B�B�B�B�B�B�B�B�B�B�,B�B�B�OB�B�B�B�B�B�B�B�B�B�B�B��B�B�vB�B��B�(B�B�B�B�B�B�.B�B��B�B�B��B��B�B�B�vB�ZB�YB�B�4B��B�_B�fB�NB�OB�ZB�fB�BB�XB�[B�vB�|B�jB�B�zB�pB�VB�dB��B��B�$B߹B�VB�
BތB�oB݀B�B�xB�.BܒBۧB�BڊB�HB�lBڄBلB�JB�LB�fBوB�IB�+B�FB�BB�B�B�B�B�B�B�B�B�B�B�'B�B�B�B�B�B�B�B�B�B�B�B�B��B� B�B�B�B�B�B�B�B�(B�B�.B�SB�pB�!B�B�B�B�$B�&B�1B�2B�"B�EB�=B�5B� B� B�DB�tB��B��B��B��B�B�B�B��BԝB�#B�1B�B�>B�B�'B�IB�[B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�CB�B��B��B�B��B��B��B��B��B�B�B�B��B��B��B�B��B� B��B��B��B��B� B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��BϸB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B�uB�{B�bB�\B�\B�PB�JB�JB�7B�1B�=B�B}�B� B�B}�B� Bz�Bw�Bv�Bw�Bu�Bu�Bs�Br�Br�Bq�Bo�Bn�Bn�Bn�Bm�Bl�BjBjBk�BhsBffBhsBhsBhsBffBe`BffBffBcTBdZBcTBcTBbNBbNBcTBcTBbNBbNB`BBaHBbNB`BBaHBcTBaHBaHBaHBaHBaHBaHBaHBaHBaHBaHB`BB`BB`BBaHBaHBaHBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BB`BB_;B_;B_;B`BB`BB`BB_;B_;B_;B^5B_;BaHB`BB`BBaHB`BB_;B_;B_;B_;B_;B_;B_;B_;B^5B^5B_;B`BBaHBbNBbNBbNBcTBbNBcTBcTBcTBbNBbNBaHB`BB`BBbNB`BBcTBcTBaHB`BB_;B_;BaHBcTBaHBaHB_;B_;BbNB_;B_;B_;B^5B`BB^5B]/B[#B_;B[#B[#B_;B[#B\)BZBZB[#BZB\)BbNBaHBiyBo�Bm�Bl�Bm�Bl�Bm�Bm�Bm�Bm�Bn�Bn�Bn�Bn�Bm�Bn�Bn�Bo�Bo�Bm�Bn�Bn�Bn�Bn�Bp�Bp�Br�Br�Bp�Bp�Bq�Bs�Bs�Bt�Bs�Bs�Bs�Bs�Bs�Bs�Br�Bs�Bs�Bs�Bs�Bt�Bu�Bu�Bu�Bu�Bv�Bw�Bw�By�B|�B�B�B�B�%B�+B�1B�7B�7B�=B�=B�=B�=B�JB�DB�JB�PB�PB�JB�PB�JB�JB�DB�JB�PB�DB�DB�PB�DB�7B�=B�1B�JB�JB�=B�=B�=B�7B�=B�=B�PB�7B�=B�DB�7B�7B�1B�=B�1B�=B�7B�JB�DB�=B�+B�7B�1B�DB�DB�PB�PB�VB�VB�VB�bB�hB�oB�oB�oB�{G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                       <#�<#ף<#�C<#�<#׎<#ޫ<$�<#��<#�&<$��<.ț<IJg<J�<>W<2��<(\,<$�-<&�9<'�e<(�<=�<L�t<P�r<(�<'d<-D�<+><843<4��<<Ex<,d}<)�5<'�e<,��<:6�<C#�<.��<(;B<%��<''�<E˔<%s<$�<#ۮ<#ٛ<#�<%�M<)�0<%�!<,��<$�Q<#��<#�<#�c<#ۮ<#�{<#�D<#�<$]h<#�<#�]<#ߜ<#��<$/%<$5w<$�<'J�<$}�<#��<#�{<#��<#�<<#�<#�<$�q<$��<#�<$�k<$}<#�I<#��<#�i<#׺<#�<$=<#�<$�<$:�<$�<#�<#ף<#�^<#��<#�a<$T�<$�<#�5<#��<#�<#��<$�L<$+<#��<$� <#�m<#�i<#�<<$/<#�&<#�0<#��<$�<#�<<#��<#��<#�<#��<#�e<#�r<#�<$"2<$Z<#��<#�<#�$<#�J<#��<$v�<%ȧ<#�<$��<$@|<#�)<#�Q<$F<#�N<#��<#�<#�&<#�"<$��<#�<$��<#��<%�d<%,#<$Sa<#�<#�g<$ K<#�D<#�J<#��<#�X<#�<#��<#�<$H�<#�N<#�o<#�)<#�$<#�{<#�8<#��<$�Q<$+<$
<#ٛ<#ף<#׺<#�8<#�N<#�<#�o<#��<#�<#��<#�l<$
�<$r<#׎<#��<#��<% �<&<%�~<#�<$3U<#׺<$
�<#ٛ<#��<#��<#��<#�m<$.<#�M<#�D<#��<#�<$�<$T�<#��<#��<#�I<#��<#�<#��<#�<#�N<#�m<#��<#��<$<#��<$F<$	<$:�<&Z�<#��<%�<%�b<#��<$��<#�e<#�o<#�D<#�<#׎<#�
<#�I<#ף<#�<#��<#��<#��<#�)<#ا<#��<#ا<#�<$G<$�<$$<(<$�<$� <$�<#�Q<#�&<#�!<$�k<$�w<$�<#��<#�<#�e<#��<#�E<#�<#��<#�C<#�8<#ڑ<#��<#�m<$�<$:�<$a<$W<#�<%��<$ʾ<%<�<#�<#�&<$"2<#�^<#��<#��<$� <$��<$U�<)�<'*�<$k�<#��<$�Q<$f�<$x+<#��<#�Q<#�N<#�J<#�c<$
<#��<#�i<#��<#�U<#�<#�<#�<$7�<$4e<$�<#��<#��<%��<$-<$.<#�<$��<$��<$_�<$�j<$.<$�<$MO<$<$�<$Y�<#��<#��<#�Q<#�4<#��<#�	<$�7<#��<#�c<#�&<#�<#�!<#�<#�e<#��<$r<$F9<$
�<$@|<$U�<$+<#�<$�<$I�<$Ş<#��<$�<#�D<#�<$O�<#�E<#�+<#ܯ<#��<$&<#��<#�<#�<#��<$ <$'<#��<$�<$�<$*<$T�<%��<#�<#�{<#��<#��<#�a<#�E<#�&<#�<#׺<#��<#�i<#�<#ߜ<#ٛ<#�J<#ޫ<#�<#�5<#�<#�e<#�{<#ף<#�<$�Q<#��<#��<#׺<$)
<%$<$:�<#��<#�^<#��<#�+<#�^<#�<#�<#׎<#��<#��<#�<$ �<$O�<$i&<$B�<$	�<#�D<#�
<#ٛ<#�J<#��<#�^<#׎<#ܯ<#�<#�]<#�M<#��<#�i<#��<#�Q<#�<#�U<#�&<#��<#�H<$<<$�<#�<#�(<#�<#��<#�<#��<#�<#�o<#�^<#��<#ߜ<#��<$�<$9�<#��<#��<#ٛ<#��<#ٛ<#��<#׎<#ا<#�+<#�r<#�l<#�<#�+<#ۮ<#��<$ �<#�]<#��<#�^<#�<#�<#��<#ڑ<#�<$ <#��<$
�<#�a<#ܯ<#�D<#��<#�*<#��<#ף<#�<#�<#׎<#�r<#��<#�l<#�!<$�<#�N<$�<#��<#�<#�<<#��<#��<#��<#׺<#�{<#�<#ޫ<#��<#�r<#�+<#ߜ<#�<#��<#׺<#�i<#�
<#��<#�<#�<#��<#�*<#�<#�U<#��<#׺<#�<#��<#�8<#��<#ף<$�<$r<$ K<#�J<$�<$L<#؄<#�<#�J<#��<#��<#�D<#�8<#�D<#�8<#�I<#�^<#��<#�e<#�I<#�o<#�<$4e<$%<#�$<#�M<%��<#�<$P�<#�<$!><#�<#�&<$ K<$/<#�<#�<#�o<#�
<#�<#�4<$o�<$3U<#��<#��<#�
<#�<#�{<#��<#�{<#�X<#׎<#��<#�<#ڑ<#�<#ޫ<#��<#ף<#�o<$ K<$�<$|d<$v<#�]<$aD<#�"<#�<#�<#�<#�M<$�j<#��<$.<$|d<#�H<#ܯ<#�<#�g<#�H<#��<#ߜ<#�<#��<#��<#�D<#��<#ޫ<#�{<#�<#�<#׎<#׎<#�{<#�
<#�<#�<#�<#ا<#�<#�<#ף<#�i<#�<#�{<#�
<#�<#�<#�i<#�
<#�<#�<#��<#�<#�<#�<#�X<#�
<#ף<#׎<#��<#׺<#��<#�M<#��<#ا<#�<#׎<#�<#�<#�o<#ۮ<#��<#��<#��<#��<$�-<#ڑ<#�<#��<$�<#ٛ<#�*<#�<#׺<#�o<#�r<#�{<#�&<$*<#�^<#�&<#�8<#�<#�o<#ߜ<#�m<#��<#�<#�{<#�X<#�<#ף<#�<#�<#׎<#�<#�<#�<<#�{<#�<#�
<#��<#�<#�<#�<#׎<#�
<#�<#ף<#��<#�"<#�+<#׺<#�i<#�c<#�<#��<#�I<#�
<#�<#ߜ<#ޫ<#�^<#�{<#׎<#׎<#�r<#׺<#�D<#��<#׎<#�X<#�i<#�]<#׎<#�
<#ا<#�
<#�
<#��<#�<#�<#�<#�<#�<#׎<#�<#�I<#�{<#�<#�
<#�<#�<#�X<#��<#�<#ۮ<#�<#�i<#׎<#�
<#ף<#׎<#�<#�X<#׎<#�<#�{<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                