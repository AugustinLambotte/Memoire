CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  \   N_CALIB       	N_HISTORY             
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
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  O�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  TD   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g$   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  k�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ~�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �T   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Ӭ   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �d   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �D   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   <   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   	<   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   <   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T <   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   0   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   L   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    T   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        t   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        |   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �Argo profile    3.1 1.2 19500101000000  20181106012856  20230209090617  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               N   NAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�t&F��B@�t&F��B11  @�t&O�`@�t&O�`@N��GG�@N��GG��;�	� ��;�	� �11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�z�@   @@  @�  @�G�@�G�@�  A   A��A#�
A?\)A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B(�B(�B   B'�
B0(�B8(�B?�
BG33BO�
BW�B`  Bh  Bo�Bw�B�
B�(�B�{B��B�{B�  B��
B��
B��
B��
B��B��B�  B�(�B�  B��
B��B�  B�(�B�{B��B�  B��B�  B�{B�  B��B��
B��
B��
B�  B�(�C {C  C  C��C  C
{C  C��C��C��C�C��C  C
=C
=C
=C {C"  C#�C&  C(
=C*  C+�C.  C0
=C2  C3��C6
=C8�C:  C;��C>
=C?��CA��CD  CE��CG��CI��CL  CN  CP  CR
=CS�CU��CX  CY��C[�C^  C`  Ca�Cc�Ce�Cg��Cj
=Cl  Cm�Co��Cr  Ct{Cv{Cx{Cz
=C|  C~
=C�
=C�C�
=C�  C�  C���C���C�C�  C�C�C���C���C�C�C�C�C���C���C���C�  C�
=C�  C���C���C�  C���C���C���C�  C���C�  C�C�  C���C���C�  C�C�C�C�C�C�C�  C���C���C���C�  C�C���C���C�  C�  C�  C�C�  C���C�C�  C���C���C�  C�  C�C�
=C�
=C�
=C�
=C�  C��C���C�
=C�
=C�C�C�C�C�  C���C���C���C�  C�  C�
=C�
=C�  C���C���C�C�
=C�C�  C�C�C�C�  C���C�C�  C�  C�  C���C�  C�  C���C�C�C�  C���C���C�  C�C�  C���C���C���C���C���C���C���C���C�  C�
=C�C�  C���C��C���D   D � D�D}qD��D��DD� D��D}qD�D� D  D� D�D��DD�D�qD	}qD
  D
}qD  D��D�qD� D  D� D�D� D�D�D�D}qD��Dz�D�RD}qD  D}qD��D� D  Dz�D  D��D�D}qD��Dz�D  D�D  Dz�D  D��DD��DD�D  Dz�D�RDxRD�RD � D!�D!��D"�D"�D#D#� D#��D$xRD$�RD%xRD%�RD&� D'D'��D(  D(� D(�qD)z�D*  D*��D+D+� D,  D,��D,��D-� D.�D.��D/D/� D/��D0}qD0�qD1}qD1��D2xRD2��D3z�D3�qD4� D5�D5��D6�D6}qD7  D7}qD7��D8� D8�qD9z�D9��D:xRD:�qD;��D<D<� D<�qD=}qD=��D>��D?�D?� D@�D@�DA�DA��DB  DB�DC�DC� DC�qDD}qDD�qDEz�DE��DF��DG  DG}qDH�DH� DH�qDI� DJDJ��DK�DK� DK�qDL�DM�DM}qDN�DN� DN��DO}qDP�DP�DQ�DQ}qDR  DR��DS  DS� DT�DT�DU  DUz�DU�qDV��DW  DW}qDX�DX�DY�DY� DZDZ��D[�D[� D[�qD\�D]D]z�D]��D^� D_D_��D_�qD`}qD`�qDa}qDb�Db� Db�qDc}qDc��Dd}qDe�De�Df�Df�Df��Dgz�Dh  Dh� Di  Di��Dj�Dj� Dk�Dk��Dl�Dl��Dl�qDm}qDn  Dn}qDo  Do� Do�qDpz�Dp�qDq� Dr  Dr��Ds  Ds}qDt  Dt��Du�Du� Du�qDvz�Dv�qDw� Dw�qDx��Dy�Dy}qDz  Dz��D{  D{}qD{��D|� D}�D}��D~  D~}qD�D� D�HD�AHD��HD�D��D�@ D�}qD�� D�HD�AHD�� D��HD��D�@ D�� D�D�HD�@ D�~�D�� D��D�@ D�� D�� D��qD�=qD�}qD��qD���D�>�D�~�D���D���D�@ D�� D���D��D�AHD�� D�� D��qD�>�D��HD���D�  D�@ D�~�D���D��qD�=qD�~�D�� D���D�>�D�~�D�� D�  D�@ D�~�D�� D�HD�>�D�}qD���D���D�>�D�� D�� D�HD�@ D�~�D�� D�HD�AHD���D��HD��D�@ D�~�D�� D�  D�>�D�~�D��HD�  D�>�D�� D�˅?k�?u?�\)?��R?���?�p�?Ǯ?�
=?�G�?��@   @�@�@��@�R@&ff@.{@333@:�H@B�\@J=q@Tz�@\(�@aG�@h��@n{@xQ�@�  @��
@��@��@�\)@��@�@���@�(�@��\@��@���@��@�\)@�z�@�Q�@�(�@��R@�G�@�ff@�=q@�\)@�33@�
=@ٙ�@�p�@�G�@��@���@�{@��@�z�@�Q�@�(�A   A�A�
AA�A	��A�Ap�A\)AG�A33A�A
=A��A�HA��A�RA ��A"�\A$z�A&ffA(��A*�HA,��A.�RA0��A2�\A3�
A5A8Q�A:=qA;�A=p�A@  AA�AC33AE�AG�AI��AJ�HAL��AN�RAQG�AR�\ATz�AVffAX��AZ=qA\(�A^{A`��Aa�Ac33AeAhQ�Ai��Ak�Amp�Ao\)AqG�As33Au�Aw
=Ax��Az�HA|��A~�RA�Q�A�G�A��\A�33A�(�A��A�ffA�
=A�  A�G�A�=qA��HA��
A��A�{A��RA��A���A���A��HA��
A���A�p�A�ffA��A���A�G�A�=qA��A�z�A��A�{A��A�Q�A���A��A�33A�(�A���A�A�
=A��A���A���A��\A��A�(�A�p�A�ffA�
=A�  A�G�A�=qA��HA��
A��A�{A��RA��A���A��A�=qA��A�z�A�p�A�{A�\)A�Q�A���A��A�33A��
A���A�AƸRA�\)A�Q�Aə�A�=qA�33A�z�A�p�A�{A�
=A�  A���A��A��HA��
A�z�A�p�AָRA�\)A�Q�A�G�A�=qA�33A��
A��A�{A޸RA�  A���AᙚA�\A��
A���A��A�ffA�\)A�Q�A�G�A�=qA��HA��
A���A�A�ffA�A��A�G�A�=qA�A�(�A���A�{A�
=A�  A���A��A��\A��A�z�A�p�A�{A�
=B   B z�B ��Bp�BB=qB�RB33B�B(�B��B�Bp�B�BffB�HB\)B�
B(�B��B	G�B	B
{B
�\B33B\)B  Bz�B��BG�B�BffB�RB33B�B(�B��B�B��B�B�\B
=B\)B�
BQ�B��BG�BB=qB�\B33B�B  Bz�B�B��B�B�\B
=B\)B�
Bz�B��Bp�B�BffB�HB33B�
B Q�B ��B!�B!B"=qB"�\B#33B#�B$  B$z�B%�B%p�B%�B&�\B'
=B'\)B'�
B(z�B(��B)G�B)�B*=qB*�RB+33B+�B,(�B,��B-G�B-��B.{B.�RB/33B/�B0  B0��B0��B1��B2{B2�\B3
=B3�B4  B4z�B5�B5p�B5�B6�\B7
=B7\)B7�
B8z�B8��B9G�B9�B:ffB:�RB;\)B;�
B<Q�B<��B=G�B=B>=qB>�RB?33B?�B@  B@��BA�BA��BB{BB�RBC
=BC�BD(�BD��BD��BEp�BF{BF�\BF�HBG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             ?�z�@   @@  @�  @�G�@�G�@�  A   A��A#�
A?\)A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B(�B(�B   B'�
B0(�B8(�B?�
BG33BO�
BW�B`  Bh  Bo�Bw�B�
B�(�B�{B��B�{B�  B��
B��
B��
B��
B��B��B�  B�(�B�  B��
B��B�  B�(�B�{B��B�  B��B�  B�{B�  B��B��
B��
B��
B�  B�(�C {C  C  C��C  C
{C  C��C��C��C�C��C  C
=C
=C
=C {C"  C#�C&  C(
=C*  C+�C.  C0
=C2  C3��C6
=C8�C:  C;��C>
=C?��CA��CD  CE��CG��CI��CL  CN  CP  CR
=CS�CU��CX  CY��C[�C^  C`  Ca�Cc�Ce�Cg��Cj
=Cl  Cm�Co��Cr  Ct{Cv{Cx{Cz
=C|  C~
=C�
=C�C�
=C�  C�  C���C���C�C�  C�C�C���C���C�C�C�C�C���C���C���C�  C�
=C�  C���C���C�  C���C���C���C�  C���C�  C�C�  C���C���C�  C�C�C�C�C�C�C�  C���C���C���C�  C�C���C���C�  C�  C�  C�C�  C���C�C�  C���C���C�  C�  C�C�
=C�
=C�
=C�
=C�  C��C���C�
=C�
=C�C�C�C�C�  C���C���C���C�  C�  C�
=C�
=C�  C���C���C�C�
=C�C�  C�C�C�C�  C���C�C�  C�  C�  C���C�  C�  C���C�C�C�  C���C���C�  C�C�  C���C���C���C���C���C���C���C���C�  C�
=C�C�  C���C��C���D   D � D�D}qD��D��DD� D��D}qD�D� D  D� D�D��DD�D�qD	}qD
  D
}qD  D��D�qD� D  D� D�D� D�D�D�D}qD��Dz�D�RD}qD  D}qD��D� D  Dz�D  D��D�D}qD��Dz�D  D�D  Dz�D  D��DD��DD�D  Dz�D�RDxRD�RD � D!�D!��D"�D"�D#D#� D#��D$xRD$�RD%xRD%�RD&� D'D'��D(  D(� D(�qD)z�D*  D*��D+D+� D,  D,��D,��D-� D.�D.��D/D/� D/��D0}qD0�qD1}qD1��D2xRD2��D3z�D3�qD4� D5�D5��D6�D6}qD7  D7}qD7��D8� D8�qD9z�D9��D:xRD:�qD;��D<D<� D<�qD=}qD=��D>��D?�D?� D@�D@�DA�DA��DB  DB�DC�DC� DC�qDD}qDD�qDEz�DE��DF��DG  DG}qDH�DH� DH�qDI� DJDJ��DK�DK� DK�qDL�DM�DM}qDN�DN� DN��DO}qDP�DP�DQ�DQ}qDR  DR��DS  DS� DT�DT�DU  DUz�DU�qDV��DW  DW}qDX�DX�DY�DY� DZDZ��D[�D[� D[�qD\�D]D]z�D]��D^� D_D_��D_�qD`}qD`�qDa}qDb�Db� Db�qDc}qDc��Dd}qDe�De�Df�Df�Df��Dgz�Dh  Dh� Di  Di��Dj�Dj� Dk�Dk��Dl�Dl��Dl�qDm}qDn  Dn}qDo  Do� Do�qDpz�Dp�qDq� Dr  Dr��Ds  Ds}qDt  Dt��Du�Du� Du�qDvz�Dv�qDw� Dw�qDx��Dy�Dy}qDz  Dz��D{  D{}qD{��D|� D}�D}��D~  D~}qD�D� D�HD�AHD��HD�D��D�@ D�}qD�� D�HD�AHD�� D��HD��D�@ D�� D�D�HD�@ D�~�D�� D��D�@ D�� D�� D��qD�=qD�}qD��qD���D�>�D�~�D���D���D�@ D�� D���D��D�AHD�� D�� D��qD�>�D��HD���D�  D�@ D�~�D���D��qD�=qD�~�D�� D���D�>�D�~�D�� D�  D�@ D�~�D�� D�HD�>�D�}qD���D���D�>�D�� D�� D�HD�@ D�~�D�� D�HD�AHD���D��HD��D�@ D�~�D�� D�  D�>�D�~�D��HD�  D�>�D�� D�˅?k�?u?�\)?��R?���?�p�?Ǯ?�
=?�G�?��@   @�@�@��@�R@&ff@.{@333@:�H@B�\@J=q@Tz�@\(�@aG�@h��@n{@xQ�@�  @��
@��@��@�\)@��@�@���@�(�@��\@��@���@��@�\)@�z�@�Q�@�(�@��R@�G�@�ff@�=q@�\)@�33@�
=@ٙ�@�p�@�G�@��@���@�{@��@�z�@�Q�@�(�A   A�A�
AA�A	��A�Ap�A\)AG�A33A�A
=A��A�HA��A�RA ��A"�\A$z�A&ffA(��A*�HA,��A.�RA0��A2�\A3�
A5A8Q�A:=qA;�A=p�A@  AA�AC33AE�AG�AI��AJ�HAL��AN�RAQG�AR�\ATz�AVffAX��AZ=qA\(�A^{A`��Aa�Ac33AeAhQ�Ai��Ak�Amp�Ao\)AqG�As33Au�Aw
=Ax��Az�HA|��A~�RA�Q�A�G�A��\A�33A�(�A��A�ffA�
=A�  A�G�A�=qA��HA��
A��A�{A��RA��A���A���A��HA��
A���A�p�A�ffA��A���A�G�A�=qA��A�z�A��A�{A��A�Q�A���A��A�33A�(�A���A�A�
=A��A���A���A��\A��A�(�A�p�A�ffA�
=A�  A�G�A�=qA��HA��
A��A�{A��RA��A���A��A�=qA��A�z�A�p�A�{A�\)A�Q�A���A��A�33A��
A���A�AƸRA�\)A�Q�Aə�A�=qA�33A�z�A�p�A�{A�
=A�  A���A��A��HA��
A�z�A�p�AָRA�\)A�Q�A�G�A�=qA�33A��
A��A�{A޸RA�  A���AᙚA�\A��
A���A��A�ffA�\)A�Q�A�G�A�=qA��HA��
A���A�A�ffA�A��A�G�A�=qA�A�(�A���A�{A�
=A�  A���A��A��\A��A�z�A�p�A�{A�
=B   B z�B ��Bp�BB=qB�RB33B�B(�B��B�Bp�B�BffB�HB\)B�
B(�B��B	G�B	B
{B
�\B33B\)B  Bz�B��BG�B�BffB�RB33B�B(�B��B�B��B�B�\B
=B\)B�
BQ�B��BG�BB=qB�\B33B�B  Bz�B�B��B�B�\B
=B\)B�
Bz�B��Bp�B�BffB�HB33B�
B Q�B ��B!�B!B"=qB"�\B#33B#�B$  B$z�B%�B%p�B%�B&�\B'
=B'\)B'�
B(z�B(��B)G�B)�B*=qB*�RB+33B+�B,(�B,��B-G�B-��B.{B.�RB/33B/�B0  B0��B0��B1��B2{B2�\B3
=B3�B4  B4z�B5�B5p�B5�B6�\B7
=B7\)B7�
B8z�B8��B9G�B9�B:ffB:�RB;\)B;�
B<Q�B<��B=G�B=B>=qB>�RB?33B?�B@  B@��BA�BA��BB{BB�RBC
=BC�BD(�BD��BD��BEp�BF{BF�\BF�HBG�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�AoAffAdZA+AVA��A�yA�AȴA�9A�\A�Av�An�AZA1A��A�jAn�A$�A9XA �AXA~�AjA`BA�A
��A	�A	��A	�A	Av�A��AQ�A9XA��AV@���@�~�@�  @��H@�v�@���@��j@�9X@�;d@�+@�M�@��/@�@�9@��/@��@��`@�Ĝ@��@��@��@���@��`@��@�@�`B@���@�J@��@�$�@�-@�=q@�E�@�E�@��@��T@�-@�7@�X@�Ĝ@�j@�A�@� �@�b@��/@���@���@���@��/@��@�%@���@���@��`@��/@��D@�@�A�@� �@���@@@�;d@�"�@��@���@�"�@�C�@�"�@�33@�;d@�;d@�\)@�@�\@�-@���@�?}@�%@��/@��@�G�@�7L@�7L@�V@���@��`@���@�%@�j@�u@�r�@�Z@�bN@�Q�@�I�@�Q�@�Q�@���@�S�@��@�@�@���@��@�ȴ@�!@�v�@�M�@�5?@��@��@�{@�J@�J@�{@�{@���@�^@��@�G�@��/@�@蛦@�u@�u@�@�bN@�1'@畁@�33@���@�ff@��@�@噚@噚@噚@噚@�7@�?}@���@���@���@���@���@��`@���@�Ĝ@�u@�I�@� �@�b@�1@��m@��;@��
@㝲@㕁@�t�@�dZ@�dZ@�dZ@�\)@�K�@�"�@���@�-@��@�@���@���@�@�@���@��@�^@���@�p�@��@���@��/@���@���@�j@�@���@���@��u@�I�@�
=@޸R@ޟ�@އ+@�M�@�@݁@�Ĝ@ܛ�@�Z@�1'@� �@�  @�1@�ƨ@�o@�1@�@ݺ^@�p�@�G�@ܴ9@�I�@�(�@��
@۝�@�S�@��@ڗ�@�V@��T@ٙ�@�x�@�X@�%@�1@��;@׶F@׍P@�t�@�|�@ם�@׮@�|�@�
=@��H@ָR@�n�@�5?@��#@�@ՙ�@�p�@�`B@�/@���@Լj@Դ9@ԃ@�r�@Ӯ@�ȴ@�v�@�=q@��T@�/@���@��`@мj@Гu@�r�@�j@�j@�Z@�1'@�|�@�^5@��@��@��T@���@�x�@�G�@��@�1'@�C�@��@��@ȓu@�Q�@��@��
@ǍP@�S�@�C�@�o@��@Ɨ�@�V@�5?@�p�@�I�@öF@Õ�@Õ�@Å@�t�@�\)@��@¸R@�-@���@���@�p�@��@� �@���@�+@�ȴ@��\@�^5@�-@�J@�@���@�@�X@��/@�Q�@�  @��
@���@�dZ@���@���@�S�@�"�@��y@���@�7L@���@���@��@�I�@�  @��@�l�@�C�@�33@�
=@�ȴ@�M�@�$�@���@���@��u@�t�@���@���@�V@�@��@��u@��@���@�-@���@��@��#@���@���@��^@���@���@�G�@�&�@��/@���@��9@��@��P@��R@��!@���@�~�@��T@�V@��`@��9@��u@��@�r�@�Z@�I�@�A�@�A�@� �@�  @��
@���@�"�@�@���@�ȴ@��R@��!@���@�~�@�n�@�$�@��7@��@�V@���@��/@�Ĝ@��@��D@��@�j@�bN@�Q�@��@���@���@��w@���@�l�@�33@��@��@�
=@�@���@��@��@��@��y@��@��@��@��@��@��y@��H@��y@��y@��y@��y@��H@��H@���@���@��@���@��!@��\@��+@�~�@�v�@�n�@�n�@�^5@�^5@�ff@�ff@�^5@�E�@�5?@�$�@�{@�{@�{@�{@�J@�{@���@��#@��-@���@���@���@���@��h@��7@��7@��7@��7@�x�@�hs@�O�@�?}@�&�@��@���@��@���@��9@��u@�r�@� �@��;@��w@��@�+@�o@�@��@��H@���@��\@��\@��+@��+@�~�@�v�@�^5@�=q@�$�@���@�O�@�G�@�/@��@�%@���@��@���@��D@��D@��@�bN@�I�@�9X@�1@���@�ƨ@��w@���@��P@��@�|�@�l�@�;d@�"�@�
=@���@��@���@��+@�ff@�n�@�n�@�n�@�n�@�ff@�ff@�^5@�V@�M�@�J@�O�@��`@��/@��@�A�@�(�@��m@��;@���@��w@��F@�t�@�;d@�@��yA+A+A"�AA��A��A�A�yA��A�9A��A1'A�A��A��A�PAt�AdZA`BA`BA\)AXAO�A;dA+A+A+A&�A&�A&�A"�A"�A�A�A�AVA
=A%A%AA%AA��AA��A��A��A��A�A�A�yA�A�yA�yA�yA�`A�HA�/A�A�A�A�A�A��A�A��A��A��AȴAȴAȴAȴAĜA��A��A�jA�jA�jA�jA�RA�RA�!A�A��A��A��A��A��A��A�uA��A��A�uA�DA�DA�DA�DA�+A�+A�+A�+A�+A�+A�+A�+A�+A�DA�+A�A�A�A�+A�A�+A�A�A�A�A�A~�Az�Av�A~�A~�Az�Av�Az�Az�Av�Av�Av�Av�Av�Ar�Ar�Ar�Ar�Av�An�Ar�Ar�An�Ar�An�An�An�An�An�An�An�AjAjAjAffA^5AZA^5AffAffAffA^5A^5AZAZA^5AVAM�AQ�AM�AI�AM�AM�AI�A9XA1'A-A1'A$�A�AbA��A�mA�#AƨA�-A��A��A�7AhsA?}A�A%A��A�A�yA�`A�HA�HA�/A�A�A�A��A��A��A��A��A��A��AȴAȴAȴAĜA��A�jA�9A�A��A��A��A��A��A��A�uA�\A�+A�A~�Av�An�AjA^5AZAVAM�AI�AE�A=qA5?A-A(�A$�A$�A �A �A�A �A �A$�A �A �A$�A(�A(�A-A1'A5?A5?A5?A5?A5?A9XA9XA=qA9XA5?A-A1'A5?A9XA=qA=qAA�AA�AA�AA�A=qA=qA9XA9XA5?A(�A �A�A�AJA1A1A  A�mA��A�wA�-A��A��A��A�hA�PAp�A?}A/A/A+A&�A"�A�A��A��A��A�A�!A�A��A��A��A�\A�An�AbNAbNAZAZAQ�AE�A=qA{A��A�#AA�-A��A�A`BA33A��A  A�A$�AJA�A�A�wA�wAƨA��A�A\)AXA\)AXA\)A`BA\)AS�A�A�A%A�/A�+Az�A��AQ�AJA��AhsA��A~�A=qA1'A�A�#A��A�^A�-A�A��A��A�hA�AdZA"�A
�A
�yA
��A
��A
ȴA
�!A
~�A
r�A
VA
-A
 �A
�A
{A
JA
A
A	��A	�A	�TA	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	A	ƨA	�#A	�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             AoAffAdZA+AVA��A�yA�AȴA�9A�\A�Av�An�AZA1A��A�jAn�A$�A9XA �AXA~�AjA`BA�A
��A	�A	��A	�A	Av�A��AQ�A9XA��AV@���@�~�@�  @��H@�v�@���@��j@�9X@�;d@�+@�M�@��/@�@�9@��/@��@��`@�Ĝ@��@��@��@���@��`@��@�@�`B@���@�J@��@�$�@�-@�=q@�E�@�E�@��@��T@�-@�7@�X@�Ĝ@�j@�A�@� �@�b@��/@���@���@���@��/@��@�%@���@���@��`@��/@��D@�@�A�@� �@���@@@�;d@�"�@��@���@�"�@�C�@�"�@�33@�;d@�;d@�\)@�@�\@�-@���@�?}@�%@��/@��@�G�@�7L@�7L@�V@���@��`@���@�%@�j@�u@�r�@�Z@�bN@�Q�@�I�@�Q�@�Q�@���@�S�@��@�@�@���@��@�ȴ@�!@�v�@�M�@�5?@��@��@�{@�J@�J@�{@�{@���@�^@��@�G�@��/@�@蛦@�u@�u@�@�bN@�1'@畁@�33@���@�ff@��@�@噚@噚@噚@噚@�7@�?}@���@���@���@���@���@��`@���@�Ĝ@�u@�I�@� �@�b@�1@��m@��;@��
@㝲@㕁@�t�@�dZ@�dZ@�dZ@�\)@�K�@�"�@���@�-@��@�@���@���@�@�@���@��@�^@���@�p�@��@���@��/@���@���@�j@�@���@���@��u@�I�@�
=@޸R@ޟ�@އ+@�M�@�@݁@�Ĝ@ܛ�@�Z@�1'@� �@�  @�1@�ƨ@�o@�1@�@ݺ^@�p�@�G�@ܴ9@�I�@�(�@��
@۝�@�S�@��@ڗ�@�V@��T@ٙ�@�x�@�X@�%@�1@��;@׶F@׍P@�t�@�|�@ם�@׮@�|�@�
=@��H@ָR@�n�@�5?@��#@�@ՙ�@�p�@�`B@�/@���@Լj@Դ9@ԃ@�r�@Ӯ@�ȴ@�v�@�=q@��T@�/@���@��`@мj@Гu@�r�@�j@�j@�Z@�1'@�|�@�^5@��@��@��T@���@�x�@�G�@��@�1'@�C�@��@��@ȓu@�Q�@��@��
@ǍP@�S�@�C�@�o@��@Ɨ�@�V@�5?@�p�@�I�@öF@Õ�@Õ�@Å@�t�@�\)@��@¸R@�-@���@���@�p�@��@� �@���@�+@�ȴ@��\@�^5@�-@�J@�@���@�@�X@��/@�Q�@�  @��
@���@�dZ@���@���@�S�@�"�@��y@���@�7L@���@���@��@�I�@�  @��@�l�@�C�@�33@�
=@�ȴ@�M�@�$�@���@���@��u@�t�@���@���@�V@�@��@��u@��@���@�-@���@��@��#@���@���@��^@���@���@�G�@�&�@��/@���@��9@��@��P@��R@��!@���@�~�@��T@�V@��`@��9@��u@��@�r�@�Z@�I�@�A�@�A�@� �@�  @��
@���@�"�@�@���@�ȴ@��R@��!@���@�~�@�n�@�$�@��7@��@�V@���@��/@�Ĝ@��@��D@��@�j@�bN@�Q�@��@���@���@��w@���@�l�@�33@��@��@�
=@�@���@��@��@��@��y@��@��@��@��@��@��y@��H@��y@��y@��y@��y@��H@��H@���@���@��@���@��!@��\@��+@�~�@�v�@�n�@�n�@�^5@�^5@�ff@�ff@�^5@�E�@�5?@�$�@�{@�{@�{@�{@�J@�{@���@��#@��-@���@���@���@���@��h@��7@��7@��7@��7@�x�@�hs@�O�@�?}@�&�@��@���@��@���@��9@��u@�r�@� �@��;@��w@��@�+@�o@�@��@��H@���@��\@��\@��+@��+@�~�@�v�@�^5@�=q@�$�@���@�O�@�G�@�/@��@�%@���@��@���@��D@��D@��@�bN@�I�@�9X@�1@���@�ƨ@��w@���@��P@��@�|�@�l�@�;d@�"�@�
=@���@��@���@��+@�ff@�n�@�n�@�n�@�n�@�ff@�ff@�^5@�V@�M�@�J@�O�@��`@��/@��@�A�@�(�@��m@��;@���@��w@��F@�t�@�;d@�@��yA+A+A"�AA��A��A�A�yA��A�9A��A1'A�A��A��A�PAt�AdZA`BA`BA\)AXAO�A;dA+A+A+A&�A&�A&�A"�A"�A�A�A�AVA
=A%A%AA%AA��AA��A��A��A��A�A�A�yA�A�yA�yA�yA�`A�HA�/A�A�A�A�A�A��A�A��A��A��AȴAȴAȴAȴAĜA��A��A�jA�jA�jA�jA�RA�RA�!A�A��A��A��A��A��A��A�uA��A��A�uA�DA�DA�DA�DA�+A�+A�+A�+A�+A�+A�+A�+A�+A�DA�+A�A�A�A�+A�A�+A�A�A�A�A�A~�Az�Av�A~�A~�Az�Av�Az�Az�Av�Av�Av�Av�Av�Ar�Ar�Ar�Ar�Av�An�Ar�Ar�An�Ar�An�An�An�An�An�An�An�AjAjAjAffA^5AZA^5AffAffAffA^5A^5AZAZA^5AVAM�AQ�AM�AI�AM�AM�AI�A9XA1'A-A1'A$�A�AbA��A�mA�#AƨA�-A��A��A�7AhsA?}A�A%A��A�A�yA�`A�HA�HA�/A�A�A�A��A��A��A��A��A��A��AȴAȴAȴAĜA��A�jA�9A�A��A��A��A��A��A��A�uA�\A�+A�A~�Av�An�AjA^5AZAVAM�AI�AE�A=qA5?A-A(�A$�A$�A �A �A�A �A �A$�A �A �A$�A(�A(�A-A1'A5?A5?A5?A5?A5?A9XA9XA=qA9XA5?A-A1'A5?A9XA=qA=qAA�AA�AA�AA�A=qA=qA9XA9XA5?A(�A �A�A�AJA1A1A  A�mA��A�wA�-A��A��A��A�hA�PAp�A?}A/A/A+A&�A"�A�A��A��A��A�A�!A�A��A��A��A�\A�An�AbNAbNAZAZAQ�AE�A=qA{A��A�#AA�-A��A�A`BA33A��A  A�A$�AJA�A�A�wA�wAƨA��A�A\)AXA\)AXA\)A`BA\)AS�A�A�A%A�/A�+Az�A��AQ�AJA��AhsA��A~�A=qA1'A�A�#A��A�^A�-A�A��A��A�hA�AdZA"�A
�A
�yA
��A
��A
ȴA
�!A
~�A
r�A
VA
-A
 �A
�A
{A
JA
A
A	��A	�A	�TA	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	��A	A	ƨA	�#A	�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�BǮBȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�BB�mB�B�B��B��B��B��B��B��B��B��B�B��B��B��B+BDBJB\BVBPB\BVBDBDBJBJBbB�B�B�B�B�B�B�B�B$�B+B-B0!B5?B<jBA�BC�BD�BE�BE�BF�BF�BH�BI�BI�BI�BH�BH�BG�BF�BF�BF�BF�BL�BN�BN�BN�BN�BN�BN�BN�BN�BN�BO�BN�BP�BP�BP�BO�BN�BN�BK�BK�BK�BJ�BM�BN�BN�BN�BN�BN�BO�BN�BM�BK�BH�BF�BE�BF�BG�BI�BI�BJ�BI�BI�BI�BI�BJ�BJ�BI�BI�BI�BH�BH�BH�BH�BG�BF�BE�BE�BD�BD�BD�BC�BB�BB�BA�BA�BB�BB�BB�BB�BB�BB�BB�BB�BB�BA�BA�B@�B@�B@�B@�B@�B?}B?}B?}B>wB=qB<jB;dB9XB7LB6FB6FB7LB7LB8RB8RB8RB7LB7LB7LB8RB8RB8RB8RB8RB7LB7LB6FB6FB6FB6FB5?B5?B5?B49B49B5?B5?B5?B5?B49B33B2-B0!B1'B0!B0!B1'B1'B1'B1'B1'B1'B0!B0!B.B.B.B.B.B.B.B.B-B-B,B(�B'�B'�B&�B%�B#�B"�B�B�B�B�B�B�B�B�B�B$�B2-B5?B5?B6FB49B49B33B49B49B33B33B33B2-B33B33B49B49B2-B0!B/B/B/B0!B2-B5?B7LB7LB8RB8RB8RB8RB8RB8RB8RB8RB7LB7LB7LB7LB7LB7LB7LB6FB6FB6FB6FB6FB6FB7LB7LB7LB7LB7LB7LB7LB7LB7LB6FB6FB6FB6FB6FB6FB5?B5?B5?B49B33B2-B0!B+B)�B(�B'�B'�B'�B(�B(�B)�B+B+B+B+B+B,B,B,B,B,B,B,B,B,B+B+B)�B)�B)�B(�B(�B(�B'�B'�B&�B&�B&�B%�B%�B$�B#�B"�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B{B{B{B{B{BuBuBoBuBuBoBoBhB\BPBDB	7B+B+B+B+B+B%B%B%B%BB%BBBBBBBBBB  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�yB�sB�sB�sB�sBǮBǮBǮBǮBǮBƨBƨBƨBǮBŢBɺB��BȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�
B�B�B�B�B�)B�/B�)B�/B�5B�5B�;B�;B�;B�HB�HB�HB�HB�NB�HB�HB�HB�ZB�ZB�ZB�`B�`B�`B�fB�`B�`B�sB�yB�sB�sB�sB�sB�mB�yB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�B�B�B�B�yB�yB�yB�yB�B�B  BB��B��B��BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�B��BB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             B̌B͠B�B�pB�SB�?B�DB�FB�MB�fB�B�B�B�-B��BϾB͢B̱B͡BѸB�qB�wB�B��B��B�B��B��B�_B��B�oB��B�^B��B�wB��BjB
�BtB�BBB�B B�B�ByB�BXB�BaBaB^B�B�B�B�B�B$mB+B-�B/�B4/B;�BA)BC}BD�BE�BE�BF�BF�BH�BJ
BJBI�BIBI�BH8BF�BF�BF�BE{BL�BOBN�BN�BN�BN�BN�BOBN�BN�BPNBN�BQ\BQBQbBP,BN�BOUBK�BK�BL@BJDBM�BOBN�BN�BN�BN�BP`BO�BNqBL�BIQBGBE�BFTBGjBI�BI�BJ�BI�BI�BI�BI�BK,BJ�BI�BI�BI�BH�BH�BH�BH�BHtBGbBE�BE�BD�BD�BD�BC�BB�BB�BA�BA�BB�BB�BB�BB�BB�BB�BB�BB�BB�BA�BBBA%B@�B@�B@�B@�B?�B?�B?�B?dB>B<�B<AB:B7�B6�B6JB7IB7PB8oB8�B8�B7QB7LB7PB8QB8xB8mB8oB8�B7�B7�B6bB6UB6wB6SB5RB5�B5MB4pB4WB5?B5?B5LB5VB4|B3�B3B08B1HB00B0"B1B1(B13B1CB1lB1B0�B0�B.LB.FB.-B.B..B.0B.,B.B-"B-�B-�B)vB(B(B'@B&YB$�B#�B
BB�B�B�B�BB�B6B"ZB2 B5�B5�B7B4�B4oB3�B4�B4�B3�B3�B3�B2�B3�B3kB4nB4�B3�B0pB/XB/ZB/EB0B1�B5+B7�B7�B8�B8�B8�B8�B8�B8}B8�B8�B7hB7�B7�B7tB7ZB7�B7oB7qB7�B6�B6�B6�B7\B7�B7qB7�B7�B7{B7XB7OB7eB7�B7lB7�B6�B6�B6dB6hB5�B5�B5�B5�B4�B2�B2�B+�B*gB)OB(ZB(hB(GB)B)FB*VB+`B+`B+FB,5B,�B,�B,=B,B,"B,$B,3B,pB,�B,�B+XB+RB*�B++B*�B)�B)�B)�B(RB(;B'4B'B&�B%�B&=B%�B$�B#�B"HB!B B B�B�BB�B�B2BoB%B�B�B�B�B�B�B�B�B�B�B,B�B�B�B"B1B.B�B�BJBUB*B`B;B�BsB9BMB8B6B?B@B4B�B]B�B4B_BB�BSB"B!B@B �B�=B�2B�6B�B�B�B�B�B��B��B�B�B�-B�OB��B�B��B�-B��B��B��B�B��B�TB��B�xB��B��B��B��B��B��B��B��B��B��B�B��B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�:B�'B��B�B�@B��B��B��B��B�B��B�B�B�B�B�B��B��B��B�rB�$B�B��B�B��B�B�B��B�B�B�B��B��B�B��B��B�B�B��B�B�B�B�B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B��B�B�$B�B��B�%B�B��B�B�B�B�B��B��B��B�B�vBǮBǮBǮBǮBǮBƨBƨBƨBǮBŢBɺB��BȴB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�
B�B�B�B�B�)B�/B�)B�/B�5B�5B�;B�;B�;B�HB�HB�HB�HB�NB�HB�HB�HB�ZB�ZB�ZB�`B�`B�`B�fB�`B�`B�sB�yB�sB�sB�sB�sB�mB�yB�yB�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�B�B�B�B�yB�yB�yB�yB�B�B  BB��B��B��BB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�B��BB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                             <5g�<5�{<$�<$�<$p<#�(<#�)<#�<#�<$p<#��<#��<#��<#��<$�J<*w<$aD<$}�<$_�<#�J<#�H<'�.<)k�<Au�<[��<E˔<<e<+'�<#��<$o�<%�d<&'<3~�<O�<$�<.p<9��<-D�<'�e<-Yc<&\<$@|<%F<%��<A�S<%��<#ߜ<%(<'<#��<#��<#��<#�<#��<#�E<#�N<#�<#�r<#�N<#��<$�<#��<$�7<$<<#�(<#��<#ף<#�$<#�o<#�i<#�&<#�N<#�<#�!<#�U<#�<$f�<$<#�&<#�<#�X<$�<#�<#�<#�<#׺<#��<#��<#��<#��<#�c<#�C<#��<#�I<$<<#��<$�<#�!<#؄<$�<#�^<#�$<$�<$�<#�E<#��<#��<#׎<#�<#ۮ<$	�<$3U<$#(<$O�<$"2<#�m<#��<#�<#�&<#�<#�0<#��<#ۮ<#؄<#��<#ף<#��<#��<#��<#�8<#�i<#�<#ף<#׎<#�X<$N�<$@|<#��<#��<#�<#�i<#ܯ<#�<#�<#�<#�&<#�l<#�8<#�<#ף<#׎<#�
<#׎<#�<#�+<#�m<#��<$<$'<#�!<#�<#�{<#�<#�<#��<#�<$�<$�<#�m<$k�<$.<#�<#�U<#�<#�<#�<#ٛ<#��<#�N<#�<#�
<#�<#�<#�r<#�D<#ٛ<#�<#�a<#�N<#�o<#׺<#�^<#׎<#�$<#�<#ף<#�E<#��<#�
<#�
<#׎<#ا<#�<$k<$e.<#ا<#�]<#׺<#�<#��<#�<#�{<#�o<#�<#��<$�<$/<#��<#ޫ<#��<#�<#�<#�o<#��<#�
<#�C<$
�<&h�<$	<#�8<#�r<#�"<$�<$U�<$�2<#�<#�(<#��<#��<#��<#�0<#��<$� <%�<(�\<#�8<#��<#�<$b�<$(<#��<$	�<#�<#�H<$H�<#��<#��<$<<#�g<#��<#ߜ<$	�<%��<#�<#�e<#�&<#�l<#׺<#�<#�C<#�4<$2G<#��<#�U<#��<#��<$.<#ܯ<#��<#��<#�o<#�<$
�<#��<#ף<#�<#��<$��<%K:<$.<#��<$�<$�;<#�<#�8<#��<#�<#��<#�{<#�<#��<#�M<$��<%͍<#�)<#��<#��<#ڑ<$G<#�l<#�M<%��<%y<$�<(�T<$z�<#��<#�5<#�Q<$�<#�"<#ٛ<#�<#��<#��<#��<#�&<$��<&�<$r�<#ߜ<#�<#�<#�o<#ܯ<#�	<$�<$N�<#�<#�<$�<$�J<$j|<$Sa<$-<$k<#�W<#�4<#�4<#��<#�i<#�<#��<$,<$O�<$m,<$�<#��<#�N<#�(<$T�<*��<(�_<#��<#�	<%�~<$z�<$#(<#��<#ۮ<#�m<$ �<$/<#��<#��<#�*<#�l<#�g<$6�<#�<$�<#�<&�<&,f<$?[<$�<$G<$r�<$��<$g�<$�w<&��<$F9<#��<#ף<#ڑ<#׎<#��<#�<#�D<#׺<$f<#��<#��<#�c<#�<$��<$i&<%<#�o<#�D<#�N<$ub<% <#�&<#��<#�^<#�o<#��<#�*<#ٛ<#׎<#�<#ޫ<#ޫ<#��<#�	<$3U<#ߜ<#�<#�<#��<#�X<#��<#��<#��<$�<$��<$.<#�*<#ۮ<#�8<#�8<#ڑ<#�<#׎<#�r<#�i<#��<#�<#�E<#�<#��<#�J<#�<#�<#�r<#�<#؄<#�{<#�i<#�i<#�<#�<#�{<#�X<#�<#�<#�<#�<#�{<#�{<#ף<#�<#�<#�<#ף<#�<#�c<#�<#�{<#׎<#��<#�<#ף<#�i<#�i<#�0<#�<#�c<#�<#׎<#�
<#��<#�8<#�<#�D<#�<#�<#�<#�<#�{<#�I<#��<#�^<#�<#ף<#�
<#��<#�<#�D<#ף<#�<#�<#�<#ا<#��<#��<#ا<#�]<#�]<#ۮ<#�<#�8<#�J<#��<#�J<$v<#��<#�J<#��<$<#�8<#��<#��<#�o<#��<#��<#�
<#�X<#�<#�{<#׺<#��<#ߜ<#�^<$P�<$	<#�<#��<#�<#ڑ<#�X<#��<#�l<#��<#�<#�{<#��<#�r<#�D<#�<#�<#��<#��<#�l<#�8<#��<#��<#��<#�<#�+<#�+<#�o<#��<#�E<$ <#��<#�<<#�
<#�
<#�
<#�{<#�<#ף<#ף<#�$<#�N<$ح<$$<#�c<#�!<$+<#�8<#�)<#׺<#�D<#�D<#��<#�Q<#�<#�"<#��<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285620181106012856QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                