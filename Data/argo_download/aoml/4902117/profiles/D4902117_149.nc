CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92020-07-05T02:00:44Z creation; 2022-08-01T16:26:22Z DMQC;      
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
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  \<   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  d    PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  ��   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �$   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ʴ   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  Ҙ   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �(   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     � �   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 � 9,   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     � A   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 � `�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     � h�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   �t   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   �t   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �t   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �t   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   �h   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��Argo profile    3.1 1.2 19500101000000  20200705020044  20220801122622  4902117 4902117 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               �   �AA  AOAO6730                            6730                            2C  2C  DD  S2A                             S2A                             7359                            7359                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�#����@�#����11  @�#����P@�#����P@N&���[B@N&���[B�?�)
�D��?�)
�D�11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�=q@   @@  @�G�@�G�@�G�@�G�A ��A��A$z�A@��A`��A�  A�  A�  A��A�  A�Q�A�  A�A��B�
B  B(�B (�B(  B0  B7�
B?�
BG33BO�
BX  B`(�Bh(�Bp(�BxQ�B�(�B�  B��
B�  B�  B�{B�{B�  B�(�B�  B��
B��B�{B�  B��B�  B�(�B�{B��B��
B��B��B�{B�(�B�  B��B�B��B�  B��B�B��C 
=C��C  C  C�C
  C�C{C{C{C{C
=C��C�C  C{C 
=C"  C#��C%�C(  C*
=C,  C.
=C0  C2  C3��C5�C7�HC9��C<  C>  C@  CB
=CD
=CE��CG��CJ  CL  CN  CP
=CR
=CT
=CV
=CW��CY��C\
=C]��C`  Cb
=Cc��Ce��Ch
=Ci��Cl  Cn
=Cp  Cr  Cs��Cu�Cw�HCz  C|{C~
=C�C�
=C���C���C���C�  C�  C�C�
=C�  C���C��C��C��C���C���C���C���C�  C�  C�C�C�C�  C�  C�C���C��C���C�C�
=C�C�  C���C���C�C�  C���C�C�C���C�C�  C���C�C�  C���C�C�
=C�
=C�C���C���C���C�  C�  C�  C�  C���C���C���C�  C�C���C�C�
=C�  C�  C�
=C�  C���C���C��C���C���C���C�C���C�  C�  C���C���C�  C���C���C�C�
=C�C�C�C�C�
=C���C���C�C�  C���C���C���C���C���C�C�C�C�C���C�  C�  C�
=C�  C�  C�C�  C���C�C�  C���C�  C�  C�C�
=C�  C�  C�C���C�C�
=C�  D   D ��D�D}qD�qD��D  D}qD  D�DD��D�qD}qD�qDz�D  D�D	D	�D
D
�DD��D  D��DD� D�qD��DD�DD�D�qDz�D�qD� D  D��D  D}qD�qD� D�D}qD�qD}qD��D��D�qDz�D�qD��DD��D  D� D  D� D  D� D  D��D�qD }qD!  D!��D"�D"��D#�D#z�D#��D$��D%�D%��D&�D&}qD&�qD'}qD(  D(z�D(�qD)}qD)�qD*��D+�D+��D,�D,� D-�D-� D-�qD.� D/  D/� D0  D0� D0�qD1}qD2  D2z�D2�qD3}qD3�qD4}qD4�qD5}qD6  D6��D7�D7��D8�D8�D9�D9��D:�D:��D;�D;}qD;�qD<}qD<�qD=� D=��D>}qD?  D?� D@D@�D@�qDA}qDB  DB� DC  DCz�DC�qDD� DE�DE� DF  DF� DF�qDG}qDG��DH��DI  DI� DJ  DJ}qDK  DK� DL  DL� DM�DM��DN  DN� DO  DO}qDO�qDP� DQ�DQ}qDR  DR��DR�qDS}qDS�qDT}qDU  DU��DV  DV}qDV�qDW}qDX  DX� DY  DY� DZ�DZ}qDZ�qD[� D\�D\�D]  D]� D^  D^� D_�D_}qD`  D`��D`�qDa}qDa�qDb� DcDc��Dc�qDd� De  De� Df�Df}qDg  Dg�Dh  Dh��DiDi� Dj  Dj� Dj�qDk� Dk��Dl}qDl�qDm� DnDn�Do�Do}qDo��Dp}qDq  Dq}qDq��Drz�Dr��Ds� DtDt��Du�Du� Du�RDv}qDw�Dw� Dx  Dx}qDx�qDy}qDy��Dz� D{D{��D|�D|��D}D}��D~  D~z�D~��Dz�D�qD�@ D��HD�� D�  D�>�D�� D�D�  D�=qD�� D���D�  D�B�D���D��HD��D�B�D���D�� D���D�@ D��HD�D��D�@ D�~�D��HD�  D�>�D�}qD��qD�HD�B�D�~�D���D���D�>�D��HD�� D�  D�@ D��HD�� D�  D�@ D��HD��HD���D�=qD�~�D��HD��D�AHD�~�D�� D�  D�@ D�~�D�� D�HD�@ D�� D��HD�HD�@ D�~�D���D���D�@ D��HD���D�HD�B�D�� D�� D�HD�>�D�~�D�� D���D�=qD�}qD��qD��)D�=qD�}qD��)D�  D�AHD�� D��HD�HD�@ D��HD�� D���D�=qD�}qD��qD��)D�=qD�~�D���D���D�=qD�}qD���D���D�>�D�� D�� D�HD�AHD��HD�D��D�AHD�~�D��HD�HD�AHD��HD�� D���D�AHD�~�D��)D��qD�>�D�� D��HD�  D�>�D�~�D���D���D�@ D�� D�� D�  D�=qD�� D�� D�  D�@ D�~�D�� D�  D�>�D�� D���D���D�>�D�~�D���D���D�AHD��HD�� D�HD�B�D�� D��qD��qD�>�D�� D�� D�  D�@ D��HD���D���D�AHD�~�D��qD���D�>�D�~�D��HD�  D�=qD�~�D���D���D�>�D�~�D�� D�HD�B�D�� D��qD�  D�AHD�~�D�� D�  D�@ D�� D�� D���D�>�D���D��HD��qD�>�D�� D���D�HD�AHD�� D���D���D�>�D�~�D��qD��qD�<)D�|)D�� D��D�B�D��HD��HD��D�B�D���D�D��D�AHD�� D�� D���D�<)D�|)D��)D��qD�=qD�}qD���D�HD�B�D���D��HD�HD�B�D���D�� D��)D�=qD�}qD��qD��qD�=qD�� D��HD��D�B�D�~�D��)D���D�@ D��HD��HD���D�<)D�~�D�� D��D�C�DÃ�D���D��D�AHD�~�D��HD��D�>�D�|)DŽqD�  D�B�Dƀ DƽqD���D�AHDǀ DǾ�D�HD�@ D�~�D��HD�  D�>�DɁHD�� D��qD�AHDʀ DʽqD�  D�@ Dˀ D��HD��D�B�D́HD��HD�  D�@ D�}qD;�D�HD�AHD΀ D�� D�HD�>�Dπ D�D�  D�=qD�}qDо�D�HD�@ Dр D�D�  D�>�DҀ DҾ�D���D�AHDӁHDӾ�D���D�@ DԀ D�� D���D�>�DՀ Dվ�D�  D�AHDր DֽqD�  D�@ DׁHD��HD���D�@ D؀ D��HD��D�@ D�}qDٽqD���D�AHD�~�Dڼ)D��qD�>�DہHD�D�  D�=qD�}qD�� D�HD�@ D�|)Dݼ)D��qD�>�Dހ D�� D�HD�B�D߀ D߽qD���D�@ D�~�D�qD�  D�@ D�}qD�� D�  D�AHD₏D�� D���D�@ D�~�D㾸D�HD�AHD� D�D�HD�@ D�HD�� D��qD�>�D� D�D�  D�>�D� D�D�HD�@ D�HD�� D��qD�@ D� D��HD��D�@ D�}qD꾸D���D�>�D�HD�D�HD�>�D�~�D�� D�HD�@ D�HD���D���D�@ D�~�D�� D��D�AHD�~�D�qD�  D�B�D���D��HD�  D�>�D�~�D�� D�  D�@ D�D�� D�  D�AHD�HD��HD�  D�>�D� D�� D��qD�@ D���D�D��D�@ D�}qD���D���D�>�D�� D�D�HD�AHD��HD�� D���D�>�D�~�D���D�  D�B�D��HD�� D�  D�=q?W
=?k�?��?�z�?��
?�33?\?���?�(�?�?�@�\@
=q@�@
=@(�@#�
@(��@0��@5@=p�@E�@J=q@Q�@W
=@^�R@fff@n{@u@z�H@�G�@��
@��@��@�\)@�33@�
=@��H@��R@�G�@��@���@���@�\)@�33@�
=@���@�p�@�G�@��@Ǯ@˅@�\)@�33@�
=@ٙ�@�p�@�G�@��@���@�@�\)@�33@�@���@�p�A ��A�\A�
AA�A	��A�A��A�RA��A�\Az�AA�A��A�Ap�A\)A!G�A#33A$z�A&ffA(Q�A*=qA,(�A.{A0  A1�A333A5�A7
=A8��A:�HA<��A>�RA@��AB�\ADz�AFffAG�AI��AK�AMp�AO\)AQG�AR�\ATz�AVffAXQ�AZ=qA\(�A^{A`  AaG�Ac33Ae�Ag
=Ah��Aj�HAl(�An{Ap  Aq�As33Au�Aw
=Ax��Az�HA|(�A~{A�  A���A��A��\A��A�z�A�p�A�ffA�
=A�  A���A���A��\A��A�z�A�p�A�ffA�
=A�  A���A���A��\A��A�z�A�p�A�{A�
=A�  A���A���A��\A��A�z�A�p�A�{A�\)A�  A���A��A��HA��
A�z�A�p�A�ffA�\)A�  A���A��A��HA��A�z�A�p�A�ffA�
=A�  A�G�A�=qA��HA��
A���A�A��RA��A���A�G�A�=qA�33A�(�A��A�A��RA��A���A���A\A�33A�(�A��A�AƸRAǮAȣ�Aə�Aʏ\A�33A�(�A��A�{AθRA�  AУ�Aљ�Aҏ\AӅA�(�A��A�{A�
=A�  A���Aٙ�Aڏ\AۅA�z�A�p�A�{A�
=A�  A���A��A�\A�A�z�A�p�A�ffA�\)A�  A�G�A��A��HA��
A���A�A�RA�A��A�A�\A�A�z�A�p�A�ffA�
=A�  A�G�A��A�33A��
A���A�A��RA��B Q�B ��BG�BB=qB�RB33B�B  B��B��Bp�B�BffB�HB\)B�
BQ�B��B	G�B	��B
=qB
�\B33B�B  Bz�B��Bp�B�BffB�HB\)B�
BQ�B��BG�B��B{B�\B
=B�B  BQ�B��BG�BB{B�RB
=B�B  Bz�B��BG�BB=qB�RB
=B�B  Bz�B��Bp�BB=qB�RB33B�B (�B ��B ��B!p�B!�B"ffB"�RB#33B#�B$(�B$��B$��B%p�B%�B&=qB&�RB'33B'�B(  B(z�B(��B)G�B)B*=qB*�RB+
=B+�B,  B,z�B,��B-G�B-B.=qB.�\B/
=B/\)B/�
B0Q�B0��B1�B1��B2{B2ffB2�HB3\)B3�B4(�B4��B4��B5p�B5�B6=qB6�RB7
=B7�B7�
B8Q�B8��B9�B9��B9�B:ffB:�RB;33B;�B<  B<z�B<��B=G�B=��B>{B>ffB>�HB?\)B?�
B@(�B@��B@��BAp�BA�BBffBB�HBC33BC�BD(�BD��BE�BEp�BE�BFffBF�HBG\)BG�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ?�=q@   @@  @�G�@�G�@�G�@�G�A ��A��A$z�A@��A`��A�  A�  A�  A��A�  A�Q�A�  A�A��B�
B  B(�B (�B(  B0  B7�
B?�
BG33BO�
BX  B`(�Bh(�Bp(�BxQ�B�(�B�  B��
B�  B�  B�{B�{B�  B�(�B�  B��
B��B�{B�  B��B�  B�(�B�{B��B��
B��B��B�{B�(�B�  B��B�B��B�  B��B�B��C 
=C��C  C  C�C
  C�C{C{C{C{C
=C��C�C  C{C 
=C"  C#��C%�C(  C*
=C,  C.
=C0  C2  C3��C5�C7�HC9��C<  C>  C@  CB
=CD
=CE��CG��CJ  CL  CN  CP
=CR
=CT
=CV
=CW��CY��C\
=C]��C`  Cb
=Cc��Ce��Ch
=Ci��Cl  Cn
=Cp  Cr  Cs��Cu�Cw�HCz  C|{C~
=C�C�
=C���C���C���C�  C�  C�C�
=C�  C���C��C��C��C���C���C���C���C�  C�  C�C�C�C�  C�  C�C���C��C���C�C�
=C�C�  C���C���C�C�  C���C�C�C���C�C�  C���C�C�  C���C�C�
=C�
=C�C���C���C���C�  C�  C�  C�  C���C���C���C�  C�C���C�C�
=C�  C�  C�
=C�  C���C���C��C���C���C���C�C���C�  C�  C���C���C�  C���C���C�C�
=C�C�C�C�C�
=C���C���C�C�  C���C���C���C���C���C�C�C�C�C���C�  C�  C�
=C�  C�  C�C�  C���C�C�  C���C�  C�  C�C�
=C�  C�  C�C���C�C�
=C�  D   D ��D�D}qD�qD��D  D}qD  D�DD��D�qD}qD�qDz�D  D�D	D	�D
D
�DD��D  D��DD� D�qD��DD�DD�D�qDz�D�qD� D  D��D  D}qD�qD� D�D}qD�qD}qD��D��D�qDz�D�qD��DD��D  D� D  D� D  D� D  D��D�qD }qD!  D!��D"�D"��D#�D#z�D#��D$��D%�D%��D&�D&}qD&�qD'}qD(  D(z�D(�qD)}qD)�qD*��D+�D+��D,�D,� D-�D-� D-�qD.� D/  D/� D0  D0� D0�qD1}qD2  D2z�D2�qD3}qD3�qD4}qD4�qD5}qD6  D6��D7�D7��D8�D8�D9�D9��D:�D:��D;�D;}qD;�qD<}qD<�qD=� D=��D>}qD?  D?� D@D@�D@�qDA}qDB  DB� DC  DCz�DC�qDD� DE�DE� DF  DF� DF�qDG}qDG��DH��DI  DI� DJ  DJ}qDK  DK� DL  DL� DM�DM��DN  DN� DO  DO}qDO�qDP� DQ�DQ}qDR  DR��DR�qDS}qDS�qDT}qDU  DU��DV  DV}qDV�qDW}qDX  DX� DY  DY� DZ�DZ}qDZ�qD[� D\�D\�D]  D]� D^  D^� D_�D_}qD`  D`��D`�qDa}qDa�qDb� DcDc��Dc�qDd� De  De� Df�Df}qDg  Dg�Dh  Dh��DiDi� Dj  Dj� Dj�qDk� Dk��Dl}qDl�qDm� DnDn�Do�Do}qDo��Dp}qDq  Dq}qDq��Drz�Dr��Ds� DtDt��Du�Du� Du�RDv}qDw�Dw� Dx  Dx}qDx�qDy}qDy��Dz� D{D{��D|�D|��D}D}��D~  D~z�D~��Dz�D�qD�@ D��HD�� D�  D�>�D�� D�D�  D�=qD�� D���D�  D�B�D���D��HD��D�B�D���D�� D���D�@ D��HD�D��D�@ D�~�D��HD�  D�>�D�}qD��qD�HD�B�D�~�D���D���D�>�D��HD�� D�  D�@ D��HD�� D�  D�@ D��HD��HD���D�=qD�~�D��HD��D�AHD�~�D�� D�  D�@ D�~�D�� D�HD�@ D�� D��HD�HD�@ D�~�D���D���D�@ D��HD���D�HD�B�D�� D�� D�HD�>�D�~�D�� D���D�=qD�}qD��qD��)D�=qD�}qD��)D�  D�AHD�� D��HD�HD�@ D��HD�� D���D�=qD�}qD��qD��)D�=qD�~�D���D���D�=qD�}qD���D���D�>�D�� D�� D�HD�AHD��HD�D��D�AHD�~�D��HD�HD�AHD��HD�� D���D�AHD�~�D��)D��qD�>�D�� D��HD�  D�>�D�~�D���D���D�@ D�� D�� D�  D�=qD�� D�� D�  D�@ D�~�D�� D�  D�>�D�� D���D���D�>�D�~�D���D���D�AHD��HD�� D�HD�B�D�� D��qD��qD�>�D�� D�� D�  D�@ D��HD���D���D�AHD�~�D��qD���D�>�D�~�D��HD�  D�=qD�~�D���D���D�>�D�~�D�� D�HD�B�D�� D��qD�  D�AHD�~�D�� D�  D�@ D�� D�� D���D�>�D���D��HD��qD�>�D�� D���D�HD�AHD�� D���D���D�>�D�~�D��qD��qD�<)D�|)D�� D��D�B�D��HD��HD��D�B�D���D�D��D�AHD�� D�� D���D�<)D�|)D��)D��qD�=qD�}qD���D�HD�B�D���D��HD�HD�B�D���D�� D��)D�=qD�}qD��qD��qD�=qD�� D��HD��D�B�D�~�D��)D���D�@ D��HD��HD���D�<)D�~�D�� D��D�C�DÃ�D���D��D�AHD�~�D��HD��D�>�D�|)DŽqD�  D�B�Dƀ DƽqD���D�AHDǀ DǾ�D�HD�@ D�~�D��HD�  D�>�DɁHD�� D��qD�AHDʀ DʽqD�  D�@ Dˀ D��HD��D�B�D́HD��HD�  D�@ D�}qD;�D�HD�AHD΀ D�� D�HD�>�Dπ D�D�  D�=qD�}qDо�D�HD�@ Dр D�D�  D�>�DҀ DҾ�D���D�AHDӁHDӾ�D���D�@ DԀ D�� D���D�>�DՀ Dվ�D�  D�AHDր DֽqD�  D�@ DׁHD��HD���D�@ D؀ D��HD��D�@ D�}qDٽqD���D�AHD�~�Dڼ)D��qD�>�DہHD�D�  D�=qD�}qD�� D�HD�@ D�|)Dݼ)D��qD�>�Dހ D�� D�HD�B�D߀ D߽qD���D�@ D�~�D�qD�  D�@ D�}qD�� D�  D�AHD₏D�� D���D�@ D�~�D㾸D�HD�AHD� D�D�HD�@ D�HD�� D��qD�>�D� D�D�  D�>�D� D�D�HD�@ D�HD�� D��qD�@ D� D��HD��D�@ D�}qD꾸D���D�>�D�HD�D�HD�>�D�~�D�� D�HD�@ D�HD���D���D�@ D�~�D�� D��D�AHD�~�D�qD�  D�B�D���D��HD�  D�>�D�~�D�� D�  D�@ D�D�� D�  D�AHD�HD��HD�  D�>�D� D�� D��qD�@ D���D�D��D�@ D�}qD���D���D�>�D�� D�D�HD�AHD��HD�� D���D�>�D�~�D���D�  D�B�D��HD�� D�  D�=q?W
=?k�?��?�z�?��
?�33?\?���?�(�?�?�@�\@
=q@�@
=@(�@#�
@(��@0��@5@=p�@E�@J=q@Q�@W
=@^�R@fff@n{@u@z�H@�G�@��
@��@��@�\)@�33@�
=@��H@��R@�G�@��@���@���@�\)@�33@�
=@���@�p�@�G�@��@Ǯ@˅@�\)@�33@�
=@ٙ�@�p�@�G�@��@���@�@�\)@�33@�@���@�p�A ��A�\A�
AA�A	��A�A��A�RA��A�\Az�AA�A��A�Ap�A\)A!G�A#33A$z�A&ffA(Q�A*=qA,(�A.{A0  A1�A333A5�A7
=A8��A:�HA<��A>�RA@��AB�\ADz�AFffAG�AI��AK�AMp�AO\)AQG�AR�\ATz�AVffAXQ�AZ=qA\(�A^{A`  AaG�Ac33Ae�Ag
=Ah��Aj�HAl(�An{Ap  Aq�As33Au�Aw
=Ax��Az�HA|(�A~{A�  A���A��A��\A��A�z�A�p�A�ffA�
=A�  A���A���A��\A��A�z�A�p�A�ffA�
=A�  A���A���A��\A��A�z�A�p�A�{A�
=A�  A���A���A��\A��A�z�A�p�A�{A�\)A�  A���A��A��HA��
A�z�A�p�A�ffA�\)A�  A���A��A��HA��A�z�A�p�A�ffA�
=A�  A�G�A�=qA��HA��
A���A�A��RA��A���A�G�A�=qA�33A�(�A��A�A��RA��A���A���A\A�33A�(�A��A�AƸRAǮAȣ�Aə�Aʏ\A�33A�(�A��A�{AθRA�  AУ�Aљ�Aҏ\AӅA�(�A��A�{A�
=A�  A���Aٙ�Aڏ\AۅA�z�A�p�A�{A�
=A�  A���A��A�\A�A�z�A�p�A�ffA�\)A�  A�G�A��A��HA��
A���A�A�RA�A��A�A�\A�A�z�A�p�A�ffA�
=A�  A�G�A��A�33A��
A���A�A��RA��B Q�B ��BG�BB=qB�RB33B�B  B��B��Bp�B�BffB�HB\)B�
BQ�B��B	G�B	��B
=qB
�\B33B�B  Bz�B��Bp�B�BffB�HB\)B�
BQ�B��BG�B��B{B�\B
=B�B  BQ�B��BG�BB{B�RB
=B�B  Bz�B��BG�BB=qB�RB
=B�B  Bz�B��Bp�BB=qB�RB33B�B (�B ��B ��B!p�B!�B"ffB"�RB#33B#�B$(�B$��B$��B%p�B%�B&=qB&�RB'33B'�B(  B(z�B(��B)G�B)B*=qB*�RB+
=B+�B,  B,z�B,��B-G�B-B.=qB.�\B/
=B/\)B/�
B0Q�B0��B1�B1��B2{B2ffB2�HB3\)B3�B4(�B4��B4��B5p�B5�B6=qB6�RB7
=B7�B7�
B8Q�B8��B9�B9��B9�B:ffB:�RB;33B;�B<  B<z�B<��B=G�B=��B>{B>ffB>�HB?\)B?�
B@(�B@��B@��BAp�BA�BBffBB�HBC33BC�BD(�BD��BE�BEp�BE�BFffBF�HBG\)BG�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�A=qA$�A�PA
�\A	G�A�Ax�A"�A�A�\A(�A�A�HA��A+A�^AC�A%A ��A I�A �@��;@���@�t�@�"�@��!@� �@�@��!@�r�@�bN@�+@�/@�Z@Ώ\@�@��@��#@�X@Ý�@�@���@���@�K�@�$�@���@���@��H@�E�@���@� �@�$�@���@�  @�@�v�@�^5@��T@��h@���@�`B@��@��`@��j@��u@���@� �@�ȴ@���@��@�-@�
=@��@��@��R@�$�@��@��`@�Ĝ@�Q�@���@��@��\@�^5@�V@�V@�E�@�$�@���@��h@�`B@��@���@��@��/@�Ĝ@�z�@�  @��@��;@��w@���@��F@�C�@��R@��\@�v�@�-@��T@��h@�x�@�G�@��@��@��7@�hs@�O�@��!@�ff@�{@���@�X@��@�%@��@���@��@�J@��+@�;d@��
@�1'@���@�&�@�x�@�ff@�o@�l�@�o@��R@�J@�V@��D@��u@��;@��@��T@�ƨ@�hs@�  @��H@�M�@���@�/@�Ĝ@���@���@�(�@���@�K�@��y@���@��R@���@�ff@�=q@�@��@��@��@��7@�7L@���@� �@��;@��w@���@��@�dZ@�+@��H@��\@�~�@�v�@�E�@�5?@��@�J@���@�$�@���@��@��
@��7@��@�9X@���@��@��+@�ȴ@�"�@�;d@�\)@�l�@�
=@��\@�ff@���@�O�@�Ĝ@�Z@�9X@�9X@�1'@��;@�+@�V@�I�@� �@��
@��w@��F@���@�t�@�"�@���@�l�@���@�=q@�?}@���@���@��@�1'@�1'@�I�@�j@�j@�bN@�A�@��@��m@��m@��;@���@��P@�\)@��@�?}@���@���@���@�7L@���@���@�C�@�7L@�ȴ@��@���@�-@��#@�^5@��@���@�v�@��-@���@��D@�bN@�r�@��9@��9@��@�1'@��m@��+@���@�G�@���@�bN@�A�@��@��@��@�b@�  @�@�@�@�  @�(�@�j@�z�@���@��u@�%@�p�@�hs@��@��`@��`@��@�p�@��-@��-@��-@���@���@���@��^@�@�@��h@�p�@�?}@�%@���@�Ĝ@�Ĝ@��/@��9@�A�@��@��@�b@�1@�b@�1'@��9@���@���@��j@+@}��@}�@}O�@}V@|�@|�j@|j@{�
@{33@z��@z=q@zJ@y�@y�^@yhs@y�^@y�#@y�^@yhs@x��@x��@x�@xb@x �@w�;@w�@w��@w�P@w|�@w\)@wK�@wK�@w\)@w\)@wl�@w|�@w��@w�@w��@w��@w�P@w;d@v�y@vȴ@w�@wK�@w+@v��@v�y@v�y@v�@v��@u�@u�@t�@t�@t�/@t��@t�j@tj@t(�@t�D@tj@tI�@s�
@sdZ@sdZ@r�\@q�@q�^@qx�@qX@qG�@qhs@qhs@qhs@qhs@qhs@qX@qG�@p��@p �@p  @o��@pĜ@p�`@p��@p��@p  @o�w@o��@o�@p1'@o�@o�;@o�@n�y@l�/@l�/@n5?@nff@nȴ@o�@o;d@o�P@o�P@o+@nȴ@nff@n$�@m�T@m�h@m�@lI�@k��@j�H@j^5@j=q@j�@iX@h��@h�9@h��@h�u@hQ�@g�@g�;@g�@g��@gK�@g�@f�R@f��@f�+@fff@f$�@f@e�@e�@e��@e?}@d�/@d��@dZ@dZ@dZ@dI�@d�@ct�@c��@dj@c��@c�@cS�@co@co@c"�@b�H@bn�@bn�@b=q@bM�@b-@a�^@a��@aG�@a�@aG�@`�9@`Q�@`A�@`A�@`A�@`A�@`1'@_�@_�@^�R@^��@^�+@^ff@^ff@^E�@]�T@]�-@]�-@]�@]p�@]?}@]?}@]?}@]?}@]O�@]O�@]?}@]?}@]?}@]?}@]?}@]�-@^E�@^�+@^v�@^�+@^�@_+@_
=@_
=@^�y@^��@^$�@^{@^E�@^$�@]�h@]?}@\�@\�@]V@]V@]�@]/@]?}@]`B@]�h@]�-@^@^$�@^$�@^{@^@^$�@^{@]��@\��@\�j@\��@\��@]`B@]@]�T@]��@]�-@]��@]��@]��@^{@^$�@^5?@]�@]�-@]?}@\��@\�@\�D@\(�@\1@\9X@\j@\I�@\9X@[ƨ@[��@[ƨ@[��@\�@[��@\�D@^��@_
=@^�@^��@^��@^�+@^ff@^�+@_�w@`r�@`�9@`��@`�9@`Ĝ@`��@`bN@_�@]?}@\9X@\1@\(�@\9X@\�D@\�/@\�@]`B@^��@_;d@_\)@_|�@_�@_�;@`bN@`Ĝ@`�u@`Q�@`b@_�@`b@`Q�@`��@`Ĝ@`�9@`��@`Ĝ@`Q�@`bN@`�`@a&�@a�@aG�@ax�@a��@a�^@a�#@b�@bn�@b��@cdZ@c�m@c�m@c�m@c�
@c��@d1@d1@d1@d�@d(�@d9X@d9X@dZ@d�D@d�j@d��@d�/@d�@d��@d��@d9X@c33@bn�@bM�@b�@a��@a��@a�^@a��@a��@ax�@ax�@a�7@a��@a��@a�@a��@bJ@b-@b^5@b�!@b�@dz�@fE�@fE�@f$�@f$�@f5?@f��@f5?@f$�@fff@fE�@fE�@e�@e��@e�h@e�@e�h@e@eV@e��@e�-@d�@d�@d�D@d�@ep�@e��@fv�@fff@f�@g|�@g�;@hb@h1'@h �@g�;@g�@gl�@g+@f�@f�@g
=@g
=@g\)@g;d@f��@f�R@f��@f�+@fff@fff@fff@fff@fff@fv�@fv�@fv�@f�+@fv�@fff@fff@fff@fv�@fv�@fff@fV@f5?@f$�@f5?@f$�@f{@f@e�@e��@e@e@e��@e�h@e�@eO�@eV@d��@d�/@d�@dz�@dj@dI�@d�@d1@c��@c��@c��@c�m@c�@c�@ct�@ct�@ct�@cS�@cC�@cC�@c"�@co@c@b�@b�H@b��@b��@b��@b�!@b��@b��@b��@b�\@bn�@bn�@b^5@bM�@bM�@b=q@b�@a��@a�@a�@a��@ax�@ahs@aX@aX@aX@a&�@`�`@`��@`Ĝ@`�9@`r�@`r�@`bN@`A�@`b@_�;@_��@_��@_�w@_��@_|�@_;d@_�@_
=@^��@^�R@^�R@^�R@^�R@^��@^��@^��@^��@^V@^@]�-@]�h@]�@]O�@]�@\��@\�@\��@\��@\�D@\j@\(�@\1@[��@[��@[��@[�m@[�
@[�
@[�
@[�F@[�F@[�@[S�@[C�@[33@[o@[o@[@Z�@Z�H@Z��@Z��@Z�!@Z��@Z�\@Z^5@Z^5@ZM�@Z-@Y��@Y�^@Y��@Y��@Y�7@Yx�@Y7L@X��@XĜ@X��@X�u@X�u@X�@X�@X�@XQ�@W�;@W��@W\)@W+@V�y@V��@V��@V�@Vȴ@V��@V��@V��@Vv�@VE�@V5?@V{@U�@U�@U��@U��@U?}@U�@UV@T��@T�@T�D@Tz�@Tz�@Tz�@Tj@Tj@Tj@TZ@TZ@TZ@T�@Sƨ@S��@SC�@R�@R�H@R��@R�\@R~�@R=q@RJ@RJ@RJ@RJ@RJ@Q��@Q�#@Q��@Q��@Q�^@Q�^@Q��@Q�7@Q&�@Q%@P��@P��@P��@Q%@Q%@Q%@P��@P��@P��@PĜ@PĜ@P��@P�u@P�@Pr�@PbN@PA�@P �@P  @O�;@O��@O��@O��@O��@O�w@O�w@O�w@O��@O��@O�@O��AI�AE�AA�A=qA5?A1'A-A$�A(�A$�A�A(�A(�A �A �A�A1A�A��A�^A��A|�AS�A?}A+AoA
��A
�HA
��A
�9A
��A
�A
bNA
5?A
  A	ƨA	�-A	�PA	p�A	S�A	/A	%A��A��A~�AbNA=qA(�A{AA�A�;AƨA�-A��A�hA�At�AhsA`BAXAO�AG�A?}A33A&�A�A�AoAoAVA
=A%AA��A��A�A�`A�AȴA�9A�A�A��A��A��A�uA�DA�Az�Ar�AffAbNAVAQ�AM�AE�A9XA1'A(�A �A�A�A�A�A{A{AbAJA1A  A��A�A�TA�
AƨA�FA��A�Al�AXAS�AC�A7LA+A"�A�AoA
=A%AAAA��A��A��A�A�yA�HA�A��A��A�!A��A�DAv�AbNAE�A(�AbA��A�mA��A�wA�A��A��A�hA�PA�7A�7A�A�A�A�A|�Ax�Ax�At�Ap�AhsA\)AO�A?}A+AoA��A�/A�RA�DAVA(�AJA�A�TA��AƨA�wA�^A�-A��A��A��A��A��A�hA�PA�A|�At�Al�AdZA\)AS�AK�AC�A?}A?}A7LA7LA/A/A/A/A/A+A+A&�A&�A"�A�A�A�A�AoA
=A
=AA ��A ��A �A �yA �`A �/A ��A ��A ȴA ĜA ��A �jA �RA �9A �A ��A ��A ��A ��A �DA �A z�A jA bNA ZA Q�A M�A M�A M�A M�A M�A I�A I�A E�A E�A E�A A�A A�A =qA =qA 5?A 1'A -A $�A $�A  �A �A �A {A {A {A bA bA JA JA 1A 1A A   @���@��@��@��@��m@��;@��;@��;@��;@��
@���@���@�ƨ@��w@��w@��F@��F@��F@��@��@��@��@��@��@���@��@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��P@��@�|�@�t�@�l�@�\)@�\)@�K�@�C�@�;d@�33@�33@�+@�+@�"�@�+@�+@�+@�+@�"�@�+@�"�@�"�@��@�o@�o@�
=@�
=@�
=@�@�@���@���@��@��@��H@��@���@���@��R@���@��\@�v�@�V@�5?@�{@��@���@��-@��h@�p�@�O�@�&�@��`@���@�I�@��@��@�l�@�+@��@���@���@���@��+@�v�@�v�@�ff@�V@�M�@�E�@�=q@�-@�$�@�J@��@��#@�@��^@��-@���@���@��h@��@�p�@�X@�7L@��@���@���@��u@�Q�@��@�l�@�ȴ@��@�`B@��@�D@�1@�\)@�ȴ@�n�@�5?@�@���@�h@�X@��@���@�r�@���G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          A=qA$�A�PA
�\A	G�A�Ax�A"�A�A�\A(�A�A�HA��A+A�^AC�A%A ��A I�A �@��;@���@�t�@�"�@��!@� �@�@��!@�r�@�bN@�+@�/@�Z@Ώ\@�@��@��#@�X@Ý�@�@���@���@�K�@�$�@���@���@��H@�E�@���@� �@�$�@���@�  @�@�v�@�^5@��T@��h@���@�`B@��@��`@��j@��u@���@� �@�ȴ@���@��@�-@�
=@��@��@��R@�$�@��@��`@�Ĝ@�Q�@���@��@��\@�^5@�V@�V@�E�@�$�@���@��h@�`B@��@���@��@��/@�Ĝ@�z�@�  @��@��;@��w@���@��F@�C�@��R@��\@�v�@�-@��T@��h@�x�@�G�@��@��@��7@�hs@�O�@��!@�ff@�{@���@�X@��@�%@��@���@��@�J@��+@�;d@��
@�1'@���@�&�@�x�@�ff@�o@�l�@�o@��R@�J@�V@��D@��u@��;@��@��T@�ƨ@�hs@�  @��H@�M�@���@�/@�Ĝ@���@���@�(�@���@�K�@��y@���@��R@���@�ff@�=q@�@��@��@��@��7@�7L@���@� �@��;@��w@���@��@�dZ@�+@��H@��\@�~�@�v�@�E�@�5?@��@�J@���@�$�@���@��@��
@��7@��@�9X@���@��@��+@�ȴ@�"�@�;d@�\)@�l�@�
=@��\@�ff@���@�O�@�Ĝ@�Z@�9X@�9X@�1'@��;@�+@�V@�I�@� �@��
@��w@��F@���@�t�@�"�@���@�l�@���@�=q@�?}@���@���@��@�1'@�1'@�I�@�j@�j@�bN@�A�@��@��m@��m@��;@���@��P@�\)@��@�?}@���@���@���@�7L@���@���@�C�@�7L@�ȴ@��@���@�-@��#@�^5@��@���@�v�@��-@���@��D@�bN@�r�@��9@��9@��@�1'@��m@��+@���@�G�@���@�bN@�A�@��@��@��@�b@�  @�@�@�@�  @�(�@�j@�z�@���@��u@�%@�p�@�hs@��@��`@��`@��@�p�@��-@��-@��-@���@���@���@��^@�@�@��h@�p�@�?}@�%@���@�Ĝ@�Ĝ@��/@��9@�A�@��@��@�b@�1@�b@�1'@��9@���@���@��j@+@}��@}�@}O�@}V@|�@|�j@|j@{�
@{33@z��@z=q@zJ@y�@y�^@yhs@y�^@y�#@y�^@yhs@x��@x��@x�@xb@x �@w�;@w�@w��@w�P@w|�@w\)@wK�@wK�@w\)@w\)@wl�@w|�@w��@w�@w��@w��@w�P@w;d@v�y@vȴ@w�@wK�@w+@v��@v�y@v�y@v�@v��@u�@u�@t�@t�@t�/@t��@t�j@tj@t(�@t�D@tj@tI�@s�
@sdZ@sdZ@r�\@q�@q�^@qx�@qX@qG�@qhs@qhs@qhs@qhs@qhs@qX@qG�@p��@p �@p  @o��@pĜ@p�`@p��@p��@p  @o�w@o��@o�@p1'@o�@o�;@o�@n�y@l�/@l�/@n5?@nff@nȴ@o�@o;d@o�P@o�P@o+@nȴ@nff@n$�@m�T@m�h@m�@lI�@k��@j�H@j^5@j=q@j�@iX@h��@h�9@h��@h�u@hQ�@g�@g�;@g�@g��@gK�@g�@f�R@f��@f�+@fff@f$�@f@e�@e�@e��@e?}@d�/@d��@dZ@dZ@dZ@dI�@d�@ct�@c��@dj@c��@c�@cS�@co@co@c"�@b�H@bn�@bn�@b=q@bM�@b-@a�^@a��@aG�@a�@aG�@`�9@`Q�@`A�@`A�@`A�@`A�@`1'@_�@_�@^�R@^��@^�+@^ff@^ff@^E�@]�T@]�-@]�-@]�@]p�@]?}@]?}@]?}@]?}@]O�@]O�@]?}@]?}@]?}@]?}@]?}@]�-@^E�@^�+@^v�@^�+@^�@_+@_
=@_
=@^�y@^��@^$�@^{@^E�@^$�@]�h@]?}@\�@\�@]V@]V@]�@]/@]?}@]`B@]�h@]�-@^@^$�@^$�@^{@^@^$�@^{@]��@\��@\�j@\��@\��@]`B@]@]�T@]��@]�-@]��@]��@]��@^{@^$�@^5?@]�@]�-@]?}@\��@\�@\�D@\(�@\1@\9X@\j@\I�@\9X@[ƨ@[��@[ƨ@[��@\�@[��@\�D@^��@_
=@^�@^��@^��@^�+@^ff@^�+@_�w@`r�@`�9@`��@`�9@`Ĝ@`��@`bN@_�@]?}@\9X@\1@\(�@\9X@\�D@\�/@\�@]`B@^��@_;d@_\)@_|�@_�@_�;@`bN@`Ĝ@`�u@`Q�@`b@_�@`b@`Q�@`��@`Ĝ@`�9@`��@`Ĝ@`Q�@`bN@`�`@a&�@a�@aG�@ax�@a��@a�^@a�#@b�@bn�@b��@cdZ@c�m@c�m@c�m@c�
@c��@d1@d1@d1@d�@d(�@d9X@d9X@dZ@d�D@d�j@d��@d�/@d�@d��@d��@d9X@c33@bn�@bM�@b�@a��@a��@a�^@a��@a��@ax�@ax�@a�7@a��@a��@a�@a��@bJ@b-@b^5@b�!@b�@dz�@fE�@fE�@f$�@f$�@f5?@f��@f5?@f$�@fff@fE�@fE�@e�@e��@e�h@e�@e�h@e@eV@e��@e�-@d�@d�@d�D@d�@ep�@e��@fv�@fff@f�@g|�@g�;@hb@h1'@h �@g�;@g�@gl�@g+@f�@f�@g
=@g
=@g\)@g;d@f��@f�R@f��@f�+@fff@fff@fff@fff@fff@fv�@fv�@fv�@f�+@fv�@fff@fff@fff@fv�@fv�@fff@fV@f5?@f$�@f5?@f$�@f{@f@e�@e��@e@e@e��@e�h@e�@eO�@eV@d��@d�/@d�@dz�@dj@dI�@d�@d1@c��@c��@c��@c�m@c�@c�@ct�@ct�@ct�@cS�@cC�@cC�@c"�@co@c@b�@b�H@b��@b��@b��@b�!@b��@b��@b��@b�\@bn�@bn�@b^5@bM�@bM�@b=q@b�@a��@a�@a�@a��@ax�@ahs@aX@aX@aX@a&�@`�`@`��@`Ĝ@`�9@`r�@`r�@`bN@`A�@`b@_�;@_��@_��@_�w@_��@_|�@_;d@_�@_
=@^��@^�R@^�R@^�R@^�R@^��@^��@^��@^��@^V@^@]�-@]�h@]�@]O�@]�@\��@\�@\��@\��@\�D@\j@\(�@\1@[��@[��@[��@[�m@[�
@[�
@[�
@[�F@[�F@[�@[S�@[C�@[33@[o@[o@[@Z�@Z�H@Z��@Z��@Z�!@Z��@Z�\@Z^5@Z^5@ZM�@Z-@Y��@Y�^@Y��@Y��@Y�7@Yx�@Y7L@X��@XĜ@X��@X�u@X�u@X�@X�@X�@XQ�@W�;@W��@W\)@W+@V�y@V��@V��@V�@Vȴ@V��@V��@V��@Vv�@VE�@V5?@V{@U�@U�@U��@U��@U?}@U�@UV@T��@T�@T�D@Tz�@Tz�@Tz�@Tj@Tj@Tj@TZ@TZ@TZ@T�@Sƨ@S��@SC�@R�@R�H@R��@R�\@R~�@R=q@RJ@RJ@RJ@RJ@RJ@Q��@Q�#@Q��@Q��@Q�^@Q�^@Q��@Q�7@Q&�@Q%@P��@P��@P��@Q%@Q%@Q%@P��@P��@P��@PĜ@PĜ@P��@P�u@P�@Pr�@PbN@PA�@P �@P  @O�;@O��@O��@O��@O��@O�w@O�w@O�w@O��@O��@O�@O��AI�AE�AA�A=qA5?A1'A-A$�A(�A$�A�A(�A(�A �A �A�A1A�A��A�^A��A|�AS�A?}A+AoA
��A
�HA
��A
�9A
��A
�A
bNA
5?A
  A	ƨA	�-A	�PA	p�A	S�A	/A	%A��A��A~�AbNA=qA(�A{AA�A�;AƨA�-A��A�hA�At�AhsA`BAXAO�AG�A?}A33A&�A�A�AoAoAVA
=A%AA��A��A�A�`A�AȴA�9A�A�A��A��A��A�uA�DA�Az�Ar�AffAbNAVAQ�AM�AE�A9XA1'A(�A �A�A�A�A�A{A{AbAJA1A  A��A�A�TA�
AƨA�FA��A�Al�AXAS�AC�A7LA+A"�A�AoA
=A%AAAA��A��A��A�A�yA�HA�A��A��A�!A��A�DAv�AbNAE�A(�AbA��A�mA��A�wA�A��A��A�hA�PA�7A�7A�A�A�A�A|�Ax�Ax�At�Ap�AhsA\)AO�A?}A+AoA��A�/A�RA�DAVA(�AJA�A�TA��AƨA�wA�^A�-A��A��A��A��A��A�hA�PA�A|�At�Al�AdZA\)AS�AK�AC�A?}A?}A7LA7LA/A/A/A/A/A+A+A&�A&�A"�A�A�A�A�AoA
=A
=AA ��A ��A �A �yA �`A �/A ��A ��A ȴA ĜA ��A �jA �RA �9A �A ��A ��A ��A ��A �DA �A z�A jA bNA ZA Q�A M�A M�A M�A M�A M�A I�A I�A E�A E�A E�A A�A A�A =qA =qA 5?A 1'A -A $�A $�A  �A �A �A {A {A {A bA bA JA JA 1A 1A A   @���@��@��@��@��m@��;@��;@��;@��;@��
@���@���@�ƨ@��w@��w@��F@��F@��F@��@��@��@��@��@��@���@��@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��P@��@�|�@�t�@�l�@�\)@�\)@�K�@�C�@�;d@�33@�33@�+@�+@�"�@�+@�+@�+@�+@�"�@�+@�"�@�"�@��@�o@�o@�
=@�
=@�
=@�@�@���@���@��@��@��H@��@���@���@��R@���@��\@�v�@�V@�5?@�{@��@���@��-@��h@�p�@�O�@�&�@��`@���@�I�@��@��@�l�@�+@��@���@���@���@��+@�v�@�v�@�ff@�V@�M�@�E�@�=q@�-@�$�@�J@��@��#@�@��^@��-@���@���@��h@��@�p�@�X@�7L@��@���@���@��u@�Q�@��@�l�@�ȴ@��@�`B@��@�D@�1@�\)@�ȴ@�n�@�5?@�@���@�h@�X@��@���@�r�@���G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�Bl�BjBs�Bt�Bv�Bl�BcTB_;B`BB`BB^5BcTB`BBdZBo�B|�B|�B~�B�B�B�B�B�B�B�B�+B�hB�7B��B��B�B�LB�}B�'B�B�3B�RB�9B�-B�RB�RB�wB�}B�qB�FB�FB�wB�dB�^B�dB�jB�RB�B�B�B��B��B��B��B��B�B��B��B�B�B�B�LB�B�B�-B�B��B�LB�3B�9B�3B�9B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B��B��B�B��B��B�B�B�B�^B�XB�RB�RB�XB�?B�FB�FB�FB�?B�dB�wB��BǮBǮB��B��B��B�
B�#B�NB�NB�NB�;B�5B�B��B�B��B��B��B�wB�LB�'B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�}BƨBǮBɺB��B�B�B�#B�)B�5B�;B�#B�)B�)B�B�B��B��B��B��B�B�B��BȴBƨBƨBŢBƨBƨBɺBɺBȴBŢB��B��BǮBÖBÖBŢBBBBĜBĜBĜBŢBÖBÖBĜBĜBĜBŢBÖB��BɺB��B��B��B��B��B��B��BȴB��B�'B�!B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B�hB�hB�bB�bB�bB�bB�bB�bB�\B�bB�bB�bB�oB�uB�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B�{B��B��B�{B�uB�uB�{B�uB�uB�uB�oB�oB�uB�uB�oB�oB�oB�uB�oB�hB�hB�bB�bB�hB�bB�oB�\B�\B�uB�oB�hB�hB�bB�bB�bB�bB�\B�\B�VB�VB�\B�VB�VB�VB�PB�PB�\B�DB�DB�DB�DB�DB�DB�JB�=B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�1B�=B�=B�7B�=B�DB�JB�DB�JB�JB�DB�7B�7B�DB�=B�7B�1B�+B�+B�1B�1B�1B�1B�1B�7B�7B�=B�DB�DB�DB�DB�DB�JB�JB�7B�1B�7B�1B�7B�DB�JB�JB�DB�DB�DB�DB�JB�JB�PB�PB�JB�JB�=B�DB�=B�=B�1B�1B�7B�7B�7B�7B�1B�1B�7B�=B�=B�B�bB�oB�oB�hB�bB�hB�bB�\B�hB�{B��B��B��B��B��B��B��B�hB�VB�DB�JB�JB�JB�PB�VB�PB�oB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�'B�!B�!B�B�'B�-B�!B�-B�-B�-B�-B�'B�'B�!B�!B�-B�!B�'B�-B�-B�B�B�B�'B�'B�9B�?B�9B�LB�RB�XB�XB�^B�^B�^B�XB�XB�RB�LB�XB�XB�XB�dB�^B�^B�XB�^B�XB�XB�XB�XB�^B�^B�^B�^B�^B�^B�^B�dB�dB�dB�dB�jB�dB�jB�dB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�dB�dB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�qB�qB�qB�qB�qB�qB�qB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�wB�}B�}B�}B�}B�}B�}B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BBBBBBBBBBBBBBBÖBÖBBBÖBÖBÖBÖBÖBBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBĜBÖBÖBÖBĜBĜBÖBÖBĜBĜBĜBĜBĜBÖBÖBÖBÖBÖBÖBĜBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBĜBÖBÖBÖBÖBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBÖBBÖBÖBÖBÖBk�Bl�Bm�Bl�Bl�Bm�Bk�BiyBk�Bl�BcTBgmBn�BgmBl�Bp�Bo�Bz�Bt�Br�Bt�B�Bp�Bm�Bq�Bp�Bp�Bm�Bn�Br�Bp�Bs�B{�B~�B�Bp�Bt�Bx�Bp�Bu�Bw�B{�B{�By�Bq�Bt�Bk�Bk�BiyBgmBiyBm�BgmBiyBffBdZBcTBaHB`BBaHBaHB_;B`BBaHBbNBaHB`BB^5B]/B]/B]/B^5B^5B^5B^5B_;BaHBbNBdZBgmBaHB]/B^5B_;B`BB`BB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB_;B]/B]/B\)B\)B\)B\)B]/B^5B^5B_;B`BBaHBbNBaHBbNBe`BgmBiyBhsBffBaHBaHBcTBcTB_;B_;B`BBaHB]/B^5B\)B\)B\)B]/B^5B^5B^5B`BB_;B`BBbNBdZBe`BffBffBgmBk�Bk�BjBiyBgmBgmBgmBffBe`BcTBbNB`BB_;B_;B`BB`BB`BB`BBaHBbNBbNBbNBdZBffBhsBiyBl�Bo�Br�Bv�Bx�B{�B�B�7B�+B�B�B� B� B}�B{�Bz�B{�B|�B}�B{�Bz�Bz�B{�Bz�B|�B}�B|�B}�B~�B}�B}�B}�B~�B}�B{�B{�B|�B}�B|�Bz�B{�B{�B{�B{�B{�B{�B{�B}�B}�B|�B}�B|�B}�B~�B}�B� B� B�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�+B�1B�+B�+B�%B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�B�%B�+B�%B�B�%B�%B�B�B�B�B�B�B�B�B�B�%B�%B�B�%B�%B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�%B�%B�%B�+B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�+B�+B�7B�=B�=B�7B�JB�JB�JB�DB�DB�JB�VB�oB��B��B��B��B��B��B�{B�uB�\B�DB�7B�7B�7B�7B�7B�7B�1B�1B�1B�7B�=B�DB�PB�PB�1B�+B�%B�+B�1B�1B�7B�DB�DB�PB�JB�VB�uB��B��B��B�B�jB�'B��B��B��B�9B�9B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Bn	Bn�By�B|B|�Bp
BeXB`�Bb Ba�B`3BeXBc}Bf�Bs�B~6B}�B�6B�B��B��B�sB�mB��B�B�B��B��B� B�-B�6B�GB̤B�gB�jB��B��B�PB��B�(B�0B�(B��B�B��B��B�B�dB��B��B�`B��B�HB��B��B�1B��B�lB��B��B��B�B�>B�BB�B�	B��B��B�MB�B��B�2B�YB�|B�B�9B�B�mB��B�2B��B��B�OB�B��B�B�2B�tB�`B�OB�hB�B�B�B�!B�tB��B�B�B�B��B�"B��B��B�PB�/B�zB�vB��B�,B�LB�?B��B�hB�3B�B�KB��B��B��B�RB��B�wB�iB�0B�iB�fB��B�jB��B��BƓB�wB�JBΌB��BڣB��B��B�[B��B��B� B�.BیBқB�(B�]B��B�B�B��B��B��B�B��B��B�B��B�pB�B��B��B�(B�B�!B��B��B��B�mB�gB��B��B�>B�B�B��B��B�B�=B�JB��B��B�B��B��B��B��B��B�B�wB�HB�KB��B��B��B�YB�.B�sB�|B��B��B�B��B��B�B�B��B��BثB�1B��B�
B�B�LB�CB�B�	B�B��BŴB��B� B�8B��B�B�SB�ZB�aB�kBëB��B�JBB�gB�dBĔBĪB��B�BèBÙBĦBĻB�B��B��B�AB��B�:B��BӑBѓBԉB��B�)B�hB��B��B�KB�|B� B� B��B��B�/B�IB�RB��B��B�XB��B�B�MB�VB��B�B�,B��B��B��B��B�fB�dB�oB�|B�oB�_B�_B�RB� B�B�MB�AB�yB��B��B��B�B� B��B�AB�B�EB��B��B��B��B��B��B��B��B�B�B�B�(B�!B��B��B��B�	B�tB�B��B��B��B��B��B��B�!B��B��B��B�B��B��B��B��B��B��B�(B�+B�B�B��B��B��B��B�UB�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�zB��B��B��B��B��B��B��B�LB�mB��B��B��B��B��B��B��B�oB��B��B�B�B��B�[B�6B��B��B��B��B��B��B��B��B��B��B��B��B�RB��B��B��B��B��B��B�/B��B��B�sB��B��B��B��B�[B�EB��B��B�yB�^B�uB��B��B��B�B�B�B��B��B��B�B�[B� B�}B�B��B��B�$B��B��B��B��B��B��B��B�pB��B��B��B��B��B��B��B��B��B�zB�uB��B��B��B��B��B�fB�dB�pB��B��B��B�B��B��B��B��B�iB�ZB��B��B�bB��B�QB�vB��B�lB��B�zB�5B��B��B�^B�DB�BB�DB�TB��B��B��B�DB�NB�LB�2B�LB�B�VB�.B�QB�9B�LB�*B�'B�%B�B�%B�2B�'B�%B�%B�B��B��B��B�DB�/B��B�B�ZB�JB�cB��B��B�PB�B�YB��B��B�|B�0B�B�-B�&B�$B�$B�B�B�B��B�"B�CB�QB�PB�.B�VB��B��B�JB�,B�0B��B��B�&B�XB�fB�PB�AB�B�B�>B�DB��B��B��B�{B�UB��B��B�QB�B�B�JB�HB��B�[B�B�B� B�RB��B��B� B��B��B�_B��B�zB�;B�tB��B�HB��B�|B�B��B��B�*B�UB�EB�}B�9B�8B�B�B�=B��B�0B�)B�[B�eB�WB�YB�B�FB��B��B��B��B�uB�XB�1B��B��B��B��B��B��B�2B�lB��B��B��B��B��B��B�B�vB�nB�IB�fB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�"B�sB��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�5B�'B�B��B��B�3B��B�=B�0B�oB�HB�]B�1B�B�B��B��B�DB��B�^B�8B��B��B��B��B�7B��B��B��B�,B�>B�hB��B��B��B��B��B�RB�.B�TB�B�uB��B��B�xB�eB�tB�YB�VB�VB�XB�QB�]B�^B�SB�lB�lB�bB�cB�WB�fB�uB�uB��B�rB�aB�rB�vB�xB�yB��B�xB�mB��B�xB�vB��B��B�vB��B��B��B�uB��B��B�wB�uB�lB�jB�yB��B�kB�vB�lB�lB��B�~B�rB��B��B�}B�}B��B��B�xB�vB��B��B�uB�uB��B��B�vB��B��B�|B��B��B��B��B�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B¿B­BBB��BBBBBBBB��B��B��BìB£B¹BýB��BðBâBÑB¢BïB��BñBàB×BÖBãBãBÕB×BïBÛBüBÿBàBäBìBÖBàBáBâBãBáBáBâBèBľBĜBæBðB��B��BĨBÛBðBĪB��B��B��BĳBãBØBãBÖBÖB��B��B��B��BýB��BÊBÖBíBäBðBÔBÖBüBĽBéBóBĳBÚBñB��B��BķBèB��BðBïBáBØB×BäBÕB×BàBÛBÙB��B��BöB��B��BèBðBýBçB��BýB×BÕBÕB×BåBòBâBÕBáBÕBäBôB��BðBâBÖBÕBËB×B×BáBÕBîBàBÕBòBâBåBâBâBðBòBðBïBâBÖBÖBÖBãB×BÐBBÖBïBãBØBk�Bl�Bm�Bl�Bl�Bm�Bk�BiyBk�Bl�BcTBgmBn�BgmBl�Bp�Bo�Bz�Bt�Br�Bt�B�Bp�Bm�Bq�Bp�Bp�Bm�Bn�Br�Bp�Bs�B{�B~�B�Bp�Bt�Bx�Bp�Bu�Bw�B{�B{�By�Bq�Bt�Bk�Bk�BiyBgmBiyBm�BgmBiyBffBdZBcTBaHB`BBaHBaHB_;B`BBaHBbNBaHB`BB^5B]/B]/B]/B^5B^5B^5B^5B_;BaHBbNBdZBgmBaHB]/B^5B_;B`BB`BB`BB`BBaHB`BB`BB`BB`BB`BB`BB`BB`BB`BBaHB_;B]/B]/B\)B\)B\)B\)B]/B^5B^5B_;B`BBaHBbNBaHBbNBe`BgmBiyBhsBffBaHBaHBcTBcTB_;B_;B`BBaHB]/B^5B\)B\)B\)B]/B^5B^5B^5B`BB_;B`BBbNBdZBe`BffBffBgmBk�Bk�BjBiyBgmBgmBgmBffBe`BcTBbNB`BB_;B_;B`BB`BB`BB`BBaHBbNBbNBbNBdZBffBhsBiyBl�Bo�Br�Bv�Bx�B{�B�B�7B�+B�B�B� B� B}�B{�Bz�B{�B|�B}�B{�Bz�Bz�B{�Bz�B|�B}�B|�B}�B~�B}�B}�B}�B~�B}�B{�B{�B|�B}�B|�Bz�B{�B{�B{�B{�B{�B{�B{�B}�B}�B|�B}�B|�B}�B~�B}�B� B� B�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�+B�1B�+B�+B�%B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�B�%B�+B�%B�B�%B�%B�B�B�B�B�B�B�B�B�B�%B�%B�B�%B�%B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�+B�%B�%B�%B�+B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�+B�+B�7B�=B�=B�7B�JB�JB�JB�DB�DB�JB�VB�oB��B��B��B��B��B��B�{B�uB�\B�DB�7B�7B�7B�7B�7B�7B�1B�1B�1B�7B�=B�DB�PB�PB�1B�+B�%B�+B�1B�1B�7B�DB�DB�PB�JB�VB�uB��B��B��B�B�jB�'B��B��B��B�9B�9B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          <%�<1j7<<�k<I�.<>�]<- 6<&��<%{@<&1�<%`�<&�n<&��<+̺<'uq<1|�<%<$P�<%�<$��<$�<$<<#�<#�<$�<$|d</y<+�<:F<iSz<u��<[>�<���<�YC<+��<8�B<,�u<3�6<$��<);-<*�<*$,<47a<(g?<%��<#�	<'�<(�<$��<$�`<(,�<*�F<'��<%�<%��<$g�<#��<$7�<#�N<#�<$H�<$0.<#��<#�U<#�N<#�<$��<&e<#�<#�4<$g�<%0<#�<#׎<#�M<$o�<$�<$i&<#�5<$?[<$��<$o�<$XX<#�<#ף<#�0<#�<#�<$�<#��<#�5<#��<#��<#ا<#�D<#ܯ<$v<$/%<#�o<#��<#�<#�<#ޫ<$F9<$\"<#��<#�8<$�<#��<$v<#�<#�<#�e<#�&<$�<#�8<#�<&/<#�4<$�<#��<$��<#�<#��<#��<#؄<$k�<$f�<$e.<$�J<$q@<$7�<$ʾ<#�W<$-<%4L<$��<$	<$<<$'<$�J<%͍<$F9<#�o<$�!<%p<&y<+�<-�z<(�<&!�<$��<$^�<$�t<$&<#�N<#�D<$6�<$!><$_�<$k<#�<#��<#ۮ<#�m<#�<#��<#ٛ<#�{<#�i<$�<$.<$��<$@|<#�H<#�N<#��<#��<#��<#��<#��<$�<#�]<#�$<#��<#��<#ۮ<#ٛ<#ا<#��<$3U<#��<%�d<)7,<+'�<$x+<$4e<%8j<%��<$<<$�<#��<#��<#׺<$#(<$7�<#��<$�<$<<$}�<$ K<#��<#�<#�<$�<%<+�<$��<#�<#�a<#�r<#�<#�<#�<$v<#ߜ<$3U<#�<(}�<%�M<$C�<#�c<#�&<$-<#�<#��<#ܯ<#�<<#ף<#�N<#�N<#�<#�<#�X<#��<#�	<#��<$(<'�B<$]h<#�<#�<$$<$3U<%͍<$��<,��<.9l<,��<&1�<'uq<$�<$=<$6�<#�<$�<$��<%@�<$�<#��<#�<#�(<#�<#�<$r<$'<'7�<%G<$k<$��<#�<#�U<#�e<#�<#�<#׎<#�<#׎<#�<#�<#��<#�U<#��<#ڑ<#�J<#�<$8�<$'<#�<#�H<#��<#�
<#�<$<<#�g<#�<#�<#�{<#ף<#�<#ۮ<#�0<#�<#�<#�^<#�<#�m<#�<#�$<#�<<#�o<#�&<$(<#�<#�<#�{<#�i<#ף<#�<$r�<%��<#�o<&�<&)�<$��<#�U<#�r<#ޫ<#��<#�r<#�<$�<$<<#�W<#�(<#�<#�D<#�l<#��<#��<#�o<#�*<#�U<#�<#��<#�&<#�<#�<#�^<#ۮ<#׎<#׎<#�i<#؄<#׺<#�<#�X<#�<#�i<#�X<#��<#ף<#׺<#�<#��<#�e<#�e<#ا<#��<#�*<#ا<#��<#׺<#�<#�C<#ޫ<$	�<$'<#��<#�<#ף<#�{<#�{<#��<#�<#�e<#�$<#��<#�N<#�N<#׺<$.<$v<#ޫ<#�J<#��<#׎<#�<#�<#�<#�<#�<#׎<#��<#�l<$%<#�D<#ٛ<$2G<#�o<#�<#�^<$<<#ߜ<#�D<#�&<#�l<#ۮ<#�<#��<$.<%m�<#ף<$}�<#�8<#�<#�e<#�*<#��<#�<#�<#�<#�!<#��<#�<#��<#��<$(<#�<$\"<#��<#�*<#�*<$'<#��<#�]<#�{<#��<#��<#�4<#��<#�<#��<#�<#�8<#�<#�<#�<#��<#�<#��<#�i<#�&<#�<#�<#�<#��<#�<#�<#�<#�<<#ޫ<#�a<#�(<#�<#��<#�<#��<#ߜ<#�0<#�<<#�J<#�<#�&<#�8<#�<#�*<#�<#؄<#�<#��<#�D<$ �<#�<#�<#�
<#�<#�
<#��<#�<$Z<#�<#�$<#ٛ<#�D<#�<#�D<#�<#�8<#�<#�r<#ף<#�]<#�<#�<#�
<#�i<#�
<#׎<#�<#�
<#�
<#�&<#��<#�H<#��<#�0<#ף<#��<#�<#؄<#�
<#��<#�&<#��<#�{<#�D<#ڑ<#�W<#�<#�<#�<#؄<#�<#�i<#׎<#׎<#�<#��<#��<#�<#�D<#�<#׎<#�{<#؄<#�<#�"<$ K<#�$<#�<#�0<#��<#�<#��<#ף<#�o<#�{<#�<#��<#�<#�{<#�&<#�<#��<#��<#�^<#��<#��<#�!<#�C<#�o<#��<#�$<#��<#�<#��<#�*<#�]<#ا<#�c<$�<%S�<#�W<#�D<#�<#�I<#��<#�<#ۮ<$z�<$'<#��<#�I<#��<#�<<#�D<#�E<$.<&L�<$k�<#ۮ<#�i<#�<#�U<#�<#�$<#��<$�h<#��<#�*<#؄<#�l<#��<#�Q<#��<#�*<#�J<#�E<#ا<#��<#ߜ<#�W<#�<#�X<#�C<#؄<#�4<#��<#�	<#�<#�<#ڑ<#ڑ<#��<#�D<#�D<#��<#�U<#�<$�<#�<#׎<#�<#�I<#�<#ף<#�<#�<#׎<#�i<#��<#�
<#��<#ڑ<#ۮ<#��<#�i<#�<<#�{<#��<$<<$c�<$�<#�]<#��<#�D<#��<#�<#׺<#��<#��<#�<#׎<#�<#ا<#�c<#�<#׺<#�D<#ۮ<#�<#�<%�<$��<#ܯ<#ף<#�&<#�<#�<#�N<#�&<#�8<#��<#�<#�U<#�D<#��<#�X<#�X<#��<#�a<$	�<#ٛ<$f<#�^<#ٛ<#ۮ<$<<#�<#�a<#�<#��<$<<#�<#�r<#�<#��<#ޫ<#��<#�<#��<#�<#�
<#��<#�<#�<#ٛ<#�^<#��<#�<#׎<#؄<#�<#�<#�<#�
<#׎<#�<#�
<#�i<#ף<#ף<#�<#�<#׎<#�<#��<#�i<#ڑ<#�<<#�<#�<<#�{<#ף<#׺<#�<#ף<#�<#��<#ף<#�{<#�8<#�<#�{<#�<#�r<#ڑ<#��<#ڑ<#ڑ<#׎<#�i<#�<#�
<#׺<#�l<#�<#�{<#�<#�<#�<#׎<#�<#ا<#�<#�{<#�{<#�<#��<#�<#�<#ף<#׎<#�<#�<#׎<#��<#�<#�i<#��<#�<#׎<#�D<#�o<#׎<#�<#ٛ<#�e<#��<#׺<#�<#�<#ܯ<#�^<#׺<#ף<#�$<#�<#�<#�<#ا<#�r<#�8<#ף<#�<#ף<#�o<#�o<#�<#��<#׺<#��<#��<#�<#�<#�<#��<#�<#�<#ף<#ߜ<#�l<#��<#؄<#�C<#�l<#ۮ<#�<#�<#�{<#�<#�$<#��<#�J<#�D<#�X<#�<#�
<#׎<#׎<#�<#�<#��<#�<#�r<#�+<#�X<#ף<#؄<#�
<#�X<#�i<#�{<#׎<#�i<#�i<#�{<#�<#ڑ<#�
<#��<#�<#�l<#��<#�{<#�<#�<#ף<#�J<#��<#�l<#ا<#׎<#�<#׎<#�
<#�
<#�8<#�<#��<#ߜ<#ۮ<#�<#�{<#�
<#ا<#ף<#�<#�<#�
<#�r<#�]<#�$<#ٛ<#ا<#�<#�D<#��<#�!<#�D<#�<#��<#�<#��<#�i<#�<#�<#ף<#�<#�<#�X<#�<#�<#��<#�U<#�*<#�<#�e<#�<#�<#ۮ<#��<#�J<#ۮ<#�<#�<#�<#�<#׺<#�o<#�{<#�<#�i<#�<#ף<#��<#�4<#�<#�{<#�
<#�<#�i<#�<#�<#�i<#�<#��<#�X<#�<#�o<#�{<#׺<#�{<#�{<#�<#�o<#�<#��<#�{<#�
<#�
<#�
<#׎<#�<#�&<#�i<#�
<#��<#׎<#�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.004(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.004(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202208010000002022080100000020220801000000202208010000002022080100000020220801000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2020070502004420200705020044QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2020070502004420200705020044QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2021011500000020210115000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2022080100000020220801000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2021V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2021V03; CTD_for_DMQC_2021V02                     2022080100000020220801000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                