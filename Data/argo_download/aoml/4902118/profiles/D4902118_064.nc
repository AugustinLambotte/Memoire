CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       $Woods Hole Oceanographic Institution   source        
Argo float     history       92018-11-06T01:28:33Z creation; 2022-08-04T19:12:10Z DMQC;      
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
resolution        =���   axis      Z        p  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  \   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  c�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �h   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  �D   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p  ��   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �$   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p  �    TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �p   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p  �L   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p �   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 � 8,   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p @   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 � _x   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     p gT   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   �$   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   �$   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �$   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �$   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   �x   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   �4   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �<   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        �\   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        �d   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       �l   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181106012833  20220804151210  4902118 4902118 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               @   @AA  AOAO6731                            6731                            2C  2C  DD  S2A                             S2A                             7360                            7360                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�P�N�^@�P�N�^11  @�P��W@@�P��W@@OY����@OY�����B���l��B���l�11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  AA  AA  ?�z�@�\@B�\@�  @��R@��R@�  A   A\)A#33A?\)A_\)A\)A�  A�  A��A��A�  A�  A�  A��B  B(�B(�B�
B(  B0(�B8  B@  BG33BO�
BX  B`  Bh(�BpQ�Bxz�B�  B��
B��
B��
B�B�  B�  B��
B��
B�B�B��
B��B��B��
B��
B�B��
B�B�  B�  B�  B�(�B�  B��B��B�  B�{B�  B��
B�{B�=qC {C  C
=C
=C  C
  C��C��C�C  C{C
=C
=C  C��C  C {C"
=C$
=C%��C(  C*
=C+�C.  C0
=C2
=C4{C6�C8
=C9��C<{C>  C?�CB  CD  CE��CH
=CJ  CK��CM��CO�CR
=CT  CV  CX
=CZ{C\
=C^  C`{Cb{Cd  Cf  Ch{Cj
=Cl  Cn
=Co��Cq�HCs�Cu��Cx{Cz  C{��C~  C�C���C��C��C���C�C�  C���C���C�C�\C�C�  C�  C�C�
=C���C��C���C���C���C���C�  C�C�\C�C�  C�
=C�  C���C�C�\C�C�  C���C���C�C�  C���C�C�  C���C�C�  C���C�
=C�C�  C���C���C�C�  C���C���C���C�
=C�  C���C���C�  C�C���C�  C�
=C�C���C���C�
=C�C�  C�  C�  C�  C���C���C���C�  C�  C�  C�  C���C���C�  C�  C�  C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C�  C�  C�  C���C���C���C���C�  C�  C�  C���C���C�  C���C���C���C���C���C���C�  C�C�C�  C�C�
=C�C�  C�  C���D }qD �qD}qD  D� D  D� D  D��D�D�DD� D��D}qD�qD� D�qD	}qD
�D
��D
��Dz�D�qD� D�D� D�qD}qD�D��D  D}qD  D� D�D�D�D� D��D� D  D��D�D� D��D}qD�D� D�qD}qD  D�D�D� D�D��D�qD��DD}qD  D�D D ��D!D!�D"D"�D#�D#� D#�qD$z�D$��D%z�D%��D&z�D&��D'z�D'�qD(}qD(��D)z�D)�qD*� D+  D+��D,�D,��D-D-�D.�D.� D.�RD/z�D/��D0}qD1  D1��D2D2� D2��D3z�D4  D4� D4�RD5}qD6�D6� D7  D7��D8  D8}qD8�qD9xRD9�RD:xRD:��D;z�D;�qD<z�D<��D=� D>D>�D?D?�D@D@�DA�DA� DA�RDB� DCDC�DDDD�DEDE�DFDF�DGDG�DH�DH��DI�DI�DJ�DJ�DK�DK��DL�DLz�DL��DM}qDM�qDN}qDN��DO}qDP�DP��DQDQ� DQ�qDR��DS  DS� DT�DT� DT��DU� DVDV� DV�RDW}qDX  DX�DY  DYxRDY��DZ� D[D[� D[�RD\}qD]  D]��D^D^��D_D_� D_��D`}qDa  Da� Db�Db�DcDc�DdDd��Dd�qDez�De�RDfxRDf��Dgz�Dh  Dh��Di  Di}qDi��Dj� DkDk�DlDl��Dm  Dm� Dn  Dn� Do  Do� Do�qDpz�Dp��Dqz�Dq�qDrz�Dr��Ds}qDs�qDt� Du�Du��DvDv�Dw�Dw��Dw�qDxxRDx�qDy}qDy��Dz}qD{  D{� D{�qD|}qD|�qD}� D~  D~}qD  D�D�  D�>�D��HD��HD�  D�@ D�� D�� D�HD�AHD�� D��)D��qD�>�D�� D�� D���D�@ D��HD���D��qD�@ D���D��HD�  D�B�D��HD���D���D�>�D�}qD���D�  D�@ D�~�D�� D�HD�@ D�� D�� D���D�B�D���D���D���D�AHD��HD��HD�HD�>�D�� D�� D�  D�>�D�~�D���D�HD�B�D�� D���D��qD�@ D���D��HD�HD�@ D�� D�� D�  D�AHD���D��HD�  D�@ D��HD���D��)D�<)D�}qD�� D�  D�<)D�}qD��qD���D�>�D�~�D�� D�HD�>�D�~�D�� D��qD�>�D�� D��HD�  D�@ D�~�D�� D���D�>�D�~�D�� D�  D�=qD�~�D���D�  D�AHD�~�D���D�  D�AHD��HD���D���D�@ D�� D��HD���D�@ D�� D��HD��D�AHD�~�D��qD���D�AHD�� D��HD��D�AHD�� D���D��)D�>�D��HD��HD�  D�>�D�~�D��HD��D�B�D��HD�� D���D�=qD�~�D�� D��qD�>�D��HD�� D���D�>�D�}qD��)D�  D�B�D��HD���D��qD�@ D���D��HD�  D�>�D�}qD�� D�HD�@ D�}qD��qD�  D�C�D���D�D�HD�AHD�� D���D��qD�@ D���D�D�HD�@ D�� D���D���D�=qD�� D�D�HD�@ D�~�D���D�  D�@ D��HD��HD�  D�=qD��HD���D��D�@ D�� D��HD�HD�B�D�~�D��)D��qD�=qD�}qD��qD��qD�=qD�|)D��qD��qD�>�D�}qD��qD���D�>�D�~�D���D�  D�B�D���D�� D��qD�=qD�~�D�� D�  D�AHD���D�� D��qD�AHD�~�D��)D���D�AHD��HD��HD���D�>�D�~�D�� D�HD�@ D�}qD��qD���D�@ D��HD�� D���D�AHD D½qD�  D�>�D�}qD�� D��D�AHD�~�D��HD���D�<)D�~�D�� D���D�<)D�~�D�D�HD�AHDǂ�D�D��D�B�DȂ�D���D�  D�>�Dɀ D��HD�HD�AHDʀ D�� D�  D�AHD˂�D�� D��qD�=qD�|)D�� D�  D�=qD�}qD;�D�  D�B�D΀ DνqD���D�>�DρHD���D��D�@ D�|)DнqD��qD�@ DсHD�D��D�B�D҂�D���D�  D�=qD�~�D��HD�  D�>�DԁHD�� D���D�>�D�}qD�� D���D�=qDցHD�� D�  D�@ D�~�D׾�D�  D�@ D�}qDؽqD��qD�>�Dـ D��HD�HD�AHDځHD��HD���D�@ Dۂ�D��HD�  D�>�D܁HD���D�HD�>�D݀ D�� D�HD�B�Dހ D��HD��D�AHD�~�D߽qD���D�AHD�~�D�qD���D�@ D�}qDᾸD�  D�>�D�~�D�qD�HD�B�D�HD�� D�  D�>�D� D��HD�  D�>�D� D��HD�HD�B�D� D�� D�HD�@ D�~�D�� D�  D�@ D�HD��HD�HD�>�D� D��HD�HD�@ D�}qD�� D�HD�@ D�~�D�qD�  D�B�D�HD�� D�HD�B�D� D�� D���D�@ D�HD�� D�  D�@ D� D��HD�  D�@ D�� D��HD�HD�AHD�D�� D���D�@ D�HD�� D��qD�>�D�~�D�D�HD�@ D�~�D��HD�HD�>�D�~�D�� D�  D�@ D�� D�� D�HD�AHD�� D���D�  D�AHD�~�D�� D�HD�@ D��HD��HD�  D�E?k�?u?�=q?��R?�{?�Q�?��?�ff?�?��H@�@\)@�@(�@&ff@.{@5@:�H@@  @J=q@Q�@W
=@\(�@fff@n{@u@}p�@�G�@��@���@�{@���@�@�Q�@�(�@�G�@��@���@���@���@�z�@�Q�@�(�@�  @��
@Ǯ@˅@У�@�z�@�Q�@�p�@�  @��
@�@�@�\)@�33@�Q�@�(�A   A�\Az�AffAQ�A	��A(�A{A��A�\Az�A
=AQ�A�HA(�A{A ��A"�\A$z�A&ffA(Q�A*=qA,(�A.{A0  A1�A3�
A6ffA8Q�A:�HA<��A>�RAAG�AB�\ADz�AG
=AH��AK�AL��AN�RAP��AQ�AS�
AUAW�AY��A[�A]p�A_\)AaG�Ac33Ae�Ag
=Ah��Aj�HAl��An�RAp��Ar�\Atz�AvffAx��Az�HA|��A~�RA�Q�A�G�A�=qA�33A�(�A��A�{A�\)A�  A�G�A��A��\A��A�z�A��A�{A�
=A��A���A���A�=qA��A�(�A���A�p�A�ffA�
=A�  A���A���A��\A�33A�(�A���A�A�ffA�\)A�Q�A���A��A��HA��
A���A�p�A�ffA�
=A�  A���A���A��\A�33A�(�A��A�{A�
=A�  A���A���A�=qA�33A�(�A���A�A�ffA�\)A�  A���A���A��\A�33A�(�A���A�A�ffA�\)A�  A���A���A\AÅA�(�A��A�{AƸRAǮA�Q�A�G�A�=qA��HA��
A���A�A�ffA�\)A�Q�A�G�A�=qA�33A��
A���A�AָRA�\)A�Q�A�G�A�=qA�33A�(�A�p�A�{A�\)A�Q�A���A�=qA�33A�(�A��A�{A�RA�  A���A��A��HA��
A���A�A�ffA�A�Q�A�G�A�=qA�33A�z�A�p�A�ffA�\)A���A���A��\A�33A�(�A��A�A�
=B   B ��B�BB=qB�RB\)B�
BQ�B��B�B��B=qB�RB33B�
BQ�B��B	p�B	�B
ffB
�HB\)B�
BQ�B��BG�B�BffB
=B�B  Bz�B�B��B{BffB�HB\)B�
BQ�B��Bp�B�B�\B
=B�B  Bz�B��B��B{BffB
=B\)B�
BQ�B��Bp�B�B�\B
=B�B   B z�B ��B!G�B!�B"ffB"�HB#�B$  B$z�B$��B%p�B%�B&ffB&�HB'\)B(  B(��B)�B)��B*{B*�\B+
=B+�B,  B,z�B,��B-��B.{B.�\B/
=B/�B0(�B0z�B0��B1p�B1�B2ffB2�RB3\)B3�
B4Q�B4��B5p�B5�B6ffB6�HB7\)B7�
B8Q�B8��B9G�B9B:=qB:�RB;33B;�
B<Q�B<��B=p�B=�B>=qB>�RB?33B?�B@Q�B@��BAG�BABBffBB�HBC\)BC�
BDQ�BD��BEG�BEBFffBF�HBG33G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ?�z�@�\@B�\@�  @��R@��R@�  A   A\)A#33A?\)A_\)A\)A�  A�  A��A��A�  A�  A�  A��B  B(�B(�B�
B(  B0(�B8  B@  BG33BO�
BX  B`  Bh(�BpQ�Bxz�B�  B��
B��
B��
B�B�  B�  B��
B��
B�B�B��
B��B��B��
B��
B�B��
B�B�  B�  B�  B�(�B�  B��B��B�  B�{B�  B��
B�{B�=qC {C  C
=C
=C  C
  C��C��C�C  C{C
=C
=C  C��C  C {C"
=C$
=C%��C(  C*
=C+�C.  C0
=C2
=C4{C6�C8
=C9��C<{C>  C?�CB  CD  CE��CH
=CJ  CK��CM��CO�CR
=CT  CV  CX
=CZ{C\
=C^  C`{Cb{Cd  Cf  Ch{Cj
=Cl  Cn
=Co��Cq�HCs�Cu��Cx{Cz  C{��C~  C�C���C��C��C���C�C�  C���C���C�C�\C�C�  C�  C�C�
=C���C��C���C���C���C���C�  C�C�\C�C�  C�
=C�  C���C�C�\C�C�  C���C���C�C�  C���C�C�  C���C�C�  C���C�
=C�C�  C���C���C�C�  C���C���C���C�
=C�  C���C���C�  C�C���C�  C�
=C�C���C���C�
=C�C�  C�  C�  C�  C���C���C���C�  C�  C�  C�  C���C���C�  C�  C�  C���C���C���C���C���C���C���C���C���C���C���C���C���C���C���C�  C�  C�  C���C���C���C���C�  C�  C�  C���C���C�  C���C���C���C���C���C���C�  C�C�C�  C�C�
=C�C�  C�  C���D }qD �qD}qD  D� D  D� D  D��D�D�DD� D��D}qD�qD� D�qD	}qD
�D
��D
��Dz�D�qD� D�D� D�qD}qD�D��D  D}qD  D� D�D�D�D� D��D� D  D��D�D� D��D}qD�D� D�qD}qD  D�D�D� D�D��D�qD��DD}qD  D�D D ��D!D!�D"D"�D#�D#� D#�qD$z�D$��D%z�D%��D&z�D&��D'z�D'�qD(}qD(��D)z�D)�qD*� D+  D+��D,�D,��D-D-�D.�D.� D.�RD/z�D/��D0}qD1  D1��D2D2� D2��D3z�D4  D4� D4�RD5}qD6�D6� D7  D7��D8  D8}qD8�qD9xRD9�RD:xRD:��D;z�D;�qD<z�D<��D=� D>D>�D?D?�D@D@�DA�DA� DA�RDB� DCDC�DDDD�DEDE�DFDF�DGDG�DH�DH��DI�DI�DJ�DJ�DK�DK��DL�DLz�DL��DM}qDM�qDN}qDN��DO}qDP�DP��DQDQ� DQ�qDR��DS  DS� DT�DT� DT��DU� DVDV� DV�RDW}qDX  DX�DY  DYxRDY��DZ� D[D[� D[�RD\}qD]  D]��D^D^��D_D_� D_��D`}qDa  Da� Db�Db�DcDc�DdDd��Dd�qDez�De�RDfxRDf��Dgz�Dh  Dh��Di  Di}qDi��Dj� DkDk�DlDl��Dm  Dm� Dn  Dn� Do  Do� Do�qDpz�Dp��Dqz�Dq�qDrz�Dr��Ds}qDs�qDt� Du�Du��DvDv�Dw�Dw��Dw�qDxxRDx�qDy}qDy��Dz}qD{  D{� D{�qD|}qD|�qD}� D~  D~}qD  D�D�  D�>�D��HD��HD�  D�@ D�� D�� D�HD�AHD�� D��)D��qD�>�D�� D�� D���D�@ D��HD���D��qD�@ D���D��HD�  D�B�D��HD���D���D�>�D�}qD���D�  D�@ D�~�D�� D�HD�@ D�� D�� D���D�B�D���D���D���D�AHD��HD��HD�HD�>�D�� D�� D�  D�>�D�~�D���D�HD�B�D�� D���D��qD�@ D���D��HD�HD�@ D�� D�� D�  D�AHD���D��HD�  D�@ D��HD���D��)D�<)D�}qD�� D�  D�<)D�}qD��qD���D�>�D�~�D�� D�HD�>�D�~�D�� D��qD�>�D�� D��HD�  D�@ D�~�D�� D���D�>�D�~�D�� D�  D�=qD�~�D���D�  D�AHD�~�D���D�  D�AHD��HD���D���D�@ D�� D��HD���D�@ D�� D��HD��D�AHD�~�D��qD���D�AHD�� D��HD��D�AHD�� D���D��)D�>�D��HD��HD�  D�>�D�~�D��HD��D�B�D��HD�� D���D�=qD�~�D�� D��qD�>�D��HD�� D���D�>�D�}qD��)D�  D�B�D��HD���D��qD�@ D���D��HD�  D�>�D�}qD�� D�HD�@ D�}qD��qD�  D�C�D���D�D�HD�AHD�� D���D��qD�@ D���D�D�HD�@ D�� D���D���D�=qD�� D�D�HD�@ D�~�D���D�  D�@ D��HD��HD�  D�=qD��HD���D��D�@ D�� D��HD�HD�B�D�~�D��)D��qD�=qD�}qD��qD��qD�=qD�|)D��qD��qD�>�D�}qD��qD���D�>�D�~�D���D�  D�B�D���D�� D��qD�=qD�~�D�� D�  D�AHD���D�� D��qD�AHD�~�D��)D���D�AHD��HD��HD���D�>�D�~�D�� D�HD�@ D�}qD��qD���D�@ D��HD�� D���D�AHD D½qD�  D�>�D�}qD�� D��D�AHD�~�D��HD���D�<)D�~�D�� D���D�<)D�~�D�D�HD�AHDǂ�D�D��D�B�DȂ�D���D�  D�>�Dɀ D��HD�HD�AHDʀ D�� D�  D�AHD˂�D�� D��qD�=qD�|)D�� D�  D�=qD�}qD;�D�  D�B�D΀ DνqD���D�>�DρHD���D��D�@ D�|)DнqD��qD�@ DсHD�D��D�B�D҂�D���D�  D�=qD�~�D��HD�  D�>�DԁHD�� D���D�>�D�}qD�� D���D�=qDցHD�� D�  D�@ D�~�D׾�D�  D�@ D�}qDؽqD��qD�>�Dـ D��HD�HD�AHDځHD��HD���D�@ Dۂ�D��HD�  D�>�D܁HD���D�HD�>�D݀ D�� D�HD�B�Dހ D��HD��D�AHD�~�D߽qD���D�AHD�~�D�qD���D�@ D�}qDᾸD�  D�>�D�~�D�qD�HD�B�D�HD�� D�  D�>�D� D��HD�  D�>�D� D��HD�HD�B�D� D�� D�HD�@ D�~�D�� D�  D�@ D�HD��HD�HD�>�D� D��HD�HD�@ D�}qD�� D�HD�@ D�~�D�qD�  D�B�D�HD�� D�HD�B�D� D�� D���D�@ D�HD�� D�  D�@ D� D��HD�  D�@ D�� D��HD�HD�AHD�D�� D���D�@ D�HD�� D��qD�>�D�~�D�D�HD�@ D�~�D��HD�HD�>�D�~�D�� D�  D�@ D�� D�� D�HD�AHD�� D���D�  D�AHD�~�D�� D�HD�@ D��HD��HD�  D�E?k�?u?�=q?��R?�{?�Q�?��?�ff?�?��H@�@\)@�@(�@&ff@.{@5@:�H@@  @J=q@Q�@W
=@\(�@fff@n{@u@}p�@�G�@��@���@�{@���@�@�Q�@�(�@�G�@��@���@���@���@�z�@�Q�@�(�@�  @��
@Ǯ@˅@У�@�z�@�Q�@�p�@�  @��
@�@�@�\)@�33@�Q�@�(�A   A�\Az�AffAQ�A	��A(�A{A��A�\Az�A
=AQ�A�HA(�A{A ��A"�\A$z�A&ffA(Q�A*=qA,(�A.{A0  A1�A3�
A6ffA8Q�A:�HA<��A>�RAAG�AB�\ADz�AG
=AH��AK�AL��AN�RAP��AQ�AS�
AUAW�AY��A[�A]p�A_\)AaG�Ac33Ae�Ag
=Ah��Aj�HAl��An�RAp��Ar�\Atz�AvffAx��Az�HA|��A~�RA�Q�A�G�A�=qA�33A�(�A��A�{A�\)A�  A�G�A��A��\A��A�z�A��A�{A�
=A��A���A���A�=qA��A�(�A���A�p�A�ffA�
=A�  A���A���A��\A�33A�(�A���A�A�ffA�\)A�Q�A���A��A��HA��
A���A�p�A�ffA�
=A�  A���A���A��\A�33A�(�A��A�{A�
=A�  A���A���A�=qA�33A�(�A���A�A�ffA�\)A�  A���A���A��\A�33A�(�A���A�A�ffA�\)A�  A���A���A\AÅA�(�A��A�{AƸRAǮA�Q�A�G�A�=qA��HA��
A���A�A�ffA�\)A�Q�A�G�A�=qA�33A��
A���A�AָRA�\)A�Q�A�G�A�=qA�33A�(�A�p�A�{A�\)A�Q�A���A�=qA�33A�(�A��A�{A�RA�  A���A��A��HA��
A���A�A�ffA�A�Q�A�G�A�=qA�33A�z�A�p�A�ffA�\)A���A���A��\A�33A�(�A��A�A�
=B   B ��B�BB=qB�RB\)B�
BQ�B��B�B��B=qB�RB33B�
BQ�B��B	p�B	�B
ffB
�HB\)B�
BQ�B��BG�B�BffB
=B�B  Bz�B�B��B{BffB�HB\)B�
BQ�B��Bp�B�B�\B
=B�B  Bz�B��B��B{BffB
=B\)B�
BQ�B��Bp�B�B�\B
=B�B   B z�B ��B!G�B!�B"ffB"�HB#�B$  B$z�B$��B%p�B%�B&ffB&�HB'\)B(  B(��B)�B)��B*{B*�\B+
=B+�B,  B,z�B,��B-��B.{B.�\B/
=B/�B0(�B0z�B0��B1p�B1�B2ffB2�RB3\)B3�
B4Q�B4��B5p�B5�B6ffB6�HB7\)B7�
B8Q�B8��B9G�B9B:=qB:�RB;33B;�
B<Q�B<��B=p�B=�B>=qB>�RB?33B?�B@Q�B@��BAG�BABBffBB�HBC\)BC�
BDQ�BD��BEG�BEBFffBF�HBG33G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@���@�@�@�@�@���@��@��@��@��y@��y@��@���@�@���@�@�
=@�
=@�o@�"�@�o@��y@��y@�@��H@��y@���@���@���@��\@��\@�~�@�~�@��+@�v�@�n�@�n�@�^5@�V@�V@�V@�M�@�^5@�v�@���@��+@�n�@�n�@�n�@�n�@�E�@��T@���@�@���@���@���@���@���@���@���@���@���@��-@��h@�p�@�X@�`B@�O�@�V@���@���@���@���@���@��u@��u@��u@��u@��D@��@�r�@�bN@�Q�@�1'@��@�1@�b@�@��@�;@�;@��@�;@�w@�w@�w@��@�P@|�@|�@l�@l�@l�@\)@\)@K�@K�@K�@+@
=@
=@~��@~�y@~�@~ȴ@~��@~��@~��@~��@~��@~��@~��@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~ȴ@~ȴ@~ȴ@~ȴ@~�@~�@~�@~�@~�@~�@~�@~�@~�y@~�y@~�y@~�y@~�y@~��@~��@~��@~��@~��@
=@
=@
=@
=@
=@
=@
=@
=@
=@
=@�@
=@
=@�@�@�@+@�@
=@
=@
=@
=@~�y@~��@�@
=@~��@~@}@}��@}�h@}�h@}p�@}p�@}�@}�@}�@}�@}�@}�h@}�-@}��@~@~5?@~�+@~��@~�R@~ȴ@
=@
=@~��@~�@~�@~�@~ȴ@~ȴ@~ȴ@~�+@}�@}��@}��@}�h@}�-@}��@}�@}?}@|�@|(�@|�@{��@|9X@{�m@{�@{�@{dZ@y��@x�9@x1'@w�@w�w@w�P@w|�@w\)@w;d@w+@v�@v��@vff@v$�@u��@u�h@u�@up�@uO�@u?}@u/@u/@uV@t�@t�j@tj@tZ@t9X@t(�@t�@s��@s��@s�m@s�
@s�
@sƨ@s�F@s��@s��@sdZ@s33@r�\@r-@rJ@q�@q�@q�^@q�7@q&�@pQ�@p �@pb@o��@o�;@o��@o�@o��@o��@o�@o��@o\)@oK�@o\)@oK�@o;d@o�@o�@o+@o
=@nȴ@n�R@n��@n��@nff@nV@n@m�-@m�-@m��@m��@mp�@mO�@m`B@mO�@mV@l�/@l�/@l�j@l�j@l�j@l��@l�j@l��@lz�@lj@lI�@l9X@l�@k��@k�m@k�
@kƨ@k�F@k��@k�@kt�@kt�@kS�@k"�@ko@j�H@j��@j��@j�\@j�\@jn�@jM�@j=q@j-@j�@jJ@i��@i��@i�@i�#@i��@i��@i�7@ix�@ihs@iX@iG�@i&�@h�`@h��@h�@h�@hr�@hr�@hr�@hr�@hr�@hr�@h�@hr�@hbN@hbN@hr�@hbN@hA�@hA�@h1'@h �@g�@g��@g�w@g�@g�@g�P@g|�@g�P@g|�@gl�@g\)@g\)@gK�@g+@g�@g
=@f�R@f��@f��@f�+@f�+@fv�@f�+@fv�@fv�@fv�@fv�@fv�@fff@fV@fV@fff@fV@fE�@f$�@f@f@e�@e�@e�@e�T@e�@e�@e�@f@e�@e�T@e�-@e��@e�h@e�@e�@e?}@eV@eV@d�/@d��@d�/@d��@d�j@d�D@dz�@dz�@dz�@dj@dZ@dI�@d9X@d�@d�@d�@d�@d�@d�@d1@d1@c��@c�m@c�
@c�F@c��@ct�@cS�@cC�@c33@co@co@c@c@b�H@b�!@b��@b�\@b�\@b�\@b~�@b=q@b�@a��@a��@ahs@aX@`��@`Ĝ@`�9@`��@`�u@`�@`�@`�@`bN@`Q�@`A�@`A�@`1'@` �@`b@_�@_��@_�@_�@_�@_�@_�@_�P@_|�@_�P@_|�@_�P@_�P@_|�@_l�@_l�@_l�@_l�@_l�@_l�@_l�@_K�@_K�@_K�@_K�@_K�@_�@^�y@^�@^�@^ȴ@^��@^��@^��@^�+@^v�@^v�@^ff@^V@^V@^E�@^V@^V@^E�@^{@^@]�@]�T@]�T@]�T@]�@]�T@]�T@]�T@]�@]�@]�@]��@]@]�-@]�h@]p�@]p�@]`B@]?}@]?}@]?}@]?}@]?}@]/@]/@]?}@]?}@]?}@]?}@]�@]/@]/@]/@]/@]?}@]?}@]/@]?}@]?}@]/@]V@\��@]V@]V@]V@\��@\�@\�@\�@\�@\�/@\�/@\�/@\�/@\�@\�@\�@\�@\�@\�/@\��@\�j@\�j@\�j@\�D@\�D@\z�@\j@\j@\z�@\j@\I�@\9X@\9X@\9X@\(�@\(�@\(�@\(�@\(�@\(�@\�@\1@[��@[��@[�
@[�
@[�
@[ƨ@[�
@[�
@[�F@[��@[��@[��@[��@[��@[��@[��@[�@[t�@[dZ@[S�@[C�@[C�@[S�@[33@[33@["�@[o@[@[@Z�@Z�@Z�@Z�@Z�H@Z��@Z��@Z��@Z��@Z�H@Z�H@Z��@Z��@Z��@Z��@Z�H@Z�H@Z�H@Z�H@Z�@Z�@Z�@[@[33@[C�@[C�@[S�@[dZ@[dZ@[dZ@[dZ@[t�@[�@[��@[ƨ@[�
@[�m@[�m@[��@\1@\9X@\Z@\��@\�/@\�@]V@]V@]V@]/@]�h@]@^{@^ff@^v�@^�@_
=@_+@_l�@_�@_�;@` �@`A�@`A�@`bN@`�@`��@`��@a�@a�@a&�@aG�@aX@a�7@a�^@b�@b~�@b�\@b��@b��@c@c"�@cS�@cdZ@ct�@c�@c�@c��@c��@cƨ@c�m@dZ@dz�@d��@eV@e?}@e`B@e�@e@f@f{@fE�@fv�@f�+@f��@f�R@fȴ@f�@fȴ@fȴ@fȴ@fȴ@f�y@g
=@g�@g;d@gK�@g\)@g|�@g�P@h  @h�u@hĜ@h�9@h�9@h��@h�9@h��@h�`@h�`@i%@i7L@ihs@i��@ix�@ix�@i7L@h��@h�9@h��@h�`@h��@h��@i%@i�7@i�@i�@i�@i�#@i�#@i��@i��@i�@i��@i�@i��@i�^@i�^@i��@ihs@i&�@h��@hĜ@h�9@h��@h��@hr�@hbN@hQ�@hA�@hA�@hA�@hA�@hQ�@hQ�@hQ�@hA�@h �@h  @g�;@g��@g�w@g��@g�P@g|�@g|�@g|�@gl�@g\)@gK�@g;d@g�@g
=@f�y@fȴ@f�R@fȴ@fȴ@f�R@f�R@f��@f��@fV@f5?@f$�@f$�@f$�@f$�@f$�@f{@e@e�@ep�@ep�@ep�@eO�@e�@d��@d�@d��@d�@d�D@dz�@dj@dZ@d9X@d(�@d1@c�m@c�
@cƨ@c�F@c��@c�F@c�F@c��@cdZ@cC�@c33@c33@co@b�@b�!@b�\@bn�@bM�@b�@a��@a��@a��@a�#@a��@aG�@a&�@`��@`�u@`�@`Q�@`1'@` �@`b@`  @_�@_��@_��@_��@_�P@_l�@_;d@_�@_
=@^�@^�+@^v�@^ff@^5?@^{@]�@]�T@]@]��@]�@]`B@]/@\�@\�@\�D@\(�@[�F@[�@[t�@[dZ@[C�@[o@[@Z��@Z��@Z�!@Z�\@Z�\@Z~�@Z^5@Z=q@ZJ@Y��@Y��@Y��@Y�7@Yx�@YX@Y7L@X��@X��@X��@X��@X��@X��@X�u@X�@Xr�@X1'@X �@Xb@Xb@X  @W�@W�;@W��@W�P@WK�@W+@W
=@V��@V��@V�y@V��@Vff@VV@V5?@U@U�-@U��@U�@U`B@U/@U�@UV@UV@T��@Tj@��@��@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@���@���@�@���@��@���@���@���@��@��@��@��@��@��@��@��@��y@��@��@��@��@��@��@��@��@��@��@���@���@��@���@��@��@��@��y@��y@��y@��H@��H@��H@��y@��y@��y@��H@��y@��@��@��y@��@��@��y@��y@��@��@��y@��y@��@��@���@���@���@���@���@��@��@���@��@��@��@��@��@��@��@��@��@��@��y@��@��@���@�@�@�@�@�@���@���@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�
=@�
=@�@�@�@�@�
=@�@�@�@���@�@���@���@���@�@�@�@�@�@�@���@���@�@���@���@���@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�
=@�@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�o@�
=@�o@�o@�o@�o@�o@�o@��@��@�"�@��@��@��@��@��@��@��@��@��@��@�"�@�"�@�"�@�"�@�"�@�"�@��@�"�@�"�@�"�@�"�@�"�@�+@�"�@�"�@�"�@�"�@�+@�33@�+@�"�@��@�@���@���@��@��@���@���@�@�@�@�@�
=@���@��H@���@�ȴ@���@���@���@��@��@��H@��@��H@��y@��y@��@���@��@��H@��H@��y@��@��H@���@���@���@���@�@�@�
=@�@�@���@�@�@�@��@���@�@���@���@�@���@���@���@��@��H@��@��@��@��@���@��H@��H@��H@��H@��y@��@��@��@��@���@���@��@���@�@���@���@��@��@��@��H@�ȴ@�ȴ@�ȴ@���@��R@��!@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��\@��\@��\@��\@���@���@���@��\G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   @���@�@�@�@�@���@��@��@��@��y@��y@��@���@�@���@�@�
=@�
=@�o@�"�@�o@��y@��y@�@��H@��y@���@���@���@��\@��\@�~�@�~�@��+@�v�@�n�@�n�@�^5@�V@�V@�V@�M�@�^5@�v�@���@��+@�n�@�n�@�n�@�n�@�E�@��T@���@�@���@���@���@���@���@���@���@���@���@��-@��h@�p�@�X@�`B@�O�@�V@���@���@���@���@���@��u@��u@��u@��u@��D@��@�r�@�bN@�Q�@�1'@��@�1@�b@�@��@�;@�;@��@�;@�w@�w@�w@��@�P@|�@|�@l�@l�@l�@\)@\)@K�@K�@K�@+@
=@
=@~��@~�y@~�@~ȴ@~��@~��@~��@~��@~��@~��@~��@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~�R@~ȴ@~ȴ@~ȴ@~ȴ@~�@~�@~�@~�@~�@~�@~�@~�@~�y@~�y@~�y@~�y@~�y@~��@~��@~��@~��@~��@
=@
=@
=@
=@
=@
=@
=@
=@
=@
=@�@
=@
=@�@�@�@+@�@
=@
=@
=@
=@~�y@~��@�@
=@~��@~@}@}��@}�h@}�h@}p�@}p�@}�@}�@}�@}�@}�@}�h@}�-@}��@~@~5?@~�+@~��@~�R@~ȴ@
=@
=@~��@~�@~�@~�@~ȴ@~ȴ@~ȴ@~�+@}�@}��@}��@}�h@}�-@}��@}�@}?}@|�@|(�@|�@{��@|9X@{�m@{�@{�@{dZ@y��@x�9@x1'@w�@w�w@w�P@w|�@w\)@w;d@w+@v�@v��@vff@v$�@u��@u�h@u�@up�@uO�@u?}@u/@u/@uV@t�@t�j@tj@tZ@t9X@t(�@t�@s��@s��@s�m@s�
@s�
@sƨ@s�F@s��@s��@sdZ@s33@r�\@r-@rJ@q�@q�@q�^@q�7@q&�@pQ�@p �@pb@o��@o�;@o��@o�@o��@o��@o�@o��@o\)@oK�@o\)@oK�@o;d@o�@o�@o+@o
=@nȴ@n�R@n��@n��@nff@nV@n@m�-@m�-@m��@m��@mp�@mO�@m`B@mO�@mV@l�/@l�/@l�j@l�j@l�j@l��@l�j@l��@lz�@lj@lI�@l9X@l�@k��@k�m@k�
@kƨ@k�F@k��@k�@kt�@kt�@kS�@k"�@ko@j�H@j��@j��@j�\@j�\@jn�@jM�@j=q@j-@j�@jJ@i��@i��@i�@i�#@i��@i��@i�7@ix�@ihs@iX@iG�@i&�@h�`@h��@h�@h�@hr�@hr�@hr�@hr�@hr�@hr�@h�@hr�@hbN@hbN@hr�@hbN@hA�@hA�@h1'@h �@g�@g��@g�w@g�@g�@g�P@g|�@g�P@g|�@gl�@g\)@g\)@gK�@g+@g�@g
=@f�R@f��@f��@f�+@f�+@fv�@f�+@fv�@fv�@fv�@fv�@fv�@fff@fV@fV@fff@fV@fE�@f$�@f@f@e�@e�@e�@e�T@e�@e�@e�@f@e�@e�T@e�-@e��@e�h@e�@e�@e?}@eV@eV@d�/@d��@d�/@d��@d�j@d�D@dz�@dz�@dz�@dj@dZ@dI�@d9X@d�@d�@d�@d�@d�@d�@d1@d1@c��@c�m@c�
@c�F@c��@ct�@cS�@cC�@c33@co@co@c@c@b�H@b�!@b��@b�\@b�\@b�\@b~�@b=q@b�@a��@a��@ahs@aX@`��@`Ĝ@`�9@`��@`�u@`�@`�@`�@`bN@`Q�@`A�@`A�@`1'@` �@`b@_�@_��@_�@_�@_�@_�@_�@_�P@_|�@_�P@_|�@_�P@_�P@_|�@_l�@_l�@_l�@_l�@_l�@_l�@_l�@_K�@_K�@_K�@_K�@_K�@_�@^�y@^�@^�@^ȴ@^��@^��@^��@^�+@^v�@^v�@^ff@^V@^V@^E�@^V@^V@^E�@^{@^@]�@]�T@]�T@]�T@]�@]�T@]�T@]�T@]�@]�@]�@]��@]@]�-@]�h@]p�@]p�@]`B@]?}@]?}@]?}@]?}@]?}@]/@]/@]?}@]?}@]?}@]?}@]�@]/@]/@]/@]/@]?}@]?}@]/@]?}@]?}@]/@]V@\��@]V@]V@]V@\��@\�@\�@\�@\�@\�/@\�/@\�/@\�/@\�@\�@\�@\�@\�@\�/@\��@\�j@\�j@\�j@\�D@\�D@\z�@\j@\j@\z�@\j@\I�@\9X@\9X@\9X@\(�@\(�@\(�@\(�@\(�@\(�@\�@\1@[��@[��@[�
@[�
@[�
@[ƨ@[�
@[�
@[�F@[��@[��@[��@[��@[��@[��@[��@[�@[t�@[dZ@[S�@[C�@[C�@[S�@[33@[33@["�@[o@[@[@Z�@Z�@Z�@Z�@Z�H@Z��@Z��@Z��@Z��@Z�H@Z�H@Z��@Z��@Z��@Z��@Z�H@Z�H@Z�H@Z�H@Z�@Z�@Z�@[@[33@[C�@[C�@[S�@[dZ@[dZ@[dZ@[dZ@[t�@[�@[��@[ƨ@[�
@[�m@[�m@[��@\1@\9X@\Z@\��@\�/@\�@]V@]V@]V@]/@]�h@]@^{@^ff@^v�@^�@_
=@_+@_l�@_�@_�;@` �@`A�@`A�@`bN@`�@`��@`��@a�@a�@a&�@aG�@aX@a�7@a�^@b�@b~�@b�\@b��@b��@c@c"�@cS�@cdZ@ct�@c�@c�@c��@c��@cƨ@c�m@dZ@dz�@d��@eV@e?}@e`B@e�@e@f@f{@fE�@fv�@f�+@f��@f�R@fȴ@f�@fȴ@fȴ@fȴ@fȴ@f�y@g
=@g�@g;d@gK�@g\)@g|�@g�P@h  @h�u@hĜ@h�9@h�9@h��@h�9@h��@h�`@h�`@i%@i7L@ihs@i��@ix�@ix�@i7L@h��@h�9@h��@h�`@h��@h��@i%@i�7@i�@i�@i�@i�#@i�#@i��@i��@i�@i��@i�@i��@i�^@i�^@i��@ihs@i&�@h��@hĜ@h�9@h��@h��@hr�@hbN@hQ�@hA�@hA�@hA�@hA�@hQ�@hQ�@hQ�@hA�@h �@h  @g�;@g��@g�w@g��@g�P@g|�@g|�@g|�@gl�@g\)@gK�@g;d@g�@g
=@f�y@fȴ@f�R@fȴ@fȴ@f�R@f�R@f��@f��@fV@f5?@f$�@f$�@f$�@f$�@f$�@f{@e@e�@ep�@ep�@ep�@eO�@e�@d��@d�@d��@d�@d�D@dz�@dj@dZ@d9X@d(�@d1@c�m@c�
@cƨ@c�F@c��@c�F@c�F@c��@cdZ@cC�@c33@c33@co@b�@b�!@b�\@bn�@bM�@b�@a��@a��@a��@a�#@a��@aG�@a&�@`��@`�u@`�@`Q�@`1'@` �@`b@`  @_�@_��@_��@_��@_�P@_l�@_;d@_�@_
=@^�@^�+@^v�@^ff@^5?@^{@]�@]�T@]@]��@]�@]`B@]/@\�@\�@\�D@\(�@[�F@[�@[t�@[dZ@[C�@[o@[@Z��@Z��@Z�!@Z�\@Z�\@Z~�@Z^5@Z=q@ZJ@Y��@Y��@Y��@Y�7@Yx�@YX@Y7L@X��@X��@X��@X��@X��@X��@X�u@X�@Xr�@X1'@X �@Xb@Xb@X  @W�@W�;@W��@W�P@WK�@W+@W
=@V��@V��@V�y@V��@Vff@VV@V5?@U@U�-@U��@U�@U`B@U/@U�@UV@UV@T��@Tj@��@��@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@���@���@�@���@��@���@���@���@��@��@��@��@��@��@��@��@��y@��@��@��@��@��@��@��@��@��@��@���@���@��@���@��@��@��@��y@��y@��y@��H@��H@��H@��y@��y@��y@��H@��y@��@��@��y@��@��@��y@��y@��@��@��y@��y@��@��@���@���@���@���@���@��@��@���@��@��@��@��@��@��@��@��@��@��@��y@��@��@���@�@�@�@�@�@���@���@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�
=@�
=@�@�@�@�@�
=@�@�@�@���@�@���@���@���@�@�@�@�@�@�@���@���@�@���@���@���@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@�
=@�@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�
=@�o@�
=@�o@�o@�o@�o@�o@�o@��@��@�"�@��@��@��@��@��@��@��@��@��@��@�"�@�"�@�"�@�"�@�"�@�"�@��@�"�@�"�@�"�@�"�@�"�@�+@�"�@�"�@�"�@�"�@�+@�33@�+@�"�@��@�@���@���@��@��@���@���@�@�@�@�@�
=@���@��H@���@�ȴ@���@���@���@��@��@��H@��@��H@��y@��y@��@���@��@��H@��H@��y@��@��H@���@���@���@���@�@�@�
=@�@�@���@�@�@�@��@���@�@���@���@�@���@���@���@��@��H@��@��@��@��@���@��H@��H@��H@��H@��y@��@��@��@��@���@���@��@���@�@���@���@��@��@��@��H@�ȴ@�ȴ@�ȴ@���@��R@��!@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@���@��\@��\@��\@��\@���@���@���@��\G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�1B�=B�=B�=B�DB�=B�=B�=B�JB�PB�PB�PB�PB�PB�PB�VB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�\B�VB�VB�VB�\B�VB�\B�\B�\B�\B�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B~�B~�B� B~�B~�B~�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�By�By�By�By�By�By�By�By�By�By�By�By�By�By�By�By�Bx�Bx�By�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{�B{�B|�B|�B|�B|�B|�B}�B~�B~�B� B�B�B�B�B�B�B�B�B�B�%B�%B�+B�1B�1B�7B�=B�=B�DB�JB�JB�JB�PB�PB�VB�VB�VB�\B�\B�\B�bB�hB�oB�uB�uB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�!B�!B�'B�'B�'B�'B�'B�!B�B�B�!B�'B�'B�-B�3B�?B�?B�?B�FB�FB�FB�FB�LB�LB�LB�LB�FB�FB�LB�LB�FB�FB�FB�FB�FB�FB�FB�FB�FB�LB�LB�LB�RB�LB�RB�LB�RB�RB�RB�RB�RB�RB�RB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�^B�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�RB�RB�RB�RB�RB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�XB�RB�RB�RB�RB�XB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�RB�XB�RB�XB�RB�RB�XB�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�B�B�%B�%B�B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�B�%B�B�B�%B�%B�%B�B�B�%B�B�%B�%B�%B�%B�%B�B�B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�B�B�%B�B�%B�B�+B�+B�%B�B�%B�%B�B�B�%B�%B�%B�B�%B�B�%B�%B�%B�B�%B�B�B�%B�B�%B�B�%B�B�B�%B�%B�B�%B�B�%B�+B�B�B�%B�B�B�%B�B�+B�B�%B�%B�%B�%B�%B�B�%B�B�%B�%B�B�B�%B�B�B�B�B�%B�B�B�+B�B�B�B�B�B�B�%B�B�B�+B�B�%B�B�B�B�B�B�B�B�B�B�B�B�%B�B�B�B�B�B�B�B�B�%B�%B�%B�B�B�%B�B�%B�B�B�%B�B�B�B�B�B�B�B�%B�%B�%B�B�B�B�B�B�B�B�B�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   B�B�$B�&B�*B�>B�:B�&B�*B�4B�%B�B�B�B�/B�B�B�$B�B�B�AB�_B�%B�B�VB�B��B�.B�B�-B�B�6B� B�B�8B�,B�#B�:B�,B� B�"B�-B�	B�B�B�PB�bB�CB�?B�@B��B��B�sB�[B�DB�SB�PB�PB�SB�RB�FB�OB�^B�wB��B��B�uB�HB�oB��B��B�SB�PB�cB��B�qB�[B�\B�ZB�bB�fB�qB�vB�uB��B��B�wB�XB��B�~B�YB�iB�rB�]B�B�hB�iB��B�tB�tB�jB�sB�lB�lB�zB�lB�zB�oB�oB��B��B�qB�zB�|B�}B�~B��B�tB�tB�tB�uB�tB�vB�iB�tB�vB�tB�tB�tB�tB�tB�vB�fB�vB�uB�qB�fB�tB�vB�wB�vB�tB�tB�tB�hB�wB�vB�vB�tB�hB�vB�tB�vB�tB�kB�vB�tB�tB�vB�wB�tB�wB�vB�vB�pB��B�|B�nB�xB�{B�nB��B��B��B��B��B��B��B��B��B��B�(B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B�"B�`B�&B��B��B��B��B��B�B�VB�EB��B��B��B�B� B��B��B��B��B�!B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�%B��B��B��B��B��B��B��B�:B��B��B��B�~B��B��B��B��B�|B��B��B��B�}B��B��B��B��B�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�mB��B��B��B�B��B�|B�{B�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�tB��B��B�~B��B��B��B�vB�iB��B��B�xB�vB�vB�vB�vB�kB�sB�wB�vB�~B�~B�pB�pB�oB�qB�~B��B��B�rB�]B�kB�]B�[B�]B�]B�\B�PB�jB�gB�]B�QB�kB�vB�[B�eB�cB�~B�pB�cB�cB�VB�nB�^B�EB�\B�\B�^B�QB�^B�kB�\B�`B��B�fB�LB�VB�MB�XB�=B�WB�KB�KB�HB�KB�ZB�YB�HB�>B�QB�TB�`B�]B�GB�OB�FB�GB�QB�8B�BB�EB�8B�OB�SB�hB�OB�KB�KB�BB�tB�aB�;B�\B�FB�-B�DB�EB�\B�FB�8B�8B�BB�=B�>B�?B�LB�3B�5B�3B�3B�5B�BB�3B�>B�@B�>B�HB�FB�DB�FB�<B�6B�@B�%B�0B�*B�@B�NB�0B�,B�!B�!B�0B�RB�;B�]B�AB�?B�*B�sB� B�B�B�B�B�B�B�$B�B�B�B�B�B�B�!B�"B�B� B�B�B�B�B�B�B�B�B�B�B�B� B� B� B� B�B�B�B~�B~�B�B~�B#B"BB~�BBBB}�B~B~B}�B~B~B}�B~B}�B}�B~B~B~B~B~B}�B}�B}�B~B}�B}�B}�B}�B}�B~B~B}�B}B}B|�B|�B}B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B}B|�B|�B|�B|�B|�B|�B|�B|�B|�B|�B}B|�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B|B{�B{�B{�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�By�By�By�By�By�By�By�By�By�By�By�By�By�By�By�By�Bx�Bx�By�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�Bx�By�By�By�By�By�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{�B{�B|�B|�B|�B|�B|�B}�B~�B~�B�B��B��B�B�B��B��B��B��B��B�B��B�B�B��B�B�B�B�-B�IB�0B�.B�5B�B�9B�UB�KB�AB�LB�7B�6B�B�$B�bB�eB�[B�HB�_B�\B�{B�|B�{B��B�|B��B�bB�rB�=B�|B�bB�rB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B�(B�B��B�B� B��B��B�B�@B�*B�\B�\B�RB�'B��B�-B�$B�B��B��B�6B�AB�LB�FB�SB�FB�-B�<B�YB�gB�UB�LB�cB�pB�~B�qB�jB�TB�RB�IB�lB�TB�SB�TB�KB�LB�NB�@B�MB�QB�[B�kB�lB�lB�aB�cB�oB�`B�dB�VB�XB�dB�fB�dB�fB�rB�iB�sB�oB�bB�IB�VB�aB�XB�dB�dB��B�pB�bB�VB�VB�VB�ZB�fB��B��B�fB�XB�XB�qB�|B�pB�gB�tB�rB�rB�eB�gB�eB�pB�eB�rB�pB�eB�eB�eB�bB�HB�WB�fB��B�rB�bB�VB�nB�sB��B�sB�pB�rB�~B�pB�XB�VB�rB��B��B�tB��B��B�dB�~B�rB�cB�cB�cB�dB�qB�VB�XB��B�tB�~B�qB�fB�B��B�fB�fB�{B�oB�qB�fB�qB�qB�rB�rB�{B��B��B�pB��B��B�B�bB�bB�rB�zB�dB�{B�YB�rB�rB�XB�bB�qB�rB��B��B�qB�XB�nB�dB�nB�rB��B�oB�zB�RB�UB�TB�_B�`B�bB��B�aB�]B�UB�`B�`B�aB�aB��B��B�oB�jB�aB�QB�`B��B��B�bB�rB��B�dB�`B�pB�lB�zB�_B�`B�YB��B��B�aB�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�B�B�%B�%B�B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�B�%B�%B�%B�%B�B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�%B�B�%B�B�B�%B�%B�%B�B�B�%B�B�%B�%B�%B�%B�%B�B�B�%B�%B�B�%B�%B�%B�%B�%B�%B�B�%B�B�B�%B�B�%B�B�+B�+B�%B�B�%B�%B�B�B�%B�%B�%B�B�%B�B�%B�%B�%B�B�%B�B�B�%B�B�%B�B�%B�B�B�%B�%B�B�%B�B�%B�+B�B�B�%B�B�B�%B�B�+B�B�%B�%B�%B�%B�%B�B�%B�B�%B�%B�B�B�%B�B�B�B�B�%B�B�B�+B�B�B�B�B�B�B�%B�B�B�+B�B�%B�B�B�B�B�B�B�B�B�B�B�B�%B�B�B�B�B�B�B�B�B�%B�%B�%B�B�B�%B�B�%B�B�B�%B�B�B�B�B�B�B�B�%B�%B�%B�B�B�B�B�B�B�B�B�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   <#؄<#�<#�<#�<#��<#�c<#�<#�<#׺<#�
<#�i<#ף<#׎<#�X<#ף<#�i<#�<#׺<#��<#�o<#�N<#�
<#��<#�E<#�X<#�N<#׺<#�
<#ף<#�
<#ا<#�<#�I<#��<#׎<#�<#�D<#׎<#�<#�<#ף<#�o<#��<#�<#�$<#�8<#�&<#�<#�<#�M<$$<#�+<#�i<#�{<#�<#�
<#�
<#�<#�<#�X<#�<#ף<#ۮ<#��<#��<#�8<#�<<#��<#��<#�"<#�<#�
<#�$<#�<#�D<#�<#�
<#�<#�{<#��<#�c<#�*<#��<#ޫ<#�r<#�D<#�X<#�+<#�o<#׺<#�<#�X<#�i<#ا<#�
<#�<#�D<#�{<#�{<#�<#�i<#�<#�<#�<#�<#�i<#�
<#�
<#ا<#�D<#�<#�i<#׎<#ף<#׺<#�C<#�<#�<#�<#�
<#�<#�<#�{<#�<#�<#�<#�<#�<#�<#�<#�<#׺<#�<#�
<#�<#׺<#�<#�<#�<#�<#�<#�<#�<#׎<#�<#�<#�<#�<#׎<#�<#�<#�<#�<#�X<#�<#�<#�<#�<#�<#�<#�<#�<#�<#�i<#�I<#�<#׎<#�<#�
<#׎<#�I<#׺<#�
<#�<#�<#�C<#�I<#�<#׎<#�<$
�<#�J<#�D<#׺<#�
<#ا<#�<#��<#�<#�<#�<#�<#��<#��<#�<#�+<#��<#�e<#׺<#�D<#��<#ߜ<#�
<#�i<#�o<#�<#�<#�i<#�<#�<#ޫ<$<<#�e<#�<#�{<#�c<#ף<#�D<#ߜ<#�a<#�	<#��<#ا<#��<#��<#�M<#�<#��<$��<$��<#��<#��<#��<#��<#ף<#��<#�D<#��<#�<#��<#�l<#��<#�&<#ޫ<#��<#��<#�D<#��<#ף<#�<#��<#�D<#ۮ<#�<#�{<#��<#׎<#�{<#�c<#�<#�{<#�{<#�<#�X<#�<#��<#�<#ۮ<#��<$<<#�<#�D<#؄<#�&<#ۮ<#�l<#�<$,<#��<#ף<#��<#׺<#�I<#�D<#�X<#�<#�i<#׺<#ޫ<#�{<#�X<#׎<#��<#��<#�<#�X<#�D<#�^<#�$<#��<#׺<#��<#��<#�&<#�<#�<#׺<#�<#ܯ<#�<#ף<#׺<#��<#��<#�<#��<#�<#�
<#�{<#׎<#�D<#�<#�{<#�<#�<#�o<#�<#��<#��<#ף<#׎<#׺<#�<#׎<#�<#�D<#��<#׺<#�8<#�<#�c<#ף<#�<#��<#��<#��<#ף<#ף<#ף<#ף<#�<#�i<#׺<#ף<#�o<#�o<#ף<#ף<#׎<#׺<#�o<#��<#�J<#؄<#�<#׺<#�<#�<#�<#�<#�
<#�{<#ף<#�i<#�<#�i<#׺<#�<#�<#׺<#׎<#��<#�<#׎<#׎<#�
<#��<#�<<#�i<#�{<#�{<#ף<#�<#ף<#�D<#�{<#��<#�&<#�o<#�<#�{<#�<#ף<#׎<#׎<#�<#�<#�<#�<#��<#׺<#�<#�&<#�0<#��<#�o<#��<#�<#�i<#�<#�<#׎<#�{<#�<#�<#�{<#�i<#׺<#��<#�<#ף<#ף<#�<#�E<#�l<#�<#�8<#׺<#�X<#׎<#ף<#�8<#׺<#�<#�<#��<#�{<#׎<#ף<#�D<#�<#�<#�<#�<#�<#��<#�<#׎<#׺<#׎<#ٛ<#�D<#��<#�D<#��<#�i<#�D<#�
<#�i<#�<#�D<#�+<#�i<#׎<#�<#�<#��<#��<#�o<#�&<#��<#�r<#��<#�<#�C<#׺<#׺<#׎<#׎<#�
<#�
<#��<#ף<#�{<#�<#�{<#ף<#ף<#�D<#�o<#�o<#�
<#�<#�<#�<#�D<#�{<#�i<#�{<#�i<#�<#ף<#׺<#�
<#�
<#�
<#�
<#�<#�<#؄<#�<#�<#�<#�<#�+<#��<#׎<#�
<#ף<#�<#�X<#�&<#׺<#׺<#�<#��<#׺<#�<#׎<#�<<#�<#�<#ܯ<#��<#׺<#��<#�<#�<#�I<#׺<#�<#�<#�X<#�<#�<#��<#ף<#�i<#�*<#�<#�<#ף<#�<#�<#�
<#�<#�<#�{<#�
<#ף<#�
<#�
<#�<#��<#�{<#�<#�
<#�<#ף<#�<#�i<#�i<#�<#ף<#��<#�i<#�X<#�<#�<#׎<#�i<#�<#�<#�<#׎<#�<#�<#�<#�i<#�<#�<#�<#�<#ף<#׺<#׎<#�
<#�<#�8<#�<#׎<#�i<#�<#�I<#׎<#�D<#ף<#�<#�<#�{<#�<#�<#�<#�
<#�<#׺<#�{<#�{<#�<#�D<#�<#�<#׎<#׎<#�<#�D<#׎<#�<#�<#�<#�<#�<#ף<#׺<#�{<#׺<#�I<#ף<#�<#�{<#�o<#�<#�{<#׎<#ף<#�<#�i<#�<#�<#�<#׎<#׎<#�{<#�<#�<#�D<#�<#�{<#�<#�<#�<#ף<#�<#�<#�<#׎<#�<#�
<#׺<#ۮ<#�i<#�<#�i<#�i<#�<#�
<#�<#�X<#�{<#�D<#ا<#��<#�i<#�<#�i<#׎<#ۮ<#�*<#ߜ<#ޫ<#��<#�c<#�<#�
<#��<#�<#�l<#�<#�e<#�<#�l<#�+<#�o<#��<#�J<#ܯ<#��<#ا<#�<#�<#�o<#�D<#�<#ٛ<#�<#�i<#�D<#��<#�8<#��<#�<#�4<#�$<#��<#�<#��<#�o<#�8<#�{<#�i<#�{<#�<#�i<#�<#ܯ<#�D<#�<#��<#��<#�J<#��<#��<#؄<#ޫ<#ޫ<#�C<#�8<#��<#��<#׺<#�<#׎<#׺<#�i<#�<#�<#�<#ا<#��<#�i<#�D<#�{<#��<#�<#ף<#�"<#��<#��<#�X<#�<#�$<#�<#�i<#�8<#�<#�D<#ۮ<#��<#�8<#��<#�<#ߜ<#ߜ<#�^<#�{<#ޫ<#�{<#�<#��<#��<#��<#�I<#�<#׎<#�
<#׎<#�
<#��<#��<#׎<#�D<#�I<#�&<#ٛ<#��<#ޫ<#ܯ<#��<#ף<#�{<#�<#�r<#ף<#׎<#ף<#�<#�
<#�<#�<#�<#�<#׺<#��<#�<#�<#׺<#��<#ٛ<#ף<#�{<#�<#�
<#�{<#ף<#�{<#ף<#�<#��<#�D<#ا<#�X<#׺<#�<#�I<#�
<#�{<#�{<#ޫ<#��<#�X<#�<#�<#�<#�<#ף<#�<#�J<#ף<#�
<#�
<#��<#��<#��<#׺<#�o<#�C<#�<#׎<#׺<#׎<#��<#׎<#�<#��<#׎<#׎<#׎<#�X<#��<#�<#ף<#�<#�<#�X<#�<#؄<#�D<#ߜ<#�D<#��<#�<#�r<#��<#�
<#�<#�<#�J<#�<#�o<#��<#��<#�{<#�r<#�<#�i<#�i<#�i<#�{<#��<#�<#�
<#�<#�o<#�r<#��<#ף<#ۮ<#�N<#ף<#ף<#��<#ا<#��<#ף<#��<#��<#�<#�<#��<#��<#ߜ<#��<#�!<#�5<#�8<#�X<#�X<#�<#ڑ<#�{<#��<#�<#�<#�<#�
<#�X<#��<#�<#�+<#ޫ<#��<#�
<#؄<#�{<#؄<#�*<#��<#ٛ<#��<#�&<#�<#�<#׎<#ף<#��<#ߜ<#׺<#�i<#�<#ף<#ף<#׺<#׺<#��<#ߜ<#ٛ<#��<#׺<#�<#ף<#�J<#��<#��<#�*<#�<#�<#ף<#��<#�<#��<#�0<#ף<#�<#�<#�4<#�I<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.008(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.008(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202208040000002022080400000020220804000000202208040000002022080400000020220804000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601283320181106012833QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601283320181106012833QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2021011500000020210115000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2022080100000020220801000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2021V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2021V03; CTD_for_DMQC_2021V02                     2022080400000020220804000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                