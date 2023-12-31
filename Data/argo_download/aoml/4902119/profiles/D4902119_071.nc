CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  f   N_CALIB       	N_HISTORY             
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
resolution        =���   axis      Z        0  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  O�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     0  T�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  g�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     0  l�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  �   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  ��   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �,   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  ��   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �(   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     0  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  ` $   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                   �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                   �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                   �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                   �   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                   �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                   �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                   �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  � �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                   x   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                   �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar        �   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�       �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �Argo profile    3.1 1.2 19500101000000  20181106012855  20230209090616  4902119 4902119 US ARGO PROJECT                                                 US ARGO PROJECT                                                 BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         BRECK OWENS, STEVEN JAYNE, P.E. ROBBINS                         PRES            TEMP            PSAL            PRES            TEMP            PSAL               G   GAA  AOAO6732                            6732                            2C  2C  DD  S2A                             S2A                             7365                            7365                            SBE602 ARM_v2.0_xmsg_ve         SBE602 ARM_v2.0_xmsg_ve         854 854 @�b��M@�b��M11  @�b����@�b����@N�*�)s�@N�*�)s��;Q֌i/�;Q֌i/11  GPS     GPS     Primary sampling: averaged [nominal 2 dbar binned data sampled at 0.5 Hz from a SBE41CP]                                                                                                                                                                        Near-surface sampling: discrete, pumped [data sampled at 1.0Hz from the same SBE41CP]                                                                                                                                                                                 AA  BA  BA  ?�z�@   @=p�@z�H@��R@�  @�  A   A  A%�A@��A`  A�  A��A�\)A�\)A��A�  A�Q�A��B   B�
B(�B(�B   B(  B0  B7�
B?�BG�BPQ�BW�B_\)Bg�
Bp(�Bw�
B�
B��
B��B��B��
B��B��B�  B�  B��
B�  B�{B��
B��B�{B��
B�B��B��
B��B�(�B�{B�{B�{B�{B��B�B��B�(�B�(�B�(�B�(�C 
=C
=C
=C  C  C	�
C�C�C
=C��C  C
=C�C�C�C
=C�C"  C${C%��C'�HC)�
C,  C.{C0
=C1��C3�HC6  C8�C:  C;�HC=�C@
=CB(�CD(�CF�CH  CI�HCL
=CN  CO��CR�CT�CU��CW��CZ
=C\  C]��C_�HCa�Cc�Cf  Ch�Cj{Cl
=Cm��Cp
=Cr
=Cs��Cu�HCw��Cy��C|
=C~{C�  C�C�C�C�\C���C���C�C�\C�C�  C�C�  C�C�
=C���C�  C���C���C���C���C�C�  C��C��C���C�C�C���C���C���C�
=C�\C�
=C�  C�  C���C�
=C�  C���C�  C�  C�  C�C�
=C�
=C�  C���C�  C�
=C�  C���C�C�
=C�C�C�C���C���C���C���C���C���C���C�  C�  C���C�  C�C�  C�  C���C���C�C�C�C�\C�C�  C�C�C�C���C���C�  C�C�  C���C�C�C���C�
=C�
=C���C���C���C���C�  C���C�  C�
=C�C�  C���C��C��C���C���C���C���C���C���C�  C�  C���C�  C���C�C�\C�C���C�C�  C��HC���C�  C�  C���D   D � D�D��D�D}qD�RDz�D�qDz�D�qD� D�D�D��Dz�D�D��D	  D	�D
  D
}qD
�qD}qD�qD��D  D� D�D� DD��D  D� D�D�DD� D  D}qD�qD}qD�D�D  Dz�D�qD� D�D��D  DxRD��Dz�D��D}qD�qD}qD��D� D�D}qD  D}qD   D ��D �qD!z�D"  D"� D"�qD#� D$D$�D%  D%}qD&  D&}qD&�qD'� D(�D(�D)  D)}qD)�qD*}qD*�qD+� D,�D,}qD-  D-��D.  D.��D/�D/}qD0  D0�D1�D1�D2�D2}qD2��D3}qD4�D4� D4�qD5}qD6�D6��D7  D7� D8  D8��D9D9�D9�qD:� D;D;��D<�D<}qD=�D=� D=��D>}qD?�D?�D@�D@z�D@��DA}qDB�DB��DC�DC�DD  DD��DE  DE}qDFDF�DF�qDG}qDH  DHu�DH�RDI}qDJ�DJ}qDK  DK�DK�qDLz�DL��DMxRDM�RDN� DO  DO� DPDP� DQ  DQ}qDQ�qDRz�DR�qDS��DT�DT�DU�DU�DU�qDV��DWDW��DX  DX��DY  DY��DZDZ��DZ��D[}qD\  D\��D\�qD]z�D]��D^}qD^�qD_}qD`�D`��D`��Da}qDa�qDb}qDc�Dc�Dd  Dd� De  De}qDf  Df}qDg�Dgz�Dg��DhxRDi  Di� Di�RDj� Dk�Dk�Dk��Dl� Dm�Dm� Dm�qDn� Dn�qDo��Dp�Dp� Dq�Dq�Dr  Drz�Dr�qDs��Ds�qDt� Du�Du��Dv�Dv��Dw  Dw�Dx�Dx��Dx�qDy��Dz�Dz��D{�D{��D|�D|� D}  D}��D~�D~� D~�qDz�D�qD�@ D�� D��HD�HD�@ D�~�D��qD�HD�@ D�� D���D���D�AHD�� D�D��D�@ D�~�D��)D��qD�=qD�|)D���D�HD�@ D���D���D���D�@ D�� D��HD�HD�>�D�}qD��qD�  D�AHD�� D���D��)D�=qD��D���D��qD�@ D�}qD��)D���D�@ D�~�D�D���D�<)D�~�D���D��D�B�D��HD��HD�D�B�D�}qD��)D��qD�=qD�~�D��HD�  D�AHD��HD���D�  D�@ D��HD�� D��qD�=qD�~�D���D��qD�>�D�~�D��qD���D�>�D�� D���D���D�AHD��HD��HD�HD�@ D�}qD���D��D�Ff?W
=?u?�  ?�\)?�{?���?�33?Ǯ?�
=?�(�?�@   @
=q@\)@z�@�R@&ff@.{@5@0��@J=q@Tz�@\(�@fff@fff@n{@n{@�  @��
@���@���@���@�z�@�Q�@��H@�G�@��@���@���@�
=@�
=@�(�@�  @��
@�ff@�=q@У�@�z�@�Q�@��H@޸R@��
@���@���@�\)@�z�@�Q�@�(�AG�A ��A�
AffA��A
�HA
�HA{A  AG�A33A�A�A=qA(�A{A   A!�A#�
A&ffA)��A*�HA,(�A.{A0  A1G�A333A5�A7�A:=qA;�A=p�AC�
AA�AC�
AE�AHQ�AJ=qAMp�AL��AN�RAP��AS33AUAUAW�AZ�HA\(�A]p�A`��Ab�\Ac�
AeAhQ�Ai��Ak�Amp�Ap  Aq�As�
Au�Aw�Ax��Az�HA}p�A~�RA�Q�A���A��\A�33A�z�A�p�A�{A��RA�Q�A�G�A��A�33A��A�p�A�ffA�\)A�Q�A��A��\A��A�(�A�z�A�{A�
=A�Q�A�G�A�=qA�33A�z�A��RA�A�\)A���A�G�A�=qA��A��A�A�ffA�Q�A���A���A��\A��
A���A�A�
=A�\)A���A�=qA�33A�(�A��A�ffA�
=A���A�G�A�=qA��
A�(�A��A�
=A�
=A�Q�A���A��A�33A�(�A��A�{A�\)A�Q�A���A�=qA�33A�(�A��A�{A�
=A�\)A��A�=qA��HAӅA���A�p�A�ffA׮Aأ�A�G�A�=qAۅA�z�A��A�ffA�\)A�  A���A��A�\A���A���A�\)A�
=A�  A�G�A��A�33A�(�A��A�{A�\)A�Q�A���A�(�A�33A��
A���A�{A�
=A��A���A�=qA��HA�{A���A�{A�
=B (�B ��B�B��B{B�\B�HB
=B(�B��B�BB{B�RB33B�B(�B��B	�B	��B
ffB
�RB33B�
BQ�B��BG�B�B=qB�HB�B(�BQ�BG�B��B�B
=B
=B\)Bz�B�B��B��B{B�\B33B�B��B��B�B��B{B�\B
=B�BQ�B��B�BB=qB
=B33B�
B (�B z�B!G�B!B!G�B"�HB#33B#�
B$Q�B$��B%G�B%�B&=qB&�RB'\)B(  B(Q�B(��B)p�B)�B*ffB+33B+\)B+�B,Q�B,��B-�B-B.ffB.�RB/33B/�
B0Q�B0��B1p�B1B2=qB2�HB3\)B3�B4Q�B4��B5�B5��B6{B6�RB7\)B8  B7�
B9�B9p�B9�B:ffB:�RB;�B<  B<  B=�B=p�B=�B>�\B>ffB?�B@(�B@z�BA�BA��BB{BB�RBC\)BC�BDQ�BD��BEG�BEBFffBF�HBG�BG�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            ?�z�@   @=p�@z�H@��R@�  @�  A   A  A%�A@��A`  A�  A��A�\)A�\)A��A�  A�Q�A��B   B�
B(�B(�B   B(  B0  B7�
B?�BG�BPQ�BW�B_\)Bg�
Bp(�Bw�
B�
B��
B��B��B��
B��B��B�  B�  B��
B�  B�{B��
B��B�{B��
B�B��B��
B��B�(�B�{B�{B�{B�{B��B�B��B�(�B�(�B�(�B�(�C 
=C
=C
=C  C  C	�
C�C�C
=C��C  C
=C�C�C�C
=C�C"  C${C%��C'�HC)�
C,  C.{C0
=C1��C3�HC6  C8�C:  C;�HC=�C@
=CB(�CD(�CF�CH  CI�HCL
=CN  CO��CR�CT�CU��CW��CZ
=C\  C]��C_�HCa�Cc�Cf  Ch�Cj{Cl
=Cm��Cp
=Cr
=Cs��Cu�HCw��Cy��C|
=C~{C�  C�C�C�C�\C���C���C�C�\C�C�  C�C�  C�C�
=C���C�  C���C���C���C���C�C�  C��C��C���C�C�C���C���C���C�
=C�\C�
=C�  C�  C���C�
=C�  C���C�  C�  C�  C�C�
=C�
=C�  C���C�  C�
=C�  C���C�C�
=C�C�C�C���C���C���C���C���C���C���C�  C�  C���C�  C�C�  C�  C���C���C�C�C�C�\C�C�  C�C�C�C���C���C�  C�C�  C���C�C�C���C�
=C�
=C���C���C���C���C�  C���C�  C�
=C�C�  C���C��C��C���C���C���C���C���C���C�  C�  C���C�  C���C�C�\C�C���C�C�  C��HC���C�  C�  C���D   D � D�D��D�D}qD�RDz�D�qDz�D�qD� D�D�D��Dz�D�D��D	  D	�D
  D
}qD
�qD}qD�qD��D  D� D�D� DD��D  D� D�D�DD� D  D}qD�qD}qD�D�D  Dz�D�qD� D�D��D  DxRD��Dz�D��D}qD�qD}qD��D� D�D}qD  D}qD   D ��D �qD!z�D"  D"� D"�qD#� D$D$�D%  D%}qD&  D&}qD&�qD'� D(�D(�D)  D)}qD)�qD*}qD*�qD+� D,�D,}qD-  D-��D.  D.��D/�D/}qD0  D0�D1�D1�D2�D2}qD2��D3}qD4�D4� D4�qD5}qD6�D6��D7  D7� D8  D8��D9D9�D9�qD:� D;D;��D<�D<}qD=�D=� D=��D>}qD?�D?�D@�D@z�D@��DA}qDB�DB��DC�DC�DD  DD��DE  DE}qDFDF�DF�qDG}qDH  DHu�DH�RDI}qDJ�DJ}qDK  DK�DK�qDLz�DL��DMxRDM�RDN� DO  DO� DPDP� DQ  DQ}qDQ�qDRz�DR�qDS��DT�DT�DU�DU�DU�qDV��DWDW��DX  DX��DY  DY��DZDZ��DZ��D[}qD\  D\��D\�qD]z�D]��D^}qD^�qD_}qD`�D`��D`��Da}qDa�qDb}qDc�Dc�Dd  Dd� De  De}qDf  Df}qDg�Dgz�Dg��DhxRDi  Di� Di�RDj� Dk�Dk�Dk��Dl� Dm�Dm� Dm�qDn� Dn�qDo��Dp�Dp� Dq�Dq�Dr  Drz�Dr�qDs��Ds�qDt� Du�Du��Dv�Dv��Dw  Dw�Dx�Dx��Dx�qDy��Dz�Dz��D{�D{��D|�D|� D}  D}��D~�D~� D~�qDz�D�qD�@ D�� D��HD�HD�@ D�~�D��qD�HD�@ D�� D���D���D�AHD�� D�D��D�@ D�~�D��)D��qD�=qD�|)D���D�HD�@ D���D���D���D�@ D�� D��HD�HD�>�D�}qD��qD�  D�AHD�� D���D��)D�=qD��D���D��qD�@ D�}qD��)D���D�@ D�~�D�D���D�<)D�~�D���D��D�B�D��HD��HD�D�B�D�}qD��)D��qD�=qD�~�D��HD�  D�AHD��HD���D�  D�@ D��HD�� D��qD�=qD�~�D���D��qD�>�D�~�D��qD���D�>�D�� D���D���D�AHD��HD��HD�HD�@ D�}qD���D��D�Ff?W
=?u?�  ?�\)?�{?���?�33?Ǯ?�
=?�(�?�@   @
=q@\)@z�@�R@&ff@.{@5@0��@J=q@Tz�@\(�@fff@fff@n{@n{@�  @��
@���@���@���@�z�@�Q�@��H@�G�@��@���@���@�
=@�
=@�(�@�  @��
@�ff@�=q@У�@�z�@�Q�@��H@޸R@��
@���@���@�\)@�z�@�Q�@�(�AG�A ��A�
AffA��A
�HA
�HA{A  AG�A33A�A�A=qA(�A{A   A!�A#�
A&ffA)��A*�HA,(�A.{A0  A1G�A333A5�A7�A:=qA;�A=p�AC�
AA�AC�
AE�AHQ�AJ=qAMp�AL��AN�RAP��AS33AUAUAW�AZ�HA\(�A]p�A`��Ab�\Ac�
AeAhQ�Ai��Ak�Amp�Ap  Aq�As�
Au�Aw�Ax��Az�HA}p�A~�RA�Q�A���A��\A�33A�z�A�p�A�{A��RA�Q�A�G�A��A�33A��A�p�A�ffA�\)A�Q�A��A��\A��A�(�A�z�A�{A�
=A�Q�A�G�A�=qA�33A�z�A��RA�A�\)A���A�G�A�=qA��A��A�A�ffA�Q�A���A���A��\A��
A���A�A�
=A�\)A���A�=qA�33A�(�A��A�ffA�
=A���A�G�A�=qA��
A�(�A��A�
=A�
=A�Q�A���A��A�33A�(�A��A�{A�\)A�Q�A���A�=qA�33A�(�A��A�{A�
=A�\)A��A�=qA��HAӅA���A�p�A�ffA׮Aأ�A�G�A�=qAۅA�z�A��A�ffA�\)A�  A���A��A�\A���A���A�\)A�
=A�  A�G�A��A�33A�(�A��A�{A�\)A�Q�A���A�(�A�33A��
A���A�{A�
=A��A���A�=qA��HA�{A���A�{A�
=B (�B ��B�B��B{B�\B�HB
=B(�B��B�BB{B�RB33B�B(�B��B	�B	��B
ffB
�RB33B�
BQ�B��BG�B�B=qB�HB�B(�BQ�BG�B��B�B
=B
=B\)Bz�B�B��B��B{B�\B33B�B��B��B�B��B{B�\B
=B�BQ�B��B�BB=qB
=B33B�
B (�B z�B!G�B!B!G�B"�HB#33B#�
B$Q�B$��B%G�B%�B&=qB&�RB'\)B(  B(Q�B(��B)p�B)�B*ffB+33B+\)B+�B,Q�B,��B-�B-B.ffB.�RB/33B/�
B0Q�B0��B1p�B1B2=qB2�HB3\)B3�B4Q�B4��B5�B5��B6{B6�RB7\)B8  B7�
B9�B9p�B9�B:ffB:�RB;�B<  B<  B=�B=p�B=�B>�\B>ffB?�B@(�B@z�BA�BA��BB{BB�RBC\)BC�BDQ�BD��BEG�BEBFffBF�HBG�BG�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@��@�@�X@�X@�hs@�X@�&�@��@�&�@�?}@�%@��`@�j@�  @��;@��
@��
@���@�w@��@��@睲@�t�@�S�@�+@�o@�E�@�%@�9@�r�@�1'@�1'@�1'@�(�@�1'@�1'@�1'@�(�@� �@�b@���@��;@���@���@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�F@�@�K�@�C�@�C�@�;d@�33@�33@��@�@�@�
=@�@���@��@���@��@���@���@���@��@��y@��y@��H@���@�ȴ@�ȴ@���@���@���@���@���@�R@���@�!@�+@�~�@���@��@�!@�;d@㕁@�@�@�|�@�j@��@�@�{@��@��@��@�{@���@��#@��#@��#@���@���@���@�@�^@�^@��@�7@�x�@�p�@�p�@�hs@�?}@�/@�7L@�7L@���@�j@�9@�@��@��@�z�@�Q�@�9X@� �@�  @��@�|�@�dZ@�K�@��H@�~�@�~�@�~�@�n�@�n�@�n�@�ff@�^5@�=q@��@��@�@�/@��`@��`@��`@��/@���@���@��@��@�@��@�bN@�9X@�1'@� �@� �@��@��@�b@�b@�1@���@��
@�|�@�|�@�\)@�K�@��@��H@�M�@��@݉7@�O�@�?}@�&�@�%@�Ĝ@܃@�9X@�9X@�9X@� �@�  @��@�ƨ@ۍP@�|�@�33@ڧ�@�E�@�{@���@�X@�%@��`@�Ĝ@ؼj@؛�@�r�@�Q�@�A�@� �@� �@� �@��@��@��@�1@���@��m@ץ�@�|�@ׅ@׍P@׍P@�S�@�
=@�"�@�+@�"�@�"�@�o@�
=@�@�@���@��H@��@�ȴ@ָR@֗�@�^5@�@��@���@��T@��T@��T@��@��@��@���@�@�J@��@���@Ցh@�x�@�hs@�hs@�?}@�7L@�&�@�%@��@���@Դ9@Լj@ԓu@�z�@�j@�A�@�9X@�1'@�(�@�(�@�  @���@�ƨ@ӶF@Ӯ@ӍP@�+@��y@��y@��H@Ұ!@ҧ�@җ�@�M�@�-@�-@�@�p�@���@�1@θR@�@�7L@�9X@� �@��@��@˥�@˕�@˅@�t�@�l�@�S�@�K�@�C�@�C�@�K�@�+@���@�
=@�
=@��@�"�@�"�@�"�@�o@ʇ+@ɩ�@�bN@���@��;@�\)@�V@���@ź^@ŉ7@�hs@�7L@ļj@�(�@å�@�o@��H@°!@�@�ff@���@��@��j@��u@��D@�j@� �@��F@�M�@���@�?}@�r�@���@�"�@��!@���@���@��h@�X@�G�@��@�&�@��@���@��@��u@��D@�Q�@��@��@�p�@�X@�&�@��/@�j@�1'@�dZ@��\@�ff@��@�&�@�A�@��@��H@�{@�/@���@��@�I�@�9X@�1'@�1'@�1'@�9X@�9X@�1'@��@��@�G�@�r�@�b@��
@�t�@�33@�
=@�v�@��@��^@���@���@��h@���@�/@��@��u@�V@�%@�1@�K�@�C�@��\@�5?@�@��@��D@��
@���@��w@���@�  @�  @�(�@�;d@�hs@��#@�(�@�M�@��@���@���@���@�`B@�V@��-@��@�`B@�O�@�/@��@���@�b@�C�@���@�{@���@�V@���@�1'@��
@��F@�l�@��@��H@��y@���@�v�@��+@���@�ȴ@�ȴ@��!@��R@�ȴ@���@��!@���@���@��!@��\@�M�@�=q@��@�J@�J@�J@���@��@��@��@��#@���@�@��^@��^@��-@��h@��@�G�@�&�@�%@��`@���@��@�z�@�Z@�I�@�b@��F@�|�@�t�@�l�@�dZ@�S�@�C�@�;d@�"�@��H@���@��+@�n�@�E�@�5?@���@��@�G�@���@���@���@��@��@��@��/@��`@���@��@��u@��u@��u@��u@��u@��D@��@�bN@�Q�@�Z@�Q�@�Z@�bN@�Z@�Z@�Z@�Q�@�A�@�1'@�9X@�9X@�1'@�1'@� �@� �@�(�@��@��@�(�@�b@���@���@���@�|�@�|�@�|�@�|�@�|�@��@��@��@��@��@��@��@��@��@��@��P@��@��@��@�t�@�|�@�|�@��@��@��@��P@��P@�dZ@�K�@�S�@�\)@�K�@�C�@�;d@�33@�C�@�;d@�33@�33@�33@�33@��@��@��@��@��@��@��@��@��@��@�h@�hs@�x�@�X@�X@�X@�O�@�X@�X@�X@�O�@�X@�O�@�X@�O�@�X@�X@�X@�O�@�X@�`B@�`B@�`B@�`B@�hs@�hs@�hs@�hs@�p�@�hs@�hs@�hs@�`B@�`B@�G�@�&�@��@�V@�V@�&�@�?}@��@�?}@�G�@�7L@��@�V@��@�V@�V@��@�%@�V@��@��@�&�@��@��@�/@�/@�/@�7L@�7L@�7L@�?}@�O�@�7L@�7L@�/@�?}@�G�@�O�@�G�@�?}@�&�@�/@�V@���@���@�%@�V@�%@���@���@���@�Ĝ@���@��@�V@��@���@���@���@��/@���@��`@��`@��@��/@��`@���@���@��@�%@�%@���@��`@�j@�j@�z�@�z�@�@�D@�z�@�Z@�Q�@�Z@�I�@�A�@�(�@��@�b@�b@��@�1@�  @�  @�  @���@�  @��@��m@���@���@��@�  @�  @���@��m@��;@��m@��@��m@��@��m@��m@��;@��;@��
@��
@��;@��;@��
@��
@��
@��
@��
@��
@��
@���@��
@���@���@��
@��
@���@���@��
@��
@��
@��
@���@���@��
@��
@��
@��
@��
@��
@��
@��
@���@���@��
@��
@��
@��
@���@���@���@���@�ƨ@�ƨ@���@�ƨ@�ƨ@�ƨ@���@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@���@�ƨ@�w@�w@�w@�ƨ@�ƨ@�w@�w@�w@�F@�@�@��@��@�@��@��@��@��@��@睲@睲@��@��@��@�@��@�@��@��@�@��@��@睲@睲@睲@睲@睲@��@��@睲@睲@��@��@��@��@��@睲@畁@睲@睲@睲@畁@畁@睲@畁@畁@畁@畁@畁@畁@�P@畁@�P@�@�@�|�@�dZ@�dZ@�dZ@�dZ@�\)@�\)@�dZ@�\)@�\)@�\)@�\)@�\)@�\)@�\)@�\)@�S�@�K�@�S�@�S�@�S�@�K�@�K�@�K�@�C�@�C�@�C�@�C�@�33@�;d@�C�@�33@�+@�"�@��@��@�"�@��@��@��@��@��@��@�o@��@�o@�o@�o@�o@�o@��@��@��@�o@�o@�
=@�@�@�@�@��y@��H@��@���@�\@�$�@��T@�^@��@�7@�x�@�p�@�p�@�hs@�X@�/@�&�@��@�%@�%@���@��`@��`@��/@��`@��/@���@��/@���@���@�Ĝ@�j@�9@�j@�9@�9@�@�@�@�@�@�@�@�@�@��@�u@�@�@�z�@�z�@�r�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            @��@�@�X@�X@�hs@�X@�&�@��@�&�@�?}@�%@��`@�j@�  @��;@��
@��
@���@�w@��@��@睲@�t�@�S�@�+@�o@�E�@�%@�9@�r�@�1'@�1'@�1'@�(�@�1'@�1'@�1'@�(�@� �@�b@���@��;@���@���@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�F@�@�K�@�C�@�C�@�;d@�33@�33@��@�@�@�
=@�@���@��@���@��@���@���@���@��@��y@��y@��H@���@�ȴ@�ȴ@���@���@���@���@���@�R@���@�!@�+@�~�@���@��@�!@�;d@㕁@�@�@�|�@�j@��@�@�{@��@��@��@�{@���@��#@��#@��#@���@���@���@�@�^@�^@��@�7@�x�@�p�@�p�@�hs@�?}@�/@�7L@�7L@���@�j@�9@�@��@��@�z�@�Q�@�9X@� �@�  @��@�|�@�dZ@�K�@��H@�~�@�~�@�~�@�n�@�n�@�n�@�ff@�^5@�=q@��@��@�@�/@��`@��`@��`@��/@���@���@��@��@�@��@�bN@�9X@�1'@� �@� �@��@��@�b@�b@�1@���@��
@�|�@�|�@�\)@�K�@��@��H@�M�@��@݉7@�O�@�?}@�&�@�%@�Ĝ@܃@�9X@�9X@�9X@� �@�  @��@�ƨ@ۍP@�|�@�33@ڧ�@�E�@�{@���@�X@�%@��`@�Ĝ@ؼj@؛�@�r�@�Q�@�A�@� �@� �@� �@��@��@��@�1@���@��m@ץ�@�|�@ׅ@׍P@׍P@�S�@�
=@�"�@�+@�"�@�"�@�o@�
=@�@�@���@��H@��@�ȴ@ָR@֗�@�^5@�@��@���@��T@��T@��T@��@��@��@���@�@�J@��@���@Ցh@�x�@�hs@�hs@�?}@�7L@�&�@�%@��@���@Դ9@Լj@ԓu@�z�@�j@�A�@�9X@�1'@�(�@�(�@�  @���@�ƨ@ӶF@Ӯ@ӍP@�+@��y@��y@��H@Ұ!@ҧ�@җ�@�M�@�-@�-@�@�p�@���@�1@θR@�@�7L@�9X@� �@��@��@˥�@˕�@˅@�t�@�l�@�S�@�K�@�C�@�C�@�K�@�+@���@�
=@�
=@��@�"�@�"�@�"�@�o@ʇ+@ɩ�@�bN@���@��;@�\)@�V@���@ź^@ŉ7@�hs@�7L@ļj@�(�@å�@�o@��H@°!@�@�ff@���@��@��j@��u@��D@�j@� �@��F@�M�@���@�?}@�r�@���@�"�@��!@���@���@��h@�X@�G�@��@�&�@��@���@��@��u@��D@�Q�@��@��@�p�@�X@�&�@��/@�j@�1'@�dZ@��\@�ff@��@�&�@�A�@��@��H@�{@�/@���@��@�I�@�9X@�1'@�1'@�1'@�9X@�9X@�1'@��@��@�G�@�r�@�b@��
@�t�@�33@�
=@�v�@��@��^@���@���@��h@���@�/@��@��u@�V@�%@�1@�K�@�C�@��\@�5?@�@��@��D@��
@���@��w@���@�  @�  @�(�@�;d@�hsG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@�/@��@���@�b@�C�@���@�{@���@�V@���@�1'@��
@��F@�l�@��@��H@��y@���@�v�@��+@���@�ȴ@�ȴ@��!@��R@�ȴ@���@��!@���@���@��!@��\@�M�@�=q@��@�J@�J@�J@���@��@��@��@��#@���@�@��^@��^@��-@��h@��@�G�@�&�@�%@��`@���@��@�z�@�Z@�I�@�b@��F@�|�@�t�@�l�@�dZ@�S�@�C�@�;d@�"�@��H@���@��+@�n�@�E�@�5?@���@��@�G�@���@���@���@��@��@��@��/@��`@���@��@��u@��u@��u@��u@��u@��D@��@�bN@�Q�@�Z@�Q�@�Z@�bN@�Z@�Z@�Z@�Q�@�A�@�1'@�9X@�9X@�1'@�1'@� �@� �@�(�@��@��@�(�@�b@���@���@���@�|�@�|�@�|�@�|�@�|�@��@��@��@��@��@��@��@��@��@��@��P@��@��@��@�t�@�|�@�|�@��@��@��@��P@��P@�dZ@�K�@�S�@�\)@�K�@�C�@�;d@�33@�C�@�;d@�33@�33@�33@�33@��@��@��@��@��@��@��@��@��@��@�h@�hs@�x�@�X@�X@�X@�O�@�X@�X@�X@�O�@�X@�O�@�X@�O�@�X@�X@�X@�O�@�X@�`B@�`B@�`B@�`B@�hs@�hs@�hs@�hs@�p�@�hs@�hs@�hs@�`B@�`B@�G�@�&�@��@�V@�V@�&�@�?}@��@�?}@�G�@�7L@��@�V@��@�V@�V@��@�%@�V@��@��@�&�@��@��@�/@�/@�/@�7L@�7L@�7L@�?}@�O�@�7L@�7L@�/@�?}@�G�@�O�@�G�@�?}@�&�@�/@�V@���@���@�%@�V@�%@���@���@���@�Ĝ@���@��@�V@��@���@���@���@��/@���@��`@��`@��@��/@��`@���@���@��@�%@�%@���@��`@�j@�j@�z�@�z�@�@�D@�z�@�Z@�Q�@�Z@�I�@�A�@�(�@��@�b@�b@��@�1@�  @�  @�  @���@�  @��@��m@���@���@��@�  @�  @���@��m@��;@��m@��@��m@��@��m@��m@��;@��;@��
@��
@��;@��;@��
@��
@��
@��
@��
@��
@��
@���@��
@���@���@��
@��
@���@���@��
@��
@��
@��
@���@���@��
@��
@��
@��
@��
@��
@��
@��
@���@���@��
@��
@��
@��
@���@���@���@���@�ƨ@�ƨ@���@�ƨ@�ƨ@�ƨ@���@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@���@�ƨ@�w@�w@�w@�ƨ@�ƨ@�w@�w@�w@�F@�@�@��@��@�@��@��@��@��@��@睲@睲@��@��@��@�@��@�@��@��@�@��@��@睲@睲@睲@睲@睲@��@��@睲@睲@��@��@��@��@��@睲@畁@睲@睲@睲@畁@畁@睲@畁@畁@畁@畁@畁@畁@�P@畁@�P@�@�@�|�@�dZ@�dZ@�dZ@�dZ@�\)@�\)@�dZ@�\)@�\)@�\)@�\)@�\)@�\)@�\)@�\)@�S�@�K�@�S�@�S�@�S�@�K�@�K�@�K�@�C�@�C�@�C�@�C�@�33@�;d@�C�@�33@�+@�"�@��@��@�"�@��@��@��@��@��@��@�o@��@�o@�o@�o@�o@�o@��@��@��@�o@�o@�
=@�@�@�@�@��y@��H@��@���@�\@�$�@��T@�^@��@�7@�x�@�p�@�p�@�hs@�X@�/@�&�@��@�%@�%@���@��`@��`@��/@��`@��/@���@��/@���@���@�Ĝ@�j@�9@�j@�9@�9@�@�@�@�@�@�@�@�@�@��@�u@�@�@�z�@�z�@�r�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B"�B$�B$�B$�B%�B+B33B5?B6FB6FB6FB6FB6FB7LB7LB7LB7LB7LB7LB7LB7LB7LB6FB6FB6FB6FB6FB6FB6FB6FB6FB6FB5?B5?B5?B5?B5?B5?B5?B49B49B49B33B33B2-B2-B2-B1'B1'B0!B0!B0!B0!B0!B0!B0!B0!B/B/B.B.B.B.B/B/B/B0!B/B/B/B/B/B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B0!B/B/B.B,B+B(�B(�B'�B'�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B&�B%�B$�B#�B#�B"�B"�B �B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B �B �B!�B!�B!�B!�B!�B"�B"�B"�B"�B"�B"�B"�B#�B$�B%�B&�B)�B+B,B-B-B-B/B0!B0!B1'B2-B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B49B5?B49B49B5?B49B49B5?B5?B49B49B49B33B2-B0!B/B.B-B-B-B,B,B,B-B-B-B.B/B/B0!B0!B2-B5?B5?B5?B5?B49B49B49B49B6FB5?B7LB6FB6FB5?B49B49B49B49B33B33B2-B1'B0!B0!B0!B0!B0!B/B/B/B/B/B/B.B.B-B-B,B,B+B)�B)�B(�B'�B'�B'�B&�B'�B&�B'�B&�B&�B&�B%�B%�B$�B#�B"�B!�B �B �B�B�B�B�B�B�B�B�B�B�B�B�B{BuBuBuBuBuBuBuBoBoBhBbB\BVBVBPBPBJBJBDBDB
=B
=B
=B
=B
=B	7B	7B1B1B1B+B1B%BBBBBBB  B  B  B  B  B  B��B��B��B��B�B�B�B�B��B��B��B�B�B�B�B�yB�yB�sB�mB�fB�`B�ZB�ZB�TB�NB�NB�NB�NB�HB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�BB�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�;B�5B�5B�5B�5B�5B�5B�/B�/B�/B�/B�/B�/B�/B�/B�/B�/B�/B�/B�/B�/B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�)B�#B�#B�)B�)B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            B�B�B�BpB�B B�B\BdB�B�BBB,B�B�B�B�B�B�B�B�B�B�B�B�B�BwB$BB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BOB�B�B�B"@B$�B%B$�B$]B)B2�B5 B68B6FB6IB6UB6sB7wB7KB7LB7VB7NB7XB7XB7VB7LB6bB6{B6`B6RB6HB6XB6�B6`B69B6IB5�B5�B5LB5MB5MB5CB5{B4~B4_B4^B3jB3�B2vB2SB2RB1�B1�B0&B0 B07B0!B0 B0.B0.B0UB/�B/�B.PB.�B.�B.B/B/*B/+B0+B/^B/HB.�B/VB/OB0VB0.B06B0!B0)B0 B0,B0#B0,B08B0RB0�B0&B0UB0<B/jB/uB.�B,�B+�B)OB)B(B($B'UB'NB'UB&�B&�B'B'B'B''B'BB'B&UB%�B$lB$&B#>B#�B!DB�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�BB�B�B�B�BBB�B �B �B!�B!�B!�B!�B!�B"�B"�B"�B"�B"�B#B#/B$[B$�B& B&�B)�B+B+�B-B-B-B/B0B0IB1eB2�B4YB4PB4<B4tB4DB4SB4jB4^B4hB4\B4/B4sB4]B4SB4xB4HB4EB4EB4<B4vB4�B4HB4UB4JB4pB4�B4�B5@B4LB4�B5JB4ZB4�B5nB5BB4�B5B5	B4�B43B1FB0ZB/�B-BB-!B-PB,wB,"B,$B-(B-B-6B.$B/&B/"B0B0ZB2tB5$B5>B5(B50B4:B49B4YB5B7�B79B7�B6�B7B6�B5
B4ZB4�B4nB3�B3�B3B1�B0�B0pB0mB0@B0�B0ZB0B/lB/XB/,B/NB.�B.�B/ B-�B,�B-<B+�B+B*�B*
B(xB(B(BB'
B(&B&�B(
B'�B'B&�B%�B&:B%�B%,B%DB"B!B!8B tB $B�B�B�BjB�B�B�B�B�B�BB�B�B�B�BwBuBhBsB�BB�B�B�B�B�B�B�B�B$B�B
�B
YB
IB
SB
6B	�B	�B�B�BKB�B	OBIB:B�B�B|B�BB B B��B��B��B��B rB��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B��B�B�KB�B�qB�!B��B�B��B�=B��B�B��B��B�^B�<B�B�|B�&B��B�4B�?B�eB�4B�.B��B�B�JB� B�[B�wB�B�\B�jB�UB�<B�>B�SB�JB�>B�IB�TB�JB�TB�JB�;B�GB�hB�QBޏB�fB�hB�jB�ZB�dB�zB�eB�MB݆BݺB݊B�>B�>B�>B�IB�IB�<B�XBܔBܑB�TB�NB�hB�HB܊B��B܂BܚB�(B�4B�6B�*B�)B�AB�B�QB�VB�KB�'B�#B�(B�&B�.B�1B�TB�<B�B�.B�B�B�-B�#B�&B�4B�=B�8B�B�#B�1B�%B�<B�#B�B�;B�#B�B�NBۅB�_B�5B�SB�$B�$B�#B�#B�B�"B�#B�$B�$B�#B�%B�"B�"B�"B�B�0B�$B�$B�>B�B�#B�B�#B�#B�B�&B�`B�IB�B�B�<B�1B�/B�0B�	B�1B�1B�#B�#B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444444441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                            <#�a<#�W<#�<#ا<#�&<$�<#��<#ܯ<#��<#�<#��<$A�<$*<#ߜ<#ף<#�<#׎<#�D<#�8<#�
<#��<#�<#��<#�<#�E<%:<&��<$.<#��<#�"<#�<#�<#�{<#�I<#�<#�<#��<#׎<#�<#ۮ<#�8<#�<#�<#�i<#�
<#�<#�<#�
<#�<#�<#�<#�D<#�M<#�<#׎<#�<#��<#׺<#�<#ۮ<#��<#�<#�X<#�{<#�{<#�{<#�X<#�I<#׎<#�<#�
<#�{<#׺<#�<#ף<#��<#׎<#�<#ף<#�<#�<#�<#�
<#�i<#�<#��<#��<#�<<#�<#��<#�C<$W<$'<#�+<#��<#�<%��<&��<$&<#��<#ף<#�
<#�<#׺<#�8<#ܯ<#�<#�
<#�X<#�<#�{<#�{<#�X<#�
<#�o<#ߜ<#�<#�{<#�<#�<#�e<#�<#׎<#�<#��<#��<#׎<#ף<#ף<#�<#�<#�<#�r<#�8<#�E<$<<#�M<#�r<#�8<$!><$
<#�<#�<#؄<#�
<#�<#׎<#׎<#�J<$v<#��<#�<$�<#�N<#�0<#�
<#׺<#��<#�X<#�<#�8<#��<#�<#�J<#ߜ<#׎<#�c<#�
<#�<<#�<#�i<#�<#�i<#ا<#�^<$�<#�<#�J<#�D<#�<#��<$q@<$�<$k<#�5<#�<#�8<#�J<#��<#�)<#��<#�
<#�<#�+<#��<#�o<#��<#�5<#�o<#��<$\"<$�<#�<#�H<$9�<$<<#��<#�<#�<#�^<#�e<#��<#��<#ޫ<#�<#�<#�<<#�<#�<#ا<#��<#�*<#�5<#�N<#�i<#ף<#�<#�<#�)<#�<#׎<#�{<#�<#�o<#׺<#�{<#�<#׎<#�8<#׺<#��<#�o<#�E<#��<$.<#��<#�e<#�D<#�I<#�<#�{<#�<#�{<#�0<#��<#�{<#��<#��<#�5<#�*<#ا<#�<#�<#�i<#�<#�^<#�8<#��<#��<#�X<#�N<#��<#�<#�&<#׺<#�{<#�{<#�<#�e<#�<#׺<#�o<#��<#�E<$F<#��<#�<#�$<#�<#�i<#�]<#�a<#��<#�<#�<$n�<$Z�<%}�<'<$�7<%s<%it<#�J<#�$<#�U<#��<#�<#�o<#�<#׺<#��<#��<#�i<#�0<#�{<#��<#�l<#�D<#�<#ا<#׺<#�<#�
<#�*<$k�<%U�<&�A<$Z<#�N<$aD<%��<$\"<#�]<#�<#ߜ<#�"<$Gd<$v�<$T�<$k�<#�<#�<#��<#�<%s<$��<#�<#�e<#��<#��<$<<$Gd<''�<$Y�<$��<$��<$��<$�3<$E<$��<$r<#��<#�<#�]<#��<#�<<#�<$8�<#�+<#��<#�$<#�"<$(<%8j<(v<#�<#�<#�a<$<<#��<%0<%
�<#�<$L<$�L<%B�<$�q<$��<%MY<%K:<$�<$v<#�<#�D<#�{<#�<#�
<#׎<#�<#ا<$��<'�|<%2?<%<$*<#�N<$'<#�	<#��<$o�<$G<$#(<#�o<#�{<#؄<#�0<$}<#�	<$�<$0.<#�<%��<$��<#��<$ʾ<$�<$+<#�N<%��<$��<#׺<#؄<#��<#�l<#�0<#��<%�J<)Ɩ<$�3G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�<#�N<$f<$m,<%s<$��<$O�<$ <$�L<$p<$� <$f<#��<$�<$9�<#�o<#�&<#��<#�N<#�o<#�N<#ף<#�<#��<#ף<#�C<#�<#��<#�<<#ڑ<#��<#ߜ<#��<#�]<#��<#�<#�<#�<#��<#׺<#�<#ף<#��<#׺<#��<#׺<#�
<#�{<#��<#�o<#��<#�^<#��<#ߜ<#ܯ<#ߜ<#�4<#��<#��<#�"<$�<#�N<#׺<#׺<#׺<#�<#�<#׎<#��<#��<#�	<#ܯ<#�8<#�&<#��<#�<$F9<#�5<#��<#�<#�i<#׎<#�<#�
<#��<#�X<#��<#�8<#��<#�<#�&<#�<#�<#�i<#ף<#�^<#��<#�X<#�i<#�i<#�{<#�X<#�
<#�<#��<#�<#�c<#ף<#�
<#ף<#�<#��<#�
<#�X<#��<#�
<#�c<#ܯ<#�W<#�<#�<#�<#�<#�<#�
<#�
<#ף<#�<#�
<#�<#�<#�
<#�<#�<#�<#�<#�X<#׎<#�<#�<#�D<#�I<#�
<#ף<#�
<#�
<#ף<#�<#�e<#�r<#�X<#�i<#��<#ף<#�{<#׎<#�<#ף<#ף<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = CTM_ADJ_PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                              PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJ = PSAL, multiplicative adjustment term r = 1, no additional adjustment necessary.                                                                                                                                                                      None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment;OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                   None                                                                                                                                                                                                                                                            None                                                                                                                                                                                                                                                            OW: r =1(+/-0), vertically averaged dS =0.003(+/-0.001),                                                                                                                                                                                                        SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT.; No significant drift detected in conductivity                                                                                                                                   SOLO-W floats auto-correct mild pressure drift by zeroing the pressure sensor while on the surface.  Additional correction was unnecessary in DMQC;      PRES_ADJ_ERR: SBE sensor accuracy + resolution error                                                   No significant temperature drift detected;  TEMP_ADJ_ERR: SBE sensor accuracy + resolution error                                                                                                                                                                No thermal mass adjustment on non-primary profiles.; No significant drift detected in conductivity                                                                                                                                                              202302090000002023020900000020230209000000202302090000002023020900000020230209000000AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285520181106012855QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�5F03E           703E            AO  AO  ARGQARGQQCPLQCPL                                                                                                                                        2018110601285520181106012855QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�0               0               WHOIWHOIARSQARSQWHQCWHQCV0.5V0.5                                                                                                                                2020010700000020200107000000QC  QC                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOI    ARSQ    WHQC    V0.5                                                                                                                                    20200709000000              CF                                      G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARSQARSQCTM CTM V1.0V1.0                                                                                                                                2023020700000020230207000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                WHOIWHOIARCAARCAOWC OWC V2.0V2.0ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     ARGO_for_DMQC_2022V03; CTD_for_DMQC_2021V02                     2023020900000020230209000000IP  IP                                  G�O�G�O�G�O�G�O�G�O�G�O�                                