CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-07-26T08:51:11Z creation; 2020-10-16T13:43:57Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_030d      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    7<   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    7L   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    7P   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    7T   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7d   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7t   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  �  7�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  �  8   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  `  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9    DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                  @  9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9D   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9L   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                  @  9P   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                  @  9�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                  @  9�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    :   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           :   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    :(   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            :,   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           :<   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           :L   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    :\   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    :`   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    :p   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :t   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    :x   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    :|   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <|   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %8.2f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  L   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %8.2f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  P   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %8.2f      FORTRAN_format        F7.1   
resolution        =���     �  c�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  s   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  ��   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  �   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �L   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  �4   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     �  д   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �    HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  �  �    HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                     ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                     ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  �L   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T  �Argo profile    3.1 1.2 19500101000000  20190726085111  20201016134357  3901988 3901988 Argo-Norway                                                     Argo-Norway                                                     Kjell Arne Mork                                                 Kjell Arne Mork                                                 PRES            TEMP            PSAL            PRES            TEMP            PSAL               ,   ,AA  IFIF                                                                2C  2B  DA  APEX                            APEX                            8415                            8415                            2.10.1                          2.10.1                          846 846 @�Z7x�Hp@�Z7x�Hp11  @�Z7���@�Z7���@Q���Y��@Q���Y���.2���m�.2���m11  GPS     GPS     AA  AA  AA  Primary sampling: mixed [averaged from 4.42 dbar to 997.72 dbar {1 sample/second, 2 dbar average from 2100 dbar to 4 dbar}; discrete otherwise]                                                                                                                 Secondary sampling: discrete [1 sample, 10 dbar interval from 2100 dbar to 4 dbar]                                                                                                                                                                                    @�p�@�
=@��A�A?�A`Q�A�(�A�(�A�Q�A�=qA�  A�  A�(�A�  A��B  B��B�B�B(  B0
=B7�HB?��BG�HBO��BW��B`  Bh
=Bp  Bw��B�HB�B�{B�  B��B���B���B�  B�B�
=B�B��B���B�B�\B�B���B�B�B���B���B�
=B�#�B��B�  B��B��B�
=B�\B�  B�B�{C \C�CC�RC��C
  C�qC�RC�qC�C�C
=C�C��C  C�C �C!�qC#�qC%��C(  C*�C+�qC-�RC/�RC2  C4  C5�qC7��C9�qC<  C>C@  CA�qCD�CF  CG�qCI�qCK��CM��CO�qCQ��CS��CU�RCW�RCY��C\  C^  C_�qCa�qCd�Cf�Ch�Ci��Ck�RCn  Cp�Cq�qCtCv  Cw��Cy��C|  C}��C��C��qC��qC�  C��qC��qC��qC���C��C�HC��)C��qC���C�  C�  C�HC���C���C��C��C�HC�  C�  C�HC�C�  C��qC�  C�  C�C�HC�  C�  C���C�  C�HC���C���C���C�C��C���C��qC��)C��)C���C��qC�HC��C��C�HC��qC�HC�C��C��qC��qC��C�  C�  C��C�fC�  C���C�  C��)C��qC�  C��C�C��C�C��qC���C���C���C��RC���C���C��)C�  C�HC���C�  C�  C���C�  C��C�HC��)C���C�HC���C�  C��qC���C�HC��qC�HC�C��C�HC��C�  C�  C��C��C�HC�  C�  C�HC�HC�  C��qC�  C��C��C�HC�HC���C��)C�  C���C���C���C��C��C��D HD �HD�D�3D3D�HD  D~�D��D}qD  D�HD�D�HD��D~D�\D\D�D	~D	��D
~�D  D�HD  D\D��D\D  D\D �D��D �D� DHD��D�\D��D�D� D�\D~D�qD� D �D� D�\D~�D�qD|�D�qD\DHD�HD �D��D3D��D�\D�HD�D��DHD��D �D � D �qD!\D"  D"� D#�D#� D#�\D$� D$�D%~D%�\D&� D&�D'\D( �D(�HD) �D)��D*HD*�HD+HD+� D+�qD,~D,��D-\D.  D.��D/�D/��D0�D0��D1�D1� D2 �D2��D2�D3\D4 �D4\D4�qD5~�D5�D6~D7HD7�HD8  D8��D9 �D9~�D9��D:\D:�\D;\D< �D<~D=  D=�HD=�\D>� D?�D?� D@HD@�3DAHDA\DB �DB��DCHDC� DC�\DD\DD�\DE\DF�DF��DF�\DG��DG�\DH�HDI�DI�HDJ �DJ�HDK  DK��DL  DL~�DMHDM�HDN �DN�HDO �DO�HDP�DP��DQHDQ��DR �DR��DR��DS\DTHDT� DT�\DU\DU��DV��DV�DW}qDX �DX�HDYHDY��DZHDZ� D[ �D[�HD\ �D\��D] �D]� D^�D^��D_ �D_\D_�\D`� D`�\Da\Da�\Db� DcHDc� Dc�\Dd\Dd�De~�Df �Df��Dg  Dg~�Dg�Dh��DiHDi� Di��Dj~�Dj�\Dk�HDl�Dl��Dm  Dm\Dn  Dn��Dn��Do~�Dp  Dp��Dq�Dq�HDr  Dr\Dr�\Ds\Ds�\Dt\Du �Du��Dv �Dv��DwHDw�3DxHDx� Dx��DynDz� A�A��HA�33B�BE�Bn��B�\)B�\)B�p�B�z�B�u�B��)Cp�C�\C=qC��C)��C3�3C=� CGp�CQk�C[W
Ceu�CoxRCy�{C��qC�]qC���C�� C��fC��fC��qC��C��\C�� C��=C���C�ffC�ٚCǫ�C���Cѹ�C���C��C��C�ǮC��fC�
C���C��C��qD��D��D׮D	u�D�3DaHD޸Dk�D�Dj�D��DJ=DVfD"t{D$��D'l�D)�fD,c�D.�D1`�D3�qD6e�D8�D;[�D=�{D@^DB� DD��DG�DJaHDL�DOqHDQ�3DTH�DV�DYG
D[h�D^r�D`�3Dcr�De��Dhp�Dj�=Dmg�Do��Dq�\G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    @�
>@���A ��A z�A@z�Aa�A��\A��\A��RA���A�ffA�ffA��\A�ffB (�B33B(�B�B �B(33B0=pB8{B@(�BH{BP(�BX(�B`33Bh=pBp33Bx(�B�
>B��B�.B��B�B�{B�{B��B��B�#�B��B�
>B�{B��B�(�B��B�]B��B��B�]B�{B�#�B�=qB�34B��B�B�
>B�#�B�(�B��B��B�.C )C{C�CC�C
�C
>CC
>C{C{C
C\C�C�C{C \C"
>C$
>C&�C(�C*\C,
>C.C0C2�C4�C6
>C8�C:
>C<�C>�C@�CB
>CD\CF�CH
>CJ
>CL�CN�CP
>CR�CT�CVCXCZ�C\�C^�C`
>Cb
>Cd\Cf{Ch\Cj�ClCn�Cp\Cr
>Ct�Cv�Cx�Cz�C|�C~�C��C��C��C�fC��C��C��C�C�
=C��C��C��C�C�fC�fC��C�C�C��C�
=C��C�fC�fC��C��C�fC��C�fC�fC��C��C�fC�fC�C�fC��C�C�GC�C��C�
=C�C��C��C��C�GC��C��C�
=C�
=C��C��C��C��C��C��C��C��C�fC�fC��C��C�fC�C�fC��C��C�fC��C��C�\C��C��C�C�C�GC���C�  C�GC��C�fC��C�C�fC�fC�C�fC��C��C��C�C��C�C�fC��C�C��C��C��C��C��C��C��C�fC�fC��C��C��C�fC�fC��C��C�fC��C�fC�
=C�
=C��C��C�C��C�fC�C�  C�C��C�
=C�
=D {D �{D�D�fDfD�{D3D��D  D��D3D�{DD�{D�D�GD�D��D	GD	�GD
�D
��D3D�{D3D��D�D��D3D��D�D�D�D�3D{D��D�D��DD�3D�D�GD �D�3D�D�3D�D��D �D� D �D��D{D�{D�D�DfD��D�D�{D�D�D{D��D �D �3D! �D!��D"3D"�3D#D#�3D$�D$�3D%GD%�GD&�D&�3D'GD'��D(�D(�{D)�D)��D*{D*�{D+{D+�3D, �D,�GD-�D-��D.3D.��D/�D/�D0D0�D1�D1�3D2�D2��D3GD3��D4�D4��D5 �D5��D6GD6�GD7{D7�{D83D8��D9�D9��D:  D:��D;�D;��D<�D<�GD=3D=�{D>�D>�3D?D?�3D@{D@�fDA{DA��DB�DB��DC{DC�3DD�DD��DE�DE��DFDF��DG�DG��DH�DH�{DIDI�{DJ�DJ�{DK3DK��DL3DL��DM{DM�{DN�DN�{DO�DO�{DP�DP�DQ{DQ��DR�DR��DS�DS��DT{DT�3DU�DU��DV�DV��DWGDW��DX�DX�{DY{DY�DZ{DZ�3D[�D[�{D\�D\��D]�D]�3D^D^�D_�D_��D`�D`�3Da�Da��Db�Db�3Dc{Dc�3Dd�Dd��DeGDe��Df�Df��Dg3Dg��DhGDh��Di{Di�3Dj�Dj��Dk�Dk�{DlDl��Dm3Dm��Dn3Dn��Do�Do��Dp3Dp�Dq�Dq�{Dr3Dr��Ds�Ds��Dt�Dt��Du�Du��Dv�Dv��Dw{Dw�fDx{Dx�3Dy�DyqGDz�3AffA��A��
B=qBF  Bo�B��B��B���Bǣ�Bٞ�B�C�C��CQ�C� C)�C3ǮC=�{CG�CQ� C[k�Ce�=Co��Cy��C�ǮC�g�C��C��=C���C�ФC���C���C���C��=C��{C��
C�p�C���Cǵ�C��
C���C��C��RC��\C���C��C��HC��C��RC�ǮD��D��D��D	z�D�RDffD��Dp�D�Dp D�DO\D[�D"y�D$��D'q�D)�D,h�D.��D1e�D3�D6j�D8�D;`�D=�D@c3DB�DD��DG�DJffDL�=DOvfDQ�RDTNDV�=DYL)D[nD^w�D`�RDcw�De��Dhu�Dj�\Dml�Do��Dr{G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�?��?�>?��??�Z?	?�?�?�#?�#?0U?�?�?#:?0U?u?�?�'?#�?��?�?�X?�g?,�?�?@O??�?q?خ?�r?4n??H?|�?N�?|�?��?��?%?A ?o ?G�?A ?4�?�?��?!-??�?�r?c ?\�?!�? �?�?�)?�)?�?�?�-?c�?V?j?�?��?j�?)_?�?
�,?
��?
��?
W�?
?	�?	�C?	�^?	�~?��?�K?��?�?��?��?�0?��?��?�P?��?�k?�4?1�?1�?�?�?�?�?
=?��?�"?�?�?�'?�b?_�??4?+k?1�?4?
�?4?
�?4?YK?��?��?S&?�*>�<6>���>���>㕁>��>ӶF>���>��>��e>��m>�H�>��>���>��>�Q>��4>��j>�H�>��>���>�h
>�\�>��k>��>���>���>�E9>���>o4�>e��>Z��>W��>V>TɆ>U2a>P��>Kx>He�>?H�>4�>4��>%��>#�f>e,>	�'>�o>M>�]>%=���=�e�=��=��B=�u�=֡b=��d=ޞ=�i�=�(�=���=�&=��=�g8=��=�y�=�
�=��8=���=���=�=�0U=��	=���=�}�=��+=���=��t=��s=�L�=�Ft=��=��=��n=w��=o��=b�=J��=@�=5��='�='�="�x=	=�= �I<��E<�<�<�9X<�S<�YK<y	l<<j;�IR;��4;���<�<(�U<_<o ;���;�)_;��.;��';�o;�t�;��;��;��4:ѷ9ѷ�*dû#�
�ѷ��{�������4׻�`B��������-���P��k���p�E�r{����'������u��6z��W��ڹ������~�������R��u�&��2-�3�}�5s�3g��2�Ž9#��?H��Gy��I�̽L�ͽW���["ѽ]���\�?�]���\�?�`�o����S&�����YK��������������_p��1'���q����������qv��3���Xy�����9X�����}�����
���y��Ǯ�Ʌ��ƨ�Τ��՛=��خ��������ی�����ս�D��iD��A ��곽�˒�Fܾ
#:��M��ϫ����Q�҉�!�.�!|�$���%��&L0�'�(�U�+6z�+���-\��0��2GE�4S��4�ؾ5�7���9	l�9�Z�<���>($�?�$�@  �A�ܾD���GE9�H��K)_�L�_�N_�QN<�R�ξUϫ�W���Z6�\w��^5?�`�e�b��c�
�f�Ծg�0�i�z�kC�k��m\��p:��q&�s��us�wKǾy	l�{�6���.�������\���
���o���������+Ծ�m]��1����]���U��W��������h���r���2���{������&龒a|�����}�������㾙�о�ھ�ߤ���[���徠4n�� \���\��{J���]���/������?��RT��₾���=q������h����]̾�ƨ��5���.���.���� ��ξ����ᾴ꾴����%F�������+������Dо�y>��rG��kQ��xl�����/쾼���(���PH��]d���Ѿ�푾�Vm��!-��[����.������	��|���S���ܾ�&���MӾx�©���a��+Ծ��T�Ƨ����RT�ǆ¾ǻ0��$��>B��K^��K^��r����U��^�����ʙ1��ں��C��)_��6z��]̾˟V�ˬq���þ�M��~(��;�͑h�Ίr���?C�?��?�Z?k�?Ft?�d?˒?�I?u�?�?Ov?�?q??
J�?	rG?�k?�0?�"?�?��?_�?O?�4>��?>� �>�/�>�#:>���>gl�>Uf�>3�>	Q�=��=�D�=��=�y�=�g�=�a�=��=we�=(�U=:*<�;��.;�{�;�YK;k����4׻��$�h����}����49X�@7�["ѽ\]d�����ԕ��Xy�����c����Dg���Dо%F�+�þ4���:�ھC���L�_�Y���c�f�k���u%F��:���������D���O��6���	���澪J�������E��9X���u������PH��Ĝ���ܾA��E9��r���k�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    ?��?�>?��??�Z?	?�?�?�#?�#?0U?�?�?#:?0U?u?�?�'?#�?��?�?�X?�g?,�?�?@O??�?q?خ?�r?4n??H?|�?N�?|�?��?��?%?A ?o ?G�?A ?4�?�?��?!-??�?�r?c ?\�?!�? �?�?�)?�)?�?�?�-?c�?V?j?�?��?j�?)_?�?
�,?
��?
��?
W�?
?	�?	�C?	�^?	�~?��?�K?��?�?��?��?�0?��?��?�P?��?�k?�4?1�?1�?�?�?�?�?
=?��?�"?�?�?�'?�b?_�??4?+k?1�?4?
�?4?
�?4?YK?��?��?S&?�*>�<6>���>���>㕁>��>ӶF>���>��>��e>��m>�H�>��>���>��>�Q>��4>��j>�H�>��>���>�h
>�\�>��k>��>���>���>�E9>���>o4�>e��>Z��>W��>V>TɆ>U2a>P��>Kx>He�>?H�>4�>4��>%��>#�f>e,>	�'>�o>M>�]>%=���=�e�=��=��B=�u�=֡b=��d=ޞ=�i�=�(�=���=�&=��=�g8=��=�y�=�
�=��8=���=���=�=�0U=��	=���=�}�=��+=���=��t=��s=�L�=�Ft=��=��=��n=w��=o��=b�=J��=@�=5��='�='�="�x=	=�= �I<��E<�<�<�9X<�S<�YK<y	l<<j;�IR;��4;���<�<(�U<_<o ;���;�)_;��.;��';�o;�t�;��;��;��4:ѷ9ѷ�*dû#�
�ѷ��{�������4׻�`B��������-���P��k���p�E�r{����'������u��6z��W��ڹ������~�������R��u�&��2-�3�}�5s�3g��2�Ž9#��?H��Gy��I�̽L�ͽW���["ѽ]���\�?�]���\�?�`�o����S&�����YK��������������_p��1'���q����������qv��3���Xy�����9X�����}�����
���y��Ǯ�Ʌ��ƨ�Τ��՛=��خ��������ی�����ս�D��iD��A ��곽�˒�Fܾ
#:��M��ϫ����Q�҉�!�.�!|�$���%��&L0�'�(�U�+6z�+���-\��0��2GE�4S��4�ؾ5�7���9	l�9�Z�<���>($�?�$�@  �A�ܾD���GE9�H��K)_�L�_�N_�QN<�R�ξUϫ�W���Z6�\w��^5?�`�e�b��c�
�f�Ծg�0�i�z�kC�k��m\��p:��q&�s��us�wKǾy	l�{�6���.�������\���
���o���������+Ծ�m]��1����]���U��W��������h���r���2���{������&龒a|�����}�������㾙�о�ھ�ߤ���[���徠4n�� \���\��{J���]���/������?��RT��₾���=q������h����]̾�ƨ��5���.���.���� ��ξ����ᾴ꾴����%F�������+������Dо�y>��rG��kQ��xl�����/쾼���(���PH��]d���Ѿ�푾�Vm��!-��[����.������	��|���S���ܾ�&���MӾx�©���a��+Ծ��T�Ƨ����RT�ǆ¾ǻ0��$��>B��K^��K^��r����U��^�����ʙ1��ں��C��)_��6z��]̾˟V�ˬq���þ�M��~(��;�͑h�Ίr���?C�?��?�Z?k�?Ft?�d?˒?�I?u�?�?Ov?�?q??
J�?	rG?�k?�0?�"?�?��?_�?O?�4>��?>� �>�/�>�#:>���>gl�>Uf�>3�>	Q�=��=�D�=��=�y�=�g�=�a�=��=we�=(�U=:*<�;��.;�{�;�YK;k����4׻��$�h����}����49X�@7�["ѽ\]d�����ԕ��Xy�����c����Dg���Dо%F�+�þ4���:�ھC���L�_�Y���c�f�k���u%F��:���������D���O��6���	���澪J�������E��9X���u������PH��Ĝ���ܾA��E9��r���k�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B|B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{JBz�Bz�Bz�BzDBy�By�Bx�By>BzBz�BzDBzxBzDBz�Bz�Bz�Bz�Bz�B{Bz�Bz�Bz�Bz�B{dB{B{0B{B{B{B{B{0B{B{Bz�B{B{B{B{Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{JB}�BHB.B�UB��B��B�DB��B�B�<B�(B��B�.B�}B�HB�HB��B��B��B� B��B��B�oB��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�oB��B��B��B��B�4B�hB�}B�}B�HB�}B�}B�bB�bB�bB�HB�bB�HB�HB�HB�B�}B�HB�}B��B��B��B� B� B� B�4B�NB�:B��B��B��B��B�B�@B�[B�@B�uB�[B�uB�[B�@B��B��B��B��B��B�aB�aB�FB�{B�FB�aB�FB��B�FB�,B�,B�B�B�B�aB�,B�{B�2B�gB��B��B�9B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�?B�?B�YB�sB�sB�sB��B��B��B��B��B�+B�EB�_B�yB��B��B��B�1B��B��B��B��B��B��B�B�B�B�QB�kB��B��B��B��B��B��B��B�#B�=B�=B�=B�=B�=B�=B�WB�qB��B��B��B��B��B�)B�]B��B��B��B��B��B��B��B�B�/B�dB�dB��B��B��B��B��B��B��B��B��B��B�B�B�5B�OB�OB�OB�jB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B�B�B�B�B�!B�!B�;B�;B�VB�VB�VB��B�pB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�'B�B�'B�B�B�B�B�B�'B�'B�'B�BB�BB�BB�BB�BB�'B�'B�\B�'B�BB�BB�'B�'B�'B�B�B�B�B�'B��B��B�B�B�'B�'B�'B�BB�BB�BB�'B�'B�'B�'B�BB�BB�BB�\B�BB�BB�BB�BB�'B�BB�BB�BB�\B�BB�BB�BB�\B�BB�BB�\B�\B�\B�\B�\B�\B�BB�BB�vB�vB�vB�\B�vB�vB�vB�vB�vB�\B�\B�\B�\B�\B�\B�\B�\B�\B�BB�\B�\B�\B�BB�'B�BB�'B�BB�\B�\B�BB�BB�BB�BB�BB�\B�BB�'B�BB�'B�'B�'B�B�;B{JB{0B|B{Bz�By$BzDBzxBzBz�Bz�B{JB{�Bz�BzDBzxBz�BzBz�BzDBz�BzBy�B��B�PB��B�@B��B�TB�B� B��B�HB��B�\B�}B��B�B��B��B�{B��B��B��B��B�mB��B�?B��B�?B��B�yB��B�B��B��B�B��B�#B��B��B�/B�B�B��B�B��B�B�B�5B��B��B�B�B��B��B�VB�B��B��B��B�BB�B��B��B��B��B��B��B��B��B��B�vB�'B�\B��B�\G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{�B{BBz�Bz�Bz�Bz?By�By�Bx�By7Bz	Bz�Bz>BznBz?Bz�Bz�Bz�Bz�Bz�B{Bz�Bz�Bz�Bz�B{\B{B{*B{B{B{B{B{)B{B{Bz�B{B{B{B{Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�Bz�B{DB}�BAB%B�NB�B��B�=B��B��B�8B� B��B�(B�vB�?B�BB��B��B��B��B��B�}B�fB��B��B��B�zB�zB��B��B��B��B��B��B��B��B��B��B��B�gB��B��B��B�yB�,B�]B�uB�wB�>B�vB�vB�]B�[B�[B�AB�YB�AB�?B�AB�B�xB�?B�vB��B��B��B��B��B��B�-B�GB�3B��B��B��B��B�B�8B�TB�8B�mB�TB�oB�VB�8B��B��B��B��B��B�YB�[B�=B�rB�?B�[B�=B��B�;B�"B�#B�
B�
B�	B�[B�%B�qB�+B�aB��B��B�2B��B��B�B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B� B�6B�6B�SB�mB�lB�lB��B��B��B��B��B�%B�>B�VB�pB��B��B��B�)B��B��B��B��B��B��B��B�B�B�KB�bB�~B��B��B��B��B��B��B�B�5B�5B�5B�4B�5B�5B�OB�hB��B��B��B��B��B�"B�TB��B��B��B��B��B��B��B�B�&B�\B�[B��B��B��B��B��B��B��B��B��B��B��B��B�,B�GB�GB�GB�dB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B� B��B��B��B��B��B��B��B��B� B�B�B�5B�2B�PB�MB�RB��B�fB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�!B�B�!B�B�B�B�B�B�"B�B�"B�=B�9B�:B�=B�;B�B�"B�WB�B�:B�;B�!B�B� B�B�B�	B�B� B��B��B�B�B�"B�"B� B�;B�;B�:B� B� B�B�B�;B�:B�9B�UB�:B�<B�;B�:B� B�<B�:B�;B�TB�:B�;B�<B�VB�<B�;B�UB�UB�UB�VB�XB�TB�<B�;B�mB�mB�lB�VB�oB�oB�qB�oB�mB�VB�UB�XB�UB�VB�VB�VB�UB�VB�:B�UB�UB�UB�;B�!B�<B�B�:B�VB�VB�:B�<B�:B�:B�<B�UB�9B� B�;B� B� B� B�B�5B{JB{0B|B{Bz�By$BzDBzxBzBz�Bz�B{JB{�Bz�BzDBzxBz�BzBz�BzDBz�BzBy�B��B�PB��B�@B��B�TB�B� B��B�HB��B�\B�}B��B�B��B��B�{B��B��B��B��B�mB��B�?B��B�?B��B�yB��B�B��B��B�B��B�#B��B��B�/B�B�B��B�B��B�B�B�5B��B��B�B�B��B��B�VB�B��B��B��B�BB�B��B��B��B��B��B��B��B��B��B�vB�'B�\B��B�\G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                                    <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     PRES_ADJUSTED = PRES - Surface Pressure                                                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            Surface pressure = -0.05 dbar                                                                                                                                                                                                                                   none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Surface Pressure = -0.08 dbar                                                                                                                                                                                                                                   Not applicable                                                                                                                                                                                                                                                  Not applicable                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              Pressure adjusted in real time by using pressure offset at the sea surface                                                                                                                                                                                      No adjustment performed (values duplicated)                                                                                                                                                                                                                     No adjustment performed (values duplicated)                                                                                                                                                                                                                     202010161343572020101613435720201016134357201907260851112019072608511120190726085111IF  IF  ARFMARFMCODACODA030d030d                                                                                                                                2019072608511120190726085111                                        G�O�G�O�G�O�G�O�G�O�G�O�                                IF  IF  ARGQARGQCOQCCOQC4.3 4.3                                                                                                                                 2019072608523020190726085230QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�000000000008FB7E000000000008FB7EIF  IF  ARGQARGQCOQCCOQC4.3 4.3                                                                                                                                 2019072608523020190726085230QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�00000000000000000000000000000000GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2016V01 + ARGO climatology                                                                 20190821171504              IP      PSAL                            @�p�G�O�Dz� G�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology                                                                 20201016134357              IP      PSAL                            @�p�G�O�Dz� G�O�G�O�G�O�                                