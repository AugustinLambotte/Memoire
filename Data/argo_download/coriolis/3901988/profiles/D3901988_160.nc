CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2019-07-26T08:51:19Z creation; 2020-10-16T13:44:02Z last update (BSH ARSQ software)    
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
resolution        =���   axis      Z        X  <�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  K�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %8.2f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        X  O�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %8.2f      FORTRAN_format        F7.1   
resolution        =���     X  b�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  r<   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  �l   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ��   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  ��   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  ��   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �L   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  �$   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  �|   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.4f     FORTRAN_format        F9.3   
resolution        :�o     X  �T   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �`   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �h   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �p   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �x   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  �  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �    HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                     �$   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �D   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �L   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �T   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                     �\   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  ݬ   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T  �Argo profile    3.1 1.2 19500101000000  20190726085119  20201016134402  3901988 3901988 Argo-Norway                                                     Argo-Norway                                                     Kjell Arne Mork                                                 Kjell Arne Mork                                                 PRES            TEMP            PSAL            PRES            TEMP            PSAL               �   �AA  IFIF                                                                2C  2B  DA  APEX                            APEX                            8415                            8415                            2.10.1                          2.10.1                          846 846 @��]���}@��]���}11  @��^�`@��^�`@Q.��NP�@Q.��NP��0���cI{�0���cI{11  GPS     GPS     AA  BA  BA  Primary sampling: mixed [averaged from 4.16 dbar to 981.63 dbar {1 sample/second, 2 dbar average from 2100 dbar to 4 dbar}; discrete otherwise]                                                                                                                 Secondary sampling: discrete [1 sample, 10 dbar interval from 2100 dbar to 4 dbar]                                                                                                                                                                                    @��@�G�A ��A\)A?\)A_�A�
A��\A�{A�  A�(�A�{A�Q�AA��B�HB��B�HB��B(
=B/�B8  B?��BH
=BO��BW��B_��Bh{Bo�Bw�
B�B�B��B�
=B���B���B�
=B�B�\B�  B���B�
=B�\B�B���B���B��B���B���B���B�  B���B���B���B���B�\B�
=B���B�B�  B�
=B�C �CC  C�C
=C
�C�qC��C��C�C�CCC�C�C��C�RC!��C$  C&�C'��C*  C,
=C.  C0  C2  C3��C5�qC7�qC9��C;��C=�qC?�qCA�RCC��CF  CH�CI�qCK�qCM�qCP  CQ�qCS�qCU��CX  CZC\�C^�C`  CbCd�Ce��Cg�qCj�Cl�CnCo��Cr�Ct
=Cv�Cx�Cz�C|  C~  C�  C�HC���C��qC��qC��qC���C�HC�C�fC��C���C��)C�HC��C�C��C�C��C���C��)C�  C��C�HC���C��qC�C�C��C���C��qC��)C��)C���C���C��)C�HC�HC���C���C�  C��C�HC��C��C�HC�  C���C���C��C��C��C��C���C��)C���C��)C�  C��C�  C��)C���C���C�  C�HC��C��C��C���C�  C�  C���C��C���C��)C�  C��qC���C�HC�fC��C���C�HC�HC�HC��C���C��C��C��C��qC���C��)C��qC��qC���C��qC�HC�HC���C���C���C��)C��C�  C���C�HC��qC���C���C�  C��C��C�HC��C�  C���C��C�HC�HC���C��qC�  C��)C���C��qC�HC�D   D ~�D ��D~D��D� D�\D~�D �D~�D��D\D��D� DHD�HD�\D~�D��D	}qD	��D
�HDHD��D��D� D  D\D �D�HD  D~�D�\D\D��D��D �D��D  D��D�\D}qD��D~�D�\D\D��D~�D��D\D�\D~�D�qD|�D��D� D�\D|�D�qD~D�\D��D�\D��D �D ��D!�D!�HD!�\D"~D"�\D#��D$ �D$\D$�\D%\D%�\D&\D' �D'��D( �D(\D(�D)� D*HD*\D*�qD+\D,�D,��D-HD-\D-�\D.��D.�\D/� D/��D0}qD0��D1\D1��D2� D2�D3}qD3�D4|�D4�D5}qD5�\D6��D7�D7�{D8�D8� D8�\D9\D9�\D:�HD;�D;��D<HD<~�D<�\D=� D> �D>�HD>�\D?\D@HD@��D@��DA\DA�\DB� DCHDC�HDD �DD��DE �DE\DF  DF��DG �DG��DH  DH~�DH��DI~�DI��DJ� DJ�\DK� DK�\DL��DM�DM��DN�DN�HDO �DO\DP �DP� DQ  DQ��DR�DR�HDS�DS� DT �DT��DU  DU� DV  DV~�DW  DW~�DW�DX\DX�\DY��DZ  DZ��D[�D[��D\HD\\D\��D]�HD]��D^}qD^�\D_� D_�\D`\DaHDa\Db  Db\Db�\Dc\Dc��Dd� De  De\De��Df~Dg �Dg�HDg��Dh~�Di �Di� Di�\Dj~�Dk  Dk�HDlHDl��Dl�Dm~�Dn�Dn�HDn�\Do~Dp �Dp�HDq �Dq��Dr �Dr~�Dr��Ds~�Ds��Dt}qDu �DuhRDv�RA\)A��
A�BG�BD�
Bf33B�\)B��fB�.B��
B��B�L�B��qC�{C� CC)��C3��C=\)CExRCN�C[ffCe��Co0�Cy��C��qC�ٚC��=C�g�C���C��C���C�޸C���C�ФC�C���C��)C¬�CǥC̨�C��RC֯\Cۗ
C��fC�qC��HC�L�C���C��qC���DmqD��D�\D	Z�D��D^�D�qDQ�D޸Dl)D׮D<)D�D"G
D$��D'S�D)ٚD,K�D.� D1k�D3׮D6@ D8�qD:� D=ۅD@QHDB�{DEL)DG�fDJa�DLθDO\)DQФDT^DV� DY1HD[�D]�D`��Dce�De��DhQHDj�RDl�HG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              @��G@�
=@�
=A=qA>=qA^�\A~�RA�  A��A�p�A���AυA�A�
>A�\)B��B�B��B�B'B/��B7�RB?�BGBO�BW�B_�Bg��Bo��Bw�\B��B��HB���B��fB���B���B��fB��HB��B��)B���B��fB��B��HB���B��
B���B��
B��
B��
B��)B���B��
B���B���B��B��fB��
B��HB��)B��fB��HB��GC�3C�C��C�QC	��C�C��C��C�C�C�3C�3C��C��C��C�fC!��C#�C%�C'��C)�C+�QC-�C/�C1�C3��C5�C7�C9��C;��C=�C?�CA�fCC��CE�CG�CI�CK�CM�CO�CQ�CS�CU��CW�CY�3C[�C]�C_�Ca�3Cc�Ce��Cg�Ci�Ck��Cm�3Co��Cq�Cs�QCu��Cw��Cy�C{�C}�C�C��RC���C��{C��{C��{C���C��RC��)C��pC���C���C��3C��RC���C��)C���C��)C���C���C��3C��
C���C��RC���C��{C��)C��)C���C���C��{C��3C��3C���C���C��3C��RC��RC���C���C��
C���C��RC���C���C��RC��
C���C���C���C���C���C���C���C��3C���C��3C��
C���C��
C��3C���C���C��
C��RC���C���C���C���C��
C��
C���C���C���C��3C��
C��{C���C��RC��pC���C���C��RC��RC��RC���C���C���C���C���C��{C���C��3C��{C��{C���C��{C��RC��RC���C���C��C��3C���C��
C���C��RC��{C��C���C��
C���C���C��RC���C��
C���C���C��RC��RC���C��{C��
C��3C���C��{C��RC��)C��
D z=D �=Dy�D�=D{�D��Dz=D�)Dz=D�=Dz�D�=D{�D��D|�D��Dz=D�=D	x�D	�=D
|�D
��D|)D�=D{�D��Dz�D�)D|�D��Dz=D��Dz�D�=D|)D�)D|)D��D|)D��Dx�D�=Dz=D��Dz�D�=Dz=D�=Dz�D��Dz=D��DxRD�=D{�D��DxRD��Dy�D��D|)D��D|)D�qD |)D �qD!|�D!��D"y�D"��D#}qD#�)D$z�D$��D%z�D%��D&z�D&�)D'~D'�)D(z�D(��D){�D)��D*z�D*��D+z�D+�D,}qD,��D-z�D-��D.|)D.��D/{�D/�=D0x�D0�=D1z�D1�=D2{�D2��D3x�D3��D4xRD4��D5x�D5��D6|)D6�qD7� D7�D8{�D8��D9z�D9��D:|�D:�D;~D;��D<z=D<��D={�D=�)D>|�D>��D?z�D?��D@|)D@�=DAz�DA��DB{�DB��DC|�DC�)DD|)DD�)DEz�DE��DF|)DF�)DG|)DG��DHz=DH�=DIz=DI�=DJ{�DJ��DK{�DK��DL|)DL�qDM|)DM�qDN|�DN�)DOz�DO�)DP{�DP��DQ~DQ�qDR|�DR�DS{�DS�)DT}qDT��DU{�DU��DVz=DV��DWz=DW��DXz�DX��DY|)DY��DZ|)DZ�qD[}qD[��D\z�D\�=D]|�D]�=D^x�D^��D_{�D_��D`z�D`��Daz�Da��Dbz�Db��Dcz�Dc�=Dd{�Dd��Dez�De�=Dfy�Df�)Dg|�Dg�=Dhz=Dh�)Di{�Di��Djz=Dj��Dk|�Dk��Dl|)Dl��Dmz=Dm�Dn|�Dn��Doy�Do�)Dp|�Dp�)Dq}qDq�)Drz=Dr�=Dsz=Ds�=Dtx�Dt�)Duc�Dv��A  A�(�A�{Bp�BD  Be\)B��B�z�B�B�k�Bـ B��HB�Q�C^�CJ=C��C)� C3�C=&fCEB�CN�RC[0�Ces3Cn��Cy�C���C���C��\C�L�C��C���C���C���C��C���C���C���C��HC�CǊ=C̎CѽqC֔{C�|)CૅC傏C�fC�1�C���C���C�q�D` D�{D��D	MqD�{DQHD� DD{D�HD^�D�=D.�DФD"9�D$�{D'FfD)�)D,>D.D1^D3�=D62�D8� D:��D=�D@C�DB�
DE>�DG��DJT{DL�HDON�DQ�3DTP�DVr�DY#�D[�=D]�D`�DcXRDe�{DhC�Dj��Dl��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�?�/�?�J#?��8?��?�z?�??�;d?��f?��E?��?��]?�$�?��M?��?�O?���?ή}?ɏ�?���?��'?��?�($?�&�?���?�a|?yDg?X��?Co?:C�?8��?7��?6YK?6��?4!?0[�?-��?%?��?	�^?�u?
��?/�?��?	�Z>�5?>���>��v>��>�A >��;>�N<>��">��}>�>��Y>�kQ>��>>�s>��4>�=q>�C�>�1'>��>��>��]>��>�S&>���>��>��0>��>���>�#:>��B>~($>y��>|�$>���>��M>�>���>��>�a|>��>�ϫ>׍P>��>�a>C�&>#n/>
ں=�xl=���>,�>��>@O>��>��>!�S>*�h>2�X>`Ĝ>m\�>n.�>m(�>q�3>z^5>��7>s�}>z��>w1�>{�>z��>s�F>s33>w�k>w1�>x�5>|�>~v�>�o >�2�>��]>���>��>��>��>��L>���>���>��1>�^�>��^>�^�>�r�>�>B>�0U>��>�~(>��>� i>�Ɇ>��}>��>�%�>�K^>���>�(>�P�>���>���>���>�}V>���>��r>�!�>��M>�\�>��>�<�>���>�j�>�6z>�P�>�x>���>���>���>�"h>�<�>�/�>���>�1>�q>���>�P�>�(>���>�*0>�l�>�+>�8>�_p>�$�>�?>�>B>��>���>���>�>�>�j�>�ں>�*0>��]>���>���>��S>��I>�˒>�5?>�PH>�}�>��>���>��#>��+>��b>�Y�>�Z�>��>�Ɇ>�Y>�s>���>�Y�>��N>�C�>��>�$>�1'>�	�>��>���>�n/>��x>��e>��.>�A�>v�}>hی>fff>b�>`��>[��>W��>V�+>T��>P��>Mq>D2�>;/�>0�|>,"h>&ff>%F>!a�>�>��>��>�*>M>1'>�>�\=��=��"=�:�=�(�=� �=���=��=�W�=�E9=�g8=�6=��[=�K�=��=�H�=��z=�Xy=���=��L=�	=�(=�_p=��=�M=�%=L��=HK^=Em]=?�[=:�=7�=5s�=2-=%�= 'R=\)<��><�T�<�w�<��N<��r<r{�<D��<*d�<(�U<�;�T�;�IR;�$;K)_;#�
:��4:�҉:�o�ѷ���'�����r�2��>��c�Ӽ�\)��+������zx��#�����������>��/��G������Ӽѷ�͞���҉��s��Z���J��YK�������"3��ATʽLc�N��Np;�N��@��JW��]���X���`�m(��f�Խo���l���lVֽo4׽jJ��t���we��y	l�v+k�|�$�~\������%�r{��n���d%��ZQ�]/�f1��ixսkC�k��f�Խe��lVֽ�+���.���3��Y��,��՛=��j��"ѽ�/��ᱽ��۾	޾	k��	�̾	޾	�'�
��5������녾��%F�*0U�,V־..��-�D�.c �6�F�D�ݾG�K�L��O��V�b�Xy>�\(��_��ezx�i��i^��l<��o��oO�r�X�q���tmƾxl"�|�[��'R��%��@���&����x��a��a��9�������S����Ӿ�������}V���ξ�.����W���3���N��{���L������X��֡��ϫ�������"���~���P��J����~���쾙e,���о�c��J���=پ�X�������~��c��J����쾛�Q���������6���$���Ѿ�6������㾝�������۾��������}���p���}���}����H������}����d��曾�˒?��?�x?���?��?�q�?G�+?71�?J#?W?>��;>�9X>��4>��o>�ی>�Dg>�I�>��o>E��=�~�>�t>hXy>o��>tm�>z��>�Xy>���>���>�"h>��b>�)_>�.�>���>�6z>�I�>�q>��>���>�k�>�#:>���>�<6>�>�K�>�~�>���>��I>bh
>U��>,�>�Q>E9=�A =��=�U2=���=��n=@7=+6z<��.<9#�;r{�:�o�,1��w����>��s��xսIR�O�ν\�?�n���v�F���&�Y�>�k������W?�	�^�~(�%���-\��N��eF�p���}p���h
������������2a���������J���C-��/G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              ?�/�?�J#?��8?��?�z?�??�;d?��f?��E?��?��]?�$�?��M?��?�O?���?ή}?ɏ�?���?��'?��?�($?�&�?���?�a|?yDg?X��?Co?:C�?8��?7��?6YK?6��?4!?0[�?-��?%?��?	�^?�u?
��?/�?��?	�Z>�5?>���>��v>��>�A >��;>�N<>��">��}>�>��Y>�kQ>��>>�s>��4>�=q>�C�>�1'>��>��>��]>��>�S&>���>��>��0>��>���>�#:>��B>~($>y��>|�$>���>��M>�>���>��>�a|>��>�ϫ>׍P>��>�a>C�&>#n/>
ں=�xl=���>,�>��>@O>��>��>!�S>*�h>2�X>`Ĝ>m\�>n.�>m(�>q�3>z^5>��7>s�}>z��>w1�>{�>z��>s�F>s33>w�k>w1�>x�5>|�>~v�>�o >�2�>��]>���>��>��>��>��L>���>���>��1>�^�>��^>�^�>�r�>�>B>�0U>��>�~(>��>� i>�Ɇ>��}>��>�%�>�K^>���>�(>�P�>���>���>���>�}V>���>��r>�!�>��M>�\�>��>�<�>���>�j�>�6z>�P�>�x>���>���>���>�"h>�<�>�/�>���>�1>�q>���>�P�>�(>���>�*0>�l�>�+>�8>�_p>�$�>�?>�>B>��>���>���>�>�>�j�>�ں>�*0>��]>���>���>��S>��I>�˒>�5?>�PH>�}�>��>���>��#>��+>��b>�Y�>�Z�>��>�Ɇ>�Y>�s>���>�Y�>��N>�C�>��>�$>�1'>�	�>��>���>�n/>��x>��e>��.>�A�>v�}>hی>fff>b�>`��>[��>W��>V�+>T��>P��>Mq>D2�>;/�>0�|>,"h>&ff>%F>!a�>�>��>��>�*>M>1'>�>�\=��=��"=�:�=�(�=� �=���=��=�W�=�E9=�g8=�6=��[=�K�=��=�H�=��z=�Xy=���=��L=�	=�(=�_p=��=�M=�%=L��=HK^=Em]=?�[=:�=7�=5s�=2-=%�= 'R=\)<��><�T�<�w�<��N<��r<r{�<D��<*d�<(�U<�;�T�;�IR;�$;K)_;#�
:��4:�҉:�o�ѷ���'�����r�2��>��c�Ӽ�\)��+������zx��#�����������>��/��G������Ӽѷ�͞���҉��s��Z���J��YK�������"3��ATʽLc�N��Np;�N��@��JW��]���X���`�m(��f�Խo���l���lVֽo4׽jJ��t���we��y	l�v+k�|�$�~\������%�r{��n���d%��ZQ�]/�f1��ixսkC�k��f�Խe��lVֽ�+���.���3��Y��,��՛=��j��"ѽ�/��ᱽ��۾	޾	k��	�̾	޾	�'�
��5������녾��%F�*0U�,V־..��-�D�.c �6�F�D�ݾG�K�L��O��V�b�Xy>�\(��_��ezx�i��i^��l<��o��oO�r�X�q���tmƾxl"�|�[��'R��%��@���&����x��a��a��9�������S����Ӿ�������}V���ξ�.����W���3���N��{���L������X��֡G�O�G�O����"���~���P��J����~���쾙e,���о�c��J���=پ�X�������~��c��J����쾛�Q���������6���$���Ѿ�6������㾝�������۾��������}���p���}���}����H������}����d��曾�˒?��?�x?���?��?�q�?G�+?71�?J#?W?>��;>�9X>��4>��o>�ی>�Dg>�I�>��o>E��=�~�>�t>hXy>o��>tm�>z��>�Xy>���>���>�"h>��b>�)_>�.�>���>�6z>�I�>�q>��>���>�k�>�#:>���>�<6>�>�K�>�~�>���>��I>bh
>U��>,�>�Q>E9=�A =��=�U2=���=��n=@7=+6z<��.<9#�;r{�:�o�,1��w����>��s��xսIR�O�ν\�?�n���v�F���&�Y�>�k������W?�	�^�~(�%���-\��N��eF�p���}p���h
������������2a���������J���C-��/G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B
ߊB
ߊB
�pB
�5B
�B
��B
�B
�B
��B
�B
�6B
�B
��B
�B
�xB�B3�B+�B-�B&�B$�B"hB$�B1B�B�B�B
�B�B�B�BKBfB
�B�BB	�B	BmBBmB�B	�B	�BDB�BHB�B4BNB4B:B�B�B@BaBB�B�B�B#BIB;B!B$�B(
B)�B*�B,�B-]B.}B/�B0�B1[B3MB6�B=<BA�BF?BIBK�BNVBQ�BW�B[=B]~B_pBY�BU�BU�BT{BT,BT�BUgBWsBX�BYBZ7B[�B]�B]�Ba�BdZBd�Bd�Be,Be�BgBf�Bg�Bg�Bh>Bh�BhXBhsBh�Bi�BjKBj�BkBk�BlBl�BncBn}Bo Bo�Bo�Bp!BpUBp�BqvBq�BrBrBrBr|Bs3Bs�BtBt�Bu?BvzBw2BxRBy�BzxBz�B{0B{0B{0B{�B|PB|PB|PB|�B|�B|�B|�B|�B}qB}�B}�B}�B~B~BB~BB~(B~wB~�B~�B~�B~�B�B�iB��B��B��B��B��B��B��B��B��B�oB��B�aB�gB�9B��B�B�?B�B��B��B�B��B��B�7B��B�#B�#B�rB��B�)B�xB��B��B�JB�PB��B�B�B�(B�HB�hB��B��B�B�aB�FB�FB�,B�,B�{B��B��B��B�FB�gB�MB�gB�MB��B��B��B��B��B��B��B�B��B��B��B�$B��B�YB�$B�
B�?B��B��B��B��B��B�B�B��B��B��B��B��B�yB��B��B��B��B��B��B��B��B��B�B�B��B�1B�1B�KB��B�KB��B�B�B�B�B��B�B��B��B�B��B�B�7B�QB�7B�7B�kB�QB�kB�QB�kB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�#B�	B�	B��B�=B�	B�=B�#B�#B�=B�=B�WB��B�=B�	B�WB�	B�	B�#B��B�qB�=B�qB�WB�B��B��B��B��B��B�)B��B��B��B�CB�CB�]B�CB�)B�)B�B�CB�CB�]B�)B�]B�CB�]B�)B�xB��B�B�)B�]B�]B�CB�]B�]B�)B��B��B��B��B�B�CB�~B�B��B�dB��B�dB�dB�~B�/B�/B�IB�IB�IB�IB�dB�IB��B��B��B��B��B��B��B�dB�B�5B�B�jB�B�jB�B�OB�OB��B��B�OB��B��B�5B��B�jB��B�jB��B��B��B��B��B��B��B��B��B��B��B��B��B�!B��B��B��B�!B��B��B��B�;B��B�B��B��B
�VB��B�!B�B�!B�B�B�B�B�B�B�B�B�B��B��B��B��B�B��B��B��B�B��B��B��B��B�!B�B��B��B��B��B��B��B��B��B��B��B��B��B��B
�pB
�tB
��B:�BIB
rB	B�B�B�BBB�B 'B+�B0�BF?BX�BUgBUMBX�BbNBh$Bg�Bj�Bn�BpUBp�BtBv�Bz*B{B{�B}�B~B�B��B�B��B�B��B�DB�xB�pB�@B��B��B��B��B�YB�mB��B��B��B�1B�qB��B�B��B��B��B��B�=B��B��B�#B��B�=B��B�B�#B�/B�]B��B�B�B��B��B��B�IB�5B�IB��B��B��B��B�B�pB�;B��B�B��B�'B�jB�B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              B
ߗB
ߘB
�}B
�BB
�B
�B
�B
��B
��B
�B
�CB
�B
�B
�B
��B�B3�B+�B-�B&�B$�B"uB%B>B�B	B�B
�B�B
B�BXBtB
�B�B'B	�B	*B{BB{B�B	�B	�BNB�BQB�B@B[B@BEB�B�BKBmB)B�B�B�B3BWBHB! B$�B(B)�B*�B,�B-jB.�B/�B0�B1hB3VB6�B=JBA�BFIBI(BLBNbBQ�BW�B[KB]�B_|BY�BU�BU�BT�BT9BT�BUrBW}BX�BY�BZCB[�B]�B]�BbBdbBd�Bd�Be7Be�Bg+Bf�Bg�Bg�BhIBh�BhdBh�BiBi�BjXBj�BkBk�BlBl�BnlBn�BoBo�Bo�Bp,BpaBp�BqBq�BrBrBrBr�Bs>Bs�Bt+Bt�BuLBv�Bw>Bx\By�Bz�Bz�B{:B{<B{=B{�B|YB|[B|ZB|�B|�B|�B|�B|�B}|B}�B}�B~B~B~NB~NB~4B~�B~�B~�BBB�B�vB��B��B��B��B��B��B��B��B��B�|B��B�lB�tB�BB��B�B�OB�B��B��B�!B��B��B�BB��B�.B�.B�}B�B�4B��B��B�B�XB�ZB��B�B�B�3B�VB�sB��B�B�B�kB�NB�NB�9B�9B��B��B��B��B�OB�rB�ZB�sB�ZB��B��B��B��B��B��B��B�B��B��B��B�/B��B�cB�,B�B�IB��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B� B� B��B�<B�;B�WB��B�SB��B��B��B��B��B��B��B��B��B��B��B�(B�CB�_B�CB�EB�wB�^B�uB�^B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B��B�-B�B�B��B�HB�B�GB�,B�.B�IB�GB�aB��B�FB�B�aB�B�B�.B��B�|B�IB�}B�aB�B��B��B��B��B��B�4B��B��B��B�NB�NB�iB�MB�4B�4B�B�MB�OB�iB�4B�gB�NB�iB�2B��B��B�B�3B�fB�gB�OB�gB�fB�5B��B��B��B��B�!B�NB��B�B��B�oB��B�oB�oB��B�8B�;B�UB�SB�SB�SB�nB�SB�B��B��B��B��B��B��B�oB�%B�>B�B�tB�'B�vB�&B�YB�YB��B��B�[B��B��B�AB��B�uB��B�wB��B��B��B��B��B��B��B��B��B��B��B��B��B�+B��B��B��B�+B��B��B��B�EB��B�B��G�O�G�O�B��B�,B�B�+B�B�B�B�B�B�B�B�B�B��B��B��B��B�B��B��B��B�B��B��B��B��B�+B�B��B��B��B��B��B��B��B��B��B��B��B��B��B
�pB
�tB
��B:�BIB
rB	B�B�B�BBB�B 'B+�B0�BF?BX�BUgBUMBX�BbNBh$Bg�Bj�Bn�BpUBp�BtBv�Bz*B{B{�B}�B~B�B��B�B��B�B��B�DB�xB�pB�@B��B��B��B��B�YB�mB��B��B��B�1B�qB��B�B��B��B��B��B�=B��B��B�#B��B�=B��B�B�#B�/B�]B��B�B�B��B��B��B�IB�5B�IB��B��B��B��B�B�pB�;B��B�B��B�'B�jB�B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111441111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111                                                                                                                                                                                                                                                                                                                                                                                                              <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     PRES_ADJUSTED = PRES - Surface Pressure                                                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL                                                                                                                                                                                                                                            Surface pressure = 0.07 dbar                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Surface Pressure = 0.21 dbar                                                                                                                                                                                                                                    Not applicable                                                                                                                                                                                                                                                  Not applicable                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              Pressure adjusted in real time by using pressure offset at the sea surface                                                                                                                                                                                      No adjustment performed (values duplicated)                                                                                                                                                                                                                     No adjustment performed (values duplicated)                                                                                                                                                                                                                     202010161344022020101613440220201016134402201907260851192019072608511920190726085119IF  IF  ARFMARFMCODACODA030d030d                                                                                                                                2019072608511920190726085119                                        G�O�G�O�G�O�G�O�G�O�G�O�                                IF  IF  ARGQARGQCOQCCOQC4.3 4.3                                                                                                                                 2019072608532920190726085329QCP$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�000000000008FB7E000000000008FB7EIF  IF  ARGQARGQCOQCCOQC4.3 4.3                                                                                                                                 2019072608532920190726085329QCF$QCF$                                G�O�G�O�G�O�G�O�G�O�G�O�00000000000040000000000000000000GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2016V01 + ARGO climatology                                                                 20190821171509              IP      PSAL                            @��G�O�Dv�RG�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO climatology                                                                 20201016134402              IP      PSAL                            @��G�O�Dv�RG�O�G�O�G�O�                                