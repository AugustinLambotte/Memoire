CDF      
      	DATE_TIME         STRING2       STRING4       STRING8       STRING16      STRING32       STRING64   @   	STRING256         N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY                     <   	DATA_TYPE                  comment       	Data type      
_FillValue                    0�   FORMAT_VERSION                 comment       File format version    
_FillValue                    0�   HANDBOOK_VERSION               comment       Data handbook version      
_FillValue                    0�   REFERENCE_DATE_TIME                 comment       !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    1    PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    1   PROJECT_NAME                  comment       Name of the project    
_FillValue                  @  1   PI_NAME                   comment       "Name of the principal investigator     
_FillValue                  @  1X   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  1�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       <0..N, 0 : launch cycle (if exists), 1 : first complete cycle   
_FillValue         ��        1�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    1�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    1�   DATE_CREATION                   comment       Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    1�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    1�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     1�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    2   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    2   INST_REFERENCE                    	long_name         Instrument type    conventions       Brand, type, serial number     
_FillValue                  @  2   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    2\   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            2`   JULD_QC                	long_name         Quality on Date and Time   conventions       Argo reference table 2     
_FillValue                    2h   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~            2l   LATITUDE               	long_name         &Latitude of the station, best estimate     units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�             2t   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�             2|   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    2�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    2�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    2�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    2�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    2�   PRES         
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  2�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  :P   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  <@   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  C�   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  E�   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  M�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  UL   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  W<   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ^�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  `�   PSAL         
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  h�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  pH   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  r8   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  y�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  {�   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   CALIBRATION_DATE            	             
_FillValue                  ,  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �<   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �L   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �P   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �`   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �d   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �h   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �lArgo profile    2.2 1.2 19500101000000  4901407 Irminger Sea, Argo equivalent                                   Amy Bower                                                       PRES            TEMP            PSAL              A   AO  20121107090145  20130503154401  4476_5267_260                   2C  D   APEXIR_SBE_3479                                                 846 @�ia  1   @�i5� @O�n��P�;�7KƧ�1   GPS     A   A   A   @��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���C�fC  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DXy�DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dmy�Dm��Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� DvfDv��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @&ff@�ff@�ffA33A#33AC33Ac33A���A���A���A���A���Aљ�AᙚA�B ��B��B��B��B ��B(��B0��B8��B@��BH��BP��BX��B`��Bh��Bp��Bx��B�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffC �C�C33C33C33C
33C33C33C33C33C33C33C33C33C33C33C 33C"33C$33C&33C(33C*33C,33C.33C033C233C433C633C833C:33C<33C>33C@33CB33CD33CF33CH33CJ33CL33CN33CP33CR33CT33CV33CX33CZ33C\33C^33C`33Cb33Cd33Cf33Ch33Cj33Cl33Cn33Cp33Cr33Ct33Cv33Cx33Cz33C|33C~33C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV�DV��DW�DW��DX�DX�fDY�DY��DZ�DZ��D[�D[��D\�D\��D]�D]��D^�D^��D_�D_��D`�D`��Da�Da��Db�Db��Dc�Dc��Dd�Dd��De�De��Df�Df��Dg�Dg��Dh�Dh��Di�Di��Dj�Dj��Dk�Dk��Dl�Dl��Dm�Dm�fDnfDn��Do�Do��Dp�Dp��Dq�Dq��Dr�Dr��Ds�Ds��Dt�Dt��Du�Du��Dv3Dv�f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A	hsA	`BA	hsA	p�A	|�A	��A	��A	�A	A	��A	A	��A	��A	ƨA	��A	A	�^A	�^A	�FA	��A	��A	��A	�A	�-A	�wA	�
A	�
A	�TA
  A
9XA
A�A
A�A
A�A
A�A
9XA
5?A
5?A
5?A
5?A
1'A
1'A
1'A
5?A
1'A
1'A
-A
$�A��A Q�@�1'@�S�@��#@��@�G�@��@�bN@�;d@�@���@�@�hs@��@��@�?}@�j@띲@�^5@�@�  @��@�b@�F@�+@��T@��@��@��@�j@�x�@�j@�(�@�+@��@�ȴ@�
=@�@��y@��@���@�R@�!@���@���@�7L@�`B@��T@��T@��#@��@���@��@�^@�x�@�p�@�O�@�?}@�?}@�&�@�7L@���@�bN@�I�@�A�@�z�@���@��y@��@�w@�|�@�l�@�V@�O�@�&�@��`@�9@�@�bN@�z�@��@���@���@���@�%@�/@�hs@�@�-@�@��T@���@�@�@�@�-@�-@�7@�/@��@�V@���@�j@�Q�@�b@�w@�@�|�@�\)@�P@���@�Q�@�j@��@�u@蛦@�u@�bN@�b@�b@��m@睲@�P@�dZ@�o@�\@��@���@�-@�@�@�%@���@�j@�D@� �@�b@��@�o@�J@�^@�7@�z�@�ƨ@�l�@߶F@�(�@� �@�ƨ@�o@ް!@�J@ݡ�@�p�@�&�@�Ĝ@��@�{@ف@�X@�&�@��@�Ĝ@��/@���@���@���@�~�@�p�@�&�@��@�V@�bN@��@��
@Ӯ@�\)@�
=@ҏ\@�hs@�%@�x�@��T@��@�@�p�@���@�j@�I�@��
@Ώ\@�=q@�X@���@�t�@ˮ@˝�@�ȴ@��@�7L@Ȭ@�j@�r�@ȓu@ȃ@�(�@���@ǝ�@�(�@�hs@ȋD@�+@őh@Ł@�-@��@�+@�K�@�K�@�S�@�C�@�33@���@��y@�"�@�o@ƸR@�E�@�=q@��#@�hs@��@ă@��@Õ�@�;d@��@¸R@¸R@¸R@°!@+@�ff@�=q@��@�@�J@���@��7@�Ĝ@��`@��/@�Z@�\)@�\)@�C�@�o@��@��#@�p�@�V@�b@�33@�
=@��!@��h@��@�&�@�%@�Z@���@��P@�@�V@�E�@�=q@�5?@���@�x�@�O�@���@��D@� �@���@�;d@��!@�J@�@���@�G�@��`@���@�bN@� �@��F@��@�S�@�
=@�^5@�x�@��@��@��@���@�|�@��@�ff@��@��T@��@��-@��@�7L@��9@�9X@�(�@�1@���@�t�@�dZ@�K�@�C�@�;d@��@���@�~�@�E�@�-@��@��^@��h@�X@���@��@��@�(�@�b@�  @��;@��@�l�@�C�@��H@�ȴ@��\@�E�@�J@��@��@��^@���@�x�@�G�@��@�V@�%@��9@�Z@�I�@�1'@�b@�  @��
@���@�K�@�+@�
=@�@���@��+@�$�@��T@���@�`B@�/@���@�Ĝ@�r�@�bN@�r�@�9X@�1@��;@��w@��P@�l�@�K�@��H@��+@�E�@�-@�$�@��@��@�hs@�O�@�/@��@��@���@���@���@��@��u@�Q�@��@�b@���@�;d@��@�o@�o@�@��@���@�n�@�E�@��-@��@�%@��/@��j@�1'@��;@��w@��F@��@�|�@�+@��H@��!@��+@�^5@��@��@��#@��^@�p�@�?}@�/@�V@���@���@��@�Z1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A	hsA	`BA	hsA	p�A	|�A	��A	��A	�A	A	��A	A	��A	��A	ƨA	��A	A	�^A	�^A	�FA	��A	��A	��A	�A	�-A	�wA	�
A	�
A	�TA
  A
9XA
A�A
A�A
A�A
A�A
9XA
5?A
5?A
5?A
5?A
1'A
1'A
1'A
5?A
1'A
1'A
-A
$�A��A Q�@�1'@�S�@��#@��@�G�@��@�bN@�;d@�@���@�@�hs@��@��@�?}@�j@띲@�^5@�@�  @��@�b@�F@�+@��T@��@��@��@�j@�x�@�j@�(�@�+@��@�ȴ@�
=@�@��y@��@���@�R@�!@���@���@�7L@�`B@��T@��T@��#@��@���@��@�^@�x�@�p�@�O�@�?}@�?}@�&�@�7L@���@�bN@�I�@�A�@�z�@���@��y@��@�w@�|�@�l�@�V@�O�@�&�@��`@�9@�@�bN@�z�@��@���@���@���@�%@�/@�hs@�@�-@�@��T@���@�@�@�@�-@�-@�7@�/@��@�V@���@�j@�Q�@�b@�w@�@�|�@�\)@�P@���@�Q�@�j@��@�u@蛦@�u@�bN@�b@�b@��m@睲@�P@�dZ@�o@�\@��@���@�-@�@�@�%@���@�j@�D@� �@�b@��@�o@�J@�^@�7@�z�@�ƨ@�l�@߶F@�(�@� �@�ƨ@�o@ް!@�J@ݡ�@�p�@�&�@�Ĝ@��@�{@ف@�X@�&�@��@�Ĝ@��/@���@���@���@�~�@�p�@�&�@��@�V@�bN@��@��
@Ӯ@�\)@�
=@ҏ\@�hs@�%@�x�@��T@��@�@�p�@���@�j@�I�@��
@Ώ\@�=q@�X@���@�t�@ˮ@˝�@�ȴ@��@�7L@Ȭ@�j@�r�@ȓu@ȃ@�(�@���@ǝ�@�(�@�hs@ȋD@�+@őh@Ł@�-@��@�+@�K�@�K�@�S�@�C�@�33@���@��y@�"�@�o@ƸR@�E�@�=q@��#@�hs@��@ă@��@Õ�@�;d@��@¸R@¸R@¸R@°!@+@�ff@�=q@��@�@�J@���@��7@�Ĝ@��`@��/@�Z@�\)@�\)@�C�@�o@��@��#@�p�@�V@�b@�33@�
=@��!@��h@��@�&�@�%@�Z@���@��P@�@�V@�E�@�=q@�5?@���@�x�@�O�@���@��D@� �@���@�;d@��!@�J@�@���@�G�@��`@���@�bN@� �@��F@��@�S�@�
=@�^5@�x�@��@��@��@���@�|�@��@�ff@��@��T@��@��-@��@�7L@��9@�9X@�(�@�1@���@�t�@�dZ@�K�@�C�@�;d@��@���@�~�@�E�@�-@��@��^@��h@�X@���@��@��@�(�@�b@�  @��;@��@�l�@�C�@��H@�ȴ@��\@�E�@�J@��@��@��^@���@�x�@�G�@��@�V@�%@��9@�Z@�I�@�1'@�b@�  @��
@���@�K�@�+@�
=@�@���@��+@�$�@��T@���@�`B@�/@���@�Ĝ@�r�@�bN@�r�@�9X@�1@��;@��w@��P@�l�@�K�@��H@��+@�E�@�-@�$�@��@��@�hs@�O�@�/@��@��@���@���@���@��@��u@�Q�@��@�b@���@�;d@��@�o@�o@�@��@���@�n�@�E�@��-@��@�%@��/@��j@�1'@��;@��w@��F@��@�|�@�+@��H@��!@��+@�^5@��@��@��#@��^@�p�@�?}@�/@�V@���@���@��@�Z1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B  BB%B%B%B%B%B%B+B+B+B+B+B+B+B%B%BBB��BVB�B�B�B�B�B�B{BoBPB+BBBB%B	7B1BBB��B��B  BBBBBB1BB%BPBJBVBDBDBPBbBuB�B�B�B�B�B �B-B1'B2-B7LB8RB:^B;dB<jB<jB>wB>wB?}B?}BA�BC�BE�BG�BG�BE�BE�BG�BJ�BN�B]/BbNBdZBdZBffBr�Bt�Bs�Br�Bq�Bq�Bq�Br�Bu�Bv�Bw�Bw�Bx�By�B|�B}�B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B�B�B�%B�1B�7B�DB�JB�PB�PB�PB�JB�PB�PB�JB�JB�JB�DB�=B�1B�%B�B�%B�B�B�B�B�B�B�B~�B{�Bx�Bw�Bv�Br�Bp�Bo�Br�Bt�Bu�Bs�Bq�Bo�Bn�Bm�Bl�Bk�BhsBdZB]/BZBZBZBZB\)B]/BZBW
BVBS�BQ�BQ�BQ�BR�BR�BS�BR�BR�BQ�BP�BN�BL�BK�BO�BR�BR�BR�BR�BP�BN�BM�BK�BG�BE�BD�BC�B>wB?}B@�B<jB9XB7LB5?B5?B6FB9XB:^B9XB9XB9XB>wBD�B?}B8RB33B49B9XB?}B@�BA�BA�BB�BC�BC�BD�BD�BF�BG�BF�BF�BH�BH�BG�BG�BF�BE�BE�BD�BC�BC�BC�BC�BC�BD�BD�BC�BC�BD�BF�BG�BF�BC�BE�BE�BD�BA�BB�BC�BC�BA�B>wB<jB;dB8RB5?B49B33B0!B.B0!B1'B.B-B,B)�B)�B(�B(�B(�B'�B'�B&�B%�B%�B%�B$�B#�B!�B �B�B�B�B�B�B�B�B�B�B�B�B�B�BuBuBuBuBuBoBbB\BbBbBbBbB\BVBPBPBJBDB
=B
=B
=B
=B
=B
=B	7B	7B	7B	7B	7B	7B	7B1B1B1B1B1B1B+B+B+B+B%B%B%B%B%B%B%BBBBBBBBBBBBBBBBBBBB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�B�yB�yB�yB�sB�sB�sB�sB�sB�mB�mB�s1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B��B��B��B��B��B��B��B��B��B��B�-B��B�VB��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B�cB�B$B#B%B=B/B%B&B)B6B)B'B B3B)B3BvB�B/BIB�B�B�B`B7BgBEB�BB
eB�B�B�B�B
jB	gBB�B� B��B B�B�BBB�B	(BpBBB�B�B�BBB�BkB�BNBB�B�BgB
B,UB0�B1zB75B8^B:8B;XB<wB<�B>�B>�B?�B?�BA�BC�BE�BHBH~BE�BE�BGQBI�BL$B[�Bb.Bd�Bd[Bd/Br8Bt�BtBr�Bq�Bq�Bq�BrBu�Bv�Bw�Bw�Bx�By�B|�B}�B�B�6B�7B�B�B�B�/B�B�[B��B�6B�2B�4B�}B��B�wB�}B�#B�IB�1B��B�oB��B�B��B�ZB�BB�_B��B��B�SB��B��B�lB��B��B�B��B��B�KB�B�.B�.B�jB�2B�aB��B�.B��B�B}kByYBx/BxRBs�Bq(Bo,BrBt�BvIBt�BrLBp�Bo:Bm�Bl�Bl!Bi�BgB^BZkBZlBZyBZ\B\B^{B[�BWNBV�BU�BRjBR	BRBS�BS�BT BS6BSoBRrBQ�BP�BMaBK3BO:BR�BS?BStBS�BQ�BOBN�BM�BHBBF�BE�BExB>;B?�BA�B=�B:vB8B5�B58B6B9yB:�B9�B9�B8sB<�BE�BA�B:�B3UB3>B7�B?ZB@SBA�BA�BB�BC�BC�BD�BDLBF�BH5BGUBF�BILBIfBHpBH[BG�BF5BF)BEBC�BC�BC�BC�BC�BD�BD�BDBC�BD8BF�BHbBG�BC�BE�BFrBFBA�BB�BC�BC�BC B?)B=B<�B9�B5�B4�B4�B1B-�B0^B2&B.�B-�B,�B+B*B)B)B)�B(DB(5B'�B&RB&�B&jB%�B$�B"�B!:B B;BMB&BBBJB�B�BB�B�BZB�B�BPB:BMBHBB~BMB�B�B�B'BBqB�B�B�B
ZB
cB
LB
LB
rB
�B	�B	�B	dB	�B	�B	wB	�B�B�BtB�B\BLB_BvB�BoB�BUB�B�B�BMB6BpBWBVBpB\B8B/B�B�B1B;B@B*BNB\B�B<B7BB PB qB��B�cB�SB�jB�CB�BB�GB�aB��B��B�;B�/B�"B�B�,B�B�B�{B�gB�:B��B��B��B�B��B��B��B��B��B��B��B�B��B��B�"B�B��B�_B�ZB��B�B�B��B��B�B��B��B�B��B�[B��B��B�pB�B��B��B�B��B�B��B��B��B�B��B��B�B�B��B�B�B�B�B�yB��B��B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#�o<#�o<#ا<#�r<#�<#�<#�E<#�U<#��<#��<#�<#ۮ<$�<#ڑ<#�]<#��<#�<#�I<#�*<#ף<#�
<#؄<#׎<#ۮ<#�<#�<#�8<#�<$"2<#�]<#�<#�<#�
<#��<#�X<#�
<#�<#�<#�i<#�<#�<#�i<#�<<#�<#ף<#�"<5��<��x<Fm<%gB<'�B<*{�<.�	<$<<$o�<&R`<'�-<)��<+��<$�<$5w<#�8<#ܯ<$��<$�L<&�?<(�(<$�R<#��<#�C<$Z<&e<$�<#�&<$r<$��<#�<$t <%�l<$/<%Z2<$XX<#��<#��<#�I<#��<#��<$�<#�{<#�&<#�4<)��<$?[<#��<$8�<#ا<#�{<#�r<#�{<#׎<#�m<#��<#׎<#��<#�]<#�<#�]<#�<#��<$Z�<#ܯ<#�<<#�m<$^�<)w�<$�<#�*<#��<#�<'��<$�<#�<#��<#�4<#�<#ޫ<#��<$%<#�<#�<#׎<#�i<#��<#�<#��<#�M<#��<#��<#�8<#׎<#�<#�<#؄<#�<#�U<$
�<#ٛ<#��<#�D<#��<$ K<#��<$<<#��<#�M<#�^<#�<$(<$�<#��<#�<#؄<#�<<#׺<#�!<$�<#�I<#�&<#�	<#ڑ<#�&<$<<$O�<$3U<$�<#�r<#ף<#�I<$��<#�<#�*<#�<$Z<#ڑ<$&<$<<%��<$.<#�(<%��<$�w<$.<#��<$/%<#׺<$�<$�7<$'<$�X<$'<#�<#�a<$!><%2?<)�<$r�<#�<#�<#��<#�&<#�$<%&<%��<#�&<$�<%��<$v<#ٛ<#�D<$�<$8�<#��<#�&<$�<$�<$P�<&,f<$�<$�<$*<#؄<#�!<$
�<$=<$^�<#��<$A�<&��<$�<%G<$��<&�%<#�<#��<$�-<%@�<$��<$^�<#�Q<#�0<#ۮ<#�]<$	�<#�<$p<$v�<&1�<$��<'N(<(B�<#ڑ<$��<%K:<#��<#�<#�0<#�<<#�o<#ا<#�<#��<#�<#ٛ<$�<$2G<#��<$}<$7�<$I�<$2G<$q@<$�<$�<#��<#�<#�<#�<#׺<#�&<#�^<#�&<#��<#�<#��<#�<$9�<$ѩ<#؄<#�{<$Z�<%{@<#�o<#�r<#�4<#��<%�<$7�<$-<%�<%,#<#�<$#(<&�<$�Q<#�<#�e<$�<$n�<#�a<$n�<$�<#��<#׺<#؄<$W<#�<#�<$U�<#��<$'<$�<$<<$k�<$�w<$ �<#�<$a<$�<$ <#�5<#��<$)
<#�<#�<$�<$�J<%Z2<$f�<#��<#�M<$i&<$MO<$e.<$f�<$6�<#ڑ<#�c<#�)<#�<$ <$T�<$A�<#�]<#�J<$�<#�<#ٛ<#�r<#׺<#׺<#ߜ<#��<$<<#�"<#�8<#��<#�<#�<#�<$f<$a<#�<$�<#ܯ<#�D<#�J<#�4<#�)<#�&<$�<#�<#�N<#�a<#�N<#��<#��<#�<#��<#�E<#�<#�<#��<#؄<$
�<$�<#��<#�+<#�J<#��<#�U<#�<$�<#��<#�^<#׺<#�<#��<$k<#��<#�5<$�<#�<#�<#�5<$p<#��<#�c<#�<#�<#��<#�J<#�4<#ޫ<#�N<$%<$
<#�)<#ܯ<#��<#�$<#�<$Gd<#ܯ<#�J<#�]<#�<#��<#׺<#��<#�J<#�l<#��<#�<#�c<$.<$)
<#ޫ<#׎<#�<#�o<#�<#��<#��<#��<$t <#�<$C�<#��<#��<$e.<$<<#��<#׺<#׎<#�<$/<$�<#�<#�U<#��<#�Q<#�M<#��<#�E<#��<#�<#�*<#ޫ<#�D<#�&<$�<#��<#��PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - surface_pressure                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT                                                                                                                                                                                   surface_pressure=-0.210                                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment                                                                                                                                                                            No significant pressure drift detectedn most recent valid surface pressure                                                                                                                                                                                      No significant temperature drift detected                                                                                                                                                                                                                       No significant drift detected in conductivity                                                                                                                                                                                                                   20121107090145              20130503000000  AO  ARGQ                                                                        20121107090145  QCP$                G�O�G�O�G�O�FFBFE           AO  ARGQ                                                                        20121107090145  QCF$                G�O�G�O�G�O�0               AO  ARCAADJP                                                                    20121107090145    IP                G�O�G�O�G�O�                WHOIARSQWHQCV0.3                                                                20130426160715  QC                  G�O�G�O�G�O�                WHOIARSQ CTMV1.0                                                                20130502165016  IP                  G�O�G�O�G�O�                WHOIARSQOW  V1.1CORIOLIS CTD_for_DMQC_2012V1                                    20130503000000  IP                  G�O�G�O�G�O�                