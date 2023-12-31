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
_FillValue                 �  :p   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  <h   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  D<   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  F4   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  N   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  U�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  W�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _�   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  a�   PSAL         
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  it   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  qH   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  s@   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  {   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  }   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �   CALIBRATION_DATE            	             
_FillValue                  ,  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �<   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �@   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �D   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �H   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �L   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��Argo profile    2.2 1.2 19500101000000  4901407 Irminger Sea, Argo equivalent                                   Amy Bower                                                       PRES            TEMP            PSAL              A   AO  20121027081556  20130503154401  4476_5267_259                   2C  D   APEXIR_SBE_3479                                                 846 @�g�qf  1   @�g�u�  @Oݲ-V�;�A�7K�1   GPS     A   A   A   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�33B   B  B  B  B   B(  B0ffB8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Doy�Dp  Dp� Dq  Dq� Dr  Dr� Ds  Dsy�Ds��Dt� DufDu� Dv  Dv� Dw  Dw� Dx  Dx� Dy  Dy� DzfDz� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @L��@���@ə�A��A$��AD��Ad��A�ffA�ffA�ffA�ffA�ffA�ffA�ffA�B33B	33B33B33B!33B)33B1��B933BA33BI33BQ33BY33Ba33Bi33Bq33By33B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���Bę�Bș�B̙�BЙ�Bԙ�Bؙ�Bܙ�B���B䙚B虚B왚B�B���B���B���C L�CL�CL�CL�CL�C
L�CL�CL�CL�CL�CL�CL�CL�CL�CL�CL�C L�C"L�C$L�C&L�C(L�C*L�C,L�C.L�C0L�C2L�C4L�C6L�C8L�C:L�C<L�C>L�C@L�CBL�CDL�CFL�CHL�CJL�CLL�CNL�CPL�CRL�CTL�CVL�CXL�CZL�C\L�C^L�C`L�CbL�CdL�CfL�ChL�CjL�ClL�CnL�CpL�CrL�CtL�CvL�CxL�CzL�C|L�C~L�C�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC��C�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fC�&fD 3D �3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D	3D	�3D
3D
�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D3D�3D 3D �3D!3D!�3D"3D"�3D#3D#�3D$3D$�3D%3D%�3D&3D&�3D'3D'�3D(3D(�3D)3D)�3D*3D*�3D+3D+�3D,3D,�3D-3D-�3D.3D.�3D/3D/�3D03D0�3D13D1�3D23D2�3D33D3�3D43D4�3D53D5�3D63D6�3D73D7�3D83D8�3D93D9�3D:3D:�3D;3D;�3D<3D<�3D=3D=�3D>3D>�3D?3D?�3D@3D@�3DA3DA�3DB3DB�3DC3DC�3DD3DD�3DE3DE�3DF3DF�3DG3DG�3DH3DH�3DI3DI�3DJ3DJ�3DK3DK�3DL3DL�3DM3DM�3DN3DN�3DO3DO�3DP3DP�3DQ3DQ�3DR3DR�3DS3DS�3DT3DT�3DU3DU�3DV3DV�3DW3DW�3DX3DX�3DY3DY�3DZ3DZ�3D[3D[�3D\3D\�3D]3D]�3D^3D^�3D_3D_�3D`3D`�3Da3Da�3Db3Db�3Dc3Dc�3Dd3Dd�3De3De�3Df3Df�3Dg3Dg�3Dh3Dh�3Di3Di�3Dj3Dj�3Dk3Dk�3Dl3Dl�3Dm3Dm�3Dn3Dn�3Do3Do��Dp3Dp�3Dq3Dq�3Dr3Dr�3Ds3Ds��Dt�Dt�3Du�Du�3Dv3Dv�3Dw3Dw�3Dx3Dx�3Dy3Dy�3Dz�Dz�3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A
�A
 �A
$�A
$�A
$�A
$�A
(�A
(�A
$�A
$�A
$�A
 �A
 �A
 �A
 �A
 �A
$�A
(�A
1'A
1'A
-A
(�A
(�A
1'A
A�A
I�A
I�A
I�A
VA
bNA
n�A
n�A
~�A
�A
v�A
r�A
r�A
v�A
�DA
v�A
��A
��A
�A
{A�
A�;@�@�z�@�(�@��
@�C�@���@�1'@�R@�5?@���@��@��D@�5?@�@陚@��@�^5@�o@@�$�@��#@���@�K�@��y@�h@睲@���@�F@�@��@�R@�+@���@�n�@�J@�5?@��@�7L@ᙚ@�+@���@��
@�7L@�(�@�!@�5?@���@�%@� �@�C�@�z�@�^@�ȴ@�@�$�@�M�@�@�S�@���@�9@��`@��@�{@�=q@�^5@�E�@�-@�V@�9@��/@�@�C�@�|�@�|�@畁@�+@�o@�@��#@��@�^@�@���@�?}@�@��@�n�@�=q@�ȴ@�(�@�  @���@�x�@�u@�u@��`@��@�$�@��@�z�@��`@���@�9@��@�&�@��@�@�1'@�9X@��;@�@�P@�\)@�l�@�|�@�t�@�t�@�S�@�o@���@��H@��@��@��H@�M�@�5?@�@�@��@��@��T@�-@�hs@�?}@�@�Q�@� �@㝲@��@�@�K�@�5?@�-@�/@�@��@�I�@�  @��y@�=q@�{@�@ݲ-@�hs@��/@�A�@�t�@�;d@��@ڟ�@�J@���@�G�@؛�@��
@׶F@ם�@�l�@�K�@�o@��@��@�Q�@�1@��@�ƨ@Ӯ@ӥ�@�;d@���@���@�v�@��@љ�@щ7@�`B@�%@��/@���@д9@�A�@� �@���@�S�@���@Ο�@�^5@��@͉7@�G�@�7L@�7L@�7L@��@�Ĝ@˥�@�|�@�C�@��@��@��@ʗ�@�@�@ɉ7@�p�@�p�@�p�@�hs@�V@��/@�bN@�  @�S�@��@Ə\@�M�@���@�X@�X@�O�@�?}@���@�J@őh@�$�@��T@�X@�%@ě�@�1@ÍP@�|�@�dZ@�K�@�33@�@���@§�@�^5@�-@��@��h@�Ĝ@��@�9X@��@��@�+@���@�5?@�@��^@�hs@�hs@���@�^5@�ff@�{@��#@���@��@�X@�/@�%@���@���@��9@���@�Z@�(�@��m@���@��F@��@���@���@��+@��@��T@���@�hs@�X@�`B@�?}@�V@��D@�9X@���@��^@�V@��u@�j@�S�@�ȴ@�n�@�5?@�@���@�`B@�/@���@���@���@�r�@�I�@�(�@�S�@���@�E�@��T@�@�hs@��/@���@�I�@�1'@�  @��;@��@�;d@��@�~�@��T@�G�@���@��u@��@�Z@���@�+@���@�n�@�-@��@���@���@�X@�hs@�hs@�`B@�G�@��@�V@���@��`@�Ĝ@���@��u@��@�j@�9X@�(�@���@�;d@���@��R@���@�-@��@��^@�x�@��@��u@�r�@�bN@�bN@�1@�|�@�+@��@���@��h@�`B@�?}@�7L@�/@�/@�&�@��@��@��@��@��@��@��`@��j@�bN@�ƨ@��@�l�@�+@��@�@�ȴ@���@�v�@�V@�=q@�$�@�@��@��T@��#@��^@��-@��@�G�@�G�@�O�@�G�@�G�@�O�@�X@�X@�G�@���@��@�  @��
@��@�\)@�;d@�@���@�^5@��@��^@��@�V@���@�1'@��@�t�@�dZ@�;d@�
=@���@��R@�v�@�$�@���@��-@���@���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   A
�A
 �A
$�A
$�A
$�A
$�A
(�A
(�A
$�A
$�A
$�A
 �A
 �A
 �A
 �A
 �A
$�A
(�A
1'A
1'A
-A
(�A
(�A
1'A
A�A
I�A
I�A
I�A
VA
bNA
n�A
n�A
~�A
�A
v�A
r�A
r�A
v�A
�DA
v�A
��A
��A
�A
{A�
A�;@�@�z�@�(�@��
@�C�@���@�1'@�R@�5?@���@��@��D@�5?@�@陚@��@�^5@�o@@�$�@��#@���@�K�@��y@�h@睲@���@�F@�@��@�R@�+@���@�n�@�J@�5?@��@�7L@ᙚ@�+@���@��
@�7L@�(�@�!@�5?@���@�%@� �@�C�@�z�@�^@�ȴ@�@�$�@�M�@�@�S�@���@�9@��`@��@�{@�=q@�^5@�E�@�-@�V@�9@��/@�@�C�@�|�@�|�@畁@�+@�o@�@��#@��@�^@�@���@�?}@�@��@�n�@�=q@�ȴ@�(�@�  @���@�x�@�u@�u@��`@��@�$�@��@�z�@��`@���@�9@��@�&�@��@�@�1'@�9X@��;@�@�P@�\)@�l�@�|�@�t�@�t�@�S�@�o@���@��H@��@��@��H@�M�@�5?@�@�@��@��@��T@�-@�hs@�?}@�@�Q�@� �@㝲@��@�@�K�@�5?@�-@�/@�@��@�I�@�  @��y@�=q@�{@�@ݲ-@�hs@��/@�A�@�t�@�;d@��@ڟ�@�J@���@�G�@؛�@��
@׶F@ם�@�l�@�K�@�o@��@��@�Q�@�1@��@�ƨ@Ӯ@ӥ�@�;d@���@���@�v�@��@љ�@щ7@�`B@�%@��/@���@д9@�A�@� �@���@�S�@���@Ο�@�^5@��@͉7@�G�@�7L@�7L@�7L@��@�Ĝ@˥�@�|�@�C�@��@��@��@ʗ�@�@�@ɉ7@�p�@�p�@�p�@�hs@�V@��/@�bN@�  @�S�@��@Ə\@�M�@���@�X@�X@�O�@�?}@���@�J@őh@�$�@��T@�X@�%@ě�@�1@ÍP@�|�@�dZ@�K�@�33@�@���@§�@�^5@�-@��@��h@�Ĝ@��@�9X@��@��@�+@���@�5?@�@��^@�hs@�hs@���@�^5@�ff@�{@��#@���@��@�X@�/@�%@���@���@��9@���@�Z@�(�@��m@���@��F@��@���@���@��+@��@��T@���@�hs@�X@�`B@�?}@�V@��D@�9X@���@��^@�V@��u@�j@�S�@�ȴ@�n�@�5?@�@���@�`B@�/@���@���@���@�r�@�I�@�(�@�S�@���@�E�@��T@�@�hs@��/@���@�I�@�1'@�  @��;@��@�;d@��@�~�@��T@�G�@���@��u@��@�Z@���@�+@���@�n�@�-@��@���@���@�X@�hs@�hs@�`B@�G�@��@�V@���@��`@�Ĝ@���@��u@��@�j@�9X@�(�@���@�;d@���@��R@���@�-@��@��^@�x�@��@��u@�r�@�bN@�bN@�1@�|�@�+@��@���@��h@�`B@�?}@�7L@�/@�/@�&�@��@��@��@��@��@��@��`@��j@�bN@�ƨ@��@�l�@�+@��@�@�ȴ@���@�v�@�V@�=q@�$�@�@��@��T@��#@��^@��-@��@�G�@�G�@�O�@�G�@�G�@�O�@�X@�X@�G�@���@��@�  @��
@��@�\)@�;d@�@���@�^5@��@��^@��@�V@���@�1'@��@�t�@�dZ@�;d@�
=@���@��R@�v�@�$�@���@��-@���@���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�jB�jB�qB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�jB�dB�jB�jB�qB�qB�qB�qB�wB�}B��B��B��B��B��B��B��B��BB��BÖBɺB��B�
B�;B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�`B�HB�fB  B+B1B
=B	7BB%B1BBB��B��B��BB��B��B��B��B��BBB  BBDBuB{BDB
=B�B�B�B{BhBoB�B$�B,B.B,B-B2-B49B9XB?}BA�BF�BJ�BK�BM�BL�BJ�BG�BF�BH�BO�BW
BYBZB[#BZB[#Be`Bl�Bn�Bm�Bn�Bn�BjBcTBaHB_;B_;BbNBjBjBhsBbNB^5B^5BaHBe`BiyBo�By�B|�B}�B}�B� B�B�B� B� B� B� B~�B~�B~�B~�B�B�B�B�B�B� B� B� B� B�B}�B}�B}�B~�B�B�B�B� B~�B~�B|�Bz�By�Bx�Bw�Bz�By�Bu�Bs�Bq�Bo�Bn�Bm�Bl�BjBiyBhsBgmBgmBffBffBdZBaHB`BB_;B^5B\)BZBYBW
BT�BS�BS�BR�BQ�BP�BO�BI�BH�BG�BG�BI�BI�BI�BI�BH�BG�BG�BG�BG�BH�BH�BG�BG�BF�BF�BE�BD�BC�BB�BA�B@�B?}B>wB>wB>wB=qB=qB=qB=qB<jB:^B:^B9XB9XB9XB9XB8RB8RB8RB7LB7LB7LB7LB7LB6FB5?B49B33B1'B0!B/B.B-B,B,B,B-B1'B2-B1'B5?B5?B33B2-B0!B.B-B-B-B-B-B-B-B,B+B)�B(�B&�B$�B$�B$�B#�B"�B!�B �B�B �B �B �B"�B%�B-B.B.B/B.B.B.B.B-B-B/B0!B0!B0!B/B/B/B0!B/B/B.B.B.B-B-B-B-B-B-B-B,B)�B&�B$�B#�B"�B!�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B{BuB{BuBoBbB\B\BVBVBPBDB
=B	7B	7B1B+B+B+B+B+B+B+B%B%B%B%B%BBBBBBBBBB  B  B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�yB�yB�yB�yB�yB�yB�yB�yB�sB�sB�mB�mB�mB�m111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   B�TB�_B�pB�sB�qB�gB�sB�}B�sB�qB�vB�jB�gB�hB�gB�]B�]B�RB�iB�uB�sB�eB�NB�?B�XB�qB�nB�QB�QB�ZB�B�TB�}B��B��B��B�vB�RBµB�4B��BɎBФB�2B��B�yB��B�cB�eB��B� BZB�&B��B�VB��B "BqB�cB�B��B߶B��B�B�B
�B�BCB�B
9B-B�B��B��B�uB�B�MB��B�jB��B��B�B�B�uB�B	lByB B�B�B$BEB�B�B�B�B�B#KB+�B/JB+�B,B1�B3@B8BB?B@BE�BJ�BK�BM�BM�BK�BH4BFfBGBNBV�BYBZB[�BZ(BX�Bc�BlBoBm�Bn�BoqBl�BdFBb
B_B^lB`eBj�BlBj�Bc�B^>B]�B`.Bd�Bh:Bm_By1B}	B~B}�B�B�\B��B�}B�B��B�LB2BAB~�B~�B�B�B�>B�mB�-B�*B�B�B�B��B~*B~;B}�B~�B�=B�$B�RB�qBHB�B}xB{9Bz�By�Bw=B{=B{jBv�Bt�BroBo�Bn�BnBn4Bk�Bi�Bh�Bg�Bg�BgABgRBe�Ba�B`�B_�B_B\�BZ�BZBX+BU2BT!BT>BS)BRHBQ^BR�BJ�BI%BG�BG�BI�BI�BJZBJBH�BH{BH�BH-BG�BH�BI9BG�BG�BF�BGRBE�BD�BD�BCVBA�B@�B@!B?B>�B>�B=wB=tB=�B=�B>
B:�B:�B9�B9^B9�B9�B94B8�B8�B7mB7NB7OB7`B7�B6�B5�B4�B49B1�B0�B/�B.�B-�B,B,B,B,CB0�B2�B0aB5�B6B3�B2�B1B.�B-1B-5B-5B-9B-ZB-[B-PB,yB+UB*XB)�B(B%SB%QB%QB$FB#�B"�B!yB rB �B!=B �B"B$�B,�B.�B.rB/cB.YB.WB.WB.VB-hB-PB/B08B0�B0oB/�B/CB/EB0pB/�B/�B.SB.�B.(B-kB-pB-(B-B-AB-bB-�B,�B,uB(JB%�B$�B#'B#lB �B MB B B UB $BB B B B�B�B�B�B�BaB5B�B$B]B�BB�B�B�B�B,BB
BdBVB�B�BwB�B�B�B
�B	�B	�B�BcBxB�BB)B7BOBfBCB?BMBVBXB;B:BDBdB<B�B�B�B *B  B��B�WB�RB�YB��B��B�B��B��B�nB��B�YB�3B�bB�[B�B��B��B��B��B��B��B��B��B��B��B��B�B�B�KB��B��B�B�B��B��B�B��B��B��B��B��B��B��B�B�B��B�B��B�B�B�B�B�B�B�B�B��B�&B�mB�oB��B��B�B��B��B�2B��B��B�)B�tB�B�B�B�uB�B�B�B��B��B�B��B��B��B�B�|B�B�p111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   <#؄<#�i<#�<#�<#�
<#�X<#�<#�{<#�<#�
<#�{<#�
<#�<#�<#�<#׎<#׎<#��<#�<#�i<#�I<#�<#�o<#ܯ<#��<#�
<#�<#�*<#�r<#��<#�<#��<#�{<#�*<#��<#�<#׎<#�^<#�r<#�<$2G<#��<)�<G�l<���<X+<&�R<$
�<$v<$��<0�=<Aә<'��<$f�<$�<#�<'��<-�<4D�<'��<#��<C�P<%U�<$�<$8�<$f<)��<$�q<$8�<'F<*�</(<&�z<$��<$�<-<$�<#�M<#�N<$�<#��<$L<$'<$�<%t�<&y�<#�<-�`<%��<-�M<$<$*<$��<%Q�<$�t<&y<&�8<%�Z<#�4<$�L<#ޫ<$�e<$�<$��<$�;<#�Q<$�k<$4e<#�<#ܯ<#��<$g�<$�V<$�<#�U<%�!<&|V<#�<#�I<#�<$'<#�i<(F.<&Z�<$ <$<<#�0<#�<$f�<(��<$�b<$I�<#�&<$Y�<&��<#�<%�<'��<%K:<#�I<$�<$�<$_�<%s<'�-<$/%<#�D<#ٛ<#�!<#�<#�<$&<$�<#�<$r<#�<#��<#�l<#�D<#�C<#�0<#�&<#ޫ<#�<#ۮ<#�l<#ף<#�<#�<$I�<#��<#�l<#�<#�o<#�E<#��<#�<#��<#�<$b�<$<#�<$R'<$A�<$
<#��<%�<$b�<$U�<$MO<#�l<#�<$<<%��<$�<#�<#�*<$
�<$<<$i&<$��<$�<#��<$ <$�<$f�<#�H<$Y�<$�j<$�(<#�J<#�+<#��<#�E<#��<$�<)��<$�J<#��<#�+<#��<#�o<#�$<$%<#�<#��<$W<$b�<$<<#�8<#�&<$�<#��<#׺<#�N<$/%<#�J<#��<$�Q<$O�<#�<#�Q<$)
<$)
<#��<#�o<#�&<#�<#ߜ<$�<%��<#�!<#��<$�<#�&<#�+<#��<$r�<#�Q<#�<#�]<#�<#�<#�C<$r<#�<$@|<$&<$�<$E<#�a<#��<$H�<$><<#׎<#׺<#��<$T�<#�g<$}<$N�<#�)<$U�<$�<$/%<$t <$C�<#��<#ۮ<#ۮ<#ܯ<#�<#�!<#�U<#��<#�<#��<$"2<$�!<$�<$ <$ <#��<$T�<$T�<$9�<$9�<#�c<$�<#�<#�<$��<#�<$�<#��<#��<#�<#�<#�<#�U<#��<#�U<#��<#ا<#��<#�<#�g<#��<#�l<#�<$]h<$<<#�&<$x+<#�C<#�m<#�W<#�<#�<<#��<#�<$J�<$!><(��<%Q�<$��<$Gd<#�<%�d<$t <$f<#��<#�<$�<#��<#�<#�l<#��<#�<#�U<#�<#�<%b<$z�<$H�<$ K<#�e<$'<$Z�<#��<$<#�8<#�4<#��<#�<$=<$
<$�<$��<$�V<$�<$�<#�]<#�<$#(<%s<$-<$.<#��<#�	<#��<#�!<#�<#��<#�<#�{<#��<#��<#��<#�<#��<#�^<#��<#�o<#�D<#ܯ<#�4<#��<$c�<$#(<$)
<#�l<#�*<$A�<#�m<#��<#�)<$!><$R'<#�E<#�o<#�&<$�<$aD<$.<#��<$�<%�Z<#��<#�J<#׺<#ף<#�<#ף<#ا<#�{<#�<#�
<#�
<#׺<#�<#��<$f<$��<#�8<#�g<#�)<#�o<#�+<#��<#�<#�&<#ߜ<#��<#�+<#ޫ<#��<#��<#ף<#ޫ<#�<#�<#�<#�<#׺<#�{<#�
<#ף<#׎<#�<#ڑ<$�<$J�<$T�<#�<#��<$	<#��<#��<$Z<#�g<#��<$#(<$� <#��<#ߜ<$��<$�e<#��<#�D<#��<#�<#��<#�8<#��<$�<$�<#��<#׺<#�C<#�PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - surface_pressure                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT                                                                                                                                                                                   surface_pressure=-0.330                                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment                                                                                                                                                                            No significant pressure drift detectedn most recent valid surface pressure                                                                                                                                                                                      No significant temperature drift detected                                                                                                                                                                                                                       No significant drift detected in conductivity                                                                                                                                                                                                                   20121027081556              20130503000000  AO  ARGQ                                                                        20121027081556  QCP$                G�O�G�O�G�O�FFBFE           AO  ARGQ                                                                        20121027081556  QCF$                G�O�G�O�G�O�0               AO  ARCAADJP                                                                    20121027081556    IP                G�O�G�O�G�O�                WHOIARSQWHQCV0.3                                                                20130426160715  QC                  G�O�G�O�G�O�                WHOIARSQ CTMV1.0                                                                20130502165016  IP                  G�O�G�O�G�O�                WHOIARSQOW  V1.1CORIOLIS CTD_for_DMQC_2012V1                                    20130503000000  IP                  G�O�G�O�G�O�                