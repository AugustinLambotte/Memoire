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
_FillValue                 �  :,   PRES_ADJUSTED            
      	   	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    comment       $In situ measurement, sea surface = 0   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  <   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  C�   PRES_ADJUSTED_ERROR          
         	long_name         SEA PRESSURE   
_FillValue        G�O�   units         decibar    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  E�   TEMP         
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  M   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   TEMP_ADJUSTED            
      	   	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  V�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ^   TEMP_ADJUSTED_ERROR          
         	long_name         $SEA TEMPERATURE IN SITU ITS-90 SCALE   
_FillValue        G�O�   units         degree_Celsius     comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  _�   PSAL         
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  g�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  o   PSAL_ADJUSTED            
      	   	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    	valid_min                	valid_max         B(     comment       In situ measurement    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  q    PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  x�   PSAL_ADJUSTED_ERROR          
         	long_name         PRACTICAL SALINITY     
_FillValue        G�O�   units         psu    comment       WContains the error on the adjusted values as determined by the delayed mode QC process.    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  zt   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    �4   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �4   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �4   CALIBRATION_DATE            	             
_FillValue                  ,  �4   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    �`   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    �d   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    �h   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �l   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �p   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    ��Argo profile    2.2 1.2 19500101000000  4901407 Irminger Sea, Argo equivalent                                   Amy Bower                                                       PRES            TEMP            PSAL               A   AO  20121011120150  20130503154400  4476_5267_256                   2C  D   APEXIR_SBE_3479                                                 846 @�c�ѻ  1   @�c��� @O�vȴ9�;������1   GPS     A   A   A   @�33@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CI�fCL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� DdfDd� De  De� Df  Df� Dg  Dg� Dh  Dh� DifDi�fDj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr�f1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@�  @���AffA&ffAFffAfffA�33A�33A�33A�33A�33A�33A�33A�33B��B	��B��B��B!��B)��B1��B9��BA��BI��BQ��BY��Ba��Bi��Bq��By��B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���B���C ffCffCffCffCffC
ffCffCffCffCffCffCffCffCffCffCffC ffC"ffC$ffC&ffC(ffC*ffC,ffC.ffC0ffC2ffC4ffC6ffC8ffC:ffC<ffC>ffC@ffCBffCDffCFffCHffCJL�CLffCNffCPffCRffCTffCVffCXffCZffC\ffC^ffC`ffCbffCdffCfffChffCjffClffCnffCpffCrffCtffCvffCxffCzffC|ffC~ffC�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33C�33D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV�DV��DW�DW��DX�DX��DY�DY��DZ�DZ��D[�D[��D\�D\��D]�D]��D^�D^��D_�D_��D`�D`��Da�Da��Db�Db��Dc�Dc��Dd  Dd��De�De��Df�Df��Dg�Dg��Dh�Dh��Di  Di� Dj�Dj��Dk�Dk��Dl�Dl��Dm�Dm��Dn�Dn��Do�Do��Dp�Dp��Dq�Dq��Dr�Dr� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A&�A&�A+A+A&�A+A+A+A+A+A/A33A7LA33A/A/A33A/A/A33A;dA7LA33A+A"�A�A"�A&�A+A+A�A&�A��A��A�\AffA�7A��A`BA��A	�A;dAC�@�+@�ȴ@��@���@�o@�@�I�@��;@��@���@ى7@�j@׍P@��#@��@ԋD@�Q�@�b@���@�j@���@�p�@�&�@�%@�%@��@ԓu@�1@�Z@��/@�p�@��@�ff@ղ-@Չ7@�5?@ְ!@�5?@�{@���@���@Դ9@ԋD@�9X@�bN@ԣ�@�bN@�I�@��m@�  @��@�@���@��@���@�E�@ѡ�@�p�@��@��@���@���@�Ĝ@У�@�r�@��/@��T@�{@��T@��@�r�@У�@��@҇+@Ӿw@� �@��/@��T@׍P@�z�@���@�dZ@�t�@׍P@�\)@��@��@�z�@�+@�A�@Ϯ@�x�@д9@�j@�|�@�p�@�G�@��@�z�@�\)@���@ʸR@ʇ+@�ff@�M�@�E�@�E�@ʗ�@ʰ!@ʰ!@���@ʏ\@�r�@�G�@���@�z�@ˮ@�9X@��@��@�
=@Η�@���@��@�E�@��H@�S�@Ϯ@�l�@ΰ!@�ff@�$�@�Ĝ@��m@˝�@�dZ@�K�@�dZ@˅@��@�p�@�~�@�@���@��@��@�n�@Ə\@Ə\@Ƈ+@�ff@�-@�$�@�^5@Ƈ+@Ƈ+@Ƈ+@�^5@ŉ7@��H@+@�^5@�$�@��T@�`B@�&�@�V@��j@���@��!@�n�@��@�@��T@��^@�hs@�%@�z�@���@��@���@�M�@�`B@�1'@��;@��@��!@�x�@��@��u@�A�@� �@��@��@�(�@�9X@� �@���@���@�l�@�\)@�+@�"�@�@��H@��@���@���@��\@��+@��+@�v�@�^5@�V@�^5@�V@�$�@��@��@�J@��@��T@���@��@���@���@�Z@��m@���@�l�@�o@��!@�V@�{@��T@�X@�1'@���@�dZ@�K�@�+@�
=@���@��!@��R@��@�z�@�G�@�X@���@�n�@���@��y@�33@�
=@��!@�v�@�^5@�5?@���@�X@�X@�G�@�?}@� �@���@���@�l�@�S�@�33@�+@�"�@�@���@���@�M�@��T@�hs@�?}@�Ĝ@�r�@�A�@�  @��@��@���@�G�@��D@��@�^5@�p�@�X@�G�@�G�@���@�j@��@��@���@�l�@�S�@�C�@�33@�"�@��y@���@�5?@���@��m@���@�E�@��@���@���@��@�Q�@�  @� �@�  @�t�@�
=@���@��@���@��@���@���@�r�@�bN@�I�@��@��@���@�n�@�5?@�@��#@���@��-@��7@�X@��j@��u@��@�z�@�z�@�z�@�Q�@�1'@��@�ƨ@�K�@�o@���@�V@�5?@�$�@��T@�@��^@��^@��-@���@��-@��h@�x�@�hs@�`B@�G�@�?}@��@���@��9@�j@��m@��F@��@��@��@���@���@���@��@�+@�o@��@��@���@�^5@�=q@�5?@�{@��@��T@���@��-@���@��h@���@��7@��@�7L@��/@�Ĝ@���@���@��9@�Z@�1@��;@��;@�ƨ@���@���@���@���@��
@��;@��;@���@��w@��w@��w@��F@���@���@�dZ@��H@���@���@���@�ȴ@�ȴ@�ȴ@���@���@��R@��R@��R@��R@��R1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111A&�A&�A+A+A&�A+A+A+A+A+A/A33A7LA33A/A/A33A/A/A33A;dA7LA33A+A"�A�A"�A&�A+A+A�A&�A��A��A�\AffA�7A��A`BA��A	�A;dAC�@�+@�ȴ@��@���@�o@�@�I�@��;@��@���@ى7@�j@׍P@��#@��@ԋD@�Q�@�b@���@�j@���@�p�@�&�@�%@�%@��@ԓu@�1@�Z@��/@�p�@��@�ff@ղ-@Չ7@�5?@ְ!@�5?@�{@���@���@Դ9@ԋD@�9X@�bN@ԣ�@�bN@�I�@��m@�  @��@�@���@��@���@�E�@ѡ�@�p�@��@��@���@���@�Ĝ@У�@�r�@��/@��T@�{@��T@��@�r�@У�@��@҇+@Ӿw@� �@��/@��T@׍P@�z�@���@�dZ@�t�@׍P@�\)@��@��@�z�@�+@�A�@Ϯ@�x�@д9@�j@�|�@�p�@�G�@��@�z�@�\)@���@ʸR@ʇ+@�ff@�M�@�E�@�E�@ʗ�@ʰ!@ʰ!@���@ʏ\@�r�@�G�@���@�z�@ˮ@�9X@��@��@�
=@Η�@���@��@�E�@��H@�S�@Ϯ@�l�@ΰ!@�ff@�$�@�Ĝ@��m@˝�@�dZ@�K�@�dZ@˅@��@�p�@�~�@�@���@��@��@�n�@Ə\@Ə\@Ƈ+@�ff@�-@�$�@�^5@Ƈ+@Ƈ+@Ƈ+@�^5@ŉ7@��H@+@�^5@�$�@��T@�`B@�&�@�V@��j@���@��!@�n�@��@�@��T@��^@�hs@�%@�z�@���@��@���@�M�@�`B@�1'@��;@��@��!@�x�@��@��u@�A�@� �@��@��@�(�@�9X@� �@���@���@�l�@�\)@�+@�"�@�@��H@��@���@���@��\@��+@��+@�v�@�^5@�V@�^5@�V@�$�@��@��@�J@��@��T@���@��@���@���@�Z@��m@���@�l�@�o@��!@�V@�{@��T@�X@�1'@���@�dZ@�K�@�+@�
=@���@��!@��R@��@�z�@�G�@�X@���@�n�@���@��y@�33@�
=@��!@�v�@�^5@�5?@���@�X@�X@�G�@�?}@� �@���@���@�l�@�S�@�33@�+@�"�@�@���@���@�M�@��T@�hs@�?}@�Ĝ@�r�@�A�@�  @��@��@���@�G�@��D@��@�^5@�p�@�X@�G�@�G�@���@�j@��@��@���@�l�@�S�@�C�@�33@�"�@��y@���@�5?@���@��m@���@�E�@��@���@���@��@�Q�@�  @� �@�  @�t�@�
=@���@��@���@��@���@���@�r�@�bN@�I�@��@��@���@�n�@�5?@�@��#@���@��-@��7@�X@��j@��u@��@�z�@�z�@�z�@�Q�@�1'@��@�ƨ@�K�@�o@���@�V@�5?@�$�@��T@�@��^@��^@��-@���@��-@��h@�x�@�hs@�`B@�G�@�?}@��@���@��9@�j@��m@��F@��@��@��@���@���@���@��@�+@�o@��@��@���@�^5@�=q@�5?@�{@��@��T@���@��-@���@��h@���@��7@��@�7L@��/@�Ĝ@���@���@��9@�Z@�1@��;@��;@�ƨ@���@���@���@���@��
@��;@��;@���@��w@��w@��w@��F@���@���@�dZ@��H@���@���@���@�ȴ@�ȴ@�ȴ@���@���@��R@��R@��R@��R@��R1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B�{B��B��B��B��B��B��B��B��B��B��B��B�B�!B�'B�?B�dB�qB�}BBBĜBƨBŢBÖBÖBǮB��B��B��B��B�
B�#B�5B�HB�NB�ZB�fB�yB�yB�B�B��B��B��B��B  B%B
=B1B1B+B+B1B	7B
=BVBoBhBuBuBoBhBoBoBoBuBuBoBhBhBbBhBhBhBoBuB�B�B�B�B�B�B�B!�B+B33B8RB=qBD�BN�BT�BR�BQ�BR�BR�BR�BP�BK�BE�B>wB1'B.B7LB49B49B/B&�B%�B$�B"�B �B�B�B�B�B�B �B!�B#�B%�B&�B)�B)�B5?B9XB8RB7LB6FB<jBG�BL�BM�BK�BI�BK�BO�BT�BXBZBYBW
BVBT�BO�BL�BK�BK�BK�BK�BL�BI�BB�B9XB8RB8RB:^B;dB=qB>wB?}B>wB>wB>wB@�BB�BD�BE�BE�BD�BA�B9XB7LB8RB8RB7LB5?B5?B49B33B0!B.B-B,B,B+B+B+B)�B(�B%�B$�B#�B!�B�B�B�B�B{BhBbB\B\B\B\B\B\BbBbBbB\BbBbBhBbBhBhBhBhBhBhBhBhBoBoBoBoBuBoBuBuBuB{B{B�B�B�B�B�B{BuBoBhBbB\BVBPBDB	7B	7B	7B	7B	7B	7B
=B
=BDBPB�B�B�B �B%�B'�B(�B)�B-B-B-B,B,B,B+B+B)�B(�B&�B&�B%�B%�B%�B%�B%�B%�B$�B$�B$�B#�B#�B#�B"�B"�B!�B!�B �B�B�B�B�BoBPB	7BBBBBBBBBBBBBBBBB��B��B��B�B�B�B�B�B�yB�sB�mB�mB�fB�ZB�TB�NB�HB�BB�;B�/B�/B�)B�#B�#B�B�B�
B�
B�
B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BɺBɺBɺBɺBɺBɺBɺBɺBȴBȴBȴBȴBǮBǮBǮBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBŢBŢBŢBŢBŢBŢBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBĜBŢBŢBŢBŢBŢBŢBŢBŢBŢBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨBƨ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B�sB�iB�uB��B�jB�vB�vB�yB�vB�mB�iB�mB��B��B�wB�mB��B�uB�mB�cB��B��B��B��B��B�nB�pB�pB�|B��B�nB�B�B�LB�:B�{B�1B�UB�B��B�$B��B��B��B�PB��B�OB��B¾B�B�kB�B�B�cB�	B�"B�B�EB�B�4B��B�=B�(B�}BޘB�wB�PB�B��B�2B�B�B��B�B�4B��B�(B�BoB
�BlBkB�B�BwB	�B
 B�B�B�B BaB�B�B{B�B�B-BoB�B�B�B�BpBvB�B�B�BBkBB BeBZB�B�B)&B2|B7)B;�BB"BMoBU�BS�BQ�BR�BS>BShBR�BN BG�BB�B1�B+�B82B4�B5�B2#B'LB&,B%�B$�B!�B�BB�B�B�B �B!XB#�B%�B&�B*'B'RB4B9�B9B8iB5jB9�BFSBL�BNyBL�BI�BKABN�BTJBW�BZ}BZ)BW�BV�BWBQ7BMLBLBK�BK�BK�BM�BL	BF�B:PB8�B8*B:B:�B==B>xB?�B>�B>�B>�B@)BBQBD�BE�BE�BFBE`B:B7�B8�B8�B8B5�B5mB4�B4�B1�B.B-�B,5B,8B+HB+�B+�B*�B*5B&�B%�B$oB#:B �B,B�BBUBHB�B�B�BgBZBCBHB�B�B�B�B}B�BpB�B�BxB�B�BtBtBiB�B�B{BcB�B�B}BsB�B�B�B�BVB�B�B�B/B�B�B�BB�B�B�B1BB
B	�B	aB	kB	lB	�B
SB
0B
�BB@B}BGB�B%>B'�B(�B*AB-�B-gB-6B,RB,�B,�B+B+B*!B*�B'�B&�B&$B&	B&B%�B%�B&B%IB%$B%QB$�B$�B$B#�B#OB"B"1B!QB!BB�B�BBB
�BbB>B'B�B�B�BvB0B_B?B4B6B4BoB}B�B �B��B��B�B�*B�B�B�B�	B��B�FB�B�;B�B�B�B��B�B�B݀B��B�IB�WB�B�9B�dB׍B�gB�WB�EB�B�,B�CB�YB��B�:B�B� B��B��B�2B�'B�B�qBѧB�DBЄB�XB�B��B�8B�B��B��B��B��B��B��B��B��B��B��B��B�B��B�8B�;BʇB�B��BɹBɽB��B��BɻB��B�<B��B��B��B��B�DB��BƸB��B��BƸB��B��B��B��BƞBƿBŸB�B�.B��B��BšBĄB�&B�B��BĢB��BđBĜBĞBĜBđBđBġBźBźBţBţBŮBŻB��B��B�jB��BƫBƫBƶBƨBƫBƶBƩBƶBƨBƨBƦBƪBƥ1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�<#�{<#�
<#׎<#�i<#�<#�<#�<#�<#�<<#�{<#�<<#�0<#��<#�<#�<<#�{<#�
<#�<<#��<#�i<#�{<#��<#ا<#�i<#׎<#�i<#�i<#�<#��<#�$<$�<$/<$7�<$�<*A<U!m<a=x<I��<YXb<T"x<z��<o�,<C��<7��<V��<e��<8�<8}�<B�2<*�<-nV<,�?<&)�<%^�<(��<%k�<$�<#�<#�g<#��<$4e<$r�<$+<#��<#��<#�<#��<$�<$?[<$ �<$W<$t <$B�<$2G<$��<#�<$��<$<<$k<#�<#�N<%��<#��<#��<#�a<#�e<#�m<#�<#��<$�<#�C<% <#�E<#�{<#�C<#�<$><<$�X<#�<$�<#�M<#��<#�<<#ף<#�J<#�<$&<%�<#�!<#�<%Q�<$F9<#��<$^�<'G�<&��<$=<$�`<%��<(�H<%e<$.<$g�<#�C<#�]<#�<$�<&/<'��<'hA<1�<$aD<(�<$x+<$�<%�<*�&<#��<#�M<$�<&�<$b�<#�8<#�!<#��<#�r<#��<#�&<#�a<#�8<#�<#ۮ<#ܯ<)K?<%�<#��<$J�<$�1<$j|<)��<%D�<#ޫ<$+<$|d<#��<$�<$|d<$9�<$a<#�(<$��<$p<$v<'�<%>�<$<<#�<#��<#�8<#ٛ<$��<'��<1�{<$�Q<#�!<#��<#�<$�<#�J<#�<#�X<#��<#�<#׺<#��<#��<#�<#�<#�<%b�<.��<$A�<#�l<#�5<#��<$Gd<#��<#�<$�<%��<%it<#��<$p<#�8<#�<#��<$	�<$Z<$n�<%s<$c�<$8�<$}<%p<&W�<$�<#�<%��<&/<$o�<$�<$a<#��<#�i<#�<#��<#�C<#�8<#�U<$�<#��<#�D<#�4<#�<<#��<#�^<#��<#�r<#�e<#�{<#�{<#�<#�o<#��<#�{<#�{<#��<#�<#ף<#�<#�<#�+<#�*<#�a<$aD<#�	<#�N<#�N<$3U<$�<#�<$F<$Z<$�<#�)<#�"<$q@<&<�<$e.<#�	<#�l<#�J<#ߜ<#��<#؄<#׎<#�a<'�|<%4L<#ߜ<#��<%:<$*<#ޫ<#�H<#�<$	<#�5<#��<#�<$/%<$/%<#�<#ا<#�8<%�<$U�<#�<#��<#�r<#ޫ<#ף<#׎<#��<#��<#�l<$ <$/%<$C�<#��<$Gd<$v<#�<#��<$�<%<%D�<$�q<$�t<%��<&U"<%`�<#�<#��<#�<<$p<$|d<$Gd<#��<#��<#�!<#�r<#�D<#ٛ<#�D<#�m<#��<$L<&�<&/<&O�<$\"<$�<$<<$�<$�<$F<$/<#ۮ<#ߜ<$aD<$1:<$=<$Gd<$�<$/<$N�<#�<$'<#�*<#�J<$}�<$��<#�<$a<#�m<#�!<#��<#�D<#��<#�&<#�N<$z�<#�U<#�D<#ף<#�<#�<#�<#ߜ<#�l<$�<$I�<#�<$*<$<<#��<#�]<#�)<#�J<#��<#�<#׎<#׺<#�&<#ޫ<#ܯ<#�D<#��<#��<#��<#��<#��<$<<$p<$W<#�<#׎<#�<#�<#ף<#��<#�<#�]<$r<#�l<#�<#�<#�<$�<#�N<#��<#�J<#�&<#��<#�D<#�J<#�D<#�<#�X<#ا<#؄<$ �<$�<#��<#��<#�<#��<$<$	�<#�&<#�&<#��<#�i<#�
<#�<#�
<#�i<#�i<#�<#��<#��<#�<#�<#�{<#��<#��<#�<$P�<#�D<#�<#�<#ף<#�
<#�<#ף<#�<#ף<#�
<#�
<#�<#�<#�PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - surface_pressure                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            PSAL_ADJ corrects Conductivity Thermal Mass (CTM), Johnson et al., 2007, JAOT                                                                                                                                                                                   surface_pressure=-0.370                                                                                                                                                                                                                                         none                                                                                                                                                                                                                                                            CTM: alpha=0.141C, tau=6.89s, rise rate = 10 cm/s with error equal to the adjustment                                                                                                                                                                            No significant pressure drift detectedn most recent valid surface pressure                                                                                                                                                                                      No significant temperature drift detected                                                                                                                                                                                                                       No significant drift detected in conductivity                                                                                                                                                                                                                   20121011120150              20130503000000  AO  ARGQ                                                                        20121011120150  QCP$                G�O�G�O�G�O�FFBFE           AO  ARGQ                                                                        20121011120150  QCF$                G�O�G�O�G�O�0               AO  ARCAADJP                                                                    20121011120150    IP                G�O�G�O�G�O�                WHOIARSQWHQCV0.3                                                                20130426160714  QC                  G�O�G�O�G�O�                WHOIARSQ CTMV1.0                                                                20130502165016  IP                  G�O�G�O�G�O�                WHOIARSQOW  V1.1CORIOLIS CTD_for_DMQC_2012V1                                    20130503000000  IP                  G�O�G�O�G�O�                