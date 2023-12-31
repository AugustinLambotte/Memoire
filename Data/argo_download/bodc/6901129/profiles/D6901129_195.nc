CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
Argo float     history       06-May-2016 13:52:41Zcreation      
references        (http://www.argodatamgt.org/Documentation   comment       bThis netCDF file is generated using BODC's argoReader and netCDF writer software (argo@bodc.ac.uk)     user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    <H   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    <X   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    <\   REFERENCE_DATE_TIME                	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    <`   DATE_CREATION                  	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    <p   DATE_UPDATE                	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    <�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    <�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  <�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  <�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  =   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        =H   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    =L   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    =P   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     =T   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    =t   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    =x   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     =|   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     =�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     =�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    =�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   axis      T      
_FillValue        A.�~       
resolution        >�E�vQ�        =�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    =�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
_FillValue        A.�~       
resolution        >�E�vQ�        =�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   	valid_min         �V�        	valid_max         @V�        axis      Y      
_FillValue        @�i�            =�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    	valid_min         �f�        	valid_max         @f�        axis      X      
_FillValue        @�i�            =�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    >   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    >   VERTICAL_SAMPLING_SCHEME                   	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    >   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ?   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ?   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ?   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        axis      Z      
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  E�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  Lx   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  S$   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  T�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  V|   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  X(   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ^�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  e�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  l,   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  m�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  o�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  q0   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  w�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ~�   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  �4   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    �   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    �    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    �<   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �X   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  �   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �|   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        �   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        �$   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        �@   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �\Argo profile    3.1 1.2 19500101000000  20210225044450  20210225044450  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125482                          2C  D   APEX                            6229                            120210                          846 @�*����1   @�*����@P��Q��4�fffff1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @���@�  A   A   A@  Aa��A���A���A�ffA�33A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BW��B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C�fC  C  C  C  C   C"  C#�fC&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV�Bv�Bv�Bv�B~�B�DB��B�-B
=B&�B��B��B	aHB
!�B
�oB
ƨB+B<jB]/Be`B`BB]/B[#BVBM�BC�B8RB.B!�B�B	7BB
��B
�B
�B
��B
��B
��B
��B
�B
�B
�B
�B
�B
�B
�B
��BBB	7BJBhB�B"�B&�B/B8RB;dB@�BD�BD�BG�BF�BH�BL�BL�BK�BN�BP�BR�BW
BS�BO�BJ�BD�BE�BG�BK�BK�BK�BK�BK�BL�BM�BQ�BQ�BT�BXBXBYB[#B]/B\)B_;BbNBdZBgmBhsBhsBhsBjBjBjBjBl�Bm�Bm�Bl�BjBiyBiyBjBjBiyBiyBiyBiyBhsBiyBhsBhsBgmBhsBiyBjBjBk�Bl�Bo�Bq�Bq�Bq�Br�Bs�Bv�Bw�Bw�Bv�Bu�Bu�Bt�Bu�Bu�Bu�Bu�Bu�Bv�Bv�Bv�Bv�Bu�Bu�Bt�Bt�Bt�Bs�Bs�Br�Br�Br�Br�Br�Bt�Bu�Bv�Bv�Bu�Bu�Bu�Bt�Bt�Bs�Bs�Br�Br�Br�Bs�Bs�Bs�Bs�Bs�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bs�Bs�Bt�Bw�By�Bz�B|�B�B�%B�+B�+B�1B�=B�=B�=B�=B�=B�7B�7B�1B�7B�=B�DB�JB�JB�JB�JB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�PB�VB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�bB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�b@� �@�b@��j@���@��@�z�@�{@�l�@���@���@�I�@�I�@��y@�1@�X@�M�@�ƨ@�X@�t�@��R@���@��@���@��@�M�@���@�ff@�ff@�G�@�33@�\)@�?}@���@�+@~5?@{o@y&�@o��@k@a��@^$�@\(�@Sƨ@O\)@M�@K"�@G��@C@=�-@:^5@9��@4�@4��@2J@.{@,�@)��@)hs@'l�@%�@$Z@!��@ bN@�@"�@b@��@��@p�@t�?�b?�(�?�K�?�33?��`?�G�?���?���?�b?�S�?��?�1?��y?�$�?�n�?�  ?�;d?�/?�?��?���?���?�p�?��?��?��j?�M�?~v�?s��?r-?p �?m�h?l��?i��?hr�?f�y?c�
?Z�H?O�;?G�?A%?@  ?;"�?6E�?/�;?,1?%��?"��?�w?�?
=?9X?�?E�?�?�+?��?n�?�D?��?ff?�?
��?C�?	7L?��?J>�dZ>� �>�>��/>�A�>�(�>ؓu>��>ؓu>�b>Ձ>���>Ƨ�>�o>�E�>��F>�->�~�>�A�>��>��u>�t�>��;>�bN>���>��w>�;d>���>�z�>�bN>��^>�%>q��>\(�>M��>E��>@�>:^5>8Q�>7K�>0 �>(��>"��>�->�P>n�>�=��=�F=�"�=���=�j=� �=��
=���=�t�=�O�=�O�=�t�=�t�=���=�^5=�S�=���>I�>$�/>k�>�  >��7>�%>��7>��7>��>w��>vȴ>m�h>gl�>Q�>R�>Z�>\(�>\(�>Z�>Z�>Z�>Y�>Xb>V>S��>Q�>M��>H�9>E��>B�\>A�7>A�7>9X>1&�>/�>+>"��>�>�u>hs>bN>\)>O�>�>o=��=��#=�S�=��`=\=�Q�=���=��-=�t�=�O�=�7L=�o=u=m�h=e`B=T��=H�9=,1=�P=+<�<�`B<�/<�9X<���<e`B<49X<t�;�`B;��
;D��    �o��o�ě���`B�o�49X�T����o��t����㼣�
��j���ͼ�`B�o�C��C��\)�\)��w�8Q�8Q�@��D���D���D���D���H�9�H�9�L�ͽT���aG��ixս�o��%��o�����+��C���hs���P���㽝�-���w�����������
�������T�� Ž�9X��9X��9X��9X��E���^5��j���������
=��/��G���l���h��F���m���   �%�o�����$ݾ$ݾ$ݾ$ݾ+�	7L�I��\)�z��+�����w�$�/�&�y�&�y�&�y�')��6E��<j�@��E�˾H�9�Kƨ�R�Z��bMӾl�D�r�!�u�~�۾���������ƨ��I����;�O߾���bN��hs��녾�n���t����P������"Ѿ��㾝/��;d��Ĝ��Ĝ��Ĝ��G����徣�
��Z��l���xվ��D��KǾ���J�š˾�I���V���`��녾�녾�z��z�և+��b�ؓu�ٙ��ٙ����1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @�z�@��HAp�A!p�AAp�Ac
>A��A�Q�A��A��A��RAиRA�RA�RB \)B\)B\)B\)B \)B(\)B0\)B8\)B@\)BH\)BP\)BW��B`\)Bh\)Bp\)Bx\)B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.B�.C 
C
C
C
C
C

C
C
C
C
C
C�pC
C
C
C
C 
C"
C#�pC&
C(
C*
C,
C.
C0
C2
C4
C6
C8
C:
C<
C>
C@
CB
CD
CF
CH
CJ
CL
CN
CP
CR
CT
CV
CX
CZ
C\
C^
C`
Cb
Cd
Cf
Ch
Cj
Cl
Cn
Cp
Cr
Ct
Cv
Cx
Cz
C|
C~
C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��C��D �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ��DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP��DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV]Bv�Bv�Bu�B}�B��B�_B��B3B%B BB��B	aIB
"�B
�DB
�(B 
B4�BZ�Bg�BaZB^�B]}BX�BSOBI-B<2B5AB'�B�BzB�B
�'B
�B
�B#B
�3B
�`B
�QB
�XB
��B
�^B
�!B
�fB
��B
�"B
��B�B|BB�BBB"�B(�B1�B9�B<�BABFBFBH�BH�BI�BN�BN�BNBR5BVvBW�BX�BYRBTzBR�BK�BF�BG�BL�BM	BM%BM�BMOBM�BO�BRCBSNBU�BX_BX�BY�B\ZB]�B]�Ba<Bc�Be�Bg�Bi\Bi�BjaBj�Bj�Bj�Bj�BmBm�Bm�Bm&Bl0Bk�BkBk�Bj�BjdBjiBj�Bj;Bi�BjBi	Bh�Bh�Bh�Bi<BjZBj�Bk`BmnBo�Br�Br\Br*BsBr�Bv�Bx.BxOBw�Bv�Bv�BuBBvYBv5Bv(BvBu�Bv�Bv�BwBw�BvpBv#Bu�BuBt�BttBt�BsBsBs-BsBr�Bs�BuIBv�Bw6BvZBv+BviBu�Bu�Bt�BtdBsBr�Br�Bs�Bs�BtBtBs�Br�Br�Br�BsMBr�Br�BsDBsBs	Br�Br�Br�Br�Br�Br�Br�Bs�Bs�Bs�Bv�ByTBzB{�B}�B�B��B�8B�'B�<B�YB��B�PB��B��B�,B�0B��B�$B�DB�_B�KB�KB�WB�WB�cB�cB�fB��B��B�wB�uB�]B�TB��B��B�oB��B��B��B��B��B�eB�cB�rB��B�sB��B�tB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�1B��B��B��B��B��B��B��B�B�B��B��B�B�B�5B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�HB�B�B�IB��B��B��B��B��B��B��B��B��B��B��B��B��@�b@� �@�b@��j@���@��@�z�@�{@�l�@���@���@�I�@�I�@��y@�1@�X@�M�@�ƨ@�X@�t�@��R@���@��@���@��@�M�@���@�ff@�ff@�G�@�33@�\)@�?}@���@�+@~5?@{o@y&�@o��@k@a��@^$�@\(�@Sƨ@O\)@M�@K"�@G��@C@=�-@:^5@9��@4�@4��@2J@.{@,�@)��@)hs@'l�@%�@$Z@!��@ bN@�@"�@b@��@��@p�@t�?�b?�(�?�K�?�33?��`?�G�?���?���?�b?�S�?��?�1?��y?�$�?�n�?�  ?�;d?�/?�?��?���?���?�p�?��?��?��j?�M�?~v�?s��?r-?p �?m�h?l��?i��?hr�?f�y?c�
?Z�H?O�;?G�?A%?@  ?;"�?6E�?/�;?,1?%��?"��?�w?�?
=?9X?�?E�?�?�+?��?n�?�D?��?ff?�?
��?C�?	7L?��?J>�dZ>� �>�>��/>�A�>�(�>ؓu>��>ؓu>�b>Ձ>���>Ƨ�>�o>�E�>��F>�->�~�>�A�>��>��u>�t�>��;>�bN>���>��w>�;d>���>�z�>�bN>��^>�%>q��>\(�>M��>E��>@�>:^5>8Q�>7K�>0 �>(��>"��>�->�P>n�>�=��=�F=�"�=���=�j=� �=��
=���=�t�=�O�=�O�=�t�=�t�=���=�^5=�S�=���>I�>$�/>k�>�  >��7>�%>��7>��7>��>w��>vȴ>m�h>gl�>Q�>R�>Z�>\(�>\(�>Z�>Z�>Z�>Y�>Xb>V>S��>Q�>M��>H�9>E��>B�\>A�7>A�7>9X>1&�>/�>+>"��>�>�u>hs>bN>\)>O�>�>o=��=��#=�S�=��`=\=�Q�=���=��-=�t�=�O�=�7L=�o=u=m�h=e`B=T��=H�9=,1=�P=+<�<�`B<�/<�9X<���<e`B<49X<t�;�`B;��
;D��    �o��o�ě���`B�o�49X�T����o��t����㼣�
��j���ͼ�`B�o�C��C��\)�\)��w�8Q�8Q�@��D���D���D���D���H�9�H�9�L�ͽT���aG��ixս�o��%��o�����+��C���hs���P���㽝�-���w�����������
�������T�� Ž�9X��9X��9X��9X��E���^5��j���������
=��/��G���l���h��F���m���   �%�o�����$ݾ$ݾ$ݾ$ݾ+�	7L�I��\)�z��+�����w�$�/�&�y�&�y�&�y�')��6E��<j�@��E�˾H�9�Kƨ�R�Z��bMӾl�D�r�!�u�~�۾���������ƨ��I����;�O߾���bN��hs��녾�n���t����P������"Ѿ��㾝/��;d��Ĝ��Ĝ��Ĝ��G����徣�
��Z��l���xվ��D��KǾ���J�š˾�I���V���`��녾�녾�z��z�և+��b�ؓu�ٙ��ٙ����1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�s<#��<$O<$�P<&�<'��<A�<*��<&O�<k-�<*��<#�<$xx<'\<J�j<G9�<L˴<(L}<(�]<$�E<%�.<(E�<*+.<:7i<;</q�<H}<<��<H�<,+�<:%�<T�X<$<*Gj<8"C<(}�<%�/<B�%<.w#<@��<*?�<%�!<=uk<,>�<%�%<&P	<)&�<,�]</�<(�&<$>q<-]<#��<'/A<*T)<%�X<%ˇ<$�<%|�<%y�<$�v<&��<$�L<'\u<&l3<(@:<,�3<:ڳ<4��<&�<9:�<3�(<Qo<F��<$�<#�<$p�<%,{<%gX<&�"<%Ƚ<%
�<&�?<#��<%o^<$��<#�.<$_Y<$x�<%�<$�<%��<'<&�<%W�<#�1<$�<%�<&��<#��<#��<$R<#�e<$w<#�<#�<$-b<&42<'V <%�K<%�<#�;<$�r<$�
<%�<$Y�<%[<$�<$(/<#��<%9f<$<#��<#ع<#�j<#ر<$�z<#�N<$��<$F�<$(<#��<$��<#ؘ<#�"<$�<$�3<${#<$�{<$�<$(7<$�<#��<#��<#�<#��<#ي<#�\<$L�<$A~<#��<%�<#�<#��<$UT<$��<#��<#��<$\<#� <#�<$u <#��<#�O<$�<$)�<$�<$8�<$u�<$c�<$�<<$B<<#��<#��<#��<#�z<#ٯ<#��<#��<#�!<#�<#�L<#�<$0]<#�<#�<$&e<$ 7<#��<#�<#�<#�<#�<#��<#׋<#�1<#�A<#�<$>�<$}�<$�<$@�<$�<+��<$�%<#�0<#�;<#�<#ׅ<#�G<$m<#�H<$z<#��<$��<#׈<#�D<#�i<#צ<#ڬ<#ײ<#׮<#�'<#�F<#ۡ<#ۻ<#�;<#�<#�t<#�d<#��<#�?<#�	<#��<#�<#ۥ<#�<#�G<#�u<#�*<#��<#٥<#�J<#�S<#��<#܆<#�<#��<$<$�<#��<#�<#�\<#�W<#��<#�Q<#�n<#�\<#�9<#��<#��<#��<#�y<#�]<#�P<#�<#�<#�=<#ٓ<#�<#ߞ<#�<#�U<#ܑ<#�A<#ۼ<#��<#��<#��<#ۺ<#۟<#�D<#�X<#��<#��<#ߡ<#�<#�D<#�X<#��<#��<#�8<#�<#ۄ<#׺<#�<#��<#�<#�<#װ<#�i<#�,<#ף<#ף<#׬<#�<#ױ<#�=<#ۻ<#��<#�Q<#��<#�<#�<#�2<#�D<#��<#�<#�<#�5<#�1<#�0<#�<#״<#�<#�<#ۿ<#��<#۷<#ע<#ע<#׬<#٤<#�L<#�L<#�S<#�~<#�t<#�<#�h<#�<#�5<#�U<#�<#�a<#�0<#�A<#ۚ<#�;<#٭<#ط<#ץ<#ס<#׫<#�0<#۽<#�5<#�x<#�<#�<#ߑ<#��<#�<#��<#׮<#׫<#�;<#ܴ<$#�<#�<#��<#�<#�<#��<#��<#�2<#��<$�<#�:<#�D<$*<$.<$'X<#ۛ<#�<<#�/<#�1<#�z<#�-<#��<#�A<#�E<#�]<#��<#�<#�<#�d<#�><#�W<#��<#�d<#׫<#�G<#��<#۶<#ٖ<#�<#�t<#��<$�w<$5\<#��<#�<$6�<#�<#�s<#��<#��<#�x<#��<#�<#�<#�b<#�w<#װ<#��<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.09                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929160807201909291608072019092916081120190929160818202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@���@���@���G�O�G�O�G�O�G�O�DV�DV�DV�G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          