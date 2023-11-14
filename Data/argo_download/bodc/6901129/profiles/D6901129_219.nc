CDF       
      	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       	DATE_TIME         N_PROF        N_PARAM       N_LEVELS  f   N_CALIB       	N_HISTORY            	   title         Argo float vertical profile    institution       BODC   source        
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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  D�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  JP   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 h  O�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 h  QP   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 h  R�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  T    PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  Y�   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  _P   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 h  d�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 h  fP   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 h  g�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  i    PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  n�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  tP   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  y�   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    zH   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �H   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �H   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �H   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  �   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �0   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �L   HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  �Argo profile    3.1 1.2 19500101000000  20210225045048  20210225045048  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125505                          2C  D   APEX                            6229                            120210                          846 @�eMt>2�1   @�eMt>2�@P睲-V�5h�9Xb1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            !A   A   A   @��@y��@�  A   A!��AA��A`  A�  A�  A�  A�  A�  A�  A�  A���B   B  B  B  B ffB(  B/��B7��B@  BH  BP  BXffB`  Bh  Bp  BxffB�33B�  B�  B�  B�33B�  B���B���B�  B�  B���B���B���B�  B�33B�33B�  B�  B�  B�33B�  B���B�  B�  B�  B�33B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C�fC  C�fC�fC  C   C"  C$  C&  C(  C*  C,  C.  C/�fC2  C4  C6  C8  C9�fC<  C>�C@  CB  CD  CE�fCH  CJ  CL  CN  CP  CR  CT  CV  CX�CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl�Cn  Co�fCq�fCs�fCu�fCx  Cz  C|  C~  C�fC�  C�  C��3C�  C��C�  C�  C�  C�  C��3C�  C�  C��3C�  C�  C��C��C�  C��3C�  C�  C��3C�  C��C��C�  C�  C�  C�  C��C�  C�  C�  C��C�  C��3C�  C�  C�  C��3C�  C��C�  C�  C��C�  C�  C�  C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C��3C�  C�  C��C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C��C�  C��3C�  C��C�  C�  C��3C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C��3C�  C�  C�  C��C�  C�  C��3D � D  D� D  D� D  Dy�D  D� D  D�fD  D� DfD� D  D� D	  D	� D
fD
� D  Dy�D  D� D  Dy�D  D� D��D� DfD� D��D� D  D� D  D�fDfD� D  D� D��Dy�D��Dy�D  D�fD  D� D  Dy�D  D� D  D� D  D� D��Dy�D��Dy�D��D � D!fD!� D"  D"� D#  D#� D#��D$y�D%  D%�fD&fD&�fD'fD'� D(  D(� D)  D)� D*  D*y�D*��D+y�D,  D,�fD-  D-y�D.  D.� D/  D/� D0  D0�fD1  D1y�D1��D2� D3�B�B�B�B�B� B~�B}�B|�Bz�By�Bt�Bn�Br�Bs�Bo�Bm�Bp�Bn�Bn�Bk�BiyBhsBhsBgmBffBdZBcTBaHB_;B`BB]/BZBS�BP�BL�BH�BF�BD�BB�B?}B<jB:^B6FB6FB49B.B,B,B'�B'�B(�B&�B&�B#�B#�B�B�B�BuBoBhBhB�B�B �B.B5?B6FB9XB:^B:^B;dB=qB@�BC�BD�BG�BK�BN�BO�BO�BP�BP�BQ�BR�BS�BS�BT�BVBW
BYB\)B^5B_;B`BBaHBbNBe`BhsBiyBk�Bm�Bm�Bn�Bn�Bo�Bo�Bq�Br�Bs�Bu�Bu�Bv�Bw�Bx�By�Bz�B{�B|�B}�B}�B~�B�B�B�B�B�B�%B�%B�%B�%B�1B�1B�1B�=B�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�=B�7B�=B�=B�=B�=B�=B�=B�7B�7B�7B�=B�7B�=B�=B�=B�7B�=B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�\B�\B�\B�bB�bB�bB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�{B�{B��B��B��B��B��B�{B�{B��B��B��B��B��B��B��B�{B�{B�{B�{B��B��B��B�{B��B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�+@�"�@�33@�C�@�K�@�K�@�S�@�o@��^@��9@���@�9X@��@�"�@���@��-@���@�b@��@���@�@��@��@�G�@l�@}�@|�j@zJ@w�@t�@sS�@q��@m�T@d�/@Y��@Q&�@L�D@HĜ@D�@@Q�@<�D@7+@2-@*J@!��@!��@`B@��@O�@��@�7?�?�u?��?͑h?��h?���?�  ?�  ?kC�?Y�#?T�j?G�?@  ?5?)�^?%`B?$Z?#o?"�\?"M�?"J? Ĝ?|�?5??�?�-?�?(�?��?dZ?�H?�?�#?��?��?�#?X?Q�?ȴ?�+?ȴ?K�?K�?�P?��?��?^5?�-?�R? A�?"M�?"��?#o?#��?#��?#��?$Z?%�?$�/?$��?#��?"�\?"��?"J?!�7?!�7?!G�? �? A�?|�?v�?j?��?dZ?��?(�?(�?(�?(�?j?^5?�+?��?��?bN?�;?�?�?��?V?{?V?	x�?�>��H>�9X>�->�V>�~�>��>���>߾w>�5?>�(�>�"�>�"�>׍P>�\)>ɺ^>�$�>ě�>Õ�>��7>�X>���>�~�>�`B>��w>�"�>���>�\)>�V>���>��\>�  >vȴ>n��>bM�>O�;>I�^>A�7>@�>@�>>v�>2->.{>-V>)��>%�T>�w>�>��>��>��>��>�>n�>n�>n�>hs>\)>
=q>�>   =��m=�F==�=�l�=�S�=�/=�
==��=Ƨ�=�^5=���=��w=���=���=���=�\)=�%=ix�=D��=@�=8Q�=,1=#�
=t�=C�<��<�h<�/<ě�<ě�<���<�/<�j<��
<�t�<�C�<D��<#�
;�`B:�o��o�D����o���
��`B��`B�o�o�49X�T����t���1��1��1��9X��9X���ͽt��,1�Y��ixսu�y�#��%��o��7L��O߽�hs��hs�������㽟�w��{��9X��Q콺^5��j������
=��/��`B����h����J�
=q�\)������w�!���"��$�/�$�/�%�T�%�T�%�T�)��,1�,1�-V�.{�/��1&�1&�2-�5?}�;dZ�B�\�I�^�T���]/�`A��ixվo���r�!�t�j�~�۾�J�����hs�����b���������"Ѿ�����5?��;d��Ĝ��G���MӾ��徢�徢�徢MӾ�MӾ��徢�徢�徢MӾ�MӾ��
���/��Z��Z���
���
���
��Z���
���
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @ 6�@�H@�N{A�=A#@�AC@�Aa�=A�ӟA�ӟA�ӟA�ӟA�ӟA�ӟA�ӟA�lB i�Bi�Bi�Bi�B �5B(i�B0iB8iB@i�BHi�BPi�BX�5B`i�Bhi�Bpi�Bx�5B�hB�4�B�4�B�4�B�hB�4�B��B��B�4�B�4�B��B��B��B�4�B�hB�hB�4�B�4�B�4�B�hB�4�B��B�4�B�4�B�4�B�hB�4�B�4�B�4�B�4�B�4�B�4�C tCtCtCtCtC
tCtCtCtCtCtC �CtC �C �CtC tC"tC$tC&tC(tC*tC,tC.tC0 �C2tC4tC6tC8tC: �C<tC>4C@tCBtCDtCF �CHtCJtCLtCNtCPtCRtCTtCVtCX4CZtC\tC^tC`tCbtCdtCftChtCjtCl4CntCp �Cr �Ct �Cv �CxtCztC|tC~tC� mC�:C�:C� mC�:C�C�:C�:C�:C�:C� mC�:C�:C� mC�:C�:C�C�C�:C� mC�:C�:C� mC�:C�C�C�:C�:C�:C�:C�C�:C�:C�:C�C�:C� mC�:C�:C�:C� mC�:C�C�:C�:C�C�:C�:C�:C�:C�:C� mC� mC�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�C�:C�:C�:C� mC�:C�:C�C�C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�C�:C�:C�:C�:C�C�:C� mC�:C�C�:C�:C� mC� mC�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�:C�C�:C�:C� mC�:C�:C�:C�C�:C�:D  6D ��D�D��D�D��D�D�7D�D��D�D�D�D��DD��D�D��D	�D	��D
D
��D�D�7D�D��D�D�7D�D��D 7D��DD��D 7D��D�D��D�D�DD��D�D��D 7D�7D 7D�7D�D�D�D��D�D�7D�D��D�D��D�D��D 7D�7D 7D�7D  7D ��D!D!��D"�D"��D#�D#��D$ 7D$�7D%�D%�D&D&�D'D'��D(�D(��D)�D)��D*�D*�7D+ 7D+�7D,�D,�D-�D-�7D.�D.��D/�D/��D0�D0�D1�D1�7D2 7D2��D3 7B�B��B��B��B� B~�B~\B~�B|HBy�Bu�By�BycBtWBq^Bs$Br�Bo�BpBk�BkdBk�Bi�Bi�Bg}Be>Be[BciBabB`�B^_B]-BZBX�BR�BLCBI�BG�BE8BBOB@�B>B<5B;�B4�B1FB2oB1�B,7B,GB.B,B)�B+OB)�B&]B!B�B'B�B�B�BB�B#B.�B5nB6�B9mB:iB:oB;�B=�B@�BC�BD�BG�BK�BN�BO�BO�BQBP�BQ�BR�BS�BTBU.BVMBWBYB\B^1B_0B`4BaBbBd�Bh?Bi0Bk%BmqBm�BnBn�Bo�BoBq�Br�Bs�Bu�Bu�Bv�Bw�Bx�By�Bz�B|B|�B~B~*B\B�-B�B�B�B�!B�#B�%B�B��B��B��B��B�{B�]B�vB�GB�QB�QB�PB�uB��B�FB��B��B�tB��B��B�{B��B��B�fB�nB�VB�?B��B�B��B��B�bB�XB�tB��B� B��B��B��B��B��B��B�ZB��B��B��B��B��B��B�B��B��B�JB�?B�YB��B�qB�PB�jB�xB��B�iB�RB�EB�QB�SB�sB�kB�JB�KB�XB�dB��B��B�B�hB��B�iB�iB�jB�jB�wB�{B��B��B��B��B��B�VB�dB��B��B��B��B��B�vB��B��B��B��B��B��B��B��B��B�rB�fB�_B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��@�+@�"�@�33@�C�@�K�@�K�@�S�@�o@��^@��9@���@�9X@��@�"�@���@��-@���@�b@��@���@�@��@��@�G�@l�@}�@|�j@zJ@w�@t�@sS�@q��@m�T@d�/@Y��@Q&�@L�D@HĜ@D�@@Q�@<�D@7+@2-@*J@!��@!��@`B@��@O�@��@�7?�?�u?��?͑h?��h?���?�  ?�  ?kC�?Y�#?T�j?G�?@  ?5?)�^?%`B?$Z?#o?"�\?"M�?"J? Ĝ?|�?5??�?�-?�?(�?��?dZ?�H?�?�#?��?��?�#?X?Q�?ȴ?�+?ȴ?K�?K�?�P?��?��?^5?�-?�R? A�?"M�?"��?#o?#��?#��?#��?$Z?%�?$�/?$��?#��?"�\?"��?"J?!�7?!�7?!G�? �? A�?|�?v�?j?��?dZ?��?(�?(�?(�?(�?j?^5?�+?��?��?bN?�;?�?�?��?V?{?V?	x�?�>��H>�9X>�->�V>�~�>��>���>߾w>�5?>�(�>�"�>�"�>׍P>�\)>ɺ^>�$�>ě�>Õ�>��7>�X>���>�~�>�`B>��w>�"�>���>�\)>�V>���>��\>�  >vȴ>n��>bM�>O�;>I�^>A�7>@�>@�>>v�>2->.{>-V>)��>%�T>�w>�>��>��>��>��>�>n�>n�>n�>hs>\)>
=q>�>   =��m=�F==�=�l�=�S�=�/=�
==��=Ƨ�=�^5=���=��w=���=���=���=�\)=�%=ix�=D��=@�=8Q�=,1=#�
=t�=C�<��<�h<�/<ě�<ě�<���<�/<�j<��
<�t�<�C�<D��<#�
;�`B:�o��o�D����o���
��`B��`B�o�o�49X�T����t���1��1��1��9X��9X���ͽt��,1�Y��ixսu�y�#��%��o��7L��O߽�hs��hs�������㽟�w��{��9X��Q콺^5��j������
=��/��`B����h����J�
=q�\)������w�!���"��$�/�$�/�%�T�%�T�%�T�)��,1�,1�-V�.{�/��1&�1&�2-�5?}�;dZ�B�\�I�^�T���]/�`A��ixվo���r�!�t�j�~�۾�J�����hs�����b���������"Ѿ�����5?��;d��Ĝ��G���MӾ��徢�徢�徢MӾ�MӾ��徢�徢�徢MӾ�MӾ��
���/��Z��Z���
���
���
��Z���
���
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<#�:<#�,<#�A<#�<#נ<#�7<$�<&��<%~�<#�+<$��<u�M<DD�<$4�<&^g<:��<(+�<%u�<%�&<$,<&�B<,��<%f<'�P<$�4<$�\<'3><'��<'�<$#�<%
�<+C�<B��<P�j<=.�<-��<*��<,�^<)��<*,�<0��<.�N<=��<:e�<#�j<+�(<A�m<;k�<1�[<2L�<7/E<8�<)��<KV�<>�C<L�$<,|f<?uE<bGn<-�<$�<(ȶ<%�s<'<(*�<$h�<#�<#�Z<#�*<#�)<#�I<#��<#��<#�<#��<#�^<#��<#�{<#��<#�Y<#�Y<#��<#�F<#�|<#��<#�<#�?<#�h<#��<#�i<#�&<#�k<#�t<#�<#�<#ڈ<#�i<$<#�<#��<#�<#��<#�<#�H<#�<#�{<#׶<#؅<#�T<#�[<#��<#�<#�)<#�1<#�<#��<#ٗ<#��<#�<#��<#��<#��<#�><#ُ<#�<#�<<#��<#ח<#��<#�:<#��<$>�<#��<$6�<#�Z<#�+<#�9<#�<#ه<#٫<#�W<#��<$=t<$�n<%Is<$CO<#�4<$l<#��<#�<$W<#�1<#��<#�<#�B<#��<#��<$o<$`<#��<#ߢ<#ܧ<#��<$V{<$�<<$�<$<$�<$|<$4�<$<#�A<#��<$\Q<#��<$a<$x<$55<$�+<#�<$ �<#ن<#��<#�<$)D<#�P<#�t<#��<#��<#�<#ߑ<#ٱ<#��<#َ<#��<#�<#�?<#��<#��<#ٵ<#܆<#�<#�<#�<#��<#��<#�E<#�E<#�F<#�b<#�E<#ߴ<#�m<#�-<#��<#��<#�P<#�<#��<#�8<#�<#��<#��<$�<#��<#�<#�%<#��<#�<#��<#߱<#��<#��<#�<#�<#�<#�&<#�p<#�w<#�D<#��<#��<#�d<#�:<#�n<#��<#�u<#ـ<#و<#��<#��<#�h<#�V<#�Y<#�`<#�<#ߗ<#��<#��<#�q<#��<#�Y<$<#��<$p<#�<#�P<#�K<#��<#�z<#�L<#�]<#ܸ<#��<#��<#�<#��<#��<#�}<#�6<#٩<#�<<#�}<#�!<#�<#�t<#ޠ<#��<#��<$<#�V<#�<#��<#��<#�F<#��<#٢<#��<#��<#�h<#��<#�<#��<#��<#��<#َ<#ِ<#ٿ<#۲<#��<#ٱ<#�<#�P<#�N<#��<$g<$ 4<#�<$	h<#�><#ߩ<#�<$�<#�<#��<$�<$<#�\<#�$<#٦<#��<#�6<#�2<#�><#߈<#��<#�<#ق<#��<#��<#�<#��<#�\<#��<#��<#�<#��<#�e<#��<#�<#��<#�<#��<#��<#�Y<#�<#؋<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.10333                                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929163915201909291639152019092916391920190929163926202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@��@��@��G�O�G�O�G�O�G�O�D3�D3�D3�G�O�G�O�G� �G� �G� �G� �G� �G� �G� �                                6389758         6389758         131072                                          