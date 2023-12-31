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
resolution        ?�������     �  ?    PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  F�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  N�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  Vl   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  X`   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  ZT   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     units         decibar    	valid_min                    	valid_max         @�p        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  \H   PSAL_ADJUSTED            
      
   	long_name         Practical salinity     standard_name         sea_water_salinity     units         psu    	valid_min         @          	valid_max         @D�        conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  d   TEMP_ADJUSTED            
      
   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      units         degree_Celsius     	valid_min         �         	valid_max         @D         conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  k�   PRES_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PRES_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  s�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   standard_name         PSAL_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  u�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   standard_name         TEMP_ADJUSTED_QC   conventions       Argo reference table 2     
_FillValue                 �  w|   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PRES_ADJUSTED_ERROR    units         decibar    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %7.1f      FORTRAN_format        F7.1   
resolution        ?�������     �  yp   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         PSAL_ADJUSTED_ERROR    units         psu    conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  �4   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     standard_name         TEMP_ADJUSTED_ERROR    units         degree_Celsius     conventions       Argo reference table 2     
_FillValue        G�O�   C_format      %9.3f      FORTRAN_format        F9.3   
resolution        ?PbM���     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    source_name       	PARAMETER      conventions       Argo reference table 3     
_FillValue                  `  ��   SCIENTIFIC_CALIB_EQUATION               	             	long_name         'Calibration equation for this parameter    source_name       SCIENTIFIC_CALIB_EQUATION      
_FillValue                    �   SCIENTIFIC_CALIB_COEFFICIENT            	             	long_name         *Calibration coefficients for this equation     source_name       SCIENTIFIC_CALIB_COEFFICIENT   
_FillValue                    �   SCIENTIFIC_CALIB_COMMENT            	             	long_name         .Comment applying to this parameter calibration     source_name       SCIENTIFIC_CALIB_COMMENT   
_FillValue                    �   SCIENTIFIC_CALIB_DATE               	            	long_name         Date of calibration    source_name       SCIENTIFIC_CALIB_DATE      conventions       YYYYMMDDHHMISS     
_FillValue                  T  �   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     source_name       HISTORY_INSTITUTION    conventions       Argo reference table 4     
_FillValue                    �p   HISTORY_STEP                     	long_name         Step in data processing    source_name       HISTORY_STEP   conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    source_name       HISTORY_SOFTWARE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     source_name       HISTORY_SOFTWARE_RELEASE   conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      source_name       HISTORY_REFERENCE      conventions       Institution dependent      
_FillValue                 �  ��   HISTORY_DATE                     	long_name         #Date the history record was created    source_name       HISTORY_DATE   conventions       YYYYMMDDHHMISS     
_FillValue                  d  ��   HISTORY_ACTION                       	long_name         Action performed on data   source_name       HISTORY_ACTION     conventions       Argo reference table 7     
_FillValue                    �   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   source_name       HISTORY_PARAMETER      conventions       Argo reference table 3     
_FillValue                  p  �    HISTORY_START_PRES                    	long_name          Start pressure action applied on   source_name       HISTORY_START_PRES     units         decibar    
_FillValue        G�O�        ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    source_name       HISTORY_STOP_PRES      units         decibar    
_FillValue        G�O�        ��   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    source_name       HISTORY_PREVIOUS_VALUE     
_FillValue        G�O�        ��   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   source_name       HISTORY_QCTEST     conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                  p  ��Argo profile    3.1 1.2 19500101000000  20210225042553  20210225042553  6901129 Argo UK                                                         Jon Turton                                                      PSAL            TEMP            PRES               �A   BO  125454                          2C  D   APEX                            6229                            120210                          846 @���>��@1   @���>��@@P�\(��4i�^5?}1   GPS     Primary sampling: mixed                                                                                                                                                                                                                                            'A   A   A   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ�fDK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP�fDQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]fD]� D^  D^� D_  D_� D`  D`� DafDa� Da��Db� Dc  Dc� Dd  Dd� De  De� Df  Df� DgfDg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dky�Dk��Dl� DmfDm�fDn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Du  Du� Dv  Dv� Dw  Dw� Dx3Dx@ B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B  BhB�B�B�B)�BO�Be`Bq�Bv�Bx�Bx�Bx�Bw�Bv�Bv�Bv�Bv�Bu�Bt�Bt�Br�Br�Br�Br�Br�Br�Br�Br�Bq�Bq�Bp�Bp�Bp�Bq�Br�Bs�Bs�Bs�Bs�Bt�Bt�Bu�Bv�Bv�Bw�Bw�Bw�Bx�Bx�By�By�By�By�By�B{�B}�B}�B}�B}�B}�B~�B~�B~�B~�B~�B~�B~�B~�B~�B~�B~�B~�B~�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�+B�+B�+B�+B�+B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�%B�%B�%B�%B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�1B�1B�1B�1B�1B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�JB�PB�PB�PB�PB�PB�PB�VB�VB�VB�VB�VB�VB�VB�VB�\B�\B�\B�\B�\B�bB�bB�bB�hB�hB�oB�oB�uB�uB�uB�{B�{B�uB�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�\)@ �9@ Ĝ@ Ĝ@ ��@ ��@ ��@ Ĝ@ ��@ ��@ Ĝ@ Ĝ@ Ĝ@ �9@ �9@ Ĝ@ Ĝ@ Ĝ@ �9@ �9@ �9@ �9@ ��@   ?��w?���?�|�?�|�?�\)?��?��?�v�?�V?�5??�p�?�O�?�/?��?�1?�~�?���?�-?��?��?�7L?��?��?Ƨ�?��?��?���?j��?Z�?S��?Q�?Q&�?O�?M�h?G�?<�?9X?7�P?1��?,I�?*��?#o?X??�?V?I�??�y?Z?S�?J? �?S�?S�>�j>���>�F>� �>�->�!>��>��>��>��>>>��>� �>�&�>�->��>>��
>�Z>��T>�;d>�"�>��>Ձ>�n�>�V>���>�1'>ě�>��7>���>��m>��#>�?}>�33>�33>��!>���>��h>�1>�>���>�x�>���>�~�>�ff>���>��w>��w>�/>��>�>�z�>�t�>��`>��;>��`>��+>��+>�n�>�bN>�V>�ƨ>�C�>�C�>�C�>��^>���>�  >x��>q��>o��>gl�>Z�>Q�>F��>C��>>v�>5?}>/�>)��>&�y>%�T>#�
>!��> Ĝ>!��>#�
>%�T>$�/>$�/>$�/>$�/>!��>��>�P>�+>t�>bN>I�>
=q>+>�>�>�>J=��=��#=��#=���=���=���=�=��=�h=�x�=�x�=�l�=�S�=�;d=�"�=���=��`=ȴ9=��=��=�j=�j=�j=�j=�j=�^5=�^5=�Q�=�Q�=�^5=�^5=�Q�=�E�=�9X=�1=��=��T=���=�hs=�7L=u=e`B=Y�=H�9=@�=@�=H�9=P�`=T��=]/=}�=y�#=ix�=]/=T��=D��=@�=0 �=+<�h<�<�`B<�`B<�h=o=\)=��=#�
=,1=0 �='�=,1=,1=,1=�w=t�=\)=\)=+=+=o<�<�`B<�/<���<�j<�1<���<�t�<e`B<T��<T��<T��<D��<49X<#�
<t�;ě�;o:�o�o��o�t��#�
�49X�u��C���C���t���1�ě�������`B��h�����o�+��P�#�
�,1�<j�H�9�L�ͽY��]/�aG��aG��e`B�ixսixսu�y�#��o�����+��7L��7L��C���O߽�\)��hs��hs��t����P���㽧��{��{��-��^5����ě��������ͽ�����ͽ������
=��S���xս��F���#�1'�I��O߾\)�n��t���+���������R��w� Ĝ�!���$�/�&�y�)��-V�49X�>vɾD���E�˾E�˾M��R�W
=�Z��["Ѿ]/�_;d�`A��aG��cS��cS��e`B�fff�ixվj~��k��o���r�!�u�w�پ}󶾃�������7L���9���9��C����;�V��bN��녾�n����Ͼ������������㾜(����R���R��;d��G���ff���y������D�������������׾�&龱&龰�׾�&龱&龱��������-��9X���j��9X��?}���پ�Q쾺^5���m���۾�  ����\��$ݾ�$ݾ�$ݾǮ�ȴ9��7L��=q��C����������;���`��hs��hs��녾�n���n���n���n���n���n���t������׍P�ؓu�ٙ��������ٙ�����ڟ��ڟ�����ڟ��ڟ��ڟ��ڟ��ڟ��ڟ��ڟ��������������ٙ��ٙ��ٙ��ٙ��ٙ�������11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @J=q@�Q�@�Q�A(�A$(�AD(�Ad(�A�{A�{A�{A�{A�{A�{A�{A�{B
=B	
=B
=B
=B!
=B)
=B1
=B9
=BA
=BI
=BQ
=BY
=Ba
=Bi
=Bq
=By
=B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BąBȅB̅BЅBԅB؅B܅B��B�B�B�B��B�B��B��C B�CB�CB�CB�CB�C
B�CB�CB�CB�CB�CB�CB�CB�CB�CB�CB�C B�C"B�C$B�C&B�C(B�C*B�C,B�C.B�C0B�C2B�C4B�C6B�C8B�C:B�C<B�C>B�C@B�CBB�CDB�CFB�CHB�CJB�CLB�CNB�CPB�CRB�CTB�CVB�CXB�CZB�C\B�C^B�C`B�CbB�CdB�CfB�ChB�CjB�ClB�CnB�CpB�CrB�CtB�CvB�CxB�CzB�C|B�C~B�C�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�.C�.C�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HC�!HD �D ��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D	�D	��D
�D
��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D�D��D �D ��D!�D!��D"�D"��D#�D#��D$�D$��D%�D%��D&�D&��D'�D'��D(�D(��D)�D)��D*�D*��D+�D+��D,�D,��D-�D-��D.�D.��D/�D/��D0�D0��D1�D1��D2�D2��D3�D3��D4�D4��D5�D5��D6�D6��D7�D7��D8�D8��D9�D9��D:�D:��D;�D;��D<�D<��D=�D=��D>�D>��D?�D?��D@�D@��DA�DA��DB�DB��DC�DC��DD�DD��DE�DE��DF�DF��DG�DG��DH�DH��DI�DI��DJ�DJ�
DK�DK��DL�DL��DM�DM��DN�DN��DO�DO��DP�DP�
DQ�DQ��DR�DR��DS�DS��DT�DT��DU�DU��DV�DV��DW�DW��DX�DX��DY�DY��DZ�DZ��D[�D[��D\�D\��D]
D]��D^�D^��D_�D_��D`�D`��Da
Da��Db
>Db��Dc�Dc��Dd�Dd��De�De��Df�Df��Dg
Dg��Dh�Dh��Di�Di��Dj�Dj��Dk�Dk�>Dl
>Dl��Dm
Dm�
Dn�Dn��Do�Do��Dp�Dp��Dq�Dq��Dr�Dr��Ds�Ds��Dt�Dt��Du�Du��Dv�Dv��Dw�Dw��Dx#�DxP�B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
��B
�UB
�B
��B
��B
��B
��B
��B
��B
�B
��B
��B
�(B
��B
��B
��B
�3B
�rB
��B
��B
��B
�/BpB�B)B�B%B6WBWYBj�Bt�Bw�By8Bx�By7Bx#Bw�Bx�BwvBw(Bv�Bu�BuBt)Bt~Bs_BswBs]BsBr�BsoBr'Bq�Bp�Bp�Bp'Bq�Bs�BtBt1Bt
Bs�Bt�Bt�Bu�Bv�Bv�Bw�Bw�Bw�Bx�Bx�By�By�By�Bz�By�B{�B~�B~WB~B~aB~DB^BRB;BQBFB:BEB0BhB+B~�B
B_B�B�+B�#B�"B�B� B�B�sB�jB�]B�B�XB�|B�dB�?B�4B�VB�4B�B��B�%B��B�_B�^B�jB�<B�2B�2B�UB��B��B��B��B�IB��B��B��B��B�OB�eB��B�jB�^B�EB�-B�:B�=B�1B�B�B�B�7B�,B�+B�-B�PB��B�TB�9B�QB�RB�_B�IB�VB�JB�>B�6B�NB�]B�PB�9B�CB�8B�8B�DB�PB�SB�VB�>B�JB�VB�XB�]B�iB�_B�vB�uB�HB�_B�KB�JB�JB�KB�VB�MB�\B�PB�EB�PB�]B�^B�bB��B�pB�cB�rB��B��B��B��B��B��B�uB�[B�HB�JB�WB�LB�
B�wB��B��B��B��B��B��B��B��B�rB��B�B�tB�^B�aB�hB�yB�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B�DB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�:B�	B��B��B�"B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�.B�!B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B��B��B��B��B�;B��B��B�	B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��?�\)@ �9@ Ĝ@ Ĝ@ ��@ ��@ ��@ Ĝ@ ��@ ��@ Ĝ@ Ĝ@ Ĝ@ �9@ �9@ Ĝ@ Ĝ@ Ĝ@ �9@ �9@ �9@ �9@ ��@   ?��w?���?�|�?�|�?�\)?��?��?�v�?�V?�5??�p�?�O�?�/?��?�1?�~�?���?�-?��?��?�7L?��?��?Ƨ�?��?��?���?j��?Z�?S��?Q�?Q&�?O�?M�h?G�?<�?9X?7�P?1��?,I�?*��?#o?X??�?V?I�??�y?Z?S�?J? �?S�?S�>�j>���>�F>� �>�->�!>��>��>��>��>>>��>� �>�&�>�->��>>��
>�Z>��T>�;d>�"�>��>Ձ>�n�>�V>���>�1'>ě�>��7>���>��m>��#>�?}>�33>�33>��!>���>��h>�1>�>���>�x�>���>�~�>�ff>���>��w>��w>�/>��>�>�z�>�t�>��`>��;>��`>��+>��+>�n�>�bN>�V>�ƨ>�C�>�C�>�C�>��^>���>�  >x��>q��>o��>gl�>Z�>Q�>F��>C��>>v�>5?}>/�>)��>&�y>%�T>#�
>!��> Ĝ>!��>#�
>%�T>$�/>$�/>$�/>$�/>!��>��>�P>�+>t�>bN>I�>
=q>+>�>�>�>J=��=��#=��#=���=���=���=�=��=�h=�x�=�x�=�l�=�S�=�;d=�"�=���=��`=ȴ9=��=��=�j=�j=�j=�j=�j=�^5=�^5=�Q�=�Q�=�^5=�^5=�Q�=�E�=�9X=�1=��=��T=���=�hs=�7L=u=e`B=Y�=H�9=@�=@�=H�9=P�`=T��=]/=}�=y�#=ix�=]/=T��=D��=@�=0 �=+<�h<�<�`B<�`B<�h=o=\)=��=#�
=,1=0 �='�=,1=,1=,1=�w=t�=\)=\)=+=+=o<�<�`B<�/<���<�j<�1<���<�t�<e`B<T��<T��<T��<D��<49X<#�
<t�;ě�;o:�o�o��o�t��#�
�49X�u��C���C���t���1�ě�������`B��h�����o�+��P�#�
�,1�<j�H�9�L�ͽY��]/�aG��aG��e`B�ixսixսu�y�#��o�����+��7L��7L��C���O߽�\)��hs��hs��t����P���㽧��{��{��-��^5����ě��������ͽ�����ͽ������
=��S���xս��F���#�1'�I��O߾\)�n��t���+���������R��w� Ĝ�!���$�/�&�y�)��-V�49X�>vɾD���E�˾E�˾M��R�W
=�Z��["Ѿ]/�_;d�`A��aG��cS��cS��e`B�fff�ixվj~��k��o���r�!�u�w�پ}󶾃�������7L���9���9��C����;�V��bN��녾�n����Ͼ������������㾜(����R���R��;d��G���ff���y������D�������������׾�&龱&龰�׾�&龱&龱��������-��9X���j��9X��?}���پ�Q쾺^5���m���۾�  ����\��$ݾ�$ݾ�$ݾǮ�ȴ9��7L��=q��C����������;���`��hs��hs��녾�n���n���n���n���n���n���t������׍P�ؓu�ٙ��������ٙ�����ڟ��ڟ�����ڟ��ڟ��ڟ��ڟ��ڟ��ڟ��ڟ��������������ٙ��ٙ��ٙ��ٙ��ٙ�������11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��<$��<#��<#ۦ<#�|<#��<#��<#��<#ُ<#��<#�<#��<#�K<#�<#�Q<#�g<#��<#��<#��<#��<#��<#��<#�Q<$$�<#�Z<#�M<#�<#�<#�H<#�<#�:<#�<#ߚ<#��<#�R<#��<#ތ<#��<$ <$J�<'��<$��<$_(<-�<:�z<(�<&.0<,%U<9�<�P�<LT
<:�i<+�<%:�<$r<#�E<$<<$<%�<'�q<$am<$6<%�<$�E<$f<%�<&�<$eF<$�<$bb<$�<#��<${<$+"<#�<#�~<#��<#��<#ޢ<$��<$�<$(�<$i<#�<#� <#�<#�<#�V<#�<#ޟ<#�-<#ٯ<#׺<#�<#��<#�*<#�K<% <#��<#�n<$M<$�<#�<$m<$'<$�<$	]<#�<$�<$ =<#�\<#��<#�9<$ <#�<#ܥ<#�<$<#�<#�<#�<#�<#ߠ<#ي<#��<$�<$	�<#��<#ܽ<#��<$Q<#��<#��<#�<#�K<#�+<#�r<#��<#��<$�<#�,<#�<#��<#�G<#�G<#�><#��<$(�<$"<$
G<$0<#�1<$l<$T�<$<$4�<#��<#��<$l<$ )<#��<#��<#�7<#��<#��<#ߢ<#ْ<#�<#��<#�[<#�E<#�+<#�<#�<$<#�<#�M<#��<#�&<#�`<#��<#�<#�B<#��<#�D<#�<#�<#�"<#�y<#�<#�0<#�A<#��<#�.<#�,<#��<#�g<#��<#�0<#��<#�<#�<#��<#�<#�s<#�3<#�<#�Q<#�<#�<#�9<#�t<#��<#ߔ<#�<#ٞ<#�<#ߪ<#�<#�M<#�'<#�_<#��<#�A<$g<#��<$	+<#��<#�><#�<#�*<#��<#ץ<#��<#پ<#�y<#��<#�4<#��<#��<#�<#�f<#ߢ<#�<$*5<#�<#�<#�<#ۜ<#�H<#�!<#�<#�<#פ<#��<#�$<#�2<#٩<#�<#��<#� <#�R<#ߠ<#ܜ<#�<#�b<#߱<#�4<#�<#��<#��<#�[<#�C<#�	<#�4<#��<#�|<#�<#�%<#�\<#�@<#ߡ<#��<#�}<#�P<#��<#�9<#��<#��<#��<#�<#�z<#�<#�J<#��<#�d<#�y<#�b<#�<#ߨ<#ߙ<#ߘ<#ߜ<#�<#�h<#�<#�<#�<#�V<#�<#�<#��<#�q<#�/<#�o<#ޡ<#��<#�s<#�<#��<#��<#ߐ<#�m<#�4<#�o<#ߏ<#ߍ<#�j<#�3<#ߚ<#�<#��<#��<#��<#�<#�U<#�!<#�Z<#�<#�A<#�<#��<#�B<#ٓ<#ߑ<#��<#��<#�"<#�<#�`<#�y<#�<$3�<#��<#�<#�<#�<#��<#�z<#�<#��<#�R<#��<#��<#߆<#��<#�Z<#�Z<#�b<#�r<$	�<$(�<#�k<#�<#��<$}<#��<#�<#�S<#��<#��<#��<#ߥ<#ߣ<#�<#�M<#�<#��<#��<#ߞ<#��<#�<#�<#�Y<#��<#��<$`<$�<#�<#ٰ<#�^<#��<#�<#�<#�]<#�<#��<#�<#��<$�<#�<#��<#�$<#� <#܆<#��<#�<$)�<#�<#�<#��<#�&<#�<#߷<#�F<#��<#ߐ<#��<#ق<#�Q<#�<#�-<#�<#��<#��<#ߌ<#ٱ<#�<#�N<#�G<#��<#�l<$<#�<#�<#��<$<#��<#�1<#��<#�]<#��<#�<#�w<#��<#�Q<#��<#��<#�k<#�<#�Y<#�<#��<#��<#��<#��<#�<#��<#�<#�L<#�<#�<#�;<#��<#�u<#�U<#�?<#��<#�r<#�!<#��<#��<#��<#��<#��<#۹<#�g<#��<#��<#۷<#�g<#ۿ<#��<#��<#۸<#٢<#ۿ<#��;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oPRES            TEMP            PSAL            PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP, where dP is SURFACE PRESSURE (minus 5 dbar for Apf-5,7,8) from next cycle.                                                                                                                                                           TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt(sw_cndr(PSAL,TEMP,PRES),TEMP,PRES_ADJUSTED)                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             PSAL_ADJUSTED = PSAL - dS                                                                                                                                                                                                                                        dP=-0.26                                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                              ds=0                                                                                                                                                                                                                                                           Pressures adjusted using despiked reported SURFACE PRESSURE (1 dBar threshold) from the subsequent profile. The quoted error is 2.4 dBar.                                                                                                                       The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   Salinity adjusted for effects of pressure adjustment. The quoted error is max(0.01, 1xOW uncertainty) in PSS-78.                                                                                                                                                N/A                                                                                                                                                                                                                                                             N/A                                                                                                                                                                                                                                                             OWC(2018v01). Mapping scales LON 3.2/0.8 LAT 1/0.5 MAPSCALE_PHI 0.1/0.02. MAPSCALE_AGE 0.69/10. MAP_P_DELTA 50. Compared with CTD2019v01 and ARGO2020v01 ref. data.                                                                                             202102231454382021022411435520210223145438202102231454382021022411435520210224114355BO  BO  BO  BO  BO  BO  BO  ARGQARGQARGQARGQARGQARSQARSQRTSPPREXRTQCRTQCSCUTnullOW  1.0 2.0 2.0 2.0 2.0 null0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                 20190929151650201909291516502019092915165420190929151701202102231444402021022314543820210224114355  CV  CV  QCP$QCP$QCP$IP  IP                                                                                                                  G�O�G�O�@9��@9��@9��G�O�G�O�G�O�G�O�Dx@ Dx@ Dx@ G�O�G�O�G� G� G� G� G� G� G�                                 6389758         6389758         131072                                          