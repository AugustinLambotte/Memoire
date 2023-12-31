CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   r   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       X2020-05-11T14:45:28Z creation; 2020-05-11T14:46:11Z last update (coriolis COQC software)   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_034d      comment_dmqc_operator         DPRIMARY | https://orcid.org/0000-0002-3512-2070 | Saout-Grit, Glazeo      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8,   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8<   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8L   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8T   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  8�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    90   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    94   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     98   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9X   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9x   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        :�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  :�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  <�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  =   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  >�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  ?T   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  A   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  B�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  CX   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  E    TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  E�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  G\   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  I$   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  I�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  t  K`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     �  K�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    V�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    V�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    W    HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    W   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  W   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    WH   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    WX   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    W\   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         Wl   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         Wp   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        Wt   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    Wx   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  M�   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    M�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    P�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    S�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  V�Argo profile    3.1 1.2 19500101000000  20200511144528  20211004092454  6902726 NARVAL                                                          Camille DAUBORD                                                 PRES            TEMP            PSAL               A   IF                                  2C  D   ARVOR                           AI2600-16FR309                  5900A04                         844 @�/l}'�}1   @�/l}'�}@Sfq]���/��#��8   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 25 dbar average from 2000 dbar to 500 dbar; 10 sec sampling, 10 dbar average from 500 dbar to 50 dbar; 10 sec sampling, 1 dbar average from 50 dbar to 2.5 dbar]                                                      A���A���A���A���A�  A�  A���A���A���B   BffBffBffBffBffBffBffB   B$  B(ffB,  B0ffB3��B733B<  B@  BDffBH  B]��B���B�33B�ffB���B���B晚B�  C� C��CffC%  C/�C933CCL�CMffCW33Ca�CkffCu�fC��C��3C�� C�ٚC���C���C��fC��fC��3C�� C��fC��fC��fC�� Cř�CʦfC���C�� C٦fCޙ�C�� C�fC홚C�3C�ٚD FfD� D��D  D9�Ds3D%�3D+��D2&fD8ffD>�fDD��DK9�DQffDW��D^fDd  Djl�Dp� Dw�D}9�D��fD�� D���D��D�P D�Y�D�|�D���D���D��fD���D��D�<�D�Y�D�s3D���D���D�� D�  D��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A���A���A���A���A�  A�  A���A���A���B   BffBffBffBffBffBffBffB   B$  B(ffB,  B0ffB3��B733B<  B@  BDffBH  B]��B���B�33B�ffB���B���B晚B�  C� C��CffC%  C/�C933CCL�CMffCW33Ca�CkffCu�fC��C��3C�� C�ٚC���C���C��fC��fC��3C�� C��fC��fC��fC�� Cř�CʦfC���C�� C٦fCޙ�C�� C�fC홚C�3C�ٚD FfD� D��D  D9�Ds3D%�3D+��D2&fD8ffD>�fDD��DK9�DQffDW��D^fDd  Djl�Dp� Dw�D}9�D��fD�� D���D��D�P D�Y�D�|�D���D���D��fD���D��D�<�D�Y�D�s3D���D���D�� D�  D��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���j���e��e`B�bMӿe��hr��^�ۿq녿�񪿕�˿��9��Q쿩�����h��r���dZ���ÿ�����𿦇+��ff������
=���+���T����?}��?}���\�Vȴ�7
==#�
�<j>�ȴ?&��?��@)G�@Y��@Vff@T�D@Dj@Co@&ȴ@(�9@$9X@��@ff@�7?��H?��?�X?��?��?�X?�`B?�p�?�b?��?�`B?v?kC�?`  ?a�7?Z^5?Hr�?7�P?7��?-��?$�/?��?�??�>��>��
>��`>��>�?}>���>�t�>�+>p��>Z�>S��>&�y>J=ȴ9=���=49X<u��`B��/�@���+��-��;d�o��P�2-�B�\�V�o�����7��1'��I����`��zᾙ����w��1��KǾ\�ɺ^��bN111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  �j���e��e`B�bMӿe��hr��^�ۿq녿�񪿕�˿��9��Q쿩�����h��r���dZ���ÿ�����𿦇+��ff������
=���+���T����?}��?}���\�Vȴ�7
==#�
�<j>�ȴ?&��?��@)G�@Y��@Vff@T�D@Dj@Co@&ȴ@(�9@$9X@��@ff@�7?��H?��?�X?��?��?�X?�`B?�p�?�b?��?�`B?v?kC�?`  ?a�7?Z^5?Hr�?7�P?7��?-��?$�/?��?�??�>��>��
>��`>��>�?}>���>�t�>�+>p��>Z�>S��>&�y>J=ȴ9=���=49X<u��`B��/�@���+��-��;d�o��P�2-�B�\�V�o�����7��1'��I����`��zᾙ����w��1��KǾ\�ɺ^��bN111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oA�dZA���A�bNA���A�I�A��TA�n�A�ƨA�M�A��;A���A��`A�|�A��A��FA�I�A���A��B 5?B �B ��B ȴB�By�B�yB49By�B�DB��B�B��B�9B	B	�VB	�TB
x�BC�B��B�B�B�yB��BŢB�BB��BB�B�B�B��B�oB�+B}�B{�B{�By�Bx�Bz�B}�Bx�Bx�Bw�Bz�B{�Bw�Bu�Bv�Bw�Bu�Bs�Bt�Bv�Bv�Bv�Bv�Bx�By�Bz�B|�B}�B~�B� B�B�%B�+B�+B�7B�DB�JB�PB�VB�bB�hB�oB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  A�f�A� yA�d�A��)A�L;A���A�q"A��1A�PaA���A��A���A�}A�!KA���A�LRA��\A���B 6�B �iB �B �B B{/B��B5�B{3B��B��B�B�-B��B	oB	��B	�B
zCBEB�?BێB�	B��B�,B�B�B�kB��BڇBـB��B�B��B��BcB}VB}WB{LBzDB|RBbBzCBzEBy>B|PB}WBy=Bw2Bx5By?Bw/Bu"Bv*Bx8Bx7Bx6Bx6BzDB{JB|PB~]BcB�gB�nB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�"B�(B�/B�-B�5B�8B�8B�AB�FB�KB�QB�VB�^111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL  + Delta_S, where Delta_S is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r                                                                                                                     none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            r= 1 (+/- 9e-06) , vertically averaged dS =0.001364 (+/- 0.01)                                                                                                                                                                                                  No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              No adjustement was necessary -Calibration error is manufacturer specified accuracy                                                                                                                                                                              Salinity drift or offset detected - OWC fit is adopted. Error = maximum [statistical uncertainty, 0.01]. OWC Method, 2.0,  -CTD2021V01 & ARGO2020V03 -                                                                                                          202005120935392021100409245420211004092454  IF  ARFMCODA034d                                                                20200511144528                      G�O�G�O�G�O�                IF  ARGQCOQC4.5                                                                 20200511144611  QCP$                G�O�G�O�G�O�000000000208F35EIF  ARGQCOQC4.5                                                                 20200511144611  QCF$                G�O�G�O�G�O�0000000000000000IF  ARSQOW  1.1 CTD2018V01 & ARGO2018V01                                        20200512093539  IP  PSAL            A���D��G�O�                IF  ARSQOW  2.0 CTD2021V01 & ARGO2020V03                                        20211004092454  IP  PSAL            A���D��G�O�                