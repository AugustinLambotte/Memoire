CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   A   	N_HISTORY          N_CALIB             title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2011-04-04T15:47:44Z creation; 2018-08-16T12:50:35Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    6�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    6�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  �  7   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  �  7�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  `  8   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        8h   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    8p   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    8t   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                  @  8x   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                  @  8�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                  @  9   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                  @  9D   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    9�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            9�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           9�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           9�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    9�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    9�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ;�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          ;�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  >    PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z          >�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  @�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���       A   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       C   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  E    TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       E�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  G�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       H0   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       J8   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  L@   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       L�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  �  N�   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o       OP   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    d   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    d   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    d   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    d$   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  �  d,   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    d�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    d�   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                     d�   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         d�   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         d�   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        e    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                     e   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  `  QX   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Q�   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    W�   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ]�   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  T  c�  c�Argo profile    3.1 1.2 19500101000000  20110404154744  20180816125035  6900842 6900842 BSH                                                             BSH                                                             Holger GIESE                                                    Holger GIESE                                                    PRES            TEMP            PSAL            PRES            TEMP                                  AA  IFIF18082685                                                        2C  2B  DA  APEX                            APEX                            5037                            5037                            061609                          061609                          846 846 @��.��A�@��.��A�11  @��/�� @��/�� @P�9XbN@P�9XbN�)9XbM��)9XbM�11  ARGOS   ARGOS   AD  AD  A   Primary sampling: discrete []                                                                                                                                                                                                                                   Near-surface sampling: discrete, unpumped []                                                                                                                                                                                                                          A�33B   BI��Bq��B���B���B���B���B�ffB�  C33CffC� C �C*�C4� C>ffCG�fCR� C\ffCe��Co��Cy��C�� C��3C��3C�L�C�@ C��C�&fC�ٚC�&fC�&fC��C�33C��3C�&fC�  C��C��3C��3C�&fC��C�&fC�33C��3C��fC��D9�D	l�D� D3D"��D/fD;l�DG�3DT��D`��Dml�Dz�D�C3D��fD�ɚD���D�Ff=��;L��                �L��>L��?fff?�33@��@9��@���AffAi��A�33G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111114444443111111111                                                   A���B ��BJfgBrfgB�33B�33B�33B�33B���B�ffCffC��C�3C L�C*L�C4�3C>��CH�CR�3C\��Ce��Co��Cz  C���C��C��C�fgC�Y�C�&gC�@ C��4C�@ C�@ C�34C�L�C��C�@ C��C�&gC��C��C�@ C�&gC�@ C�L�C��C�  C�34DFgD	y�D��D  D"�gD/3D;y�DH  DT��DagDmy�Dz�D�I�D���D�� D�3D�L�    �������ͽ��ͽ��ͽ��;���=���?L��?�ff@ff@333@�ffA��Ah  A�ffG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111114444443331111111                                                   @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O����+�t��\)�+�C��49X�Y���9X��Q콾vɽ�Q��w��t�=C�=���=��=#�
<���<D���D���u������j�����1'�hs����"��)��,1�333�;dZ�D���N��T���Xb�_;d�ixվm�h�p�׾s�F�x���z�H�}󶾃o�����7L��\)���u��S���^5����V��^5��T��Ϳhs�
=��m�#���$���*=q�-O�=m�h=aG�=e`B=]/=P�`=T��=0 �=�%=��=�%=��=��=q���o�#�
�uG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111113333333111111111                                                   ���+�t��\)�+�C��49X�Y���9X��Q콾vɽ�Q��w��t�=C�=���=��=#�
<���<D���D���u������j�����1'�hs����"��)��,1�333�;dZ�D���N��T���Xb�_;d�ixվm�h�p�׾s�F�x���z�H�}󶾃o�����7L��\)���u��S���^5����V��^5��T��Ϳhs�
=��m�#���$���*=q�-O�=m�h=aG�=e`B=]/=P�`=T��=0 �=�%=��=�%=��=��=q���o�#�
�uG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�111111111111111111111111111111111111111111111111111111111111111113333333331111111                                                   ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B2-B2-B2-B33B33B33B2-B2-B6FB7LB8RB9XB@�BG�BXBp�Bw�B{�B|�B}�B� B�B�B�B�+B�1B�1B�7B�=B�DB�DB�JB�PB�PB�PB�VB�VB�VB�\B�bB�hB�oB�oB�uB�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111                                                                   B2B2B2B3B3B3B2B2B6$B7-B83B99B@cBG�BW�Bp�Bw�B{�B|�B}�B�B��B��B� B�B�B�B�B�B�&B�&B�)B�0B�0B�1B�6B�6B�6B�<B�CB�FB�PB�PB�UB�\B�[B�gB�lB�tB�yB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�11111111111111111111111111111111111111111111111111111111111111111                                                                   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�PRES            TEMP            PSAL            PRES            TEMP                            PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     PRES_ADJUSTED = PRES - Surface Pressure                                                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Surface pressure = -0.2 dbar                                                                                                                                                                                                                                    none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Surface Pressure = 0.1 dbar                                                                                                                                                                                                                                     Not applicable                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.                                                                                                              Pressure adjusted in real time by using pressure offset at the sea surface                                                                                                                                                                                      No adjustment performed (values duplicated)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     2018081612503520180816125035201808161250352018041312080720180413120807              IF  IF  ARGQARFMCOARCODA1.0 019c                                                                                                                                2011040411443520180413120807QCP$                                    G�O�G�O�G�O�G�O�G�O�G�O�DEBFC                           IF  IF  ARGQARGQCOARCOQC1.0 3.2                                                                                                                                 2011040411443520180413120854QCF$QCP$                                G�O�G�O�G�O�G�O�G�O�G�O�02000           000000000008FB7EIF  IF  CORTARGQCOOACOQC5.2 3.2 RTQCGL01                                                                                                                        2011040505572220180413120854QCP$QCF$PSAL                            G�O�G�O�G�O�G�O�G�O�G�O�                0000000000000000IF      CORT    COOA    5.2     RTQCGL01                                                                                                                        20110406044547              QCF$    PSAL                            G�O�G�O�G�O�G�O�G�O�G�O�4                               IF      CORT    COOA    5.2     RTQCGL01                                                                                                                        20110406044227              QCF$    TEMP                            G�O�G�O�G�O�G�O�G�O�G�O�4                               IF      CORT    COOA    5.2     RTQCGL01                                                                                                                        20110405055116              QCP$    TEMP                            G�O�G�O�G�O�G�O�G�O�G�O�                                IF              SCOO    1.3                                                                                                                                     20110404172256              QC                                      G�O�G�O�G�O�G�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD reference database, Version: CTD_for_DMQC_2010V2, 04/11                                                                20110712114307              IP      PSAL                            A�33G�O�D�FfG�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2012V01 + ARGO climatology                                                                 20120626104111              IP      PSAL                            A�33G�O�D�FfG�O�G�O�G�O�                                IF  IF          COFCCOFC3.0 3.0                                                                                                                                 2018061115480520180611154805                                        G�O�G�O�G�O�G�O�G�O�G�O�                                GE      ARSQ    OW      1.0     ARGO CTD ref. database: CTD_for_DMQC_2017V01 + ARGO climatology                                                                 20180816125035              IP      PSAL                            A�33G�O�D�FfG�O�G�O�G�O�                                