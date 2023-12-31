CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       �2022-01-26T15:40:01Z creation; 2022-07-24T09:13:06Z last update (coriolis COQC software); 2023-02-15T11:51:39ZZ DMQC performed on CORE variables   
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.41   Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_050k      comment_dmqc_operator         iPRIMARY | http://orcid.org/0000-0003-2516-6106 | Jan Even Øie Nilsen, Institute of Marine Research (IMR)         @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    8d   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    8t   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    8x   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    8|   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    8�   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    8�   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    8�   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  8�   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  8�   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  94   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        9d   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    9h   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    9l   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     9p   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    9�   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    9�   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     9�   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     9�   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     9�   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    9�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        ?F�l�l   
_FillValue        A.�~       axis      T      comment_on_resolution         �JULD resolution is 1 minute, except when JULD = JULD_LOCATION or when JULD = JULD_FIRST_MESSAGE (TRAJ file variable); in that case, JULD resolution is 1 second         9�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    :   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            :   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           :   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           :   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    :    POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    :$   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    :,   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    :0   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    :4   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    :8   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        ;8   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        `  ;<   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        `  C4   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  I�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     `  K,   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  Q�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  W�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  Y�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  _�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  a|   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  g�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  n<   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  o�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                 �  v4   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     `  w�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    ��   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    ��   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    ��   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         ��   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �    HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ~,   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ~\   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    �\   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    �\   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  �\Argo profile    3.1 1.2 19500101000000  20220126154001  20230215115139  6903558 Norway-BGC-Argo                                                 Kjell Arne Mork                                                 PRES            TEMP            PSAL               D   IF                                  2C  D   ARVOR_D                         AD2700-18NO003                  5608A13                         838 @���q�1   @����8�@R�bM���O�;1   GPS     A   A   A   Primary sampling: averaged [10 sec sampling, 1 dbar average from surface to 10 dbar; 10 sec sampling, 2 dbar average from 10 dbar to 1000 dbar; 10 sec sampling, 5 dbar average from 1000 dbar to 1000 dbar]                                                       AK33AfffA���A���A�  A�33A�ffA�ffA���A�  B  B
  B��B��B!��B)33B2  B9��BA��BI��BQ33BY33B`��Bh��BpffBy��B�ffB�33B�  B���B���B�ffB���B�  B���B�ffB���B�33B�  B���B���B���B�33B���Bș�B�33B���B�  B�ffB���B�33B���B���B�  B�B���B���B���C � CffCffCffCL�C
ffC� C� CffCffC33CL�C� C��C� CL�C ffC"L�C$ffC&L�C(33C*��C,� C.�C0L�C233C433C6ffC8ffC:� C<� C>� C@ffCBL�CD33CFffCH33CJffCL33CN33CPL�CRffCTffCVffCX� CZffC\L�C^ffC`� Cb� CdffCfL�Ch33Cj33ClL�CnL�CpffCr33CtL�Cv33CxL�Cz33C|ffC~L�C��C�@ C��C��C�33C�33C�33C�33C�33C�&fC�&fC�@ C�33C�@ C�&fC�&fC�&fC��C�33C�&fC��C��C�33C��C�&fC��C�@ C�L�C��C�33C�&fC�@ C�&fC��C�&fC�&fC�&fC�&fC�@ C�@ C�L�C�@ C�@ C�&fC�&fC�33C�33C�33C�&fC�&fC�&fC�&fC��C�33C�&fC�33C�33C�33C�33C�&fC�L�C�33C�&fC�33C�@ C�33C�&fC�&fC�33C�33C�@ C�&fC�@ C�@ C�33C�@ C��C�33C�&fC�&fC�@ C�33C�&fC�@ C��C�33C��C�33C�@ C�33C�33C�&fC��C�&fC�&fC�@ C��C��C�33C��C�&fC�&fC�33C�&fC��C�33C�&fC�33C�33C�&fC�&fC�&fC�&fC�33C�L�C�L�C�33C�&fC�@ C�L�C�@ C�33C�&fC�@ C�L�C�33C��C�&fD 3D �3DfD��D3D�3D  D��D3D� D  D� D3D� D  D� D3D�3D	�D	��D
�D
��D3D� D�D�3D�D� D�D� D�D�3D�D��D�D��D�D��D�D� D  D��D�D��D3D��D  D��D  D�3D�D��D  D� D  D� D�D��D�D��D�D�3D�D��D &fD ��D!3D!� D"&fD"� D#3D#��D$�D$��D%�D%��D&�D&�3D'  D'��D(fD(��D)�D)� D*  D*�3D+�D+�3D,�D,��D-&fD-�3D.fD.� D/�D/�3D0�D0��D1  D1��D2  D2�3D3�D3��D4�D4��D5  D5� D6  D6��D73D7��D83D8�3D9  D9��D:3D:�3D;�D;��D<�D<��D=  D=��D>3D>�3D?�D?�3D@3D@�3DA  DA� DB&fDB��DC�DC� DD�DD� DE  DE��DF�DF�fDG�DG�fDH�DH�3DI3DI� DJ3DJ��DK�DK��DL�DL�3DM�DM��DN3DNs3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111AK33AfffA���A���A�  A�33A�ffA�ffA���A�  B  B
  B��B��B!��B)33B2  B9��BA��BI��BQ33BY33B`��Bh��BpffBy��B�ffB�33B�  B���B���B�ffB���B�  B���B�ffB���B�33B�  B���B���B���B�33B���Bș�B�33B���B�  B�ffB���B�33B���B���B�  B�B���B���B���C � CffCffCffCL�C
ffC� C� CffCffC33CL�C� C��C� CL�C ffC"L�C$ffC&L�C(33C*��C,� C.�C0L�C233C433C6ffC8ffC:� C<� C>� C@ffCBL�CD33CFffCH33CJffCL33CN33CPL�CRffCTffCVffCX� CZffC\L�C^ffC`� Cb� CdffCfL�Ch33Cj33ClL�CnL�CpffCr33CtL�Cv33CxL�Cz33C|ffC~L�C��C�@ C��C��C�33C�33C�33C�33C�33C�&fC�&fC�@ C�33C�@ C�&fC�&fC�&fC��C�33C�&fC��C��C�33C��C�&fC��C�@ C�L�C��C�33C�&fC�@ C�&fC��C�&fC�&fC�&fC�&fC�@ C�@ C�L�C�@ C�@ C�&fC�&fC�33C�33C�33C�&fC�&fC�&fC�&fC��C�33C�&fC�33C�33C�33C�33C�&fC�L�C�33C�&fC�33C�@ C�33C�&fC�&fC�33C�33C�@ C�&fC�@ C�@ C�33C�@ C��C�33C�&fC�&fC�@ C�33C�&fC�@ C��C�33C��C�33C�@ C�33C�33C�&fC��C�&fC�&fC�@ C��C��C�33C��C�&fC�&fC�33C�&fC��C�33C�&fC�33C�33C�&fC�&fC�&fC�&fC�33C�L�C�L�C�33C�&fC�@ C�L�C�@ C�33C�&fC�@ C�L�C�33C��C�&fD 3D �3DfD��D3D�3D  D��D3D� D  D� D3D� D  D� D3D�3D	�D	��D
�D
��D3D� D�D�3D�D� D�D� D�D�3D�D��D�D��D�D��D�D� D  D��D�D��D3D��D  D��D  D�3D�D��D  D� D  D� D�D��D�D��D�D�3D�D��D &fD ��D!3D!� D"&fD"� D#3D#��D$�D$��D%�D%��D&�D&�3D'  D'��D(fD(��D)�D)� D*  D*�3D+�D+�3D,�D,��D-&fD-�3D.fD.� D/�D/�3D0�D0��D1  D1��D2  D2�3D3�D3��D4�D4��D5  D5� D6  D6��D73D7��D83D8�3D9  D9��D:3D:�3D;�D;��D<�D<��D=  D=��D>3D>�3D?�D?�3D@3D@�3DA  DA� DB&fDB��DC�DC� DD�DD� DE  DE��DF�DF�fDG�DG�fDH�DH�3DI3DI� DJ3DJ��DK�DK��DL�DL�3DM�DM��DN3DNs3111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��?%�T?$�/?#S�?"�\?"�\?"�\?"M�?"J?!G�? �?��?�-?��??}>�33>��!>e`B>F��>8Q�>1&�>/�>9X>I�^>V>_;d>]/>M��>I�^>F��>E��>B�\>D��>:^5>7K�>0 �>,1>,1>'�>#�
>(��>hs>�>\)>��>�w>$�/>%�T>%�T>2->9X>;dZ>?|�>B�\>C��>B�\>?|�>:^5>49X>,1> Ĝ>\)>�>�>   =�=�=�F=�F=��==�h=�h=�h=�=�x�=�l�=�`B=�`B=�`B=�l�=�h==��=��=�h=�`B=�;d=�"�=�
==���=��=���=���=���=���=���=�
==�;d=�S�=�`B=�l�=�F=��=��#=��=�F=�F=�F=�F=�F=�=�=���=�F=�=��==��=��=��=�h=�h=�x�=�`B=�"�=��`=��`=ě�=\=��=�j=�E�=�-=� �=�{=�{=�{=�-=�9X=� �=���=��
=��w=��=�\)=�C�=��=u=aG�=]/=Y�=P�`=D��='�=�P=t�=�P=\)=o<���<�1<�o<o;�o:�o    �o��o���
��`B�o�49X�#�
�#�
�D���u��C���t����
��1��j�ě�������`B��h���o�C��\)�t���P����w��w����w��w�'0 Ž0 Ž49X�<j�D���H�9�P�`�P�`�Y��Y��T���Y��e`B�e`B�m�h�q���u�}󶽃o��+��7L��C���\)��hs��t����P�������㽝�-�����������
���T���T���T����罩�罩�罩�罩�罬1��1��{��{��{��{�� Ž�-��-��-��-��9X��9X��9X��E���Q콸Q콸Q콺^5��^5��Q콸Q콸Q콺^5��^5��^5��j��^5��j��^5��^5��^5��Q콸Q콸Q콸Q콶E���E���E���Q콸Q콶E���E���Q콸Q콸Q콸Q콸Q콺^5��Q콺^5��Q콺^5��j��^5��j��^5��^5��vɽ�vɽ�j��vɽ������\�ě��Ƨ�ȴ9�ȴ9�������ͽ��ͽ�����������`����������
=��
=�����"ѽ�"ѽ�/��;d��;d��G���G���S���`B��l���xս�xս�xս���h��h��������ٽ��m���   �   �   �J�o�����+�1'�1'�
=q�C��C��I��I��O߾V�V�\)�bN�hs�t��z����+��P��P��u�������������R��w��w� Ĝ�!���"��"��#�
�$�/�$�/�$�/�%�T�&�y�&�y�')��)��+�+�,1�-V�-V�/��/��0 ž0 ž1&�0 ž2-�333�333�49X111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111?%�T?$�/?#S�?"�\?"�\?"�\?"M�?"J?!G�? �?��?�-?��??}>�33>��!>e`B>F��>8Q�>1&�>/�>9X>I�^>V>_;d>]/>M��>I�^>F��>E��>B�\>D��>:^5>7K�>0 �>,1>,1>'�>#�
>(��>hs>�>\)>��>�w>$�/>%�T>%�T>2->9X>;dZ>?|�>B�\>C��>B�\>?|�>:^5>49X>,1> Ĝ>\)>�>�>   =�=�=�F=�F=��==�h=�h=�h=�=�x�=�l�=�`B=�`B=�`B=�l�=�h==��=��=�h=�`B=�;d=�"�=�
==���=��=���=���=���=���=���=�
==�;d=�S�=�`B=�l�=�F=��=��#=��=�F=�F=�F=�F=�F=�=�=���=�F=�=��==��=��=��=�h=�h=�x�=�`B=�"�=��`=��`=ě�=\=��=�j=�E�=�-=� �=�{=�{=�{=�-=�9X=� �=���=��
=��w=��=�\)=�C�=��=u=aG�=]/=Y�=P�`=D��='�=�P=t�=�P=\)=o<���<�1<�o<o;�o:�o    �o��o���
��`B�o�49X�#�
�#�
�D���u��C���t����
��1��j�ě�������`B��h���o�C��\)�t���P����w��w����w��w�'0 Ž0 Ž49X�<j�D���H�9�P�`�P�`�Y��Y��T���Y��e`B�e`B�m�h�q���u�}󶽃o��+��7L��C���\)��hs��t����P�������㽝�-�����������
���T���T���T����罩�罩�罩�罩�罬1��1��{��{��{��{�� Ž�-��-��-��-��9X��9X��9X��E���Q콸Q콸Q콺^5��^5��Q콸Q콸Q콺^5��^5��^5��j��^5��j��^5��^5��^5��Q콸Q콸Q콸Q콶E���E���E���Q콸Q콶E���E���Q콸Q콸Q콸Q콸Q콺^5��Q콺^5��Q콺^5��j��^5��j��^5��^5��vɽ�vɽ�j��vɽ������\�ě��Ƨ�ȴ9�ȴ9�������ͽ��ͽ�����������`����������
=��
=�����"ѽ�"ѽ�/��;d��;d��G���G���S���`B��l���xս�xս�xս���h��h��������ٽ��m���   �   �   �J�o�����+�1'�1'�
=q�C��C��I��I��O߾V�V�\)�bN�hs�t��z����+��P��P��u�������������R��w��w� Ĝ�!���"��"��#�
�$�/�$�/�$�/�%�T�&�y�&�y�')��)��+�+�,1�-V�-V�/��/��0 ž0 ž1&�0 ž2-�333�333�49X111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
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
��B�B;dBp�Bq�By�B� B�B�B�7B�JB�VB�\B�bB�hB�hB�uB�{B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�B�B�B�B�B�B�B�!B�!B�B�!B�!B�B�B�B�!B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B�B��B�B�B��B�B��B�B��B�B�B�B�B��B��B�B�B�B�B�B�B�B�B�B�B��B�B�B��B�B�B�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111B
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
��B�B;dBp�Bq�By�B� B�B�B�7B�JB�VB�\B�bB�hB�hB�uB�{B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�B�B�B�B�B�B�B�B�!B�!B�B�!B�!B�B�B�B�!B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B��B��B�B��B�B�B��B�B��B�B��B�B�B�B�B��B��B�B�B�B�B�B�B�B�B�B�B��B�B�B��B�B�B�B��B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES.                                                                                                                                                                                                                                           TEMP_ADJUSTED = TEMP.                                                                                                                                                                                                                                           PSAL_ADJUSTED = PSAL.                                                                                                                                                                                                                                           none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            The quoted error is manufacturer specified accuracy.                                                                                                                                                                                                            The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity offset or drift detected. The quoted error is max[0.01, statistical uncertainty] in PSS-78.                                                                                                                                             202302131647062023021316470620230213164706  IF  ARFMCODA046k                                                                20220126154001                      G�O�G�O�G�O�                IF  ARGQCOQC5.6                                                                 20220126154121                      G�O�G�O�G�O�                IF  ARGQCOQC5.6                                                                 20220126154121                      G�O�G�O�G�O�                IF  ARFMCODA050d                                                                20220614090127                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220614090239                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220614090239                      G�O�G�O�G�O�                IF  ARFMCODA050f                                                                20220624085532                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220624085643                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220624085643                      G�O�G�O�G�O�                IF  ARFMCODA050h                                                                20220704090304                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220704090416                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220704090416                      G�O�G�O�G�O�                IF  ARFMCODA050k                                                                20220714090326                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220714090430                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220714090430                      G�O�G�O�G�O�                IF  ARFMCODA050k                                                                20220724091156                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220724091306  QCP$                G�O�G�O�G�O�000000000288F35EIF  ARGQCOQC5.8                                                                 20220724091306  QCF$                G�O�G�O�G�O�0000000000000000