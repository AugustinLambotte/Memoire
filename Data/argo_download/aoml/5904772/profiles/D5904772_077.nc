CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:31Z creation      
references        (http://www.argodatamgt.org/Documentation   comment       	free text      user_manual_version       3.2    Conventions       Argo-3.2 CF-1.6    featureType       trajectoryProfile      comment_dmqc_operator         (Matthew Alkire, University of Washington      @   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    6�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    6�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    6�   REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    6�   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    7   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    7   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    7,   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  74   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  7t   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  0  7�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        7�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    7�   DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    7�   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     7�   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    8   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    8   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     8   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     88   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     8X   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    8x   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~       axis      T           8|   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    8�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ      
_FillValue        A.�~            8�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           8�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           8�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    8�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    8�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    8�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        9�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    9�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    9�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        �  9�   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    A�   PRES_ADJUSTED            
      	   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  C�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    K�   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     �  M�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  U�   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    ]�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  _�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    g�   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  i�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  q�   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    y�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  {�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                    �|   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  �|   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  �t   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    ��   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  ��   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �    HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �0   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �4   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �D   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �H   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �L   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004031  20190405134403  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               MA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�k-�e�1   @�k.@ys�@O#�E����B���O�;1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    MA   B   B   @�ff@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� DmfDm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� DtٚDy��D�	�D�L�D�� D�� D���D�L�D���D���D�3D�9�D��fD��3D�3D�C3D�C3D��fD�	�D�6fD�i�D�� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@�G�A��A,��AL��Al��A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�A�Q�B(�B(�B(�B(�B#(�B+(�B3(�B;(�BC(�BK(�BS(�B[(�Bc(�Bk(�Bs(�B{(�B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{B��{BŔ{Bɔ{B͔{Bє{BՔ{Bٔ{Bݔ{B�{B�{B�{B�{B�{B��{B��{B��{C �=C�=C�=C�=C�=C
�=C�=C�=C�=C�=C�=C�=C�=C�=C�=C�=C �=C"�=C$�=C&�=C(�=C*�=C,�=C.�=C0�=C2�=C4�=C6�=C8�=C:�=C<�=C>�=C@�=CB�=CD�=CF�=CH�=CJ�=CL�=CN�=CP�=CR�=CT�=CV�=CX�=CZ�=C\�=C^�=C`�=Cb�=Cd�=Cf�=Ch�=Cj�=Cl�=Cn�=Cp�=Cr�=Ct�=Cv�=Cx�=Cz�=C|�=C~�=C�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�q�C�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eC�eD 2�D ��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D	2�D	��D
2�D
��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D2�D��D 2�D ��D!2�D!��D"2�D"��D#2�D#��D$2�D$��D%2�D%��D&2�D&��D'2�D'��D(2�D(��D)2�D)��D*2�D*��D+2�D+��D,2�D,��D-2�D-��D.2�D.��D/2�D/��D02�D0��D12�D1��D22�D2��D32�D3��D42�D4��D52�D5��D62�D6��D72�D7��D82�D8��D92�D9��D:2�D:��D;2�D;��D<2�D<��D=2�D=��D>2�D>��D?2�D?��D@2�D@��DA2�DA��DB2�DB��DC2�DC��DD2�DD��DE2�DE��DF2�DF��DG2�DG��DH2�DH��DI2�DI��DJ2�DJ��DK2�DK��DL2�DL��DM2�DM��DN2�DN��DO2�DO��DP2�DP��DQ2�DQ��DR2�DR��DS2�DS��DT2�DT��DU2�DU��DV2�DV��DW2�DW��DX2�DX��DY2�DY��DZ2�DZ��D[2�D[��D\2�D\��D]2�D]��D^2�D^��D_2�D_��D`2�D`��Da2�Da��Db2�Db��Dc2�Dc��Dd2�Dd��De2�De��Df2�Df��Dg2�Dg��Dh2�Dh��Di2�Di��Dj2�Dj��Dk2�Dk��Dl2�Dl��Dm8�Dm��Dn2�Dn��Do2�Do��Dp2�Dp��Dq2�Dq��Dr2�Dr��Ds2�Ds��Dt2�Dt��Du)Dy�\D�"�D�fD��HD��HD��D�fD���D��D�{D�R�D���D��{D�{D�\{D�\{D��D�"�D�O�D��D��H111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�S�@�O�@۾w@ڟ�@��@�V@��#@� �@�Z@�V@��@���@���@�
=@��T@��-@���@�+@�O�@�A�@�
=@�^5@���@���@���@��@�7L@��m@���@��@�V@��@���@��@��+@��T@�x�@�/@���@��/@��9@�r�@��D@�9X@���@�t�@�
=@��R@���@�-@�/@���@��@���@���@�^5@��T@�`B@�/@��@�b@�l�@�C�@�33@��H@�M�@�7L@��/@��@�I�@��
@��
@���@��P@�;d@�+@���@���@�=q@���@�O�@���@���@��9@�r�@�  @���@���@�K�@�@���@��@��-@�G�@��u@�|�@�l�@�t�@�l�@�S�@�33@��H@�M�@�=q@�E�@�-@�G�@��9@��D@�dZ@���@�1@� �@�9X@�I�@�1@�1@�  @��;@��
@��P@�;d@��@�ff@��@���@�O�@�O�@�X@�O�@�7L@�&�@���@�z�@�1'@��m@��@�S�@�o@���@�~�@�^5@�V@�V@��@��@�$�@�-@�$�@���@�hs@�hs@�hs@�O�@�V@��D@�r�@�j@�bN@�bN@�bN@�r�@�z�@�z�@��@�z�@�j@�j@�Z@�Z@�A�@�1@��@�@|�@
=@~��@~�+@~ff@~E�@~5?@~$�@~{@~{@}�@}��@}�-@}�@}�@}�@}p�@}O�@}/@}V@|��@|�@|�@|j@|1@{�m@{��@{S�@{33@{@{o@{@z�@z�H@z��@z��@z�!@z��@y��@yx�@yx�@yX@yG�@y7L@x��@x�9@xbN@xQ�@xQ�@xA�@xA�@xQ�@xbN@xbN@xQ�@xb@x  @w��@w�@vv�@u�@u�-@u`B@u�@t(�@s�m@s��@t9X@t(�@t9X@t�D@t��@t��@tI�@s��@s�m@s�
@s�
@s�F@s�@sdZ@sS�@sS�@r��@rn�@rn�@r�@q�#@q��@qx�@p��@p1'@p �@p  @o�@o�;@o��@o��@o�w@o��@o��@o�P@o|�@o�P@o�@o�;@o�@o�;@o�w@ol�@o�P@o�@n�y@n�@nȴ@n��@n�+@nff@nE�@n5?@n5?@n5?@n$�@n{@n@n@n@m��@m�h@mO�@m/@l�@l��@l�/@l�j@lz�@lI�@k��@k�m@kƨ@k��@k��@k��@k��@k��@k��@k��@kdZ@k33@ko@k@j�H@j�H@j�H@j�H@j�@j�H@j��@j��@j�\@j~�@j�\@j~�@j^5@j=q@j�@j�@jJ@i��@i�#@i�#@i�#@i��@i��@i�#@i�#@i��@i��@i�7@ix�@i�7@i�7@ihs@iG�@i7L@i7L@i&�@i�@i�@i%@h��@h�`@h��@hĜ@h�9@h��@h�u@h�@hQ�@h �@h �@h �@h1'@h �@h  @g�;@g�;@g��@g��@g��@g�w@g�@g�@g��@g��@g|�@gl�@gl�@gK�@g�@g
=@g
=@g
=@g
=@g
=@f��@f�y@f�@fȴ@fȴ@f�R@f��@fv�@fff@fV@fE�@f$�@f$�@f{@f{@f{@f@e�@e�@e�T@e�-@e�h@ep�@e`B@e`B@eO�@eO�@eO�@eO�@eO�@eO�@eO�@e?}@e?}@e/@e�@d��@d�@d�D@dj@dI�@dI�@dI�@dZ@dI�@d9X@dZ@dj@dj@dZ@dI�@d9X@c��@c��@c�m@c��@c��@c��@c�F@c��@c��@ct�@cS�@c33@c33@c33@c"�@c@b�H@bn�@b-@a��@a��@aG�@a7L@a%@a�@`��@`��@a%@`��@`bN@`Q�@`r�@`�@`A�@`b@_�@_�;@_�@_��@_�P@]�@[�@[�
@\1@_��@c�m@f@gK�@fȴ@e�@c33@`�9@]�-@Yx�@S�F@P��@M��@I�^@IG�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @�S�@�O�@۾w@ڟ�@��@�V@��#@� �@�Z@�V@��@���@���@�
=@��T@��-@���@�+@�O�@�A�@�
=@�^5@���@���@���@��@�7L@��m@���@��@�V@��@���@��@��+@��T@�x�@�/@���@��/@��9@�r�@��D@�9X@���@�t�@�
=@��R@���@�-@�/@���@��@���@���@�^5@��T@�`B@�/@��@�b@�l�@�C�@�33@��H@�M�@�7L@��/@��@�I�@��
@��
@���@��P@�;d@�+@���@���@�=q@���@�O�@���@���@��9@�r�@�  @���@���@�K�@�@���@��@��-@�G�@��u@�|�@�l�@�t�@�l�@�S�@�33@��H@�M�@�=q@�E�@�-@�G�@��9@��D@�dZ@���@�1@� �@�9X@�I�@�1@�1@�  @��;@��
@��P@�;d@��@�ff@��@���@�O�@�O�@�X@�O�@�7L@�&�@���@�z�@�1'@��m@��@�S�@�o@���@�~�@�^5@�V@�V@��@��@�$�@�-@�$�@���@�hs@�hs@�hs@�O�@�V@��D@�r�@�j@�bN@�bN@�bN@�r�@�z�@�z�@��@�z�@�j@�j@�Z@�Z@�A�@�1@��@�@|�@
=@~��@~�+@~ff@~E�@~5?@~$�@~{@~{@}�@}��@}�-@}�@}�@}�@}p�@}O�@}/@}V@|��@|�@|�@|j@|1@{�m@{��@{S�@{33@{@{o@{@z�@z�H@z��@z��@z�!@z��@y��@yx�@yx�@yX@yG�@y7L@x��@x�9@xbN@xQ�@xQ�@xA�@xA�@xQ�@xbN@xbN@xQ�@xb@x  @w��@w�@vv�@u�@u�-@u`B@u�@t(�@s�m@s��@t9X@t(�@t9X@t�D@t��@t��@tI�@s��@s�m@s�
@s�
@s�F@s�@sdZ@sS�@sS�@r��@rn�@rn�@r�@q�#@q��@qx�@p��@p1'@p �@p  @o�@o�;@o��@o��@o�w@o��@o��@o�P@o|�@o�P@o�@o�;@o�@o�;@o�w@ol�@o�P@o�@n�y@n�@nȴ@n��@n�+@nff@nE�@n5?@n5?@n5?@n$�@n{@n@n@n@m��@m�h@mO�@m/@l�@l��@l�/@l�j@lz�@lI�@k��@k�m@kƨ@k��@k��@k��@k��@k��@k��@k��@kdZ@k33@ko@k@j�H@j�H@j�H@j�H@j�@j�H@j��@j��@j�\@j~�@j�\@j~�@j^5@j=q@j�@j�@jJ@i��@i�#@i�#@i�#@i��@i��@i�#@i�#@i��@i��@i�7@ix�@i�7@i�7@ihs@iG�@i7L@i7L@i&�@i�@i�@i%@h��@h�`@h��@hĜ@h�9@h��@h�u@h�@hQ�@h �@h �@h �@h1'@h �@h  @g�;@g�;@g��@g��@g��@g�w@g�@g�@g��@g��@g|�@gl�@gl�@gK�@g�@g
=@g
=@g
=@g
=@g
=@f��@f�y@f�@fȴ@fȴ@f�R@f��@fv�@fff@fV@fE�@f$�@f$�@f{@f{@f{@f@e�@e�@e�T@e�-@e�h@ep�@e`B@e`B@eO�@eO�@eO�@eO�@eO�@eO�@eO�@e?}@e?}@e/@e�@d��@d�@d�D@dj@dI�@dI�@dI�@dZ@dI�@d9X@dZ@dj@dj@dZ@dI�@d9X@c��@c��@c�m@c��@c��@c��@c�F@c��@c��@ct�@cS�@c33@c33@c33@c"�@c@b�H@bn�@b-@a��@a��@aG�@a7L@a%@a�@`��@`��@a%@`��@`bN@`Q�@`r�@`�@`A�@`b@_�@_�;G�O�@_��@_�P@]�@[�@[�
@\1@_��@c�m@f@gK�@fȴ@e�@c33@`�9@]�-@Yx�@S�F@P��@M��@I�^@IG�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�PB�PB�PB�JB�=B�%B�+B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�VB�VB�VB�VB�\B�VB�VB�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�JB�JB�DB�DB�DB�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�+B�+B�B�B�B�B�bB��B�B�?B�dB�qB�wB�wB�}B�jB�^B�^B�RB�?B�j111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B��B��B��B��B��B��B��B�&B�:B�LB�RB�YB�WB�eB�eB�^B�_B�^B�mB�jB�lB�dB�fB�eB�dB�_B�WB�PB�DB�:B�8B�,B�,B�-B�%B�B�B�$B�,B�0B�1B�:B�IB�QB�\B�\B�\B�dB�iB�nB�oB�pB�qB�uB�vB�wB�vB�zB�zB�|B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�vB�uB�zB�yB�yB�xB�{B�{B�zB��B��B�sB�hB�aB�VB�bB�nB�tB�zB�}B��B��B�B��B��B��B�{B�uB�oB�iB�`B�fB�hB�nB�oB�oB�nB�nB�nB�hB�gB�aB�[B�]B�UB�UB�UB�UB�[B�[B�[B�[B�]B�[B�\B�UB�VB�TB�TB�MB�HB�IB�IB�FB�IB�OB�MB�UB�UB�UB�UB�UB�_B�_B�]B�]B�]B�UB�UB�UB�UB�PB�NB�RB�NB�PB�NB�RB�NB�PB�PB�NB�OB�UB�TB�VB�UB�SB�^B�[B�^B�^B�[B�[B�UB�UB�VB�VB�UB�_B�\B�\B�]B�]B�_B�]B�\B�SB�WB�VB�WB�UB�PB�PB�PB�OB�WB�WB�WB�NB�VB�XB�YB�WB�XB�XB�XB�OB�HB�IB�HB�DB�BB�>B�>B�>B�CB�JB�JB�JB�OB�RB�RB�JB�JB�JB�HB�IB�IB�IB�KB�JB�EB�CB�CB�=B�:B�6B�7B�3B�+B�+B�,B�,B�-B�+B�+B�,B�,B�-B�-B�+B�,B�,B�2B�4B�7B�2B�2B�2B�2B�2B�2B�/B�,B�,B�-B�,B�,B�,B�-B�-B�,B�,B�0B�-B�#B�(B�%B�%B� B� B�&B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B��B��B�ZB��B��B��B�B�B�B�B�B��B��B��B��B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.79 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344042019040513440420190405134404  AO  ARCAADJP                                                                    20181006004031    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004031  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004031  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134404  IP                  G�O�G�O�G�O�                