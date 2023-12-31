CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:26Z creation      
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
_FillValue                    ��   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %10.3f     FORTRAN_format        F10.3      
resolution        :�o     �  ��   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  0  ��   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    ��   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    ��   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    ��   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  ,  ��   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    ��   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    ��   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    ��   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    �    HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  �   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    �D   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    �T   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    �X   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         �h   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         �l   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        �p   HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004026  20190405134359  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               8A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�6�I��1   @�6�����@N����F�A�I�^5?1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    8A   B   B   @33@�  @�33A��A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B��B   B(  B0  B8ffB@  BH  BP  BX  B`  BhffBpffBx  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�33B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4�fD5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� DtٚDy��D���D�,�D���D��fD���D�0 D���D��3D��fD�0 D��3D���D� D�P D�s3D���D�fD�9�D�l�D�� 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @1G�@�
=@�=pA	�A'�AG�Ag�A�A�A�A�A�A�A�A�B�HB	�HB�HBz�B!�HB)�HB1�HB:G�BA�HBI�HBQ�HBY�HBa�HBjG�BrG�By�HB��B��B��qB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BнqB��B��B��B��B��B��B�#�B��B��B��B��C xRCxRCxRCxRCxRC
xRCxRCxRCxRCxRCxRCxRCxRCxRCxRCxRC xRC"xRC$xRC&xRC(xRC*xRC,xRC.xRC0xRC2xRC4xRC6xRC8xRC:xRC<xRC>xRC@xRCBxRCDxRCFxRCHxRCJxRCLxRCNxRCPxRCRxRCTxRCVxRCXxRCZxRC\xRC^xRC`xRCbxRCdxRCfxRChxRCjxRClxRCnxRCpxRCrxRCtxRCvxRCxxRCzxRC|xRC~xRC�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)C�<)D D �DD�DD�DD�DD�DD�DD�DD�DD�D	D	�D
D
�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�DD�D D �D!D!�D"D"�D#D#�D$D$�D%D%�D&D&�D'D'�D(D(�D)D)�D*D*�D+D+�D,D,�D-D-�D.D.�D/D/�D0D0�D1D1�D2D2�D3D3�D4D4�zD5D5�D6D6�D7D7�D8D8�D9D9�D:D:�D;D;�D<D<�D=D=�D>D>�D?D?�D@D@�DADA�DBDB�DCDC�DDDD�DEDE�DFDF�DGDG�DHDH�DIDI�DJDJ�DKDK�DLDL�DMDM�DNDN�DODO�DPDP�DQDQ�DRDR�DSDS�DTDT�DUDU�DVDV�DWDW�DXDX�DYDY�DZDZ�D[D[�D\D\�D]D]�D^D^�D_D_�D`D`�DaDa�DbDb�DcDc�DdDd�DeDe�DfDf�DgDg�DhDh�DiDi�DjDj�DkDk�DlDl�DmDm�DnDn�DoDo�DpDp�DqDq�DrDr�DsDs�DtDt�Dt��Dy��D��D�;�D���D��pD��D�?
D���D��=D�pD�?
D��=D���D�
D�_
Dڂ=D���D�pD�H�D�{�D��
1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��`@��@��@��`@��`@���@��@���@�%@�%@���@���@�%@�V@�V@�%@�V@��@�V@��@�&�@��@��@�&�@��@�/@�&�@�/@�/@�&�@�&�@�&�@�&�@�&�@�/@�/@�/@�/@�/@�/@�/@�7L@�7L@�7L@�7L@�&�@�bN@��@�&�@��T@���@��@���@�`B@�%@�Ĝ@� �@���@�^5@��D@��@��F@��P@�l�@�\)@�33@�33@��@��!@�=q@��^@�?}@�&�@��@�r�@��@��@��`@��9@���@�+@�v�@���@���@�hs@���@�5?@���@�p�@��@���@� �@��
@��@��
@�33@�
=@�7L@��j@���@�/@��9@�Q�@� �@�
=@��@��T@���@���@�ff@�hs@��@��@��7@�E�@�@�V@�?}@�hs@�hs@��7@���@���@��-@���@�7L@��@�r�@��w@��P@�o@�
=@�;d@�
=@�v�@��@���@���@��h@��h@�X@��@���@��/@���@��@�z�@�A�@��@�  @�w@+@~ȴ@~�+@~E�@~$�@~@}�@}��@}`B@}/@|��@|�@|Z@|1@{�
@{t�@z��@zn�@y�@y�^@y��@y��@y%@xĜ@x�@xbN@x �@w�;@w�w@wl�@v��@vȴ@v��@vff@vE�@v@u��@u?}@t�/@t��@tI�@s��@sƨ@s��@sdZ@so@s33@r�@r��@r^5@r�@q��@q�#@q�^@q��@qx�@qhs@qhs@qX@q&�@p��@p�u@pbN@pb@o�@o��@o�w@o�w@o�@o�@o��@o|�@ol�@oK�@o;d@oK�@o+@n��@nȴ@n�R@n�+@nE�@n5?@n5?@n$�@n@n@m�@m�@m�T@m��@m�h@mp�@m`B@m?}@m/@m�@m�@m�@l��@l�/@l�@l�D@l�D@lz�@lj@lI�@l9X@l�@k�m@k�
@k�F@kt�@kC�@kC�@k33@k"�@k"�@ko@j�@j��@j��@j�\@j~�@jn�@jM�@j=q@j�@jJ@i�@i�@i�#@i�#@i�^@i�7@i&�@i%@h�`@h�`@h��@h��@hĜ@h�9@h�`@i%@i%@hĜ@h�u@h�u@h��@h�u@h�u@h�@hbN@hQ�@hQ�@hQ�@hQ�@hQ�@hA�@hA�@h1'@h �@h �@hb@h  @h  @g�@g�;@g��@g�w@g�;@g��@g��@g�w@g�@g��@g��@g|�@gl�@gl�@gl�@gK�@g+@g�@g
=@f��@f��@g
=@g
=@f�y@f�R@f�y@g
=@g
=@f��@f�y@f�@f�R@f�R@f��@f�+@fff@f$�@f@f@e��@e�-@e�-@e?}@e�@e/@d�@d�@d�D@dz�@dz�@dj@dZ@d(�@d(�@d(�@d�@d1@c�m@cƨ@c��@c�@cS�@cC�@c33@c33@c33@c33@c"�@co@co@co@b��@b^5@b^5@b=q@b-@b�@bJ@a�@a��@a�7@ax�@ahs@ahs@ax�@ax�@ax�@a�7@ax�@ax�@ahs@ahs@aX@aX@a�@a�@a%@`��@`��@`�9@`�@`bN@`bN@`A�@`1'@`1'@`1'@`b@_�@_��@_��@_�w@_��@_��@_�;@`  @`A�@`r�@`1'@` �@_�;@_�w@_�@_\)@_�@_
=@_
=@^��@^ȴ@^��@^�+@^v�@^v�@^v�@^v�@^ff@^V@^V@^E�@^$�@^$�@^$�@^5?@^�+@^�+@^�+@^v�@^ff@^ff@^V@^V@^5?@^$�@^@]�@]�@]��@]��@]��@]�-@]�h@]p�@]�@]`B@]`B@]p�@]/@\�@\�@\�@\��@\�/@[��@Z�@ZJ@["�@Y%@X�@W�P@W�@Z�@_�w@eO�@h1'@i�@i�#@hb@e��@a��@_�@]�@Z�!@XA�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 @��`@��@��@��`@��`@���@��@���@�%@�%@���@���@�%@�V@�V@�%@�V@��@�V@��@�&�@��@��@�&�@��@�/@�&�@�/@�/@�&�@�&�@�&�@�&�@�&�@�/@�/@�/@�/@�/@�/@�/@�7L@�7L@�7L@�7L@�&�@�bN@��@�&�@��T@���@��@���@�`B@�%@�Ĝ@� �@���@�^5@��D@��@��F@��P@�l�@�\)@�33@�33@��@��!@�=q@��^@�?}@�&�@��@�r�@��@��@��`@��9@���@�+@�v�@���@���@�hs@���@�5?@���@�p�@��@���@� �@��
@��@��
@�33@�
=@�7L@��j@���@�/@��9@�Q�@� �@�
=@��@��T@���@���@�ff@�hs@��@��@��7@�E�@�@�V@�?}@�hs@�hs@��7@���@���@��-@���@�7L@��@�r�@��w@��P@�o@�
=@�;d@�
=@�v�@��@���@���@��h@��h@�X@��@���@��/@���@��@�z�@�A�@��@�  @�w@+@~ȴ@~�+@~E�@~$�@~@}�@}��@}`B@}/@|��@|�@|Z@|1@{�
@{t�@z��@zn�@y�@y�^@y��@y��@y%@xĜ@x�@xbN@x �@w�;@w�w@wl�@v��@vȴ@v��@vff@vE�@v@u��@u?}@t�/@t��@tI�@s��@sƨ@s��@sdZ@so@s33@r�@r��@r^5@r�@q��@q�#@q�^@q��@qx�@qhs@qhs@qX@q&�@p��@p�u@pbN@pb@o�@o��@o�w@o�w@o�@o�@o��@o|�@ol�@oK�@o;d@oK�@o+@n��@nȴ@n�R@n�+@nE�@n5?@n5?@n$�@n@n@m�@m�@m�T@m��@m�h@mp�@m`B@m?}@m/@m�@m�@m�@l��@l�/@l�@l�D@l�D@lz�@lj@lI�@l9X@l�@k�m@k�
@k�F@kt�@kC�@kC�@k33@k"�@k"�@ko@j�@j��@j��@j�\@j~�@jn�@jM�@j=q@j�@jJ@i�@i�@i�#@i�#@i�^@i�7@i&�@i%@h�`@h�`@h��@h��@hĜ@h�9@h�`@i%@i%@hĜ@h�u@h�u@h��@h�u@h�u@h�@hbN@hQ�@hQ�@hQ�@hQ�@hQ�@hA�@hA�@h1'@h �@h �@hb@h  @h  @g�@g�;@g��@g�w@g�;@g��@g��@g�w@g�@g��@g��@g|�@gl�@gl�@gl�@gK�@g+@g�@g
=@f��@f��@g
=@g
=@f�y@f�R@f�y@g
=@g
=@f��@f�y@f�@f�R@f�R@f��@f�+@fff@f$�@f@f@e��@e�-@e�-@e?}@e�@e/@d�@d�@d�D@dz�@dz�@dj@dZ@d(�@d(�@d(�@d�@d1@c�m@cƨ@c��@c�@cS�@cC�@c33@c33@c33@c33@c"�@co@co@co@b��@b^5@b^5@b=q@b-@b�@bJ@a�@a��@a�7@ax�@ahs@ahs@ax�@ax�@ax�@a�7@ax�@ax�@ahs@ahs@aX@aX@a�@a�@a%@`��@`��@`�9@`�@`bN@`bN@`A�@`1'@`1'@`1'@`b@_�@_��@_��@_�w@_��@_��@_�;@`  @`A�@`r�@`1'@` �@_�;@_�w@_�@_\)@_�@_
=@_
=@^��@^ȴ@^��@^�+@^v�@^v�@^v�@^v�@^ff@^V@^V@^E�@^$�@^$�@^$�@^5?@^�+@^�+@^�+@^v�@^ff@^ff@^V@^V@^5?@^$�@^@]�@]�@]��@]��@]��@]�-@]�h@]p�@]�@]`B@]`B@]p�@]/@\�@\�@\�@\��G�O�@[��@Z�@ZJ@["�@Y%@X�@W�P@W�@Z�@_�w@eO�@h1'@i�@i�#@hb@e��@a��@_�@]�@Z�!@XA�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B(�B(�B)�B(�B(�B(�B)�B)�B(�B)�B)�B)�B)�B)�B)�B(�B(�B(�B(�B)�B(�B(�B(�B(�B(�B(�B(�B1'B@�BL�BcTB�B�=B�VB�PB�VB�PB�PB�JB�DB�=B�7B�=B�DB�PB�VB�VB�bB�bB�bB�\B�\B�hB�oB�oB�oB�hB�oB��B��B��B��B��B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B��B��B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�JB�JB�JB�JB�JB�DB�DB�DB�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�+B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B~�B~�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B{�B{�B{�B{�B{�B{�Bx�Bx�Bv�Bz�Bu�Bs�Bq�Bq�B}�B�DB��B�B�-B�XB�dB�jB�XB�XB�dB�dB�j1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B)�B(�B(�B)�B(�B(�B(�B)�B)�B(�B)�B)�B)�B)�B)�B)�B(�B(�B(�B(�B)�B(�B(�B(�B(�B(�B(�B(�B0�B@CBL�BcB��B��B�B�B�B�B�B�
B�B��B��B��B� B�B�B�B�B� B�B�B�B�%B�+B�+B�+B�&B�-B�LB�KB�DB�FB�>B�:B�:B�AB�PB�hB�kB�eB�eB�_B�\B�^B�]B�qB�jB�kB�YB�RB�^B�hB�_B�\B�^B�LB�<B�?B�DB�YB�YB�FB�GB�MB�[B�qB�hB�dB�kB�oB�rB�zB�~B��B��B��B��B��B��B�{B�uB�pB�sB�~B�yB�tB�wB�tB�tB�sB�sB�tB�uB�uB�uB�|B�~B�{B�|B�zB�vB�sB�sB�uB�sB�sB�uB�sB�sB�uB�vB�tB�vB�wB�tB�pB�mB�qB�jB�kB�jB�kB�kB�kB�kB�iB�iB�iB�iB�gB�iB�gB�eB�eB�cB�fB�eB�eB�fB�_B�^B�[B�^B�VB�XB�XB�VB�VB�VB�XB�ZB�QB�QB�QB�RB�LB�RB�PB�PB�PB�OB�PB�PB�LB�LB�LB�HB�EB�LB�HB�HB�JB�HB�DB�EB�EB�EB�DB�EB�DB�FB�CB�BB�CB�DB�DB�BB�>B�CB�?B�DB�DB�@B�@B�?B�?B�?B�<B�?B�>B�>B�?B�>B�9B�;B�>B�?B�;B�;B�9B�9B�9B�;B�3B�5B�3B�2B�3B�3B�2B�3B�2B�*B�*B�,B�,B�,B�-B�/B�-B�)B�'B�'B�'B�'B�$B�'B�"B�!B�!B� B�"B� B�!B�!B�!B�!B� B� B�B�B� B�!B� B�B�B� B�!B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B}�B~�B~�B~�B~�B~�B~�B~�B}�B}�B}�B}�B}�B}�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B{�B{�B{�B{�B{�G�O�Bx�Bx�Bv�Bz�Bu�BsxBqjBqiB}�B�B�tB��B��B�B�'B�/B�B�B�(B�(B�,1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.47 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051343592019040513435920190405134359  AO  ARCAADJP                                                                    20181006004026    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004026  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004026  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134359  IP                  G�O�G�O�G�O�                